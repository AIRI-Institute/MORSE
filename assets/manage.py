#!/usr/bin/env python3
"""Install and verify frozen inference assets, without importing ML libraries."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile
import urllib.request
import zipfile

REPOSITORY = Path(os.environ.get("MORSE_ROOT", str(Path(__file__).resolve().parents[1]))).expanduser().resolve()
REGISTRY = Path(__file__).with_name("registry.json")


def digest(path: Path, algorithm: str = "sha256") -> str:
    result = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def cache_root() -> Path:
    return Path(os.environ.get("MORSE_ASSETS_DIR", str(REPOSITORY / ".cache/morse/assets"))).expanduser().resolve()


def safe_path(root: Path, relative: str) -> Path:
    target = root / relative
    if not target.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Asset path escapes its directory: {relative}")
    return target


def _wheel_files(asset: dict, directory: Path):
    """Yield validated pure-Python wheel members; never import package code."""
    spec = asset["runtime_package"]
    with zipfile.ZipFile(safe_path(directory, spec["wheel"])) as archive:
        wheel_metadata = archive.read(spec["dist_info"] + "/WHEEL").decode("utf-8")
        if "Root-Is-Purelib: true" not in wheel_metadata.splitlines() or "Tag: py3-none-any" not in wheel_metadata.splitlines():
            raise ValueError("Runtime preparation requires the registered pure-Python wheel")
        # Only explicitly registered, byte-identical aliases may share a physical
        # file. No general case folding or filename normalization is permitted.
        for alias, details in spec.get("identical_file_aliases", {}).items():
            canonical = details["canonical_path"]
            if alias == canonical or details["sha256"] != details["canonical_sha256"]:
                raise ValueError("Registered wheel alias must preserve identical bytes")
            for relative, expected in ((alias, details["sha256"]),
                                       (canonical, details["canonical_sha256"])):
                result = hashlib.sha256()
                with archive.open(relative) as stream:
                    for block in iter(lambda: stream.read(1024 * 1024), b""):
                        result.update(block)
                if result.hexdigest() != expected:
                    raise ValueError(f"Registered wheel alias checksum mismatch: {relative}")
        seen = set()
        for item in archive.infolist():
            relative = item.filename
            if (not relative or "\\" in relative or relative.startswith("/")
                    or ".." in Path(relative).parts or stat.S_ISLNK(item.external_attr >> 16)):
                raise ValueError(f"Unsafe wheel member: {relative}")
            if item.is_dir():
                continue
            if relative in seen or Path(relative).parts[0] not in {spec["package"], spec["dist_info"]}:
                raise ValueError(f"Unexpected or duplicate wheel member: {relative}")
            seen.add(relative)
            yield archive, item


def verify_runtime_package(asset: dict, directory: Path) -> dict:
    spec = asset["runtime_package"]
    package_root = safe_path(directory, spec["directory"])
    expected = {}
    physical_files = set()
    for archive, item in _wheel_files(asset, directory):
        result = hashlib.sha256()
        with archive.open(item) as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                result.update(block)
        expected[item.filename] = result.hexdigest()
        physical = spec.get("identical_file_aliases", {}).get(item.filename, {}).get("canonical_path", item.filename)
        physical_files.add(physical)
        path = safe_path(package_root, physical)
        if path.is_symlink() or not path.is_file() or digest(path) != expected[item.filename]:
            raise ValueError(f"Missing or invalid prepared runtime file: {item.filename}")
    actual = {str(path.relative_to(package_root)) for path in package_root.rglob("*") if path.is_file()}
    if actual != physical_files:
        raise ValueError("Prepared runtime contains files outside the registered wheel")
    return {"status": "verified", "directory": str(package_root), "files_sha256": expected,
            "identical_file_aliases": spec.get("identical_file_aliases", {})}


def prepare_runtime_package(asset: dict, directory: Path) -> None:
    spec = asset.get("runtime_package")
    if not spec:
        return
    destination = safe_path(directory, spec["directory"])
    if destination.exists():
        verify_runtime_package(asset, directory)
        return
    with tempfile.TemporaryDirectory(prefix=".runtime-package-", dir=directory) as temporary:
        stage = Path(temporary) / "package"
        stage.mkdir()
        for archive, item in _wheel_files(asset, directory):
            if item.filename in spec.get("identical_file_aliases", {}):
                continue
            target = safe_path(stage, item.filename)
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(item) as incoming, target.open("xb") as outgoing:
                shutil.copyfileobj(incoming, outgoing, length=1024 * 1024)
        # The raw wheel digest is verified before this function is called.
        # Expose only a complete tree; existing prepared trees are never replaced.
        stage.rename(destination)
    verify_runtime_package(asset, directory)


def verify_asset(asset: dict, directory: Path, *, dependency_root: Path | None = None,
                 require_runtime: bool = False) -> dict:
    verified = []
    for item in asset.get("files", []):
        path = safe_path(directory, item["path"])
        if not path.is_file():
            raise FileNotFoundError(f"Missing {asset['id']}: {path}")
        algorithm = item.get("algorithm", "sha256")
        actual = digest(path, algorithm)
        if actual != item[algorithm]:
            raise ValueError(f"Checksum mismatch for {asset['id']}/{item['path']}: expected {item[algorithm]}, got {actual}")
        verified.append({"path": item["path"], "sha256": digest(path) if algorithm != "sha256" else actual,
                         "bytes": path.stat().st_size})
    if asset.get("kind") == "ribonn-ensemble":
        manifest = json.loads((directory / "ensemble.json").read_text())
        if manifest.get("record") != "17258709" or manifest.get("n_models") != 50 or manifest.get("n_targets") != 78:
            raise ValueError("Expected the frozen official human 50-model / 78-target ensemble")
        if len(manifest.get("models", [])) != 50:
            raise ValueError("RiboNN manifest must list 50 model files")
        items = [("human/runs.csv", manifest["runs_csv_sha256"])]
        items += [(item["path"], item["sha256"]) for item in manifest["models"]]
        for relative, expected in items:
            path = safe_path(directory, relative)
            if not path.is_file() or digest(path) != expected:
                raise ValueError(f"Missing or invalid RiboNN file: {relative}")
            verified.append({"path": relative, "sha256": expected, "bytes": path.stat().st_size})
    result = {"id": asset["id"], "status": "verified", "directory": str(directory), "files": verified}
    if asset.get("runtime_package"):
        package_root = safe_path(directory, asset["runtime_package"]["directory"])
        if package_root.exists() or require_runtime:
            result["runtime_package"] = verify_runtime_package(asset, directory)
        else:
            # Historical raw-asset verification remains valid without preparing a
            # new container layout. --runtime requires the complete mounted tree.
            result["runtime_package"] = {"status": "not_prepared", "directory": str(package_root)}
    if asset.get("depends_on"):
        registry = {item["id"]: item for item in json.loads(REGISTRY.read_text())["assets"]}
        parent = dependency_root or directory.parent
        result["dependencies"] = [verify_asset(registry[name], parent / name) for name in asset["depends_on"]]
    if require_runtime and asset.get("runtime_dependencies"):
        registry = {item["id"]: item for item in json.loads(REGISTRY.read_text())["assets"]}
        parent = dependency_root or directory.parent
        result["runtime_dependencies"] = [verify_asset(registry[name], parent / name, require_runtime=True)
                                          for name in asset["runtime_dependencies"]]
    return result


def fetch(item: dict, output: Path) -> None:
    if item.get("repository_path"):
        source = safe_path(REPOSITORY, item["repository_path"])
        if not source.is_file():
            raise FileNotFoundError("Bundled profile installation requires the full MORSE checkout; "
                                    "set MORSE_ROOT to its root directory")
        shutil.copyfile(source, output)
        return
    url = item.get("url")
    if not url:
        raise ValueError(f"No public download is registered for {item['path']}; use --source with an authorized ready checkpoint")
    output.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "MORSE-inference-assets/1.0"})
    with urllib.request.urlopen(request, timeout=60) as incoming, output.open("wb") as outgoing:
        shutil.copyfileobj(incoming, outgoing, length=1024 * 1024)
    if item.get("transform") == "postponed_annotations":
        raw = output.read_bytes()
        expected = item["download_sha256"]
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError(f"Upstream code checksum mismatch: {item['path']}")
        output.write_bytes(b"from __future__ import annotations\n\n" + raw)


def install_asset(asset: dict, root: Path, source: Path | None) -> dict:
    destination = root / asset["id"]
    if destination.exists():
        verified = verify_asset(asset, destination)
        if verified.get("runtime_package", {}).get("status") == "not_prepared":
            prepare_runtime_package(asset, destination)
            return verify_asset(asset, destination, require_runtime=True)
        return verified
    root.mkdir(parents=True, exist_ok=True)
    # Installation becomes visible only after every registered digest passes.
    with tempfile.TemporaryDirectory(prefix=f".{asset['id']}-", dir=root) as temporary:
        stage = Path(temporary) / "asset"
        stage.mkdir()
        if asset.get("kind") == "ribonn-ensemble":
            if source:
                if not source.is_dir():
                    raise ValueError("RiboNN --source must be the prepared ensemble directory")
                for relative in ["weights.zip", "ensemble.json", "human/runs.csv"]:
                    incoming = safe_path(source, relative)
                    if not incoming.is_file():
                        raise FileNotFoundError(incoming)
                    target = safe_path(stage, relative)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(incoming, target)
                manifest = json.loads((stage / "ensemble.json").read_text())
                for model in manifest.get("models", []):
                    relative = model["path"]
                    target = safe_path(stage, relative)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(safe_path(source, relative), target)
            else:
                subprocess.run([sys.executable, str(REPOSITORY / "scorers/ribonn/prepare_assets.py"),
                                "--directory", str(stage)], check=True)
        else:
            for item in asset.get("files", []):
                target = safe_path(stage, item["path"])
                target.parent.mkdir(parents=True, exist_ok=True)
                if source:
                    candidates = [source] if source.is_file() and len(asset["files"]) == 1 else [
                        safe_path(source, name) for name in [item["path"], *item.get("source_names", [])]]
                    incoming = next((path for path in candidates if path.is_file()), None)
                    if incoming is None:
                        raise FileNotFoundError(f"Cannot locate {item['path']} in {source}")
                    shutil.copyfile(incoming, target)
                    if item.get("transform") == "postponed_annotations" and digest(target) == item["download_sha256"]:
                        target.write_bytes(b"from __future__ import annotations\n\n" + target.read_bytes())
                else:
                    fetch(item, target)
        verify_asset(asset, stage, dependency_root=root)
        prepare_runtime_package(asset, stage)
        receipt = verify_asset(asset, stage, dependency_root=root,
                               require_runtime=bool(asset.get("runtime_package")))
        if "runtime_package" in receipt:
            receipt["runtime_package"].pop("directory", None)
        receipt["source_kind"] = ("local_ready_asset" if source else "registered_bundled_source"
                                  if all(item.get("repository_path") for item in asset.get("files", []))
                                  else "registered_public_source")
        receipt.pop("directory", None)
        (stage / "morse-installation.json").write_text(json.dumps(receipt, indent=2) + "\n")
        stage.rename(destination)
    return verify_asset(asset, destination)


def runtime_paths(asset: dict, directory: Path, root: Path) -> dict:
    def expand(value):
        return value.replace("{asset}", str(directory)).replace("{root}", str(root))
    env = {key: expand(value)
           for key, value in asset.get("runtime_env", {}).items()}
    mounts = [{**mount, "source": expand(mount["source"])} for mount in asset.get("runtime_mounts", [])]
    return {"id": asset["id"], "directory": str(directory), "runtime_env": env,
            "depends_on": asset.get("depends_on", []),
            "runtime_dependencies": asset.get("runtime_dependencies", []),
            "runtime_mounts": mounts,
            "runtime_mounts_present": all(Path(mount["source"]).exists() for mount in mounts),
            "present": all((directory / item["path"]).is_file() for item in asset.get("files", []))
                       and all((root / dependency).is_dir() for dependency in asset.get("depends_on", []))}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["list", "install", "verify", "paths"])
    parser.add_argument("--asset", help="Registry asset ID; omit for list, paths, or verify all")
    parser.add_argument("--source", type=Path, help="Authorized ready checkpoint file or prepared asset directory")
    parser.add_argument("--root", type=Path, default=cache_root(), help="Cache root (default MORSE_ASSETS_DIR)")
    parser.add_argument("--json", action="store_true", help="Machine-readable result")
    parser.add_argument("--runtime", action="store_true", help="For verify, require the prepared package and external runtime dependencies too")
    args = parser.parse_args(argv)
    registry = json.loads(REGISTRY.read_text())
    assets = {item["id"]: item for item in registry["assets"]}
    if args.asset and args.asset not in assets:
        parser.error(f"Unknown asset {args.asset!r}; run list")
    if args.command == "install" and not args.asset:
        parser.error("install requires --asset; use list to inspect source and license status")
    if args.source and args.command != "install":
        parser.error("--source is valid only for install")
    if args.runtime and args.command != "verify":
        parser.error("--runtime is valid only for verify")
    root = args.root.expanduser().resolve()
    chosen = [assets[args.asset]] if args.asset else list(assets.values())
    result = []
    try:
        for asset in chosen:
            directory = root / asset["id"]
            if args.command == "list":
                result.append({"id": asset["id"], "paper": asset.get("paper", False),
                               "license_status": asset["license_status"],
                               "download_available": all(item.get("url") or item.get("repository_path") for item in asset.get("files", [])),
                               "installed": directory.is_dir()})
            elif args.command == "paths":
                result.append(runtime_paths(asset, directory, root))
            elif args.command == "install":
                for dependency in [*asset.get("depends_on", []), *asset.get("runtime_dependencies", [])]:
                    install_asset(assets[dependency], root, None)
                result.append(install_asset(asset, root, args.source.expanduser().resolve() if args.source else None))
            else:
                result.append(verify_asset(asset, directory, require_runtime=args.runtime))
    except (ValueError, OSError, KeyError, zipfile.BadZipFile, subprocess.CalledProcessError) as error:
        print(f"Asset error: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result[0] if args.asset else result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
