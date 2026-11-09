#!/usr/bin/env python3
"""Download and verify the upstream RiboDecode inference assets."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import urllib.request
from pathlib import Path


REPO = "ribodecode"
DEFAULT_ROOT = Path(os.environ.get("MORSE_ASSETS_DIR", str(Path(__file__).resolve().parents[2] / ".cache/morse/assets")))
UPSTREAM_COMMIT = "fcf2d50decf99ebaf0ed2cd58d1dd9976a820470"
FILES = {
    "TranslationModel-1.1.0-py3-none-any.whl": {
        "url": "https://drive.usercontent.google.com/download?id=19Yyl8uUbjQn0TAA4E52ZUGBd_pb_RyW8&export=download&confirm=t",
        "sha256": "22d4b398a399878f7573074fd1c381b2f3894bd3a62ede9cd27989dc47ac276a",
    },
    "env_file.csv": {
        "url": f"https://raw.githubusercontent.com/wangfanfff/RiboDecode/{UPSTREAM_COMMIT}/env_file.csv",
        "sha256": "9a85da49ff6d5de94cd920688d196b49af14b21e4e3fdfdda10e356edd277062",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, path: Path, expected: str) -> None:
    if path.is_file() and sha256(path) == expected:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "morse-assets/1.0"})
    with urllib.request.urlopen(request) as response, temporary.open("wb") as output:
        shutil.copyfileobj(response, output)
    actual = sha256(temporary)
    if actual != expected:
        temporary.unlink(missing_ok=True)
        raise ValueError(f"SHA-256 mismatch for {url}: {actual}")
    temporary.replace(path)


def archive_image(image: str, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".part")
    temporary.unlink(missing_ok=True)
    try:
        with gzip.open(temporary, "wb", compresslevel=6) as handle:
            process = subprocess.Popen(
                ["docker", "image", "save", image], stdout=subprocess.PIPE
            )
            assert process.stdout is not None
            shutil.copyfileobj(process.stdout, handle, length=1024 * 1024)
            process.stdout.close()
            if process.wait() != 0:
                raise subprocess.CalledProcessError(process.returncode, process.args)
        temporary.replace(output)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--docker-image")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()

    asset_dir = args.root / f"{REPO}_files"
    data_dir = asset_dir / "data"
    if not args.verify_only:
        for name, spec in FILES.items():
            fetch(spec["url"], data_dir / name, spec["sha256"])
        if args.docker_image:
            archive_image(args.docker_image, asset_dir / "docker" / f"{REPO}-1.0.0.tar.gz")

    for name, spec in FILES.items():
        path = data_dir / name
        if not path.is_file() or sha256(path) != spec["sha256"]:
            raise ValueError(f"Missing or invalid asset: {path}")
    if sum(1 for _ in (data_dir / "env_file.csv").open()) != 10553:
        raise ValueError("env_file.csv must contain 10,552 genes plus its header")

    files = []
    for folder in ("data", "docker", "misc"):
        for path in sorted((asset_dir / folder).rglob("*")):
            if path.is_file() and path.name != "assets_manifest.json":
                files.append(
                    {
                        "path": str(path.relative_to(asset_dir)),
                        "bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )
    manifest = {
        "repository": REPO,
        "upstream_commit": UPSTREAM_COMMIT,
        "asset_dir": str(asset_dir),
        "sources": FILES,
        "files": files,
    }
    manifest_path = asset_dir / "misc" / "assets_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Verified {len(files)} files; manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
