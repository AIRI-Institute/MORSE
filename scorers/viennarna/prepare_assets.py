#!/usr/bin/env python3
"""Prepare the ViennaRNA asset directory and optionally archive its image."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path


REPO = "viennarna"
DEFAULT_ROOT = Path(os.environ.get("MORSE_ASSETS_DIR", str(Path(__file__).resolve().parents[2] / ".cache/morse/assets")))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    args = parser.parse_args()
    asset_dir = args.root / f"{REPO}_files"
    for folder in ("data", "models", "docker", "misc"):
        (asset_dir / folder).mkdir(parents=True, exist_ok=True)
    if args.docker_image:
        archive_image(args.docker_image, asset_dir / "docker" / f"{REPO}-1.0.0.tar.gz")

    files = []
    for folder in ("docker", "misc"):
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
        "asset_dir": str(asset_dir),
        "viennarna_version": "2.7.2",
        "external_data_or_weights": False,
        "files": files,
    }
    manifest_path = asset_dir / "misc" / "assets_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Verified {len(files)} files; manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
