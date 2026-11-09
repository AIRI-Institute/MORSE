#!/usr/bin/env python3
"""Download and verify the pinned GENCODE inputs for UTR library building."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import urllib.request
from pathlib import Path


REPO = "utr-library-pipeline"
from asset_paths import asset_directory

RELEASE = 50
SOURCES = {
    "gencode.v50.annotation.gff3.gz": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gff3.gz",
        "sha256": "2aaf245c91ed00e80920953add6cfaffcccc876dc0aceeb6ca0c86d15875899a",
    },
    "gencode.v50.transcripts.fa.gz": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz",
        "sha256": "5a320f524d73b5793518eb19b118829033713443d0f42af20a67bb31cc06cf56",
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, help="Override MORSE_ASSETS_DIR for this download/verification")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()

    asset_dir = asset_directory(args.root)
    raw_dir = asset_dir / "data" / "raw" / f"gencode_v{RELEASE}"
    if not args.verify_only:
        for name, spec in SOURCES.items():
            fetch(spec["url"], raw_dir / name, spec["sha256"])
    for name, spec in SOURCES.items():
        path = raw_dir / name
        if not path.is_file() or sha256(path) != spec["sha256"]:
            raise ValueError(f"Missing or invalid asset: {path}")

    files = []
    for folder in ("data", "misc"):
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
        "gencode_release": RELEASE,
        "asset_dir": str(asset_dir),
        "sources": SOURCES,
        "files": files,
    }
    manifest_path = asset_dir / "misc" / "assets_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Verified {len(files)} files; manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
