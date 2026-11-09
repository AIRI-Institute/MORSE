#!/usr/bin/env python3
"""Verify bundled IRNA panels and optionally archive the Docker image."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path


REPO = "irna"
VERSION = "1.0.0"
SNAPSHOT = "2026-08-21"
DEFAULT_ROOT = Path(os.environ.get("MORSE_ASSETS_DIR", str(Path(__file__).resolve().parents[2] / ".cache/morse/assets")))
PROJECT_DIR = Path(__file__).resolve().parent
EXPECTED_MODELS = {"human": 157, "mouse": 170}


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


def verify_panel(path: Path, species: str) -> dict[str, object]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != EXPECTED_MODELS[species]:
        raise ValueError(
            f"{species} panel has {len(rows)} models, expected {EXPECTED_MODELS[species]}"
        )
    names = [row["scanmir_name"] for row in rows]
    if len(set(names)) != len(names):
        raise ValueError(f"{species} panel contains duplicate scanMiR names")
    fraction_columns = [
        name for name in rows[0] if name.startswith("fraction_total__")
    ]
    if not fraction_columns:
        raise ValueError(f"{species} panel has no expression fractions")
    return {
        "species": species,
        "path": str(path.relative_to(PROJECT_DIR)),
        "models": len(rows),
        "expression_profiles": len(fraction_columns),
        "sha256": sha256(path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--docker-image")
    args = parser.parse_args()

    panel_root = PROJECT_DIR / "panels" / SNAPSHOT
    panels = [
        verify_panel(panel_root / species / "core_975.tsv", species)
        for species in EXPECTED_MODELS
    ]
    sources_path = panel_root / "sources.json"
    sources = json.loads(sources_path.read_text())
    if sources.get("snapshot") != SNAPSHOT:
        raise ValueError("IRNA sources snapshot does not match the bundled panel")

    asset_dir = args.root / f"{REPO}_files"
    for folder in ("docker", "misc"):
        (asset_dir / folder).mkdir(parents=True, exist_ok=True)
    if args.docker_image:
        archive_image(
            args.docker_image,
            asset_dir / "docker" / f"{REPO}-{VERSION}.tar.gz",
        )

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
        "version": VERSION,
        "snapshot": SNAPSHOT,
        "asset_dir": str(asset_dir),
        "scanmir_version": "1.16.0",
        "scanmirdata_version": "1.16.0",
        "panels": panels,
        "source_manifest_sha256": sha256(sources_path),
        "files": files,
    }
    manifest_path = asset_dir / "misc" / "assets_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Verified {len(panels)} panels; manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
