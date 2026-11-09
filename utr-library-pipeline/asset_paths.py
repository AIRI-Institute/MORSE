"""Public cache paths shared by UTR acquisition and library construction."""
from __future__ import annotations
import os
from pathlib import Path


def cache_root() -> Path:
    default = Path(__file__).resolve().parents[1] / ".cache/morse/assets"
    return Path(os.environ.get("MORSE_ASSETS_DIR", default)).expanduser().resolve()


def asset_directory(root: Path | None = None) -> Path:
    if root is not None:
        return root.expanduser().resolve() / "utr-library-pipeline"
    explicit = os.environ.get("MORSE_UTR_LIBRARY_PIPELINE_FILES")
    return Path(explicit).expanduser().resolve() if explicit else cache_root() / "utr-library-pipeline"
