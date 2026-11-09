"""Maintainer utility: record the reviewed public evidence bundle's checksums.

Do not run this to hide a failed verification. Review an intentional evidence or
implementation change before refreshing its manifest and committing that change.
"""
from pathlib import Path

from .io import HERE, sha256, write_json


def is_public_bundle_file(relative):
    """Keep local unit-test sources outside the public evidence manifest."""
    return not (
        any(part in {"__pycache__", "artifacts", "results", "tests", "test", "testing", "pytest"} for part in relative.parts)
        or relative.name == "bundle_manifest.json"
        or (relative.suffix == ".py" and (relative.name.startswith("test_") or relative.name.endswith("_test.py")))
    )


def main():
    files = {}
    for path in sorted(HERE.rglob("*")):
        relative = path.relative_to(HERE)
        if not path.is_file() or not is_public_bundle_file(relative):
            continue
        files[str(relative)] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    write_json(HERE / "provenance/bundle_manifest.json", {"schema_version": 1, "files": files})
    print(f"Recorded {len(files)} reviewed paper bundle files")


if __name__ == "__main__":
    main()
