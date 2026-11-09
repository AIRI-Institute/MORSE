"""Build-time and run-time proof of the active native baseline implementation."""
from __future__ import annotations
import hashlib
import importlib.metadata
import inspect
import json
from pathlib import Path
import platform


def main():
    lock = json.loads(Path(__file__).with_name("source-lock.json").read_text())
    from vaxpress.evolution_chamber import CDSEvolutionChamber
    import RNA
    versions = {name: importlib.metadata.version(name) for name in lock["packages"]}
    for name, spec in lock["packages"].items():
        if versions[name] != spec["version"]:
            raise ValueError(f"Runtime dependency differs: {name} {versions[name]} != {spec['version']}")
    source_root = Path("/opt/VaxPress")
    for name, expected in lock["upstream"]["verified_source_files"].items():
        if hashlib.sha256((source_root / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Installed upstream snapshot differs: {name}")
    active = Path(inspect.getsourcefile(CDSEvolutionChamber))
    expected = lock["upstream"]["verified_source_files"]["vaxpress/evolution_chamber.py"]
    if hashlib.sha256(active.read_bytes()).hexdigest() != expected:
        raise ValueError("Active CDSEvolutionChamber differs from the pinned upstream source")
    if RNA.__version__ != "2.7.2":
        raise ValueError("The shared baseline requires ViennaRNA 2.7.2")
    print(json.dumps({"status": "verified", "upstream_commit": lock["upstream"]["commit"],
        "upstream_archive_sha256": lock["upstream"]["archive_sha256"], "python": platform.python_version(),
        "active_native_search_source": str(active), "packages": versions,
        "historical_image_identity_claimed": False}, indent=2))


if __name__ == "__main__":
    main()
