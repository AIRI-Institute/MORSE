"""Fetch the fixed upstream Git archive and verify its exact bytes before use."""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import urllib.request

LOCK = Path(__file__).with_name("source-lock.json")


def extract_verified(payload: bytes, destination: Path, upstream: dict) -> None:
    actual = hashlib.sha256(payload).hexdigest()
    if actual != upstream["archive_sha256"]:
        raise ValueError(f"VaxPress archive checksum mismatch: {actual}")
    if destination.exists():
        raise FileExistsError(f"Refusing to replace existing source: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".vaxpress-source-", dir=destination.parent) as temporary:
        stage = Path(temporary) / "source"
        stage.mkdir()
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
            for member in archive.getmembers():
                parts = Path(member.name).parts
                if not parts or parts[0] != upstream["archive_root"] or ".." in parts or Path(member.name).is_absolute():
                    raise ValueError(f"Unexpected pinned archive member: {member.name}")
                target = stage.joinpath(*parts[1:])
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                elif member.isfile():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with archive.extractfile(member) as source:
                        target.write_bytes(source.read())
                    target.chmod(member.mode & 0o777)
                else:
                    raise ValueError(f"Unsupported archive link/device: {member.name}")
        for name, expected in upstream["verified_source_files"].items():
            if hashlib.sha256((stage / name).read_bytes()).hexdigest() != expected:
                raise ValueError(f"Pinned source file differs: {name}")
        (stage / ".morse_source_receipt.json").write_text(json.dumps(upstream, indent=2) + "\n")
        stage.rename(destination)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--archive", type=Path, help="Use already-downloaded archive bytes, with the same immutable checksum")
    args = parser.parse_args(argv)
    upstream = json.loads(LOCK.read_text())["upstream"]
    if args.archive:
        payload = args.archive.read_bytes()
    else:
        request = urllib.request.Request(upstream["archive_url"], headers={"User-Agent": "MORSE-pinned-baseline/1"})
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = response.read()
    extract_verified(payload, args.output.resolve(), upstream)
    print(f"Verified VaxPress {upstream['commit']} at {args.output.resolve()}")


if __name__ == "__main__":
    main()
