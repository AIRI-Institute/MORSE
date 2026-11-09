"""Release packaging for the optional checksum-addressed raw observations."""
from __future__ import annotations

from pathlib import Path
import shutil
import tarfile

from .io import HERE, read_json, sha256
from .verification import verify


def manage(action,archive,raw_root):
    archive=Path(archive).resolve()
    raw_root=Path(raw_root).resolve()
    manifest=read_json(HERE/"provenance/raw_sources.json")
    if action=="pack":
        verify(raw=True,raw_root=raw_root)
        if archive.exists():raise ValueError(f"Refusing to replace release archive: {archive}")
        archive.parent.mkdir(parents=True,exist_ok=True)
        with tarfile.open(archive,"w:gz",compresslevel=6) as bundle:
            for relative in sorted(manifest["sources"]):
                bundle.add(raw_root/relative,arcname=relative,recursive=False)
    else:
        with tarfile.open(archive,"r:gz") as bundle:
            members=bundle.getmembers()
            if {m.name for m in members}!=set(manifest["sources"]):raise ValueError("Raw release archive file set differs from the frozen manifest")
            for member in members:
                relative=Path(member.name)
                if not member.isfile() or relative.is_absolute() or ".." in relative.parts:
                    raise ValueError(f"Unsafe raw archive member: {member.name}")
                destination=raw_root/relative
                if destination.exists():
                    if sha256(destination)==manifest["sources"][member.name]:continue
                    raise ValueError(f"Refusing to overwrite a different existing raw file: {destination}")
                destination.parent.mkdir(parents=True,exist_ok=True)
                with bundle.extractfile(member) as source,destination.open("xb") as target:shutil.copyfileobj(source,target)
                if sha256(destination)!=manifest["sources"][member.name]:
                    destination.unlink()
                    raise ValueError(f"Raw archive checksum mismatch: {member.name}")
        verify(raw=True,raw_root=raw_root)
    return {"status":"pass","action":action,"archive":str(archive),"archive_sha256":sha256(archive),"archive_bytes":archive.stat().st_size,"raw_files":len(manifest["sources"])}
