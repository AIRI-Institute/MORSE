# syntax=docker/dockerfile:1.7
FROM pytorch/pytorch:2.2.0-cuda12.1-cudnn8-runtime@sha256:e5c32997194edfc6fb0235223d1cd0433056045a468e7dfeaac390a3fcf8ae6d
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 PIP_NO_CACHE_DIR=1
WORKDIR /app
COPY requirements.txt /app/
RUN python -m pip install -r requirements.txt
COPY encoding.py model.py server.py /app/
COPY third_party /app/third_party
COPY notices /usr/share/morse/notices
COPY Dockerfile /usr/share/morse/sources/Dockerfile
RUN --mount=type=bind,from=assets,target=/morse-assets,ro python - <<'PY'
from pathlib import Path, PurePosixPath
import hashlib
import json
import shutil

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

manifest = json.loads(Path("/usr/share/morse/notices/image-assets.json").read_text())
assert manifest["scorer_version"] == "1.0.1" and manifest["selected_models"] == 50
for item in manifest["inference_source_and_notices"]:
    source = Path("/app/third_party/ribonn") / item["path"]
    if sha256(source) != item["sha256"]:
        raise RuntimeError("RiboNN inference source/notice SHA-256 mismatch: " + item["path"])
for item in manifest["files"]:
    relative = PurePosixPath(item["path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise RuntimeError("Unsafe RiboNN asset path")
    source = Path("/morse-assets") / item["path"]
    if not source.is_file() or source.stat().st_size != item["bytes"] or sha256(source) != item["sha256"]:
        raise RuntimeError("RiboNN ready asset SHA-256 mismatch: " + item["path"])
    target = Path("/models") / item["path"]
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)
ensemble = json.loads(Path("/models/ensemble.json").read_text())
assert ensemble["n_models"] == 50 and len(ensemble["models"]) == 50
assert len(list(Path("/models").rglob("state_dict.pth"))) == 50
PY
ENV RIBONN_MODEL_DIR=/models
LABEL org.opencontainers.image.source="https://github.com/AIRI-Institute/MORSE" \
      org.opencontainers.image.version="1.0.1"

ENTRYPOINT ["python", "-u", "/app/server.py"]

