# syntax=docker/dockerfile:1.7
FROM pytorch/pytorch:2.12.1-cuda13.0-cudnn9-runtime@sha256:72f863fa1fe13d5d87a72d00db2c85fb2d43409ee08dd26bc469de4c8a28b427

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_BREAK_SYSTEM_PACKAGES=1

RUN python -m pip install "pip==26.1.2" "setuptools==81.0.0" "wheel==0.47.0" \
    && python -m pip install \
      "numpy==2.4.6" \
      "multimolecule==0.2.0" \
    && mkdir -p /usr/share/morse/sources/multimolecule \
    && python -m pip download --no-deps \
      --dest /usr/share/morse/sources/multimolecule "multimolecule==0.2.0" \
    && python -m pip download --no-deps --no-binary=multimolecule \
      --dest /usr/share/morse/sources/multimolecule "multimolecule==0.2.0" \
    && python -m pip freeze --all > /usr/share/morse/sources/python-packages.txt

# The wheel embeds AGPL text in METADATA; its original standalone notices are
# in the exact sdist. Preserve both with verified hashes beside the runtime.
RUN python - <<'PY'
from importlib.metadata import distribution
from pathlib import Path
import hashlib
import json
import tarfile
import zipfile

dist = distribution("multimolecule")
destination = Path("/usr/share/morse/notices/multimolecule")
sources = Path("/usr/share/morse/sources/multimolecule")
expected_archives = {
    "multimolecule-0.2.0-py3-none-any.whl": "a3476171e810b7ce76d9b066b4a061dbaca8ea326673ab7df740c5cc88b1c272",
    "multimolecule-0.2.0.tar.gz": "b4a1cb6b952d4698432f46c7e1e1251a2cdf2ff2ed5c8cf8bf5e94a010198c0f",
}
archives = [sources / name for name in expected_archives]
weight_suffixes = (".pt", ".pth", ".ckpt", ".safetensors", ".onnx", ".h5", ".hdf5", ".bin", ".p", ".pkl", ".pickle")
for archive in archives:
    if not archive.is_file() or hashlib.sha256(archive.read_bytes()).hexdigest() != expected_archives[archive.name]:
        raise RuntimeError("Exact MultiMolecule archive SHA-256 mismatch: " + archive.name)
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as handle:
            names = handle.namelist()
    else:
        with tarfile.open(archive) as handle:
            names = handle.getnames()
    if any(name.lower().endswith(weight_suffixes) for name in names):
        raise RuntimeError("Unexpected checkpoint in MultiMolecule source/wheel: " + archive.name)

metadata = dist.read_text("METADATA").encode("utf-8")
metadata_sha256 = "1500df703ca5333e2e7dcc9e7ee3ab14b9a43b130acec24ad790cb7d5087189b"
if dist.version != "0.2.0" or hashlib.sha256(metadata).hexdigest() != metadata_sha256:
    raise RuntimeError("Installed MultiMolecule metadata differs from the exact licensed release")
destination.mkdir(parents=True, exist_ok=True)
(destination / "METADATA").write_bytes(metadata)
notices = {
    "license.md": "e91ab743b5267586856dd817d1fc5933ecea9432c864cd08366152ab647cd055",
    "license.header": "93fa1f0441a81b392e56dd0d3bb06f81cd5dc3fff9c908777e54b68c4bd39438",
}
with tarfile.open(sources / "multimolecule-0.2.0.tar.gz") as archive:
    for name, expected in notices.items():
        with archive.extractfile("multimolecule-0.2.0/" + name) as source:
            contents = source.read()
        if hashlib.sha256(contents).hexdigest() != expected:
            raise RuntimeError("Original MultiMolecule license notice SHA-256 mismatch: " + name)
        (destination / name).write_bytes(contents)
(destination / "provenance.json").write_text(json.dumps({
    "package": "multimolecule", "version": "0.2.0", "license": "AGPL-3.0-or-later",
    "notice_source": "Exact PyPI 0.2.0 source distribution; full license also embedded in wheel METADATA",
    "archives_sha256": expected_archives, "notices_sha256": notices,
    "installed_metadata_sha256": metadata_sha256,
}, indent=2) + "\n")
with (sources / "SHA256SUMS").open("w") as output:
    for archive in archives:
        output.write(hashlib.sha256(archive.read_bytes()).hexdigest() + "  " + archive.name + "\n")
PY

WORKDIR /app
COPY server.py /app/server.py
COPY runtime_model.py /app/runtime_model.py
COPY model_assets /app/model_assets
COPY notices /usr/share/morse/notices
COPY Dockerfile /usr/share/morse/sources/Dockerfile

# The named context is the verified rinalmo-imrnads asset directory only.
# Copy the ready checkpoint directly; no fitting data or repository cache is used.
COPY --from=assets best_model.pt /models/best_model.pt
RUN python - <<'PY_CHECKPOINT'
from pathlib import Path
import hashlib
import json

path = Path("/models/best_model.pt")
expected = "923d1e9f843d12fe29a9515a95b6a32eae286017790dfeef1c82215446958f0a"
digest = hashlib.sha256()
with path.open("rb") as source:
    for block in iter(lambda: source.read(1024 * 1024), b""):
        digest.update(block)
if digest.hexdigest() != expected:
    raise RuntimeError("RiNALMo IMrnaDS checkpoint SHA-256 mismatch")
Path("/usr/share/morse/sources/checkpoint-manifest.json").write_text(json.dumps({
    "asset": "rinalmo-imrnads", "path": str(path), "sha256": expected,
    "size_bytes": path.stat().st_size,
}, indent=2) + "\n")
PY_CHECKPOINT

ENV RINALMO_CHECKPOINT=/models/best_model.pt \
    MORSE_EXPECTED_CHECKPOINT_SHA256=923d1e9f843d12fe29a9515a95b6a32eae286017790dfeef1c82215446958f0a

LABEL org.opencontainers.image.source="https://github.com/AIRI-Institute/MORSE" \
      org.opencontainers.image.version="2.1.1"

ENTRYPOINT ["python", "/app/server.py"]
