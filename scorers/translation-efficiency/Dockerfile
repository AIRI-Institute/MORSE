# syntax=docker/dockerfile:1.7
FROM pytorch/pytorch:2.2.0-cuda12.1-cudnn8-runtime@sha256:e5c32997194edfc6fb0235223d1cd0433056045a468e7dfeaac390a3fcf8ae6d

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_BREAK_SYSTEM_PACKAGES=1

# Official release wheels avoid CUDA compilation on an emulated/no-GPU builder.
# Keep this dependency prefix identical in all three Orthrus Dockerfiles.
RUN <<'SH'
set -eu
python - <<'PY'
from pathlib import Path
import hashlib
import shutil
import sys
import urllib.request
import torch

assert sys.version_info[:2] == (3, 10), "Official wheels require CPython 3.10"
assert torch.__version__.split("+", 1)[0] == "2.2.0", "Pinned PyTorch 2.2.0 required"
assert torch.version.cuda and torch.version.cuda.split(".")[0] == "12", "CUDA 12 PyTorch required"
assert not torch._C._GLIBCXX_USE_CXX11_ABI, "Official wheels require the old C++ ABI"
wheels = (
    ("causal_conv1d-1.2.0.post2+cu122torch2.2cxx11abiFALSE-cp310-cp310-linux_x86_64.whl",
     "https://github.com/Dao-AILab/causal-conv1d/releases/download/v1.2.0.post2/causal_conv1d-1.2.0.post2%2Bcu122torch2.2cxx11abiFALSE-cp310-cp310-linux_x86_64.whl",
     "9fffdd36640e6e52e9c3d2ad6ac5483f109cecc2fdeedfc5767487893f8ccae2"),
    ("mamba_ssm-1.2.0.post1+cu122torch2.2cxx11abiFALSE-cp310-cp310-linux_x86_64.whl",
     "https://github.com/state-spaces/mamba/releases/download/v1.2.0.post1/mamba_ssm-1.2.0.post1%2Bcu122torch2.2cxx11abiFALSE-cp310-cp310-linux_x86_64.whl",
     "b264292652a34fb9dd0ce880a34a4407ba7256a3338388d056769ec29a4581c9"),
)
directory = Path("/tmp/morse-orth-runtimes")
directory.mkdir()
for name, url, expected in wheels:
    path = directory / name
    with urllib.request.urlopen(url, timeout=120) as response, path.open("wb") as output:
        shutil.copyfileobj(response, output, length=1024 * 1024)
    with path.open("rb") as source:
        digest = hashlib.sha256()
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != expected:
        raise RuntimeError("Official wheel SHA-256 mismatch: " + name)
PY
python -m pip install "pip==24.0" "setuptools==69.2.0" "wheel==0.43.0" "ninja==1.11.1.1"
python -m pip install \
    "torch==2.2.0" "triton==2.2.0" "einops==0.7.0" "packaging==23.2" \
    "huggingface-hub==0.21.4" "numpy==1.24.3" "transformers==4.38.1" \
    /tmp/morse-orth-runtimes/*.whl
python - <<'PY'
import torch
import causal_conv1d_cuda
import selective_scan_cuda
from causal_conv1d import causal_conv1d_fn
from mamba_ssm import Mamba

assert torch.__version__.split("+", 1)[0] == "2.2.0"
assert torch.version.cuda and torch.version.cuda.split(".")[0] == "12"
assert not torch._C._GLIBCXX_USE_CXX11_ABI
PY
rm -rf /tmp/morse-orth-runtimes
SH


# Only these registered inference assets enter the named build contexts.
COPY --from=backbone /config.json /orthrus_hf.py /pytorch_model.bin /backbone/
RUN printf '%s\n' \
      'a8b513b3d73356e7ac8a59516ebbf51bc03f635e6329ad32350c289f84c8a288  /backbone/config.json' \
      '35fb82d090003c9bf91be1820d5cc0e22c1092d28d9805951d3f284d728cdbe4  /backbone/orthrus_hf.py' \
      'a188c974493679d9eb5d18b5ecc3148f3634cd785e20f8303d6f313f62f00c08  /backbone/pytorch_model.bin' \
      | sha256sum -c - \
    && mkdir -p /models \
    && ln -s /backbone /models/orthrus-backbone

COPY --from=checkpoint /model_final.pt /models/model_final.pt
RUN printf '%s\n' 'c380316fe282640d563b9254e9ab8031e129f9e768c9ce04c0c7955bb0ac6bab  /models/model_final.pt' | sha256sum -c -

ENV ORTHRUS_CHECKPOINT=/models/model_final.pt \
    ORTHRUS_MODEL_SOURCE=/backbone \
    MORSE_EXPECTED_CHECKPOINT_SHA256=c380316fe282640d563b9254e9ab8031e129f9e768c9ce04c0c7955bb0ac6bab \
    HF_HUB_OFFLINE=1 \
    TRANSFORMERS_OFFLINE=1

WORKDIR /app
COPY model.py server.py portable_orthrus.py server_heldout.py heldout_encoding.py /app/
COPY notices /usr/share/morse/notices
COPY Dockerfile /usr/share/morse/sources/Dockerfile

LABEL org.opencontainers.image.source="https://github.com/AIRI-Institute/MORSE" \
      org.opencontainers.image.version="2.2.1"

ENTRYPOINT ["python", "-u", "/app/server.py"]
