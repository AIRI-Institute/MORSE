# syntax=docker/dockerfile:1.7
FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu20.04@sha256:4d0748e17319c7afde812a57943b9febe9c20237a7e168f9d18224034dff9d14

ENV DEBIAN_FRONTEND=noninteractive \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility \
    PIP_NO_CACHE_DIR=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH=/opt/ribodecode

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    python3.8 \
    python3.8-distutils \
    python3-pip \
    unzip \
    && rm -rf /var/lib/apt/lists/*

RUN python3.8 -m pip install "pip==24.0" "setuptools==69.2.0" "wheel==0.43.0"

RUN python3.8 -m pip install \
    --extra-index-url https://download.pytorch.org/whl/cu118 \
    "torch==2.0.1+cu118"

RUN python3.8 -m pip install \
    "numpy==1.24.4" \
    "pandas==2.0.3" \
    "p-tqdm==1.4.0"

WORKDIR /app
COPY server.py /app/server.py
COPY notices /usr/share/morse/notices
COPY Dockerfile /usr/share/morse/sources/Dockerfile
RUN --mount=type=bind,from=assets,target=/morse-assets,ro python3.8 - <<'PY'
from pathlib import Path, PurePosixPath
import hashlib
import json
import shutil
import subprocess
import sys

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

manifest = json.loads(Path("/usr/share/morse/notices/image-assets.json").read_text())
assert manifest["scorer_version"] == "2.1.1" and manifest["package_version"] == "1.1.0"
for item in manifest["files"]:
    relative = PurePosixPath(item["path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise RuntimeError("Unsafe TranslationModel asset path")
    source = Path("/morse-assets") / item["path"]
    if not source.is_file() or source.stat().st_size != item["bytes"] or sha256(source) != item["sha256"]:
        raise RuntimeError("TranslationModel release asset SHA-256 mismatch: " + item["path"])
wheel = Path("/morse-assets/TranslationModel-1.1.0-py3-none-any.whl")
subprocess.check_call([sys.executable, "-m", "pip", "install", "--no-deps", "--target", "/opt/ribodecode", str(wheel)])
Path("/assets").mkdir(parents=True, exist_ok=True)
shutil.copyfile("/morse-assets/env_file.csv", "/assets/env_file.csv")
source_root = Path("/usr/share/morse/sources/TranslationModel")
source_root.mkdir(parents=True, exist_ok=True)
for item in manifest["installed_package_files"]:
    installed = Path("/opt/ribodecode") / item["path"]
    if not installed.is_file() or installed.stat().st_size != item["bytes"] or sha256(installed) != item["sha256"]:
        raise RuntimeError("Installed TranslationModel source/weights differ: " + item["path"])
    # Keep all supplied preferred Python/configuration/condition and packaging
    # sources accessible beside the build recipe. The native tensor checkpoint
    # remains at its runnable package path without a second 476 MB copy.
    if item["path"] != "TranslationModel/best_model.p":
        target = source_root / item["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(installed, target)
license_path = Path("/opt/ribodecode/TranslationModel-1.1.0.dist-info/LICENSE.txt")
assert sha256(license_path) == "42baafb5997a6f54e6eb8dc186416e386cd178793d85e089d403e59f926bcd20"
shutil.copyfile("/usr/share/morse/notices/image-assets.json", source_root / "image-assets.json")
(source_root / "ASSET_LOCATIONS.json").write_text(json.dumps({
    "preferred_python_source": "/opt/ribodecode/TranslationModel/inference.py",
    "native_tensor_checkpoint": "/opt/ribodecode/TranslationModel/best_model.p",
    "model_configuration": "/opt/ribodecode/TranslationModel/model_config.json",
    "environment_template": "/assets/env_file.csv",
    "adapter_source": "/app/server.py",
    "installation_recipe": "/usr/share/morse/sources/Dockerfile",
    "license": "/opt/ribodecode/TranslationModel-1.1.0.dist-info/LICENSE.txt",
}, indent=2) + "\n")
PY

LABEL org.opencontainers.image.source="https://github.com/AIRI-Institute/MORSE" \
      org.opencontainers.image.version="2.1.1"

ENTRYPOINT ["python3.8", "/app/server.py"]
