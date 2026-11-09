FROM python:3.12-slim@sha256:e695ab189fa2d082b1c61c69a6b888a5adeca61ceeccad107b1ed2a125721a9e

LABEL org.opencontainers.image.source="https://github.com/AIRI-Institute/MORSE" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.licenses="MIT"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    CSI_PROFILE=/app/profiles/homo_sapiens_consensus_wu2019_v1.json

WORKDIR /app
COPY model.py server.py LICENSE /app/
COPY profiles /app/profiles

ENTRYPOINT ["python", "/app/server.py"]
