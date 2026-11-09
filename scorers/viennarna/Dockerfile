FROM python:3.12-slim@sha256:e695ab189fa2d082b1c61c69a6b888a5adeca61ceeccad107b1ed2a125721a9e

LABEL org.opencontainers.image.source="https://github.com/AIRI-Institute/MORSE" \
      org.opencontainers.image.version="2.2.0"

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

RUN pip install --no-cache-dir "ViennaRNA==2.7.2"

WORKDIR /app
COPY scorer.py /app/scorer.py
COPY LICENSE ViennaRNA-COPYING.txt /app/licenses/

ENTRYPOINT ["python", "/app/scorer.py"]
