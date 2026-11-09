# Build context: paper/baselines/vaxpress; use --platform linux/amd64.
FROM python:3.12-slim@sha256:e695ab189fa2d082b1c61c69a6b888a5adeca61ceeccad107b1ed2a125721a9e
LABEL org.opencontainers.image.title="MORSE native VaxPress shared-fitness baseline" \
      org.opencontainers.image.source="https://github.com/ChangLabSNU/VaxPress" \
      org.opencontainers.image.revision="dfc36638978e4c8f0925045018fd65678d3dcb3e"
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
COPY source-lock.json /baseline/source-lock.json
COPY requirements.lock.txt /baseline/requirements.lock.txt
COPY fetch_source.py /baseline/fetch_source.py
COPY verify_runtime.py /baseline/verify_runtime.py
RUN python /baseline/fetch_source.py --output /opt/VaxPress \
    && python -m pip install --only-binary=:all: --require-hashes -r /baseline/requirements.lock.txt \
    && python -m pip install --no-deps --no-build-isolation /opt/VaxPress \
    && python -m pip check \
    && python /baseline/verify_runtime.py
WORKDIR /work
ENTRYPOINT ["python3"]
CMD ["/baseline/verify_runtime.py"]
