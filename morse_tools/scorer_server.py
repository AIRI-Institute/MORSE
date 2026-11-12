"""Optional localhost HTTP adapter for one already prepared JSONL scorer."""
from __future__ import annotations

from http.server import BaseHTTPRequestHandler, HTTPServer
import json
import math
import os
import selectors
import signal
import subprocess
import time

from morse.scorers import PersistentProcessScorer
from morse.config import scorer_parameters_from_mapping

MAX_BYTES = 64 * 1024 * 1024


class _Backend(PersistentProcessScorer):
    def __init__(self, command, timeout_seconds):
        super().__init__(["/bin/bash", "-c", command], timeout_seconds=timeout_seconds)
        self.started = False

    def _ensure_process(self):
        if self.started:
            if self._process is None or self._process.poll() is not None:
                raise RuntimeError("The prepared scorer stopped; restart it explicitly")
            return self._process
        self._process = subprocess.Popen(self.command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                         stderr=None, text=True, bufsize=1, start_new_session=True)
        self.started = True
        return self._process

    def _readline_with_timeout(self, process):
        deadline = time.monotonic() + self.timeout_seconds
        output = bytearray()
        with selectors.DefaultSelector() as selector:
            selector.register(process.stdout, selectors.EVENT_READ)
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0 or not selector.select(remaining):
                    raise TimeoutError("Prepared scorer response timed out")
                chunk = os.read(process.stdout.fileno(), min(65536, MAX_BYTES + 1 - len(output)))
                if not chunk:
                    raise RuntimeError("Prepared scorer closed stdout")
                output.extend(chunk)
                if len(output) > MAX_BYTES:
                    raise RuntimeError("Prepared scorer stdout exceeds the response limit")
                if b"\n" in output:
                    line, tail = output.split(b"\n", 1)
                    if tail.strip():
                        raise RuntimeError("Prepared scorer emitted more than one JSON response")
                    try:
                        return line.decode()
                    except UnicodeDecodeError as error:
                        raise RuntimeError("Prepared scorer stdout is not UTF-8 JSON") from error

    def close(self):
        process, self._process = self._process, None
        if process is None:
            return
        if process.stdin:
            try:
                process.stdin.close()
            except OSError:
                pass
        try:
            process.wait(timeout=2)
        except subprocess.TimeoutExpired:
            pass
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(process.pid, sig)
            except ProcessLookupError:
                break
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                pass
        if process.stdout:
            process.stdout.close()


def serve(command, port, timeout_seconds=1800):
    """Explicitly start one persistent command; never construct Docker commands."""
    if not isinstance(command, str) or not command.strip() or not isinstance(port, int) or isinstance(port, bool) or not 0 < port < 65536 or not isinstance(timeout_seconds, (float, int)) or not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError("Supply a nonempty command, port1–65535 and positive timeout")
    backend = _Backend(command, timeout_seconds)

    class Handler(BaseHTTPRequestHandler):
        def reply(self, status, body):
            payload = json.dumps(body, allow_nan=False).encode()
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)

        def do_GET(self):
            if self.path != "/health":
                return self.reply(404, {"error": "Use POST /score"})
            ready = backend._process is not None and backend._process.poll() is None
            self.reply(200 if ready else 503, {"ready": ready})

        def do_POST(self):
            request_id = None
            try:
                if self.path != "/score":
                    return self.reply(404, {"error": "Use POST /score"})
                length = int(self.headers.get("Content-Length", "-1"))
                if not 0 <= length <= MAX_BYTES:
                    raise ValueError("Invalid or oversized JSON request")
                self.connection.settimeout(30)
                request = json.loads(self.rfile.read(length))
                request_id = request["request_id"]
                sequences = request["sequences"]
                if not isinstance(sequences, list) or not all(isinstance(x, str) for x in sequences):
                    raise ValueError("sequences must be a list of strings")
                parameters = scorer_parameters_from_mapping(request)
            except (ValueError, KeyError, TypeError) as error:
                return self.reply(400, {"request_id": request_id, "error": str(error)[:2000]})
            try:
                outputs = backend.score_with_context(sequences, request.get("sequence_contexts"),
                                                     parameters=parameters)
                response = {"request_id": request_id, "outputs": outputs}
                if backend.resolved_parameters is not None:
                    response["parameters"] = backend.resolved_parameters
                json.dumps(response, allow_nan=False)
                self.reply(200, response)
            except (OSError, RuntimeError, TimeoutError, ValueError) as error:
                backend.close()
                self.reply(503, {"request_id": request_id, "error": str(error)[:2000]})

    def terminate(signum, frame):
        raise KeyboardInterrupt

    previous = signal.signal(signal.SIGTERM, terminate)
    try:
        with HTTPServer(("127.0.0.1", port), Handler) as server:
            backend._ensure_process()  # Explicit utility startup; never on an HTTP request.
            try:
                server.serve_forever()
            except KeyboardInterrupt:
                pass
    finally:
        backend.close()
        signal.signal(signal.SIGTERM, previous)
    return 0
