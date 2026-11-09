from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

RUN_LOG_FILENAME = "run_log.jsonl.gz"
RUN_LOG_EVENT_TYPES = (
    "evaluation",
    "population_member",
    "scorer_timing",
    "initial_population",
)


def write_run_log(
    path: Path,
    evaluation_rows: list[dict[str, Any]],
    population_rows: list[dict[str, Any]],
    scorer_timing_rows: list[dict[str, Any]],
    initial_population_rows: list[dict[str, Any]] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=6) as handle:
        for event_type, rows in (
            ("evaluation", evaluation_rows),
            ("population_member", population_rows),
            ("scorer_timing", scorer_timing_rows),
            ("initial_population", initial_population_rows or []),
        ):
            for row in rows:
                handle.write(
                    json.dumps(
                        {"type": event_type, **row},
                        default=_json_default,
                        sort_keys=True,
                    )
                )
                handle.write("\n")


def read_run_log_tables(path: Path) -> dict[str, list[dict[str, Any]]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing run log: {path}")
    tables = {event_type: [] for event_type in RUN_LOG_EVENT_TYPES}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            event = json.loads(line)
            event_type = event.pop("type", None)
            if event_type not in tables:
                raise ValueError(f"Unsupported event type in {path}:{line_number}: {event_type!r}")
            tables[event_type].append(event)
    return tables


def _json_default(value: Any) -> Any:
    item = getattr(value, "item", None)
    if callable(item):
        return item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")
