from __future__ import annotations

import json
from itertools import islice
from typing import Any, Iterable

from morse.config import ScorerConfig
from morse.construct import VaccineConstruct
from morse.router import RoutedSequence, ScorerInputRouter
from morse.scorers import make_raw_scorer


class ObjectiveScorer:
    """Own one configured Core scorer and cache its projected requests."""

    def __init__(self, config: ScorerConfig) -> None:
        self.config = config
        self._scorer = make_raw_scorer(config)
        self._cache: dict[str, dict[str, Any]] = {}
        self._router = ScorerInputRouter()

    def score(
        self,
        constructs: list[VaccineConstruct],
        *,
        objective: str,
        batch_size: int,
    ) -> list[dict[str, Any]]:
        routed = [
            self._router.route_with_context(construct, self.config.scope)
            for construct in constructs
        ]
        keys = [
            json.dumps(
                {"sequence": item.sequence, "context": item.context()},
                sort_keys=True,
                separators=(",", ":"),
            )
            for item in routed
        ]
        misses: dict[str, RoutedSequence] = {}
        for key, item in zip(keys, routed, strict=True):
            if key not in self._cache and key not in misses:
                misses[key] = item

        for batch in _batched(list(misses.items()), batch_size):
            request_keys: list[str] = []
            batch_sequences: list[str] = []
            batch_contexts: list[dict[str, object]] = []
            for key, item in batch:
                request_keys.append(key)
                batch_sequences.append(item.sequence)
                batch_contexts.append(item.context())
            batch_outputs = self._scorer.score_with_context(
                batch_sequences, batch_contexts
            )
            if len(batch_outputs) != len(batch_sequences):
                raise RuntimeError(
                    f"Scorer {self.config.name!r} returned {len(batch_outputs)} "
                    f"outputs for {len(batch_sequences)} sequences"
                )
            for key, outputs in zip(request_keys, batch_outputs, strict=True):
                outputs = dict(outputs)
                if objective not in outputs:
                    raise RuntimeError(
                        f"Scorer {self.config.name!r} did not return objective "
                        f"{objective!r}"
                    )
                self._cache[key] = outputs

        outputs = [dict(self._cache[key]) for key in keys]
        if any(objective not in item for item in outputs):
            raise RuntimeError(
                f"Scorer {self.config.name!r} did not return objective "
                f"{objective!r}"
            )
        return outputs

    def close(self) -> None:
        close = getattr(self._scorer, "close", None)
        if callable(close):
            close()


class ScorerSessions:
    def __init__(self) -> None:
        self._sessions: dict[str, ObjectiveScorer] = {}

    def get(self, config: ScorerConfig) -> ObjectiveScorer:
        session = self._sessions.get(config.name)
        if session is None:
            session = ObjectiveScorer(config)
            self._sessions[config.name] = session
        return session

    def close(self) -> None:
        for session in self._sessions.values():
            session.close()
        self._sessions.clear()


def _batched(
    values: list[tuple[str, RoutedSequence]],
    size: int,
) -> Iterable[list[tuple[str, RoutedSequence]]]:
    iterator = iter(values)
    while batch := list(islice(iterator, size)):
        yield batch
