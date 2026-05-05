from __future__ import annotations

import os
import time
from contextlib import contextmanager
from typing import Iterator


_TIMING_ENABLED = os.environ.get("RXGAMING_TIMING", "").strip().lower() in {
    "1",
    "true",
    "yes",
    "on",
    "debug",
}


def timing_enabled() -> bool:
    return _TIMING_ENABLED


@contextmanager
def timed_block(name: str, **context: object) -> Iterator[None]:
    if not _TIMING_ENABLED:
        yield
        return

    start = time.perf_counter()
    try:
        yield
    finally:
        elapsed_ms = (time.perf_counter() - start) * 1000.0
        ordered = ", ".join(f"{key}={value}" for key, value in sorted(context.items()))
        suffix = f" ({ordered})" if ordered else ""
        print(f"[rxgaming timing] {name}: {elapsed_ms:.2f} ms{suffix}", flush=True)
