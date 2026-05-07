"""Retry helpers for recoverable failures.

`retry_once` — two-attempt wrapper for per-item operations (e.g. one
protein hitting a web API). Quick fail-fast behaviour.

`retry_with_backoff` — N-attempt wrapper with linear backoff
(initial_delay × attempt) for whole-batch operations (e.g. submitting
a 500-protein batch to a DTU server). Sleeps longer between tries
because the cost of giving up is much higher.
"""

from __future__ import annotations

import logging
import time
from typing import Any, Callable

logger = logging.getLogger(__name__)


def retry_with_backoff(
    func: Callable[[], Any],
    *,
    max_attempts: int = 3,
    initial_delay: float = 30.0,
    label: str = "operation",
    retry_on: type[BaseException] | tuple[type[BaseException], ...] = Exception,
) -> Any:
    """Run *func* with linear backoff between failures.

    Sleeps ``initial_delay × attempt`` seconds between attempts (so 30 s,
    60 s, 90 s with default settings). Returns the first successful
    result; re-raises the final exception after exhausting *max_attempts*.

    Args:
        func: Zero-argument callable. Wrap your real call in a lambda
            if it has arguments.
        max_attempts: Total tries including the first. Must be >= 1.
        initial_delay: Base sleep, scaled by attempt number.
        label: Prefix for log messages — caller's identifier (e.g. the
            batch number, or "DTU SignalP").
        retry_on: Exception class or tuple. Other exceptions propagate
            immediately without retry.
    """
    if max_attempts < 1:
        raise ValueError(f"max_attempts must be >= 1, got {max_attempts}")

    for attempt in range(1, max_attempts + 1):
        try:
            return func()
        except retry_on as e:
            if attempt < max_attempts:
                wait = initial_delay * attempt
                logger.warning(
                    "%s attempt %d/%d failed: %s. Retrying in %ss...",
                    label,
                    attempt,
                    max_attempts,
                    e,
                    wait,
                )
                time.sleep(wait)
            else:
                raise


def retry_once(func: Callable, item_id: str, delay: float = 1.0, **kwargs) -> tuple[Any, str]:
    """Execute *func*, retry once on failure, return ``(result, status)``.

    Args:
        func: Callable that takes *item_id* as its first argument plus
            any additional keyword arguments.
        item_id: Identifier for logging (e.g. a protein ID).
        delay: Seconds to wait between attempts.
        **kwargs: Extra keyword arguments forwarded to *func*.

    Returns:
        ``(result, "success")`` on success.
        ``(None, "failed: <reason>")`` after two consecutive failures.
    """
    for attempt in range(2):
        try:
            result = func(item_id, **kwargs)
            return result, "success"
        except Exception as e:
            error_msg = f"{type(e).__name__}: {str(e)[:200]}"
            if attempt == 0:
                logger.warning(
                    "%s: attempt 1 failed (%s), retrying in %ss...",
                    item_id,
                    error_msg,
                    delay,
                )
                time.sleep(delay)
            else:
                logger.error("%s: attempt 2 failed (%s), skipping", item_id, error_msg)
                return None, f"failed: {error_msg}"

    # Should never reach here, but guard against it.
    return None, "failed: unknown"
