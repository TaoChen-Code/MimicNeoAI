# nodemon_pool.py
# -*- coding: utf-8 -*-
"""
Non-daemon multiprocessing Pool.

Why:
  - The default multiprocessing Pool uses daemon worker processes.
  - Daemon processes cannot spawn child processes.
  - Some bioinformatics pipelines / wrappers (e.g., callers that spawn subprocess trees)
    may require non-daemon workers.

Usage:
  from nodemon_pool import NoDaemonPool

  with NoDaemonPool(processes=4) as pool:
      ...
"""

from __future__ import annotations

import multiprocessing
import multiprocessing.pool


class _NoDaemonProcess(multiprocessing.Process):
    """Process class that always behaves as non-daemon (allows child processes)."""

    # For Python < 3.8 compatibility pattern (property override)
    def _get_daemon(self) -> bool:
        return False

    def _set_daemon(self, value: bool) -> None:
        # ignore any attempt to set daemon=True
        return

    daemon = property(_get_daemon, _set_daemon)


# Python ≥ 3.8: override Pool's Process factory
if hasattr(multiprocessing, "get_start_method"):
    class NoDaemonPool(multiprocessing.pool.Pool):
        """A multiprocessing Pool whose worker processes are non-daemon."""

        @staticmethod
        def Process(_, *args, **kwargs):
            return _NoDaemonProcess(*args, **kwargs)

else:
    # Python < 3.8 fallback
    class NoDaemonPool(multiprocessing.pool.Pool):
        """A multiprocessing Pool whose worker processes are non-daemon."""
        Process = _NoDaemonProcess