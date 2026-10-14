"""Workers for an experiment's independent jobs, sharing the cores a run is given.

A run on a compute server sets REPOKIT_WORKERS to its cores and the linear-algebra libraries to one thread each, so
that many single-threaded workers fill the cores, as a scan of many small models does. An experiment with fewer
independent jobs than cores would then leave most of them idle. `map` runs one worker per job, up to the workers
given, and shares the cores among them as linear-algebra threads. The workers are started fresh rather than forked,
with the thread counts set before they start, because the libraries read their thread counts when they load. The
result does not depend on how many workers run it.
"""
from __future__ import annotations

import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor

THREADS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def shares(jobs: int, workers: int) -> tuple[int, int]:
    """The worker processes and the threads each takes: one process per job up to the workers given, and the cores
    shared equally among them."""
    processes = max(1, min(workers, jobs))
    return processes, max(1, workers // processes)


def map(func, jobs: list) -> list:
    """`func` over the jobs, in order: in this process when one worker is given, else in fresh workers, each with its
    share of the cores' threads. `func` must be importable from its module, as a worker imports it."""
    workers = int(os.environ.get("REPOKIT_WORKERS", "1"))
    if workers <= 1:
        return [func(job) for job in jobs]
    processes, threads = shares(len(jobs), workers)
    saved = {name: os.environ.get(name) for name in THREADS}
    try:
        for name in THREADS:
            os.environ[name] = str(threads)
        with ProcessPoolExecutor(processes, mp_context=multiprocessing.get_context("spawn")) as pool:
            return list(pool.map(func, jobs))
    finally:
        for name, value in saved.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
