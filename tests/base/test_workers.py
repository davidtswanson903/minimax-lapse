"""An experiment's independent jobs share the cores a run is given: one worker per job, the cores' threads shared
among them, the result in order and the same however many workers run it."""
import math
import os

from repokit import workers


def test_the_cores_are_shared_among_one_worker_per_job():
    assert workers.shares(6, 24) == (6, 4)
    assert workers.shares(232, 24) == (24, 1)
    assert workers.shares(1, 24) == (1, 24)
    assert workers.shares(6, 1) == (1, 1)


def test_the_jobs_come_back_in_order_with_one_worker_or_several(monkeypatch):
    monkeypatch.setenv("REPOKIT_WORKERS", "1")
    serial = workers.map(math.sqrt, [1.0, 4.0, 9.0])
    monkeypatch.setenv("REPOKIT_WORKERS", "2")
    before = os.environ.get("OPENBLAS_NUM_THREADS")
    parallel = workers.map(math.sqrt, [1.0, 4.0, 9.0])
    assert serial == parallel == [1.0, 2.0, 3.0]
    assert os.environ.get("OPENBLAS_NUM_THREADS") == before, "the threads are restored after the workers end"
