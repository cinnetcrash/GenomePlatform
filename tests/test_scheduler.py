import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

import pytest
import scheduler
from config import MAX_THREADS


@pytest.fixture(autouse=True)
def clean_scheduler():
    scheduler._active_jobs.clear()
    yield
    scheduler._active_jobs.clear()


def test_no_active_jobs_uses_max_threads():
    assert scheduler.threads_for_stage() == MAX_THREADS


def test_two_jobs_halves_threads():
    scheduler.job_start("job1")
    scheduler.job_start("job2")
    expected = max(1, MAX_THREADS // 2)
    assert scheduler.threads_for_stage() == expected


def test_floor_is_one():
    for i in range(MAX_THREADS + 10):
        scheduler.job_start(f"job{i}")
    assert scheduler.threads_for_stage() == 1


def test_job_done_removes_from_active():
    scheduler.job_start("jobA")
    scheduler.job_done("jobA")
    assert scheduler.active_count() == 0
    assert scheduler.threads_for_stage() == MAX_THREADS


def test_job_done_unknown_id_is_noop():
    scheduler.job_done("does-not-exist")
    assert scheduler.active_count() == 0


def test_active_count():
    scheduler.job_start("x")
    scheduler.job_start("y")
    assert scheduler.active_count() == 2
