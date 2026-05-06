import pytest
from errors import (
    PipelineError,
    FLYE_FAILED, MLST_FAILED, TOOL_NOT_FOUND, STAGE_TIMEOUT,
    ERROR_LABELS,
)


def test_pipeline_error_is_exception():
    with pytest.raises(PipelineError):
        raise PipelineError(
            code=FLYE_FAILED, message="Flye assembly failed",
            detail="exit code 1", stage="assembly",
        )


def test_pipeline_error_str_contains_code_and_message():
    e = PipelineError(
        code=MLST_FAILED, message="MLST typing failed",
        detail="mlst exited 1", stage="mlst",
    )
    assert "GP-401" in str(e)
    assert "MLST typing failed" in str(e)


def test_error_labels_cover_all_constants():
    codes = [
        FLYE_FAILED, MLST_FAILED, TOOL_NOT_FOUND, STAGE_TIMEOUT,
    ]
    for code in codes:
        assert code in ERROR_LABELS, f"{code} missing from ERROR_LABELS"


def test_recoverable_defaults_false():
    e = PipelineError(code=FLYE_FAILED, message="x", detail="", stage="assembly")
    assert e.recoverable is False


def test_recoverable_can_be_true():
    e = PipelineError(
        code=STAGE_TIMEOUT, message="timed out", detail="",
        stage="assembly", recoverable=True,
    )
    assert e.recoverable is True


@pytest.fixture
def tmp_db(monkeypatch, tmp_path):
    import database as db
    monkeypatch.setattr(db, "DB_PATH", tmp_path / "test.db")
    db.init_db()
    return db


def test_error_code_column_persists(tmp_db):
    db = tmp_db
    job_id = "a" * 32
    db.create_job(job_id, "test.fastq.gz", "2099-01-01T00:00:00+00:00", "iphash")
    db.update_job_status(
        job_id, "failed",
        error="[GP-301] Flye assembly failed: exit code 1",
        error_code="GP-301",
        error_detail="exit code 1",
    )
    job = db.get_job(job_id)
    assert job["error_code"] == "GP-301"
    assert job["error_detail"] == "exit code 1"
    assert job["status"] == "failed"


def test_error_code_defaults_to_none(tmp_db):
    db = tmp_db
    job_id = "b" * 32
    db.create_job(job_id, "test.fastq.gz", "2099-01-01T00:00:00+00:00", "iphash")
    job = db.get_job(job_id)
    assert job.get("error_code") is None
