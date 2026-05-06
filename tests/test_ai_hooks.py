import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

from unittest.mock import patch, MagicMock
import pytest


def _mock_message(text: str) -> MagicMock:
    msg = MagicMock()
    msg.content = [MagicMock(text=text)]
    return msg


def _reload_ai(monkeypatch, key: str):
    """Reload config then ai_interpreter so ANTHROPIC_API_KEY is picked up."""
    import importlib, config, ai_interpreter
    monkeypatch.setenv("ANTHROPIC_API_KEY", key)
    importlib.reload(config)
    importlib.reload(ai_interpreter)
    return ai_interpreter


def test_explain_error_returns_string(monkeypatch):
    ai_interpreter = _reload_ai(monkeypatch, "test-key")
    with patch("ai_interpreter.anthropic.Anthropic") as MockClient:
        MockClient.return_value.messages.create.return_value = _mock_message(
            "The Flye assembler exited with a non-zero code."
        )
        result = ai_interpreter.explain_error(
            error_code="GP-301", stage="assembly", detail="flye: error: exit code 1",
        )
    assert isinstance(result, str) and len(result) > 0


def test_explain_error_returns_none_without_api_key(monkeypatch):
    ai_interpreter = _reload_ai(monkeypatch, "")
    result = ai_interpreter.explain_error(
        error_code="GP-301", stage="assembly", detail="err"
    )
    assert result is None


def test_explain_error_returns_none_on_api_exception(monkeypatch):
    ai_interpreter = _reload_ai(monkeypatch, "test-key")
    with patch("ai_interpreter.anthropic.Anthropic") as MockClient:
        MockClient.return_value.messages.create.side_effect = Exception("API error")
        result = ai_interpreter.explain_error(
            error_code="GP-502", stage="mlst", detail="tool not found"
        )
    assert result is None


def test_interpret_assembly_qc_returns_string(monkeypatch):
    ai_interpreter = _reload_ai(monkeypatch, "test-key")
    with patch("ai_interpreter.anthropic.Anthropic") as MockClient:
        MockClient.return_value.messages.create.return_value = _mock_message(
            "Assembly N50 of 45 kb is below threshold."
        )
        result = ai_interpreter.interpret_assembly_qc(
            stats={"n50_bp": 45000, "total_contigs": 180,
                   "total_length_bp": 5100000, "gc_percent": 51.2},
            top_organism="Escherichia coli",
        )
    assert isinstance(result, str) and len(result) > 0


def test_interpret_assembly_qc_returns_none_without_key(monkeypatch):
    ai_interpreter = _reload_ai(monkeypatch, "")
    result = ai_interpreter.interpret_assembly_qc(
        stats={"n50_bp": 100000, "total_contigs": 50,
               "total_length_bp": 5000000, "gc_percent": 51.0},
        top_organism=None,
    )
    assert result is None


def test_summarise_comparison_returns_string(monkeypatch):
    ai_interpreter = _reload_ai(monkeypatch, "test-key")
    with patch("ai_interpreter.anthropic.Anthropic") as MockClient:
        MockClient.return_value.messages.create.return_value = _mock_message(
            "Samples A and B cluster tightly."
        )
        result = ai_interpreter.summarise_comparison(
            sample_names=["SampleA", "SampleB"],
            mlst_types=["ST131", "ST131"],
            amr_matrix={"SampleA": ["blaCTX-M-15"], "SampleB": ["blaCTX-M-15"]},
            closest_pair=("SampleA", "SampleB", 0.0003),
        )
    assert isinstance(result, str) and len(result) > 0


def test_summarise_comparison_returns_none_without_key(monkeypatch):
    ai_interpreter = _reload_ai(monkeypatch, "")
    result = ai_interpreter.summarise_comparison(
        sample_names=["A", "B"], mlst_types=["ST131", "ST10"],
        amr_matrix={}, closest_pair=None,
    )
    assert result is None


def test_interpret_assembly_qc_returns_none_on_api_exception(monkeypatch):
    ai_interpreter = _reload_ai(monkeypatch, "test-key")
    with patch("ai_interpreter.anthropic.Anthropic") as MockClient:
        MockClient.return_value.messages.create.side_effect = Exception("API error")
        result = ai_interpreter.interpret_assembly_qc(
            stats={"n50_bp": 100000, "total_contigs": 50,
                   "total_length_bp": 5000000, "gc_percent": 51.0},
            top_organism="E. coli",
        )
    assert result is None


def test_summarise_comparison_returns_none_on_api_exception(monkeypatch):
    ai_interpreter = _reload_ai(monkeypatch, "test-key")
    with patch("ai_interpreter.anthropic.Anthropic") as MockClient:
        MockClient.return_value.messages.create.side_effect = Exception("API error")
        result = ai_interpreter.summarise_comparison(
            sample_names=["A", "B"], mlst_types=["ST131", "ST10"],
            amr_matrix={"A": [], "B": []}, closest_pair=None,
        )
    assert result is None
