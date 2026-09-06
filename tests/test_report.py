"""Stage-5 report generator tests (drugos.report)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from drugos.pipeline import RunSpec, run_pipeline
from drugos.report.render import (
    _plain,
    _risk_color,
    _severity_color,
    render_all,
    render_html,
    render_json,
    render_markdown,
    svg_line_chart,
    write_report,
)


def test_severity_and_risk_colors() -> None:
    assert _severity_color(4) == "#f87171"
    assert _severity_color(2) == "#f5c268"
    assert _severity_color(1) == "#a78bfa"
    assert _severity_color(0) == "#86e087"
    assert _risk_color(0.9) == "#f87171"
    assert _risk_color(0.3) == "#f5c268"
    assert _risk_color(0.1) == "#86e087"


def test_plain() -> None:
    assert _plain(None) == "n/a"
    assert _plain(3.14159) == "3.142"
    assert _plain(2) == "2"
    assert _plain("abc") == "abc"


def test_svg_line_chart() -> None:
    t = np.linspace(0.0, 6.0, 61)
    svg = svg_line_chart(t, {"plasma": t**0.5}, y_label="concentration")
    assert svg.startswith("<svg")
    assert 'aria-label="concentration"' in svg
    assert "plasma" in svg

    flat = svg_line_chart(np.linspace(0.0, 6.0, 61), {"a": np.ones(61)})
    assert "a" in flat

    single_t = np.array([2.0, 2.0])
    assert svg_line_chart(single_t, {"x": np.array([1.0, 1.0])})

    with pytest.raises(ValueError):
        svg_line_chart(np.array([1.0, 2.0, 3.0]), {"x": np.array([1.0, 1.0, 1.0, 1.0])})
    with pytest.raises(ValueError):
        svg_line_chart(np.array([1.0]), {"x": np.array([1.0])})
    with pytest.raises(ValueError):
        svg_line_chart(np.array([1.0, 2.0]), {})


def test_svg_line_chart_html_escape() -> None:
    svg = svg_line_chart(np.array([0.0, 1.0]), {"<img onerror=x>": np.array([1.0, 2.0])})
    assert "<img" not in svg


def test_render_json_roundtrip(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    text = render_json(result)
    data = json.loads(text)
    assert data["clinical"]["verdict"] == result.verdict
    assert "cns" in data["organ"]
    assert data["manifest"]["name"] == "warfarin"


def test_render_markdown(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    md = render_markdown(result)
    assert md.startswith("# DrugOS pipeline report — warfarin")
    assert "Organ trajectories" in md
    assert "CNS exposure ratio" in md
    assert "Hy's Law" in md
    assert "Target engagement" in md
    assert "Peak occupancy" in md
    assert "Pathway signaling" in md
    assert "fold-change" in md
    assert "class prior only" not in md or "drug-induced" in md.lower()


def test_render_html(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    page = render_html(result)
    assert page.startswith("<!doctype html>")
    assert "DrugOS pipeline report" in page
    assert "<svg" in page
    assert "not for clinical decision-making" in page
    assert "Target engagement" in page
    assert "Pathway signaling" in page
    assert "Peak occupancy" in page


def test_render_no_pathway(fast_warfarin: RunSpec, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fast_warfarin, "include_pathway", False)
    result = run_pipeline(fast_warfarin)
    md = render_markdown(result)
    page = render_html(result)
    assert "Pathway simulation disabled for this run." in md
    assert "Pathway simulation disabled for this run." in page


def test_render_empty_panel(no_safety_panel: RunSpec) -> None:
    result = run_pipeline(no_safety_panel)
    md = render_markdown(result)
    page = render_html(result)
    assert "| — | — | — |" in md  # no engaged targets, markdown row
    assert "no engaged targets" in page  # empty engagement table, html


def test_render_pathway_without_readout(fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    assert result.pathway is not None
    result.pathway.model.readout = None
    md = render_markdown(result)
    page = render_html(result)
    assert "no declared readout node" in md
    assert "no declared readout node" in page


def test_render_all_and_write(tmp_path: Path, fast_warfarin: RunSpec) -> None:
    result = run_pipeline(fast_warfarin)
    artifacts = render_all(result)
    assert set(artifacts) == {"json", "markdown", "html"}
    paths = write_report(result, str(tmp_path))
    assert set(paths) == {"json", "markdown", "html"}
    assert (tmp_path / "report.json").exists()
    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "report.html").exists()
    assert (tmp_path / "report.md").read_text().startswith("# DrugOS")
    assert json.loads((tmp_path / "report.json").read_text())["manifest"]["dose_mg"] == 5.0
