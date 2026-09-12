"""Report generator: JSON / markdown / HTML rendering (doc/03 2.6, D20).

``render_*`` functions turn a ``pipeline.RunResult`` into the three report
artifacts the clinical layer must produce: a structured JSON contract, a
human-readable markdown summary and a self-contained HTML page with inline SVG
charts (no external assets, deep-purple dark theme reused by the playground).
"""

from __future__ import annotations

import html
import json
from typing import Any

import numpy as np

from drugos.organ.base import NDArray
from drugos.pipeline import RunResult
from drugos.version import __version__

_GREEN = "#86e087"
_AMBER = "#f5c268"
_RED = "#f87171"
_PURPLE_BG = "#171022"
_PURPLE_PANEL = "#241a36"
_PURPLE_PANEL_2 = "#2e2145"
_PURPLE_ACCENT = "#a78bfa"
_PURPLE_ACCENT_2 = "#7c3aed"
_TEXT = "#ece6f6"
_MUTED = "#a398c4"


def _severity_color(grade: int) -> str:
    if grade >= 3:
        return _RED
    if grade == 2:
        return _AMBER
    if grade == 1:
        return _PURPLE_ACCENT
    return _GREEN


def _risk_color(risk: float) -> str:
    if risk >= 0.5:
        return _RED
    if risk >= 0.3:
        return _AMBER
    return _GREEN


def render_json(result: RunResult, indent: int = 2) -> str:
    """Pretty JSON serialization of the full report contract."""
    return json.dumps(result.to_contract(), indent=indent)


# ---------------------------------------------------------------------------
# SVG charts
# ---------------------------------------------------------------------------
def svg_line_chart(
    t: NDArray,
    curves: dict[str, NDArray],
    *,
    width: int = 520,
    height: int = 170,
    y_label: str = "",
    pad: float = 0.08,
) -> str:
    """Inline SVG polyline chart of one or more series over ``t`` (hours)."""
    n = len(t)
    if n < 2 or any(len(v) != n for v in curves.values()) or not curves:
        raise ValueError("all curves must be equal length to t with >= 2 points")
    values = np.concatenate([np.asarray(v, dtype=float) for v in curves.values()])
    y_lo = float(np.min(values))
    y_hi = float(np.max(values))
    if y_hi - y_lo < 1e-12:
        y_hi = y_lo + 1.0
        y_lo -= 1.0
    span = y_hi - y_lo
    y_lo -= pad * span
    y_hi += pad * span
    xs = np.asarray(t, dtype=float)
    x_lo = float(xs[0])
    x_hi = float(xs[-1])
    if x_hi - x_lo < 1e-12:
        x_hi = x_lo + 1.0

    def point(x: float, y: float) -> str:
        px = 8 + (x - x_lo) / (x_hi - x_lo) * (width - 16)
        py = height - 18 - (y - y_lo) / (y_hi - y_lo) * (height - 40)
        return f"{px:.1f},{py:.1f}"

    def polyline(vals: NDArray) -> str:
        return " ".join(point(float(x), float(v)) for x, v in zip(xs, vals, strict=True))

    colors = ["#a78bfa", "#f5c268", "#7c3aed", "#86e087", "#c4b5fd", "#f87171"]
    parts = [
        f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
        f'role="img" aria-label="{html.escape(y_label or "chart")}">',
        "<style>.do-grid{stroke:#3a2c55;stroke-width:1}.do-axis{stroke:#a398c4}"
        ".do-lab{fill:#a398c4;font:10px ui-monospace,monospace}</style>",
    ]
    parts.append(f'<path d="M8 {height - 18} H{width - 8}" class="do-axis"/>')
    parts.append(f'<path d="M8 {height - 18} V{height - 56}" class="do-axis"/>')
    for i in range(4):
        gy = height - 18 - (i + 1) * (height - 40) / 4
        parts.append(f'<path d="M8 {gy:.1f} H{width - 8}" class="do-grid"/>')
        parts.append(
            f'<text x="4" y="{gy + 3:.1f}" class="do-lab" text-anchor="end">'
            f"{y_lo + (i + 1) * (y_hi - y_lo) / 4:.2g}</text>"
        )
    for i, (name, vals) in enumerate(curves.items()):
        color = colors[i % len(colors)]
        parts.append(
            f'<polyline points="{polyline(np.asarray(vals, dtype=float))}" '
            f'fill="none" stroke="{color}" stroke-width="1.6"/>'
        )
        last_x = float(xs[-1])
        last_y = float(np.asarray(vals, dtype=float)[-1])
        px, py = point(last_x, last_y).split(",")
        parts.append(
            f'<text x="{px}" y="{float(py) - 5:.1f}" fill="{color}" '
            'class="do-lab">' + html.escape(name) + "</text>"
        )
    parts.append(
        f'<text x="{width / 2:.0f}" y="{height - 4}" class="do-lab" '
        'text-anchor="middle">time (h)</text>'
    )
    parts.append("</svg>")
    return "".join(parts)


def _svg_risk_bar(risk: float, ci_lo: float, ci_hi: float) -> str:
    """Horizontal risk bar with a 95% credible-interval whisker."""
    x0, w = 40.0, 200.0
    cx = x0 + w * min(max(risk, 0.0), 1.0)
    lx = x0 + w * min(max(ci_lo, 0.0), 1.0)
    rx = x0 + w * min(max(ci_hi, 0.0), 1.0)
    return (
        f'<svg viewBox="0 0 {int(x0 + w + 4)} 14" xmlns="http://www.w3.org/2000/svg">'
        f'<rect x="{x0}" y="4" width="{w}" height="6" rx="3" fill="#3a2c55"/>'
        f'<line x1="{lx:.1f}" y1="7" x2="{rx:.1f}" y2="7" stroke="{_MUTED}" '
        f'stroke-width="3"/>'
        f'<circle cx="{cx:.1f}" cy="7" r="4" fill="{_risk_color(risk)}"/>'
        f'<text x="{x0 + w + 8}" y="10" fill="{_MUTED}" font-size="10">'
        f"{risk:.2f}</text></svg>"
    )


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------
def render_markdown(result: RunResult) -> str:
    """Human-readable markdown summary of a pipeline run."""
    liver = result.organ.liver
    cardiac = result.organ.cardiac
    kidney = result.organ.kidney
    bio = result.organ.biomarkers
    trust = result.to_contract()["trust"]
    lines = [
        f"# DrugOS pipeline report — {result.name}",
        "",
        f"- Dose {_plain(result.pk.dose_mg)} mg ({result.pk.route.value}) · "
        f"Cmax {_plain(result.metrics.cmax_mg_l)} mg/L · "
        f"AUC0-t {_plain(result.metrics.auc_last_mgh_l)} mg·h/L · "
        f"tmax {_plain(result.metrics.tmax_h)} h",
        "",
        "## Prediction reliability",
        "",
        f"- Regime: **{trust['reliability']['regime']}** — "
        f"{trust['reliability']['label']} (reliability "
        f"{trust['reliability']['reliability']}).",
        f"- Basis: {trust['reliability']['basis']}",
        f"- Recommended parameter-ensemble CV: {_plain(trust['reliability']['band_cv'])}",
        f"- {trust['reliability']['disclaimer']}",
        "",
        "## Target engagement",
        "",
        "| Target | Peak occupancy | Time-at-target (h) |",
        "|---|---|---|",
    ]
    for site, eng in result.panel.ranked_by_time_at_target():
        lines.append(f"| {site} | {_plain(eng.peak_occupancy)} | {_plain(eng.time_at_target_h)} |")
    if not result.panel.results:
        lines.append("| — | — | — |")

    lines += ["", "## Pathway signaling", ""]
    if result.pathway is None:
        lines.append("- Pathway simulation disabled for this run.")
    else:
        readout = result.pathway.model.readout
        if readout is not None:
            peak_fold = float(np.max(result.pathway.readout_fold_change()))
            lines.append(
                f"- **{readout}** peak fold-change {_plain(peak_fold)} vs drug-free baseline."
            )
        else:
            lines.append("- Pathway response simulated (no declared readout node).")

    lines += [
        "",
        "## Organ trajectories",
        "",
        f"- **Liver (DILI):** grade {liver.dili_grade}, ALT "
        f"{_plain(liver.peak_alt_uln)} xULN, bilirubin {_plain(liver.peak_bilirubin_uln)} "
        f"xULN, Hy's Law {'met' if liver.hy_law else 'not met'}",
        f"- **Cardiac:** ΔQTc {_plain(float(max(cardiac.delta_qtc_ms)))} ms, peak QTc "
        f"{_plain(float(max(cardiac.qtc_ms)))} ms ({cardiac.tdpr_band} band), "
        f"MAP {_plain(cardiac.map_mmhg)} mmHg",
        f"- **Kidney:** AKI grade {kidney.aki_grade}, peak Scr ratio "
        f"{_plain(kidney.peak_scr_ratio)}, min GFR {_plain(kidney.min_gfr_ml_min)} mL/min",
        "",
        "## Clinical biomarkers",
        "",
        "| Biomarker | Value | Unit | Grade | Severity | Onset (h) | Duration (h) |",
        "|---|---|---|---|---|---|---|",
    ]
    for _group, rows in bio.items():
        for row in rows:
            lines.append(
                f"| {row.name} | {_plain(row.value)} | {row.unit} | {row.grade} | "
                f"{row.severity} | {_opt(row.onset_h)} | {_opt(row.duration_h)} |"
            )
    lines += ["", "## Composite toxicity", ""]
    lines.append("| Endpoint | Risk | 95% CI | Grade | Driver | Reason |")
    lines.append("|---|---|---|---|---|---|")
    for r in result.toxicity.risks:
        driver_note = next(
            (e.note for e in r.evidence if e.kind.value == r.driver.value), "class prior only"
        )
        lines.append(
            f"| {r.label} | {_plain(r.risk)} | {_plain(r.ci_lo)}–{_plain(r.ci_hi)} | "
            f"{r.grade} | {r.driver.value} | {driver_note} |"
        )
    lines += [
        "",
        f"**Verdict:** {result.verdict}",
        "",
        "> Research-grade model output; not for clinical decision-making (doc/08).",
    ]
    agreement = trust.get("empirical_agreement")
    if agreement is not None:
        agreement_rows_md = "\n".join(
            f"| {row['label']} | {_plain(row['observed'])} | {_plain(row['predicted'])} | "
            f"{(str(row['fold_error']) if row['fold_error'] is not None else 'n/a')} | "
            f"{'within 2x' if row['within_2x'] else 'outside 2x'} |"
            for row in agreement["observations"]
        )
        lines += [
            "",
            "## Empirical agreement (observed vs predicted)",
            "",
            "| Endpoint | Observed | Predicted | Fold error | Within 2x |",
            "|---|---|---|---|---|",
            agreement_rows_md,
            "",
            f"_Policy: {agreement['policy']}_",
        ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------
def _engagement_panel_html(result: RunResult) -> str:
    """Target-engagement table panel (doc/05 5.3, occupancy profile)."""
    rows = "".join(
        f"<tr><td>{html.escape(site)}</td><td>{_plain(eng.peak_occupancy)}</td>"
        f"<td>{_plain(eng.time_at_target_h)}</td></tr>"
        for site, eng in result.panel.ranked_by_time_at_target()
    )
    if not rows:
        rows = "<tr><td colspan='3'>no engaged targets</td></tr>"
    return (
        '<section class="panel"><div class="h2">Target engagement</div>'
        "<table><tr><th>Target</th><th>Peak occupancy</th>"
        "<th>Time-at-target (h)</th></tr>" + rows + "</table></section>"
    )


def _pathway_panel_html(result: RunResult) -> str:
    """Pathway-signaling panel (doc/05 5.3, pathway activity profile)."""
    if result.pathway is None:
        body = '<div class="mut">Pathway simulation disabled for this run.</div>'
    else:
        readout = result.pathway.model.readout
        if readout is not None:
            peak_fold = float(np.max(result.pathway.readout_fold_change()))
            body = (
                "<table><tr><th>Readout node</th><th>Peak fold change vs baseline</th></tr>"
                f"<tr><td>{html.escape(readout)}</td><td>{_plain(peak_fold)}</td></tr></table>"
            )
        else:
            body = '<div class="mut">Pathway response simulated (no declared readout node).</div>'
    return '<section class="panel"><div class="h2">Pathway signaling</div>' + body + "</section>"


def render_html(result: RunResult) -> str:
    """Self-contained HTML report with inline SVG charts (dark purple theme)."""
    liver = result.organ.liver
    cardiac = result.organ.cardiac
    kidney = result.organ.kidney
    contract = result.to_contract()

    bio_rows = []
    for group, rows in result.organ.biomarkers.items():
        bio_rows.append(f'<h4 class="group">{html.escape(group)}</h4>')
        bio_rows.append(
            "<table><tr><th>Biomarker</th><th>Value</th><th>Unit</th>"
            "<th>Grade</th><th>Severity</th><th>Onset (h)</th><th>Duration (h)</th></tr>"
        )
        for row in rows:
            color = _severity_color(row.grade)
            bio_rows.append(
                f"<tr><td>{html.escape(row.name)}</td><td>{_plain(row.value)}</td>"
                f"<td>{row.unit}</td>"
                f'<td><span class="grade" style="background:{color}">{row.grade}</span></td>'
                f"<td>{row.severity}</td><td>{_opt(row.onset_h)}</td>"
                f"<td>{_opt(row.duration_h)}</td></tr>"
            )
        bio_rows.append("</table>")

    tox_rows = []
    for r in result.toxicity.risks:
        driver_note = next(
            (e.note for e in r.evidence if e.kind.value == r.driver.value), "class prior only"
        )
        tox_rows.append(
            f"<tr><td>{html.escape(r.label)}</td>"
            f"<td>{_svg_risk_bar(r.risk, r.ci_lo, r.ci_hi)}</td>"
            f"<td>{_plain(r.ci_lo)}–{_plain(r.ci_hi)}</td>"
            f'<td><span class="grade" '
            f'style="background:{_severity_color(r.grade)}">{r.grade}</span></td>'
            f"<td>{r.driver.value}</td><td>{html.escape(driver_note)}</td></tr>"
        )

    rel = contract["trust"]["reliability"]
    agreement = contract["trust"].get("empirical_agreement")
    emp_rows = ""
    if agreement is not None:
        emp_rows = (
            "<table><tr><th>Endpoint</th><th>Observed</th><th>Predicted</th>"
            "<th>Fold error</th><th>Within 2x</th></tr>"
        )
        for row in agreement["observations"]:
            emp_rows += (
                "<tr><td>" + html.escape(row["label"]) + "</td>"
                f"<td>{_plain(row['observed'])}</td><td>{_plain(row['predicted'])}</td>"
                f"<td>{(str(row['fold_error']) if row['fold_error'] is not None else 'n/a')}</td>"
                f"<td>{'within 2x' if row['within_2x'] else 'outside 2x'}</td></tr>"
            )
        emp_rows += "</table>"

    css = (
        f":root{{--bg:{_PURPLE_BG};--panel:{_PURPLE_PANEL};--panel2:{_PURPLE_PANEL_2};"
        f"--acc:{_PURPLE_ACCENT};--acc2:{_PURPLE_ACCENT_2};--text:{_TEXT};--mut:{_MUTED};}}"
        "body{background:linear-gradient(160deg,var(--bg),#100a18 70%);color:var(--text);"
        "font-family:ui-monospace,'SF Mono',Menlo,Consolas,monospace;margin:0;padding:28px}"
        ".wrap{max-width:960px;margin:0 auto}.panel{background:var(--panel);"
        "border:1px solid #3a2c55;border-radius:14px;padding:20px;margin:18px 0}"
        "h1{font-size:22px}.h2{color:var(--acc);font-size:15px;text-transform:uppercase;"
        "letter-spacing:.08em;margin:0 0 10px}table{border-collapse:collapse;width:100%}"
        "th,td{padding:7px 10px;text-align:left;border-bottom:1px solid #2e2145;font-size:13px}"
        "th{color:var(--mut);font-weight:600}.grade{display:inline-block;color:#171022;"
        "font-weight:700;border-radius:10px;min-width:22px;text-align:center;padding:2px 6px}"
        ".group{color:var(--acc2);margin:14px 0 4px;font-size:13px}.verdict{font-size:15px;"
        "font-weight:700}.mut{color:var(--mut);font-size:12px}svg{max-width:100%;height:auto}"
    )

    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        f"<title>DrugOS report — {html.escape(result.name)}</title>"
        f'<style>{css}</style></head><body><div class="wrap">'
        "<h1>DrugOS pipeline report</h1>"
        f'<div class="verdict" style="color:{_risk_color(result.toxicity.overall_risk())}">'
        f"{html.escape(result.verdict)}</div>"
        '<div class="mut">compound <b>'
        + html.escape(result.name)
        + "</b> · dose "
        + f"{_plain(result.pk.dose_mg)} mg ({result.pk.route.value}) · model {__version__}</div>"
        '<section class="panel"><div class="h2">Pharmacokinetics</div>'
        "<table><tr><th>Cmax (mg/L)</th><th>AUC0-t (mg·h/L)</th><th>AUC0-inf (mg·h/L)</th>"
        "<th>tmax (h)</th></tr>"
        f"<tr><td>{_plain(result.metrics.cmax_mg_l)}</td>"
        f"<td>{_plain(result.metrics.auc_last_mgh_l)}</td>"
        f"<td>{_plain(result.metrics.auc_inf_mgh_l)}</td>"
        f"<td>{_plain(result.metrics.tmax_h)}</td></tr></table>"
        + svg_line_chart(
            result.pk.t,
            {"plasma total (mg/L)": result.pk.plasma_total},
            y_label="plasma concentration",
        )
        + "</section>"
        + _engagement_panel_html(result)
        + _pathway_panel_html(result)
        + '<section class="panel"><div class="h2">Organ trajectories</div>'
        + svg_line_chart(
            liver.t_h,
            {"ALT (U/L)": liver.alt_u_l, "bilirubin (mg/dL)": liver.bilirubin_mg_dl},
            y_label="liver biomarkers",
        )
        + svg_line_chart(
            cardiac.t_h,
            {"QTc (ms)": cardiac.qtc_ms},
            y_label="QTc",
        )
        + svg_line_chart(
            kidney.t_h,
            {"GFR (mL/min)": kidney.gfr_ml_min},
            y_label="GFR",
        )
        + "</section>"
        '<section class="panel"><div class="h2">Clinical biomarkers</div>'
        + "".join(bio_rows)
        + "</section>"
        '<section class="panel"><div class="h2">Composite toxicity</div>'
        "<table><tr><th>Endpoint</th><th>Risk</th><th>95% CI</th><th>Grade</th>"
        "<th>Driver</th><th>Reason</th></tr>" + "".join(tox_rows) + "</table></section>"
        '<section class="panel"><div class="h2">Profile & exposure</div>'
        f"<table><tr><th>Sex</th><th>Route</th><th>Unbound Cmax (nM)</th><th>Liver free "
        f"Cmax (nM)</th><th>Kidney free Cmax (nM)</th><th>hERG IC50 (nM)</th>"
        f"<th>BSEP/mito IC50 (nM)</th></tr><tr>"
        f"<td>{result.spec.profile.sex.value}</td>"
        f"<td>{result.pk.route.value}</td>"
        f"<td>{_plain(contract['clinical']['exposure']['plasma_cmax_unbound_nm'])}</td>"
        f"<td>{_plain(contract['clinical']['exposure']['liver_cmax_free_nm'])}</td>"
        f"<td>{_plain(contract['clinical']['exposure']['kidney_cmax_free_nm'])}</td>"
        f"<td>{_plain(contract['clinical']['exposure']['qt_ic50_nm'])}</td>"
        f"<td>{_plain(contract['clinical']['exposure']['dili_ic50_nm'])}</td></tr></table>"
        "</section>"
        '<section class="panel"><div class="h2">Prediction reliability</div>'
        f'<div class="mut">regime <b>{html.escape(rel["regime"])}</b> — '
        f"{html.escape(rel['label'])} (reliability {html.escape(rel['reliability'])})</div>"
        f'<div class="mut" style="margin-top:6px">{html.escape(rel["basis"])}</div>'
        f'<div class="mut" style="margin-top:6px">recommended parameter-ensemble CV '
        f"{_plain(rel['band_cv'])}</div>"
        f'<div class="mut" style="margin-top:6px">{html.escape(rel["disclaimer"])}</div>'
        + (emp_rows or "")
        + "</section>"
        '<div class="mut">Research-grade model output; not for clinical decision-making '
        " (doc/08). Runs are not reproduced trials.</div>"
        "</div></body></html>"
    )


def render_all(result: RunResult) -> dict[str, str]:
    """Render the three report artifacts as strings."""
    return {
        "json": render_json(result),
        "markdown": render_markdown(result),
        "html": render_html(result),
    }


def write_report(result: RunResult, directory: str = ".") -> dict[str, str]:
    """Write report.json / report.md / report.html under ``directory``."""
    rendered = render_all(result)
    paths = {
        "json": f"{directory}/report.json",
        "markdown": f"{directory}/report.md",
        "html": f"{directory}/report.html",
    }
    for key, path in paths.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(rendered[key])
    return paths


# ---------------------------------------------------------------------------
# formatting helpers
# ---------------------------------------------------------------------------
def _plain(value: Any) -> str:
    if value is None:
        return "n/a"
    return f"{float(value):.4g}" if isinstance(value, (int, float)) else str(value)


def _opt(value: float | None) -> str:
    return _plain(value) if value is not None else "—"


__all__ = [
    "render_all",
    "render_html",
    "render_json",
    "render_markdown",
    "svg_line_chart",
    "write_report",
]
