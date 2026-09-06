"""Report layer: JSON / markdown / HTML rendering of pipeline runs (D20)."""

from drugos.report.render import (
    render_all,
    render_html,
    render_json,
    render_markdown,
    svg_line_chart,
    write_report,
)

__all__ = [
    "render_all",
    "render_html",
    "render_json",
    "render_markdown",
    "svg_line_chart",
    "write_report",
]
