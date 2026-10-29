"""
The HTML page (Jinja2, inline CSS, no scripts, no external resources: it works
offline from a dataset folder) and the summary tables.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from processor import analysis as A

TEMPLATES = Path(__file__).parent / "templates"


def environment() -> Environment:
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=select_autoescape(["html", "j2"]),
                      trim_blocks=True, lstrip_blocks=True)
    env.filters["pct"] = A.pct
    env.filters["pval"] = A.pval
    env.filters["num"] = A.num
    env.filters["count"] = A.count
    env.filters["pctv"] = lambda v, d=1: "–" if v is None or pd.isna(v) else f"{float(v):.{d}f}%"
    return env


def write_index(report_dir: Path, context: dict[str, Any]) -> Path:
    page = report_dir / "index.html"
    page.write_text(environment().get_template("index.html.j2").render(**context), encoding="utf-8")
    return page


def qc_table(qc: dict[str, Any]) -> pd.DataFrame:
    cols = ["sample_id", "batch", "condition", "n_events", "qc_in", "qc_out", "qc_pct", "qc_flagged", "qc_status",
            "db_in", "db_n", "db_pct", "db_out", "db_status", "kept_pct"]
    names = {"qc_in": "peacoqc_events_in", "qc_out": "peacoqc_events_out", "qc_pct": "peacoqc_pct_removed",
             "qc_flagged": "peacoqc_flagged", "qc_status": "peacoqc_status", "db_in": "cleanet_events_in",
             "db_n": "cleanet_doublets", "db_pct": "cleanet_pct_doublets", "db_out": "cleanet_events_out",
             "db_status": "cleanet_status", "kept_pct": "pct_events_kept"}
    df = pd.DataFrame(qc["rows"])
    df = df[[c for c in cols if c in df.columns]]
    return df.rename(columns=names)


def metrics_json(ctx: dict[str, Any], own: dict[str, Any]) -> dict[str, Any]:
    ov, qc, pop, da, ds, batch = (ctx.get(k) for k in ("ov", "qc", "pop", "da", "ds", "batch"))
    out: dict[str, Any] = {"app": own["app"], "appVersion": own["version"], "generatedAt": ctx["generated_at"],
                           "title": ctx["title"], "workflows": ctx["workflows"], "sectionErrors": ctx["errors"]}
    if ov:
        out["study"] = {k: ov.get(k) for k in ("n_samples", "batches", "conditions", "events", "n_markers",
                                                "n_clustering", "n_state", "dropped", "cytometer", "instruments")}
        out["study"]["observations"] = ov["observations"]
        out["study"]["caveats"] = ov["caveats"]
    if qc:
        out["qc"] = {k: qc.get(k) for k in ("peacoqc", "doublets", "observations", "caveats")}
    if pop:
        out["populations"] = {"metaclusters": pop["metaclusters"], "total_cells": pop["total_cells"],
                              "composition_means": pop.get("composition_means"),
                              "observations": pop["observations"], "caveats": pop["caveats"]}
    if da:
        out["differential_abundance"] = {k: da.get(k) for k in ("fdr", "method", "block", "comparisons",
                                                                 "observations", "caveats", "failed")}
        for c in (out["differential_abundance"].get("comparisons") or []):
            c.pop("figure", None)
    if ds:
        out["differential_state"] = {k: ds.get(k) for k in ("fdr", "method", "comparisons", "observations",
                                                             "caveats", "failed")}
    if batch:
        out["batch_effects"] = {k: batch.get(k) for k in ("labels", "runs", "scores", "reductions", "best",
                                                           "normalization", "observations", "caveats")}
    out["methods"] = [{"app": g["app"], "versions": [r["version"] for r in g["runs"]],
                       "citations": [c["doi"] for c in g["citations"]]} for g in ctx["methods"]["apps"]]
    out["inputData"] = ctx["credit"]
    return out
