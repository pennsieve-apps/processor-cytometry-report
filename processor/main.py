"""
processor-cytometry-report: the last step of both immune cytometry workflows
(A. Immune Profiling & Differential Abundance, B. Multi-batch Harmonization & QC).
It reads whatever upstream outputs it is given and writes

  OUTPUT_DIR/report/index.html          one self-contained page (no scripts, no external resources)
  OUTPUT_DIR/report/figures/*.png       its figures (own and copied from upstream steps)
  OUTPUT_DIR/report/provenance.json
  OUTPUT_DIR/tables/composition.csv     sample x label fractions, with batch and condition
  OUTPUT_DIR/tables/da_significant.csv  significant metaclusters per comparison
  OUTPUT_DIR/tables/batch_scores.csv    per batch and assessment LABEL
  OUTPUT_DIR/tables/qc.csv              events per sample through QC and doublet removal
  OUTPUT_DIR/tables/metrics.json        the numbers behind the page

Every input is optional; a section is rendered when its inputs are present. The
run fails only when none of the inputs is found.
"""
from __future__ import annotations

import json
import logging
import os
import shutil
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("processor-cytometry-report")
# The Lambda runtime configures the root logger at WARNING before basicConfig can.
log.setLevel(logging.INFO)
for noisy in ("matplotlib", "PIL", "fontTools"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "processor-cytometry-report"
OUTPUTS = ("report", "tables")
PARAMETERS = ("REPORT_TITLE", "FDR")
DEFAULT_TITLE = "Immune Cytometry Report"
REPO_ROOT = Path(__file__).resolve().parent.parent


def _now() -> datetime:
    return datetime.now(timezone.utc)


def version_of(root: Path) -> str:
    for p in (root / "VERSION", Path("/processor/VERSION")):
        if p.is_file():
            return p.read_text().strip()
    return "unknown"


def tool_versions() -> dict[str, str]:
    import jinja2
    import matplotlib
    import numpy
    import pandas

    return {"python": sys.version.split()[0], "pandas": pandas.__version__, "numpy": numpy.__version__,
            "matplotlib": matplotlib.__version__, "jinja2": jinja2.__version__}


def parse_fdr(value: str | None) -> float | None:
    if value is None or not str(value).strip():
        return None
    f = float(value)
    if not 0 < f < 1:
        raise ValueError(f"FDR must be between 0 and 1, got {value}")
    return f


def _copy(src: Path | None, fig_dir: Path, name: str) -> str | None:
    if src is None or not Path(src).is_file():
        return None
    shutil.copyfile(src, fig_dir / name)
    return name


def _section(name: str, fn, errors: list[dict[str, str]], *args):
    """Run one section; a failing section is logged, reported on the page, and the rest goes on."""
    t0 = time.monotonic()
    try:
        result = fn(*args)
        log.info("section %s: %s (%.1fs)", name, "rendered" if result else "no inputs", time.monotonic() - t0)
        return result
    except Exception as exc:  # noqa: BLE001 - one broken input must not sink the report
        log.error("section %s FAILED: %s\n%s", name, exc, traceback.format_exc())
        errors.append({"section": name, "error": f"{type(exc).__name__}: {exc}"})
        return None


def main() -> int:
    started = _now()
    t_start = time.monotonic()
    input_dir = os.environ.get("INPUT_DIR")
    output_dir = os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

    title = os.environ.get("REPORT_TITLE") or DEFAULT_TITLE
    try:
        fdr = parse_fdr(os.environ.get("FDR"))
    except ValueError as exc:
        log.error("bad FDR: %s", exc)
        return 2
    parameters = {"REPORT_TITLE": title, "FDR": os.environ.get("FDR", "")}

    from processor import analysis as A, figures as F, inputs as I, report as R

    inp = I.gather(input_dir)
    for port in I.PORTS:
        files = inp.found.get(port, [])
        log.info("input %s: %d file(s)%s", port, len(files),
                 (": " + ", ".join(inp.rel(f) for f in files[:6]) + (" …" if len(files) > 6 else "")) if files else "")
    log.info("input assessment labels: %s; normalization methods: %s; provenance files: %d",
             ", ".join(inp.assessments) or "none", ", ".join(inp.normalization) or "none", len(inp.provenances))
    if not inp.has_any:
        listing = [inp.rel(p) for p in inp.files]
        log.error("none of the expected inputs was found (samples.csv, panel.csv, peacoqc.csv, doublets.csv, "
                  "counts.csv, umap_*.png, da.csv, ds.csv, a batch-assess summary.json, cytonorm/cycombine qc); "
                  "INPUT_DIR holds %d file(s): %s", len(listing), ", ".join(listing[:40]) + (" …" if len(listing) > 40 else ""))
        return 1

    report_dir = output_dir / "report"
    tables_dir = output_dir / "tables"
    fig_dir = report_dir / "figures"
    for d in (report_dir, tables_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    errors: list[dict[str, str]] = []
    figs: dict[str, Any] = {}
    reference_level = (inp.da_summary or {}).get("reference")

    ov = _section("overview", A.overview, errors, inp)
    qc = _section("qc", A.qc, errors, inp)
    if qc:
        figs["qc"] = _section("qc figure", F.qc_bars, errors, qc["rows"], fig_dir)
        qc["plots"] = [(sid, _copy(p, fig_dir, f"peacoqc_{sid}.png")) for sid, p in qc.get("peacoqc_plots", [])]
    pop = _section("populations", A.populations, errors, inp, reference_level)
    if pop:
        figs["heatmap"] = _section("heatmap", F.marker_heatmap, errors, inp.medians, pop["clustering_markers"] or [],
                                   pop["state_markers"] or [], fig_dir)
        if not figs["heatmap"]:
            figs["heatmap"] = _copy(pop.get("flowsom_heatmap"), fig_dir, "flowsom_heatmap.png")
        figs["composition"] = _section("composition figure", F.composition_bars, errors, pop["composition"], fig_dir)
        figs["star_tree"] = _copy(pop.get("star_tree"), fig_dir, "flowsom_star_tree.png")
        figs["umap"] = [(p.stem.replace("umap_", "").replace("_", " "), _copy(p, fig_dir, p.name))
                        for p in pop["umap_figures"]]
    da = _section("differential abundance", A.differential_abundance, errors, inp, fdr)
    if da and not da.get("failed"):
        for i, comp in enumerate(da["comparisons"]):
            comp["figure"] = _section("DA figure", F.da_points, errors, inp.counts, inp.samples, comp, fig_dir,
                                      f"da_{i + 1}.png")
        figs["volcano"] = _copy(da.get("volcano"), fig_dir, "diffcyt_volcano.png")
    ds = _section("differential state", A.differential_state, errors, inp, fdr)
    if ds and not ds.get("failed"):
        ds["figure_files"] = [_copy(p, fig_dir, f"diffcyt_ds_{p.name}") for p in ds.get("figures", [])]
    batch = _section("batch effects", A.batch_effects, errors, inp)
    if batch and batch.get("labels"):
        figs["batch_scores"] = _section("batch score figure", F.batch_scores, errors, batch, fig_dir)
        comparison = "reference" if any(A.marker_emd(inp, l, "reference") is not None for l in batch["labels"]) else "batches"
        mats = {l: A.marker_emd(inp, l, comparison) for l in batch["labels"]}
        flag_sets = {}
        for l in batch["labels"]:
            f = inp.assessments[l]["flags"]
            if f is not None and {"level", "comparison", "batch", "marker"} <= set(f.columns):
                sel = f[(f["level"] == "marker") & (f["comparison"] == comparison)]
                flag_sets[l] = {(str(b), str(m)) for b, m in zip(sel["batch"], sel["marker"])}
        batch["heatmap_comparison"] = comparison
        figs["batch_emd"] = _section("batch EMD figure", F.emd_heatmaps, errors, mats, flag_sets, fig_dir,
                                     "EMD per batch and marker, " + ("each batch's reference sample vs the pooled references"
                                                                     if comparison == "reference" else
                                                                     "each batch vs the other batches (mean over conditions)"))

    app_version = version_of(REPO_ROOT)
    tools = tool_versions()
    own = {"app": APP_NAME, "version": app_version, "tools": tools, "parameters": parameters}
    meth = _section("methods", A.methods, errors, inp, own) or {"apps": []}
    credit = _section("data credit", A.data_credit, errors, inp) or {}

    workflows = []
    if pop or da or ds:
        workflows.append("Immune Profiling & Differential Abundance")
    if batch:
        workflows.append("Multi-batch Harmonization & QC")

    generated_at = started.isoformat(timespec="seconds")
    context = {
        "title": title, "generated_at": generated_at, "app_version": app_version, "tools": tools,
        "ov": ov, "qc": qc, "pop": pop, "da": da, "ds": ds, "batch": batch, "methods": meth, "credit": credit,
        "figs": figs, "errors": errors, "workflows": workflows, "fdr_override": fdr,
        "inputs_found": {k: len(v) for k, v in inp.found.items()},
    }
    R.write_index(report_dir, context)

    # Tables.
    written = []
    if pop and pop.get("composition") is not None:
        pop["composition"].to_csv(tables_dir / "composition.csv", index=False)
        written.append("composition.csv")
    if da and not da.get("failed"):
        A.da_table(da).to_csv(tables_dir / "da_significant.csv", index=False)
        written.append("da_significant.csv")
    if batch and batch.get("scores"):
        A.batch_scores_table(batch).to_csv(tables_dir / "batch_scores.csv", index=False)
        written.append("batch_scores.csv")
    if qc:
        R.qc_table(qc).to_csv(tables_dir / "qc.csv", index=False)
        written.append("qc.csv")
    metrics = R.metrics_json(context, own)
    with open(tables_dir / "metrics.json", "w") as fh:
        json.dump(A.clean(metrics), fh, indent=2)
    written.append("metrics.json")

    finished = _now()
    provenance = {
        "app": APP_NAME, "version": app_version, "tools": tools, "parameters": parameters,
        "startedAt": generated_at, "finishedAt": finished.isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - t_start, 2),
        "runtime": "lambda" if os.environ.get("AWS_LAMBDA_FUNCTION_NAME") else "container",
        "inputsFound": {k: [inp.rel(p) for p in v] for k, v in inp.found.items()},
        "assessmentLabels": list(inp.assessments), "normalizationMethods": list(inp.normalization),
        "provenanceFiles": [inp.rel(p) for p, _ in inp.provenances],
        "sections": {k: bool(v) for k, v in (("overview", ov), ("qc", qc), ("populations", pop),
                                             ("differential_abundance", da), ("differential_state", ds),
                                             ("batch_effects", batch))},
        "sectionErrors": errors,
        "outputs": {"report": sorted(str(p.relative_to(report_dir)) for p in report_dir.rglob("*") if p.is_file())
                    + ["provenance.json"], "tables": written},
    }
    with open(report_dir / "provenance.json", "w") as fh:
        json.dump(A.clean(provenance), fh, indent=2)
    with open(tables_dir / "provenance.json", "w") as fh:
        json.dump(A.clean(provenance), fh, indent=2)
    log.info("wrote %s (%d section error(s)) and %d table(s) in %.1fs", report_dir / "index.html", len(errors),
             len(written), time.monotonic() - t_start)
    return 0


if __name__ == "__main__":
    sys.exit(main())
