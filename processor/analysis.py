"""
What the report says: per section, the tables, the numbers and the observations
(plain sentences drawn from the data) and caveats (what to keep in mind before
trusting a number). Figures are drawn in figures.py from what this returns.

The text describes what the data show; it does not tell the reader what to do.
"""
from __future__ import annotations

import json
import math
import re
from typing import Any

import numpy as np
import pandas as pd

from processor.inputs import Inputs


# ---------------------------------------------------------------- formatting

def pct(frac: float | None, digits: int | None = None) -> str:
    """A 0-1 fraction as a percentage with sensible precision: 19%, 1.8%, 0.06%."""
    if frac is None or (isinstance(frac, float) and math.isnan(frac)):
        return "–"
    v = 100.0 * float(frac)
    if digits is not None:
        return f"{v:.{digits}f}%"
    a = abs(v)
    if a == 0:
        return "0%"
    if a >= 10:
        return f"{v:.0f}%"
    if a >= 1:
        return f"{v:.1f}%"
    if a >= 0.1:
        return f"{v:.2f}%"
    return f"{v:.1g}%" if a >= 0.001 else "<0.001%"


def pval(p: float | None) -> str:
    """6e-93, 0.003, 0.04."""
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "–"
    p = float(p)
    if p == 0:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.0e}".replace("e-0", "e-")
    return f"{p:.2g}"


def num(value: Any, digits: int = 2, signed: bool = False) -> str:
    if value is None:
        return "–"
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if math.isnan(f):
        return "–"
    return f"{f:{'+' if signed else ''}.{digits}f}"


def count(value: Any) -> str:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return "–"
    if math.isnan(f):
        return "–"
    return f"{int(round(f)):,}"


def truthy(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in ("true", "1", "yes", "t")


def joined(items: list[str], conj: str = "and") -> str:
    items = [str(i) for i in items]
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + f" {conj} " + items[-1]


def clean(obj: Any) -> Any:
    """JSON-safe: NaN -> None, numpy scalars -> Python."""
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return None if math.isnan(f) or math.isinf(f) else f
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def mc_name(mc: Any, label: Any) -> str:
    return f"MC{mc} {label}" if label not in (None, "", mc) and not (isinstance(label, float) and math.isnan(label)) else f"MC{mc}"


# ---------------------------------------------------------------- overview

DONOR_COLUMNS = ("donor", "donor_id", "subject", "subject_id", "patient", "patient_id", "participant", "individual")


def donor_design(inp: Inputs) -> dict[str, Any] | None:
    """Which donor each sample comes from: a donor column of samples.csv, else (for the demo dataset) the donors
    recorded from its $TUBENAME keywords. Returns per condition the donors and the batches they sit in."""
    from processor import citations as C

    s = inp.samples
    if s is None or "condition" not in s:
        return None
    col = next((c for c in s.columns if c.lower() in DONOR_COLUMNS), None)
    if col:
        donor = s[col].astype(str)
        source = f"the '{col}' column of samples.csv"
    elif "file" in s:
        donor = s["file"].astype(str).map(C.demo_donor)
        if donor.isna().any():
            return None
        source = "the demo files' $TUBENAME keywords, as recorded in this app"
    else:
        return None
    df = pd.DataFrame({"condition": s["condition"].astype(str), "donor": donor,
                       "batch": s["batch"].astype(str) if "batch" in s else "NA"})
    per = {}
    for cond, g in df.groupby("condition", sort=False):
        per[cond] = {d: sorted(gg["batch"].unique()) for d, gg in g.groupby("donor")}
    conds = list(per)
    shared = set.intersection(*[set(v) for v in per.values()]) if len(per) > 1 else set()
    n_samples = df.groupby("condition").size().to_dict()
    return {"source": source, "per_condition": per, "conditions": conds,
            "donors": sorted(df["donor"].unique()), "shared": sorted(shared),
            "condition_is_donor": len(per) > 1 and not shared,
            "n_samples": {k: int(v) for k, v in n_samples.items()},
            "batch_varies": {c: len(v) > 1 for c, v in per.items()}}


def donor_text(d: dict[str, Any]) -> str:
    parts = []
    for cond, donors in d["per_condition"].items():
        if len(donors) == 1:
            (dn, _), = donors.items()
            parts.append(f"{cond}: donor {dn} in every batch")
        else:
            parts.append(f"{cond}: " + joined([f"donor {dn} ({', '.join(b)})" for dn, b in donors.items()]))
    return "; ".join(parts)


def _cofactor(value: Any) -> str:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return ""
    return "" if math.isnan(f) else f"{f:g}"


def overview(inp: Inputs) -> dict[str, Any] | None:
    s = inp.samples
    prep = inp.prep or {}
    if s is None and inp.panel is None and not prep:
        return None
    out: dict[str, Any] = {"observations": [], "caveats": []}
    if s is not None:
        n = len(s)
        batches = sorted(s["batch"].unique()) if "batch" in s else []
        conds = sorted(s["condition"].unique()) if "condition" in s else []
        ref = (inp.da_summary or {}).get("reference") or (inp.ds_summary or {}).get("reference")
        if ref in conds:  # the reference level first
            conds = [ref] + [c for c in conds if c != ref]
        out.update(n_samples=n, batches=batches, conditions=conds)
        if batches and conds:
            tab = pd.crosstab(s["batch"], s["condition"]).reindex(index=batches, columns=conds, fill_value=0)
            out["design"] = {"columns": conds, "rows": [{"batch": b, "counts": [int(tab.loc[b, c]) for c in conds],
                                                         "total": int(tab.loc[b].sum())} for b in batches],
                             "totals": [int(tab[c].sum()) for c in conds], "total": n}
            missing = [(b, c) for b in batches for c in conds if tab.loc[b, c] == 0]
            if len(conds) > 1 and not missing:
                out["observations"].append(
                    f"{n} samples in {len(batches)} batches; every batch holds every condition "
                    f"({joined([f'{int(tab[c].sum())} {c}' for c in conds])}), so condition and batch are not confounded.")
            elif len(conds) > 1:
                out["observations"].append(f"{n} samples in {len(batches)} batches and {len(conds)} conditions.")
                out["caveats"].append(
                    "Not every batch holds every condition (" +
                    joined([f"{b} has no {c}" for b, c in missing[:6]]) + ("…" if len(missing) > 6 else "") +
                    "); a difference between conditions can partly be a difference between batches.")
            else:
                out["observations"].append(f"{n} samples in {len(batches)} batch(es), one condition ({joined(conds)}).")
        extra = [c for c in s.columns if c not in ("sample_id", "file", "batch", "condition", "n_events")]
        out["extra_columns"] = extra
        if "n_events" in s:
            ev = pd.to_numeric(s["n_events"], errors="coerce")
            out["events"] = {"total": float(ev.sum()), "median": float(ev.median()), "min": float(ev.min()),
                             "max": float(ev.max())}
            small = s[ev < 0.2 * ev.median()]
            out["observations"].append(
                f"{count(ev.sum())} events in all; per sample median {count(ev.median())} "
                f"(range {count(ev.min())}–{count(ev.max())}).")
            if len(small):
                out["caveats"].append(
                    joined([f"{r.sample_id} ({count(r.n_events)} events)" for r in small.itertuples()]) +
                    (" has" if len(small) == 1 else " have") +
                    " fewer than a fifth of the median events per sample; proportions of rare populations in "
                    + ("it" if len(small) == 1 else "them") + " rest on few cells.")
        out["samples"] = s.to_dict("records")
    if inp.panel is not None:
        p = inp.panel.copy()
        typ = p["type"] if "type" in p else pd.Series(["marker"] * len(p))
        markers = p[typ == "marker"]
        rows = []
        for r in markers.itertuples():
            rows.append({"marker": r.marker, "channel": r.channel,
                         "clustering": truthy(getattr(r, "clustering", False)),
                         "state": truthy(getattr(r, "state", False)),
                         "cofactor": _cofactor(getattr(r, "cofactor", ""))})
        out["panel"] = rows
        out["n_markers"] = len(rows)
        out["n_clustering"] = sum(r["clustering"] for r in rows)
        out["n_state"] = sum(r["state"] for r in rows)
        out["other_channels"] = [str(c) for c in p[typ != "marker"]["channel"]]
        cof = sorted({r["cofactor"] for r in rows if r["cofactor"]})
        out["observations"].append(
            f"{len(rows)} markers kept ({out['n_clustering']} used for clustering, {out['n_state']} state markers"
            + (f"; arcsinh cofactor {joined(cof)}" if cof else "") + ")" +
            (f", plus {len(out['other_channels'])} scatter/time/other channels." if out["other_channels"] else "."))
    dropped: dict[str, dict[str, Any]] = {}
    for sid, rec in (prep.get("samples") or {}).items():
        for d in rec.get("channels_dropped") or []:
            ch = d.get("channel")
            if ch and ch not in dropped:
                dropped[ch] = {"channel": ch, "marker": d.get("marker", ""), "reason": d.get("reason", "")}
    for ch in prep.get("channelsDropped") or []:
        dropped.setdefault(ch, {"channel": ch, "marker": "", "reason": ""})
    out["dropped"] = list(dropped.values())
    if dropped:
        out["observations"].append(
            f"{len(dropped)} channel(s) dropped because the panel differs between samples: " +
            joined([f"{d['channel']}" + (f" ({d['marker']})" if d["marker"] else "") for d in dropped.values()]) + ".")
    dd = donor_design(inp)
    out["donors"] = dd
    if dd:
        out["observations"].append(f"Donors (from {dd['source']}): {donor_text(dd)}.")
        if dd["condition_is_donor"]:
            out["caveats"].append(
                "No donor appears in more than one condition, so every difference between conditions is also a "
                "difference between donors; the data cannot separate the two.")
        few = {c: len(v) for c, v in dd["per_condition"].items() if len(v) < dd["n_samples"].get(c, 0)}
        if few:
            out["caveats"].append(
                "Samples within a condition are repeated measurements of few donors (" +
                joined([f"{dd['n_samples'][c]} {c} samples from {n} donor{'s' if n > 1 else ''}" for c, n in few.items()]) +
                "). Tests that treat each sample as an independent replicate measure reproducibility across "
                "batches and aliquots, not variation between people.")
        varies = [c for c, v in dd["batch_varies"].items() if v]
        if varies:
            out["caveats"].append(
                "In " + joined(varies) + " samples the donor changes with the batch, so a batch difference there is "
                "partly a donor difference.")
    out["cytometer"] = prep.get("cytometerType")
    instruments = sorted({str(r.get("cytometer")) for r in (prep.get("samples") or {}).values() if r.get("cytometer")})
    out["instruments"] = instruments
    return out


# ---------------------------------------------------------------- QC

def qc(inp: Inputs) -> dict[str, Any] | None:
    pq, db = inp.peacoqc, inp.doublets
    if pq is None and db is None:
        return None
    out: dict[str, Any] = {"observations": [], "caveats": []}
    ids: list[str] = []
    for df in (inp.samples, pq, db):
        if df is not None:
            ids += [i for i in df["sample_id"] if i not in ids]
    rows = []
    for sid in ids:
        row: dict[str, Any] = {"sample_id": sid}
        if inp.samples is not None:
            m = inp.samples[inp.samples["sample_id"] == sid]
            if len(m):
                row["batch"] = m.iloc[0].get("batch")
                row["condition"] = m.iloc[0].get("condition")
                row["n_events"] = m.iloc[0].get("n_events")
        if pq is not None:
            m = pq[pq["sample_id"] == sid]
            if len(m):
                r = m.iloc[0]
                row.update(qc_in=r.get("events_in"), qc_out=r.get("events_out"), qc_pct=r.get("pct_removed"),
                           qc_flagged=truthy(r.get("flagged")), qc_status=r.get("status"),
                           qc_notes=r.get("notes") if isinstance(r.get("notes"), str) else "")
        if db is not None:
            m = db[db["sample_id"] == sid]
            if len(m):
                r = m.iloc[0]
                row.update(db_in=r.get("events_in"), db_n=r.get("doublets"), db_pct=r.get("pct_doublets"),
                           db_out=r.get("events_out"), db_status=r.get("status"),
                           db_notes=r.get("notes") if isinstance(r.get("notes"), str) else "")
        start = row.get("n_events", row.get("qc_in", row.get("db_in")))
        end = row.get("db_out", row.get("qc_out"))
        if start is not None and end is not None and not pd.isna(start) and not pd.isna(end) and float(start) > 0:
            row["kept_pct"] = 100.0 * float(end) / float(start)
        rows.append(row)
    out["rows"] = rows
    if pq is not None:
        pct_rm = pd.to_numeric(pq["pct_removed"], errors="coerce")
        worst = pq.loc[pct_rm.idxmax()] if pct_rm.notna().any() else None
        flagged = [str(s) for s, f in zip(pq["sample_id"], pq.get("flagged", [])) if truthy(f)]
        failed = [str(s) for s, st in zip(pq["sample_id"], pq.get("status", [])) if str(st) == "failed"]
        out["peacoqc"] = {"median_pct": float(pct_rm.median()), "min_pct": float(pct_rm.min()),
                          "max_pct": float(pct_rm.max()), "flagged": flagged, "failed": failed,
                          "events_in": float(pd.to_numeric(pq["events_in"], errors="coerce").sum()),
                          "events_out": float(pd.to_numeric(pq["events_out"], errors="coerce").sum())}
        limit = None
        prov = _provenance_of(inp, "processor-peacoqc")
        if prov:
            limit = (prov.get("parameters") or {}).get("MAX_REMOVED_PCT")
        out["peacoqc"]["limit"] = limit
        out["observations"].append(
            f"PeacoQC removed a median {pct_rm.median():.1f}% of events per sample (range {pct_rm.min():.1f}–"
            f"{pct_rm.max():.1f}%" + (f", most in {worst['sample_id']}" if worst is not None else "") + ").")
        if inp.samples is not None and "batch" in inp.samples:
            byb = _with_meta(pq[["sample_id"]].assign(p=pct_rm.values), inp).groupby("batch")["p"].median()
            high = byb[(byb > 2 * pct_rm.median()) & (byb > 5)]
            if len(byb) > 1 and len(high):
                out["observations"].append(
                    joined([f"{b} (median {v:.1f}%)" for b, v in high.items()]) +
                    f" lost more events to QC than the other batches (overall median {pct_rm.median():.1f}%), "
                    "so the instability is a property of " + ("those runs" if len(high) > 1 else "that run") +
                    " rather than of single files.")
        if flagged:
            out["observations"].append(
                f"{len(flagged)} file(s) flagged" + (f" for losing more than {limit}% of events" if limit else "") +
                f": {joined(flagged)}. Flagged files were still passed on; their PeacoQC plots are below.")
        else:
            out["observations"].append("No file was flagged" + (f" (limit {limit}% removed)." if limit else "."))
        if failed:
            out["caveats"].append(f"PeacoQC failed on {joined(failed)}; those samples did not continue.")
        notes = [str(n) for n in pq.get("notes", []) if isinstance(n, str) and n.strip()]
        if notes:
            common = pd.Series(notes).value_counts()
            for note, k in common.items():
                out["caveats"].append(f"PeacoQC note on {k} sample(s): “{note.strip()}”")
        plots = []
        if inp.peacoqc_dir is not None:
            for sid in flagged or []:
                p = inp.peacoqc_dir / "plots" / f"{sid}.png"
                if p.is_file():
                    plots.append((sid, p))
        out["peacoqc_plots"] = plots[:6]
    if db is not None:
        dp = pd.to_numeric(db["pct_doublets"], errors="coerce")
        failed = [str(s) for s, st in zip(db["sample_id"], db.get("status", [])) if str(st) == "failed"]
        passed = [str(s) for s, st in zip(db["sample_id"], db.get("status", [])) if str(st) == "passed_through"]
        out["doublets"] = {"median_pct": float(dp.median()), "min_pct": float(dp.min()), "max_pct": float(dp.max()),
                           "total": float(pd.to_numeric(db["doublets"], errors="coerce").sum()),
                           "failed": failed, "passed_through": passed}
        worst = db.loc[dp.idxmax()] if dp.notna().any() else None
        out["observations"].append(
            f"Cleanet called a median {dp.median():.2f}% of events doublets (range {dp.min():.2f}–{dp.max():.2f}%"
            + (f", most in {worst['sample_id']}" if worst is not None else "") + ").")
        if passed:
            out["caveats"].append(f"Cleanet passed {joined(passed)} through unchecked (see notes).")
        if failed:
            out["caveats"].append(f"Cleanet failed on {joined(failed)}.")
    return out


def _provenance_of(inp: Inputs, app: str) -> dict[str, Any] | None:
    for _, prov in inp.provenances:
        if prov.get("app") == app:
            return prov
    return None


# ---------------------------------------------------------------- populations

def _label_map(inp: Inputs) -> dict[int, str]:
    if inp.labels is not None:
        return {int(r.metacluster): str(r.label) for r in inp.labels.itertuples()}
    if inp.counts is not None and "label" in inp.counts:
        return {int(m): str(l) for m, l in inp.counts[["metacluster", "label"]].drop_duplicates().values}
    return {}


def _with_meta(df: pd.DataFrame, inp: Inputs) -> pd.DataFrame:
    if inp.samples is None:
        df = df.copy()
        df["batch"] = "NA"
        df["condition"] = "NA"
        return df
    cols = [c for c in ("sample_id", "batch", "condition") if c in inp.samples]
    df = df.merge(inp.samples[cols], on="sample_id", how="left")
    for c in ("batch", "condition"):
        if c not in df:
            df[c] = "NA"
        df[c] = df[c].fillna("NA").astype(str)
    return df


def composition(inp: Inputs) -> pd.DataFrame | None:
    """sample x label fractions (labels summed over their metaclusters), with batch and condition."""
    if inp.counts is None:
        return None
    c = inp.counts.copy()
    if "label" not in c:
        lm = _label_map(inp)
        c["label"] = c["metacluster"].map(lambda m: lm.get(int(m), f"MC{m}"))
    tot = c.groupby("sample_id")["n"].transform("sum")
    c["frac"] = c["n"] / tot
    wide = c.pivot_table(index="sample_id", columns="label", values="frac", aggfunc="sum", fill_value=0.0)
    order = wide.mean().sort_values(ascending=False).index
    wide = wide[order].reset_index()
    wide = _with_meta(wide, inp)
    return wide[["sample_id", "batch", "condition"] + list(order)]


def populations(inp: Inputs, reference_level: str | None) -> dict[str, Any] | None:
    if inp.counts is None and inp.labels is None and inp.medians is None and not inp.umap_figures:
        return None
    out: dict[str, Any] = {"observations": [], "caveats": []}
    lm = _label_map(inp)
    prov = inp.clusters_provenance or _provenance_of(inp, "processor-flowsom") or {}
    total_cells = None
    rows = []
    if inp.labels is not None:
        lab = inp.labels
        total_cells = float(pd.to_numeric(lab.get("n_cells", pd.Series(dtype=float)), errors="coerce").sum()) or None
        for r in lab.itertuples():
            rows.append({"metacluster": int(r.metacluster), "label": str(r.label),
                         "n_cells": getattr(r, "n_cells", None), "frac": getattr(r, "frac", None),
                         "profile": getattr(r, "profile", "")})
    elif inp.counts is not None:
        g = inp.counts.groupby("metacluster")["n"].sum()
        total_cells = float(g.sum())
        for mc, n in g.items():
            rows.append({"metacluster": int(mc), "label": lm.get(int(mc), f"MC{mc}"), "n_cells": n,
                         "frac": n / total_cells, "profile": ""})
    rows.sort(key=lambda r: -(float(r["frac"]) if r["frac"] is not None else 0))
    out["metaclusters"] = rows
    out["total_cells"] = total_cells
    comp = composition(inp)
    out["composition"] = comp
    if rows:
        n_mc = len(rows)
        labels = sorted({r["label"] for r in rows})
        top = rows[:3]
        out["observations"].append(
            f"{n_mc} metaclusters with {len(labels)} distinct labels"
            + (f" over {count(total_cells)} cells" if total_cells else "") + ". Largest: " +
            joined([f"{mc_name(r['metacluster'], r['label'])} ({pct(r['frac'])})" for r in top]) + ".")
        dup = {}
        for r in rows:
            dup.setdefault(r["label"], []).append(r["metacluster"])
        multi = {k: sorted(v) for k, v in dup.items() if len(v) > 1}
        if multi:
            out["observations"].append(
                "Some labels cover several metaclusters, which differ in markers the label does not name: " +
                joined([f"{k} (MC{', MC'.join(str(m) for m in v)})" for k, v in list(multi.items())[:5]]) + ".")
        unassigned = [r for r in rows if "unassigned" in r["label"].lower()]
        if unassigned:
            out["caveats"].append(
                f"{len(unassigned)} metacluster(s) matched no labelling rule ('unassigned', "
                f"{pct(sum(float(r['frac'] or 0) for r in unassigned))} of cells).")
    # Labels are a heuristic.
    method = prov.get("method") if isinstance(prov.get("method"), dict) else {}
    if rows:
        out["caveats"].append(
            "Labels are a transparent heuristic, not a classifier: each metacluster's median is scaled per marker "
            "across metaclusters and matched against a fixed rule table"
            + (" (flowsom's provenance.json lists the rules)" if method.get("rules") else "")
            + ". Read them as names for the marker profiles in the heatmap, and check the profile before trusting one.")
    panel_markers = {str(m).upper().replace("-", "") for m in (inp.panel["marker"] if inp.panel is not None else
                                                               (inp.medians.columns if inp.medians is not None else []))}
    if rows and any("treg" in r["label"].lower() for r in rows) and panel_markers and not panel_markers & {"CD127", "IL7R", "IL7RA"}:
        out["caveats"].append(
            "The panel has no CD127, so CD25+FOXP3+ CD4 T cells cannot be told apart from recently activated "
            "conventional CD4 T cells (which also up-regulate CD25 and can express FOXP3); hence 'Treg or activated'.")
    dbl = [r for r in rows if "doublet" in r["label"].lower()]
    if dbl:
        f = sum(float(r["frac"] or 0) for r in dbl)
        text = (f"{joined([mc_name(r['metacluster'], r['label']) for r in dbl])} "
                f"({pct(f)} of cells) co-express lineage markers of two cell types")
        out["caveats"].append(text + (": doublets that escaped the doublet step, or a rare true population."
                                      if inp.doublets is not None else
                                      "; no doublet-removal output was among the inputs, so these may be doublets."))
    # Composition by condition.
    if comp is not None and comp["condition"].nunique() > 1:
        label_cols = [c for c in comp.columns if c not in ("sample_id", "batch", "condition")]
        means = comp.groupby("condition")[label_cols].mean()
        out["composition_means"] = {"conditions": list(means.index),
                                    "rows": [{"label": l, "means": [float(means.loc[c, l]) for c in means.index]}
                                             for l in label_cols]}
    elif comp is not None:
        label_cols = [c for c in comp.columns if c not in ("sample_id", "batch", "condition")]
        out["composition_means"] = {"conditions": ["all"],
                                    "rows": [{"label": l, "means": [float(comp[l].mean())]} for l in label_cols]}
    out["umap_figures"] = inp.umap_figures
    ep = inp.embedding_provenance or _provenance_of(inp, "processor-cytometry-umap") or {}
    if inp.umap_figures:
        n = ep.get("cells_embedded")
        params = ep.get("parameters") or {}
        out["umap_note"] = (f"UMAP of {count(n)} cells" if n else "UMAP") + (
            f" sampled evenly across samples, on the {len(params.get('MARKERS_used') or [])} clustering markers"
            if params.get("MARKERS_used") else "") + (
            f" (n_neighbors {params.get('N_NEIGHBORS')}, min_dist {params.get('MIN_DIST')})" if params.get("N_NEIGHBORS") else "") + \
            ". Distances between distant islands in a UMAP are not meaningful; neighbourhoods are."
    out["clustering_markers"] = prov.get("clustering_markers") or (
        [r["marker"] for r in (overview_panel(inp) or []) if r["clustering"]])
    out["state_markers"] = prov.get("state_markers") or [r["marker"] for r in (overview_panel(inp) or []) if r["state"]]
    out["star_tree"] = (inp.clusters_dir / "star_tree.png") if inp.clusters_dir and (inp.clusters_dir / "star_tree.png").is_file() else None
    out["flowsom_heatmap"] = (inp.clusters_dir / "heatmap.png") if inp.clusters_dir and (inp.clusters_dir / "heatmap.png").is_file() else None
    return out


def overview_panel(inp: Inputs) -> list[dict[str, Any]] | None:
    if inp.panel is None:
        return None
    p = inp.panel
    typ = p["type"] if "type" in p else pd.Series(["marker"] * len(p))
    return [{"marker": r.marker, "clustering": truthy(getattr(r, "clustering", False)),
             "state": truthy(getattr(r, "state", False))} for r in p[typ == "marker"].itertuples()]


# ---------------------------------------------------------------- differential abundance

def _levels(comparison: str) -> tuple[str, str]:
    m = re.match(r"^(.*) vs (.*)$", str(comparison))
    return (m.group(1), m.group(2)) if m else (str(comparison), "")


def batch_consistency(inp: Inputs, mc: int, level: str, ref: str) -> dict[str, Any] | None:
    """Per batch: the cluster's mean fraction in `level` minus in `ref`; how many batches agree with the overall sign,
    and whether one batch carries most of the difference."""
    if inp.counts is None or inp.samples is None or "batch" not in inp.samples or "condition" not in inp.samples:
        return None
    c = inp.counts
    tot = c.groupby("sample_id")["n"].sum()
    sub = c[c["metacluster"].astype(int) == int(mc)].groupby("sample_id")["n"].sum()
    frac = (sub.reindex(tot.index, fill_value=0) / tot).rename("frac").reset_index()
    frac = _with_meta(frac, inp)
    per = frac.groupby(["batch", "condition"])["frac"].mean().unstack()
    if level not in per or ref not in per:
        return None
    diff = (per[level] - per[ref]).dropna()
    if diff.empty:
        return None
    overall = float(diff.mean())
    sign = np.sign(overall)
    agree = int((np.sign(diff) == sign).sum())
    share = float(diff.abs().max() / diff.abs().sum()) if diff.abs().sum() > 0 else 0.0
    top_batch = str(diff.abs().idxmax())
    # Leave the top batch out: does the direction survive?
    rest = diff.drop(top_batch)
    survives = bool(len(rest) == 0 or np.sign(rest.mean()) == sign)
    return {"n_batches": int(len(diff)), "agree": agree, "share_top": share, "top_batch": top_batch,
            "survives_without_top": survives, "per_batch": {str(k): float(v) for k, v in diff.items()}}


def differential_abundance(inp: Inputs, fdr: float | None) -> dict[str, Any] | None:
    if inp.da is None:
        return None
    da = inp.da.copy()
    summ = inp.da_summary or {}
    if summ.get("status") == "failed":
        return {"failed": summ.get("error") or "diffcyt's DA test failed", "observations": [], "caveats": []}
    threshold = fdr if fdr is not None else float(summ.get("fdr") or 0.05)
    da["p_adj"] = pd.to_numeric(da["p_adj"], errors="coerce")
    da["logFC"] = pd.to_numeric(da["logFC"], errors="coerce")
    da["sig"] = da["p_adj"] < threshold
    out: dict[str, Any] = {"observations": [], "caveats": [], "fdr": threshold, "method": summ.get("method"),
                           "block": summ.get("block"), "comparisons": []}
    for comp, g in da.groupby("comparison", sort=False):
        level, ref = _levels(comp)
        tested = g[g["tested"].map(truthy)] if "tested" in g else g
        sig = g[g["sig"]].sort_values("p_adj")
        rows = []
        for r in sig.itertuples():
            fl = getattr(r, f"mean_frac_{level}", None) if f"mean_frac_{level}" in g else None
            fr = getattr(r, f"mean_frac_{ref}", None) if f"mean_frac_{ref}" in g else None
            cons = batch_consistency(inp, r.metacluster, level, ref)
            row = {"metacluster": int(r.metacluster), "label": str(r.label), "logFC": float(r.logFC),
                   "p_adj": float(r.p_adj), "frac_level": fl, "frac_ref": fr, "consistency": cons, "notes": []}
            if cons:
                if cons["agree"] < cons["n_batches"]:
                    row["notes"].append(f"same direction in {cons['agree']} of {cons['n_batches']} batches")
                if cons["n_batches"] >= 3 and (cons["share_top"] > 0.5 or not cons["survives_without_top"]):
                    row["notes"].append(f"{cons['top_batch']} carries {cons['share_top'] * 100:.0f}% of the difference"
                                        + ("" if cons["survives_without_top"] else "; without it the direction reverses"))
            if fl is not None and fr is not None and max(fl, fr) < 0.005:
                row["notes"].append("rare in both groups (<0.5% of cells)")
            if "doublet" in str(r.label).lower():
                row["notes"].append("doublet-like cluster")
            rows.append(row)
        up = int((sig["logFC"] > 0).sum())
        down = int((sig["logFC"] < 0).sum())
        entry = {"comparison": comp, "level": level, "reference": ref, "tested": int(len(tested)),
                 "clusters": int(len(g)), "significant": rows, "up": up, "down": down}
        out["comparisons"].append(entry)
        out["observations"].append(
            f"{comp}: {len(rows)} of {len(tested)} tested metaclusters differ in abundance at FDR < {threshold:g} "
            f"({up} more abundant in {level}, {down} less).")
        for row in rows[:4]:
            if row["frac_level"] is not None and row["frac_ref"] is not None:
                verb = "rose" if row["frac_level"] > row["frac_ref"] else "fell"
                cons = row["consistency"]
                tail = ""
                if cons:
                    tail = (f"; in all {cons['n_batches']} batches" if cons["agree"] == cons["n_batches"]
                            else f"; in {cons['agree']} of {cons['n_batches']} batches")
                out["observations"].append(
                    f"{mc_name(row['metacluster'], row['label'])} {verb} from {pct(row['frac_ref'])} to "
                    f"{pct(row['frac_level'])} of cells (log2FC {row['logFC']:+.2f}, FDR {pval(row['p_adj'])}{tail}).")
        disagree = [r for r in rows if r["consistency"] and r["consistency"]["agree"] < r["consistency"]["n_batches"]]
        if disagree:
            out["caveats"].append(
                "Not every batch goes the same way for " +
                joined([f"{mc_name(r['metacluster'], r['label'])} ({r['consistency']['agree']} of "
                        f"{r['consistency']['n_batches']})" for r in disagree[:6]]) +
                ("…" if len(disagree) > 6 else "") + "; the grey lines in the figure show which.")
        dominated = [r for r in rows if r["consistency"] and r["consistency"]["n_batches"] >= 3
                     and r["consistency"]["share_top"] > 0.5]
        if dominated:
            out["caveats"].append(
                "One batch carries most of the difference for " +
                joined([f"{mc_name(r['metacluster'], r['label'])} ({r['consistency']['top_batch']}, "
                        f"{r['consistency']['share_top'] * 100:.0f}%)" for r in dominated[:6]]) +
                ("…" if len(dominated) > 6 else "") +
                "; these results rest largely on that batch's samples.")
    dd = donor_design(inp)
    if dd and dd["condition_is_donor"]:
        out["caveats"].append(
            "The conditions come from different donors (" + donor_text(dd) + "), so these differences mix the "
            "effect of the condition with differences between donors, and the samples in a condition are "
            "aliquots and batches of the same donor(s) rather than independent people. The FDRs are therefore "
            "much smaller than a study with one sample per donor would give.")
    big_up = [r for c in out["comparisons"] for r in c["significant"]
              if r["frac_level"] is not None and r["frac_ref"] is not None and r["frac_level"] - r["frac_ref"] > 0.1]
    if big_up:
        r = big_up[0]
        out["caveats"].append(
            f"Proportions are compositional: {mc_name(r['metacluster'], r['label'])} grew by "
            f"{(r['frac_level'] - r['frac_ref']) * 100:.0f} percentage points, which by itself lowers every other "
            "population's share. A fall in proportion is not necessarily a fall in absolute numbers.")
    if summ.get("method") == "edgeR" or summ.get("logFC"):
        out["caveats"].append(
            "The log2 fold change is the " + str(summ.get("logFC") or summ.get("method")) +
            (f", adjusted for {summ.get('block')} as a {summ.get('block_effect') or 'fixed'} effect" if summ.get("block") else "") +
            "; the percentages beside it are plain means of the per-sample proportions, so the two need not match exactly.")
    out["volcano"] = (inp.da_dir / "volcano.png") if inp.da_dir and (inp.da_dir / "volcano.png").is_file() else None
    return out


def da_table(da_section: dict[str, Any] | None) -> pd.DataFrame:
    rows = []
    for c in (da_section or {}).get("comparisons", []):
        for r in c["significant"]:
            cons = r["consistency"] or {}
            rows.append({"comparison": c["comparison"], "metacluster": r["metacluster"], "label": r["label"],
                         "logFC": r["logFC"], "p_adj": r["p_adj"], f"mean_frac_reference": r["frac_ref"],
                         "mean_frac_level": r["frac_level"], "reference": c["reference"], "level": c["level"],
                         "batches_same_direction": cons.get("agree"), "n_batches": cons.get("n_batches"),
                         "notes": "; ".join(r["notes"])})
    return pd.DataFrame(rows, columns=["comparison", "metacluster", "label", "logFC", "p_adj", "mean_frac_reference",
                                       "mean_frac_level", "reference", "level", "batches_same_direction", "n_batches",
                                       "notes"])


# ---------------------------------------------------------------- differential state

def differential_state(inp: Inputs, fdr: float | None, top_n: int = 15) -> dict[str, Any] | None:
    if inp.ds is None:
        return None
    ds = inp.ds.copy()
    summ = inp.ds_summary or {}
    if summ.get("status") == "failed":
        return {"failed": summ.get("error") or "diffcyt's DS test failed", "observations": [], "caveats": []}
    threshold = fdr if fdr is not None else float(summ.get("fdr") or 0.05)
    ds["p_adj"] = pd.to_numeric(ds["p_adj"], errors="coerce")
    ds["logFC"] = pd.to_numeric(ds["logFC"], errors="coerce")
    ds["sig"] = ds["p_adj"] < threshold
    out: dict[str, Any] = {"observations": [], "caveats": [], "fdr": threshold, "method": summ.get("method"),
                           "comparisons": []}
    for comp, g in ds.groupby("comparison", sort=False):
        level, ref = _levels(comp)
        tested = g[g["tested"].map(truthy)] if "tested" in g else g
        sig = g[g["sig"]]
        by_p = sig.sort_values(["p_adj", "logFC"], ascending=[True, False]).head(top_n)
        by_effect = sig.reindex(sig["logFC"].abs().sort_values(ascending=False).index).head(top_n)

        def rows_of(df):
            res = []
            for r in df.itertuples():
                res.append({"metacluster": int(r.metacluster), "label": str(r.label), "marker": str(r.marker),
                            "logFC": float(r.logFC), "p_adj": float(r.p_adj),
                            "ref": getattr(r, f"mean_median_{ref}", None) if f"mean_median_{ref}" in g else None,
                            "level": getattr(r, f"mean_median_{level}", None) if f"mean_median_{level}" in g else None})
            return res
        per_marker = sig.groupby("marker").size().sort_values(ascending=False)
        entry = {"comparison": comp, "level": level, "reference": ref, "tested": int(len(tested)),
                 "significant": int(len(sig)), "up": int((sig["logFC"] > 0).sum()), "down": int((sig["logFC"] < 0).sum()),
                 "top_p": rows_of(by_p), "top_effect": rows_of(by_effect),
                 "per_marker": [{"marker": m, "n": int(n)} for m, n in per_marker.items()]}
        out["comparisons"].append(entry)
        out["observations"].append(
            f"{comp}: {len(sig)} of {len(tested)} tested cluster × marker pairs differ at FDR < {threshold:g} "
            f"({entry['up']} higher in {level}, {entry['down']} lower).")
        if len(by_effect):
            r = entry["top_effect"][0]
            out["observations"].append(
                f"Largest shift: {r['marker']} in {mc_name(r['metacluster'], r['label'])}, "
                f"{num(r['ref'])} → {num(r['level'])} (median arcsinh; difference {r['logFC']:+.2f}, FDR {pval(r['p_adj'])}).")
        if len(per_marker):
            out["observations"].append(
                "Markers with the most significant clusters: " +
                joined([f"{m} ({n})" for m, n in list(per_marker.items())[:5]]) + ".")
        if len(tested) and len(sig) > 0.5 * len(tested):
            out["caveats"].append(
                f"More than half of the pairs tested are significant in {comp}; with this many hits the effect size "
                "(the difference in median expression) says more than the p-value. The table is ranked both ways.")
    dd = donor_design(inp)
    if dd and dd["condition_is_donor"]:
        out["caveats"].append("As for abundance, the conditions come from different donors; a shift in expression "
                              "mixes the condition with the donor.")
    out["caveats"].append(
        "diffcyt's logFC for state markers is a difference in median arcsinh-transformed expression "
        "(level minus reference), not a log fold change. At high intensity arcsinh behaves like a natural log, so a "
        "difference of 1 there is roughly a 2.7-fold change in signal; near zero it is much less.")
    out["figures"] = [p for p in ((inp.ds_dir / "logfc.png") if inp.ds_dir else None,
                                  (inp.ds_dir / "heatmap.png") if inp.ds_dir else None) if p is not None and p.is_file()]
    return out


# ---------------------------------------------------------------- batch effects

LABEL_ORDER = ["raw", "uncorrected", "peacoqc", "cytonorm", "cycombine"]


def _order_labels(labels: list[str]) -> list[str]:
    return sorted(labels, key=lambda l: (LABEL_ORDER.index(l) if l in LABEL_ORDER else len(LABEL_ORDER), l))


def batch_effects(inp: Inputs) -> dict[str, Any] | None:
    if not inp.assessments and not inp.normalization:
        return None
    out: dict[str, Any] = {"observations": [], "caveats": [], "runs": [], "batches": [], "scores": []}
    labels = _order_labels(list(inp.assessments))
    out["labels"] = labels
    batches: list[str] = []
    for label in labels:
        a = inp.assessments[label]
        s = a["summary"]
        comps = s.get("comparisons") or {}
        run = {"label": label, "median_emd": s.get("median_emd"), "median_emd_source": s.get("median_emd_source"),
               "ref_median": (comps.get("reference") or {}).get("median_emd"),
               "ref_flagged": (comps.get("reference") or {}).get("n_flagged"),
               "bat_median": (comps.get("batches") or {}).get("median_emd"),
               "bat_flagged": (comps.get("batches") or {}).get("n_flagged"),
               "flagged_batches": list(s.get("flagged_batches") or []), "n_batches": s.get("n_batches"),
               "n_samples": s.get("n_samples"), "reference": s.get("reference") or {},
               "parameters": s.get("parameters") or {}}
        out["runs"].append(run)
        for b, v in (s.get("per_batch") or {}).items():
            if b not in batches:
                batches.append(b)
            out["scores"].append({"label": label, "batch": b, "score": v.get("score"),
                                  "score_source": v.get("score_source"),
                                  "mean_emd_reference": v.get("mean_emd_reference"),
                                  "mean_emd_batches": v.get("mean_emd_batches"),
                                  "n_flagged_markers": v.get("n_flagged_markers"),
                                  # markers assessed (a marker flagged only as missing from the panel is left out)
                                  "flagged_markers": [m for m in (v.get("flagged_markers") or [])
                                                      if not s.get("markers") or m in s["markers"]],
                                  "flagged": bool(v.get("flagged"))})
    by_cond: dict[str, dict[str, float]] = {}
    for label in labels:
        eb = inp.assessments[label]["emd_batches"]
        if eb is not None and {"condition", "emd"} <= set(eb.columns):
            for c, v in eb.groupby(eb["condition"].astype(str))["emd"].median().items():
                by_cond.setdefault(c, {})[label] = float(v)
    out["by_condition"] = [{"condition": c, "cells": [v.get(l) for l in labels]} for c, v in by_cond.items()]
    out["batches"] = sorted(batches)
    score = {(r["label"], r["batch"]): r for r in out["scores"]}
    out["score_rows"] = [{"batch": b, "cells": [score.get((l, b)) for l in labels]} for b in out["batches"]]
    base = next((l for l in labels if l in ("raw", "uncorrected")), None)
    corrected = [l for l in labels if l != base]
    runs = {r["label"]: r for r in out["runs"]}
    if base:
        r0 = runs[base]
        out["observations"].append(
            f"Before correction ({base}): median EMD {num(r0['median_emd'], 3)}"
            + (f" ({r0['median_emd_source']} comparison)" if r0["median_emd_source"] else "") + "; "
            + (f"flagged batches: {joined(r0['flagged_batches'])}." if r0["flagged_batches"] else "no batch flagged."))
        worst = max((score[(base, b)] for b in out["batches"] if (base, b) in score),
                    key=lambda r: r["score"] or 0, default=None)
        if worst and worst["score"] is not None:
            out["observations"].append(
                f"The most distinct batch before correction is {worst['batch']} (score {num(worst['score'], 3)}"
                + (f"; {len(worst['flagged_markers'])} flagged markers"
                   + (", including " if len(worst["flagged_markers"]) > 6 else ": ")
                   + joined(worst["flagged_markers"][:6]) if worst["flagged_markers"] else "") + ").")
    reductions = []
    for l in corrected:
        r = runs[l]
        if base and runs[base]["median_emd"] and r["median_emd"] is not None:
            red = 1 - r["median_emd"] / runs[base]["median_emd"]
            red_b = (1 - r["bat_median"] / runs[base]["bat_median"]) if (r["bat_median"] is not None and runs[base]["bat_median"]) else None
            reductions.append((l, red, red_b))
            out["observations"].append(
                f"After {l}: median EMD {num(r['median_emd'], 3)}, {red * 100:.0f}% lower than {base}"
                + (f"; batch-vs-batch median {num(runs[base]['bat_median'], 3)} → {num(r['bat_median'], 3)} "
                   f"({red_b * 100:.0f}% lower)" if red_b is not None else "") + ".")
            remaining = sorted([(b, score[(l, b)]["score"]) for b in out["batches"] if (l, b) in score
                                and score[(l, b)]["score"] is not None], key=lambda x: -x[1])
            if len(remaining) >= 2:
                b, v = remaining[0]
                others = [x[1] for x in remaining[1:]]
                if v > 1.5 * float(np.median(others)):
                    out["observations"].append(
                        f"{b} remains the most distinct batch after {l} (score {num(v, 3)}, other batches "
                        f"{num(min(others), 3)}–{num(max(others), 3)}).")
    dd0 = donor_design(inp)
    if base and len(by_cond) > 1 and corrected:
        for c, v in by_cond.items():
            if base not in v:
                continue
            best_c = min((l for l in corrected if l in v), key=lambda l: v[l], default=None)
            if best_c is None:
                continue
            same = dd0 and not dd0["batch_varies"].get(c, True)
            out["observations"].append(
                f"Batch vs batch within {c} samples" + (" (one donor in every batch)" if same else
                                                        " (donor changes with batch)" if dd0 else "") + ": median EMD "
                + num(v[base], 3) + f" ({base}) → " + ", ".join(f"{num(v[l], 3)} ({l})" for l in corrected if l in v)
                + f"; lowest after {best_c}.")
        if dd0:
            steady = [c for c in by_cond if not dd0["batch_varies"].get(c, True)]
            varies = [c for c in by_cond if dd0["batch_varies"].get(c, False)]
            if steady and varies:
                bs = {c: min((l for l in corrected if l in by_cond[c]), key=lambda l: by_cond[c][l]) for c in steady + varies}
                if any(bs[a] != bs[b] for a in steady for b in varies):
                    out["observations"].append(
                        "The corrections rank differently by condition: where only the batch differs ("
                        + joined(steady) + f") {bs[steady[0]]} came out lowest; where donor and batch change together ("
                        + joined(varies) + f") {bs[varies[0]]} did, and part of that reduction may be donor differences removed rather than batch effect.")
    out["reductions"] = [{"label": l, "reduction": red, "reduction_batches": rb} for l, red, rb in reductions]
    if len(reductions) >= 2:
        best = max(reductions, key=lambda x: x[1])
        best_b = max([x for x in reductions if x[2] is not None], key=lambda x: x[2], default=None)
        text = (f"On these numbers {best[0]} brought the batches closest together "
                f"({best[1] * 100:.0f}% lower median EMD than {base}")
        if best_b and best_b[0] == best[0]:
            text += ", and also in the batch-vs-batch comparison)."
        elif best_b:
            text += (f"); in the batch-vs-batch comparison, which does not use the reference samples alone, "
                     f"{best_b[0]} did ({best_b[2] * 100:.0f}% lower). The two comparisons rank the corrections differently.")
        else:
            text += ")."
        out["observations"].append(text)
        out["best"] = best[0]
    elif len(reductions) == 1:
        out["best"] = reductions[0][0] if reductions[0][1] > 0 else None
    # Markers still flagged after every correction.
    if corrected:
        sets = []
        for l in corrected:
            f = inp.assessments[l]["flags"]
            if f is not None and "marker" in f:
                keep = (f["level"] == "marker") & f["comparison"].isin(["reference", "batches"])
                sets.append(set(f[keep]["marker"].dropna().astype(str)))
        common = set.intersection(*sets) if sets else set()
        if common:
            out["observations"].append(
                "Markers flagged in at least one batch after every correction: " + joined(sorted(common)[:12]) + ".")
    # Caveats.
    if len(labels) > 1:
        out["caveats"].append(
            "Flags are relative to each run's own spread (median + k·MAD across batches, with a floor), so a batch "
            "that is still flagged after correction stands out among the corrected batches; it need not have a large "
            "shift. Compare the EMD values, which are on one scale in the figure below.")
    ref_vals = {json.dumps(run["reference"].get("values")) for run in out["runs"] if run["reference"]}
    if "cytonorm" in labels:
        cn = inp.normalization.get("cytonorm", {}).get("provenance") or _provenance_of(inp, "processor-cytonorm") or {}
        mode = str((cn.get("parameters") or {}).get("MODE") or "")
        if mode != "noControls":
            out["caveats"].append(
                "CytoNorm learns its correction from the reference (control) samples of each batch"
                + (f" ({joined(runs['cytonorm']['reference'].get('values') or [])})" if runs['cytonorm']['reference'].get('values') else "")
                + "; the reference comparison measures those same samples, so it shows how well the model fits its "
                "training data. The batch-vs-batch comparison (all samples within a condition) and CytoNorm's "
                "held-out replicates (below, when present) are the fairer checks.")
    if "cycombine" in labels:
        cy = inp.normalization.get("cycombine", {}).get("provenance") or _provenance_of(inp, "processor-cycombine") or {}
        if cy.get("subsampled"):
            mx = (cy.get("parameters") or {}).get("MAX_CELLS_PER_SAMPLE")
            out["caveats"].append(f"cyCombine's output is a subsample (up to {count(mx)} cells per sample), so its "
                                  "assessment covers fewer cells than the others.")
        if cy.get("covariate_used"):
            out["caveats"].append(
                f"cyCombine was given the '{(cy.get('parameters') or {}).get('COVARIATE_COLUMN')}' column as a "
                "covariate, so differences between conditions were protected from the correction.")
    dd = donor_design(inp)
    if dd and out["runs"]:
        ref_info = out["runs"][0]["reference"] or {}
        ref_col, ref_values = ref_info.get("column"), ref_info.get("values")
        if isinstance(ref_values, str):
            ref_values = [ref_values]
        varies = [c for c, v in dd["batch_varies"].items() if v]
        steady = [c for c, v in dd["batch_varies"].items() if not v]
        if steady:
            out["observations"].append(
                "The " + joined(steady) + " samples are the same donor in every batch, so differences between "
                "batches there are technical; " + ("the reference samples" + (f" ({joined(ref_values)})" if ref_values else "")
                + " are among them, which makes the reference comparison a clean measure of the batch effect."
                if ref_values else "they are the cleanest measure of the batch effect."))
        if varies:
            out["caveats"].append(
                "In the batch-vs-batch comparison the " + joined(varies) + " samples change donor with the batch ("
                + "; ".join(joined([f"{dn} in {', '.join(b)}" for dn, b in dd["per_condition"][c].items()]) for c in varies)
                + "), so part of their EMD is a donor difference that no batch correction should remove.")
    if len(ref_vals) > 1:
        out["caveats"].append("The assessments used different reference samples; their reference EMDs are not "
                              "directly comparable.")
    out["normalization"] = normalization(inp)
    out["heatmap_marker_rows"] = {}
    return out


def normalization(inp: Inputs) -> list[dict[str, Any]]:
    res = []
    for method in sorted(inp.normalization):
        e = inp.normalization[method]
        item: dict[str, Any] = {"method": method, "rows": [], "notes": []}
        ev = e.get("evaluation.csv")
        if ev is not None and "metric" in ev:
            vals = dict(zip(ev["metric"].astype(str), ev["value"]))
            item["metrics"] = vals
            if "emd_reduction" in vals:
                item["notes"].append(
                    f"cyCombine's own evaluation: EMD reduction {float(vals['emd_reduction']):.2f} (1 = all batch "
                    f"distance removed, within re-clustered cell populations), MAD score {float(vals.get('mad_score', float('nan'))):.3f} "
                    "(small = within-batch spread, i.e. biology, kept).")
        em = e.get("emd_by_marker.csv")
        if em is not None and {"marker", "emd_uncorrected", "emd_corrected"} <= set(em.columns):
            em = em.copy()
            em["reduction"] = 1 - em["emd_corrected"] / em["emd_uncorrected"].replace(0, np.nan)
            em = em.sort_values("emd_uncorrected", ascending=False)
            item["by_marker"] = [{"marker": r.marker, "before": r.emd_uncorrected, "after": r.emd_corrected,
                                  "reduction": r.reduction} for r in em.itertuples()]
            item["by_marker_columns"] = ("EMD before", "EMD after")
        rs_mc = e.get("replicate_spread_by_metacluster.csv")
        rs = e.get("replicate_spread.csv")
        spread = rs_mc if rs_mc is not None and {"group", "role", "sd_before", "sd_after"} <= set(rs_mc.columns) else rs
        if spread is not None and {"group", "role", "sd_before", "sd_after"} <= set(spread.columns):
            per_metacluster = spread is rs_mc
            item["spread_source"] = ("SD across batches of the per-batch marker median within each metacluster, "
                                     "mean over metaclusters × markers" if per_metacluster else
                                     "SD across batches of the per-batch marker quantiles, mean over markers × quantiles")
            g = spread.groupby(["group", "role"], sort=False)[["sd_before", "sd_after"]].mean().reset_index()
            g["reduction"] = 1 - g["sd_after"] / g["sd_before"]
            g["order"] = g["role"].map({"training": 0, "held-out": 1}).fillna(2)
            g = g.sort_values(["order", "group"])
            item["groups"] = [{"group": r.group, "role": r.role, "before": r.sd_before, "after": r.sd_after,
                               "reduction": r.reduction} for r in g.itertuples()]
            train = g[g["role"] == "training"]
            held = g[g["role"] == "held-out"]
            def change(r):
                return (f"{r.group} {r.reduction * 100:.0f}% lower" if r.reduction >= 0
                        else f"{r.group} {-r.reduction * 100:.0f}% higher")
            parts = []
            if len(train):
                parts.append("its training reference(s) " + joined([change(r) for r in train.itertuples()]))
            if len(held):
                parts.append("samples it never saw " + joined([change(r) for r in held.itertuples()]))
            item["notes"].append("Spread across batches after CytoNorm: " + "; ".join(parts) + ".")
            if len(train) and len(held):
                best_held = held.sort_values("reduction", ascending=False).iloc[0]
                weak = held[held["reduction"] < 0.5 * best_held["reduction"]]
                if len(weak):
                    item["notes"].append(
                        f"The held-out samples differ: {best_held['group']} improved much more than "
                        + joined(list(weak["group"])) + ". Samples that are not replicates of the training "
                        "references (another donor, another condition) are not expected to converge.")
            # Per marker, on the held-out group that improved most (the honest check).
            src = rs if rs is not None and {"group", "marker", "sd_before", "sd_after"} <= set(rs.columns) else spread
            if "marker" in src.columns:
                pick = held.sort_values("reduction", ascending=False)["group"].iloc[0] if len(held) else g["group"].iloc[0]
                pm = src[src["group"] == pick].groupby("marker")[["sd_before", "sd_after"]].mean()
                pm["reduction"] = 1 - pm["sd_after"] / pm["sd_before"]
                pm = pm.sort_values("sd_before", ascending=False)
                item["by_marker"] = [{"marker": m, "before": r.sd_before, "after": r.sd_after, "reduction": r.reduction}
                                     for m, r in pm.iterrows()]
                item["by_marker_columns"] = (f"SD across batches before ({pick})", "after")
                neg = pm[pm["reduction"] <= 0]
                if len(neg):
                    item["notes"].append(f"For {pick}, spread did not fall for " + joined(list(neg.index[:8])) + ".")
        cv = e.get("metacluster_cv.csv")
        if cv is not None and "cv" in cv:
            mx = float(pd.to_numeric(cv["cv"], errors="coerce").max())
            item["max_cv"] = mx
            if mx > 1.5:
                item["notes"].append(
                    f"The largest metacluster CV across batches is {mx:.2f}; above ~1.5–2 the CytoNorm authors read "
                    "this as the clustering being driven by the batch effect, which undermines per-cluster normalisation.")
            else:
                item["notes"].append(f"Metacluster CVs across batches are at most {mx:.2f} (below the ~1.5 warning level).")
        res.append(item)
    return res


def batch_scores_table(batch: dict[str, Any] | None) -> pd.DataFrame:
    rows = []
    for r in (batch or {}).get("scores", []):
        rows.append({**r, "flagged_markers": " ".join(r["flagged_markers"])})
    return pd.DataFrame(rows, columns=["label", "batch", "score", "score_source", "mean_emd_reference",
                                       "mean_emd_batches", "n_flagged_markers", "flagged_markers", "flagged"])


def marker_emd(inp: Inputs, label: str, comparison: str = "reference") -> pd.DataFrame | None:
    """batch x marker EMD for one assessment (reference comparison, else batches averaged over conditions)."""
    a = inp.assessments.get(label)
    if not a:
        return None
    df = a["emd_reference"] if comparison == "reference" else a["emd_batches"]
    if df is None or df.empty or not {"batch", "marker", "emd"} <= set(df.columns):
        return None
    return df.pivot_table(index="batch", columns="marker", values="emd", aggfunc="mean")


# ---------------------------------------------------------------- methods & data credit

def methods(inp: Inputs, own: dict[str, Any]) -> dict[str, Any]:
    from processor import citations as C

    groups: dict[str, dict[str, Any]] = {}
    for path, prov in inp.provenances:
        app = str(prov.get("app"))
        version = prov.get("version") or prov.get("appVersion") or ""
        tools = prov.get("tools") or prov.get("tool_versions") or {}
        if not tools and prov.get("tool"):
            tools = {prov.get("tool"): prov.get("toolVersion")}
        params = prov.get("parameters") or {}
        g = groups.setdefault(app, {"app": app, "runs": [], "step": C.APPS.get(app, {}).get("step", ""),
                                    "citations": [C.reference(k) for k in C.APPS.get(app, {}).get("cite", [])]})
        key = json.dumps([version, params], sort_keys=True, default=str)
        for run in g["runs"]:
            if run["key"] == key:
                run["files"].append(inp.rel(path))
                for k, v in tools.items():
                    run["tools"].setdefault(k, v)
                break
        else:
            g["runs"].append({"key": key, "version": version, "tools": dict(tools),
                              "parameters": params, "files": [inp.rel(path)],
                              "method": prov.get("method") if isinstance(prov.get("method"), str) else None})
    groups[own["app"]] = {"app": own["app"], "step": C.APPS[own["app"]]["step"],
                          "citations": [C.reference(k) for k in C.APPS[own["app"]]["cite"]],
                          "runs": [{"key": "", "version": own["version"], "tools": own["tools"],
                                    "parameters": own["parameters"], "files": [], "method": None}]}
    ordered = sorted(groups.values(), key=lambda g: (C.ORDER.index(g["app"]) if g["app"] in C.ORDER else 99, g["app"]))
    for g in ordered:
        for run in g["runs"]:
            run["parameters_text"] = ", ".join(f"{k}={_param(v)}" for k, v in run["parameters"].items())
            run["tools_text"] = ", ".join(f"{k} {v}" for k, v in run["tools"].items() if v)
    return {"apps": ordered}


def _param(v: Any) -> str:
    if isinstance(v, (list, dict)):
        return json.dumps(v, separators=(",", ":"))
    if v is None:
        return "null"
    return str(v)


def data_credit(inp: Inputs) -> dict[str, Any]:
    from processor import citations as C

    prep_prov = _provenance_of(inp, "processor-fcs-prep") or {}
    files: list[str] = []
    found = prep_prov.get("inputsFound") or {}
    if isinstance(found, dict):
        files = [str(f) for f in (found.get("fcs") or [])]
    sample_sheet = found.get("sampleSheet") if isinstance(found, dict) else None
    if not files and inp.samples is not None and "file" in inp.samples:
        files = [str(f) for f in inp.samples["file"]]
    if not files and inp.prep:
        files = [str(r.get("file") or r.get("input")) for r in (inp.prep.get("samples") or {}).values()]
    names = {f.rsplit("/", 1)[-1] for f in files}
    demo = C.demo_file_names()
    is_demo = bool(names) and len(names & demo) >= max(1, len(names) // 2)
    return {"files": sorted(names), "n_files": len(names), "sample_sheet": sample_sheet,
            "demo": C.DEMO_DATASET if is_demo else None, "n_demo": len(names & demo)}
