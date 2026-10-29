"""
A small synthetic set of upstream outputs, following the workflow contract, for
tests and `make run`: 3 batches x 2 conditions x 2 samples, 5 metaclusters, 6
markers, and every input port (metadata, qc, doublets, clusters, embedding, da,
ds, assessment raw/cytonorm/cycombine, normalization cytonorm/cycombine).

  python -m tests.make_fixture data/fixture --layout both [--workflow all|a|b]

Layouts: linked (INPUT_DIR/<port>/...) and flat (everything under arbitrary,
nested folder names, found by file name).
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

BATCHES = ["b1", "b2", "b3"]
CONDITIONS = ["unstim", "stim"]
MARKERS = [("CD3", True, False), ("CD4", True, False), ("CD8", True, False), ("CD56", True, False),
           ("CD25", False, True), ("FOXP3", True, False), ("GZMB", False, True)]
MCS = [(1, "CD4 T", 0.40), (2, "CD8 T", 0.25), (3, "NK", 0.10), (4, "CD4 T CD25+FOXP3+ (Treg or activated)", 0.20),
       (5, "T/B doublet-like", 0.05)]
TINY_PNG = bytes.fromhex(
    "89504e470d0a1a0a0000000d49484452000000010000000108060000001f15c489"
    "0000000d49444154789c6360000002000154a24f5d0000000049454e44ae426082")

WORKFLOW_A = ["metadata", "qc", "doublets", "clusters", "embedding", "da", "ds"]
WORKFLOW_B = ["metadata", "qc", "assessment", "normalization"]


def _prov(app: str, **extra) -> dict:
    return {"app": app, "version": "0.1.0", "tools": {"R": "4.5.0"}, "parameters": {"SEED": 42}, **extra}


def _write(path: Path, content) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(content, pd.DataFrame):
        content.to_csv(path, index=False)
    elif isinstance(content, (dict, list)):
        path.write_text(json.dumps(content, indent=1))
    elif isinstance(content, bytes):
        path.write_bytes(content)
    else:
        path.write_text(str(content))


def port_files() -> dict[str, dict[str, object]]:
    """{port: {relative path inside the port: content}}"""
    rng = np.random.default_rng(1)
    samples = []
    for b in BATCHES:
        for c in CONDITIONS:
            for k in (1, 2):
                sid = f"{b}_{c}{k}"
                samples.append({"sample_id": sid, "file": f"{sid}.fcs", "batch": b, "condition": c,
                                "n_events": 30000 if sid != "b1_unstim1" else 4000, "control": f"{c}{k}",
                                "donor": "D1" if c == "unstim" else ("D3" if b == "b3" else "D2")})
    samples = pd.DataFrame(samples)
    panel = pd.DataFrame([{"channel": "FSC-A", "marker": "FSC-A", "type": "scatter", "cofactor": "",
                           "clustering": "false", "state": "false"}] +
                         [{"channel": f"CH{i}-A", "marker": m, "type": "marker", "cofactor": 6000,
                           "clustering": str(cl).lower(), "state": str(st).lower()} for i, (m, cl, st) in enumerate(MARKERS)])
    prep = {"app": "processor-fcs-prep", "cytometerType": "spectral", "channelsDropped": ["CH9-A"],
            "samples": {r.sample_id: {"file": r.file, "cytometer": "Aurora",
                                      "channels_dropped": [{"channel": "CH9-A", "marker": "CD28",
                                                            "reason": "missing in 4 of 12 samples"}]}
                        for r in samples.itertuples()}}
    prep_prov = _prov("processor-fcs-prep", inputsFound={"fcs": list(samples["file"]), "sampleSheet": "samples.csv"})

    pq = pd.DataFrame({"sample_id": samples["sample_id"], "events_in": samples["n_events"],
                       "pct_removed": [2.0, 1.5, 3.0, 2.5, 12.0, 14.0, 11.0, 35.0, 1.0, 2.0, 1.5, 0.5],
                       "status": "ok", "notes": ""})
    pq["events_out"] = (pq["events_in"] * (1 - pq["pct_removed"] / 100)).round().astype(int)
    pq["flagged"] = np.where(pq["pct_removed"] > 30, "true", "false")
    pq.loc[0, "notes"] = "There are not enough bins for a robust isolation tree analysis."
    db = pd.DataFrame({"sample_id": samples["sample_id"], "events_in": pq["events_out"],
                       "pct_doublets": np.round(rng.uniform(0.5, 3, len(samples)), 3), "status": "ok", "notes": ""})
    db["doublets"] = (db["events_in"] * db["pct_doublets"] / 100).round().astype(int)
    db["events_out"] = db["events_in"] - db["doublets"]

    counts = []
    for r in samples.itertuples():
        n = int(r.n_events)
        fr = np.array([f for _, _, f in MCS])
        if r.condition == "stim":
            fr = fr * np.array([0.8, 1.0, 3.0, 1.0, 1.0])            # NK up in every batch
            if r.batch == "b2":
                fr = fr * np.array([1.0, 3.0, 1.0, 1.0, 1.0])        # CD8 up only in b2
        fr = fr / fr.sum()
        for (mc, label, _), f in zip(MCS, fr):
            counts.append({"sample_id": r.sample_id, "metacluster": mc, "label": label, "n": int(n * f), "frac": f})
    counts = pd.DataFrame(counts)
    tot = counts.groupby("metacluster")["n"].sum()
    labels = pd.DataFrame([{"metacluster": mc, "label": label, "n_cells": int(tot[mc]), "frac": tot[mc] / tot.sum(),
                            "rule": mc, "profile": "CD3+ CD4+"} for mc, label, _ in MCS])
    medians = pd.DataFrame([{"metacluster": mc, "label": label, "n_cells": int(tot[mc]),
                             **{m: round(float(rng.uniform(0, 3)), 3) for m, _, _ in MARKERS}} for mc, label, _ in MCS])
    flowsom_prov = _prov("processor-flowsom", output="clusters",
                         clustering_markers=[m for m, c, _ in MARKERS if c], state_markers=[m for m, _, s in MARKERS if s],
                         method={"labels": "heuristic", "rules": [{"label": "CD4 T"}]})

    # diffcyt-like tables computed from the counts
    comp = counts.merge(samples[["sample_id", "condition"]], on="sample_id")
    comp["f"] = comp["n"] / comp.groupby("sample_id")["n"].transform("sum")
    mean = comp.groupby(["metacluster", "condition"])["f"].mean().unstack()
    da = pd.DataFrame({"metacluster": mean.index, "label": [l for _, l, _ in MCS], "comparison": "stim vs unstim",
                       "logFC": np.log2(mean["stim"] / mean["unstim"]).round(4).values, "logCPM": 15.0,
                       "p_val": [0.02, 0.003, 1e-12, 0.6, 0.9], "p_adj": [0.03, 0.006, 5e-12, 0.7, 0.9],
                       "tested": "true", "mean_frac_unstim": mean["unstim"].values, "mean_frac_stim": mean["stim"].values})
    da["significant"] = np.where(da["p_adj"] < 0.05, "true", "false")
    da_summary = {"status": "ok", "method": "edgeR", "fdr": 0.05, "reference": "unstim", "block": "batch",
                  "block_effect": "fixed", "logFC": "edgeR model coefficient (log2 fold change)"}
    ds_rows = []
    for mc, label, _ in MCS:
        for m in ("CD25", "GZMB"):
            d = 1.2 if (mc == 4 and m == "CD25") else 0.1 * mc
            p = 1e-6 if d > 1 else (0.01 if mc < 3 else 0.5)
            ds_rows.append({"metacluster": mc, "label": label, "marker": m, "comparison": "stim vs unstim",
                            "logFC": d, "p_val": p, "p_adj": p * 2, "tested": "true",
                            "significant": "true" if p * 2 < 0.05 else "false",
                            "mean_median_unstim": 0.5, "mean_median_stim": 0.5 + d})
    ds = pd.DataFrame(ds_rows)
    ds_summary = {"status": "ok", "method": "limma", "fdr": 0.05, "reference": "unstim"}
    diffcyt_prov = _prov("processor-diffcyt")

    def assessment(label: str, scale: float) -> dict[str, object]:
        mk = [m for m, _, _ in MARKERS]
        ref = pd.DataFrame([{"batch": b, "marker": m, "emd": round(scale * (0.3 if b == "b1" else 0.1) * (1 + i / 10), 4),
                             "reference": "unstim1", "threshold": 0.2 * scale, "flagged": "false", "n_samples": 1,
                             "n_cells": 1000} for b in BATCHES for i, m in enumerate(mk)])
        ref["flagged"] = np.where(ref["emd"] > ref["threshold"], "true", "false")
        bat = pd.DataFrame([{"batch": b, "condition": c, "marker": m, "emd": round(scale * 0.15, 4), "threshold": 1,
                             "flagged": "false", "n_samples": 2, "n_cells": 2000}
                            for b in BATCHES for c in CONDITIONS for m in mk])
        flags = ref[ref["flagged"] == "true"].assign(level="marker", comparison="reference", condition="",
                                                     reason="reference sample differs")
        flags = pd.concat([flags, pd.DataFrame([{"level": "marker", "comparison": "missing", "batch": "b3",
                                                 "marker": "CD28", "reason": "marker only in 2 of 3 batches"}])])
        flags = flags[["level", "comparison", "batch", "condition", "reference", "marker", "emd", "threshold", "reason"]]
        per = {}
        for b in BATCHES:
            sub = ref[ref["batch"] == b]
            fm = list(sub[sub["flagged"] == "true"]["marker"])
            per[b] = {"score": float(sub["emd"].mean()), "score_source": "reference",
                      "mean_emd_reference": float(sub["emd"].mean()), "mean_emd_batches": scale * 0.15,
                      "n_flagged_markers": len(fm), "flagged_markers": fm, "flagged": len(fm) > 1}
        summary = {"label": label, "n_batches": 3, "n_markers": len(mk), "n_samples": 12,
                   "median_emd": float(ref["emd"].median()), "median_emd_source": "reference",
                   "comparisons": {"reference": {"n": len(ref), "median_emd": float(ref["emd"].median()),
                                                 "n_flagged": int((ref["flagged"] == "true").sum())},
                                   "batches": {"n": len(bat), "median_emd": float(bat["emd"].median()), "n_flagged": 0}},
                   "per_batch": per, "flagged_batches": [b for b, v in per.items() if v["flagged"]],
                   "reference": {"column": "control", "values": ["unstim1"], "n_samples": 3},
                   "parameters": {"LABEL": label}}
        return {f"{label}/summary.json": summary, f"{label}/emd_reference.csv": ref, f"{label}/emd_batches.csv": bat,
                f"{label}/flags.csv": flags, f"{label}/heatmap_reference.png": TINY_PNG,
                f"{label}/provenance.json": _prov("processor-cytometry-batch-assess", parameters={"LABEL": label})}

    spread = pd.DataFrame([{"group": g, "role": role, "n_batches": 3, "marker": m, "channel": "x", "quantile": q,
                            "sd_before": 0.2, "sd_after": 0.2 * (0.4 if role == "training" else 0.7), "reduction": 0}
                           for g, role in (("unstim1", "training"), ("unstim2", "held-out"))
                           for m, _, _ in MARKERS for q in (0.05, 0.5, 0.95)])
    spread["reduction"] = 1 - spread["sd_after"] / spread["sd_before"]
    normalization = {
        "cytonorm/replicate_spread.csv": spread,
        "cytonorm/metacluster_cv.csv": pd.DataFrame({"metacluster": [1, 2], "cv": [0.3, 0.8], "pct_b1": [50, 50]}),
        "cytonorm/provenance.json": _prov("processor-cytonorm", parameters={"MODE": "controls", "SEED": 42}),
        "cycombine/evaluation.csv": pd.DataFrame({"metric": ["emd_reduction", "mad_score"], "value": [0.7, 0.02]}),
        "cycombine/emd_by_marker.csv": pd.DataFrame({"marker": [m for m, _, _ in MARKERS], "emd_uncorrected": 2.0,
                                                     "emd_corrected": 0.6}),
        "cycombine/provenance.json": _prov("processor-cycombine", tool_versions={"cyCombine": "0.3.0"}, subsampled=True,
                                           covariate_used=True,
                                           parameters={"COVARIATE_COLUMN": "condition", "MAX_CELLS_PER_SAMPLE": 50000}),
    }
    emb = gzip.compress(b"sample_id,metacluster,label,umap1,umap2\nb1_unstim1,1,CD4 T,0.1,0.2\n", mtime=0)
    return {
        "metadata": {"samples.csv": samples, "panel.csv": panel, "prep.json": prep, "provenance.json": prep_prov},
        "qc": {"peacoqc.csv": pq, "plots/b2_stim2.png": TINY_PNG,
               "provenance.json": _prov("processor-peacoqc", parameters={"MAX_REMOVED_PCT": 30})},
        "doublets": {"doublets.csv": db, "plots/b1_unstim1.png": TINY_PNG, "provenance.json": _prov("processor-cleanet")},
        "clusters": {"counts.csv": counts, "labels.csv": labels, "medians.csv": medians, "heatmap.png": TINY_PNG,
                     "star_tree.png": TINY_PNG, "provenance.json": flowsom_prov},
        "embedding": {"embedding.csv.gz": emb, "umap_by_label.png": TINY_PNG, "umap_by_condition.png": TINY_PNG,
                      "provenance.json": _prov("processor-cytometry-umap", cells_embedded=1200,
                                               parameters={"N_NEIGHBORS": 15, "MIN_DIST": 0.3, "MARKERS_used": ["CD3"]})},
        "da": {"da.csv": da, "summary.json": da_summary, "volcano.png": TINY_PNG, "provenance.json": diffcyt_prov},
        "ds": {"ds.csv": ds, "summary.json": ds_summary, "logfc.png": TINY_PNG, "provenance.json": diffcyt_prov},
        "assessment": {**assessment("raw", 1.0), **assessment("cytonorm", 0.3), **assessment("cycombine", 0.5)},
        "normalization": normalization,
    }


FLAT_FOLDERS = {"metadata": "prep/metadata", "qc": "peacoqc/qc", "doublets": "cleanet/qc", "clusters": "flowsom/clusters",
                "embedding": "umap", "da": "diffcyt/da", "ds": "diffcyt/ds", "assessment": "batch-assess",
                "normalization": "norm-qc"}


def make_fixture(root: Path, layout: str = "both", workflow: str = "all") -> dict[str, Path]:
    root = Path(root)
    ports = {"all": WORKFLOW_A + WORKFLOW_B, "a": WORKFLOW_A, "b": WORKFLOW_B}[workflow]
    files = port_files()
    out = {}
    for lay in (["linked", "flat"] if layout == "both" else [layout]):
        base = root / lay
        if base.exists():
            shutil.rmtree(base)
        for port in dict.fromkeys(ports):
            folder = base / (port if lay == "linked" else FLAT_FOLDERS[port])
            for rel, content in files[port].items():
                _write(folder / rel, content)
        out[lay] = base
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--layout", default="both", choices=["both", "linked", "flat"])
    ap.add_argument("--workflow", default="all", choices=["all", "a", "b"])
    a = ap.parse_args()
    for k, v in make_fixture(Path(a.root), a.layout, a.workflow).items():
        print(k, v)
