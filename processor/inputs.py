"""
Finding the upstream outputs.

Linked layout: each input port arrives in INPUT_DIR/<port>/ (placement: subfolder).
Flat layout (no links, an older compute node, files dropped together): everything
somewhere under INPUT_DIR, found by file name. Either way a file is recognised by
its name and by the files next to it (e.g. a summary.json beside da.csv is diffcyt's
DA summary; one with a per_batch key is a batch-assess summary). INPUT_DIR entries
are often symlinks, so the walk follows them.

Several links can feed one port (assessment, normalization): their folders arrive
side by side, e.g. assessment/raw/, assessment/cytonorm/, normalization/cytonorm/,
normalization/cycombine/. A single normalization source may also sit at the root of
its port.
"""
from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

log = logging.getLogger("processor-cytometry-report")

PORTS = ["metadata", "qc", "doublets", "clusters", "embedding", "da", "ds", "assessment", "normalization"]

# Files that identify a normalization method's qc folder.
CYTONORM_FILES = ("replicate_spread.csv", "metacluster_cv.csv", "reference_quantiles.csv")
CYCOMBINE_FILES = ("evaluation.csv", "emd_by_marker.csv")


def walk(root: Path) -> list[Path]:
    """Every file under root, following symlinked folders, in a stable order."""
    out: list[Path] = []
    seen: set[str] = set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:  # a symlink loop
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames.sort()
        for name in sorted(filenames):
            p = Path(dirpath) / name
            if p.is_file():
                out.append(p)
    return out


def read_json(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.is_file():
        return None
    try:
        with open(path) as fh:
            value = json.load(fh)
        return value if isinstance(value, dict) else None
    except (OSError, ValueError) as exc:
        log.warning("cannot read %s: %s", path, exc)
        return None


def read_csv(path: Path | None, **kwargs) -> pd.DataFrame | None:
    if path is None or not path.is_file():
        return None
    try:
        return pd.read_csv(path, **kwargs)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        log.warning("cannot read %s: %s", path, exc)
        return None


def sibling(path: Path, name: str) -> Path | None:
    p = path.parent / name
    return p if p.is_file() else None


@dataclass
class Inputs:
    input_dir: Path
    files: list[Path] = field(default_factory=list)
    # metadata
    samples: pd.DataFrame | None = None
    panel: pd.DataFrame | None = None
    prep: dict[str, Any] | None = None
    # qc (peacoqc) and doublets (cleanet)
    peacoqc: pd.DataFrame | None = None
    peacoqc_dir: Path | None = None
    doublets: pd.DataFrame | None = None
    doublets_dir: Path | None = None
    # clusters (flowsom)
    counts: pd.DataFrame | None = None
    labels: pd.DataFrame | None = None
    medians: pd.DataFrame | None = None
    clusters_dir: Path | None = None
    clusters_provenance: dict[str, Any] | None = None
    # embedding (umap)
    umap_figures: list[Path] = field(default_factory=list)
    embedding_provenance: dict[str, Any] | None = None
    # diffcyt
    da: pd.DataFrame | None = None
    da_summary: dict[str, Any] | None = None
    da_dir: Path | None = None
    ds: pd.DataFrame | None = None
    ds_summary: dict[str, Any] | None = None
    ds_dir: Path | None = None
    # batch-assess, one entry per LABEL
    assessments: dict[str, dict[str, Any]] = field(default_factory=dict)
    # cytonorm / cycombine qc, one entry per method
    normalization: dict[str, dict[str, Any]] = field(default_factory=dict)
    # every provenance.json found: (path, content)
    provenances: list[tuple[Path, dict[str, Any]]] = field(default_factory=list)
    found: dict[str, list[Path]] = field(default_factory=dict)

    def rel(self, path: Path | None) -> str | None:
        if path is None:
            return None
        try:
            return str(Path(path).relative_to(self.input_dir))
        except ValueError:
            return str(path)

    @property
    def has_any(self) -> bool:
        return any(self.found.get(p) for p in PORTS)


def _candidates(inp: Inputs, port: str, name: str) -> list[Path]:
    """Files called `name`: in the port's own folder when it exists and has one, else anywhere."""
    port_dir = inp.input_dir / port
    if port_dir.is_dir():
        inside = [f for f in inp.files if f.name == name and _under(f, port_dir)]
        if inside:
            return inside
    return [f for f in inp.files if f.name == name]


def _under(path: Path, folder: Path) -> bool:
    try:
        path.relative_to(folder)
        return True
    except ValueError:
        return False


def _has_columns(df: pd.DataFrame | None, *cols: str) -> bool:
    return df is not None and all(c in df.columns for c in cols)


def _first_table(inp: Inputs, port: str, name: str, *cols: str) -> tuple[pd.DataFrame | None, Path | None]:
    for path in _candidates(inp, port, name):
        df = read_csv(path)
        if _has_columns(df, *cols):
            return df, path
        log.warning("%s: not the expected table (needs %s); ignored", inp.rel(path), ", ".join(cols))
    return None, None


def gather(input_dir: Path) -> Inputs:
    inp = Inputs(input_dir=input_dir, files=walk(input_dir))
    found = {p: [] for p in PORTS}

    # metadata (fcs-prep)
    inp.samples, p = _first_table(inp, "metadata", "samples.csv", "sample_id")
    if p:
        inp.samples["sample_id"] = inp.samples["sample_id"].astype(str)
        for col in ("batch", "condition"):
            if col in inp.samples.columns:
                inp.samples[col] = inp.samples[col].astype(str)
        found["metadata"].append(p)
    inp.panel, p = _first_table(inp, "metadata", "panel.csv", "channel", "marker")
    if p:
        found["metadata"].append(p)
    for path in _candidates(inp, "metadata", "prep.json"):
        inp.prep = read_json(path)
        if inp.prep is not None:
            found["metadata"].append(path)
            break

    # qc (peacoqc)
    inp.peacoqc, p = _first_table(inp, "qc", "peacoqc.csv", "sample_id")
    if p:
        inp.peacoqc["sample_id"] = inp.peacoqc["sample_id"].astype(str)
        inp.peacoqc_dir = p.parent
        found["qc"].append(p)

    # doublets (cleanet)
    inp.doublets, p = _first_table(inp, "doublets", "doublets.csv", "sample_id")
    if p:
        inp.doublets["sample_id"] = inp.doublets["sample_id"].astype(str)
        inp.doublets_dir = p.parent
        found["doublets"].append(p)

    # clusters (flowsom)
    inp.counts, p = _first_table(inp, "clusters", "counts.csv", "sample_id", "metacluster", "n")
    if p:
        inp.counts["sample_id"] = inp.counts["sample_id"].astype(str)
        inp.clusters_dir = p.parent
        found["clusters"].append(p)
        inp.clusters_provenance = read_json(sibling(p, "provenance.json"))
    inp.labels, p2 = _first_table(inp, "clusters", "labels.csv", "metacluster", "label")
    if p2:
        found["clusters"].append(p2)
        inp.clusters_dir = inp.clusters_dir or p2.parent
    inp.medians, p3 = _first_table(inp, "clusters", "medians.csv", "metacluster")
    if p3:
        found["clusters"].append(p3)
        inp.clusters_dir = inp.clusters_dir or p3.parent
    if inp.clusters_dir is not None:
        for name in ("heatmap.png", "star_tree.png"):
            if (inp.clusters_dir / name).is_file():
                found["clusters"].append(inp.clusters_dir / name)

    # embedding (umap): its figures; the embedding itself is not needed
    umap_dirs = {f.parent for f in _candidates(inp, "embedding", "embedding.csv.gz")}
    port_dir = inp.input_dir / "embedding"
    figs = [f for f in inp.files if f.name.startswith("umap_") and f.suffix == ".png"
            and (f.parent in umap_dirs or (port_dir.is_dir() and _under(f, port_dir)))]
    if not figs and not port_dir.is_dir():
        figs = [f for f in inp.files if f.name.startswith("umap_") and f.suffix == ".png"]
    inp.umap_figures = figs
    found["embedding"] = sorted({*figs, *(d / "embedding.csv.gz" for d in umap_dirs)})
    for d in sorted(umap_dirs | {f.parent for f in figs}):
        inp.embedding_provenance = inp.embedding_provenance or read_json(d / "provenance.json")

    # diffcyt da / ds
    inp.da, p = _first_table(inp, "da", "da.csv", "metacluster", "logFC", "p_adj")
    if p:
        inp.da_dir = p.parent
        inp.da_summary = read_json(sibling(p, "summary.json"))
        found["da"].append(p)
    inp.ds, p = _first_table(inp, "ds", "ds.csv", "metacluster", "marker", "logFC", "p_adj")
    if p:
        inp.ds_dir = p.parent
        inp.ds_summary = read_json(sibling(p, "summary.json"))
        found["ds"].append(p)

    # assessment (batch-assess, one folder per LABEL)
    for path in _candidates(inp, "assessment", "summary.json"):
        summary = read_json(path)
        if not summary or "per_batch" not in summary:
            continue
        label = str(summary.get("label") or path.parent.name)
        if label in inp.assessments:
            log.warning("assessment %s found twice (%s and %s); keeping the first", label,
                        inp.rel(inp.assessments[label]["dir"]), inp.rel(path.parent))
            continue
        d = path.parent
        inp.assessments[label] = {
            "label": label, "dir": d, "summary": summary,
            "emd_reference": read_csv(sibling(path, "emd_reference.csv")),
            "emd_batches": read_csv(sibling(path, "emd_batches.csv")),
            "flags": read_csv(sibling(path, "flags.csv")),
            "provenance": read_json(sibling(path, "provenance.json")),
        }
        found["assessment"].append(path)

    # normalization (cytonorm / cycombine qc): method-named subfolders, or at the root
    norm_dirs: dict[Path, str] = {}
    port_dir = inp.input_dir / "normalization"
    for f in inp.files:
        if port_dir.is_dir() and not _under(f, port_dir):
            continue
        if f.name in CYTONORM_FILES:
            norm_dirs.setdefault(f.parent, "cytonorm")
        elif f.name in CYCOMBINE_FILES:
            norm_dirs.setdefault(f.parent, "cycombine")
    for d, guess in sorted(norm_dirs.items()):
        prov = read_json(d / "provenance.json")
        app = str((prov or {}).get("app") or "")
        method = "cytonorm" if "cytonorm" in app else "cycombine" if "cycombine" in app else guess
        if method in inp.normalization:
            log.warning("normalization %s found twice (%s and %s); keeping the first", method,
                        inp.rel(inp.normalization[method]["dir"]), inp.rel(d))
            continue
        entry: dict[str, Any] = {"method": method, "dir": d, "provenance": prov}
        for name in CYTONORM_FILES + CYCOMBINE_FILES + ("emd.csv", "mad.csv", "replicate_spread_by_metacluster.csv"):
            if (d / name).is_file():
                entry[name] = read_csv(d / name)
                found["normalization"].append(d / name)
        inp.normalization[method] = entry

    # every provenance.json (for methods & citations)
    for f in inp.files:
        if f.name == "provenance.json":
            prov = read_json(f)
            if prov and prov.get("app"):
                inp.provenances.append((f, prov))

    inp.found = found
    return inp
