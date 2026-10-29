"""
The report's own figures (matplotlib, PNG). One look across them: the light chart
surface, thin marks, recessive axes, a fixed categorical order (never cycled), a
single-hue sequential ramp for magnitudes. Every figure returns its file name, or
None when it has nothing to draw.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
OTHER = "#c3c2b7"
CRITICAL = "#d03b3b"
# Categorical slots, in order (validated: adjacent pairs pass CVD and normal-vision separation).
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
# Sequential blue, light -> dark.
SEQ = LinearSegmentedColormap.from_list("seq", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DPI = 110

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.titlecolor": INK, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
})


def _save(fig, out_dir: Path, name: str) -> str:
    fig.savefig(out_dir / name, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return name


def _short(text: str, n: int = 38) -> str:
    return text if len(text) <= n else text[: n - 1] + "…"


# ---------------------------------------------------------------- QC

def qc_bars(rows: list[dict[str, Any]], out_dir: Path) -> str | None:
    have_qc = any(r.get("qc_pct") is not None and not pd.isna(r.get("qc_pct")) for r in rows)
    have_db = any(r.get("db_pct") is not None and not pd.isna(r.get("db_pct")) for r in rows)
    if not (have_qc or have_db):
        return None
    panels = [(k, t) for k, t, have in (("qc_pct", "Removed by PeacoQC (%)", have_qc),
                                         ("db_pct", "Called doublets by Cleanet (%)", have_db)) if have]
    n = len(rows)
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels) + 1.2, 0.22 * n + 1.0), sharey=True,
                             squeeze=False)
    y = np.arange(n)
    for ax, (key, title) in zip(axes[0], panels):
        vals = np.array([float(r.get(key)) if r.get(key) is not None and not pd.isna(r.get(key)) else np.nan
                         for r in rows])
        flagged = [bool(r.get("qc_flagged")) and key == "qc_pct" for r in rows]
        colors = [CRITICAL if f else SLOTS[0] for f in flagged]
        ax.barh(y, np.nan_to_num(vals), color=colors, height=0.7, edgecolor=SURFACE, linewidth=1)
        for yi, v, f in zip(y, vals, flagged):
            if not np.isnan(v):
                ax.text(v, yi, f" {v:.1f}" + (" ▲ flagged" if f else ""), va="center", fontsize=7,
                        color=CRITICAL if f else INK2)
        ax.set_title(title, loc="left")
        ax.xaxis.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_xlim(0, (np.nanmax(vals) if np.isfinite(np.nanmax(vals)) else 1) * 1.35 + 0.1)
    axes[0][0].set_yticks(y, [r["sample_id"] for r in rows], fontsize=7)
    axes[0][0].invert_yaxis()
    return _save(fig, out_dir, "qc_removed.png")


# ---------------------------------------------------------------- populations

def marker_heatmap(medians: pd.DataFrame | None, clustering: list[str], state: list[str], out_dir: Path) -> str | None:
    if medians is None or medians.empty:
        return None
    skip = {"metacluster", "label", "n_cells", "frac"}
    markers = [c for c in medians.columns if c not in skip]
    order = [m for m in clustering if m in markers] + [m for m in state if m in markers and m not in clustering]
    order += [m for m in markers if m not in order]
    if not order:
        return None
    m = medians[order].apply(pd.to_numeric, errors="coerce")
    lo, hi = m.min(), m.max()
    scaled = (m - lo) / (hi - lo).replace(0, np.nan)
    rows = [f"MC{r.metacluster} " + _short(str(getattr(r, "label", "")), 34) for r in medians.itertuples()]
    fig, ax = plt.subplots(figsize=(0.32 * len(order) + 3.2, 0.26 * len(rows) + 1.6))
    ax.imshow(scaled.fillna(0).values, aspect="auto", cmap=SEQ, vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(order)), order, rotation=90, fontsize=8)
    ax.set_yticks(range(len(rows)), rows, fontsize=8)
    for i, mk in enumerate(order):
        if mk in clustering:
            ax.get_xticklabels()[i].set_fontweight("bold")
    n_cl = len([x for x in order if x in clustering])
    if 0 < n_cl < len(order):
        ax.axvline(n_cl - 0.5, color=INK, linewidth=1.2)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("Median expression per metacluster (scaled 0–1 per marker)", loc="left")
    sm = plt.cm.ScalarMappable(cmap=SEQ, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01)
    cb.set_ticks([0, 1], labels=["lowest", "highest"])
    cb.outline.set_visible(False)
    return _save(fig, out_dir, "marker_heatmap.png")


def composition_bars(comp: pd.DataFrame | None, out_dir: Path, top: int = 7) -> str | None:
    if comp is None or comp.empty:
        return None
    labels = [c for c in comp.columns if c not in ("sample_id", "batch", "condition")]
    shown = labels[:top]
    rest = labels[top:]
    df = comp.sort_values(["condition", "batch", "sample_id"]).reset_index(drop=True)
    data = df[shown].copy()
    if rest:
        data["Other"] = df[rest].sum(axis=1)
    colors = SLOTS[: len(shown)] + ([OTHER] if rest else [])
    n = len(df)
    fig, ax = plt.subplots(figsize=(max(6.0, 0.3 * n + 3.5), 4.2))
    x = np.arange(n)
    bottom = np.zeros(n)
    for col, color in zip(data.columns, colors):
        vals = data[col].values * 100
        ax.bar(x, vals, bottom=bottom, color=color, width=0.82, edgecolor=SURFACE, linewidth=1.2,
               label=_short(col, 40))
        bottom += vals
    ax.set_xticks(x, df["sample_id"], rotation=90, fontsize=7)
    ax.set_ylabel("% of cells")
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.6, n - 0.4)
    # Condition groups under the axis.
    conds = df["condition"].tolist()
    start = 0
    for i in range(1, n + 1):
        if i == n or conds[i] != conds[start]:
            if df["condition"].nunique() > 1:
                ax.annotate(conds[start], xy=((start + i - 1) / 2, 101), ha="center", va="bottom", fontsize=8,
                            color=INK, fontweight="bold", annotation_clip=False)
                if i < n:
                    ax.axvline(i - 0.5, color=INK, linewidth=1)
            start = i
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=7.5, reverse=True)
    ax.set_title("Composition per sample" + (", grouped by condition" if df["condition"].nunique() > 1 else ""),
                 loc="left", pad=16)
    return _save(fig, out_dir, "composition.png")


# ---------------------------------------------------------------- differential abundance

def da_points(counts: pd.DataFrame | None, samples: pd.DataFrame | None, comparison: dict[str, Any], out_dir: Path,
              name: str, max_panels: int = 12) -> str | None:
    rows = comparison["significant"][:max_panels]
    if counts is None or samples is None or not rows or "condition" not in samples:
        return None
    level, ref = comparison["level"], comparison["reference"]
    tot = counts.groupby("sample_id")["n"].sum()
    meta = samples.set_index("sample_id")
    ncol = min(4, len(rows))
    nrow = math.ceil(len(rows) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.3 * ncol + 0.4, 2.35 * nrow + 0.5), squeeze=False)
    rng = np.random.default_rng(0)
    for ax, r in zip(axes.flat, rows):
        sub = counts[counts["metacluster"].astype(int) == r["metacluster"]].groupby("sample_id")["n"].sum()
        frac = (sub.reindex(tot.index, fill_value=0) / tot * 100).rename("pct").to_frame()
        frac = frac.join(meta[[c for c in ("batch", "condition") if c in meta]], how="left")
        frac = frac[frac["condition"].isin([ref, level])]
        xs = {ref: 0, level: 1}
        if "batch" in frac:
            bm = frac.groupby(["batch", "condition"])["pct"].mean().unstack()
            for b, v in bm.iterrows():
                if ref in v and level in v and not (pd.isna(v[ref]) or pd.isna(v[level])):
                    ax.plot([0, 1], [v[ref], v[level]], color=OTHER, linewidth=0.9, zorder=1)
        for cond, color in ((ref, SLOTS[0]), (level, SLOTS[1])):
            vals = frac[frac["condition"] == cond]["pct"].values
            jitter = rng.uniform(-0.09, 0.09, len(vals))
            ax.scatter(xs[cond] + jitter, vals, s=22, color=color, edgecolor=SURFACE, linewidth=0.8, zorder=3)
            if len(vals):
                ax.hlines(np.mean(vals), xs[cond] - 0.22, xs[cond] + 0.22, color=INK, linewidth=1.4, zorder=4)
        ax.set_xticks([0, 1], [_short(ref, 14), _short(level, 14)], fontsize=7.5)
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylim(bottom=0)
        ax.yaxis.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_title(f"MC{r['metacluster']} {_short(r['label'], 24)}\nlog2FC {r['logFC']:+.2f}, FDR {r['p_adj']:.0e}"
                     .replace("e-0", "e-"), fontsize=8, loc="left")
        ax.tick_params(labelsize=7)
    for ax in list(axes.flat)[len(rows):]:
        ax.set_visible(False)
    axes[0][0].set_ylabel("% of cells")
    fig.suptitle("Significant metaclusters: % of cells per sample (bar = mean; grey lines join batch means)",
                 x=0.01, ha="left", fontsize=9, color=INK2)
    fig.tight_layout()
    return _save(fig, out_dir, name)


# ---------------------------------------------------------------- batch effects

def batch_scores(batch: dict[str, Any], out_dir: Path) -> str | None:
    labels = batch.get("labels") or []
    rows = batch.get("score_rows") or []
    if not labels or not rows:
        return None
    nb, nl = len(rows), len(labels)
    fig, ax = plt.subplots(figsize=(max(5.5, 0.9 * nb * max(1, nl / 2) + 2), 3.4))
    width = 0.8 / nl
    x = np.arange(nb)
    for j, label in enumerate(labels):
        vals = [(r["cells"][j] or {}).get("score") for r in rows]
        vals = np.array([np.nan if v is None else float(v) for v in vals])
        flagged = [bool((r["cells"][j] or {}).get("flagged")) for r in rows]
        pos = x - 0.4 + width * (j + 0.5)
        ax.bar(pos, np.nan_to_num(vals), width=width, color=SLOTS[j % len(SLOTS)], edgecolor=SURFACE,
               linewidth=1.2, label=label)
        for p, v, f in zip(pos, vals, flagged):
            if f and not np.isnan(v):
                ax.text(p, v, "▲", ha="center", va="bottom", fontsize=7, color=CRITICAL)
    ax.set_xticks(x, [r["batch"] for r in rows])
    ax.set_ylabel("batch score (mean EMD)")
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(ncols=nl, loc="upper right", fontsize=8)
    ax.set_title("Per-batch score before and after correction (one scale; ▲ = flagged in that run)", loc="left")
    return _save(fig, out_dir, "batch_scores.png")


def emd_heatmaps(mats: dict[str, pd.DataFrame], flags: dict[str, set[tuple[str, str]]], out_dir: Path,
                 title: str) -> str | None:
    mats = {k: v for k, v in mats.items() if v is not None and not v.empty}
    if not mats:
        return None
    markers: list[str] = []
    batches: list[str] = []
    for m in mats.values():
        markers += [c for c in m.columns if c not in markers]
        batches += [b for b in m.index if b not in batches]
    batches = sorted(batches)
    first = next(iter(mats.values()))
    order = list(first.mean().sort_values(ascending=False).index)
    markers = order + [m for m in markers if m not in order]
    vmax = max(float(np.nanmax(m.values)) for m in mats.values())
    n = len(mats)
    fig, axes = plt.subplots(n, 1, figsize=(0.3 * len(markers) + 2.4, (0.25 * len(batches) + 0.8) * n + 0.6),
                             squeeze=False)
    for ax, (label, m) in zip(axes[:, 0], mats.items()):
        m = m.reindex(index=batches, columns=markers)
        ax.imshow(m.values, aspect="auto", cmap=SEQ, vmin=0, vmax=vmax, interpolation="nearest")
        for (b, mk) in flags.get(label, set()):
            if b in batches and mk in markers:
                i, j = batches.index(b), markers.index(mk)
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=CRITICAL, linewidth=1.2))
        ax.set_yticks(range(len(batches)), batches, fontsize=8)
        ax.set_xticks(range(len(markers)), markers if ax is axes[-1, 0] else [""] * len(markers), rotation=90,
                      fontsize=8)
        ax.set_title(label, loc="left", fontsize=9)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
    sm = plt.cm.ScalarMappable(cmap=SEQ, norm=plt.Normalize(0, vmax))
    cb = fig.colorbar(sm, ax=axes[:, 0].tolist(), fraction=0.02, pad=0.01)
    cb.set_label("EMD (one scale for all runs)")
    cb.outline.set_visible(False)
    fig.suptitle(title, x=0.01, ha="left", fontsize=10, fontweight="bold", color=INK)
    return _save(fig, out_dir, "batch_emd_heatmaps.png")
