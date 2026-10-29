# processor-cytometry-report

The last step of both immune cytometry workflows, and the only one whose output is for a
person:

- **A. Immune Profiling & Differential Abundance:** fcs-prep → peacoqc → cleanet → flowsom →
  {cytometry-umap, diffcyt} → **cytometry-report**
- **B. Multi-batch Harmonization & QC:** fcs-prep → peacoqc → {batch-assess `LABEL=raw`,
  cytonorm → batch-assess `LABEL=cytonorm`, cycombine → batch-assess `LABEL=cycombine`} →
  **cytometry-report**

It reads whatever upstream outputs it is given and writes one HTML page an immunologist can read
and check, plus summary tables. **Every input is optional**: the report detects what it has and
renders the matching sections (both workflows' inputs together render both). It works offline:
no scripts, no external resources, figures next to the page.

The text is written as **observations from the data** ("NK CD56dim CD16+ rose from 1.8% to 19%
of cells (log2FC +3.47, FDR 6e-93; in all 6 batches)") with the **caveats** that come with them
(labels being heuristic, a result carried by one batch, donors confounded with conditions, a
correction scored on its own training samples). It does not give advice.

## Sections

| Section | Shown when | What it says |
|---|---|---|
| Overview | `metadata` | samples per batch × condition (and whether they are confounded), events per sample, donors (a `donor`/`subject`/`patient` column of `samples.csv`, carried through by fcs-prep from the sample sheet), panel with clustering/state flags and cofactors, dropped channels with reasons |
| Quality control | `qc` and/or `doublets` | events in/removed per sample (PeacoQC, Cleanet), flagged files with their PeacoQC plots, batches that lost more than the others, PeacoQC notes |
| Populations | `clusters` (+ `embedding`) | labelled metaclusters with their marker profiles, a median-expression heatmap scaled per marker, composition per sample (stacked bars by condition) and per condition, UMAP figures; caveats on heuristic labels, labels that span several metaclusters, a Treg label without CD127 in the panel, doublet-like clusters |
| Differential abundance | `da` (+ `clusters`, `metadata`) | significant metaclusters with mean % per condition, log2FC and FDR; for each, how many batches go the same way and whether one batch carries most of the difference; a points-per-sample figure with batch means joined; compositional caveat |
| Differential state | `ds` | counts per comparison, top cluster × marker hits ranked by effect and by FDR, markers with the most hits, diffcyt's figures; what the "logFC" of a state marker means |
| Batch effects | `assessment` (+ `normalization`) | median EMD before/after each correction (reference and batch-vs-batch comparisons, and batch-vs-batch per condition), per-batch scores and an EMD heatmap on **one colour scale** across runs (batch-assess's own heatmaps scale per run), flagged batches and markers, which correction brought the batches closest together on which comparison; CytoNorm's training vs held-out replicate spread and metacluster CVs, cyCombine's EMD reduction and MAD score |
| Methods & citations | any `provenance.json` among the inputs | per app: version, tool versions, parameters, and the papers it rests on (table in `processor/citations.py`) |
| Input data | `metadata` | the original FCS files and sample sheet; a data credit for the demo dataset when its files are among the inputs |

### How the batch-consistency check works

For a significant metacluster, the report computes per batch the difference between the condition
means of the cluster's share of cells (level − reference). *Batches agreeing* counts the batches
whose difference has the overall sign; a batch *carries* the result when its difference is more
than half of the summed absolute differences (with 3+ batches), and the report says whether the
direction survives leaving that batch out. On the demo data this catches, e.g., MC20 (CD16+ HLA-DR+,
2.0% → 0.006%), where batch01 carries 99% of the difference.

## Ports

Inputs (all `placement: subfolder`, none required). Linked from a workflow, each arrives in
`INPUT_DIR/<port>/`; without links (everything side by side, possibly nested, symlinks followed)
the files are found by name anywhere under `INPUT_DIR`, and recognised by their neighbours
(a `summary.json` with `per_batch` is a batch-assess summary; one next to `da.csv` is diffcyt's).

| Input | Files | From |
|---|---|---|
| `metadata` | `samples.csv`, `panel.csv`, `prep.json`, `provenance.json` | fcs-prep |
| `qc` | `peacoqc.csv`, `plots/<sample_id>.png`, `provenance.json` | peacoqc (qc output) |
| `doublets` | `doublets.csv`, `plots/`, `provenance.json` | cleanet (qc output) |
| `clusters` | `counts.csv`, `labels.csv`, `medians.csv`, `star_tree.png`, `provenance.json` | flowsom (clusters output) |
| `embedding` | `umap_*.png`, `provenance.json` | cytometry-umap |
| `da` | `da.csv`, `summary.json`, `volcano.png`, `provenance.json` | diffcyt (da output) |
| `ds` | `ds.csv`, `summary.json`, `logfc.png`, `heatmap.png`, `provenance.json` | diffcyt (ds output) |
| `assessment` | `<LABEL>/summary.json`, `emd_reference.csv`, `emd_batches.csv`, `flags.csv`, `provenance.json` | cytometry-batch-assess, one link per LABEL |
| `normalization` | `cytonorm/` (`replicate_spread.csv`, `replicate_spread_by_metacluster.csv`, `metacluster_cv.csv`, `provenance.json`) and/or `cycombine/` (`evaluation.csv`, `emd_by_marker.csv`, `provenance.json`) | cytonorm / cycombine (qc outputs) |

Several links feed `assessment` and `normalization`: their folders arrive side by side. The
normalization apps write their qc into a method-named folder (`qc/cytonorm/`, `qc/cycombine/`),
so the two never collide; a single source with its files at the root of the port works too
(the method is read from its `provenance.json`, else from its file names).

Outputs (two, so each in `OUTPUT_DIR/<name>/`):

| Output | Files |
|---|---|
| `report` | `index.html`, `figures/*.png` (the report's own figures and copies of upstream ones), `provenance.json` |
| `tables` | `composition.csv`, `da_significant.csv`, `batch_scores.csv`, `qc.csv`, `metrics.json`, `provenance.json` |

`composition.csv`: `sample_id, batch, condition`, then one column per label (fraction of the
sample's cells, labels summed over their metaclusters, largest first). `da_significant.csv`:
`comparison, metacluster, label, logFC, p_adj, mean_frac_reference, mean_frac_level, reference,
level, batches_same_direction, n_batches, notes`. `batch_scores.csv`: `label, batch, score,
score_source, mean_emd_reference, mean_emd_batches, n_flagged_markers, flagged_markers, flagged`.
`qc.csv`: per sample `n_events`, PeacoQC events in/out, % removed, flagged, Cleanet doublets and
events out, % kept. `metrics.json`: every number, observation and caveat of the page, per section.

A table is written only when its section has inputs. A section whose inputs are present but
unreadable is reported on the page and in `provenance.json` (`sectionErrors`); the rest is still
rendered. The step fails (exit 1) only when none of the inputs is found; exit 2 on a bad parameter.

## Parameters

| Name | Default | What it does |
|---|---|---|
| `REPORT_TITLE` | `Immune Cytometry Report` | Title of the page. |
| `FDR` | (empty) | Threshold for what the report counts and lists as significant in DA and DS, e.g. `0.01`. Empty: diffcyt's own (its `FDR`, from `summary.json`). Display only; diffcyt's tables are not changed. |

## Results on the demo data

All upstream outputs of the demo study (24 files, 6 batches; `data/work/` of
workflow-immune-cytometry: prep, peacoqc, flowsom, umap, diffcyt with
`REFERENCE_LEVEL=unstimulated`, batch-assess raw/cytonorm/cycombine, cytonorm and cycombine qc;
no cleanet output yet). Some of what the page says:

- Donors: unstimulated tcu1/tcu2 are one donor (181C) in every batch (the sample sheet's `donor` column); the stimulated samples are
  donor 198C (batches 02, 03, 04, 06) or 199C (01, 05). So stimulation and donor coincide, and the
  DA/DS FDRs measure reproducibility across batches and aliquots, not variation between people.
- QC: PeacoQC removed a median 1.7% per sample (0–18.9%); batch03 lost more than the others.
- DA (stimulated vs unstimulated): 13 of 20 metaclusters at FDR < 0.05; NK CD56dim CD16+ 1.8% →
  19% (FDR 6e-93, all 6 batches); DN T (gd-like) 4.7% → 0.69%; MC20 CD16+ HLA-DR+ carried 99% by
  batch01.
- Batch effects: median reference EMD 0.135 raw → 0.018 CytoNorm (−86%) → 0.063 cyCombine (−53%);
  batch-vs-batch within unstimulated (one donor) 0.158 → 0.044 / 0.066, within stimulated (donor
  changes with batch) 0.140 → 0.143 / 0.040. CytoNorm cut the spread of its training reference tcu1
  by 84% and of the held-out tcu2 by 49%; batch01 stays the most distinct batch after either.

## Runtime

`cpu: 1024`, `memory: 2048`, `computeTypes: [standard, lambda]` (the image serves both; see
`processor/entry.py` and `processor/handler.py`). Measured with `docker run --memory 2g --cpus 1`
under amd64 emulation on an M-series Mac, peak from `/sys/fs/cgroup/memory.peak`:

| Run | Peak memory | Wall time (idle host) | Wall time (host busy with another build) |
|---|---|---|---|
| all inputs (both workflows) | 164-166 MB | 5.6 s | 22-30 s |
| workflow A inputs only | 158-162 MB | 4.0 s | 18-20 s |
| workflow B inputs only | 137 MB | 2.8 s | 12-15 s |

The report never reads events: only the per-sample/per-cluster tables upstream steps write, so
memory does not grow with the number of events. Image: 532 MB on disk, 123 MB compressed (`python:3.12.11-slim-bookworm`).

## Running it

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python           # synthetic fixture, both layouts, per workflow, HTML snapshot
make run PYTHON=.venv/bin/python            # Docker over data/fixture/linked/ -> data/output/
INPUT_HOST=/path/to/inputs make run         # Docker over your own folder (one subfolder per port)
```

`tests/snapshots/index.html` is the page rendered from the fixture; after an intended change to
the page, `UPDATE_SNAPSHOTS=1 make test` rewrites it (look at the diff first).

## Limitations

- Labels and caveats are only as good as the upstream outputs: the report re-derives proportions
  and batch consistency from `counts.csv`, but takes test results (logFC, FDR) from diffcyt.
- Donor confounding is detected from a donor-like column of `samples.csv` (`donor`, `subject`,
  `patient`, …), which fcs-prep carries through from the sample sheet; without one the report
  cannot see it.
- The batch comparison relies on batch-assess runs over the same samples and reference; runs with
  different references are flagged as not directly comparable.
- Figures are static PNGs (no hover); every figure has a table or a CSV beside it.
- The page is light-themed only.

## Licences

Our code: Apache-2.0 (`LICENSE`). Bundled components and their licences: `THIRD_PARTY.md`.

## Citations

The report cites, per step, the papers in `processor/citations.py` (shown in its Methods
section). Every PMID was checked with PubMed E-utilities `esummary` (title, first author, year,
journal, DOI) on 2026-10-28; papers not in PubMed against Crossref (`api.crossref.org/works/<DOI>`)
and arXiv/Zenodo DOIs against DataCite. The report's own:

- Nowicka M, Krieg C, Crowell HL, et al. CyTOF workflow: differential discovery in
  high-throughput high-dimensional cytometry datasets. *F1000Res* 6:748 (2017).
  doi:10.12688/f1000research.11622.3. PMID 28663787.
- Liechti T, Weber LM, Ashhurst TM, et al. An updated guide for the perplexed: cytometry in the
  high-dimensional era. *Nat Immunol* 22(10):1190-1197 (2021). doi:10.1038/s41590-021-01006-z.
  PMID 34489590.
- Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput Sci Eng* 9(3):90-95 (2007).
  doi:10.1109/MCSE.2007.55. (not in PubMed; checked on Crossref)

Per step (in the report when that step's provenance is among the inputs): FlowKit PMID 34804056,
Bendall 2011 PMID 21551058, cyCombine PMID 35361793 (fcs-prep); PeacoQC PMID 34549881, flowCore
PMID 19358741 (peacoqc); Cleanet PMID 41065042 (cleanet); FlowSOM PMID 25573116, FlowSOM protocol
PMID 34172973, ConsensusClusterPlus PMID 20427518 (flowsom); UMAP Becht 2019 PMID 30531897 and
McInnes 2018 doi:10.48550/arXiv.1802.03426 (DataCite) (umap); diffcyt PMID 31098416, edgeR PMIDs
19910308 and 22287627, limma PMID 25605792, voom PMID 24485249, lme4 doi:10.18637/jss.v067.i01
(Crossref), multcomp PMID 18481363 (diffcyt); CytoNorm PMID 31633883, CytoNorm 2.0 PMID 39871681
(cytonorm); ComBat PMID 16632515, sva PMID 22257669, kohonen doi:10.18637/jss.v087.i07 (Crossref)
(cycombine); EMD Orlova 2016 PMID 27008164, CytoBatchFlagR PMID 41968341, Rubner 2000
doi:10.1023/A:1026543900054 (Crossref), SciPy PMID 32015543 (batch-assess).

Demo data: Koenig Z, Andorf S, Tilburgs T, Eswar S. CytoBatchFlagR: Spectral flow cytometry
dataset. Zenodo (2025). doi:10.5281/zenodo.15388817 (CC BY 4.0; checked on DataCite).
