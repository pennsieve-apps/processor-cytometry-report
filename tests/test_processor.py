"""End-to-end: the synthetic fixture through the processor, in both input layouts, per workflow, plus the
Lambda handler and an HTML snapshot."""
from __future__ import annotations

import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pytest

from processor import handler, inputs
from processor import main as app
from tests.make_fixture import make_fixture

SNAPSHOT = Path(__file__).parent / "snapshots" / "index.html"
FIXED_TIME = datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc)


@pytest.fixture(scope="session")
def fixtures(tmp_path_factory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("fixture")
    out = make_fixture(root / "all", "both", "all")
    out["a"] = make_fixture(root / "a", "linked", "a")["linked"]
    out["b"] = make_fixture(root / "b", "linked", "b")["linked"]
    return out


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for name, value in params.items():
        monkeypatch.setenv(name, value)
    return app.main()


def page(out: Path) -> str:
    return (out / "report" / "index.html").read_text()


def assert_self_contained(out: Path) -> None:
    html = page(out)
    assert "<script" not in html.lower()
    assert not re.search(r'(src|href)="(?!https://doi\.org/|#|figures/)', html), "only local figures and DOI links"
    assert "<link" not in html and "@import" not in html
    for src in re.findall(r'<img src="([^"]+)"', html):
        assert (out / "report" / src).is_file(), src


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_all_inputs(fixtures, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fixtures[layout], out, monkeypatch) == 0
    html = page(out)
    for anchor in ("overview", "qc", "populations", "da", "ds", "batch", "methods", "data"):
        assert f'id="{anchor}"' in html, anchor
    assert_self_contained(out)

    comp = pd.read_csv(out / "tables" / "composition.csv")
    assert len(comp) == 12 and {"sample_id", "batch", "condition", "NK", "CD4 T"} <= set(comp.columns)
    labels = [c for c in comp.columns if c not in ("sample_id", "batch", "condition")]
    assert comp[labels].sum(axis=1).round(6).eq(1).all()

    da = pd.read_csv(out / "tables" / "da_significant.csv")
    assert list(da["metacluster"]) == [3, 2, 1]  # by FDR
    cd8 = da[da["metacluster"] == 2].iloc[0]
    assert "b2 carries" in cd8["notes"]  # the CD8 rise is only in batch b2
    assert da[da["metacluster"] == 3].iloc[0]["batches_same_direction"] == 3

    scores = pd.read_csv(out / "tables" / "batch_scores.csv")
    assert set(scores["label"]) == {"raw", "cytonorm", "cycombine"} and len(scores) == 9
    m = json.loads((out / "tables" / "metrics.json").read_text())
    assert m["batch_effects"]["labels"] == ["raw", "cytonorm", "cycombine"]
    assert m["batch_effects"]["best"] == "cytonorm"
    assert {n["method"] for n in m["batch_effects"]["normalization"]} == {"cytonorm", "cycombine"}
    assert m["inputData"]["n_files"] == 12 and m["inputData"]["demo"] is None
    qc = pd.read_csv(out / "tables" / "qc.csv")
    assert qc.set_index("sample_id").loc["b2_stim2", "peacoqc_flagged"]

    study = " ".join(m["study"]["caveats"])
    assert "No donor appears in more than one condition" in study and "donor changes with the batch" in study
    assert any("different donors" in c for c in m["differential_abundance"]["caveats"])
    assert any("change donor with the batch" in c for c in m["batch_effects"]["caveats"])
    obs = " ".join(m["populations"]["caveats"])
    assert "no CD127" in obs and "doublet" in obs
    assert any("NK rose from" in o for o in m["differential_abundance"]["observations"])
    assert any("b2 (median" in o for o in m["qc"]["observations"])  # a batch losing more to QC
    # Methods: every app's provenance, with citations; this report's own entry last.
    apps = [a["app"] for a in m["methods"]]
    assert apps[0] == "processor-fcs-prep" and apps[-1] == "processor-cytometry-report"
    assert "processor-cytometry-batch-assess" in apps and "PMID 31098416" in html
    prov = json.loads((out / "report" / "provenance.json").read_text())
    assert set(prov["assessmentLabels"]) == {"raw", "cytonorm", "cycombine"}
    assert all(prov["sections"].values()) and not prov["sectionErrors"]


def test_workflow_a_only(fixtures, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixtures["a"], out, monkeypatch) == 0
    html = page(out)
    assert 'id="da"' in html and 'id="populations"' in html and 'id="batch"' not in html
    assert not (out / "tables" / "batch_scores.csv").exists()
    assert "Immune Profiling" in html and "Harmonization" not in html


def test_workflow_b_only(fixtures, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixtures["b"], out, monkeypatch) == 0
    html = page(out)
    assert 'id="batch"' in html and 'id="da"' not in html and 'id="populations"' not in html
    assert (out / "tables" / "batch_scores.csv").exists() and not (out / "tables" / "composition.csv").exists()
    assert "held-out replicates" in html and "EMD reduction 0.70" in html
    assert_self_contained(out)


def test_single_normalization_at_port_root(fixtures, tmp_path):
    src = tmp_path / "in"
    shutil.copytree(fixtures["b"], src)
    shutil.rmtree(src / "normalization" / "cytonorm")
    for f in (src / "normalization" / "cycombine").iterdir():
        f.rename(src / "normalization" / f.name)
    found = inputs.gather(src)
    assert list(found.normalization) == ["cycombine"]
    assert found.normalization["cycombine"]["dir"] == src / "normalization"


def test_fdr_parameter_and_title(fixtures, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixtures["linked"], out, monkeypatch, FDR="0.01", REPORT_TITLE="My Study") == 0
    da = pd.read_csv(out / "tables" / "da_significant.csv")
    assert list(da["metacluster"]) == [3, 2]
    assert "<title>My Study</title>" in page(out)


def test_bad_fdr(fixtures, tmp_path, monkeypatch):
    assert run(fixtures["linked"], tmp_path / "out", monkeypatch, FDR="2") == 2


def test_no_inputs(tmp_path, monkeypatch):
    (tmp_path / "in" / "other").mkdir(parents=True)
    (tmp_path / "in" / "other" / "readme.txt").write_text("nothing here")
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


def test_symlinked_inputs(fixtures, tmp_path, monkeypatch):
    src = tmp_path / "in"
    src.mkdir()
    for port in fixtures["linked"].iterdir():
        os.symlink(port, src / port.name)
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    assert 'id="batch"' in page(out) and 'id="da"' in page(out)


def test_handler(fixtures, tmp_path, monkeypatch):
    out = tmp_path / "out"
    monkeypatch.setenv("FDR", "0.5")  # left over from a warm start: must be cleared
    result = handler.handler({"inputDir": str(fixtures["linked"]), "outputDir": str(out),
                              "REPORT_TITLE": "Lambda run"}, None)
    assert result["status"] == "ok" and "FDR" not in os.environ
    assert "<title>Lambda run</title>" in page(out)
    with pytest.raises(RuntimeError):
        handler.handler({"inputDir": str(tmp_path / "missing"), "outputDir": str(out)}, None)


def test_html_snapshot(fixtures, tmp_path, monkeypatch):
    """The page for the fixture, byte for byte. After an intended change: UPDATE_SNAPSHOTS=1 make test."""
    monkeypatch.setattr(app, "_now", lambda: FIXED_TIME)
    monkeypatch.setattr(app, "tool_versions", lambda: {"python": "3.12", "pandas": "x", "numpy": "x",
                                                        "matplotlib": "x", "jinja2": "x"})
    out = tmp_path / "out"
    assert run(fixtures["linked"], out, monkeypatch) == 0
    html = page(out)
    if os.environ.get("UPDATE_SNAPSHOTS") or not SNAPSHOT.exists():
        SNAPSHOT.parent.mkdir(parents=True, exist_ok=True)
        SNAPSHOT.write_text(html)
    assert html == SNAPSHOT.read_text(), "index.html changed; inspect the diff, then UPDATE_SNAPSHOTS=1 make test"
