"""
The papers behind each step, keyed by app name. Every PMID was checked against
PubMed E-utilities (esummary: title, first author, year, journal, DOI); papers not
in PubMed against Crossref (journal articles) or DataCite (arXiv, Zenodo). See the
README's Citations section for how.

The methods section lists, for every app whose provenance.json is among the
inputs, its version, parameters and these citations.
"""
from __future__ import annotations

from typing import Any

P: dict[str, dict[str, str]] = {
    "flowkit": dict(authors="White S, Quinn J, Enzor J, et al.", title="FlowKit: A Python Toolkit for Integrated Manual and Automated Cytometry Analysis Workflows", journal="Front Immunol", year="2021", where="12:768541", doi="10.3389/fimmu.2021.768541", pmid="34804056"),
    "bendall": dict(authors="Bendall SC, Simonds EF, Qiu P, et al.", title="Single-cell mass cytometry of differential immune and drug responses across a human hematopoietic continuum", journal="Science", year="2011", where="332(6030):687-696", doi="10.1126/science.1198704", pmid="21551058"),
    "cycombine": dict(authors="Pedersen CB, Dam SH, Barnkob MB, et al.", title="cyCombine allows for robust integration of single-cell cytometry datasets within and across technologies", journal="Nat Commun", year="2022", where="13(1):1698", doi="10.1038/s41467-022-29383-5", pmid="35361793"),
    "nowicka": dict(authors="Nowicka M, Krieg C, Crowell HL, et al.", title="CyTOF workflow: differential discovery in high-throughput high-dimensional cytometry datasets", journal="F1000Res", year="2017", where="6:748", doi="10.12688/f1000research.11622.3", pmid="28663787"),
    "peacoqc": dict(authors="Emmaneel A, Quintelier K, Sichien D, et al.", title="PeacoQC: Peak-based selection of high quality cytometry data", journal="Cytometry A", year="2022", where="101(4):325-338", doi="10.1002/cyto.a.24501", pmid="34549881"),
    "flowcore": dict(authors="Hahne F, LeMeur N, Brinkman RR, et al.", title="flowCore: a Bioconductor package for high throughput flow cytometry", journal="BMC Bioinformatics", year="2009", where="10:106", doi="10.1186/1471-2105-10-106", pmid="19358741"),
    "cleanet": dict(authors="Ionita M, McKeague ML, Painter MM, et al.", title="Cleanet: Robust Doublet Detection in Cytometry Data Based on Protein Expression Patterns", journal="Cytometry A", year="2025", where="107(11):716-729", doi="10.1002/cyto.a.24961", pmid="41065042"),
    "flowsom": dict(authors="Van Gassen S, Callebaut B, Van Helden MJ, et al.", title="FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data", journal="Cytometry A", year="2015", where="87(7):636-645", doi="10.1002/cyto.a.22625", pmid="25573116"),
    "flowsom_protocol": dict(authors="Quintelier K, Couckuyt A, Emmaneel A, Aerts J, Saeys Y, Van Gassen S", title="Analyzing high-dimensional cytometry data using FlowSOM", journal="Nat Protoc", year="2021", where="16(8):3775-3801", doi="10.1038/s41596-021-00550-0", pmid="34172973"),
    "consensusclusterplus": dict(authors="Wilkerson MD, Hayes DN", title="ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking", journal="Bioinformatics", year="2010", where="26(12):1572-1573", doi="10.1093/bioinformatics/btq170", pmid="20427518"),
    "umap_becht": dict(authors="Becht E, McInnes L, Healy J, et al.", title="Dimensionality reduction for visualizing single-cell data using UMAP", journal="Nat Biotechnol", year="2019", where="37(1):38-44", doi="10.1038/nbt.4314", pmid="30531897"),
    "umap": dict(authors="McInnes L, Healy J, Melville J", title="UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction", journal="arXiv:1802.03426", year="2018", where="", doi="10.48550/arXiv.1802.03426", pmid=""),
    "diffcyt": dict(authors="Weber LM, Nowicka M, Soneson C, Robinson MD", title="diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering", journal="Commun Biol", year="2019", where="2:183", doi="10.1038/s42003-019-0415-5", pmid="31098416"),
    "edger": dict(authors="Robinson MD, McCarthy DJ, Smyth GK", title="edgeR: a Bioconductor package for differential expression analysis of digital gene expression data", journal="Bioinformatics", year="2010", where="26(1):139-140", doi="10.1093/bioinformatics/btp616", pmid="19910308"),
    "edger_glm": dict(authors="McCarthy DJ, Chen Y, Smyth GK", title="Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation", journal="Nucleic Acids Res", year="2012", where="40(10):4288-4297", doi="10.1093/nar/gks042", pmid="22287627"),
    "limma": dict(authors="Ritchie ME, Phipson B, Wu D, et al.", title="limma powers differential expression analyses for RNA-sequencing and microarray studies", journal="Nucleic Acids Res", year="2015", where="43(7):e47", doi="10.1093/nar/gkv007", pmid="25605792"),
    "voom": dict(authors="Law CW, Chen Y, Shi W, Smyth GK", title="voom: Precision weights unlock linear model analysis tools for RNA-seq read counts", journal="Genome Biol", year="2014", where="15(2):R29", doi="10.1186/gb-2014-15-2-r29", pmid="24485249"),
    "lme4": dict(authors="Bates D, Mächler M, Bolker B, Walker S", title="Fitting Linear Mixed-Effects Models Using lme4", journal="J Stat Softw", year="2015", where="67(1)", doi="10.18637/jss.v067.i01", pmid=""),
    "multcomp": dict(authors="Hothorn T, Bretz F, Westfall P", title="Simultaneous inference in general parametric models", journal="Biom J", year="2008", where="50(3):346-363", doi="10.1002/bimj.200810425", pmid="18481363"),
    "cytonorm": dict(authors="Van Gassen S, Gaudilliere B, Angst MS, Saeys Y, Aghaeepour N", title="CytoNorm: A Normalization Algorithm for Cytometry Data", journal="Cytometry A", year="2020", where="97(3):268-278", doi="10.1002/cyto.a.23904", pmid="31633883"),
    "cytonorm2": dict(authors="Quintelier KLA, Willemsen M, Bosteels V, Aerts JGJV, Saeys Y, Van Gassen S", title="CytoNorm 2.0: A flexible normalization framework for cytometry data without requiring dedicated controls", journal="Cytometry A", year="2025", where="107(2):69-87", doi="10.1002/cyto.a.24910", pmid="39871681"),
    "combat": dict(authors="Johnson WE, Li C, Rabinovic A", title="Adjusting batch effects in microarray expression data using empirical Bayes methods", journal="Biostatistics", year="2007", where="8(1):118-127", doi="10.1093/biostatistics/kxj037", pmid="16632515"),
    "sva": dict(authors="Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD", title="The sva package for removing batch effects and other unwanted variation in high-throughput experiments", journal="Bioinformatics", year="2012", where="28(6):882-883", doi="10.1093/bioinformatics/bts034", pmid="22257669"),
    "kohonen": dict(authors="Wehrens R, Kruisselbrink J", title="Flexible Self-Organizing Maps in kohonen 3.0", journal="J Stat Softw", year="2018", where="87(7)", doi="10.18637/jss.v087.i07", pmid=""),
    "emd_orlova": dict(authors="Orlova DY, Zimmerman N, Meehan S, et al.", title="Earth Mover's Distance (EMD): A True Metric for Comparing Biomarker Expression Levels in Cell Populations", journal="PLoS One", year="2016", where="11(3):e0151859", doi="10.1371/journal.pone.0151859", pmid="27008164"),
    "cytobatchflagr": dict(authors="Eswar S, Koenig ZT, Tursi AR, Cobeña-Reyes J, Tilburgs T, Andorf S", title="CytoBatchFlagR: A Comprehensive Framework to Objectively Assess High-Parameter Cytometry Data for Batch Effects", journal="Cytometry A", year="2026", where="109(4):241-256", doi="10.1002/cyto.a.70024", pmid="41968341"),
    "emd_rubner": dict(authors="Rubner Y, Tomasi C, Guibas LJ", title="The Earth Mover's Distance as a Metric for Image Retrieval", journal="Int J Comput Vis", year="2000", where="40(2):99-121", doi="10.1023/A:1026543900054", pmid=""),
    "scipy": dict(authors="Virtanen P, Gommers R, Oliphant TE, et al.", title="SciPy 1.0: fundamental algorithms for scientific computing in Python", journal="Nat Methods", year="2020", where="17(3):261-272", doi="10.1038/s41592-019-0686-2", pmid="32015543"),
    "liechti": dict(authors="Liechti T, Weber LM, Ashhurst TM, et al.", title="An updated guide for the perplexed: cytometry in the high-dimensional era", journal="Nat Immunol", year="2021", where="22(10):1190-1197", doi="10.1038/s41590-021-01006-z", pmid="34489590"),
    "matplotlib": dict(authors="Hunter JD", title="Matplotlib: A 2D Graphics Environment", journal="Comput Sci Eng", year="2007", where="9(3):90-95", doi="10.1109/MCSE.2007.55", pmid=""),
}

# What each step is, in a word, and the papers it rests on (most important first).
APPS: dict[str, dict[str, Any]] = {
    "processor-fcs-prep": {"step": "FCS import, panel harmonisation, arcsinh transform",
                           "cite": ["flowkit", "bendall", "cycombine", "nowicka"]},
    "processor-peacoqc": {"step": "Signal-stability QC (PeacoQC)", "cite": ["peacoqc", "flowcore"]},
    "processor-cleanet": {"step": "Doublet detection (Cleanet)", "cite": ["cleanet", "flowcore"]},
    "processor-flowsom": {"step": "Clustering (FlowSOM) and heuristic labels",
                          "cite": ["flowsom", "flowsom_protocol", "consensusclusterplus", "flowcore"]},
    "processor-cytometry-umap": {"step": "UMAP embedding", "cite": ["umap_becht", "umap"]},
    "processor-diffcyt": {"step": "Differential abundance and state (diffcyt)",
                          "cite": ["diffcyt", "nowicka", "edger", "edger_glm", "limma", "voom", "lme4", "multcomp"]},
    "processor-cytonorm": {"step": "Batch normalisation (CytoNorm)",
                           "cite": ["cytonorm", "cytonorm2", "flowsom", "flowcore"]},
    "processor-cycombine": {"step": "Batch correction (cyCombine)",
                            "cite": ["cycombine", "combat", "sva", "kohonen", "flowcore"]},
    "processor-cytometry-batch-assess": {"step": "Batch-effect assessment (EMD)",
                                         "cite": ["emd_orlova", "cycombine", "cytobatchflagr", "emd_rubner", "scipy"]},
    "processor-cytometry-report": {"step": "This report", "cite": ["nowicka", "liechti", "matplotlib"]},
}

# Pipeline order for the methods section.
ORDER = list(APPS)

DEMO_DATASET = {
    "citation": "Koenig Z, Andorf S, Tilburgs T, Eswar S. CytoBatchFlagR: Spectral flow cytometry dataset. Zenodo (2025).",
    "doi": "10.5281/zenodo.15388817",
    "licence": "CC BY 4.0",
    "description": ("Cytek Aurora spectral flow cytometry of healthy adult peripheral blood leukocytes, "
                    "pre-gated CD45+CD14- live single cells; 6 batches x 4 aliquots (2 unstimulated, "
                    "2 stimulated with anti-CD3/CD28 + IL-2 for 7 days)."),
}


def demo_file_names() -> set[str]:
    """The 24 FCS file names of the demo dataset (B01-B06 x tcs3, tcs4, tcu1, tcu2)."""
    return {f"B0{b}_{c}_CD45+CD14-.fcs" for b in range(1, 7) for c in ("tcs3", "tcs4", "tcu1", "tcu2")}


def reference(key: str) -> dict[str, str]:
    ref = dict(P[key])
    ref["key"] = key
    ref["url"] = f"https://doi.org/{ref['doi']}"
    return ref
