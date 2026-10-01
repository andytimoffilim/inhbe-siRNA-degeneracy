# inhbe-siRNA-degeneracy

Data and scripts for the cautionary note:

> **A Cautionary Note on Sequence-Only Models of RNA:**
> **Degenerate Embeddings When Base Sequences Repeat**

## Summary

Sequence-only models of RNA (e.g., RNA-FM) consume a string over
`{A, U, G, C}`. When a dataset contains multiple samples with identical
base sequences but different measured outcomes, any such model produces
identical embeddings for them, and prediction is necessarily constant
on each equivalence class.

We state this as a simple pre-screening rule, illustrate it with a small
INHBE dataset from patent WO 2026/055502 A3, and give a practical
checklist for when sequence-only models are applicable.

## The rule

Let `X` be a dataset and `f : X -> Y` a function to be predicted.
Suppose a model class `M` factors through a map that depends only on a
projection `π : X -> Σ*` onto strings over an alphabet `Σ`. For
sequence-only RNA models, `Σ = {A, U, G, C}` and `π` extracts the base
string, discarding chemical modifications. Then:

> If `x₁, x₂ ∈ X` satisfy `π(x₁) = π(x₂)` but `f(x₁) ≠ f(x₂)`,
> no model in `M` can fit `f` exactly. The prediction is constant on
> the equivalence class `{x : π(x) = π(x₁)}`.

**Practical form.** Let `N` be the number of samples and `K` the number
of distinct base sequences. If `K < N`, sequence-only predictions are
constant on each equivalence class. Furthermore, check whether `f`
varies within equivalence classes `{x : π(x) = const}`. If it does, no
sequence-only model can fit that variation, regardless of dataset size.

## Repository structure
.
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
├── data/
│ ├── NM_031479.5.fasta # INHBE, Homo sapiens
│ ├── NM_008382.3.fasta # INHBE, Mus musculus
│ ├── inhbe_golden.csv # 6 duplexes with IC50 and sites
│ ├── rna_fm_pairwise_results.csv # RNA-FM distances and similarities
│ └── inhbe_blast.tsv # BLAST hits
├── scripts/
│ ├── ocr/ # OCR pipeline for the source patent
│ ├── parsing/ # HELM/IC50 parsing, dataset build
│ ├── analysis/ # RNA-FM pairwise, BLAST, chemistry
│ └── utils/ # helper scripts
├── ocr_text/ # raw OCR for key pages
├── paper/ # manuscript source
└── results/ # output directory

text

## Data

### Source patent

**WO 2026/055502 A3** (Wave Life Sciences, PCT/US2025/045166) —
*Inhibin subunit beta E (INHBE) gene-targeting double stranded
oligonucleotides, compositions, and methods of using the same*.

The published document is a 147-page scanned PDF without a text layer.
It is **not** redistributed here. Download it from
[WIPO Patentscope](https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2026055502).

### Reference sequences

| File | Description | Length |
|------|-------------|--------|
| `data/NM_031479.5.fasta` | INHBE, *Homo sapiens* | 2460 nt |
| `data/NM_008382.3.fasta` | INHBE, *Mus musculus* | 2211 nt |

### Final dataset

`data/inhbe_golden.csv` — 6 duplexes with base sequences, binding sites,
and human-hepatocyte IC50 values.

| DSR           | Sense ID    | Antisense ID | Site      | IC50 (nM) |
|---------------|-------------|--------------|-----------|-----------|
| DSR-0104292   | SSR-0108051 | SSR-0108323  | 1107-1127 | 0.94      |
| DSR-0104290   | SSR-0108051 | SSR-0108321  | 1107-1127 | 0.74      |
| DSR-0104522   | SSR-0108051 | SSR-0108682  | 1107-1127 | 0.54      |
| DSR-0104523   | SSR-0108051 | SSR-0108697  | 1107-1127 | 0.55      |
| DSR-0104455   | SSR-0108001 | SSR-0108301  | 1107-1127 | 0.55      |
| DSR-0104212   | SSR-0108001 | SSR-0108302  | 1107-1127 | 0.90      |

**Key property:** four of six duplexes share one base sequence; two share
another. RNA-FM pairwise distances therefore take only two distinct
values per metric.

## Reproducibility

### Requirements

```bash
pip install -r requirements.txt
External tools (must be in PATH):

Tesseract OCR — for OCR

Poppler — for pdf2image

BLAST+ — for blastn, makeblastdb

Pipeline
1. OCR the source patent (requires the PDF, placed in pdf/):

bash
python scripts/ocr/ocr_find_examples.py
python scripts/ocr/ocr_high_quality.py
python scripts/ocr/ocr_expand.py
2. Parse HELM and IC50 data, build the golden dataset:

bash
python scripts/parsing/find_all_ic50.py
python scripts/parsing/parse_ic50_v2.py
python scripts/parsing/parse_and_build_golden.py
3. Run the analysis:

bash
python scripts/analysis/rna_fm_pairwise_inhbe.py
python scripts/analysis/check_embedding_variance.py
python scripts/analysis/chem_features.py
python scripts/analysis/blast_inhbe.py
Notes
Each script computes REPO_ROOT = Path(__file__).resolve().parents[2]
and resolves all internal paths relative to the repository root.

RNA-FM weights are downloaded automatically by the fm package on
first run.

The PDF of the source patent is not included; download it separately
and place it in pdf/.

Diagnostic checklist
Before applying a sequence-only model to a new dataset:

Count distinct base sequences K and samples N. If K < N, some
degeneracy is present.

For each equivalence class of identical sequences, check whether the
target variable varies within the class. If it does, no sequence-only
model can fit that variation.

If degeneracy is extensive, consider (a) chemistry-aware models that
take modification annotations as input, or (b) hybrid architectures
that concatenate sequence embeddings with explicit chemical
descriptors.

If degeneracy is not present but the model still underperforms,
consider other factors (target-site accessibility, mRNA secondary
structure, off-target effects) before concluding that the approach is
uninformative.

License
Code: MIT License (see LICENSE).

Data: CC BY 4.0.

Citation
If you use this dataset or the diagnostic rule, please cite:

text
Timofeev, A.V., Anufriev, A. (2026).
A Cautionary Note on Sequence-Only Models of RNA:
Degenerate Embeddings When Base Sequences Repeat.
Preprint.
Contact
Andrey V. Timofeev — corresponding author

LLC "AI Center for SCO+ Countries", St. Petersburg, Russia