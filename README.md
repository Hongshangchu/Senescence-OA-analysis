# Senescence programs in OA: main pipeline, GWAS-locus ORA and an
# independent single-cell validation (GSE220243)

Analysis code accompanying:

> **Single-cell and bulk transcriptomics reveal convergent
> senescence-associated remodeling in osteoarthritic cartilage** (in
> preparation; manuscript formatted for *Arthritis Research & Therapy*)

This repository is the GitHub-ready, three-script refactor of the two
original analysis scripts cited in the manuscript's data-availability
statement — `run_all11.R (pipeline v7.6)` and
`run_all12_GWAS_senescence_ora.R` — plus the new independent-cohort
validation of the SenMayo enrichment in FC/EffC chondrocytes
(GSE220243).

## Why three scripts instead of one?

The two original scripts are **integrated but deliberately not merged**:

- `01_pipeline_main.R` (≈1,150 lines) is the full discovery pipeline: it
  needs several GB of GEO downloads (GSE169454, GSE114007, GSE104782),
  20+ Bioconductor/CRAN packages and hours of runtime, and its headline
  SenMayo statistic is a **GSEA** of the ranked disease transcriptome.
- `02_ORA_GWAS_senescence.R` and `03_validation_GSE220243.R` are
  self-contained: all reference gene sets ship in `data/ref/`, and the
  validation reads a user-downloaded public object. They run in seconds.

Merging them would force anyone reproducing the ORA/validation to first
install the full heavy environment and download gigabytes of data, so
the three share the same repository conventions instead: one
`scripts/00_utils.R` helper file (logger, `save_fig` with PDF + 600-dpi
TIFF output, download-with-retry, SenMayo loading), one `output/`
layout (`output/tables`, `output/figures`, `output/logs`) and one
`data/README.md` input specification.

| Script | Analysis | Output | Runtime / data |
| --- | --- | --- | --- |
| `scripts/01_pipeline_main.R` | Full discovery pipeline (modules 1–15: preprocessing + Harmony, pseudobulk DESeq2, bulk age × disease interaction, GSE104782 replication, GO/GSEA enrichment, donor-level SenMayo scoring, cis-MR, CellChat, Monocle3 pseudotime, coloc) | `output/tables/Table1–6, TableS1–S12_*.csv`; `output/figures/Fig1–10, FigS0–S13_*.{pdf,tiff}` | Heavy: GEO downloads + Bioconductor env; hours |
| `scripts/02_ORA_GWAS_senescence.R` | SenMayo / GO:0090398 over-representation among Boer 2021 OA GWAS loci; cross-cohort direction-consistency test; optional MAGMA input, proteomics check and summary figure | `output/tables/TableS13–S16_*.csv`, `output/figures/FigS10_*.pdf` | Base R only; seconds |
| `scripts/03_validation_GSE220243.R` | Independent-cohort validation: SenMayo module enrichment in FC-like / EffC-like chondrocytes in GSE220243 | `output/GSE220243/*.csv`, `output/figures/GSE220243/FigS11–S12_*.pdf` | Seurat (+ optional harmony/ggplot2); minutes |

## Repository layout

```
senescence-OA-analysis/
├── scripts/
│   ├── 00_utils.R                     # shared helpers (logger, save_fig,
│   │                                  #  dl_with_retry, load_senmayo, ora_test…)
│   ├── 01_pipeline_main.R             # full discovery pipeline (run_all11 v7.6)
│   ├── 02_ORA_GWAS_senescence.R       # GWAS-locus ORA (run_all12)
│   └── 03_validation_GSE220243.R      # independent-cohort validation
├── data/
│   ├── README.md                      # input specifications
│   ├── ref/                           # committed reference gene sets (see below)
│   └── GSE220243/                     # user-downloaded dataset (gitignored)
├── output/                            # tables + figures + logs (gitignored,
│                                      #  created at runtime)
└── .gitignore
```

## Requirements

**Script 01** (the heavy one):

- R >= 4.2 (cairo backend for TIFF output), Bioconductor: `GEOquery`,
  `DESeq2`, `limma`, `ComplexHeatmap`, `clusterProfiler`,
  `org.Hs.eg.db`, `Biobase`, `monocle3`, `CellChat` (>= 1.6), `coloc`,
  `TxDb.Hsapiens.UCSC.hg38.knownGene`.
- CRAN: `Seurat` (>= v5), `harmony`, `tidyverse`, `data.table`,
  `readxl`, `patchwork`, `EnhancedVolcano`, `enrichplot`, `pROC`,
  `TwoSampleMR`, `MRInstruments`, `ieugwasr`.
- OpenGWAS JWT: set `OPENGWAS_JWT` in `~/.Renviron` (MR clumping);
  never commit the token.

**Script 02**: base R only (optional `ggplot2` for the figure).

**Script 03**: `Seurat` (>= v4; v5 recommended), optional `harmony`
(batch correction by donor) and `ggplot2` (figures).

## Quick start

```sh
# 0) shared helpers are sourced automatically by each script

# 1) Full discovery pipeline (needs GEO downloads; see data/README.md)
Rscript scripts/01_pipeline_main.R

# 2) GWAS-locus ORA (no external data; all reference sets shipped)
Rscript scripts/02_ORA_GWAS_senescence.R

# 3) Independent-cohort validation (download GSE220243 first)
Rscript scripts/03_validation_GSE220243.R
```

All three scripts locate the repository root automatically, so they can
be run from the repo root or from `scripts/`. All outputs carry a
timestamped log under `output/logs/`.

## Reference data (`data/ref/`, committed)

| Set | Source | Size |
| --- | --- | --- |
| SenMayo | Saul et al., *Nat Commun* 2022 (125 human genes) | 125 |
| GO:0090398 cellular senescence | QuickGO, human, NOT-qualified excluded | 99 |
| Boer 2021 OA GWAS lead-locus genes (all / protein-coding) | Boer et al., *Cell* 2021 (103 lead SNPs); gene mapping by GWAS Catalog / Ensembl genomic context | 659 / 298 |
| ORA background | Ensembl BioMart protein-coding genes with an HGNC symbol | 19,473 |

## Two different SenMayo analyses — do not conflate them

- **GSEA of the ranked disease transcriptome** (main pipeline,
  module 8 → `TableS7`, `Fig7`): the strong, manuscript-headline
  signal — **NES = 2.64, adjusted P = 2.99 × 10⁻¹³**, 34 leading-edge
  genes.
- **ORA of Boer 2021 OA GWAS loci** (script 02 → `TableS13`): a
  supporting, weak-signal analysis. With the protein-coding background
  (n = 19,473) and the 298 protein-coding locus genes, SenMayo shows
  only a marginal overlap (2/124 genes; OR ≈ 1.06, one-sided Fisher
  P = 0.57, BH adj. P = 0.57) while GO:0090398 is nominally enriched
  (4/83; OR ≈ 3.29, P = 0.039, BH adj. P = 0.077).

The weak GWAS-locus ORA is consistent with the manuscript's central
claim that aging–OA linkage operates at the level of convergent
programs rather than gene-level overlap.

## Validation cohort (GSE220243)

GSE220243 is an independent human single-cell cohort deposited by
**Kim Y, Jung G, Han S. *CARTILAGE* 2026,
doi:10.1177/19476035261481078** (72,616 chondrocytes), who identified
superficial-zone chondrocytes as a catabolic and senescent focus in OA
cartilage. The GEO series contains **12 cartilage samples
(`CartNorm1`–`CartNorm6`, `CartOA1`–`CartOA6`)** plus meniscus samples
(`Men*`, kept under `data/GSE220243/meniscus_raw/` and not used here);
the validation uses only the cartilage samples. Because
the public object may not carry the authors' cell-type labels, the
validation uses **annotation-free primary tests** that mirror the
discovery-cohort hypothesis:

- cluster-level Spearman correlation between mean SenMayo and mean
  FC-module (`COL1A1`, `COL1A2`) / EffC-module (`ADAMTS5`, `CASP3`)
  scores;
- cell-level one-sided Wilcoxon enrichment of SenMayo in top-tercile
  FC-like / EffC-like cells (descriptive: cells are not independent
  within a donor);
- if annotations are present, the label-based test and a donor-level
  OA vs normal SenMayo contrast (median fold change + AUC) are added.
  When the object lacks a condition column but sample names encode the
  disease (e.g. `CartNorm*` / `CartOA*`), the script derives it
  automatically so the disease contrast runs with 6 vs 6 donors.

## Notes

- The MAGMA module (script 02, module 4) requires genome-wide summary
  statistics that Boer 2021 did not deposit in GWAS Catalog; download
  instructions are in the script header and `data/README.md`.
- Numbers in the manuscript supplementary tables (Table S13–S16,
  Fig S10) are produced by script 02; the headline SenMayo statistics
  (NES = 2.64, adjusted P = 2.99 × 10⁻¹³) come from script 01's GSEA
  module. Keep `data/ref/` untouched for reproducibility.
- The refactor preserves all module parameters, the fixed RNG seed
  (`set.seed(123)`), output file names and the per-figure PDF + 600-dpi
  TIFF dual output of the original scripts; results remain directly
  comparable to the manuscript.

## License

All rights reserved. Contact the corresponding author for reuse.
