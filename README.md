# scDARI: discrete Dynamical Enhancer-driven gene Regulatory network Inference from paired scRNA-seq and scATAC-seq data

## Project Overview

scDARI is an R package that reconstructs **activating enhancer-driven gene regulatory networks (eGRNs)** — ordered transcription factor (TF)–enhancer–target gene (TG) cascades — from paired single-cell RNA sequencing (scRNA-seq) and single-cell assay for transposase-accessible chromatin sequencing (scATAC-seq) data.

The workflow proceeds as follows: cells are ordered along pseudotime and adaptively aggregated into metacells until every retained element jointly meets a non-zero-observation floor across the three analysis series (TF expression, TFE; TG activity, TGA; TG expression, TGE); JASPAR2024 motif scanning and ChIP-Atlas priors define a cell-type-specific initial GRN (iGRN); three lagged discrete time-evolution equations, grounded in Hill-type activation and homeostatic negative feedback, are fitted per branch by a hybrid genetic algorithm (GA) and gradient ascent; branch results are integrated through a four-level optimal-prediction selection; and final TF–peak–TG triples pass an eleven-feature-adjusted inverse-quantile filter with two connectivity rescue mechanisms.

Key features:

- Paired-modality integration that links chromatin accessibility to gene expression at single-cell resolution
- Pseudotime ordering with adaptive metacell aggregation that mitigates single-cell sparsity while preserving branch resolution
- Initial networks from JASPAR2024 position weight matrices (PWMs) refined by ChIP-Atlas ChIP-seq priors
- A discrete dynamical model with 23 trainable parameters per TF–TG pair, fitted by hybrid GA global search plus gradient-ascent local refinement, with multiple early-stopping criteria
- Four-level hierarchical selection of the optimal prediction per pair and cell-count-weighted integration across branches
- An eleven-feature-adjusted inverse-quantile peGRN filter with two connectivity rescue mechanisms
- Two confidence levels of the final network: `loose` (comprehensive) and `strict` (conserved across all branches)

## Installation

### Prerequisites

- **R 4.4**: scDARI is developed and tested under R 4.4.
- **Operating system**: Linux and WSL2 are recommended because parallel execution (`mclapply`) is supported; on Windows, set `ncores = 1` throughout (parallel execution is unavailable).
- **igraph 2.0.3**: the Monocle2 trajectory step requires this specific version; later igraph releases are incompatible.

### Installing dependencies

```r
# igraph 2.0.3 must be installed from the CRAN archive
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz",
  repos = NULL, type = "source"
)

# CRAN packages
install.packages(c("dplyr", "ggplot2", "Signac", "SeuratObject", "Seurat", "VGAM",
                   "ggrepel", "minpack.lm", "data.table", "tibble", "tidyr", "pROC",
                   "DDRTree", "GA", "psych", "qs", "devtools", "rbioapi", "peakRAM"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomeInfoDb", "TFBSTools", "JASPAR2024", "IRanges", "ChIPseeker",
                       "Biobase", "BiocGenerics", "monocle", "Biostrings", "S4Vectors",
                       "DOSE", "org.Hs.eg.db"))

# Genome annotation packages — choose the set matching your reference build
# hg38 / GRCh38:
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "EnsDb.Hsapiens.v86",
                       "BSgenome.Hsapiens.UCSC.hg38"))
# hg19 / GRCh37 (optional):
# BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene",
#                        "EnsDb.Hsapiens.v75",
#                        "BSgenome.Hsapiens.UCSC.hg19"))
```

### Installing scDARI

```r
devtools::install_github("hth-buaa/scDARI")
library(scDARI)
packageVersion("scDARI")  # 2.0
```

## Quick Start

### Data requirements

1. **Cell-matched scRNA-seq and scATAC-seq Seurat objects** with identical cell barcodes. For data from separate experiments, establish one-to-one cell correspondence and unify identifiers first. The RNA assay must be a Seurat v5 `Assay5` object (convert v4 `Assay` objects with `SeuratObject::CreateAssay5Object`); the ATAC assay is a Signac `ChromatinAssay`.
2. **Cell type labels**: a `seurat_annotations` character column in the metadata of both objects.
3. **Fragment files**: a `_fragments.tsv.gz` file with its `.tbi` tabix index (see Part 0 if you need to generate them from the count matrix).

### Workflow at a glance

| Stage | Main functions | Output |
|---|---|---|
| Part 0: data preparation | manual / Signac | matched Seurat objects, indexed fragments |
| Part 1: preprocessing | `get_TFs_from_priors`, `identify_TGs`, `filter_peak_annotations`, `get_expression_and_activity_matrix` | TF list, refined peak–TG annotation, TFE/TGA/TGE matrices |
| Part 2: pseudotime & grouping | `get_interest_cell_type_data`, `order_pseudotime_and_divide_branches`, `get_genes_pseudotime_info`, `cell_grouping` | branches, pseudotime series, adaptive metacells |
| Part 3: iGRN | `get_iGRN_by_TFBS_pwm_by_JASPAR2024`, `filter_iGRN_by_ChIP_Atlas_priors`, `get_cell_type_tau`, `get_branch_iGRN` | global and branch iGRNs, threshold tau |
| Part 4: model training | `select_training_set`, `set_init_params`, `model_train` | trained per-pair models (23 parameters each) |
| Part 5: GRN & peGRN | `get_data_for_inferring_GRN`, `infer_GRN`, `select_best_prediction`, `filter_GRN`, `infer_eGRN`, `filter_eGRN` | pGRN and final loose/strict peGRN |

The fully executable walkthrough (with every argument, expected output and runtime log) is provided in `scDARI_guide_line.R` / `scDARI_guide_line.txt`. It uses the public `pbmcMultiome` dataset (monocyte lineage) as the running example.

## Detailed Workflow

### Part 0: Data Loading and Preprocessing

- **Example data**: the tutorial installs `pbmcMultiome` via SeuratData (or a local install of `pbmcMultiome.SeuratData_0.1.4.tar.gz`); the raw objects contain 36,601 genes / 108,377 peaks across 11,909 matched cells.
- **Assay upgrade**: convert RNA `Assay` (v4) to `Assay5` (v5); the ATAC assay remains a `ChromatinAssay`.
- **Canonical chromosomes**: remove peaks on non-standard contigs with `GenomeInfoDb::standardChromosomes`.
- **Cell type selection**: subset both modalities to the cell types of interest; when a trajectory spans several related labels (e.g., CD14 Mono and CD16 Mono), store the originals in `original_seurat_annotations` and set a unified label (`"Mono"`) in `seurat_annotations`.
- **Optional downsampling**: cells can be ranked by the combined non-zero frequency of RNA and ATAC and truncated (the tutorial keeps the top 3,000); remove non-coding (`AC`/`AL`/`LINC`-prefixed) and mitochondrial genes and keep the top 10,000 variable genes (VST); TF-IDF-normalize peaks and keep the top 10,000 most variable. Excessive cell downsampling reduces metacell count and gene coverage.
- **Fragment files**: download the official `_fragments.tsv.gz` when available; otherwise export the ATAC count matrix in triplet form to an uncompressed BED-like TSV, then `bgzip -i file.tsv` and `tabix -p bed file.tsv.gz` (Linux).

### Part 1: Data Preprocessing

#### 1.1 Identifying TFs — `get_TFs_from_priors`

Integrates candidate TFs from JASPAR2024 PWMs (all ten collections by default) and ChIP-Atlas target-gene priors for the chosen `genome_assembly` (e.g., `"hg38"`), applies six expression-quality criteria (absolute and genome-wide percentile thresholds for non-zero frequency, non-zero mean, and overall mean), and validates candidates against ChIP-seq priors. Key arguments: `distance_from_TSS` (1/5/10 kb), `prior_min_TG_number`, `prior_top_n_keep`, `prior_top_q_keep`, and `excluded_cell_types_in_priors` (set to matched cell types to exclude self-priors and prevent circular evaluation). Prior files are downloaded per TF with a status table (`download_information.rds`); rerun the function until every TF downloads successfully.

#### 1.2 Annotating peaks and identifying TGs — `identify_TGs` then `filter_peak_annotations`

- `identify_TGs` annotates peaks with ChIPseeker using a `TxDb` object and an `org.*.eg.db` annotation package. `fragment_annotation_range` defines the TSS-centered search region (default 500,000 bp; broadly defined enhancers include promoter-annotated peaks), and `flankDistance` (default 10,000 bp) controls flanking-gene lookup. Only genes detected in the scRNA-seq object are retained.
- `filter_peak_annotations` refines peak–TG links with an inverse-distance quantile curve (`distance_points` / `min_quantile_points`, tutorial anchors (1, 0.35), (500, 0.55), (500,000, 1.55)), gene-expression quality weights (default 0.7/0.2/0.1), peak-accessibility scales, and `one_peak_at_least_per_gene = TRUE` as a per-gene fallback. Repeated runs use disk caching.

#### 1.3 Expression and activity matrices — `get_expression_and_activity_matrix`

Builds the TF/TG expression matrices from RNA and the TG activity matrix with `Signac::GeneActivity` (`flanking_region`, default 10,000 bp, requiring the indexed fragment file), returning aligned `scRNAseq_for_GRN` / `scATACseq_for_GRN` objects and filtered TF/TG lists.

### Part 2: Pseudotime Inference, Branching, and Cell Grouping

#### 2.1 Cell-type extraction and trajectory — `get_interest_cell_type_data`, `order_pseudotime_and_divide_branches`

`get_interest_cell_type_data` aligns both modalities for each target cell type and applies zero-rate QC. `order_pseudotime_and_divide_branches` fits a Monocle2 DDRTree trajectory (`auto_param_selection = TRUE` by default; set it to `FALSE` to tune DDRTree parameters manually, e.g., `ncenter`) and interactively requests, from the exported trajectory plots, the number of branches and the integer state sequence of each branch. Pseudotime is a computational proxy for temporal ordering, not physical time.

#### 2.2 Pseudotime series and adaptive grouping — `get_genes_pseudotime_info`, `cell_grouping`

- `get_genes_pseudotime_info` builds per-branch TFE/TGA/TGE matrices ordered by pseudotime. Gene filtering uses importance-score quantiles `alpha_TF` / `alpha_TG` (set to 0 to disable importance filtering), cell filtering uses `alpha_cell`, and optional per-series non-zero-ratio floors are supplied through `sparsity_cut`.
- `cell_grouping` expands contiguous pseudotime windows until every retained gene accumulates the required number of non-zero observations in all three series. The minimum-count target is set by fitting `y = c·ln(ax + b)` through three anchor points (`points_x_for_fitting_nc_and_nmin`, `points_y_for_fitting_nc_and_nmin`); anchors can be derived from the mean branch cell count. **At least 8 metacells per branch are recommended** (5–7 are acceptable thanks to regularization); avoid more than ~20 groups. The function returns aggregated non-zero means and pseudotime endpoints for dynamic modeling.

### Part 3: Constructing Initial GRNs (iGRNs)

1. **PWM scanning** — `get_iGRN_by_TFBS_pwm_by_JASPAR2024` scans accessible peaks with JASPAR2024 PWMs (`min_score_for_matchPWM = "0%"` keeps all candidates; stringent filtering happens later) and returns a TG × TF matrix of normalized scores (0 to ~0.95).
2. **Prior refinement** — `filter_iGRN_by_ChIP_Atlas_priors` applies prior-weighted correction, a dual global/per-TF threshold (calibrated by `expected_peak_number_per_TF` and the dataset's peak-to-TG ratio), an absolute score floor, and target-gene rescue, then re-normalizes while preserving zeros; all-zero rows/columns are dropped.
3. **Binarization threshold** — `get_cell_type_tau` computes the cell-type-specific threshold tau; edges above tau are regulatory.
4. **Branch iGRNs** — `get_branch_iGRN` restricts the global iGRN to each branch's retained TFs and TGs.

### Part 4: Training the Predictive Models (per branch)

1. **Training set** — `select_training_set` keeps the strongest partners per regulon (`top_n_for_max_theta_in_each_regulon`, typically 1), optionally applies mutual-top refinement (`refine`, `keep_all_TF`), and appends the grouped TFE/TGA/TGE features. Training pairs serve only as kinetic templates; inference later covers every candidate pair.
2. **Parameter initialization** — `set_init_params` sets initial values and bounds for all **23 trainable parameters** per pair, with baseline expression initialized adaptively from training-set means.
3. **Training** — `model_train` runs, per epoch, a GA global search (`popSize_ga`, `maxiter_ga`, `pcrossover_ga`, `pmutation_ga`) followed by gradient-ascent refinement (`iterations_grad`), with early stopping after `early_stop` stagnant epochs and convergence tolerances `eps_theta` / `eps_loss`; `max_epochs` bounds the total. Low-quality training pairs can be removed afterward (the tutorial drops pairs with `theta_i < 0.011`). Stop causes and per-pair epoch/time statistics are returned for diagnostics. Regulons are trained in parallel on Linux.

### Part 5: GRN Prediction and peGRN Construction

1. **Inference data** — `get_data_for_inferring_GRN` attaches grouped time-series features to every candidate pair per branch (and a global deduplicated pair list).
2. **Prediction** — `infer_GRN` predicts `theta_p` and fitting loss for all pairs, either in one pass or branch-by-branch (the latter avoids accumulated memory use on large datasets).
3. **Four-level optimal selection** — `select_best_prediction` chooses, per TF–TG pair, the best of the branch's trained models under dual quantile thresholds, then integrates branches by cell-count-weighted averaging.
4. **pGRN filtering** — `filter_GRN` applies dual thresholds: the loss cutoff is `min`(loss quantile, `fit_loss_upper`) and the strength cutoff is `max`(strength quantile, `theta_p_lower`, typically tau). The **loose pGRN** retains pairs present in any branch; the **strict pGRN** only those present in every branch.
5. **peGRN assembly** — `infer_eGRN` joins PWM matches, the pGRN, and peak–TG annotations into TF–peak–TG triples annotated with 25 multi-omics features (TF–peak, TF–TG, and per-component peak/TG/TF features with normalized ranks).
6. **peGRN filtering** — `filter_eGRN` fits an inverse quantile curve (threshold decreases as predicted strength rises; tutorial anchors (0, 5.79), (0.5, 4.66), (1, 4.23)), adjusts it per edge with eleven feature-specific factors (one penalty and ten relaxations, with `factor_balance` / `factor_range`), and applies two sequential rescues under fixed quality floors: peak–TG rescue (`rescue_peakTG`, at least `min_TF_per_peakTG` TF per peak–TG pair) and TF rescue (`rescue_TF`, at least `min_peakTG_per_TF` pairs per TF). Outputs include the full `loose_filtered` / `strict_filtered` tables (features, IQ/AQ scores, rescue and keep flags) and the three-column `loose_filtered_simplified` / `strict_filtered_simplified` edge lists (TF, peak, TG).

## Output Description

1. **Prepared data**: matched Seurat v5 objects, indexed fragments, TFE/TGA/TGE matrices.
2. **iGRNs**: global and per-branch TG × TF score matrices and the binarization threshold tau.
3. **Trained models**: per-pair optimal parameters, predicted/observed strengths, fitting loss, epochs, training time, and stop causes.
4. **pGRNs**: loose and strict TF–TG prediction tables with per-branch and integrated values.
5. **peGRNs**: loose and strict TF–peak–TG edge tables, both full (25 features, adjustment factors, filtering flags) and simplified (TF, peak, TG).

All functions are modular and independent; save every intermediate result with `qs::qsave` / `saveRDS` so any step can be rerun or resumed separately.

## Memory Management

- scDARI uses parallel (`mclapply`) execution; set `ncores = parallel::detectCores() - 1` on Linux and `ncores = 1` on Windows. Lower `ncores` if RAM is exhausted.
- Parallel workers may not release memory immediately; restart the R session and reload saved intermediates if needed.
- For large inference jobs, prefer branch-by-branch inference (Way 2 in Section 5.2) and use checkpoints (`save.image`) between branches.

## Common Issues and Solutions

1. **ChIP-Atlas prior download failures**: inspect `download_information.rds` for rows marked `"URL not found or download error"`, download them manually from the recorded URLs to the recorded destinations, and rerun `get_TFs_from_priors` (completed files are reused).
2. **Example dataset download failures**: extend `options(timeout)` and retry, or install `pbmcMultiome.SeuratData_0.1.4.tar.gz` locally with `install.packages(..., repos = NULL, type = "source")`.
3. **Too few (or zero) metacells, or failed grouping**: relax the anchor targets in `points_x/y_for_fitting_nc_and_nmin`; in `get_genes_pseudotime_info`, raise non-zero-ratio weights or increase `alpha_TF` / `alpha_TG` to retain fewer, better-supported genes. Reverse these changes when groups are too numerous.
4. **Training convergence problems**: increase `max_epochs` or `early_stop`; enlarge `popSize_ga` / `maxiter_ga` to improve global search; adjust post-training removal thresholds rather than forcing poorly fitted pairs.
5. **Out-of-memory errors**: reduce `ncores`, downsample cells/genes/peaks, or switch to branch-by-branch inference; use `qs` for fast large-object I/O.

## Citation

If you use scDARI in your research, please cite:

[... citation to be added upon publication ...]

## Related Links

1. R package: https://github.com/hth-buaa/scDARI
2. Full step-by-step tutorial: https://github.com/hth-buaa/scDARI/blob/main/scDARI_guide_line.R (also available as `scDARI_guide_line.txt`)
3. Paper code: https://github.com/hth-buaa/scDARI_paper_code
