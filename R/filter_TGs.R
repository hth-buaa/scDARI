#' @title Filter Target Genes by scRNA-seq Expression Quality Thresholds
#'
#' @description Compute per-gene expression quality metrics and their genome-wide percentile ranks from a scRNA-seq \code{counts} layer,
#' apply three absolute thresholds and three genome-wide quantile thresholds simultaneously to retain qualified target genes,
#' and return a sorted character vector of retained gene symbols.
#' Built-in disk caching accelerates repeated runs.
#'
#' @param scRNAseq Seurat object containing single-cell RNA-seq data with a \code{counts} layer,
#' used for target gene expression quantification and quality screening.
#' @param forced_feature_calculation Logical.
#' If \code{TRUE}, force recomputation of gene expression features and overwrite the cached RDS file.
#' Defaults to \code{FALSE}: the cached \code{TG_features_all.rds} is loaded when available to skip the computationally expensive feature calculation step.
#' @param min_TG_nonzero_freq Numeric between 0 and 1. Minimum absolute non-zero cell frequency. Defaults to \code{0.05}.
#' @param min_TG_nonzero_mean Numeric. Minimum absolute mean expression in non-zero cells. Defaults to \code{0} (no filter).
#' @param min_TG_overall_mean Numeric. Minimum absolute mean expression across all cells. Defaults to \code{0} (no filter).
#' @param min_TG_nonzero_freq_quantile Numeric between 0 and 1. Minimum genome-wide quantile of non-zero cell frequency. Defaults to \code{0} (no filter).
#' @param min_TG_nonzero_mean_quantile Numeric between 0 and 1. Minimum genome-wide quantile of mean expression in non-zero cells. Defaults to \code{0} (no filter).
#' @param min_TG_overall_mean_quantile Numeric between 0 and 1. Minimum genome-wide quantile of overall mean expression. Defaults to \code{0} (no filter).
#'
#' @returns Sorted character vector of target gene symbols that pass all expression quality thresholds.
#'
#' @details
#' The function implements a 3-step target gene filtering pipeline:
#' \enumerate{
#'   \item \strong{Feature quantification with cache}
#'   \itemize{
#'     \item Extract the raw gene-by-cell count matrix from the scRNA-seq Seurat object.
#'     \item Compute three per-gene metrics from raw counts: non-zero cell frequency,
#'       mean expression in non-zero cells, and overall mean expression across all cells.
#'     \item Convert each metric to genome-wide percentile ranks (0 = lowest, 1 = highest)
#'       against all genes, prefix all columns with \code{TG_}, and cache the feature table to disk.
#'     \item When \code{forced_feature_calculation = FALSE} and the cached RDS is present,
#'       load it directly to skip feature recomputation.
#'   }
#'   \item \strong{Expression quality filtering}
#'   \itemize{
#'     \item Retain genes satisfying all three absolute thresholds
#'       (\code{min_TG_nonzero_freq}, \code{min_TG_nonzero_mean}, \code{min_TG_overall_mean})
#'       and all three quantile thresholds (\code{min_TG_nonzero_freq_quantile},
#'       \code{min_TG_nonzero_mean_quantile}, \code{min_TG_overall_mean_quantile}) simultaneously.
#'   }
#'   \item \strong{Output}
#'   \itemize{
#'     \item Sort the retained target gene symbols and return them as a character vector.
#'   }
#' }
#'
#' @note
#' \itemize{
#'   \item Feature cache files (e.g. \code{TG_features_all.rds}) are stored under
#'     \code{./1.2 Data Preprocessing - Filter TGs/}, which is created automatically if absent.
#'   \item Percentile ranks are computed against all genes in the scRNA-seq count matrix,
#'     providing a genome-wide baseline for quality assessment.
#'   \item The input Seurat object must contain a \code{counts} layer for raw metric calculation.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load scRNA-seq data
#' scRNAseq = readRDS("./scRNAseq.rds")
#'
#' # Filter target genes with default parameters
#' TGs = filter_TGs(
#'   scRNAseq = scRNAseq,
#'   forced_feature_calculation = FALSE,
#'   min_TG_nonzero_freq = 0.05
#' )
#'
#' saveRDS(TGs, file = "./1.2 Data Preprocessing - Filter TGs/TGs.rds")
#' }
filter_TGs = function(
    scRNAseq,
    forced_feature_calculation = FALSE,
    min_TG_nonzero_freq = 0.05,
    min_TG_nonzero_mean = 0,
    min_TG_overall_mean = 0,
    min_TG_nonzero_freq_quantile = 0,
    min_TG_nonzero_mean_quantile = 0,
    min_TG_overall_mean_quantile = 0
) {
  ### Start.
  t_start = base::Sys.time()
  message("Run: Filtering peak annotation results. Time: ", t_start, ".")
  cache_dir = "./1.2 Data Preprocessing - Filter TGs/"
  if (!base::dir.exists(cache_dir)) base::dir.create(cache_dir, recursive = TRUE)

  ### Compute or load gene expression features (with disk cache).
  # Recompute features when forced by the user or when the cached table is missing.
  if (
    forced_feature_calculation ||
    !base::file.exists(base::file.path(cache_dir, "TG_features_all.rds"))
  ) {
    # Extract the raw gene-by-cell count matrix.
    message("Extract raw gene-by-cell count matrix.")
    gene_counts_all = Seurat::GetAssayData(scRNAseq, layer = "counts")

    # Compute three per-gene statistics: non-zero frequency, mean in non-zero cells, overall mean.
    message("Compute per-gene statistics: non-zero frequency, mean of non-zero values, overall mean.")
    n_cells_rna = base::ncol(gene_counts_all)
    gene_row_sums = Matrix::rowSums(gene_counts_all)
    gene_n_nonzero = Matrix::rowSums(gene_counts_all != 0)
    TG_features_all = base::data.frame(
      nonzero_freq = gene_n_nonzero / n_cells_rna,
      nonzero_mean = base::ifelse(gene_n_nonzero == 0, 0, gene_row_sums / gene_n_nonzero),
      overall_mean = gene_row_sums / n_cells_rna,
      row.names = base::rownames(gene_counts_all),
      check.names = FALSE
    )
    base::rm(gene_counts_all, gene_row_sums, gene_n_nonzero, n_cells_rna)

    # Convert the statistics to genome-wide percentile ranks, prefix columns with TG_, and cache.
    message("Add percentile ranks for each TG feature (0 = lowest, 1 = highest).")
    TG_features_all[base::paste0(base::names(TG_features_all), "_quantile")] =
      base::lapply(TG_features_all, function(x) base::rank(x) / base::length(x))
    base::colnames(TG_features_all) = base::paste0("TG_", base::colnames(TG_features_all))
    base::saveRDS(TG_features_all, base::file.path(cache_dir, "TG_features_all.rds"))
  } else {
    # Load the cached feature table to skip feature recomputation.
    message("Read feature data of genes and peaks.")
    TG_features_all = base::readRDS(base::file.path(cache_dir, "TG_features_all.rds"))
  }

  ### Filter target genes by expression quality.
  # Keep a gene only when it passes all three absolute and three quantile thresholds.
  message("Apply expression quality filters.")
  TG_features_all$keep = base::with(
    TG_features_all,
    TG_nonzero_freq >= min_TG_nonzero_freq &
      TG_nonzero_mean >= min_TG_nonzero_mean &
      TG_overall_mean >= min_TG_overall_mean &
      TG_nonzero_freq_quantile >= min_TG_nonzero_freq_quantile &
      TG_nonzero_mean_quantile >= min_TG_nonzero_mean_quantile &
      TG_overall_mean_quantile >= min_TG_overall_mean_quantile
  )

  # Sort and collect the target gene symbols that passed every threshold.
  message("Extract TGs passing all quality thresholds.")
  TGs_filtered_by_features = base::sort(base::rownames(TG_features_all[TG_features_all$keep, ]))

  ### End.
  t_end = base::Sys.time()
  message("Finish: Filtering peak annotation results. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(TGs_filtered_by_features)
}
