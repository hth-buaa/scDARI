#' @title Construct Cell-Type-Specific Initial GRNs (iGRNs) from kNN Mutual Information and CLR
#'
#' @description For each cell type,
#' this function estimates the pairwise mutual information (MI) between target genes (TGs, rows) and transcription factors (TFs, columns) from log-normalized scRNA-seq expression,
#' using the Kraskov k-nearest-neighbor MI estimator (averaged over several noise-perturbed replicates).
#' It then applies the Context Likelihood of Relatedness (CLR) background correction and rank-normalizes the scores per TF,
#' producing an initial TG-by-TF association matrix for downstream Gene Regulatory Network (GRN) modeling.
#'
#' @param interest_cell_type_data List. Output from \code{\link{get_interest_cell_type_data}},
#' providing the per-cell-type Seurat object \code{interest_scRNAseq_for_GRN}.
#' @param interest_cell_type_genes_pseudotime_info List. Output from \code{\link{get_genes_pseudotime_info}},
#' providing the retained TF/TG sets (\code{Branches_TFs}, \code{Branches_TGs}) and the branch cell tables (\code{Branches}).
#' @param replicate_time Positive integer.
#' Number of kNN-MI replicates averaged to stabilize the estimate against the random tie-breaking noise.
#' Defaults to \code{10}.
#' @param k Positive integer. Number of nearest neighbors to consider to estimate the mutual information used by the Kraskov MI estimator,
#' passed to \code{\link[parmigene]{knnmi.cross}}.
#' It must be smaller than the number of cells (samples) in the pooled cell-type matrix.
#' Defaults to \code{3}.
#' @param noise Non-negative numeric.
#' Magnitude of the random perturbation noise added before nearest-neighbor distance computation to break ties,
#' passed to \code{\link[parmigene]{knnmi.cross}}.
#' Defaults to \code{1e-12}.
#' @param seed Single integer passed to \code{\link[base]{set.seed}} before the MI replicates to ensure reproducibility,
#' or \code{NULL} to leave the random number generator untouched.
#' Defaults to \code{123}.
#'
#' @return A named list with one element per cell type. Each element is a data frame used as an association matrix:
#' \itemize{
#'   \item Rows: target gene symbols (TGs)
#'   \item Columns: transcription factor symbols (TFs)
#'   \item Values: per-TF rank-normalized CLR scores in (0, 1] for supported edges;
#'     \code{0} indicates no supported TG-TF association
#' }
#'
#' @details
#' For each cell type the workflow is:
#' \enumerate{
#'   \item Pools the TF/TG sets retained across all branches and the cells assigned to the branches,
#'     and subsets the cell-type Seurat object to the genes and cells that are present.
#'   \item Log-normalizes the subset (an idempotent safeguard when the input object has no normalized layer) and extracts the \code{data} layer;
#'     log transformation compresses the heavy count tail and reduces ties,
#'     while MI is invariant to invertible monotonic transformations.
#'   \item Estimates the cross-MI matrix TG-by-TF with \code{\link[parmigene]{knnmi.cross}} and averages \code{replicate_time} noise-perturbed replicates.
#'   \item Applies CLR: z-scores each TG against its TF background (rows) and each TF against its TG background (columns),
#'     truncates negative z-scores to zero,
#'     and combines the two sides as \eqn{\sqrt{z_{TG}^{2} + z_{TF}^{2}}} (Faith et al., 2007).
#'   \item Rank-normalizes the positive scores independently within each TF column to (0, 1] (\code{rank/n} with average ties),
#'     leaving unsupported edges at \code{0}.
#' }
#'
#' @note
#' \itemize{
#'   \item Set a random seed (e.g. \code{set.seed(123)}) before calling this function for reproducibility,
#'     because the kNN estimator uses random tie-breaking noise in every replicate.
#'   \item Rank normalization is performed per TF (column): the strongest TG of every TF receives 1,
#'     so absolute scores are not directly comparable across TFs.
#'   \item The Seurat object must contain a normalized \code{data} layer (normalization is run upstream).
#'   \item This function uses scRNA-seq expression only;
#'     \code{\link[parmigene]{knnmi.cross}} is a compiled estimator.
#'   \item Required dependencies: \code{parmigene} and \code{Seurat} (with \code{Matrix} support).
#' }
#'
#' @note
#' \itemize{
#'   \item Pass \code{seed} (default 123) for reproducibility,
#'     because the kNN estimator adds random tie-breaking noise in every replicate.
#'   \item Rank normalization is per TF (column);
#'     the highest positive score approaches 1 (it is strictly 1 only when the maximum is unique, because ties use average ranks).
#'     Absolute scores are therefore not directly comparable across TFs.
#'   \item Pooled cells come from the branch tables (\code{Branches});
#'     they may include cells removed by the per-branch \code{alpha_cell} filter upstream,
#'     while the gene sets come from the retained branches.
#'   \item This function uses scRNA-seq expression only and runs cell types sequentially;
#'     \code{\link[parmigene]{knnmi.cross}} is a compiled estimator.
#'   \item Required dependencies: \code{parmigene} and \code{Seurat}.
#' }
#'
#' @references
#' Sales, G. & Romualdi, C. (2011). parmigene: a parallel R package for mutual information
#' identification of gene networks. \emph{Bioinformatics}, 27(13), 1876-1877.
#'
#' Faith, J. J. et al. (2007). Large-scale mapping and validation of Escherichia coli
#' transcriptional regulation from a compendium of expression profiles. \emph{PLoS Biology}, 5(1), e8.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_data = base::readRDS(base::paste0(
#'   "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'   "/interest_cell_type_data.rds"
#' ))
#' interest_cell_type_genes_pseudotime_info = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#'
#' interest_cell_type_iGRN = calculate_iGRN_old(
#'   interest_cell_type_data = interest_cell_type_data,
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   replicate_time = 10,
#'   k = 3,
#'   noise = 1e-12,
#'   seed = 123
#' )
#'
#' base::saveRDS(
#'   interest_cell_type_iGRN,
#'   file = "./3 get iGRN/interest_cell_type_iGRN.rds"
#' )
#' }
calculate_iGRN_old = function(
    interest_cell_type_data,
    interest_cell_type_genes_pseudotime_info,
    replicate_time = 10,
    k = 3,
    noise = 1e-12,
    seed = 123
) {
  ### Validate parameters.
  if (!(base::is.numeric(replicate_time) && base::length(replicate_time) == 1 &&
        replicate_time >= 1 && replicate_time == base::round(replicate_time))) {
    base::stop("replicate_time must be a single positive integer.")
  }
  if (!(base::is.numeric(k) && base::length(k) == 1 && k >= 1 && k == base::round(k))) {
    base::stop("k must be a single positive integer.")
  }
  if (!(base::is.numeric(noise) && base::length(noise) == 1 && noise >= 0)) {
    base::stop("noise must be a single non-negative number.")
  }
  if (!(base::is.null(seed) ||
        (base::is.numeric(seed) && base::length(seed) == 1 && seed == base::round(seed)))) {
    base::stop("seed must be a single integer or NULL.")
  }
  # The two inputs must describe the same set of cell types.
  cell_types = base::names(interest_cell_type_genes_pseudotime_info)
  missing_in_data = base::setdiff(cell_types, base::names(interest_cell_type_data))
  if (base::length(missing_in_data) > 0) {
    base::stop("interest_cell_type_data is missing cell type(s): ",
               base::paste(missing_in_data, collapse = ", "), ".")
  }
  for (cell_type in cell_types) {
    if (base::is.null(interest_cell_type_data[[cell_type]][["interest_scRNAseq_for_GRN"]])) {
      base::stop("interest_cell_type_data[['", cell_type,
                 "']] must contain 'interest_scRNAseq_for_GRN'.")
    }
    info = interest_cell_type_genes_pseudotime_info[[cell_type]]
    if (base::is.null(info[["Branches_TFs"]]) || base::is.null(info[["Branches_TGs"]]) ||
        base::is.null(info[["Branches"]])) {
      base::stop("interest_cell_type_genes_pseudotime_info[['", cell_type,
                 "']] must contain Branches_TFs, Branches_TGs and Branches.")
    }
  }

  ### Start.
  t_start = base::Sys.time()
  message("Run: Constructing Initial Gene Regulatory Networks (iGRNs) from mutual information. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Define shared helpers.
  # CLR background correction for a rectangular TG (rows) x TF (cols) MI matrix.
  # Each TG is z-scored across its TF background and each TF across its TG background;
  # negative z-scores contribute nothing, and the two sides are combined in quadrature.
  CLR_rectangular = function(M) {
    z_row = base::t(base::suppressWarnings(base::scale(base::t(M)))) # Row (TG) background.
    z_col = base::suppressWarnings(base::scale(M))                  # Column (TF) background.
    z_row[base::is.na(z_row)] = 0 # Zero-variance background contributes nothing on that side.
    z_col[base::is.na(z_col)] = 0 # Zero-variance background contributes nothing on that side.
    base::sqrt(base::pmax(z_row, 0)^2 + base::pmax(z_col, 0)^2)
  }

  # Per-column (per-TF) empirical rank normalization of positive values to (0, 1];
  # zero/NA entries stay at 0.
  rank01_col = function(M) {
    out = base::apply(M, MARGIN = 2, function(Col) {
      ok = !base::is.na(Col) & Col > 0
      n_ok = base::sum(ok)
      if (n_ok == 0) {
        return(base::rep(0, base::length(Col)))
      }
      Col[ok] = base::rank(Col[ok], ties.method = "average") / n_ok
      Col[!ok] = 0
      return(Col)
    })
    base::dimnames(out) = base::dimnames(M)
    return(base::as.data.frame(out))
  }

  if (!base::is.null(seed)) {base::set.seed(seed)}

  ### Build the iGRN of each cell type.
  interest_cell_type_iGRN_all = base::list()
  for (cell_type in cell_types) {
    message("Constructing iGRN for cell type ", cell_type, ".")

    # Pool the retained TF/TG sets across branches and the cells assigned to branches.
    info = interest_cell_type_genes_pseudotime_info[[cell_type]]
    TFs = base::sort(base::unique(base::Reduce(base::union, info[["Branches_TFs"]])))
    TGs = base::sort(base::unique(base::Reduce(base::union, info[["Branches_TGs"]])))
    genes = base::sort(base::union(TFs, TGs))
    cells = base::unique(base::unlist(base::lapply(info[["Branches"]], base::rownames)))
    base::rm(info)

    # Restrict to genes/cells actually present in the cell-type Seurat object.
    scRNA_obj = interest_cell_type_data[[cell_type]][["interest_scRNAseq_for_GRN"]]
    keep_genes = base::intersect(genes, base::rownames(scRNA_obj))
    keep_cells = base::intersect(cells, base::colnames(scRNA_obj))
    dropped_genes = base::length(genes) - base::length(keep_genes)
    dropped_cells = base::length(cells) - base::length(keep_cells)
    if (dropped_genes > 0 || dropped_cells > 0) {
      message("[", cell_type, "] Dropped ", dropped_genes, " gene(s) and ",
              dropped_cells, " cell(s) absent from the Seurat object.")
    }
    TFs_keep = base::intersect(TFs, keep_genes)
    TGs_keep = base::intersect(TGs, keep_genes)
    if (base::length(TFs_keep) == 0 || base::length(TGs_keep) == 0 || base::length(keep_cells) == 0) {
      base::stop("[", cell_type, "] no usable TF, TG or cell left to estimate mutual information.")
    }
    interest_scRNAseq = scRNA_obj[keep_genes, keep_cells]
    base::rm(scRNA_obj)

    # Guarantee a log-normalized 'data' layer (Seurat v5 'layer', v4 'slot' fallback).
    interest_scRNAseq = Seurat::NormalizeData(interest_scRNAseq, verbose = FALSE)
    expr = base::tryCatch(
      base::as.matrix(SeuratObject::GetAssayData(interest_scRNAseq, layer = "data")),
      error = function(e) base::as.matrix(SeuratObject::GetAssayData(interest_scRNAseq, slot = "data"))
    )
    # Kraskov takes neighbors among cells (samples), so k must be below the cell count.
    if (k >= base::ncol(expr)) {
      base::stop("[", cell_type, "] k = ", k, " must be smaller than the number of cells (",
                 base::ncol(expr), ").")
    }

    # Kraskov kNN cross-MI between TG rows and TF columns, averaged over noise replicates.
    message("[", cell_type, "] Estimating TG-TF mutual information (",
            base::length(TGs_keep), " TGs x ", base::length(TFs_keep), " TFs x ",
            base::ncol(expr), " cells; ", replicate_time, " replicates, k = ", k, ").")
    mi_replicates = base::replicate(
      n = replicate_time,
      expr = parmigene::knnmi.cross(
        mat1 = expr[TGs_keep, , drop = FALSE],
        mat2 = expr[TFs_keep, , drop = FALSE],
        k = k,
        noise = noise
      ),
      simplify = FALSE
    )
    MI_cross_avg = base::Reduce("+", mi_replicates) / replicate_time

    # Guarantee TG rownames and TF colnames regardless of the estimator's returned dimnames.
    base::rownames(MI_cross_avg) = TGs_keep
    base::colnames(MI_cross_avg) = TFs_keep

    # CLR background correction, then per-TF rank normalization to (0, 1].
    message("[", cell_type, "] Applying CLR correction and rank normalization.")
    CLR_cross = CLR_rectangular(MI_cross_avg)
    interest_cell_type_iGRN_all[[cell_type]] = rank01_col(CLR_cross)
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Constructing Initial Gene Regulatory Networks (iGRNs) from mutual information. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_iGRN_all)
}
