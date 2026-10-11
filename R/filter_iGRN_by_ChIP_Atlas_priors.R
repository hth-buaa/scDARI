#' @title Refine Initial GRNs with ChIP-Atlas Priors and Dual-Threshold QC
#'
#' @description Refines cell-type-specific initial gene regulatory networks (iGRNs) with preprocessed ChIP-Atlas ChIP-seq priors.
#' It performs prior-weighted edge correction, dual-threshold quality control with target gene (TG) rescue, and a zero-preserving value transformation,
#' balancing network sparsity with numerical stability for downstream dynamic modeling.
#'
#' @param interest_cell_type_iGRN Named list. Output of \code{\link{calculate_iGRN}}; each element
#' is a target gene (rows) by TF (columns) association matrix (data frame or matrix).
#' @param min_scale Positive numeric, typically in (0, 1], default \code{1}. Multiplier applied to
#' the positive minimum when building the min-max lower bound, used both before and after prior
#' correction. With \code{min_scale = 1} the weakest positive edge of each TF is mapped to zero;
#' values below 1 keep it positive.
#' @param theta_i_top_quantile Numeric in (0, 1]. Global quantile cutoff; edges whose corrected
#' score is at or above the \code{1 - theta_i_top_quantile} quantile are retained.
#' Defaults to \code{0.15} (global top 15%).
#' @param theta_i_top_n Positive integer. Per-TF top-N cutoff; the top \code{theta_i_top_n} targets
#' of each TF are retained regardless of the global threshold. Defaults to \code{50}.
#' @param theta_i_absolute_min Non-negative numeric. Score floor applied after the square-root
#' transform (before re-normalization); values below it are set to zero. Defaults to \code{0.05}.
#' @param TG_rescue_number Non-negative integer. Minimum number of TF regulators guaranteed for
#' each target gene after dual-threshold QC. TGs with fewer non-zero edges have their top pre-QC
#' (post-correction) edges restored. Internally capped at the number of TFs. Set \code{0} to disable.
#' Defaults to \code{1}.
#' @param ncores Positive integer. Number of cores for parallel processing via
#'   \code{\link[parallel]{mclapply}}.
#' Defaults to 1; set to 1 on Windows.
#'
#' @return A named list of the same structure as the input, containing refined iGRN data frames for
#' each cell type. Each matrix has been prior-weighted, quality-controlled and transformed:
#' \itemize{
#'   \item Non-zero edges lie approximately in (0.1, 0.92] (the upper bound follows from the
#'     1.1-fold re-normalization ceiling followed by \code{0.9*x + 0.1}).
#'   \item Zero edges remain exactly zero, preserving the sparse network structure.
#'   \item All-zero rows (target genes) and columns (TFs) are removed; a cell type with no
#'     prior-supported edge is returned as a 0 x 0 data frame.
#' }
#'
#' @details
#' The pipeline is applied independently to each cell type, in three stages. Cell types are
#' independent and are processed in parallel when a fork-capable platform is available. All
#' per-TF and per-TG operations are evaluated as vectorized matrix operations, which gives the
#' same numerical results as the element-wise implementation while being substantially faster.
#'
#' \enumerate{
#'   \item \strong{Prior-weighted edge correction}
#'     \itemize{
#'       \item Loads the preprocessed ChIP-Atlas priors from the standard output path of
#'         \code{\link{get_TFs_from_priors}}.
#'       \item For each TF, min-max normalizes the incoming scores with the
#'         \code{min_scale}-adjusted lower bound, flooring negative values to zero.
#'       \item Multiplies the normalized scores by ChIP-seq prior intensity; TFs (or TGs) without
#'         prior support receive zero intensity.
#'       \item Applies a square-root transform and floors values below \code{theta_i_absolute_min}.
#'       \item Re-normalizes with the same \code{min_scale} offset and a 1.1-fold maximum expansion.
#'     }
#'
#'   \item \strong{Dual-threshold QC with target gene rescue}
#'     \itemize{
#'       \item Computes a global threshold as the \code{1 - theta_i_top_quantile} quantile of all
#'         positive corrected scores.
#'       \item For each TF, keeps targets passing either the per-TF top-\code{theta_i_top_n} rank or
#'         the global quantile threshold; all other edges are zeroed.
#'       \item If \code{TG_rescue_number > 0}, restores the top pre-QC edges of any TG left with
#'         fewer than that many non-zero regulators.
#'       \item Removes all-zero rows and columns.
#'     }
#'
#'   \item \strong{Zero-preserving minimum value correction}
#'     Applies \code{0.9*x + 0.1} and resets values \eqn{\le 0.1} to exactly zero, shifting
#'     non-zero edges above 0.1 while preserving sparsity.
#' }
#'
#' The prior file is read from
#' \code{../1.1 Data Preprocessing - Get TFs From Database/ChIP_Atlas_Priors/simplified_processed_priors.rds}
#' after switching into the \code{3 get iGRN/} working directory.
#'
#' @note
#' \itemize{
#'   \item Depends on the preprocessed ChIP-Atlas priors from \code{\link{get_TFs_from_priors}};
#'     \code{simplified_processed_priors.rds} must exist before running this function.
#'   \item The function creates (if missing) and enters \code{3 get iGRN/}, restoring the original
#'     working directory on exit.
#'   \item Input validation covers list structure, parameter ranges and prior-file existence;
#'     invalid inputs raise informative errors.
#'   \item TFs present in the iGRN but absent from the prior table receive zero intensity and are
#'     removed after QC.
#'   \item The per-TF top-N rule keeps a minimal target set per TF, while TG rescue prevents target
#'     genes from losing all regulators during aggressive sparsification.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
#' ncores = if (.Platform$OS.type == "unix") {
#'   base::max(1L, base::min(
#'     parallel::detectCores(logical = FALSE) - 1L,
#'     base::length(interest_cell_type_iGRN)
#'   ))
#' } else {
#'   1L
#' }
#' interest_cell_type_iGRN = filter_iGRN_by_ChIP_Atlas_priors(
#'   interest_cell_type_iGRN = interest_cell_type_iGRN,
#'   min_scale = 0.2,
#'   theta_i_top_quantile = 0.15,
#'   theta_i_top_n = 50,
#'   theta_i_absolute_min = 0.05,
#'   TG_rescue_number = 1,
#'   ncores = ncores
#' )
#'
#' base::saveRDS(
#'   interest_cell_type_iGRN,
#'   file = "./3 get iGRN/interest_cell_type_iGRN_refined.rds"
#' )
#' }
filter_iGRN_by_ChIP_Atlas_priors = function(
    interest_cell_type_iGRN,
    min_scale = 1,
    theta_i_top_quantile = 0.15,
    theta_i_top_n = 50,
    theta_i_absolute_min = 0.05,
    TG_rescue_number = 1,
    ncores = 1
) {
  ### Validate inputs.
  if (!base::is.list(interest_cell_type_iGRN) || base::length(interest_cell_type_iGRN) == 0) {
    base::stop("'interest_cell_type_iGRN' must be a non-empty named list of iGRN matrices.")
  }
  for (cell_type in base::names(interest_cell_type_iGRN)) {
    if (!base::is.matrix(interest_cell_type_iGRN[[cell_type]]) &&
        !base::is.data.frame(interest_cell_type_iGRN[[cell_type]])) {
      base::stop(base::sprintf("Element '%s' of 'interest_cell_type_iGRN' must be a matrix or data frame.", cell_type))
    }
  }
  if (theta_i_top_quantile <= 0 || theta_i_top_quantile > 1) {
    base::stop("'theta_i_top_quantile' must be in the range (0, 1].")
  }
  if (!base::is.numeric(theta_i_top_n) || base::length(theta_i_top_n) != 1 ||
      theta_i_top_n < 1 || theta_i_top_n != base::round(theta_i_top_n)) {
    base::stop("'theta_i_top_n' must be a single positive integer.")
  }
  if (!base::is.numeric(theta_i_absolute_min) || base::length(theta_i_absolute_min) != 1 ||
      theta_i_absolute_min < 0) {
    base::stop("'theta_i_absolute_min' must be a single non-negative number.")
  }
  if (!base::is.numeric(TG_rescue_number) || base::length(TG_rescue_number) != 1 ||
      TG_rescue_number < 0 || TG_rescue_number != base::round(TG_rescue_number)) {
    base::stop("'TG_rescue_number' must be a single non-negative integer.")
  }
  if (!base::is.numeric(min_scale) || base::length(min_scale) != 1 ||
      base::is.na(min_scale) || min_scale <= 0) {
    base::stop("'min_scale' must be a single positive number.")
  }
  if (min_scale > 1) {
    base::warning("'min_scale' > 1; the normalized lower bound will exceed the observed minimum.")
  }

  ### Start.
  t_start = base::Sys.time()
  message("Run: Filtering the global iGRN by ChIP-Atlas priors for each cell type. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)

  # Restore the original working directory even if the function exits early on error.
  base::on.exit(base::setwd(original_dir), add = TRUE)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Read the preprocessed ChIP-Atlas priors once for all cell types.
  priors_path = base::file.path(
    "..", "1.1 Data Preprocessing - Get TFs From Database",
    "ChIP_Atlas_Priors", "simplified_processed_priors.rds"
  )
  if (!base::file.exists(priors_path)) {
    base::stop("Preprocessed prior file not found: ", priors_path,
               ". Please run get_TFs_from_priors() first.")
  }
  message("Reading ChIP-Atlas priors.")
  priors_simplify = base::readRDS(priors_path)
  if (!base::is.null(priors_simplify) && !base::is.list(priors_simplify)) {
    base::stop("The prior file must contain a list keyed by TF (see get_TFs_from_priors()).")
  }
  rescue_n = TG_rescue_number # Keep the requested rescue count constant across cell types.
  top_n = base::as.integer(theta_i_top_n)

  # Refine one cell type using vectorized matrix operations (same math as the column/row loops).
  refine_cell_type = function(igrn) {
    if (!base::is.data.frame(igrn)) {
      igrn = base::as.data.frame(igrn, optional = TRUE)
    }
    M = base::as.matrix(igrn)
    base::storage.mode(M) = "double"
    nTG = base::nrow(M)
    nTF = base::ncol(M)
    if (nTG == 0L || nTF == 0L) {
      return(igrn[FALSE, FALSE, drop = FALSE])
    }
    tg_names = base::rownames(M)
    tf_names = base::colnames(M)

    # Build the TG x TF ChIP-seq prior-intensity matrix (missing entries become zero).
    P = base::matrix(0, nrow = nTG, ncol = nTF,
                     dimnames = base::list(tg_names, tf_names))
    for (tf in tf_names) {
      pr = priors_simplify[[tf]]
      if (!base::is.null(pr)) {
        vals = pr[tg_names, "intensity"]
        vals[base::is.na(vals)] = 0
        P[, tf] = vals
      }
    }

    # Stage 1: prior-weighted correction, vectorized across all TFs.
    # Per-TF min-max normalization of the incoming scores (non-positive entries ignored).
    mx_before = base::apply(M, 2L, base::max)
    M_pos = M
    M_pos[M_pos <= 0] = NA_real_
    mn_before = base::suppressWarnings(base::apply(M_pos, 2L, base::min, na.rm = TRUE))
    max_before = mx_before + 0.0001
    min_before = mn_before * min_scale
    norm_before = base::sweep(M, 2L, min_before, "-")
    norm_before = base::sweep(norm_before, 2L, max_before - min_before, "/")
    norm_before[norm_before < 0] = 0
    norm_before[, mx_before == 0] = 0
    norm_before[!base::is.finite(norm_before)] = 0

    # Multiply by prior intensity, square-root transform, and floor weak values.
    theta_after = norm_before * P
    theta_after = theta_after ^ 0.5
    theta_after[theta_after < theta_i_absolute_min] = 0

    # Re-normalize each TF with the same offset and a 1.1-fold ceiling.
    mx_after = base::apply(theta_after, 2L, base::max)
    A_pos = theta_after
    A_pos[A_pos <= 0] = NA_real_
    mn_after = base::suppressWarnings(base::apply(A_pos, 2L, base::min, na.rm = TRUE))
    max_after = mx_after + 0.0001
    min_after = mn_after * min_scale
    M2 = base::sweep(theta_after, 2L, min_after, "-")
    M2 = base::sweep(M2, 2L, max_after * 1.1 - min_after, "/")
    M2[M2 < 0] = 0
    M2[, mx_after == 0] = 0
    M2[!base::is.finite(M2)] = 0

    # Without any prior-supported edge the matrix cannot be thresholded; return 0 x 0.
    if (!base::any(M2 > 0)) {
      message("No prior-supported edge; returning a 0 x 0 matrix.")
      return(igrn[FALSE, FALSE, drop = FALSE])
    }

    # Stage 2a: dual-threshold QC (per-TF top-N union global quantile).
    positive_theta = M2[M2 > 0]
    min_theta = stats::quantile(x = positive_theta, probs = 1 - theta_i_top_quantile)

    # Per-TF cutoff equals the top-N-th largest positive value (ties included, matching
    # rank(..., ties.method = "min")); TFs with fewer than N positives keep all of them.
    boundary = base::apply(M2, 2L, function(v) {
      vp = v[v > 0]
      nv = base::length(vp)
      if (nv == 0L) {
        Inf
      } else if (nv >= top_n) {
        base::sort(vp, decreasing = TRUE)[top_n]
      } else {
        base::min(vp)
      }
    })
    cutoff = base::pmin(min_theta, boundary)
    keep = base::sweep(M2, 2L, cutoff, ">=")
    Mqc = M2
    Mqc[!keep] = 0

    # Pre-compute each TG's top pre-QC (post-correction) edges for the rescue step.
    k_rescue = base::min(rescue_n, nTF)
    if (rescue_n > 0) {
      if (k_rescue == 1L) {
        # Vectorized argmax; ties resolve to the first TF, matching order(-v)[1].
        top_idx = base::matrix(base::max.col(M2, ties.method = "first"), ncol = 1L)
      } else {
        top_idx = base::t(base::apply(M2, 1L, function(v) base::order(-v)[base::seq_len(k_rescue)]))
      }
      row_rep = base::rep(base::seq_len(nTG), each = k_rescue)
      top_col = base::as.vector(base::t(top_idx))
      top_val = base::matrix(M2[base::cbind(row_rep, top_col)],
                             nrow = nTG, ncol = k_rescue, byrow = TRUE)

      # Stage 2b: rescue TGs left with fewer than k_rescue non-zero regulators.
      n_nonzero = base::rowSums(Mqc > 0)
      rescue_rows = base::which(n_nonzero < k_rescue)
      if (base::length(rescue_rows)) {
        ri = base::rep(rescue_rows, each = k_rescue)
        ci = base::as.vector(base::t(top_idx[rescue_rows, , drop = FALSE]))
        vi = base::as.vector(base::t(top_val[rescue_rows, , drop = FALSE]))
        Mqc[base::cbind(ri, ci)] = vi
      }
    }

    # Drop all-zero rows (TGs) and columns (TFs).
    keep_row = base::rowSums(Mqc != 0) > 0
    keep_col = base::colSums(Mqc != 0) > 0
    Mqc = Mqc[keep_row, keep_col, drop = FALSE]

    # Stage 3: shift non-zero values above 0.1 while resetting exact zeros.
    Mqc = Mqc * 0.9 + 0.1
    Mqc[Mqc <= 0.1] = 0
    base::as.data.frame(Mqc)
  }

  ### Refine every cell type (in parallel on fork-capable platforms).
  refined = parallel::mclapply(
    interest_cell_type_iGRN, refine_cell_type, mc.cores = ncores
  )
  base::names(refined) = base::names(interest_cell_type_iGRN)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Filtering the global iGRN by ChIP-Atlas priors for each cell type. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(refined)
}
