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
#' The pipeline is applied independently to each cell type, in three stages.
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
#'
#' interest_cell_type_iGRN = filter_iGRN_by_ChIP_Atlas_priors(
#'   interest_cell_type_iGRN = interest_cell_type_iGRN,
#'   min_scale = 0.2,
#'   theta_i_top_quantile = 0.15,
#'   theta_i_top_n = 50,
#'   theta_i_absolute_min = 0.05,
#'   TG_rescue_number = 1
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
    TG_rescue_number = 1
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
  if (!base::is.numeric(min_scale) || base::length(min_scale) == 1 && min_scale <= 0) {
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
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Read the preprocessed ChIP-Atlas priors.
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

  ### Refine each cell type.
  for (cell_type in base::names(interest_cell_type_iGRN)) {
    message("[", cell_type, "] Correcting iGRN weights with prior intensity.")
    igrn = interest_cell_type_iGRN[[cell_type]]
    if (!base::is.data.frame(igrn)) {
      igrn = base::as.data.frame(igrn, optional = TRUE)
    }
    nTF = base::ncol(igrn)

    # Stage 1: prior-weighted correction, one TF column at a time.
    for (tf in base::colnames(igrn)) {
      tf_tgs = igrn[, tf, drop = FALSE]
      base::colnames(tf_tgs)[1] = "theta_i_before_fix"

      # Min-max normalize the incoming scores with the min_scale-adjusted lower bound.
      if (base::max(tf_tgs$theta_i_before_fix) == 0) {
        tf_tgs$norm_theta_i_before_fix = 0
      } else {
        max_before = base::max(tf_tgs$theta_i_before_fix) + 0.0001
        min_before = base::min(tf_tgs$theta_i_before_fix[tf_tgs$theta_i_before_fix > 0]) * min_scale
        tf_tgs$norm_theta_i_before_fix =
          (tf_tgs$theta_i_before_fix - min_before) / (max_before - min_before)
        tf_tgs$norm_theta_i_before_fix[tf_tgs$norm_theta_i_before_fix < 0] = 0
      }

      # Multiply by ChIP-seq prior intensity (missing TF/TG entries receive zero).
      if (!base::is.null(priors_simplify[[tf]])) {
        tf_tgs$prior_intensity = priors_simplify[[tf]][base::rownames(tf_tgs), "intensity"]
        tf_tgs$prior_intensity[base::is.na(tf_tgs$prior_intensity)] = 0
      } else {
        tf_tgs$prior_intensity = 0
      }

      # Square-root transform, then floor weak values.
      tf_tgs$theta_i_after_fix = tf_tgs$norm_theta_i_before_fix * tf_tgs$prior_intensity
      tf_tgs$theta_i_after_fix = tf_tgs$theta_i_after_fix ^ 0.5
      tf_tgs$theta_i_after_fix[tf_tgs$theta_i_after_fix < theta_i_absolute_min] = 0

      # Re-normalize with the same offset and a 1.1-fold ceiling so values stay below 1.
      if (base::max(tf_tgs$theta_i_after_fix) > 0) {
        max_after = base::max(tf_tgs$theta_i_after_fix) + 0.0001
        min_after = base::min(tf_tgs$theta_i_after_fix[tf_tgs$theta_i_after_fix > 0]) * min_scale
        tf_tgs$norm_theta_i_after_fix =
          (tf_tgs$theta_i_after_fix - min_after) / (max_after * 1.1 - min_after)
        tf_tgs$norm_theta_i_after_fix[tf_tgs$norm_theta_i_after_fix < 0] = 0
      } else {
        tf_tgs$norm_theta_i_after_fix = 0
      }

      igrn[, tf] = tf_tgs[base::rownames(igrn), "norm_theta_i_after_fix"]
    }

    # Without any prior-supported edge the matrix cannot be thresholded; return 0 x 0.
    positive_theta = base::unlist(igrn)
    positive_theta = positive_theta[positive_theta > 0]
    if (base::length(positive_theta) == 0) {
      message("[", cell_type, "] No prior-supported edge; returning a 0 x 0 matrix.")
      interest_cell_type_iGRN[[cell_type]] = igrn[FALSE, FALSE, drop = FALSE]
      next
    }

    # Stage 2a: dual-threshold QC (per-TF top-N union global quantile).
    message("[", cell_type, "] QC: keep per-TF top-", theta_i_top_n,
            " or global top ", theta_i_top_quantile * 100, "%.")
    min_theta = stats::quantile(x = positive_theta, probs = 1 - theta_i_top_quantile)
    TGs = base::rownames(igrn)

    # Pre-compute each TG's top pre-QC (post-correction) edges for the rescue step.
    k_rescue = base::min(rescue_n, nTF)
    if (rescue_n > 0) {
      top_tf_cols = base::paste0("top_TF", base::seq_len(k_rescue))
      top_theta_cols = base::paste0("top_theta_i", base::seq_len(k_rescue))
      TG_rescue = base::data.frame(
        base::matrix(NA_character_, nrow = base::length(TGs), ncol = k_rescue,
                     dimnames = base::list(TGs, top_tf_cols)),
        base::matrix(0, nrow = base::length(TGs), ncol = k_rescue,
                     dimnames = base::list(TGs, top_theta_cols)),
        check.names = FALSE
      )
      for (tg in TGs) {
        top_TFs = base::colnames(igrn)[base::order(-base::as.numeric(igrn[tg, ]))][base::seq_len(k_rescue)]
        TG_rescue[tg, top_tf_cols] = top_TFs
        TG_rescue[tg, top_theta_cols] = base::as.numeric(igrn[tg, top_TFs])
      }
    }

    # Per TF: retain edges passing either the top-N rank or the global quantile threshold.
    for (tf in base::colnames(igrn)) {
      tf_tgs_fixed = igrn[, tf, drop = FALSE]
      tf_tgs_fixed = tf_tgs_fixed[tf_tgs_fixed[, 1] > 0, , drop = FALSE]
      TGs_positive = base::rownames(tf_tgs_fixed)
      tf_tgs_fixed$rank = base::rank(-tf_tgs_fixed[, 1], ties.method = "min")
      TGs_by_top = TGs_positive[tf_tgs_fixed$rank <= theta_i_top_n]
      TGs_by_min = TGs_positive[tf_tgs_fixed[, 1] >= min_theta]
      TGs_remaining = base::union(TGs_by_top, TGs_by_min)
      igrn[base::setdiff(TGs, TGs_remaining), tf] = 0
    }

    # Stage 2b: rescue TGs left with fewer than k_rescue non-zero regulators.
    if (rescue_n > 0) {
      message("[", cell_type, "] QC: rescuing TGs with fewer than ", k_rescue, " non-zero TF edge(s).")
      for (tg in TGs) {
        if (base::sum(base::as.numeric(igrn[tg, ]) > 0) < k_rescue) {
          for (rn in base::seq_len(k_rescue)) {
            igrn[tg, TG_rescue[tg, base::paste0("top_TF", rn)]] =
              TG_rescue[tg, base::paste0("top_theta_i", rn)]
          }
        }
      }
    }

    # Drop all-zero rows (TGs) and columns (TFs).
    igrn = igrn[base::rowSums(igrn != 0) > 0, base::colSums(igrn != 0) > 0, drop = FALSE]

    # Stage 3: shift non-zero values above 0.1 while resetting exact zeros.
    igrn = igrn * 0.9 + 0.1
    igrn[igrn <= 0.1] = 0

    interest_cell_type_iGRN[[cell_type]] = igrn
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Filtering the global iGRN by ChIP-Atlas priors for each cell type. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_iGRN)
}
