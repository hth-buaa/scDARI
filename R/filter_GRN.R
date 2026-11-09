#' @title Filter GRN by Feature Thresholding, Importance Scoring and TF Rescue
#'
#' @description Converts the continuous predictions into a final edge set without
#' fitting a global threshold curve. Each normalized feature is first aligned to
#' a common direction, truncated at a per-feature quantile and rescaled over its
#' positive values. The prepared features are combined into a single importance
#' score, whose high quantile defines the retained edges, and a TF rescue step
#' guarantees a minimum number of targets per regulator. Loose and strict
#' filtered sets are returned.
#'
#' @param interest_cell_type_pGRN List. Output of
#'   \code{\link{select_best_prediction}}, named by cell type and containing
#'   \code{pred_summary}.
#' @param feat_cut_off_threshold Named list of nine numeric values in
#'   \code{[0, 1]}. Names are the features used throughout the function.
#'   Each value is a quantile probability: for that feature, pairs below
#'   the corresponding quantile are set to zero. Larger values remove more
#'   low-ranked pairs. Closely related to \code{feat_scale_min}, which rescales
#'   the values that survive this cutoff, and to \code{feat_weights},
#'   whose names must match these names.
#' @param feat_scale_min Named list of nine numeric values in \code{[0, 1]}.
#'   After truncation, the positive values of a feature are min-max rescaled
#'   onto \code{[feat_scale_min, 1]}, so this value is the smallest score
#'   retained by a positive feature. A value of \code{0} keeps the natural
#'   minimum at zero, whereas \code{1} collapses all positive values to one.
#'   Related to \code{feat_cut_off_threshold}, which determines which values
#'   are positive before rescaling.
#' @param feat_weights Named list of nine numeric values in \code{[0, 1]},
#'   which sums to one. Each value is the weight of its feature in the weighted
#'   sum that forms the importance score. The weights should be normalized to
#'   sum to one so the score stays on the feature scale; the function does not
#'   enforce this, and a non-unit sum only rescales the importance score without
#'   changing its ranking. Names must match those of
#'   \code{feat_cut_off_threshold}.
#' @param IS_quan_cut_off Numeric in \code{[0, 1]}, default \code{0.02}.
#'   Quantile probability of the importance score: pairs whose score is below
#'   this quantile are rejected. Larger values retain fewer edges. Drives
#'   \code{filter_keep} and interacts with the TF rescue step.
#' @param rescue_TF Logical, default \code{TRUE}. Whether to run the TF rescue
#'   step after the importance-score cutoff. When \code{FALSE},
#'   \code{min_TG_per_TF} and \code{TF_rescue_quantile_threshold} are unused.
#' @param min_TG_per_TF Positive integer, default \code{1}. Minimum number of
#'   targets kept for each TF after rescue. Only relevant when
#'   \code{rescue_TF = TRUE}.
#' @param TF_rescue_quantile_threshold Numeric in \code{[0, 1]}, default
#'   \code{0.00}. Pooled low quantile of \code{FG.theta_p_norm.R} among rescue
#'   candidates; candidates below it are rejected. Only relevant when
#'   \code{rescue_TF = TRUE}.
#'
#' @return The input list with the following per cell type:
#' \itemize{
#'   \item \code{pred_summary}: the original frame augmented with one
#'     \code{<feature>.C} prepared column per feature, the importance score
#'     \code{IS}, and the logical columns \code{filter_keep},
#'     \code{rescue_keep} and \code{keep}.
#'   \item \code{pred_summary_loose}: rows with \code{keep = TRUE}.
#'   \item \code{pred_summary_strict}: complete cases of the loose set,
#'     obtained with \code{\link[stats]{na.omit}}.
#' }
#'
#' @details For each cell type the function proceeds as follows.
#' \enumerate{
#'   \item Aligns feature directions: \code{FG.Loss_norm.R} is the only
#'     lower-is-better feature and is transformed to \code{1 - value}; the other
#'     features are used unchanged.
#'   \item For each feature, computes its cutoff quantile and sets the lower
#'     tail to zero.
#'   \item Min-max rescales the surviving positive values onto
#'     \code{[feat_scale_min, 1]}, producing the \code{.C} columns.
#'   \item Forms the importance score as a weighted sum of the \code{.C}
#'     columns.
#'   \item Sets \code{filter_keep} for pairs at or above the importance-score
#'     quantile.
#'   \item If rescue is enabled, ranks pairs within each TF by a fixed priority
#'     order and rescues top candidates until the TF has
#'     \code{min_TG_per_TF} kept targets, rejecting candidates below the pooled
#'     low quantile.
#'   \item Builds the loose set from
#'     \code{keep = filter_keep | rescue_keep}, and the strict set by removing
#'     rows containing missing values.
#' }
#'
#' @note
#' \itemize{
#'   \item The function enters the \code{5 Infer GRN} folder and restores the
#'     original working directory on exit, including on error.
#'   \item \code{feat_cut_off_threshold}, \code{feat_scale_min} and
#'     \code{feat_weights} must use identical feature names.
#'   \item All quantiles use the default sample quantile method
#'     (\code{stats::quantile}, type 7).
#'   \item The feature preparation is vectorized over a single numeric matrix
#'     and the rescue step uses integer tabulation, avoiding repeated whole-table
#'     copies while preserving the same numerical output.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = readRDS(
#'   "./5 Infer GRN/interest_cell_type_pGRN.rds"
#' )
#' interest_cell_type_pGRN = filter_GRN(interest_cell_type_pGRN)
#' saveRDS(interest_cell_type_pGRN,
#'         file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
filter_GRN = function(
    interest_cell_type_pGRN,
    feat_cut_off_threshold = base::list(`FG.theta_p_norm.R` = 0.48,
                                        `FG.Loss_norm.R` = 0.63,
                                        `FG.theta_i.R` = 0.67,
                                        `F.NZER.R` = 0.24,
                                        `F.NZEA.R` = 0.32,
                                        `F.EA.R` = 0.77,
                                        `G.NZER.R` = 0.23,
                                        `G.NZEA.R` = 0.33,
                                        `G.EA.R` = 0.21),
    feat_scale_min = base::list(`FG.theta_p_norm.R` = 0.22,
                                `FG.Loss_norm.R` = 0.55,
                                `FG.theta_i.R` = 0.26,
                                `F.NZER.R` = 0.64,
                                `F.NZEA.R` = 0.75,
                                `F.EA.R` = 0.39,
                                `G.NZER.R` = 0.56,
                                `G.NZEA.R` = 0.43,
                                `G.EA.R` = 0.31),
    feat_weights = base::list(`FG.theta_p_norm.R` = 0.25,
                              `FG.Loss_norm.R` = 0.00,
                              `FG.theta_i.R` = 0.26,
                              `F.NZER.R` = 0.15,
                              `F.NZEA.R` = 0.11,
                              `F.EA.R` = 0.09,
                              `G.NZER.R` = 0.03,
                              `G.NZEA.R` = 0.07,
                              `G.EA.R` = 0.04),
    IS_quan_cut_off = 0.02,
    rescue_TF = TRUE,
    min_TG_per_TF = 1,
    TF_rescue_quantile_threshold = 0.00
) {
  t_start = Sys.time()
  message("Run: Filtering GRN. Time: ", t_start, ".")

  original_dir = getwd()
  new_folder = "5 Infer GRN"
  if (!dir.exists(new_folder)) dir.create(new_folder, recursive = TRUE)
  setwd(new_folder)
  on.exit(setwd(original_dir), add = TRUE)

  # Feature names are defined once from the cutoff list.
  feat_names = names(feat_cut_off_threshold)

  for (cell_type in names(interest_cell_type_pGRN)) {
    message("Filtering GRN of cell type ", cell_type, ".")
    GRN = interest_cell_type_pGRN[[cell_type]][["pred_summary"]]
    n = nrow(GRN)

    ### Prepare the nine features on one numeric matrix
    # Copy the raw feature columns into a single matrix to avoid repeated
    # data.frame column extraction and whole-table copies.
    X = as.matrix(GRN[, feat_names, drop = FALSE])

    for (j in seq_along(feat_names)) {
      feat = feat_names[j]

      # Loss is lower-is-better, so reverse it to share the larger-is-better
      # direction used by every other feature.
      v = if (feat == "FG.Loss_norm.R") 1 - X[, j] else X[, j]

      # Zero the lower tail at the per-feature quantile cutoff.
      cut_value = stats::quantile(v, probs = feat_cut_off_threshold[[feat]])
      v[v < cut_value] = 0  # strict less-than keeps values equal to the cutoff

      # Rescale the positive values to [feat_scale_min, 1].
      smin = feat_scale_min[[feat]]
      pos = v > 0  # positive mask after truncation
      pv = v[pos]
      rng = range(pv)
      pv = (pv - rng[1]) / (rng[2] - rng[1] + 0.0000000001) * (1 - smin) + smin
      v[pos] = pv
      X[, j] = v
    }
    colnames(X) = paste0(feat_names, ".C")  # prepared column names

    ### Compute the importance score and apply its cutoff
    # Weighted sum of the prepared columns in one matrix product.
    w = unlist(feat_weights)[feat_names]
    IS = as.numeric(X %*% w)

    # Keep pairs at or above the high importance-score quantile.
    IS_cut_off = unname(stats::quantile(IS, probs = IS_quan_cut_off))
    filter_keep = IS >= IS_cut_off

    ### Rescue underrepresented TFs
    rescue_keep = logical(n)  # all FALSE until a rescue is needed

    if (rescue_TF) {
      message("Rescuing TFs to keep at least ", min_TG_per_TF, " target each.")

      # Fixed priority order shared with the original implementation.
      ord = order(
        -GRN$FG.theta_p_norm.R, -GRN$FG.theta_i, GRN$FG.Loss_norm.R,
        -GRN$G.EA, -GRN$F.EA,
        -GRN$G.NZEA, -GRN$F.NZEA,
        -GRN$G.NZER, -GRN$F.NZER
      )

      # Compact integer TF ids and within-TF ranks derived from the order.
      tf_id = as.integer(as.factor(GRN$TF))
      tftg_rank = integer(n)
      tftg_rank[ord] = stats::ave(seq_len(n), tf_id[ord], FUN = seq_along)  # within-TF rank in priority order

      # Restrict candidates to the top min_TG_per_TF rows of each TF.
      cand_idx = which(tftg_rank <= min_TG_per_TF)
      q_thresh = stats::quantile(GRN$FG.theta_p_norm.R[cand_idx],
                                 TF_rescue_quantile_threshold, names = FALSE)
      cand_idx = cand_idx[GRN$FG.theta_p_norm.R[cand_idx] >= q_thresh]

      # Count already kept targets per TF with integer tabulation.
      per_tf_count = tabulate(tf_id[filter_keep], nbins = max(tf_id))
      candidate_keep_count = per_tf_count[tf_id[cand_idx]]
      candidate_keep_count[is.na(candidate_keep_count)] = 0  # TFs currently without a kept edge

      # Rescue only while the TF is still below the required target count.
      need_idx = cand_idx[
        candidate_keep_count < min_TG_per_TF &
          tftg_rank[cand_idx] <= (min_TG_per_TF - candidate_keep_count)
      ]
      rescue_keep[need_idx] = TRUE
    }

    keep = filter_keep | rescue_keep  # union of the cutoff and rescue rules

    ### Write results back and build the loose and strict sets
    GRN[, colnames(X)] = X
    GRN$IS = IS
    GRN$filter_keep = filter_keep
    GRN$rescue_keep = rescue_keep
    GRN$keep = keep

    interest_cell_type_pGRN[[cell_type]][["pred_summary"]] = GRN
    pred_summary_loose = GRN[GRN$keep, ]
    interest_cell_type_pGRN[[cell_type]][["pred_summary_loose"]] =
      pred_summary_loose
    interest_cell_type_pGRN[[cell_type]][["pred_summary_strict"]] =
      stats::na.omit(pred_summary_loose)
  }

  t_end = Sys.time()
  message("Finish: Filtering GRN. Running time: ", format(t_end - t_start))
  return(interest_cell_type_pGRN)
}
