#' @title Filter GRN by an Adjusted Quantile Threshold with TF Rescue
#'
#' @description Converts the continuous predictions into a final edge set with
#' an adaptive threshold rather than one fixed global cutoff. A decreasing
#' initial threshold curve is fitted from anchor points, then adjusted per
#' pair by relaxation and penalty factors built from normalized rank features.
#' Pairs whose normalized strength rank reaches the adjusted threshold are
#' kept, and a TF rescue step ensures each TF retains a minimum number of
#' targets. Loose and strict filtered sets are produced.
#'
#' @param interest_cell_type_pGRN List. Output of
#'   \code{\link{select_best_prediction}}, named by cell type and containing
#'   \code{pred_summary}.
#' @param theta_p_points Numeric vector, default \code{c(0, 0.5, 1)}.
#'   Predicted-strength anchor points for fitting the initial threshold curve.
#' @param initial_quantile_points Numeric vector, default
#'   \code{c(1.5, 0.9, 0.6)}. Initial-threshold anchor values aligned
#'   positionally to \code{theta_p_points}.
#' @param theta_p_quantile_fit_start Named list, default
#'   \code{list(a = 0.333333, b = 0.333333, c = -0.5)}. Starting values
#'   \code{a}, \code{b}, \code{c} for the nonlinear least-squares fit.
#' @param factor_balance Named list. Center value \code{fb} of each
#'   adjustment feature. Names must match those of \code{factor_range} and
#'   \code{factor_attribute}.
#' @param factor_range Named list. Strength \code{fr} of each adjustment
#'   feature.
#' @param factor_attribute Named list. Direction of each feature, either
#'   \code{"relaxation"} or \code{"penalty"}.
#' @param rescue_TF Logical, default \code{TRUE}. Whether to run the TF
#'   rescue step.
#' @param min_TG_per_TF Positive integer, default \code{1}. Minimum number
#'   of targets kept per TF after rescue.
#' @param TF_rescue_quantile_threshold Numeric in \code{[0, 1]}, default
#'   \code{0.1}. Pooled low quantile of \code{FG.theta_p_norm.R}; rescue
#'   candidates below it are rejected.
#'
#' @return The input list with the following per cell type:
#' \itemize{
#'   \item \code{pred_summary}: the original frame augmented with the initial
#'     threshold \code{IT}, one \code{<feature>-F} column per feature, the
#'     adjusted threshold \code{AT}, and the logical columns
#'     \code{filter_keep}, \code{rescue_keep} and \code{keep}.
#'   \item \code{pred_summary_loose}: rows with \code{keep = TRUE}.
#'   \item \code{pred_summary_strict}: complete cases of the loose set,
#'     obtained with \code{\link[stats]{na.omit}}.
#' }
#'
#' @details For each cell type the function proceeds as follows.
#' \enumerate{
#'   \item Fits the initial threshold curve
#'     \eqn{\mathrm{IT} = 1/(a\,\theta_p + b) + c} from the anchor points
#'     with \eqn{a > 0}, and evaluates \code{IT} from \code{FG.theta_p}.
#'   \item For each feature, computes its factor as
#'     \eqn{1 + f_r(f_b - v)} for relaxation and
#'     \eqn{1 + f_r(v - f_b)} for penalty.
#'   \item Forms the adjusted threshold as the product of \code{IT} and all
#'     factors on the log scale, clamps non-positive values, and clamps
#'     \code{AT} to \code{[0, 1]}.
#'   \item Sets \code{filter_keep} when
#'     \code{FG.theta_p_norm.R >= AT}.
#'   \item If rescue is enabled, ranks pairs within each TF by a priority
#'     order and rescues top candidates until the TF has
#'     \code{min_TG_per_TF} kept targets, rejecting candidates below the
#'     pooled low quantile.
#'   \item Builds the loose set from \code{keep = filter_keep | rescue_keep},
#'     and the strict set by removing rows containing missing values.
#' }
#'
#' @note
#' \itemize{
#'   \item The function enters the \code{5 Infer GRN} folder and restores the
#'     original working directory on exit, including on error.
#'   \item \code{factor_balance}, \code{factor_range} and
#'     \code{factor_attribute} must use identical names.
#'   \item \code{IT} may exceed 1 for weak predicted strengths; after
#'     clamping those pairs face \code{AT = 1} and are filtered out.
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
    theta_p_points = c(0, 0.5, 1),
    initial_quantile_points = c(1.5, 0.9, 0.6),
    theta_p_quantile_fit_start = list(a = 0.333333, b = 0.333333, c = -0.5),
    factor_balance = list(
      FG.Loss_norm.R = 0.2, FG.theta_i.R = 0.2,
      F.NZER.R = 0.2, F.NZEA.R = 0.2, F.EA.R = 0.2,
      G.NZER.R = 0.2, G.NZEA.R = 0.2, G.EA.R = 0.2
    ),
    factor_range = list(
      FG.Loss_norm.R = 0.25, FG.theta_i.R = 0.25,
      F.NZER.R = 0.05, F.NZEA.R = 0.05, F.EA.R = 0.05,
      G.NZER.R = 0.15, G.NZEA.R = 0.15, G.EA.R = 0.15
    ),
    factor_attribute = list(
      FG.Loss_norm.R = "penalty", FG.theta_i.R = "relaxation",
      F.NZER.R = "relaxation", F.NZEA.R = "relaxation", F.EA.R = "relaxation",
      G.NZER.R = "relaxation", G.NZEA.R = "relaxation", G.EA.R = "relaxation"
    ),
    rescue_TF = TRUE,
    min_TG_per_TF = 1,
    TF_rescue_quantile_threshold = 0.1
) {
  t_start = Sys.time()
  message("Run: Filtering GRN. Time: ", t_start, ".")

  original_dir = getwd()
  new_folder = "5 Infer GRN"
  if (!dir.exists(new_folder)) dir.create(new_folder, recursive = TRUE)
  setwd(new_folder)
  on.exit(setwd(original_dir), add = TRUE)

  ### Fit the inverse threshold curve IT = 1/(a theta + b) + c
  fit = minpack.lm::nlsLM(
    initial_quantile_points ~ 1 / (theta_p_points * coef_a + coef_b) + coef_c,
    start = list(coef_a = theta_p_quantile_fit_start$a,
                 coef_b = theta_p_quantile_fit_start$b,
                 coef_c = theta_p_quantile_fit_start$c),
    lower = c(coef_a = 1e-9, coef_b = -Inf, coef_c = -Inf),
    upper = c(coef_a = Inf, coef_b = Inf, coef_c = Inf)
  )
  coef_a = as.numeric(stats::coef(fit)["coef_a"])
  coef_b = as.numeric(stats::coef(fit)["coef_b"])
  coef_c = as.numeric(stats::coef(fit)["coef_c"])
  rm(fit)
  message("Fitted coefficients: a=", round(coef_a, 6),
          ", b=", round(coef_b, 6), ", c=", round(coef_c, 6), ".")

  for (cell_type in names(interest_cell_type_pGRN)) {
    message("Filtering GRN of cell type ", cell_type, ".")
    GRN = interest_cell_type_pGRN[[cell_type]][["pred_summary"]]

    GRN$IT = 1 / (GRN$FG.theta_p * coef_a + coef_b) + coef_c

    factor_features = names(factor_balance)
    for (ff in factor_features) {
      fb = factor_balance[[ff]]; fr = factor_range[[ff]]
      fa = (-1)^as.numeric(factor_attribute[[ff]] == "penalty")
      GRN[, paste0(ff, "-F")] = (1 + fa * fb * fr) - fa * GRN[, ff] * fr
    }

    at_cols = c("IT", paste0(factor_features, "-F"))
    at_matrix = as.matrix(GRN[, at_cols])
    at_matrix[is.na(at_matrix) | at_matrix <= 0] = .Machine$double.eps
    GRN$AT = exp(rowSums(log(at_matrix)))
    GRN$AT = pmin(pmax(GRN$AT, 0), 1)

    GRN$filter_keep = GRN$FG.theta_p_norm.R >= GRN$AT
    GRN$filter_keep[is.na(GRN$filter_keep)] = FALSE

    GRN$rescue_keep = FALSE
    if (rescue_TF) {
      message("Rescuing TFs to keep at least ", min_TG_per_TF, " target each.")
      ord = order(
        -GRN$FG.theta_p_norm.R, -GRN$FG.theta_i.R, GRN$FG.Loss_norm.R,
        -GRN$G.EA.R, -GRN$F.EA.R,
        -GRN$G.NZEA.R, -GRN$F.NZEA.R,
        -GRN$G.NZER.R, -GRN$F.NZER.R
      )
      GRN_ordered = GRN[ord, ]
      GRN_ordered$TFTG_rank = stats::ave(
        seq_len(nrow(GRN_ordered)), GRN_ordered$TF, FUN = seq_along)

      candidates = GRN_ordered[GRN_ordered$TFTG_rank <= min_TG_per_TF, ]
      q_thresh = stats::quantile(candidates$FG.theta_p_norm.R,
                                 TF_rescue_quantile_threshold, names = FALSE)
      candidates = candidates[candidates$FG.theta_p_norm.R >= q_thresh, ]

      keep_count = tapply(GRN$filter_keep, GRN$TF, sum)
      candidate_keep_count = keep_count[candidates$TF]
      candidate_keep_count[is.na(candidate_keep_count)] = 0
      need_rescue = candidate_keep_count < min_TG_per_TF &
        candidates$TFTG_rank <= (min_TG_per_TF - candidate_keep_count)
      GRN[rownames(candidates)[need_rescue], "rescue_keep"] = TRUE
    }
    GRN$keep = GRN$filter_keep | GRN$rescue_keep

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
