#' @title Filter GRN by an Adjusted Quantile Threshold with TF Rescue
#'
#' @description Filters inferred GRN results with an adaptive threshold.
#' An initial threshold curve decreasing with predicted strength is fitted from anchor points,
#' then adjusted per pair by relaxation and penalty factors built from normalized rank features.
#' Pairs whose normalized strength rank reaches the adjusted threshold are kept,
#' and a TF rescue step ensures each TF retains a minimum number of targets.
#' Loose and strict filtered sets are produced.
#'
#' @param interest_cell_type_pGRN List. Output of \code{\link{select_best_prediction}}.
#' @param theta_p_points Numeric vector.
#' Predicted-strength anchor points for fitting the initial threshold curve.
#' Defaults to \code{c(0, 0.5, 1)}.
#' @param initial_quantile_points Numeric vector.
#' Initial-threshold anchor values aligned to \code{theta_p_points}.
#' Defaults to \code{c(1.5, 0.9, 0.6)}.
#' @param theta_p_quantile_fit_start Named list.
#' Starting values \code{a}, \code{b}, \code{c} for the nonlinear fit.
#' Defaults to \code{list(a = 0.333333, b = 0.333333, c = -0.5)}.
#' @param factor_balance Named list.
#' Center value \code{fb} of each adjustment feature.
#' Names must match those of \code{factor_range} and \code{factor_attribute}.
#' @param factor_range Named list. Strength \code{fr} of each adjustment feature.
#' @param factor_attribute Named list. Direction of each feature,
#' either \code{"relaxation"} or \code{"penalty"}.
#' @param rescue_TF Logical. Whether to run the TF rescue step.
#' Defaults to TRUE.
#' @param min_TG_per_TF Positive integer.
#' Minimum number of targets kept per TF after rescue.
#' Defaults to 1.
#' @param TF_rescue_quantile_threshold Numeric between 0 and 1.
#' Pooled low quantile of \code{FG.theta_p_norm.R};
#' rescue candidates below it are rejected.
#' Defaults to 0.1.
#'
#' @return The input list with the following per cell type:
#' \itemize{
#'   \item \code{pred_summary}:
#'   the original frame augmented with the initial threshold \code{IT},
#'   one \code{<feature>-F} column per feature, the adjusted threshold \code{AT},
#'   and the logical columns \code{filter_keep}, \code{rescue_keep} and \code{keep}.
#'   \item \code{pred_summary_loose}: rows with \code{keep = TRUE}.
#'   \item \code{pred_summary_strict}:
#'   complete cases of the loose set, obtained with \code{\link[stats]{na.omit}}.
#' }
#'
#' @details For each cell type, the function proceeds as follows:
#' \enumerate{
#'   \item Fits the initial threshold curve \eqn{\mathrm{IT} = 1/(a\,\theta_p + b) + c} from the anchor points,
#'   with \eqn{a > 0}, and evaluates \code{IT} from \code{FG.theta_p}.
#'   \item For each feature,
#'   computes its factor as \eqn{1 + f_r(f_b - v)} for relaxation and \eqn{1 + f_r(v - f_b)} for penalty.
#'   \item Forms the adjusted threshold as the product of \code{IT} and all factors,
#'   computed on the log scale with non-positive values clamped, and clamps \code{AT} to \[0, 1\].
#'   \item Sets \code{filter_keep} when \code{FG.theta_p_norm.R >= AT}.
#'   \item If rescue is enabled, ranks pairs within each TF by a priority order,
#'   and rescues top candidates until the TF has \code{min_TG_per_TF} kept targets,
#'   rejecting candidates below the pooled low quantile.
#'   \item Builds the loose set from \code{keep = filter_keep | rescue_keep},
#'   and the strict set by removing rows containing missing values.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and switches to the \code{5 Infer GRN} directory,
#'   and restores the original working directory on exit.
#'   \item \code{factor_balance}, \code{factor_range} and \code{factor_attribute} must use identical names.
#'   \item \code{IT} may exceed 1 for weak predicted strengths;
#'   after clamping those pairs face \code{AT = 1} and are filtered out.
#'   \item The loose set is suitable for exploratory network construction,
#'   and the strict complete-case set for high-confidence downstream analysis.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' interest_cell_type_pGRN = filter_GRN(
#'   interest_cell_type_pGRN = interest_cell_type_pGRN
#' )
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
filter_GRN = function(
    interest_cell_type_pGRN,
    theta_p_points = base::c(0, 0.5, 1),
    initial_quantile_points = base::c(1.5, 0.9, 0.6),
    theta_p_quantile_fit_start = base::list(a = 0.333333, b = 0.333333, c = -0.5),
    factor_balance = base::list(
      FG.Loss_norm.R = 0.2, FG.theta_i.R = 0.2,
      F.NZER.R = 0.2, F.NZEA.R = 0.2, F.EA.R = 0.2,
      G.NZER.R = 0.2, G.NZEA.R = 0.2, G.EA.R = 0.2
    ),
    factor_range = base::list(
      FG.Loss_norm.R = 0.25, FG.theta_i.R = 0.25,
      F.NZER.R = 0.05, F.NZEA.R = 0.05, F.EA.R = 0.05,
      G.NZER.R = 0.15, G.NZEA.R = 0.15, G.EA.R = 0.15
    ),
    factor_attribute = base::list(
      FG.Loss_norm.R = "penalty", FG.theta_i.R = "relaxation",
      F.NZER.R = "relaxation", F.NZEA.R = "relaxation", F.EA.R = "relaxation",
      G.NZER.R = "relaxation", G.NZEA.R = "relaxation", G.EA.R = "relaxation"
    ),
    rescue_TF = TRUE,
    min_TG_per_TF = 1,
    TF_rescue_quantile_threshold = 0.1
) {
  ### Start.
  t_start = base::Sys.time()
  message("Run: Filtering GRN. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "5 Infer GRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  # Restore the working directory on normal exit or error.
  base::on.exit({
    base::setwd(original_dir)
    message("The working directory has been restored to: ", base::getwd(), ".")
  }, add = TRUE)

  ### Fit the inverse initial-threshold curve, IT = 1/(a theta_p + b) + c.
  message("Fitting the initial threshold curve.")
  fit = minpack.lm::nlsLM(
    initial_quantile_points ~ 1 / (theta_p_points * coef_a + coef_b) + coef_c,
    start = base::list(coef_a = theta_p_quantile_fit_start$a,
                       coef_b = theta_p_quantile_fit_start$b,
                       coef_c = theta_p_quantile_fit_start$c),
    lower = base::c(coef_a = 1e-9, coef_b = -Inf, coef_c = -Inf),
    upper = base::c(coef_a = Inf, coef_b = Inf, coef_c = Inf)
  )
  coef_a = base::as.numeric(stats::coef(fit)["coef_a"])
  coef_b = base::as.numeric(stats::coef(fit)["coef_b"])
  coef_c = base::as.numeric(stats::coef(fit)["coef_c"])
  base::rm(fit)
  message("Fitted coefficients: a=", base::round(coef_a, 6),
          ", b=", base::round(coef_b, 6),
          ", c=", base::round(coef_c, 6), ".")

  ### Process.
  for (cell_type in base::names(interest_cell_type_pGRN)) {
    message("Filtering GRN of cell type ", cell_type, ".")
    GRN = interest_cell_type_pGRN[[cell_type]][["pred_summary"]]

    # Initial threshold from predicted strength.
    GRN$IT = 1 / (GRN$FG.theta_p * coef_a + coef_b) + coef_c

    # Feature-specific adjustment factors.
    factor_features = base::names(factor_balance)
    for (ff in factor_features) {
      message("Computing ", factor_attribute[[ff]], " factor for ", ff, ".")
      fb = factor_balance[[ff]]
      fr = factor_range[[ff]]
      # Relaxation uses fa = 1; penalty uses fa = -1.
      fa = (-1)^base::as.numeric(factor_attribute[[ff]] == "penalty")
      GRN[, base::paste0(ff, "-F")] =
        (1 + fa * fb * fr) - fa * GRN[, ff] * fr
    }

    # Adjusted threshold as the product of IT and factors on the log scale.
    at_cols = base::c("IT", base::paste0(factor_features, "-F"))
    at_matrix = base::as.matrix(GRN[, at_cols])
    at_matrix[base::is.na(at_matrix) | at_matrix <= 0] = .Machine$double.eps
    GRN$AT = base::exp(base::rowSums(base::log(at_matrix)))
    GRN$AT = base::pmin(base::pmax(GRN$AT, 0), 1)

    # Primary threshold filtering; missing comparisons are treated as FALSE.
    GRN$filter_keep = GRN$FG.theta_p_norm.R >= GRN$AT
    GRN$filter_keep[base::is.na(GRN$filter_keep)] = FALSE

    # TF rescue: keep a minimum number of targets per TF.
    GRN$rescue_keep = FALSE
    if (rescue_TF) {
      message("Rescuing TFs to keep at least ", min_TG_per_TF, " target each.")
      # Priority order used to rank candidates within each TF.
      ord = base::order(
        -GRN$FG.theta_p_norm.R, -GRN$FG.theta_i.R, GRN$FG.Loss_norm.R,
        -GRN$G.EA.R, -GRN$F.EA.R,
        -GRN$G.NZEA.R, -GRN$F.NZEA.R,
        -GRN$G.NZER.R, -GRN$F.NZER.R
      )
      GRN_ordered = GRN[ord, ]
      GRN_ordered$TFTG_rank = stats::ave(
        base::seq_len(base::nrow(GRN_ordered)),
        GRN_ordered$TF, FUN = base::seq_along
      )

      # Top candidates per TF, restricted by a pooled low-quantile quality cut.
      candidates = GRN_ordered[GRN_ordered$TFTG_rank <= min_TG_per_TF, ]
      q_thresh = stats::quantile(candidates$FG.theta_p_norm.R,
                                 TF_rescue_quantile_threshold, names = FALSE)
      candidates = candidates[candidates$FG.theta_p_norm.R >= q_thresh, ]

      # Rescue only as many top candidates as needed to fill the gap.
      keep_count = base::tapply(GRN$filter_keep, GRN$TF, sum)
      candidate_keep_count = keep_count[candidates$TF]
      candidate_keep_count[base::is.na(candidate_keep_count)] = 0
      need_rescue = candidate_keep_count < min_TG_per_TF &
        candidates$TFTG_rank <= (min_TG_per_TF - candidate_keep_count)
      GRN[base::rownames(candidates)[need_rescue], "rescue_keep"] = TRUE
    }
    GRN$keep = GRN$filter_keep | GRN$rescue_keep

    # Store the augmented summary.
    interest_cell_type_pGRN[[cell_type]][["pred_summary"]] = GRN

    # Loose set: pairs kept by threshold or rescue.
    pred_summary_loose = GRN[GRN$keep, ]
    interest_cell_type_pGRN[[cell_type]][["pred_summary_loose"]] = pred_summary_loose

    # Strict set: complete cases of the loose set.
    interest_cell_type_pGRN[[cell_type]][["pred_summary_strict"]] =
      stats::na.omit(pred_summary_loose)
  }

  ### End.
  t_end = base::Sys.time()
  message("Finish: Filtering GRN. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start), ".")
  return(interest_cell_type_pGRN)
}
