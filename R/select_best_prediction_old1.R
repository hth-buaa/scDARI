#' @title Select Best GRN Predictions with Multi-Level Filtering
#'
#' @description Post-processes inferred GRN results to select high-confidence predictions.
#' For each candidate pair within a branch, the optimal prediction is chosen
#' from multiple trained models by a four-level hierarchical strategy combining
#' quantile thresholds and absolute bounds.
#' Branch-level results are then integrated into a cell-type summary using
#' cell-count weighted averaging.
#'
#' @param interest_cell_type_pGRN List. Output of \code{\link{infer_GRN}}.
#' @param interest_cell_type_tau Named numeric vector.
#' Output of \code{\link{get_cell_type_tau}}, the cell-type-specific binarization thresholds.
#' @param interest_cell_type_genes_pseudotime_info List.
#' Output of \code{\link{get_genes_pseudotime_info}}, containing branch cell counts for weighted integration.
#' @param fit_loss_quantile_threshold Numeric between 0 and 1.
#' Quantile threshold for loss-based filtering. Defaults to 0.2.
#' @param theta_p_quantile_threshold Numeric between 0 and 1.
#' Quantile threshold for regulatory-strength filtering. Defaults to 0.8.
#' @param fit_loss_upper Positive numeric. Absolute upper bound of acceptable loss;
#' 999 applies no filtering by default and is intended for debugging,
#' while quantile thresholds are recommended for routine use. Defaults to 999.
#' @param theta_p_lower Numeric between 0 and 1.
#' Absolute lower bound of acceptable regulatory strength;
#' 0 applies no filtering by default and is intended for debugging. Defaults to 0.
#' @param each Positive integer. Maximum number of pairs per parallel chunk.
#' Defaults to 1000000.
#' @param ncores Positive integer.
#' Number of cores for parallel processing via \code{\link[parallel]{mclapply}}.
#' Defaults to 1; set to 1 on Windows.
#'
#' @return The input list \code{interest_cell_type_pGRN} with added elements,
#' organized by cell type:
#' \itemize{
#'   \item For each branch, a \code{best_result} frame with the selected
#'     regulatory strength \code{theta_p}, its fitting loss and the selected
#'     model name for each pair.
#'   \item A cell-type-level \code{pred_summary} frame integrating branch
#'     values and cell-count weighted averages of \code{theta_p} and
#'     \code{fit_loss}.
#' }
#'
#' @details For each cell type, the function proceeds as follows:
#' \enumerate{
#'   \item Computes branch thresholds by combining quantile statistics with
#'     absolute bounds, giving two loss thresholds, tight and standard, and one
#'     minimum strength threshold.
#'   \item Splits pairs into chunks for memory control and parallel processing.
#'   \item For each pair, applies four-level selection with sequential fallback:
#'     \enumerate{
#'       \item Below the tight loss threshold, take the maximum strength.
#'       \item Otherwise below the standard loss threshold, take the maximum
#'         strength.
#'       \item Otherwise above the minimum strength, take the minimum loss.
#'       \item Otherwise take the minimum loss directly.
#'     }
#'   \item Removes boundary-invalid strengths below 0.011 and infinite losses.
#'   \item Integrates branch results with cell-count weighted averages to form
#'     the cell-type summary.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and switches to the \code{5 Infer GRN}
#'     directory, and restores the original working directory on exit.
#'   \item Reduce \code{each} if out-of-memory errors occur.
#'   \item New elements are appended to a copy of the input list, which is
#'     returned.
#'   \item A fixed threshold of 0.011 removes boundary-invalid strengths after
#'     selection.
#'   \item Parallel processing relies on \code{\link[parallel]{mclapply}};
#'     set \code{ncores = 1} on Windows.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
#' interest_cell_type_genes_pseudotime_info = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#' ncores = parallel::detectCores() - 1 # Linux
#' # ncores = 1 # Windows
#' interest_cell_type_pGRN = select_best_prediction_old1(
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   interest_cell_type_tau = interest_cell_type_tau,
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   fit_loss_quantile_threshold = 0.2,
#'   theta_p_quantile_threshold = 0.8,
#'   fit_loss_upper = 999,
#'   theta_p_lower = 0,
#'   each = 1000,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
select_best_prediction_old1 = function(
    interest_cell_type_pGRN,
    interest_cell_type_tau,
    interest_cell_type_genes_pseudotime_info,
    fit_loss_quantile_threshold = 0.2,
    theta_p_quantile_threshold = 0.8,
    fit_loss_upper = 999,
    theta_p_lower = 0,
    each = 1000000,
    ncores = 1
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

  ### Process.
  for (cell_type in base::names(interest_cell_type_pGRN)) {
    interest_cell_type_pGRN[[cell_type]][["pred_summary"]] = base::list()
    branch_index = base::seq_along(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]])

    for (n in branch_index) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      branch_result = interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]]
      re_pred = base::rownames(branch_result[["theta_p"]])
      N_pred_all = base::length(re_pred)
      message("Branch ", n, " has ", N_pred_all, " candidate pairs.")

      # Split pairs into chunks.
      message("Partitioning data.")
      n_part = base::ceiling(N_pred_all / each)
      re_pred_list = base::list()
      for (s in base::seq_len(n_part)) {
        if (s != n_part) {
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):(each * s)]
        } else {
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):N_pred_all]
        }
      }

      # Per-chunk inputs.
      pGRN_list = base::list()
      for (s in base::seq_len(n_part)) {
        chunk_rows = re_pred_list[[s]]
        pGRN_list[[s]] = base::list(
          theta_p = branch_result[["theta_p"]][chunk_rows, , drop = FALSE],
          loss = branch_result[["loss"]][chunk_rows, , drop = FALSE]
        )
      }

      # Loss thresholds: tight quantile and standard quantile, bounded above.
      message("Computing branch filtering thresholds.")
      all_loss = base::unlist(branch_result[["loss"]])
      loss_max_s1 = base::min(
        stats::quantile(all_loss, probs = fit_loss_quantile_threshold / 2, na.rm = TRUE),
        fit_loss_upper
      )
      loss_max_s2 = base::min(
        stats::quantile(all_loss, probs = fit_loss_quantile_threshold, na.rm = TRUE),
        fit_loss_upper
      )
      # Minimum strength: high quantile, bounded below by tau and the absolute bound.
      all_theta = base::unlist(branch_result[["theta_p"]])
      theta_p_min = base::max(
        stats::quantile(all_theta, probs = theta_p_quantile_threshold, na.rm = TRUE),
        interest_cell_type_tau[cell_type],
        theta_p_lower
      )

      # Select the best prediction for every pair in one chunk.
      process_data_slice = function(pGRN,
                                    loss_max_s1,
                                    loss_max_s2,
                                    theta_p_min) {
        theta_frame = pGRN[["theta_p"]]
        loss_frame = pGRN[["loss"]]
        for (re in base::rownames(theta_frame)) {
          theta_vec = base::unlist(theta_frame[re, , drop = FALSE])
          loss_vec = base::unlist(loss_frame[re, , drop = FALSE])
          best_result = base::data.frame(base::matrix(
            NA, nrow = 1, ncol = 3,
            dimnames = base::list(re, c("theta_p", "loss", "fit"))
          ))

          # Four-level selection with sequential fallback.
          eligible1 = base::which(loss_vec <= loss_max_s1)
          eligible2 = base::which(loss_vec <= loss_max_s2)
          eligible3 = base::which(theta_vec >= theta_p_min)
          if (base::length(eligible1) > 0) {
            pick = eligible1[base::which.max(theta_vec[eligible1])]
          } else if (base::length(eligible2) > 0) {
            pick = eligible2[base::which.max(theta_vec[eligible2])]
          } else if (base::length(eligible3) > 0) {
            pick = eligible3[base::which.min(loss_vec[eligible3])]
          } else {
            pick = base::which.min(loss_vec)
          }

          best_result[re, "theta_p"] = theta_vec[pick]
          best_result[re, "loss"] = loss_vec[pick]
          best_result[re, "fit"] = base::names(loss_vec)[pick]
          pGRN[[re]] = base::list(best_result = best_result)
        }

        # Merge pair results into a chunk-level frame.
        pair_results = pGRN[3:base::length(pGRN)]
        pGRN[["best_result"]] = base::as.data.frame(data.table::rbindlist(base::lapply(
          pair_results, function(x) x[["best_result"]]
        )))
        base::rownames(pGRN[["best_result"]]) = base::names(pair_results)

        # Remove boundary-invalid strengths and infinite losses.
        pGRN[["best_result"]] = pGRN[["best_result"]][
          pGRN[["best_result"]]$theta_p >= 0.011 &
            !base::is.infinite(pGRN[["best_result"]]$loss), , drop = FALSE]

        return(pGRN[["best_result"]])
      }

      prediction_list = parallel::mclapply(
        X = base::seq_len(n_part),
        FUN = function(s) {
          process_data_slice(pGRN = pGRN_list[[s]],
                             loss_max_s1 = loss_max_s1,
                             loss_max_s2 = loss_max_s2,
                             theta_p_min = theta_p_min)
        },
        mc.cores = ncores
      )
      base::gc()

      # Bind chunks and suffix the columns with the branch index.
      message("Summarizing best_result of this branch.")
      branch_best = base::do.call(base::rbind, prediction_list)
      base::colnames(branch_best) = base::paste0(base::colnames(branch_best), "_branch", n)
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]] = branch_best
    }

    # Build the integrated cell-type summary.
    message("Integrating branch results.")
    branch_blocks = interest_cell_type_pGRN[[cell_type]][base::paste0("branch", branch_index)]
    rows = base::sort(base::unique(base::Reduce(base::c, base::lapply(
      branch_blocks, function(x) base::rownames(x[["best_result"]])
    ))))
    cols = base::c(
      "TF", "TG",
      base::Reduce(base::c, base::lapply(branch_blocks, function(x) base::colnames(x[["best_result"]])))
    )
    pred_summary = base::data.frame(base::matrix(
      NA, nrow = base::length(rows), ncol = base::length(cols),
      dimnames = base::list(rows, cols)
    ))
    pred_summary[, c("TF", "TG")] = base::do.call(
      base::rbind,
      base::strsplit(base::rownames(pred_summary), "_to_")
    )

    for (n in branch_index) {
      branch_best = branch_blocks[[n]][["best_result"]]
      pred_summary[base::rownames(branch_best), base::colnames(branch_best)] = branch_best
    }

    n_branch = base::length(branch_index)
    weights = base::unlist(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]])
    if (n_branch > 1) {
      pred_summary$theta_p = base::apply(
        pred_summary[, base::grep("^theta_p_branch([0-9]{1,2})$", base::names(pred_summary), value = TRUE)],
        1, function(x) {
          valid = base::which(!base::is.na(x))
          if (base::length(valid) == 0) 0
          else base::sum(x[valid] * weights[valid]) / base::sum(weights[valid])
        })
      pred_summary$fit_loss = base::apply(
        pred_summary[, base::grep("^loss_branch([0-9]{1,2})$", base::names(pred_summary), value = TRUE)],
        1, function(x) {
          valid = base::which(!base::is.na(x))
          if (base::length(valid) == 0) Inf
          else base::sum(x[valid] * weights[valid]) / base::sum(weights[valid])
        })
    } else if (n_branch == 1) {
      pred_summary$theta_p = pred_summary$theta_p_branch1
      pred_summary$fit_loss = pred_summary$loss_branch1
    } else {
      base::stop("The number of branches is less than 1.")
    }

    interest_cell_type_pGRN[[cell_type]][["pred_summary"]] = pred_summary
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Filtering GRN. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_pGRN)
}
