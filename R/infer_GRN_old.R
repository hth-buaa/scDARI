#' @title Infer Gene Regulatory Networks with Trained Models
#'
#' @description Predicts the optimal regulatory strength and corresponding loss
#' for every candidate TF-TG pair using pre-trained scDARI models.
#' Processing can be restricted to selected cell types or branches, chunked
#' parallel computation is supported, and branch-level results are saved
#' automatically.
#'
#' @param interest_cell_type_branch_data_for_GRN_infer List. Output of
#'   \code{\link{get_data_for_inferring_GRN}}, containing branch-level pairs
#'   with group expression features.
#' @param interest_cell_type_tau Named numeric vector. Output of
#'   \code{\link{get_cell_type_tau}}, the cell-type-specific binarization
#'   thresholds.
#' @param interest_cell_type_branch_model_train Nested list. Output of
#'   \code{\link{model_train}}, containing the trained model details for each
#'   branch.
#' @param interest_cell_type_group List. Output of \code{\link{cell_grouping}},
#'   containing pseudotime interval lengths for each branch.
#' @param select_cell_type Character vector or NULL. Cell types to process;
#' NULL processes all.
#' Defaults to NULL.
#' @param select_branch Numeric vector or NULL. Branch indices to process;
#' NULL processes all.
#' Defaults to NULL.
#' @param each Positive integer. Maximum number of pairs per parallel chunk,
#' controlling peak memory.
#' Defaults to 1000000.
#' @param ncores Positive integer. Number of cores for parallel processing via
#'   \code{\link[parallel]{mclapply}}.
#' Defaults to 1; set to 1 on Windows.
#'
#' @return Nested list named \code{interest_cell_type_pGRN}, organized by cell
#' type and branch, each branch containing two data frames:
#' \itemize{
#'   \item \code{theta_p}: rows are candidate pairs and columns are trained
#'     models, values are the optimized regulatory strengths.
#'   \item \code{loss}: same rows and columns, values are the minimum
#'     prediction losses at the optimized strengths.
#' }
#' Each branch result is also saved as an RDS file in the \code{5 Infer GRN}
#' directory.
#'
#' @details For each selected branch, inference proceeds as follows:
#' \enumerate{
#'   \item Splits candidate pairs into chunks of at most \code{each} pairs.
#'   \item Processes chunks in parallel: for each pair, every trained model is
#'     applied and Brent optimization finds the strength within 0.01-0.99 that
#'     minimizes prediction loss.
#'   \item Binds chunk results into the branch-level \code{theta_p} and
#'     \code{loss} frames.
#'   \item Saves the branch result as an RDS file.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and switches to the \code{5 Infer GRN}
#'     directory, and restores the original working directory on exit.
#'   \item Reduce \code{each} if out-of-memory errors occur.
#'   \item Use \code{select_cell_type} and \code{select_branch} to process data
#'     in batches, which supports interrupted recovery.
#'   \item Parallel processing relies on \code{\link[parallel]{mclapply}};
#'     set \code{ncores = 1} on Windows.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_data_for_GRN_infer = base::readRDS(
#'   "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds"
#' )
#' interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
#' interest_cell_type_branch_model_train = base::readRDS(
#'   "./4 Build Prediction Model/interest_cell_type_branch_model_train.rds"
#' )
#' interest_cell_type_group = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' each = 1000
#' ncores = parallel::detectCores() - 1 # Linux
#' # ncores = 1 # Windows
#'
#' # Way 1: single run, when memory is sufficient.
#' interest_cell_type_pGRN = infer_GRN_old(
#'   interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
#'   interest_cell_type_tau = interest_cell_type_tau,
#'   interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
#'   interest_cell_type_group = interest_cell_type_group,
#'   select_cell_type = NULL,
#'   select_branch = NULL,
#'   each = each,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#'
#' # Way 2: branch-by-branch, more stable for large data.
#' interest_cell_type_pGRN = base::list()
#' for (cell_type in base::names(interest_cell_type_branch_data_for_GRN_infer)) {
#'   interest_cell_type_pGRN[[cell_type]] = base::list()
#'   for (n in 1:(base::length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)) {
#'     interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]] = infer_GRN_old(
#'       interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
#'       interest_cell_type_tau = interest_cell_type_tau,
#'       interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
#'       interest_cell_type_group = interest_cell_type_group,
#'       select_cell_type = cell_type,
#'       select_branch = n,
#'       each = each,
#'       ncores = ncores
#'     )[[cell_type]][[base::paste0("branch", n)]]
#'   }
#' }
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
infer_GRN_old = function(
    interest_cell_type_branch_data_for_GRN_infer,
    interest_cell_type_tau,
    interest_cell_type_branch_model_train,
    interest_cell_type_group,
    select_cell_type = NULL,
    select_branch = NULL,
    each = 1000000,
    ncores = 1
) {
  ### Start.
  t_start = base::Sys.time()
  message("Run: Inferring GRN. Time: ", t_start, ".")
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

  interest_cell_type_pGRN = base::list()
  if (base::is.null(select_cell_type)) {
    range_cell_type = base::names(interest_cell_type_branch_data_for_GRN_infer)
  } else {
    range_cell_type = select_cell_type
  }

  for (cell_type in range_cell_type) {
    interest_cell_type_pGRN[[cell_type]] = base::list()
    tau = base::as.numeric(interest_cell_type_tau[cell_type])
    if (base::is.null(select_branch)) {
      range_branch = base::seq_len(
        base::length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)
    } else {
      range_branch = select_branch
    }

    for (n in range_branch) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      state_data = interest_cell_type_branch_data_for_GRN_infer[[cell_type]][[base::paste0("branch", n)]]
      re_pred = base::rownames(state_data)
      N_pred_all = base::length(re_pred)
      message(N_pred_all, " candidate pairs will be predicted.")

      re_train = base::sort(base::names(
        interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]]))
      message(base::length(re_train), " trained models will be applied.")

      # Split candidate pairs into chunks.
      message("Partitioning data.")
      n_part = base::ceiling(N_pred_all / each)
      idx = base::max(base::ceiling(n_part / 50), 2)
      re_pred_list = base::list()
      for (s in base::seq_len(n_part)) {
        if (s %% idx == 1) {
          message("Partitioning progress: ", s, " / ", n_part, ", ",
                  base::round(s / n_part * 100), "%.")
        }
        if (s != n_part) {
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):(each * s)]
        } else {
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):N_pred_all]
        }
      }

      # Shared feature columns and parameters for the branch.
      TFE_col = base::grep("TFE_", base::colnames(state_data), value = TRUE)
      TGE_col = base::grep("TGE_", base::colnames(state_data), value = TRUE)
      Length = base::length(TFE_col)
      train_params = base::lapply(
        interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]],
        function(sublist) sublist[["best_params"]]
      )
      Tslot_K = interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]]

      # Process one chunk of candidate pairs.
      process_data_slice = function(pGRN = base::list(theta_p = NA, loss = NA),
                                    RE_PRED,
                                    ncols = base::length(re_train),
                                    STATE_DATA,
                                    TFE_col,
                                    TGE_col,
                                    Length,
                                    re_train,
                                    train_params,
                                    tau,
                                    Tslot_K) {
        for (re in RE_PRED) {
          theta_p = base::data.frame(base::matrix(
            0, nrow = 1, ncol = ncols, dimnames = base::list(re, re_train)
          ))
          loss = base::data.frame(base::matrix(
            Inf, nrow = 1, ncol = ncols, dimnames = base::list(re, re_train)
          ))

          TFE = base::as.numeric(STATE_DATA[re, TFE_col])
          TGE = base::as.numeric(STATE_DATA[re, TGE_col])

          # Apply each trained model.
          for (re_tr in re_train) {
            pa_best = train_params[[re_tr]]
            bad_loss = 1e12

            r2 = base::as.numeric(pa_best["r2"])
            a1 = base::as.numeric(pa_best["alpha1"])
            a2 = base::as.numeric(pa_best["alpha2"])
            b1 = base::as.numeric(pa_best["beta1"])
            b2 = base::as.numeric(pa_best["beta2"])
            s = base::as.numeric(pa_best["s"])
            u1i = base::as.numeric(pa_best["u1i"])
            u1j = base::as.numeric(pa_best["u1j"])
            u2i = base::as.numeric(pa_best["u2i"])
            u2j = base::as.numeric(pa_best["u2j"])
            v_tilde = base::as.numeric(pa_best["v_tilde"])

            vi = scDARI::Hill_cal(x = TFE[-Length], Dissociation_Constant = u1i,
                                  Hill_Coefficient = u2i)
            vj = scDARI::Hill_cal(x = TGE[-Length], Dissociation_Constant = u1j,
                                  Hill_Coefficient = u2j)
            vij = vi - vj
            # Clamp saturation away from the endpoints for S.
            vj_safe = base::pmin(base::pmax(vj, 1e-12), 1 - 1e-12)

            # Theta-independent boundary factor and feedback strength.
            rho_vec = TFE[-Length] / (TFE[-Length] + 0.01)
            S_vec = base::sapply(vj_safe, function(v) {
              scDARI::S_cal(s = s, v = v, v_tilde = v_tilde)
            })
            par0 = base::min(base::max(tau, 0.01 + 1e-6), 0.99 - 1e-6)

            # Optimize the regulatory strength by Brent's method.
            predict_result = stats::optim(par = par0, fn = function(theta) {
              R = scDARI::R_cal(
                theta_TF_TG = theta, tau = tau, r1 = 200, r2 = r2,
                r3 = 1, r4 = 1000, r5 = 1
              )
              if (!base::is.finite(R)) return(bad_loss)
              MSE = base::mean(base::sapply(2:Length, function(K) {
                K_1 = K - 1
                Tslot_K_K_1 = Tslot_K[K_1]
                error_f1 = TGE[K_1] + R * a1 * vij[K_1] * Tslot_K_K_1 + b1 - TGE[K]
                error_f2 = TFE[K_1] + rho_vec[K_1] * R * a2 * S_vec[K_1] *
                  Tslot_K_K_1 + b2 - TFE[K]
                error_f1 ^ 2 + error_f2 ^ 2
              }))
              if (!base::is.finite(MSE)) return(bad_loss)
              return(MSE)
            }, lower = 0.01, upper = 0.99, method = "Brent")

            theta_p[re, re_tr] = predict_result$par
            loss[re, re_tr] = predict_result$value
          }

          # Drop make.unique columns containing a literal ".".
          loss = loss[!base::grepl(".", base::names(loss), fixed = TRUE)]
          theta_p = theta_p[!base::grepl(".", base::names(theta_p), fixed = TRUE)]
          pGRN[[re]] = base::list(theta_p = theta_p, loss = loss)
        }

        # Merge pair results into chunk-level frames.
        pair_results = pGRN[3:base::length(pGRN)]
        pGRN[["theta_p"]] = base::as.data.frame(data.table::rbindlist(base::lapply(
          pair_results, function(x) x[["theta_p"]]
        )))
        base::rownames(pGRN[["theta_p"]]) = base::names(pair_results)
        pGRN[["loss"]] = base::as.data.frame(data.table::rbindlist(base::lapply(
          pair_results, function(x) x[["loss"]]
        )))
        base::rownames(pGRN[["loss"]]) = base::names(pair_results)
        return(pGRN)
      }

      # Evaluate chunks in parallel.
      prediction_list = parallel::mclapply(
        X = base::seq_len(n_part),
        FUN = function(s) {
          process_data_slice(
            RE_PRED = re_pred_list[[s]],
            STATE_DATA = state_data[re_pred_list[[s]], -(1:3), drop = FALSE],
            TFE_col = TFE_col,
            TGE_col = TGE_col,
            Length = Length,
            re_train = re_train,
            train_params = train_params,
            tau = tau,
            Tslot_K = Tslot_K
          )
        },
        mc.cores = ncores
      )
      base::gc()

      # Bind chunks into branch-level frames.
      message("Summarizing theta_p of this branch.")
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["theta_p"]] =
        base::do.call(base::rbind, base::lapply(prediction_list, function(x) x[["theta_p"]]))
      message("Summarizing loss of this branch.")
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["loss"]] =
        base::do.call(base::rbind, base::lapply(prediction_list, function(x) x[["loss"]]))

      # Save the branch result.
      message("Saving the predictive results of this branch.")
      base::saveRDS(
        interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]],
        file = base::paste0("./interest_cell_type_pGRN.",
                            base::gsub(" ", "_", cell_type), ".branch", n, ".rds")
      )
    }
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Inferring GRN. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_pGRN)
}

