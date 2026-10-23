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
#' @param select_branch Numeric vector or NULL. Branch indices to process.
#' NULL processes all.
#' Defaults to NULL.
#' @param each Positive integer. Maximum number of pairs per parallel chunk,
#' controlling peak memory.
#' Defaults to 1000000.
#' @param ncores Positive integer. Number of cores for parallel processing via
#'   \code{\link[parallel]{mclapply}}.
#' Defaults to 1; set to 1 on Windows.
#' @param save_check_point logical. Defaults to \code{TRUE} for saving the
#'   check point after inferring per branch.
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
#' @details For a fixed pair and model, the prediction loss is a quadratic
#' in the effect efficiency \eqn{R}:
#' \deqn{\mathrm{MSE}(R) = A R^2 + B R + C,}
#' built from the consecutive-group residuals. Its unconstrained minimum is
#' \eqn{R^* = -B/(2A)}. Because \eqn{R(\theta)} is a monotonic gating function
#' whose shape depends only on \eqn{\tau} and the model's \eqn{r_2}, a single
#' monotonic inverse (shared by every pair) maps \eqn{R^*} back to the optimal
#' strength \eqn{\theta}. All coefficients are computed with matrix algebra
#' over the pairs, so the per-pair Brent optimization is eliminated.
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
#' interest_cell_type_pGRN = infer_GRN(
#'   interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
#'   interest_cell_type_tau = interest_cell_type_tau,
#'   interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
#'   interest_cell_type_group = interest_cell_type_group,
#'   select_cell_type = NULL,
#'   select_branch = NULL,
#'   each = each,
#'   ncores = ncores,
#'   save_check_point = TRUE
#' )
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#'
#' # Way 2: branch-by-branch, more stable for large data.
#' interest_cell_type_pGRN = base::list()
#' for (cell_type in base::names(interest_cell_type_branch_data_for_GRN_infer)) {
#'   interest_cell_type_pGRN[[cell_type]] = base::list()
#'   for (n in 1:(base::length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)) {
#'     interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]] = infer_GRN(
#'       interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
#'       interest_cell_type_tau = interest_cell_type_tau,
#'       interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
#'       interest_cell_type_group = interest_cell_type_group,
#'       select_cell_type = cell_type,
#'       select_branch = n,
#'       each = each,
#'       ncores = ncores,
#'       save_check_point = TRUE
#'     )[[cell_type]][[base::paste0("branch", n)]]
#'   }
#' }
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
infer_GRN = function(
    interest_cell_type_branch_data_for_GRN_infer,
    interest_cell_type_tau,
    interest_cell_type_branch_model_train,
    interest_cell_type_group,
    select_cell_type = NULL,
    select_branch = NULL,
    each = 1000000,
    ncores = 1,
    save_check_point = TRUE
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
  # Restore the original working directory even if the function is interrupted.
  base::on.exit(base::setwd(original_dir), add = TRUE)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  # Vectorized Hill saturation preserving matrix dimensions. The exported
  # Hill_cal coerces its input with as.numeric, which flattens a matrix; this
  # wrapper keeps the (pairs x terms) shape. Mathematically identical.
  Hill_sat = function(x, K, n) {
    base::storage.mode(x) = "double"
    xn = x^n
    xn / (K + xn)
  }

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

      branch_models = interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]]
      re_train = base::sort(base::names(branch_models))
      n_model_all = base::length(re_train)
      # Output columns keep unique names; if duplicate model names exist, the
      # original data.frame(matrix()) makes them unique ("a","a.1") and then
      # drops the ".1" columns, while per-name assignment makes the last
      # duplicate win. Replicate that: unique columns, later duplicates win.
      re_train_cols = base::unique(re_train)
      n_model = base::length(re_train_cols)
      message(n_model_all, " trained models will be applied.")

      # Shared feature columns and interval lengths for the branch.
      TFE_col = base::grep("TFE_", base::colnames(state_data), value = TRUE)
      TGE_col = base::grep("TGE_", base::colnames(state_data), value = TRUE)
      n_group = base::length(TFE_col)
      Tslot_K = base::as.numeric(
        interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]]
      )
      # The MSE uses transitions between consecutive groups: K_1 = 1..(n_group-1).
      term_idx = base::seq_len(n_group - 1)
      Tslot_term = Tslot_K[term_idx]

      # Extract the branch feature matrices once (rows = candidate pairs).
      TFE_mat = base::as.matrix(state_data[, TFE_col, drop = FALSE])
      TGE_mat = base::as.matrix(state_data[, TGE_col, drop = FALSE])
      storage.mode(TFE_mat) = "double"
      storage.mode(TGE_mat) = "double"
      TFE_k = TFE_mat[, term_idx, drop = FALSE]
      TFE_next = TFE_mat[, term_idx + 1L, drop = FALSE]
      TGE_k = TGE_mat[, term_idx, drop = FALSE]
      TGE_next = TGE_mat[, term_idx + 1L, drop = FALSE]
      # The boundary factor rho depends only on the pair, not on the model.
      rho_mat = TFE_k / (TFE_k + 0.01)

      # Per-model scalars and the shared monotonic inverse lookup table.
      # R(theta) depends only on tau (fixed here) and the model's r2, so the
      # theta -> R mapping is shared by every candidate pair.
      model_meta = stats::setNames(base::lapply(re_train, function(re_tr) {
        pa_best = branch_models[[re_tr]][["best_params"]]
        getp = function(nm) base::as.numeric(pa_best[nm])
        r2 = getp("r2")
        theta_grid = base::seq(0.01, 0.99, length.out = 2000)
        R_grid = R_cal(
          theta_TF_TG = theta_grid, tau = tau, r1 = 200, r2 = r2,
          r3 = 1, r4 = 1000, r5 = 1
        )
        base::list(
          r2 = r2,
          a1 = getp("alpha1"), a2 = getp("alpha2"),
          b1 = getp("beta1"), b2 = getp("beta2"),
          s = getp("s"),
          u1i = getp("u1i"), u2i = getp("u2i"),
          u1j = getp("u1j"), u2j = getp("u2j"),
          v_tilde = getp("v_tilde"),
          R_lo = R_grid[1], R_hi = R_grid[base::length(R_grid)],
          theta_grid = theta_grid, R_grid = R_grid
        )
      }), re_train)

      # Split candidate pairs into chunks of at most `each`.
      each = base::as.integer(each)
      if (!(base::is.finite(each) && each > 0)) {
        base::stop("'each' must be a positive integer.")
      }
      n_part = base::ceiling(N_pred_all / each)
      chunk_idx = base::seq_len(n_part)
      chunk_bounds = base::lapply(chunk_idx, function(s) {
        lo = each * (s - 1) + 1L
        hi = base::min(each * s, N_pred_all)
        base::seq.int(lo, hi)
      })

      # Process one chunk of candidate pairs. All algebra is vectorized over
      # the pairs in the chunk; the loop runs only over the trained models.
      process_chunk = function(rows) {
        n_rows = base::length(rows)
        theta_out = base::matrix(NA_real_, nrow = n_rows, ncol = n_model)
        loss_out = base::matrix(NA_real_, nrow = n_rows, ncol = n_model)
        # Column that each sorted model writes to (duplicates share one column).
        col_for_model = base::match(re_train, re_train_cols)

        TFE_k_c = TFE_k[rows, , drop = FALSE]
        TFE_next_c = TFE_next[rows, , drop = FALSE]
        TGE_k_c = TGE_k[rows, , drop = FALSE]
        TGE_next_c = TGE_next[rows, , drop = FALSE]
        rho_c = rho_mat[rows, , drop = FALSE]
        # Interval lengths replicated per row; multiplying a matrix by a bare
        # length-8 vector would recycle down columns, so build an explicit frame.
        Tslot_c = base::matrix(
          Tslot_term, nrow = n_rows, ncol = base::length(Tslot_term), byrow = TRUE
        )

        for (j in base::seq_len(n_model_all)) {
          meta = model_meta[[j]]
          out_col = col_for_model[j]

          # Theta-independent saturations and feedback strength.
          vi = Hill_sat(TFE_k_c, meta$u1i, meta$u2i)
          vj = Hill_sat(TGE_k_c, meta$u1j, meta$u2j)
          vij = vi - vj
          vj_safe = base::pmin(base::pmax(vj, 1e-12), 1 - 1e-12)
          S_mat = S_cal(s = meta$s, v = vj_safe, v_tilde = meta$v_tilde)

          # Residuals decomposed as c + R * d for each consecutive group.
          c1 = TGE_k_c + meta$b1 - TGE_next_c
          d1 = meta$a1 * vij * Tslot_c
          c2 = TFE_k_c + meta$b2 - TFE_next_c
          d2 = rho_c * meta$a2 * S_mat * Tslot_c

          # Quadratic coefficients of the MSE in R.
          A = base::rowMeans(d1 * d1 + d2 * d2)
          B = base::rowMeans(2 * (c1 * d1 + c2 * d2))
          C = base::rowMeans(c1 * c1 + c2 * c2)

          R_lo = meta$R_lo
          R_hi = meta$R_hi
          par0 = base::min(base::max(tau, 0.01 + 1e-6), 0.99 - 1e-6)
          R_par0 = R_cal(
            theta_TF_TG = par0, tau = tau, r1 = 200, r2 = meta$r2,
            r3 = 1, r4 = 1000, r5 = 1
          )

          # Degenerate (flat-in-R) cases.
          flatA = A < 1e-12
          flatAB = flatA & base::abs(B) < 1e-12
          # Unconstrained minimizer R* = -B/(2A), guarded against zero A.
          R_star = -B / (2 * base::pmax(A, 1e-12))
          R_clamped = base::pmin(base::pmax(R_star, R_lo), R_hi)
          # Flat A but non-flat B: monotonic loss -> the sign of B picks the edge.
          edge_lo = flatA & B > 0
          edge_hi = flatA & B < 0
          R_clamped[edge_lo] = R_lo
          R_clamped[edge_hi] = R_hi
          # Fully flat: replicate the original Brent behaviour (start point).
          R_clamped[flatAB] = R_par0

          # Map R back to theta. Interior rows use the shared monotonic inverse
          # (approx for a close start, then vectorized Newton refinement).
          at_lo = R_clamped <= R_lo
          at_hi = R_clamped >= R_hi
          interior = !(at_lo | at_hi | flatAB)
          theta = base::rep(par0, n_rows)
          if (base::any(interior)) {
            R_t = R_clamped[interior]
            th = stats::approx(
              x = meta$R_grid, y = meta$theta_grid, xout = R_t,
              rule = 2, ties = "ordered"
            )$y
            h = 1e-6
            for (it in 1:3) {
              th_plus = base::pmin(th + h, 0.99)
              th_minus = base::pmax(th - h, 0.01)
              Rn = R_cal(
                theta_TF_TG = th, tau = tau, r1 = 200, r2 = meta$r2,
                r3 = 1, r4 = 1000, r5 = 1
              )
              Rp = (R_cal(
                theta_TF_TG = th_plus, tau = tau, r1 = 200, r2 = meta$r2,
                r3 = 1, r4 = 1000, r5 = 1
              ) - R_cal(
                theta_TF_TG = th_minus, tau = tau, r1 = 200, r2 = meta$r2,
                r3 = 1, r4 = 1000, r5 = 1
              )) / (th_plus - th_minus)
              Rp = base::pmax(Rp, 1e-12)
              th = th - (Rn - R_t) / Rp
              th = base::pmin(base::pmax(th, 0.01), 0.99)
            }
            theta[interior] = th
          }
          theta[at_lo] = 0.01
          theta[at_hi] = 0.99
          theta[flatAB] = par0

          # Final R and the minimum loss evaluated at the optimum.
          R_final = R_cal(
            theta_TF_TG = theta, tau = tau, r1 = 200, r2 = meta$r2,
            r3 = 1, r4 = 1000, r5 = 1
          )
          loss = A * R_final * R_final + B * R_final + C

          theta_out[, out_col] = theta
          loss_out[, out_col] = loss
        }

        theta_df = base::as.data.frame(theta_out)
        loss_df = base::as.data.frame(loss_out)
        base::names(theta_df) = re_train_cols
        base::names(loss_df) = re_train_cols
        base::rownames(theta_df) = re_pred[rows]
        base::rownames(loss_df) = re_pred[rows]
        base::list(theta_p = theta_df, loss = loss_df)
      }

      # Evaluate chunks in parallel (fork). On a single core, run serially.
      ncores = base::as.integer(ncores)
      if (!(base::is.finite(ncores) && ncores > 0)) {
        base::stop("'ncores' must be a positive integer.")
      }
      if (ncores > 1 && n_part > 1) {
        prediction_list = parallel::mclapply(
          X = chunk_bounds,
          FUN = process_chunk,
          mc.cores = base::min(ncores, n_part)
        )
      } else {
        prediction_list = base::lapply(chunk_bounds, process_chunk)
      }

      # Validate every chunk result so a failed fork cannot silently drop pairs.
      if (base::length(prediction_list) != n_part) {
        base::stop("Expected ", n_part, " chunk results but got ",
                   base::length(prediction_list), "; reduce 'ncores' or 'each'.")
      }
      for (s in chunk_idx) {
        x = prediction_list[[s]]
        if (!base::is.list(x) ||
            !base::all(c("theta_p", "loss") %in% base::names(x)) ||
            !base::is.data.frame(x[["theta_p"]]) ||
            !base::is.data.frame(x[["loss"]]) ||
            base::nrow(x[["theta_p"]]) != base::length(chunk_bounds[[s]])) {
          base::stop("Chunk ", s, " returned an invalid result; ",
                     "reduce 'ncores' or 'each' and retry.")
        }
      }
      base::gc()

      # Bind chunks into branch-level frames (chunk order is the pair order).
      message("Summarizing theta_p of this branch.")
      branch_result = base::list(
        theta_p = base::do.call(base::rbind, base::lapply(prediction_list, `[[`, "theta_p")),
        loss = base::do.call(base::rbind, base::lapply(prediction_list, `[[`, "loss"))
      )
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]] = branch_result

      # Save the branch result.
      if (save_check_point) {
        message("Saving the predictive results of this branch.")
        qs::qsave(
          branch_result,
          file = base::paste0("./interest_cell_type_pGRN.",
                              base::gsub(" ", "_", cell_type), ".branch", n, ".qsave")
        )
      }
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
