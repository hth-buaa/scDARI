#' @title Infer the Branch GRN
#'
#' @description Applies every trained model of a branch to each candidate
#' TF-TG pair and re-optimizes the regulatory strength for that pair,
#' returning the predicted strength and loss under each model. Candidate pairs
#' are split into fixed-size chunks; by default every chunk is written to disk
#' immediately and reduced to a lightweight file reference, so the large
#' pair-by-model matrices never accumulate in memory. An in-memory mode
#' reproduces the original all-in-memory behavior for small datasets.
#'
#' @param interest_cell_type_branch_data_for_GRN_infer List. Output of
#'   \code{\link{get_data_for_inferring_GRN}}, named by cell type, with one
#'   feature table per branch and a global \code{GRN} element.
#' @param interest_cell_type_tau Named numeric vector. Output of
#'   \code{\link{get_cell_type_tau}}. One tau value per cell
#'   type, used both in the gating function and to bound the optimized
#'   strength.
#' @param interest_cell_type_branch_model_train List. Output of
#'   \code{\link{model_train}}. Trained branch models, named by cell type.
#'   Each branch entry supplies the fitted parameters
#'   (\code{r2}, \code{alpha1/2}, \code{beta1/2}, \code{s}, \code{u1i/u2i},
#'   \code{u1j/u2j}, \code{v_tilde}) for every model.
#' @param interest_cell_type_group List. Output of \code{cell_group}; its
#'   \code{Branches_Tslot_k} elements provide the pseudotime interval lengths
#'   used in the residual terms.
#' @param select_cell_type Character vector or \code{NULL}, default
#'   \code{NULL}. Cell types to process. \code{NULL} means all cell types.
#' @param select_branch Integer vector or \code{NULL}, default \code{NULL}.
#'   Branch indices to process. \code{NULL} means all branches.
#' @param each Positive integer, default \code{100000}. Maximum number of
#'   candidate pairs per chunk, which directly controls peak memory. Smaller
#'   chunks use less memory but create more files. Related to \code{ncores}:
#'   size chunks so work is spread evenly across workers.
#' @param ncores Positive integer, default \code{1}. Number of workers passed
#'   to \code{\link[parallel]{mclapply}}; set to \code{1} on Windows.
#' @param save_check_point Logical, default \code{FALSE}. Retained for
#'   interface compatibility; in disk mode every chunk is already persisted.
#' @param keep_in_memory Logical, default \code{FALSE}. If \code{FALSE},
#'   chunks are streamed to disk and the branch node is a reference. If
#'   \code{TRUE}, all chunk frames are combined and held in memory, which
#'   reproduces the original behavior but may exhaust memory.
#' @param chunk_dir Character or \code{NULL}, default \code{NULL}. Directory
#'   for chunk files. \code{NULL} resolves to \code{pGRN_chunks} under the
#'   current working directory; the directory is created and normalized.
#'
#' @return A list named by cell type. Each branch node is either:
#' \itemize{
#'   \item a disk reference containing \code{.storage = "disk"},
#'     \code{chunk_files}, \code{chunk_nrow}, \code{model_cols},
#'     \code{n_model} and \code{n_pair}, when \code{keep_in_memory = FALSE};
#'     or
#'   \item the in-memory frames \code{theta_p} and \code{loss} (pairs as
#'     rows, models as columns), when \code{keep_in_memory = TRUE}.
#' }
#'
#' @details For a fixed pair and model, the prediction loss is a quadratic in
#' the effect efficiency \eqn{R}:
#' \deqn{\mathrm{MSE}(R) = A R^2 + B R + C,}
#' built from the consecutive-group residuals. The per-pair pipeline is:
#' \enumerate{
#'   \item Extracts the shared feature matrices once, including adjacent-group
#'     TF/TG features, the pseudotime intervals and the observed saturation
#'     \code{rho = TFE/(TFE + 0.01)}.
#'   \item For each model, computes the Hill-saturated TF/TG terms, applies
#'     the \code{S_cal} transform, and aggregates the per-group residuals
#'     into row-wise coefficients \code{A}, \code{B}, \code{C}.
#'   \item Forms the unconstrained minimum \eqn{R^* = -B/(2A)}, handles flat
#'     and edge cases, maps interior values through the shared monotone
#'     inverse \code{R_cal} grid, and refines them with three Newton steps.
#'   \item Converts the strength back through \code{R_cal} and evaluates
#'     \code{loss = A R R + B R + C}, using multiplication rather than
#'     exponentiation to preserve bit-identical results.
#' }
#'
#' @note
#' \itemize{
#'   \item Chunk files use \code{qs::qsave} when \pkg{qs} is available and
#'     fall back to \code{saveRDS} otherwise.
#'   \item The function enters the \code{5 Infer GRN} folder and restores the
#'     original working directory on exit, including on error.
#'   \item Reduce \code{each} if out-of-memory errors occur; use
#'     \code{select_cell_type} and \code{select_branch} to process in batches.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = infer_GRN(
#'   interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
#'   interest_cell_type_tau = interest_cell_type_tau,
#'   interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
#'   interest_cell_type_group = interest_cell_type_group,
#'   each = 100000,
#'   ncores = parallel::detectCores() - 1
#' )
#' }
infer_GRN = function(
    interest_cell_type_branch_data_for_GRN_infer,
    interest_cell_type_tau,
    interest_cell_type_branch_model_train,
    interest_cell_type_group,
    select_cell_type = NULL,
    select_branch = NULL,
    each = 100000,
    ncores = 1,
    save_check_point = FALSE,
    keep_in_memory = FALSE,
    chunk_dir = NULL
) {
  t_start = Sys.time()
  message("Run: Inferring branch GRN. Time: ", t_start, ".")

  original_dir = getwd()
  new_folder = "5 Infer GRN"
  if (!dir.exists(new_folder)) dir.create(new_folder, recursive = TRUE)
  setwd(new_folder)
  on.exit(setwd(original_dir), add = TRUE)

  ### Chunk writer with an rds fallback

  .write_qs = function(x, file_base) {
    if (requireNamespace("qs", quietly = TRUE)) {
      f = file_base
      qs::qsave(x, file = f)
    } else {
      f = paste0(file_base, ".rds")
      saveRDS(x, file = f)
    }
    f
  }

  # Absolute root directory for chunk files, immune to later setwd calls.
  if (is.null(chunk_dir)) chunk_dir = file.path(getwd(), "pGRN_chunks")
  if (!dir.exists(chunk_dir)) dir.create(chunk_dir, recursive = TRUE)
  chunk_dir = normalizePath(chunk_dir, mustWork = TRUE)

  ### Shape-preserving Hill saturation

  # Keep the matrix shape instead of flattening it with as.numeric.
  Hill_sat = function(x, K, n) {
    storage.mode(x) = "double"
    xn = x^n
    xn / (K + xn)
  }

  interest_cell_type_pGRN = list()
  range_cell_type = if (is.null(select_cell_type))
    names(interest_cell_type_branch_data_for_GRN_infer) else select_cell_type

  for (cell_type in range_cell_type) {
    interest_cell_type_pGRN[[cell_type]] = list()
    tau = as.numeric(interest_cell_type_tau[cell_type])
    range_branch = if (is.null(select_branch))
      seq_len(length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1) else
        select_branch

    # File-name-safe cell type tag.
    ct_tag = gsub("[^A-Za-z0-9]", "_", cell_type)

    for (n in range_branch) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      state_data = interest_cell_type_branch_data_for_GRN_infer[[cell_type]][[
        paste0("branch", n)]]
      re_pred = rownames(state_data)
      N_pred_all = length(re_pred)
      message(N_pred_all, " candidate pairs will be predicted.")

      branch_models = interest_cell_type_branch_model_train[[cell_type]][[
        paste0("branch", n)]]
      re_train = sort(names(branch_models))
      re_train_cols = unique(re_train)
      n_model = length(re_train_cols)
      n_model_all = length(re_train)
      message(n_model_all, " trained models will be applied.")

      ### Shared feature matrices, extracted once

      TFE_col = grep("TFE_", colnames(state_data), value = TRUE)
      TGE_col = grep("TGE_", colnames(state_data), value = TRUE)
      n_group = length(TFE_col)
      Tslot_K = as.numeric(
        interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]])
      term_idx = seq_len(n_group - 1)
      Tslot_term = Tslot_K[term_idx]

      TFE_mat = as.matrix(state_data[, TFE_col, drop = FALSE])
      TGE_mat = as.matrix(state_data[, TGE_col, drop = FALSE])
      storage.mode(TFE_mat) = "double"; storage.mode(TGE_mat) = "double"
      TFE_k = TFE_mat[, term_idx, drop = FALSE]
      TFE_next = TFE_mat[, term_idx + 1L, drop = FALSE]
      TGE_k = TGE_mat[, term_idx, drop = FALSE]
      TGE_next = TGE_mat[, term_idx + 1L, drop = FALSE]
      rho_mat = TFE_k / (TFE_k + 0.01)

      ### Scalar parameters and the shared monotone inverse per model

      model_meta = stats::setNames(lapply(re_train, function(re_tr) {
        pa_best = branch_models[[re_tr]][["best_params"]]
        getp = function(nm) as.numeric(pa_best[nm])
        r2 = getp("r2")
        theta_grid = seq(0.01, 0.99, length.out = 2000)
        R_grid = scDARI::R_cal(theta_grid, tau = tau, r1 = 200, r2 = r2,
                               r3 = 1, r4 = 1000, r5 = 1)
        list(r2 = r2,
             a1 = getp("alpha1"), a2 = getp("alpha2"),
             b1 = getp("beta1"), b2 = getp("beta2"),
             s = getp("s"),
             u1i = getp("u1i"), u2i = getp("u2i"),
             u1j = getp("u1j"), u2j = getp("u2j"),
             v_tilde = getp("v_tilde"),
             R_lo = R_grid[1], R_hi = R_grid[length(R_grid)],
             theta_grid = theta_grid, R_grid = R_grid)
      }), re_train)

      each = as.integer(each)
      if (!(is.finite(each) && each > 0)) stop("'each' must be positive.")
      n_part = ceiling(N_pred_all / each)
      chunk_bounds = lapply(seq_len(n_part), function(s) {
        lo = each * (s - 1) + 1L
        hi = min(each * s, N_pred_all)
        seq.int(lo, hi)
      })

      ### Full computation for one chunk, identical to the original math

      process_chunk = function(rows) {
        n_rows = length(rows)
        theta_out = matrix(NA_real_, n_rows, n_model)
        loss_out = matrix(NA_real_, n_rows, n_model)
        col_for_model = match(re_train, re_train_cols)

        TFE_k_c = TFE_k[rows, , drop = FALSE]
        TFE_next_c = TFE_next[rows, , drop = FALSE]
        TGE_k_c = TGE_k[rows, , drop = FALSE]
        TGE_next_c = TGE_next[rows, , drop = FALSE]
        rho_c = rho_mat[rows, , drop = FALSE]
        Tslot_c = matrix(Tslot_term, n_rows, length(Tslot_term), byrow = TRUE)

        for (j in seq_len(n_model_all)) {
          meta = model_meta[[j]]
          out_col = col_for_model[j]

          vi = Hill_sat(TFE_k_c, meta$u1i, meta$u2i)
          vj = Hill_sat(TGE_k_c, meta$u1j, meta$u2j)
          vij = vi - vj
          vj_safe = pmin(pmax(vj, 1e-12), 1 - 1e-12)
          S_mat = scDARI::S_cal(s = meta$s, v = vj_safe, v_tilde = meta$v_tilde)

          c1 = TGE_k_c + meta$b1 - TGE_next_c
          d1 = meta$a1 * vij * Tslot_c
          c2 = TFE_k_c + meta$b2 - TFE_next_c
          d2 = rho_c * meta$a2 * S_mat * Tslot_c

          A = rowMeans(d1 * d1 + d2 * d2)
          B = rowMeans(2 * (c1 * d1 + c2 * d2))
          C = rowMeans(c1 * c1 + c2 * c2)

          R_lo = meta$R_lo; R_hi = meta$R_hi
          par0 = min(max(tau, 0.01 + 1e-6), 0.99 - 1e-6)
          R_par0 = scDARI::R_cal(par0, tau = tau, r1 = 200, r2 = meta$r2,
                                 r3 = 1, r4 = 1000, r5 = 1)

          flatA = A < 1e-12
          flatAB = flatA & abs(B) < 1e-12
          R_star = -B / (2 * pmax(A, 1e-12))
          R_clamped = pmin(pmax(R_star, R_lo), R_hi)
          edge_lo = flatA & B > 0
          edge_hi = flatA & B < 0
          R_clamped[edge_lo] = R_lo
          R_clamped[edge_hi] = R_hi
          R_clamped[flatAB] = R_par0

          at_lo = R_clamped <= R_lo
          at_hi = R_clamped >= R_hi
          interior = !(at_lo | at_hi | flatAB)
          theta = rep(par0, n_rows)
          if (any(interior)) {
            R_t = R_clamped[interior]
            th = stats::approx(meta$R_grid, meta$theta_grid, xout = R_t,
                               rule = 2, ties = "ordered")$y
            h = 1e-6
            for (it in 1:3) {
              th_plus = pmin(th + h, 0.99)
              th_minus = pmax(th - h, 0.01)
              Rn = scDARI::R_cal(th, tau = tau, r1 = 200, r2 = meta$r2,
                                 r3 = 1, r4 = 1000, r5 = 1)
              Rp = (scDARI::R_cal(th_plus, tau = tau, r1 = 200, r2 = meta$r2,
                                  r3 = 1, r4 = 1000, r5 = 1) -
                      scDARI::R_cal(th_minus, tau = tau, r1 = 200, r2 = meta$r2,
                                    r3 = 1, r4 = 1000, r5 = 1)) /
                (th_plus - th_minus)
              Rp = pmax(Rp, 1e-12)
              th = th - (Rn - R_t) / Rp
              th = pmin(pmax(th, 0.01), 0.99)
            }
            theta[interior] = th
          }
          theta[at_lo] = 0.01
          theta[at_hi] = 0.99
          theta[flatAB] = par0

          R_final = scDARI::R_cal(theta, tau = tau, r1 = 200, r2 = meta$r2,
                                  r3 = 1, r4 = 1000, r5 = 1)

          # Use R * R rather than R^2: mathematically equal but avoiding a
          # possible one-ULP floating-point difference from exponentiation.
          loss = A * R_final * R_final + B * R_final + C

          theta_out[, out_col] = theta
          loss_out[, out_col] = loss
        }

        theta_df = as.data.frame(theta_out)
        loss_df = as.data.frame(loss_out)
        names(theta_df) = re_train_cols
        names(loss_df) = re_train_cols
        rownames(theta_df) = re_pred[rows]
        rownames(loss_df) = re_pred[rows]
        list(theta_p = theta_df, loss = loss_df)
      }

      ncores = as.integer(ncores)
      if (!(is.finite(ncores) && ncores > 0))
        stop("'ncores' must be a positive integer.")

      if (keep_in_memory) {

        ### Compatibility path: full in-memory rbind, the original behavior

        if (ncores > 1 && n_part > 1) {
          prediction_list = parallel::mclapply(chunk_bounds, process_chunk,
                                               mc.cores = min(ncores, n_part))
        } else {
          prediction_list = lapply(chunk_bounds, process_chunk)
        }
        branch_result = list(
          theta_p = do.call(rbind, lapply(prediction_list, `[[`, "theta_p")),
          loss = do.call(rbind, lapply(prediction_list, `[[`, "loss"))
        )
        interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]] =
          branch_result
      } else {

        ### Memory-friendly path: persist each chunk and return paths only

        # Unique logical file base for every chunk.
        chunk_base = file.path(
          chunk_dir,
          sprintf("pGRN.%s.branch%d.chunk", ct_tag, n)
        )

        # Worker: compute one chunk, persist it, return path and row count.
        process_to_disk = function(s) {
          rows = chunk_bounds[[s]]
          res = process_chunk(rows)
          f = .write_qs(res, paste0(chunk_base, s))
          list(file = f, nrow = length(rows))
        }
        if (ncores > 1 && n_part > 1) {
          meta_list = parallel::mclapply(seq_len(n_part), process_to_disk,
                                         mc.cores = min(ncores, n_part))
        } else {
          meta_list = lapply(seq_len(n_part), process_to_disk)
        }

        # Guard against failed workers, which mclapply may return as non-lists.
        if (length(meta_list) != n_part ||
            !all(vapply(meta_list, is.list, logical(1)))) {
          stop("Chunk processing failed; reduce 'each' or 'ncores'.")
        }
        branch_ref = list(
          .storage = "disk",
          chunk_files = vapply(meta_list, `[[`, "", "file"),
          chunk_nrow = vapply(meta_list, `[[`, integer(1), "nrow"),
          model_cols = re_train_cols,
          n_model = n_model,
          n_pair = N_pred_all
        )
        interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]] =
          branch_ref
        message("Branch ", n, " streamed to ", n_part, " chunk files.")
      }
    }
  }

  t_end = Sys.time()
  message("Finish: Inferring branch GRN. Running time: ",
          format(t_end - t_start))
  return(interest_cell_type_pGRN)
}
