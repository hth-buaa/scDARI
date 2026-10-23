#' @title Select and Aggregate Best GRN Predictions
#'
#' @description Post-processes inferred GRN results.
#' For each candidate pair within a branch, the trained models with the lowest
#' prediction loss are selected and their predicted strengths and losses are
#' combined by loss-derived weights, with optional down-weighting of boundary
#' strengths.
#' Branch results are integrated into a cell-type summary by cell-count weighted
#' averaging, predicted strengths and losses are normalized within each TF, and
#' observed regulatory-strength ranks together with gene expression features
#' are attached.
#'
#' @param interest_cell_type_pGRN List. Output of \code{\link{infer_GRN}}.
#' @param interest_cell_type_genes_pseudotime_info List. Output of
#'   \code{\link{get_genes_pseudotime_info}}, containing branch cell counts for
#'   weighted integration.
#' @param interest_cell_type_iGRN List. Output of \code{\link{calculate_iGRN}} or
#'   \code{\link{filter_iGRN_by_ChIP_Atlas_priors}}, used to retrieve the
#'   observed regulatory strength \code{theta_i}.
#' @param Basic_Info List returned by \code{\link{get_expression_matrix}},
#'   containing the Seurat object \code{scRNAseq_for_GRN} from which gene-level
#'   expression features are computed.
#' @param interest_cell_type_tau Named numeric vector. Output of
#'   \code{\link{get_cell_type_tau}}, used to shift normalized positive
#'   strengths above \code{tau}.
#' @param topn_max Positive integer. Maximum number of models averaged per pair.
#' Defaults to 25.
#' @param topn_min Positive integer. Minimum number of models averaged per pair.
#' Defaults to 5.
#' @param topn_ratio Numeric between 0 and 1. Number of models averaged as a
#'   proportion of available trained models.
#' Defaults to 0.1.
#' @param fix_factor Numeric in (0, 1]. Factor applied to the weights of
#'   boundary strengths, with 1 applying no down-weighting.
#' Defaults to 1.0.
#' @param min_scale Numeric in (0, 1]. Multiplier on the smallest positive value
#'   when setting the min-max lower bound, applied to strength and loss
#'   normalization.
#' Defaults to 0.99.
#' @param each Positive integer. Maximum number of pairs per parallel chunk.
#' Defaults to 1000000.
#' @param ncores Positive integer. Number of cores for parallel processing via
#'   \code{\link[parallel]{mclapply}}.
#' Defaults to 1; set to 1 on Windows.
#'
#' @return The input list with a \code{pred_summary} frame per cell type.
#'
#' @details The per-pair "select top-n lowest-loss models and loss-weight them"
#' operation is a pure row-wise reduction: the column permutation that orders
#' the loss matrix is computed with \code{\link{apply}} (C-level ordering),
#' the selected strengths and losses are gathered by linear indexing, and the
#' triangular weights and weighted means are computed by matrix algebra. No
#' per-pair R loop, per-pair list append, or duplicated chunk frames are built.
#'
#' @note
#' \itemize{
#'   \item The function creates and switches to the \code{5 Infer GRN}
#'     directory, and restores the original working directory on exit.
#'   \item \code{fix_factor} must be in (0, 1].
#'   \item Reduce \code{each} if out-of-memory errors occur.
#'   \item Parallel processing relies on \code{\link[parallel]{mclapply}};
#'     set \code{ncores = 1} on Windows.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' interest_cell_type_genes_pseudotime_info = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#' interest_cell_type_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
#' Basic_Info = base::readRDS("./1.3 Data Preprocessing - Filter Expression Matrix/Basic_Info.rds")
#' interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
#' ncores = parallel::detectCores() - 1 # Linux
#' # ncores = 1 # Windows
#' interest_cell_type_pGRN = select_best_prediction(
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   interest_cell_type_iGRN = interest_cell_type_iGRN,
#'   Basic_Info = Basic_Info,
#'   interest_cell_type_tau = interest_cell_type_tau,
#'   each = 1000,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
select_best_prediction = function(
    interest_cell_type_pGRN,
    interest_cell_type_genes_pseudotime_info,
    interest_cell_type_iGRN,
    Basic_Info,
    interest_cell_type_tau,
    topn_max = 25,
    topn_min = 5,
    topn_ratio = 0.1,
    fix_factor = 1.0,
    min_scale = 0.99,
    each = 1000000,
    ncores = 1
) {
  ### Start.
  t_start = base::Sys.time()
  message("Run: Selecting best predictions. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "5 Infer GRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  base::on.exit({
    base::setwd(original_dir)
    message("The working directory has been restored to: ", base::getwd(), ".")
  }, add = TRUE)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  # Parameter validation.
  if (!(base::is.numeric(fix_factor) && base::length(fix_factor) == 1 &&
        fix_factor > 0 && fix_factor <= 1)) {
    base::stop("'fix_factor' must be a single value in (0, 1].")
  }
  if (!(base::is.numeric(min_scale) && base::length(min_scale) == 1 &&
        min_scale > 0 && min_scale <= 1)) {
    base::stop("'min_scale' must be a single value in (0, 1].")
  }

  ### Shared gene expression features, computed once for all cell types.
  message("Computing shared gene-level expression features.")
  scRNAseq = Seurat::NormalizeData(Basic_Info$scRNAseq_for_GRN, verbose = FALSE)
  rna_data = SeuratObject::GetAssayData(scRNAseq, layer = "data")
  n_cell = base::ncol(rna_data)
  nonzero_n = Matrix::rowSums(rna_data > 0)
  row_sum = Matrix::rowSums(rna_data)
  gene_features = base::data.frame(
    nonzero_freq = nonzero_n / n_cell,
    nonzero_mean = base::ifelse(nonzero_n > 0, row_sum / nonzero_n, 0),
    overall_mean = Matrix::rowMeans(rna_data)
  )
  gene_features[paste0(colnames(gene_features), "-norm_rank")] =
    base::lapply(gene_features, function(x) {
      r = base::rank(x, ties.method = "min", na.last = "keep")
      r / base::max(r, na.rm = TRUE)
    })
  # Named-vector views for fast feature attachment.
  gf_named = base::lapply(gene_features, function(x) {
    stats::setNames(base::as.numeric(x), base::rownames(gene_features))
  })

  ### Process.
  for (cell_type in base::names(interest_cell_type_pGRN)) {
    branch_index = base::seq_along(
      interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]]
    )
    n_branch = base::length(branch_index)
    if (n_branch < 1) base::stop("The number of branches is less than 1.")
    tau = base::as.numeric(interest_cell_type_tau[cell_type])
    cell_weights = base::as.numeric(
      interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]]
    )

    # Per-branch vectorized top-n selection.
    branch_best_list = base::list()
    for (n in branch_index) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      branch_result = interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]]
      theta_p_mat = base::as.matrix(branch_result[["theta_p"]])
      loss_mat = base::as.matrix(branch_result[["loss"]])
      storage.mode(theta_p_mat) = "double"
      storage.mode(loss_mat) = "double"
      N = base::nrow(theta_p_mat); M = base::ncol(theta_p_mat)
      message("Branch ", n, " has ", N, " candidate pairs.")
      # Number of models averaged, clamped to [topn_min, topn_max] and M.
      topn = base::min(base::max(base::round(M * topn_ratio), topn_min),
                       topn_max, M)

      # Column permutation ordering each row by loss; keep only the top-n.
      # The ordering is the only R-level loop, so parallelize it across rows.
      topn_order = function(mat) {
        base::t(base::apply(mat, 1, function(x) base::order(x)[base::seq_len(topn)]))
      }
      use_parallel = base::is.numeric(ncores) && ncores > 1 &&
        .Platform$OS.type == "unix"
      if (use_parallel) {
        n_workers = base::as.integer(base::min(ncores, N))
        block_size = base::ceiling(N / n_workers)
        chunk_id = ((base::seq_len(N) - 1L) %/% block_size) + 1L
        row_chunks = base::split(base::seq_len(N), chunk_id)
        order_list = parallel::mclapply(
          X = row_chunks,
          FUN = function(idx) topn_order(loss_mat[idx, , drop = FALSE]),
          mc.cores = n_workers
        )
        # Validate worker results (no failure, one integer row per input row).
        if (base::length(order_list) != base::length(row_chunks) ||
            base::any(vapply(order_list, function(x) !base::is.matrix(x), logical(1)))) {
          base::stop("Parallel ordering failed: a worker did not return a matrix.")
        }
        order_idx = base::do.call(base::rbind, order_list)
        if (base::nrow(order_idx) != N || base::ncol(order_idx) != topn) {
          base::stop("Parallel ordering failed: row count mismatch.")
        }
      } else {
        order_idx = topn_order(loss_mat)
      }
      ri = base::rep(base::seq_len(N), topn)
      ci = base::as.vector(order_idx)
      loss_top = base::matrix(loss_mat[cbind(ri, ci)], N, topn)
      theta_top = base::matrix(theta_p_mat[cbind(ri, ci)], N, topn)

      # Triangular weights from loss; equal weights when all selected losses tie.
      lmin = loss_top[, 1]
      lmax = loss_top[, topn]
      rng = lmax - lmin
      w = base::matrix(1, N, topn)
      tied = rng <= 0
      if (base::any(!tied)) {
        w[!tied, ] = (lmax[!tied] - loss_top[!tied, ]) / rng[!tied]
      }
      # Down-weight boundary strengths.
      boundary = theta_top < 0.0101 | theta_top > 0.9899
      w[boundary] = w[boundary] * fix_factor
      w_sum = base::rowSums(w)
      # Guard against a zero total weight (cannot happen with fix_factor > 0).
      w_sum = base::pmax(w_sum, .Machine$double.eps)
      theta_best = base::rowSums(theta_top * w) / w_sum
      loss_best = base::rowSums(loss_top * w) / w_sum

      # Remove boundary-invalid strengths and infinite losses.
      keep = theta_best >= 0.011 & base::is.finite(loss_best)
      branch_best = base::data.frame(
        theta_p = theta_best[keep], loss = loss_best[keep],
        row.names = base::rownames(theta_p_mat)[keep]
      )
      base::colnames(branch_best) = base::paste0(base::colnames(branch_best),
                                                 "_branch", n)
      interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]][["best_result"]] =
        branch_best
      branch_best_list[[n]] = branch_best
      base::rm(order_idx, loss_top, theta_top, w)
      base::gc(verbose = FALSE)
    }

    ### Build the integrated cell-type summary.
    message("Integrating branch results.")
    all_rows = base::sort(base::unique(base::unlist(base::lapply(
      branch_best_list, base::rownames
    ))))
    Ns = base::length(all_rows)
    theta_cols = base::paste0("theta_p_branch", branch_index)
    loss_cols = base::paste0("loss_branch", branch_index)
    cols = base::c(
      "TF", "TG", theta_cols, loss_cols,
      "FG.theta_p", "FG.Loss",
      "FG.theta_p_norm", "FG.theta_p_norm.R", "FG.Loss_norm", "FG.Loss_norm.R",
      "FG.theta_i", "FG.theta_i.R",
      "F.NZER", "F.NZER.R", "F.NZEA", "F.NZEA.R", "F.EA", "F.EA.R",
      "G.NZER", "G.NZER.R", "G.NZEA", "G.NZEA.R", "G.EA", "G.EA.R"
    )
    pred_summary = base::data.frame(base::matrix(
      NA_real_, nrow = Ns, ncol = base::length(cols),
      dimnames = base::list(all_rows, cols)
    ))
    pair_split = base::strsplit(all_rows, "_to_")
    pred_summary$TF = base::vapply(pair_split, `[`, "", 1L)
    pred_summary$TG = base::vapply(pair_split, `[`, "", 2L)
    # Fill branch columns via named-vector matching (fast, vectorized).
    for (b in branch_index) {
      bb = branch_best_list[[b]]
      th = stats::setNames(bb[[1]], base::rownames(bb))
      ls = stats::setNames(bb[[2]], base::rownames(bb))
      pred_summary[[theta_cols[b]]] = th[all_rows]
      pred_summary[[loss_cols[b]]] = ls[all_rows]
    }

    # Integrate strength and loss by branch cell counts.
    theta_B = base::as.matrix(pred_summary[, theta_cols, drop = FALSE])
    loss_B = base::as.matrix(pred_summary[, loss_cols, drop = FALSE])
    W = base::matrix(cell_weights, nrow = Ns, ncol = n_branch, byrow = TRUE)
    mask_t = !base::is.na(theta_B); theta_B[!mask_t] = 0
    mask_l = !base::is.na(loss_B); loss_B[!mask_l] = 0
    den_t = base::rowSums(W * mask_t)
    den_l = base::rowSums(W * mask_l)
    pred_summary$FG.theta_p = base::ifelse(
      den_t > 0, base::rowSums(theta_B * W * mask_t) / base::pmax(den_t, .Machine$double.eps), 0
    )
    pred_summary$FG.Loss = base::ifelse(
      den_l > 0, base::rowSums(loss_B * W * mask_l) / base::pmax(den_l, .Machine$double.eps), Inf
    )

    ### Within-TF normalization helper (build TG x TF matrix, scale per column).
    normalize_within_TF = function(value, tau_shift) {
      tg = pred_summary$TG; tf = pred_summary$TF
      all_tg = base::sort(base::unique(tg))
      all_tf = base::sort(base::unique(tf))
      Mx = base::matrix(0, nrow = base::length(all_tg), ncol = base::length(all_tf),
                        dimnames = base::list(all_tg, all_tf))
      Mx[cbind(base::match(tg, all_tg), base::match(tf, all_tf))] = value
      Mx_pos = Mx; Mx_pos[Mx_pos <= 0] = Inf
      col_min = base::apply(Mx_pos, 2, base::min)
      col_max = base::apply(Mx, 2, base::max)
      active = col_max > 0
      out = base::matrix(0, base::nrow(Mx), base::ncol(Mx),
                         dimnames = base::dimnames(Mx))
      if (base::any(active)) {
        min_before = col_min[active] * min_scale
        max_before = col_max[active] + 0.0001
        denom = max_before - min_before
        scaled = (Mx[, active, drop = FALSE] - base::matrix(
          min_before, nrow = base::nrow(Mx), ncol = base::sum(active), byrow = TRUE
        )) / base::matrix(denom, nrow = base::nrow(Mx), ncol = base::sum(active),
                          byrow = TRUE)
        if (tau_shift) scaled = tau + (1 - tau) * scaled
        out[, active] = scaled
      }
      out[Mx <= 0] = 0
      out[out < 0] = 0
      # Map back to the pred_summary rows.
      out[cbind(base::match(tg, base::rownames(out)),
                base::match(tf, base::colnames(out)))]
    }

    # Strengths: tau-shifted.
    message("Normalizing predicted strengths within each TF.")
    pred_summary$FG.theta_p_norm = normalize_within_TF(
      pred_summary$FG.theta_p, tau_shift = TRUE
    )
    # Losses: no shift.
    message("Normalizing losses within each TF.")
    pred_summary$FG.Loss_norm = normalize_within_TF(
      pred_summary$FG.Loss, tau_shift = FALSE
    )

    ### Global normalized ranks (ties "min") over the positive values.
    rank_over_positive = function(v, zero_rank = NA_real_) {
      out = base::rep(NA_real_, base::length(v))
      pos = !base::is.na(v) & v > 0
      if (base::any(pos)) {
        r = base::rank(v[pos], ties.method = "min") / base::sum(pos)
        out[pos] = r
      }
      if (base::length(zero_rank) == 1 && !base::is.na(zero_rank)) {
        non_pos = !pos
        out[non_pos] = zero_rank
      }
      out
    }
    pred_summary$FG.theta_p_norm.R = rank_over_positive(
      pred_summary$FG.theta_p_norm, zero_rank = 0
    )
    pred_summary$FG.Loss_norm.R = rank_over_positive(
      pred_summary$FG.Loss_norm, zero_rank = NA_real_
    )

    ### Observed theta_i from the iGRN matrix.
    message("Attaching observed theta_i and its rank.")
    iGRN = base::as.matrix(interest_cell_type_iGRN[[cell_type]])
    pred_summary$FG.theta_i = base::as.numeric(
      iGRN[cbind(base::match(pred_summary$TG, base::rownames(iGRN)),
                 base::match(pred_summary$TF, base::colnames(iGRN)))]
    )
    # theta_i rank across the whole iGRN positive values, 0 for theta_i == 0.
    pos_i = base::which(base::as.numeric(iGRN) > 0)
    theta_i_rank = base::rep(NA_real_, base::length(pred_summary$FG.theta_i))
    if (base::length(pos_i) > 0) {
      r = base::rank(iGRN[pos_i], ties.method = "min") / base::length(pos_i)
      # Build a value -> rank lookup (unique values).
      vals = iGRN[pos_i]
      lookup = base::data.frame(v = vals, r = r)
      lookup = lookup[!base::duplicated(lookup$v), ]
      theta_i_rank = lookup$r[base::match(pred_summary$FG.theta_i, lookup$v)]
    }
    theta_i_rank[pred_summary$FG.theta_i == 0] = 0
    pred_summary$FG.theta_i.R = theta_i_rank

    ### Attach the precomputed TF and TG gene features.
    message("Attaching gene-level expression features.")
    attach_features = function(gene_col, feature) {
      gf_named[[feature]][pred_summary[[gene_col]]]
    }
    pred_summary$F.NZER = attach_features("TF", "nonzero_freq")
    pred_summary$F.NZER.R = attach_features("TF", "nonzero_freq-norm_rank")
    pred_summary$G.NZER = attach_features("TG", "nonzero_freq")
    pred_summary$G.NZER.R = attach_features("TG", "nonzero_freq-norm_rank")
    pred_summary$F.NZEA = attach_features("TF", "nonzero_mean")
    pred_summary$F.NZEA.R = attach_features("TF", "nonzero_mean-norm_rank")
    pred_summary$G.NZEA = attach_features("TG", "nonzero_mean")
    pred_summary$G.NZEA.R = attach_features("TG", "nonzero_mean-norm_rank")
    pred_summary$F.EA = attach_features("TF", "overall_mean")
    pred_summary$F.EA.R = attach_features("TF", "overall_mean-norm_rank")
    pred_summary$G.EA = attach_features("TG", "overall_mean")
    pred_summary$G.EA.R = attach_features("TG", "overall_mean-norm_rank")

    interest_cell_type_pGRN[[cell_type]][["pred_summary"]] = pred_summary
  }

  ### End.
  t_end = base::Sys.time()
  message("Finish: Selecting best predictions. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start), ".")
  return(interest_cell_type_pGRN)
}
