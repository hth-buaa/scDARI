#' @importFrom data.table :=
#'
#' @title Select the Best Prediction for Each TF-TG
#'
#' @description Collapses the many per-model predictions of each TF-TG pair
#' into one strength and one loss, integrates the branch-level results by
#' branch cell count, normalizes values within each TF, and attaches the
#' observed iGRN strength and the gene expression features. Disk-backed chunk
#' files are streamed, reduced and released one at a time, so large inference
#' outputs never need to be fully loaded; in-memory branch objects are
#' supported as well.
#'
#' @param interest_cell_type_pGRN List. Output of \code{\link{infer_GRN}},
#'   named by cell type. Each branch node is either a disk reference or an
#'   in-memory object holding \code{theta_p} and \code{loss}.
#' @param interest_cell_type_genes_pseudotime_info List. Output of
#' \code{\link{get_genes_pseudotime_info}}. Per cell type, provides
#'   \code{Branches_TFE} (the number of branches) and \code{Branches_n_cell}
#'   (the cell count used to weight integration).
#' @param interest_cell_type_iGRN List. Output of \code{\link{calculate_iGRN}}
#'   or \code{\link{filter_iGRN_by_ChIP_Atlas_priors}}. Per cell type, the
#'   observed iGRN matrix (targets as rows, TFs as columns), used to attach
#'   \code{theta_i}.
#' @param Basic_Info List returned by \code{\link{get_expression_matrix}}.
#'   Contains \code{scRNAseq_for_GRN}, either a Seurat object or a mock list
#'   with a count matrix, from which gene features are computed.
#' @param interest_cell_type_tau Named numeric vector. Output of
#'   \code{\link{get_cell_type_tau}}. One tau per cell type, used as the lower
#'   shift when normalizing strengths within a TF.
#' @param topn_max Positive integer greater than or equal to \code{topn_min},
#'   default \code{1}. Upper bound on the number of models averaged per pair.
#'   Related to \code{topn_min} and \code{topn_ratio}, which jointly determine
#'   the effective count.
#' @param topn_min Positive integer less than or equal to \code{topn_max},
#'   default \code{1}. Lower bound on the number of models averaged per pair.
#' @param topn_ratio Numeric in \code{(0, 1]}, default \code{0.20}. The
#'   effective count is \code{round(n_model * topn_ratio)}, then clamped to
#'   \code{[topn_min, topn_max]} and never above the number of models.
#' @param fix_factor Numeric in \code{(0, 1]}, default \code{1.00}. Factor
#'   applied to the weights of boundary strengths (below 0.0101 or above
#'   0.9899), which reflect the optimization bounds. \code{1} means no
#'   down-weighting.
#' @param min_scale Numeric in \code{[0, 1]}, default \code{0.00}. Multiplies
#'   the smallest positive value when setting the within-TF min-max lower
#'   bound. Related to \code{max_scale}, which sets the upper bound.
#' @param max_scale Numeric \code{>= 1}, default \code{1.00}. Multiplies the
#'   maximum value (plus a small offset) when setting the upper bound.
#' @param each Positive integer, default \code{1000000}. Retained for
#'   interface compatibility; the streaming path is driven by chunk files.
#' @param ncores Positive integer, default \code{1}. Retained for interface
#'   compatibility; chunk reduction is serial to avoid copying large data.
#'
#' @return A list named by cell type. Each entry keeps the branch nodes and
#'   gains \code{pred_summary}, a data frame with one row per pair containing
#'   the per-branch values, the integrated \code{FG.theta_p} and
#'   \code{FG.Loss}, their within-TF normalizations and ranks, the observed
#'   \code{FG.theta_i} and its rank, and the \code{F.*}/\code{G.*} TF/TG
#'   expression features.
#'
#' @details The per-pair "select the top-n lowest-loss models and weight them
#' by loss" step is a pure row-wise reduction: the column permutation that
#' orders the loss matrix is computed with \code{\link{apply}} (C-level
#' ordering), the selected values are gathered by linear indexing, and the
#' triangular weights and weighted means use matrix algebra. No per-pair R
#' loop or duplicated chunk frame is built. The branch values are then
#' integrated by cell counts, normalized within each TF, and augmented with
#' the observed strength and gene features.
#'
#' @note
#' \itemize{
#'   \item The function enters the \code{5 Infer GRN} folder and restores the
#'     original working directory on exit, including on error.
#'   \item Disk chunks are read, reduced and released sequentially, so peak
#'     memory is bounded by one chunk plus the small aggregated table.
#'   \item Missing branch values are ignored when computing the weighted means.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = select_best_prediction(
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   interest_cell_type_iGRN = interest_cell_type_iGRN,
#'   Basic_Info = Basic_Info,
#'   interest_cell_type_tau = interest_cell_type_tau
#' )
#' }
select_best_prediction = function(
    interest_cell_type_pGRN,
    interest_cell_type_genes_pseudotime_info,
    interest_cell_type_iGRN,
    Basic_Info,
    interest_cell_type_tau,
    topn_max = 1,
    topn_min = 1,
    topn_ratio = 0.20,
    fix_factor = 1.00,
    min_scale = 0.00,
    max_scale = 1.00,
    each = 1000000,
    ncores = 1
) {
  t_start = Sys.time()
  message("Run: Selecting best predictions. Time: ", t_start, ".")

  original_dir = getwd()
  new_folder = "5 Infer GRN"
  if (!dir.exists(new_folder)) dir.create(new_folder, recursive = TRUE)
  setwd(new_folder)
  on.exit(setwd(original_dir), add = TRUE)

  ### Parameter validation, matching the original checks
  if (!(is.numeric(fix_factor) && length(fix_factor) == 1 &&
        fix_factor > 0 && fix_factor <= 1))
    stop("'fix_factor' must be a single value in (0, 1].")
  if (!(is.numeric(min_scale) && length(min_scale) == 1 &&
        min_scale >= 0 && min_scale <= 1))
    stop("'min_scale' must be a single value in [0, 1].")
  if (!(is.numeric(max_scale) && length(max_scale) == 1 && max_scale >= 1))
    stop("'max_scale' must be a single value not less than 1.")

  ### Chunk reader for rds and qs files
  .read_qs = function(file) {
    if (grepl("\\.rds$", file)) readRDS(file) else qs::qread(file)
  }

  ### Expression matrix from a Seurat object or a mock list
  # Real Seurat input runs NormalizeData; the mock list returns its matrix
  # directly, so both local testing and the server run work.
  get_rna_data = function(sc) {
    if (is.list(sc) && !is.null(sc$.mock_counts)) return(sc$.mock_counts)
    if (requireNamespace("Seurat", quietly = TRUE) &&
        methods::is(sc, "Seurat")) {
      s = Seurat::NormalizeData(sc, verbose = FALSE)
      return(SeuratObject::GetAssayData(s, layer = "data"))
    }
    stop("Cannot extract an expression matrix from Basic_Info.")
  }

  ### Shared gene-level expression features
  message("Computing shared gene-level expression features.")
  rna_data = get_rna_data(Basic_Info$scRNAseq_for_GRN)
  n_cell = ncol(rna_data)
  nonzero_n = Matrix::rowSums(rna_data > 0)
  row_sum = Matrix::rowSums(rna_data)
  gene_features = data.frame(
    nonzero_freq = nonzero_n / n_cell,
    nonzero_mean = ifelse(nonzero_n > 0, row_sum / nonzero_n, 0),
    overall_mean = Matrix::rowMeans(rna_data)
  )
  gene_features[paste0(colnames(gene_features), "-norm_rank")] =
    lapply(gene_features, function(x) {
      r = rank(x, ties.method = "min", na.last = "keep")
      r / max(r, na.rm = TRUE)
    })
  gf_named = lapply(gene_features, function(x)
    stats::setNames(as.numeric(x), rownames(gene_features)))

  ### Top-n selection, triangular weighting and survival for one chunk
  # Inputs are the chunk theta_p/loss frames, with pairs as rows and models as columns.
  select_chunk_topn = function(theta_df, loss_df, topn, fix_factor) {
    theta_mat = as.matrix(theta_df)
    loss_mat = as.matrix(loss_df)
    N = nrow(theta_mat); M = ncol(theta_mat)

    if (topn == 1L) {
      # max.col on -loss resolves ties to the left, equal to order(x)[1].
      ci = max.col(-loss_mat, ties.method = "first")
      ri = seq_len(N)
      theta_top = matrix(theta_mat[cbind(ri, ci)], N, 1L)
      loss_top = matrix(loss_mat[cbind(ri, ci)], N, 1L)
    } else {
      order_idx = matrix(as.vector(apply(
        loss_mat, 1, function(x) order(x)[seq_len(topn)])),
        N, topn, byrow = TRUE)
      ri = rep(seq_len(N), topn); ci = as.vector(order_idx)
      loss_top = matrix(loss_mat[cbind(ri, ci)], N, topn)
      theta_top = matrix(theta_mat[cbind(ri, ci)], N, topn)
    }

    lmin = loss_top[, 1L]; lmax = loss_top[, topn]
    rng = lmax - lmin
    w = matrix(1, N, topn)
    tied = rng <= 0
    if (any(!tied))
      w[!tied, ] = (lmax[!tied] - loss_top[!tied, ]) / rng[!tied]
    boundary = theta_top < 0.0101 | theta_top > 0.9899
    w[boundary] = w[boundary] * fix_factor
    w_sum = pmax(rowSums(w), .Machine$double.eps)
    theta_best = rowSums(theta_top * w) / w_sum
    loss_best = rowSums(loss_top * w) / w_sum

    keep = theta_best >= 0.011 & is.finite(loss_best)
    data.frame(theta_p = theta_best[keep], loss = loss_best[keep],
               row.names = rownames(theta_mat)[keep])
  }

  for (cell_type in names(interest_cell_type_pGRN)) {
    branch_index = seq_along(
      interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]])
    n_branch = length(branch_index)
    if (n_branch < 1) stop("The number of branches is less than 1.")
    tau = as.numeric(interest_cell_type_tau[cell_type])
    cell_weights = as.numeric(
      interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]])

    branch_best_list = list()
    for (n in branch_index) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      node = interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]]

      # Model count for both disk references and in-memory objects.
      M = if (!is.null(node$.storage) && node$.storage == "disk")
        node$n_model else ncol(as.matrix(node[["theta_p"]]))
      topn = min(max(round(M * topn_ratio), topn_min), topn_max, M)
      message("topn = ", topn, " among ", M, " models.")

      if (!is.null(node$.storage) && node$.storage == "disk") {
        # Streaming: read, reduce and release each chunk file.
        chunk_best = vector("list", length(node$chunk_files))
        for (s in seq_along(node$chunk_files)) {
          ch = .read_qs(node$chunk_files[s])     # list(theta_p, loss)
          chunk_best[[s]] = select_chunk_topn(
            ch$theta_p, ch$loss, topn, fix_factor)
          rm(ch); if (s %% 4L == 0L) gc(verbose = FALSE)
        }
        branch_best = do.call(rbind, chunk_best)
      } else {
        # In-memory branch object produced with keep_in_memory = TRUE.
        branch_best = select_chunk_topn(
          node[["theta_p"]], node[["loss"]], topn, fix_factor)
      }
      colnames(branch_best) = paste0(colnames(branch_best), "_branch", n)

      # Keep the small best table on the branch node but drop the large matrices.
      interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]][[
        "best_result"]] = branch_best
      branch_best_list[[n]] = branch_best
      rm(branch_best); gc(verbose = FALSE)
    }

    ### Cross-branch integration, matching the original
    message("Integrating branch results.")
    all_rows = sort(unique(unlist(lapply(branch_best_list, rownames))))
    Ns = length(all_rows)
    theta_cols = paste0("theta_p_branch", branch_index)
    loss_cols = paste0("loss_branch", branch_index)
    cols = c("TF", "TG", theta_cols, loss_cols,
             "FG.theta_p", "FG.Loss",
             "FG.theta_p_norm", "FG.theta_p_norm.R",
             "FG.Loss_norm", "FG.Loss_norm.R",
             "FG.theta_i", "FG.theta_i.R",
             "F.NZER", "F.NZER.R", "F.NZEA", "F.NZEA.R", "F.EA", "F.EA.R",
             "G.NZER", "G.NZER.R", "G.NZEA", "G.NZEA.R", "G.EA", "G.EA.R")
    pred_summary = data.frame(matrix(
      NA_real_, Ns, length(cols), dimnames = list(all_rows, cols)))
    pair_split = strsplit(all_rows, "_to_")
    pred_summary$TF = vapply(pair_split, `[`, "", 1L)
    pred_summary$TG = vapply(pair_split, `[`, "", 2L)
    for (b in branch_index) {
      bb = branch_best_list[[b]]
      pred_summary[[theta_cols[b]]] = stats::setNames(bb[[1]],
                                                      rownames(bb))[all_rows]
      pred_summary[[loss_cols[b]]] = stats::setNames(bb[[2]],
                                                     rownames(bb))[all_rows]
    }

    theta_B = as.matrix(pred_summary[, theta_cols, drop = FALSE])
    loss_B = as.matrix(pred_summary[, loss_cols, drop = FALSE])
    W = matrix(cell_weights, Ns, n_branch, byrow = TRUE)
    mask_t = !is.na(theta_B); theta_B[!mask_t] = 0
    mask_l = !is.na(loss_B); loss_B[!mask_l] = 0
    den_t = rowSums(W * mask_t)
    den_l = rowSums(W * mask_l)
    pred_summary$FG.theta_p = ifelse(
      den_t > 0, rowSums(theta_B * W * mask_t) / pmax(den_t, .Machine$double.eps), 0)
    pred_summary$FG.Loss = ifelse(
      den_l > 0, rowSums(loss_B * W * mask_l) / pmax(den_l, .Machine$double.eps), Inf)

    ### Within-TF normalization by data.table grouping
    normalize_within_TF = function(value, tau_shift) {
      dt = data.table::data.table(.id = seq_along(value),
                                  TF = pred_summary$TF, v = value)

      # Positive minimum cmin and maximum cmax per TF (cmin = Inf if none).
      # dt[, data.table::`:=`(
      dt[, `:=`(
        .cmin = { p = v[v > 0 & !is.na(v)]; if (length(p)) min(p) else Inf },
        .cmax = suppressWarnings(max(v, na.rm = TRUE))
      ), by = TF]
      dt[, .active := .cmax > 0]
      # dt[, data.table::`:=`(.active, .cmax > 0)]
      min_before = dt$.cmin * min_scale
      max_before = (dt$.cmax + 0.0001) * max_scale
      denom = max_before - min_before
      sc = (dt$v - min_before) / denom
      if (tau_shift) sc = tau + (1 - tau) * sc
      sc[dt$v <= 0 | is.na(dt$v)] = 0
      sc[sc < 0] = 0
      sc[!dt$.active] = 0
      sc
    }
    message("Normalizing predicted strengths within each TF.")
    pred_summary$FG.theta_p_norm = normalize_within_TF(
      pred_summary$FG.theta_p, tau_shift = TRUE)
    message("Normalizing losses within each TF.")
    pred_summary$FG.Loss_norm = normalize_within_TF(
      pred_summary$FG.Loss, tau_shift = FALSE)

    ### Global positive-edge ranks with ties "min"
    rank_over_positive = function(v, zero_rank = NA_real_) {
      out = rep(NA_real_, length(v))
      pos = !is.na(v) & v > 0
      if (any(pos)) out[pos] = rank(v[pos], ties.method = "min") / sum(pos)
      if (!is.na(zero_rank)) out[!pos] = zero_rank
      out
    }
    pred_summary$FG.theta_p_norm.R = rank_over_positive(
      pred_summary$FG.theta_p_norm, zero_rank = 0)
    pred_summary$FG.Loss_norm.R = rank_over_positive(
      pred_summary$FG.Loss_norm, zero_rank = NA_real_)

    ### Attach observed theta_i and its rank
    message("Attaching observed theta_i and its rank.")
    iGRN = as.matrix(interest_cell_type_iGRN[[cell_type]])
    pred_summary$FG.theta_i = as.numeric(
      iGRN[cbind(match(pred_summary$TG, rownames(iGRN)),
                 match(pred_summary$TF, colnames(iGRN)))])
    pos_i = which(as.numeric(iGRN) > 0)
    theta_i_rank = rep(NA_real_, nrow(pred_summary))
    if (length(pos_i) > 0) {
      r = rank(iGRN[pos_i], ties.method = "min") / length(pos_i)
      vals = iGRN[pos_i]
      lookup = data.frame(v = vals, r = r)
      lookup = lookup[!duplicated(lookup$v), ]
      theta_i_rank = lookup$r[match(pred_summary$FG.theta_i, lookup$v)]
    }
    theta_i_rank[pred_summary$FG.theta_i == 0] = 0
    pred_summary$FG.theta_i.R = theta_i_rank

    ### Attach gene expression features
    message("Attaching gene-level expression features.")
    attach_features = function(gene_col, feature)
      gf_named[[feature]][pred_summary[[gene_col]]]
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

  t_end = Sys.time()
  message("Finish: Selecting best predictions. Running time: ",
          format(t_end - t_start))
  return(interest_cell_type_pGRN)
}
