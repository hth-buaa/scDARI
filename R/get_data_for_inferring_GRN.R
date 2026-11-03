#' @title Prepare Data for GRN Inference
#'
#' @description Expands every non-zero TF-TG pair of the branch-specific iGRN
#' matrices and binds the same group-level TFE/TGE features used in model
#' training, producing the inference-ready datasets consumed by
#' \code{\link{infer_GRN}}. Non-zero pairs are located directly from the
#' often-sparse iGRN matrices without dense materialization, and all feature
#' binding is vectorized. Both branch-level tables and a global de-duplicated
#' pair set are returned.
#'
#' @param interest_cell_type_branch_iGRN List. Output of
#'   \code{get_branch_iGRN}, named by cell type. Each cell type holds one iGRN
#'   matrix per branch with target genes as rows and transcription factors as
#'   columns; a matrix may be a dense matrix or a \code{sparseMatrix}. Non-zero
#'   entries define the candidate pairs, and their values supply the observed
#'   \code{theta_i}.
#' @param interest_cell_type_group List. Output of \code{cell_grouping}, named
#'   by cell type. The elements \code{Branches_TFE_T} and \code{Branches_TGE_T}
#'   hold, per branch, the aggregated TF-expression and TG-expression feature
#'   matrices (genes as rows and pseudotime groups as columns) that are aligned
#'   to the candidate pairs.
#' @param n_part Positive integer, default \code{1}. Number of within-branch
#'   partitions, used only when a single branch is processed; with several
#'   branches the branch itself is the parallel unit and this value is ignored.
#'   The default \code{1} means no splitting. Related to \code{ncores}: each
#'   partition is dispatched as one parallel worker.
#' @param ncores Positive integer, default \code{1}. Number of workers passed
#'   to \code{\link[parallel]{mclapply}}. On Windows forking is unavailable, so
#'   the value is forced to \code{1}. Related to \code{n_part}, which sets how
#'   many pieces the workers receive.
#'
#' @return A list named by cell type. Each entry contains:
#' \itemize{
#'   \item \code{branch1}, \code{branch2}, ...: data frames with one row per
#'     non-zero pair, including \code{TG}, \code{TF}, \code{theta_i} and the
#'     interleaved group features \code{TFE_K}, \code{TGE_K} for every group K.
#'     Row names are \code{<TF>_to_<TG>}.
#'   \item \code{GRN}: the global pair set de-duplicated by \code{(TG, TF)} across
#'     branches, retaining only \code{TG}, \code{TF} and \code{theta_i}.
#' }
#'
#' @details For each cell type the function runs the following pipeline.
#' \enumerate{
#'   \item Extracts non-zero entries in column-major order without coercing a
#'     sparse matrix to dense: \code{Matrix::summary} for sparse matrices and
#'     \code{which(mat != 0)} for dense matrices, with \code{NA} entries
#'     discarded.
#'   \item Builds the branch table by vectorized row lookup, interleaving each
#'     selected TF's and TG's group features into \code{TFE_K}/\code{TGE_K}.
#'   \item Schedules the work: parallel across independent branches, or
#'     optionally partitions one branch into \code{n_part} pieces.
#'   \item Binds the branch pairs and de-duplicates them by \code{(TG, TF)} to
#'     form the global \code{GRN}.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and enters the \code{5 Infer GRN} folder and
#'     restores the original working directory on exit, including on error.
#'   \item The work is indexing and memory-allocation bound rather than CPU
#'     bound; the serial vectorized build is already fast, and splitting a
#'     single branch can be slightly slower because of fork and recombine
#'     overhead.
#'   \item Pairs with zero or \code{NA} strength are excluded automatically.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_iGRN = qs::qread(
#'   "./3 get iGRN/interest_cell_type_branch_iGRN.qsave"
#' )
#' interest_cell_type_group = qs::qread(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.qsave"
#' )
#' ncores = parallel::detectCores() - 1
#' interest_cell_type_branch_data_for_GRN_infer = get_data_for_inferring_GRN(
#'   interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
#'   interest_cell_type_group = interest_cell_type_group,
#'   n_part = ncores,
#'   ncores = ncores
#' )
#' qs::qsave(interest_cell_type_branch_data_for_GRN_infer,
#'           file = "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.qsave")
#' }
get_data_for_inferring_GRN = function(
    interest_cell_type_branch_iGRN,
    interest_cell_type_group,
    n_part = 1,
    ncores = 1
) {
  t_start = Sys.time()
  message("Run: Getting data for inferring GRN. Time: ", t_start, ".")

  original_dir = getwd()
  new_folder = "5 Infer GRN"
  if (!dir.exists(new_folder)) dir.create(new_folder, recursive = TRUE)
  setwd(new_folder)

  on.exit({
    setwd(original_dir)
    message("Finish: Getting data. Running time: ",
            format(Sys.time() - t_start))
  }, add = TRUE)

  is_windows = .Platform$OS.type == "windows"
  mc_cores = if (is_windows) 1L else max(1L, as.integer(ncores))
  n_part = max(1L, as.integer(n_part))

  ### Sparse-aware extraction of non-zero entries
  extract_nonzero = function(mat) {
    if (methods::is(mat, "sparseMatrix")) {

      # Read the internal column-major triplets without a dense copy.
      sm = Matrix::summary(mat)
      row = sm$i; col = sm$j; val = sm$x
    } else {
      mat = as.matrix(mat)
      nr = nrow(mat)

      # Dense path: linear indices are already in column-major order.
      lin = which(mat != 0)
      row = ((lin - 1L) %% nr) + 1L
      col = ((lin - 1L) %/% nr) + 1L
      val = mat[lin]
    }

    # Drop NA values to match the which(mat != 0) contract.
    ok = !is.na(val)
    list(row = row[ok], col = col[ok], val = val[ok])
  }

  ### Vectorized construction of one branch table
  build_branch = function(TGs_sel, TFs_sel, theta_sel, TFE, TGE) {
    n_pair = length(theta_sel)
    nK = ncol(TFE)
    tfe = TFE[TFs_sel, , drop = FALSE]
    tge = TGE[TGs_sel, , drop = FALSE]
    odd = seq(1L, 2L * nK, by = 2L)
    even = seq(2L, 2L * nK, by = 2L)
    feat = matrix(NA_real_, nrow = n_pair, ncol = 2L * nK)

    if (n_pair > 0L) {
      feat[, odd] = tfe
      feat[, even] = tge
    }
    colnames(feat) = paste0(rep(c("TFE_", "TGE_"), times = nK),
                            rep(seq_len(nK), each = 2L))
    out = data.frame(TG = TGs_sel, TF = TFs_sel, theta_i = theta_sel, feat,
                     check.names = FALSE, stringsAsFactors = FALSE)

    if (n_pair > 0L) rownames(out) = paste0(TFs_sel, "_to_", TGs_sel)
    out
  }

  interest_cell_type_branch_data_for_GRN_infer = list()

  for (cell_type in names(interest_cell_type_branch_iGRN)) {
    message("Process cell type ", cell_type, ".")
    branch_GRN = interest_cell_type_branch_iGRN[[cell_type]]
    n_branch = length(branch_GRN)
    TFE_T = interest_cell_type_group[[cell_type]][["Branches_TFE_T"]]
    TGE_T = interest_cell_type_group[[cell_type]][["Branches_TGE_T"]]

    # Build one branch through sparse extraction and vectorized binding.
    worker_serial = function(n) {
      mat = branch_GRN[[n]]
      ez = extract_nonzero(mat)

      # Read names directly from the matrix, avoiding a dense as.matrix copy.
      TGs = rownames(mat); TFs = colnames(mat)
      build_branch(TGs[ez$row], TFs[ez$col], ez$val,
                   as.matrix(TFE_T[[n]]), as.matrix(TGE_T[[n]]))
    }

    # Build one branch with optional within-branch partitioning.
    worker_within = function(n) {
      mat = branch_GRN[[n]]
      ez = extract_nonzero(mat)
      TGs = rownames(mat); TFs = colnames(mat)
      TGs_sel = TGs[ez$row]; TFs_sel = TFs[ez$col]; theta_sel = ez$val
      N = length(theta_sel)

      if (N == 0L) {
        return(build_branch(TGs_sel, TFs_sel, theta_sel,
                            as.matrix(TFE_T[[n]]), as.matrix(TGE_T[[n]])))
      }
      chunk_size = ceiling(N / min(n_part, N))
      cid = ((seq_len(N) - 1L) %/% chunk_size) + 1L
      chunks = unname(split(seq_len(N), cid))
      slices = parallel::mclapply(
        chunks,
        function(pos) build_branch(TGs_sel[pos], TFs_sel[pos], theta_sel[pos],
                                   as.matrix(TFE_T[[n]]), as.matrix(TGE_T[[n]])),
        mc.cores = mc_cores
      )

      if (!all(vapply(slices, is.data.frame, logical(1)))) {
        stop("A parallel worker failed (often OOM). Reduce n_part/ncores.")
      }
      do.call(rbind, slices)
    }

    cross_parallel = (n_branch > 1L && mc_cores > 1L)
    within_parallel = (!cross_parallel && n_branch == 1L &&
                         mc_cores > 1L && n_part > 1L)

    if (cross_parallel) {
      branch_dfs = parallel::mclapply(seq_len(n_branch), worker_serial,
                                      mc.cores = min(mc_cores, n_branch))

      if (!all(vapply(branch_dfs, is.data.frame, logical(1))))
        stop("A branch worker failed (often OOM).")
    } else if (within_parallel) {
      branch_dfs = lapply(seq_len(n_branch), worker_within)
    } else {
      branch_dfs = lapply(seq_len(n_branch), worker_serial)
    }

    for (n in seq_len(n_branch))
      message("Branch ", n, " has ", nrow(branch_dfs[[n]]),
              " non-zero candidate pairs.")

    ### Global de-duplicated pair set

    pair_frames = lapply(branch_dfs, function(df) {
      q = df[, c("TG", "TF", "theta_i"), drop = FALSE]; rownames(q) = NULL; q
    })
    df_all = do.call(rbind, pair_frames)
    df_all = df_all[!duplicated(df_all[, c("TG", "TF")]), , drop = FALSE]

    if (nrow(df_all) > 0L)
      rownames(df_all) = paste0(df_all$TF, "_to_", df_all$TG)
    names(branch_dfs) = paste0("branch", seq_len(n_branch))
    interest_cell_type_branch_data_for_GRN_infer[[cell_type]] =
      c(branch_dfs, list(GRN = df_all))
  }

  return(interest_cell_type_branch_data_for_GRN_infer)
}
