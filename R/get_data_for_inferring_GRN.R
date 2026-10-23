#' @title Prepare Data for GRN Inference
#'
#' @description Builds inference-ready datasets by combining branch-level iGRN
#' TF-TG pairs with group expression features from cell grouping.
#' Non-zero pairs are located directly from the (often sparse) branch iGRN
#' matrices via vectorized indexing, and the TFE/TGE feature binding is fully
#' vectorized (C-level matrix operations), eliminating the original row-by-row
#' R loops. When several branches are present the independent branches are
#' processed in parallel by forking; a single branch uses an optional
#' within-branch partition controlled by \code{n_part}. Both branch-level full
#' feature sets and a global de-duplicated pair set are returned.
#'
#' @param interest_cell_type_branch_iGRN List. Output of
#'   \code{\link{get_branch_iGRN}}, containing branch-specific iGRN matrices
#'   with target genes as rows and TFs as columns.
#' @param interest_cell_type_group List. Output of \code{\link{cell_grouping}},
#'   containing the aggregated TFE/TGE features for each branch.
#' @param n_part Positive integer. Number of within-branch partitions; used
#'   only for a single branch (with multiple branches the branch itself is the
#'   parallel unit).
#' Defaults to 1.
#' @param ncores Positive integer. Number of cores for parallel processing via
#'   \code{\link[parallel]{mclapply}}.
#' Defaults to 1; on Windows it is forced to 1 (no forking).
#'
#' @return Nested list organized by cell type, each entry containing:
#' \itemize{
#'   \item \code{branch1}, \code{branch2}, ...: branch data frames with one row
#'     per non-zero TF-TG pair, including \code{TG}, \code{TF}, \code{theta_i}
#'     and the group features \code{TFE_K} and \code{TGE_K} for each group K.
#'   \item \code{GRN}: global de-duplicated pair set across branches, containing
#'     only \code{TG}, \code{TF} and \code{theta_i}.
#' }
#' Pairs with \code{theta_i = 0} (or \code{NA}) are removed from both branch
#' and global sets.
#'
#' @details For each cell type, the function proceeds as follows:
#' \enumerate{
#'   \item Locates non-zero entries of each branch iGRN matrix with
#'     \code{which(mat != 0)} and converts the linear indices to TG/TF
#'     positions, yielding the non-zero pairs (with values) in column-major
#'     order in one vectorized step.
#'   \item Appends the TFE/TGE group features by vectorized matrix indexing:
#'     each selected TF's / TG's group-expression row is looked up to match the
#'     pair order.
#'   \item With multiple branches, builds the branches in parallel; with a
#'     single branch, optionally splits it into \code{n_part} partitions.
#'   \item Binds the non-zero pairs across branches and de-duplicates by the
#'     pair (TG, TF) to form the global set.
#' }
#' The pairs are generated directly from the row/column names of the branch
#' iGRN matrices, so every pair is guaranteed to belong to that branch's
#' TF/TG sets; the per-row membership test and row-by-row assignment of the
#' original implementation are therefore unnecessary and removed.
#'
#' @note
#' \itemize{
#'   \item The function creates and switches to the \code{5 Infer GRN} directory,
#'     and restores the original working directory on exit (including on error).
#'   \item Feature binding is memory-allocation/indexing bound, not CPU bound.
#'     For a single branch the vectorized serial build is already near optimal
#'     and within-branch splitting may be slightly slower from fork/recombine
#'     overhead; for multiple branches the cross-branch parallelism uses the
#'     independent branches and can speed the work up on a multi-socket server.
#'   \item Pairs with zero (or \code{NA}) regulatory strength are excluded
#'     automatically.
#'   \item The \code{GRN} element gives a global, de-duplicated view without
#'     group features.
#'   \item Parallel processing relies on \code{\link[parallel]{mclapply}};
#'     on Windows \code{ncores} is forced to 1.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_iGRN = base::readRDS(
#'   "./3 get iGRN/interest_cell_type_branch_iGRN.rds"
#' )
#' interest_cell_type_group = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' ncores = parallel::detectCores() - 1 # Linux
#' # ncores = 1 # Windows
#' n_part = ncores
#' interest_cell_type_branch_data_for_GRN_infer = get_data_for_inferring_GRN(
#'   interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
#'   interest_cell_type_group = interest_cell_type_group,
#'   n_part = n_part,
#'   ncores = ncores
#' )
#' base::saveRDS(
#'   interest_cell_type_branch_data_for_GRN_infer,
#'   file = "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds"
#' )
#' }
get_data_for_inferring_GRN = function(
    interest_cell_type_branch_iGRN,
    interest_cell_type_group,
    n_part = 1,
    ncores = 1
) {
  ### Start.
  t_start = base::Sys.time()
  base::message("Run: Getting data for inferring GRN for each branch. Time: ", t_start, ".")

  # --- Working directory (restored robustly, even on error) ---
  original_dir = base::getwd()
  new_folder = "5 Infer GRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    base::message("Folder created: ", new_folder, ".")
  } else {
    base::message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  base::message("The current working directory has been switched to: ", base::getwd(), ".")
  base::on.exit(
    {
      base::setwd(original_dir)
      base::message("The current working directory has been switched to: ", base::getwd(), ".")
      t_end = base::Sys.time()
      base::message("Finish: Getting data for inferring GRN for each branch. Time: ", t_end, ".")
      base::message("Running time: ", base::format(t_end - t_start))
    },
    add = TRUE
  )

  # --- Parallel settings: mclapply forks only on Unix-like systems ---
  is_windows = .Platform$OS.type == "windows"
  mc_cores = if (is_windows) 1L else base::max(1L, as.integer(ncores))
  n_part = base::max(1L, as.integer(n_part))

  # Contiguous, balanced index chunks covering 1:N (never more chunks than
  # rows, so it is safe for small datasets).
  make_chunks = function(N, n_part) {
    N = as.integer(N)
    if (N <= 0L) return(base::list())
    n_part = base::min(base::max(1L, n_part), N)
    chunk_size = base::ceiling(N / n_part)
    chunk_id = ((base::seq_len(N) - 1L) %/% chunk_size) + 1L
    base::unname(base::split(base::seq_len(N), chunk_id))
  }

  # Stop clearly if any forked worker failed. mclapply returns a non-data-frame
  # (a raw vector) for a child that errored or was killed (e.g. OOM); catching
  # it here prevents silently producing too few rows.
  validate_workers = function(lst, context) {
    good = base::vapply(lst, base::is.data.frame, logical(1))
    if (!base::all(good)) {
      base::stop(
        "A parallel worker failed while ", context, " (often out of memory or a ",
        "worker error). Reduce n_part/ncores or set ncores = 1. Bad job index: ",
        base::paste(base::which(!good), collapse = ","),
        call. = FALSE
      )
    }
    return(lst)
  }

  # Build a branch feature data frame from selected TG/TF/theta vectors, fully
  # vectorized. Zero-length inputs return a typed empty data frame.
  build_branch = function(TGs_sel, TFs_sel, theta_sel, TFE, TGE) {
    n_pair = base::length(theta_sel)
    nK = base::ncol(TFE)
    # Each selected TF's / TG's group-expression row is repeated to pair order.
    tfe = TFE[TFs_sel, , drop = FALSE]
    tge = TGE[TGs_sel, , drop = FALSE]
    odd = base::seq(1L, 2L * nK, by = 2L)   # TFE_K positions
    even = base::seq(2L, 2L * nK, by = 2L)  # TGE_K positions
    feat = base::matrix(NA_real_, nrow = n_pair, ncol = 2L * nK)
    if (n_pair > 0L) {
      feat[, odd] = tfe
      feat[, even] = tge
    }
    base::colnames(feat) = base::paste0(
      base::rep(c("TFE_", "TGE_"), times = nK),
      base::rep(base::seq_len(nK), each = 2L)
    )
    out = base::data.frame(
      TG = TGs_sel, TF = TFs_sel, theta_i = theta_sel, feat,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    # Note: paste0(character(0), "_to_", character(0)) returns "_to_" (length 1)
    # because of the length-1 constant, so only set row names when rows exist.
    if (n_pair > 0L) {
      base::rownames(out) = base::paste0(TFs_sel, "_to_", TGs_sel)
    }
    return(out)
  }

  # Extract the non-zero pairs of one branch iGRN matrix (column-major order).
  # Returns a list with the selected symbols/values plus the dense TFE/TGE.
  extract_pairs = function(n) {
    mat = base::as.matrix(branch_GRN[[n]])
    TGs = base::rownames(mat)
    TFs = base::colnames(mat)
    nTG = base::length(TGs)
    # Linear indices of non-zero, non-NA entries (column-major order).
    lin = base::which(mat != 0)
    row = ((lin - 1L) %% nTG) + 1L
    col = ((lin - 1L) %/% nTG) + 1L
    return(list(
      N_all = nTG * base::length(TFs),
      TGs_sel = TGs[row],
      TFs_sel = TFs[col],
      theta_sel = mat[lin],
      TFE = base::as.matrix(TFE_T[[n]]),
      TGE = base::as.matrix(TGE_T[[n]])
    ))
  }

  ### Process.
  interest_cell_type_branch_data_for_GRN_infer = base::list()

  for (cell_type in base::names(interest_cell_type_branch_iGRN)) {
    base::message("Process cell type ", cell_type, ".")
    branch_GRN = interest_cell_type_branch_iGRN[[cell_type]]
    n_branch = base::length(branch_GRN)

    TFE_T = interest_cell_type_group[[cell_type]][["Branches_TFE_T"]]
    TGE_T = interest_cell_type_group[[cell_type]][["Branches_TGE_T"]]

    # Build one complete branch serially (vectorized).
    worker_serial = function(n) {
      p = extract_pairs(n)
      return(build_branch(p$TGs_sel, p$TFs_sel, p$theta_sel, p$TFE, p$TGE))
    }

    # Build one branch using within-branch partitions (single-branch case).
    worker_within = function(n) {
      p = extract_pairs(n)
      N = base::length(p$theta_sel)
      if (N == 0L) {
        return(build_branch(p$TGs_sel, p$TFs_sel, p$theta_sel, p$TFE, p$TGE))
      }
      chunks = make_chunks(N, n_part)
      slices = parallel::mclapply(
        X = chunks,
        FUN = function(pos) build_branch(
          p$TGs_sel[pos], p$TFs_sel[pos], p$theta_sel[pos], p$TFE, p$TGE
        ),
        mc.cores = mc_cores
      )
      slices = validate_workers(slices, "binding within-branch features")
      return(base::do.call(base::rbind, slices))
    }

    # Choose the parallel granularity (never nest parallelism).
    cross_parallel = (n_branch > 1L && mc_cores > 1L)
    within_parallel = (!cross_parallel && n_branch == 1L &&
                         mc_cores > 1L && n_part > 1L)

    base::message("Selecting non-zero TF-TG pairs and adding group states (TFE, TGE).")
    if (cross_parallel) {
      base::message("Processing ", n_branch, " branches in parallel.")
      branch_dfs = parallel::mclapply(
        X = base::seq_len(n_branch),
        FUN = worker_serial,
        mc.cores = base::min(mc_cores, n_branch)
      )
      branch_dfs = validate_workers(branch_dfs, "processing branches")
    } else if (within_parallel) {
      base::message("Processing the single branch in ", n_part, " partitions.")
      branch_dfs = base::lapply(base::seq_len(n_branch), worker_within)
    } else {
      base::message("Processing branches serially.")
      branch_dfs = base::lapply(base::seq_len(n_branch), worker_serial)
    }

    for (n in base::seq_len(n_branch)) {
      base::message("Branch ", n, " has ", base::nrow(branch_dfs[[n]]),
                    " non-zero candidate pairs.")
    }

    # Global de-duplicated pair set from the non-zero branch pairs.
    pair_frames = base::lapply(branch_dfs, function(df) {
      q = df[, c("TG", "TF", "theta_i"), drop = FALSE]
      base::rownames(q) = NULL
      return(q)
    })
    df_all = base::do.call(base::rbind, pair_frames)
    df_all = df_all[!base::duplicated(df_all[, c("TG", "TF")]), , drop = FALSE]
    if (base::nrow(df_all) > 0L) {
      base::rownames(df_all) = base::paste0(df_all$TF, "_to_", df_all$TG)
    }

    base::names(branch_dfs) = base::paste0("branch", base::seq_len(n_branch))
    interest_cell_type_branch_data_for_GRN_infer[[cell_type]] =
      c(branch_dfs, base::list(GRN = df_all))
  }

  ### End (working directory + timing handled by on.exit).
  return(interest_cell_type_branch_data_for_GRN_infer)
}
