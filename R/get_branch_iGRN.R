#' @title Extract Branch-Specific Initial Gene Regulatory Networks (iGRNs)
#'
#' @description Extracts branch-specific initial gene regulatory networks (iGRNs) from the complete cell-type-specific iGRNs.
#' For each branch of each cell type,
#' the global TG-by-TF association matrix is subsetted to the transcription factors (TFs) and target genes (TGs) present in that branch,
#' without recomputing any regulatory strength.
#'
#' @param interest_cell_type_iGRN Named list. Output of \code{\link{calculate_iGRN}} or \code{\link{filter_iGRN_by_ChIP_Atlas_priors}};
#' each element is a target gene (rows) by TF (columns) association matrix or data frame.
#' @param interest_cell_type_group Named list. Output of \code{\link{cell_grouping}},
#' providing per-branch \code{Branches_TFE_T} (TFs as rows) and \code{Branches_TGE_T} (TGs as rows) matrices.
#' Its cell type names must match those in \code{interest_cell_type_iGRN}.
#' @param ncores Positive integer. Number of cores for parallel computation via \code{\link[parallel]{mclapply}}.
#' Defaults to \code{1}. Parallel processing is only supported on Unix-like systems; it is forced to \code{1} on Windows.
#'
#' @return A named list with cell types as names;
#' each element is a list of branch-specific iGRNs aligned with the branch order in \code{\link{cell_grouping}}.
#' Each branch iGRN is a data frame with target genes as rows and TFs as columns,
#' containing the regulatory strengths inherited from the complete network for matched pairs (unmatched pairs stay zero).
#' All-zero rows (target genes) and all-zero columns (TFs) are removed;
#' a branch with no supported edge is returned as a 0 x 0 data frame.
#'
#' @details
#' For each cell type and each branch, the function:
#' \enumerate{
#'   \item Reads the branch TG set from \code{Branches_TGE_T} and TF set from \code{Branches_TFE_T}.
#'   \item Initializes an all-zero subnetwork with those rows and columns.
#'   \item Copies the strengths of TG-TF pairs present in both the branch and the complete network
#'         (vectorized submatrix extraction; no strength is recalculated).
#'   \item Removes all-zero rows and columns.
#' }
#'
#' @note
#' \itemize{
#'   \item This function only extracts and filters existing values; it never recomputes scores.
#'   \item Cell type names in the two inputs must match exactly.
#'   \item Parallel processing runs at the cell type level via \code{\link[parallel]{mclapply}}.
#'   \item The function creates (if missing) and enters \code{3 get iGRN/}, restoring the original
#'     working directory on exit.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
#' interest_cell_type_group = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_branch_iGRN = get_branch_iGRN(
#'   interest_cell_type_iGRN = interest_cell_type_iGRN,
#'   interest_cell_type_group = interest_cell_type_group,
#'   ncores = ncores
#' )
#' base::saveRDS(
#'   interest_cell_type_branch_iGRN,
#'   file = "./3 get iGRN/interest_cell_type_branch_iGRN.rds"
#' )
#' }
get_branch_iGRN = function(
    interest_cell_type_iGRN,
    interest_cell_type_group,
    ncores = 1
) {
  ### Validate inputs.
  if (!base::is.list(interest_cell_type_iGRN) || base::length(interest_cell_type_iGRN) == 0 ||
      base::is.null(base::names(interest_cell_type_iGRN))) {
    base::stop("'interest_cell_type_iGRN' must be a non-empty named list of iGRN matrices.")
  }
  if (!base::is.list(interest_cell_type_group) || base::length(interest_cell_type_group) == 0 ||
      base::is.null(base::names(interest_cell_type_group))) {
    base::stop("'interest_cell_type_group' must be a non-empty named list (output of cell_grouping()).")
  }
  missing_in_group = base::setdiff(base::names(interest_cell_type_iGRN),
                                   base::names(interest_cell_type_group))
  missing_in_iGRN  = base::setdiff(base::names(interest_cell_type_group),
                                   base::names(interest_cell_type_iGRN))
  if (base::length(missing_in_group) > 0 || base::length(missing_in_iGRN) > 0) {
    base::stop("Cell type names of the two inputs must match. ",
               "Only in iGRN: ", base::paste(missing_in_group, collapse = ", "), "; ",
               "only in group: ", base::paste(missing_in_iGRN, collapse = ", "), ".")
  }
  for (cell_type in base::names(interest_cell_type_iGRN)) {
    if (!base::is.matrix(interest_cell_type_iGRN[[cell_type]]) &&
        !base::is.data.frame(interest_cell_type_iGRN[[cell_type]])) {
      base::stop(base::sprintf("iGRN element '%s' must be a matrix or data frame.", cell_type))
    }
    if (!base::all(c("Branches_TFE_T", "Branches_TGE_T") %in%
                   base::names(interest_cell_type_group[[cell_type]]))) {
      base::stop("group element '", cell_type, "' must contain Branches_TFE_T and Branches_TGE_T.")
    }
  }
  if (!base::is.numeric(ncores) || base::length(ncores) != 1 ||
      ncores < 1 || ncores != base::round(ncores)) {
    base::stop("'ncores' must be a single positive integer.")
  }
  if (.Platform$OS.type == "windows" && ncores > 1) {
    message("Parallel computing is not supported on Windows. ncores is forced to 1.")
    ncores = 1
  }

  ### Start.
  t_start = base::Sys.time()
  message("Run: Extracting branch-specific iGRNs for each cell type. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Define the per-cell-type worker.
  get_branch_iGRN_process = function(cell_type, full_iGRN, full_group) {
    full = full_iGRN[[cell_type]]
    if (!base::is.data.frame(full)) {
      full = base::as.data.frame(full, optional = TRUE)
    }
    full_TGs = base::rownames(full)
    full_TFs = base::colnames(full)

    branches_TFE = full_group[[cell_type]][["Branches_TFE_T"]] # Rows are branch TFs.
    branches_TGE = full_group[[cell_type]][["Branches_TGE_T"]] # Rows are branch TGs.
    n_branch = base::length(branches_TGE)
    if (base::length(branches_TFE) != n_branch) {
      base::stop("Cell type '", cell_type, "': numbers of TFE and TGE branches do not match.")
    }

    branch_iGRN = base::vector("list", n_branch)
    for (n in base::seq_len(n_branch)) {
      branch_TGs = base::rownames(branches_TGE[[n]])
      branch_TFs = base::rownames(branches_TFE[[n]])

      # All-zero subnetwork sized by the branch gene sets.
      sub = base::as.data.frame(base::matrix(
        0,
        nrow = base::length(branch_TGs),
        ncol = base::length(branch_TFs),
        dimnames = base::list(branch_TGs, branch_TFs)
      ))

      # Copy strengths only for TG-TF pairs present in both the branch and the complete network.
      common_TGs = base::intersect(branch_TGs, full_TGs)
      common_TFs = base::intersect(branch_TFs, full_TFs)
      if (base::length(common_TGs) > 0 && base::length(common_TFs) > 0) {
        sub[common_TGs, common_TFs] =
          base::as.matrix(full[common_TGs, common_TFs, drop = FALSE])
      }

      # Remove all-zero rows (TGs) and columns (TFs).
      keep_TG = base::rowSums(sub != 0) > 0
      keep_TF = base::colSums(sub != 0) > 0
      sub = sub[keep_TG, keep_TF, drop = FALSE]
      branch_iGRN[[n]] = sub

      if (base::nrow(sub) == 0 || base::ncol(sub) == 0) {
        message("[", cell_type, "] Branch ", n, " has no supported edge; returned as 0 x 0.")
      }
    }
    return(branch_iGRN)
  }

  ### Extract branch iGRNs in parallel at the cell type level.
  interest_cell_type_branch_iGRN = parallel::mclapply(
    X = base::names(interest_cell_type_iGRN),
    FUN = function(cell_type) {
      message("Extracting branch iGRNs for cell type ", cell_type, ".")
      get_branch_iGRN_process(
        cell_type = cell_type,
        full_iGRN = interest_cell_type_iGRN,
        full_group = interest_cell_type_group
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_branch_iGRN) = base::names(interest_cell_type_iGRN)

  # Surface child-process errors instead of returning silent 'try-error' objects.
  failed = base::vapply(
    interest_cell_type_branch_iGRN,
    function(x) methods::is(x, "try-error"),
    base::logical(1)
  )
  if (base::any(failed)) {
    base::stop("Branch iGRN extraction failed for cell type(s): ",
               base::paste(base::names(failed)[failed], collapse = ", "), ".")
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Extracting branch-specific iGRNs for each cell type. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_branch_iGRN)
}
