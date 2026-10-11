#' @ title Select Training Sets for Gene Regulatory Network Modeling
#'
#' @description Extracts high-confidence TF-TG interactions from branch-level initial GRNs,
#' to build training sets for downstream regulatory strength prediction.
#' Per-TF and per-TG top-ranked interactions are first combined into an initial set,
#' then subjected to multi-layer threshold filtering,
#' with optional mutual-top refinement and TF coverage protection,
#' after which pseudotime-resolved group expression features are appended for each branch.
#'
#' @param interest_cell_type_branch_iGRN Named list. Output of \code{\link{get_branch_iGRN}},
#' containing branch-specific iGRN matrices for each cell type.
#' @param interest_cell_type_tau Named numeric vector. Output of \code{\link{get_cell_type_tau}},
#' giving one binarization threshold per cell type;
#' its names must match the cell types.
#' @param interest_cell_type_group Named list. Output of \code{\link{cell_grouping}},
#' containing the aggregated per-branch matrices \code{Branches_TFE_T} and \code{Branches_TGE_T}.
#' @param top_n_for_max_theta_in_each_regulon Positive integer.
#' Number of top interactions retained separately for each TF (its top targets)
#' and for each TG (its top regulators);
#' the union of the two selections forms the initial training set.
#' Defaults to \code{1}.
#' @param max_theta_retention_rate_in_tops Numeric in \[0, 1\].
#' Fraction retained by \code{theta_i} within the initial training set;
#' the cutoff is the \code{1 - rate} quantile and edges are kept inclusively.
#' Defaults to \code{1} (no filtering).
#' @param max_theta_retention_rate_in_all Numeric in \[0, 1\].
#' Fraction retained by \code{theta_i} over the global candidate set;
#' the cutoff is the \code{1 - rate} quantile and edges are kept inclusively.
#' Defaults to \code{1} (no filtering).
#' @param retention_number_of_tops Positive integer.
#' Maximum number of interactions kept by absolute descending \code{theta_i} rank
#' within the initial training set.
#' Defaults to \code{1000}.
#' @param refine Logical.
#' Whether to apply mutual-top refinement.
#' Defaults to \code{TRUE}.
#' @param keep_all_TF Logical.
#' When \code{refine = TRUE},
#' whether to add back the strongest pre-refinement edge of each TF lost during refinement.
#' Defaults to \code{TRUE}.
#' Added edges keep their original filter flags,
#' and are therefore not guaranteed to survive the final triple filters.
#' @param top_ratio_in_refine Numeric in \[0, 1\].
#' Fraction retained by \code{theta_i} after mutual-top refinement;
#' only used when \code{refine = TRUE}.
#' Defaults to \code{1} (no filtering).
#' @param ncores Positive integer.
#' Number of cores for \code{\link[parallel]{mclapply}}.
#' Defaults to \code{1}, and is forced to \code{1} on Windows.
#'
#' @return A nested list \code{result[[cell_type]][[branch]]},
#' with an additional \code{all} element for each cell type:
#' \itemize{
#'   \item \code{branch1}, \code{branch2}, ...:
#'     data frames of selected pairs whose TF and TG both occur in that branch's subnetwork,
#'     with columns \code{TF}, \code{TG}, \code{theta_i},
#'     followed by the interleaved group features
#'     \code{TFE_1, TGE_1, ..., TFE_K, TGE_K},
#'     where K is the number of cell groups in that branch (K may differ across branches);
#'     feature values are within-group means of non-zero values,
#'     and may be \code{NA} when a gene has no non-zero observation in that group.
#'   \item \code{all}:
#'     the merged global training set of the cell type,
#'     containing only \code{TF}, \code{TG}, \code{theta_i} (no group features).
#' }
#' Row names are \code{<TF>_to_<TG>}.
#' Cell types or branches with no pair at or above \code{tau}
#' contain zero-row data frames with the \code{TF}, \code{TG}, \code{theta_i} columns.
#'
#' @details For each cell type the workflow is:
#' \enumerate{
#'   \item Builds all TF-TG pairs of every branch from the branch subnetwork,
#'     attaches \code{theta_i},
#'     drops pairs below the cell-type \code{tau},
#'     then merges branches and deduplicates into a global candidate set.
#'   \item Builds the initial training set as the union of,
#'     computed separately,
#'     the top \code{top_n_for_max_theta_in_each_regulon} targets of each TF
#'     and the top regulators of each TG.
#'   \item Marks three criteria on that initial set (not applied yet):
#'     the within-set quantile,
#'     the global-candidate quantile,
#'     and the absolute top-\code{retention_number_of_tops} rank.
#'   \item If \code{refine = TRUE},
#'     keeps only pairs whose TF is a top regulator of some TG
#'     and whose TG is a top target of some TF,
#'     then applies the \code{top_ratio_in_refine} quantile.
#'     If \code{keep_all_TF = TRUE},
#'     the strongest pre-refinement edge of each TF lost in this step is added back.
#'   \item Applies the three marked criteria together (all must be TRUE),
#'     and formats the surviving pairs.
#'   \item Maps the global set back to each branch,
#'     keeps pairs whose TF and TG both occur in that branch's subnetwork,
#'     and appends the interleaved \code{TFE_*}/\code{TGE_*} group features.
#' }
#'
#' @note
#' \itemize{
#'   \item Rate parameters must lie in \[0, 1\],
#'     and integer parameters must be single positive integers;
#'     the per-TF/per-TG top count is additionally clipped inside the worker
#'     (values below 1.5 become 1, otherwise it is rounded).
#'   \item Edges added back by \code{keep_all_TF} are not exempt from the final triple filters,
#'     so full TF coverage is not guaranteed under stringent cutoffs.
#'   \item Refinement can greatly reduce the edge count;
#'     tune \code{top_ratio_in_refine} to balance coverage and confidence.
#'   \item The function creates and enters \code{4 Build Prediction Model},
#'     restoring the original working directory on exit.
#'   \item Parallelization runs at the cell type level via \code{\link[parallel]{mclapply}},
#'     and is not supported on Windows.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_branch_iGRN.rds")
#' interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
#' interest_cell_type_group = base::readRDS(
#'   "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds"
#' )
#' interest_cell_type_branch_training_set = select_training_set(
#'   interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
#'   interest_cell_type_tau = interest_cell_type_tau,
#'   interest_cell_type_group = interest_cell_type_group,
#'   top_n_for_max_theta_in_each_regulon = 1,
#'   max_theta_retention_rate_in_tops = 0.8,
#'   max_theta_retention_rate_in_all = 0.8,
#'   retention_number_of_tops = 1000,
#'   refine = TRUE,
#'   keep_all_TF = TRUE,
#'   top_ratio_in_refine = 0.8,
#'   ncores = parallel::detectCores() - 1 # in Linux; ncores = 1 on Windows.
#' )
#' base::saveRDS(
#'   interest_cell_type_branch_training_set,
#'   file = "./4 Build Prediction Model/interest_cell_type_branch_training_set.rds"
#' )
#' }
select_training_set = function(
    interest_cell_type_branch_iGRN,
    interest_cell_type_tau,
    interest_cell_type_group,
    top_n_for_max_theta_in_each_regulon = 1,
    max_theta_retention_rate_in_tops = 1,
    max_theta_retention_rate_in_all = 1,
    retention_number_of_tops = 1000,
    refine = TRUE,
    keep_all_TF = TRUE,
    top_ratio_in_refine = 1,
    ncores = 1
) {
  ### Validate inputs.
  if (!base::is.list(interest_cell_type_branch_iGRN) ||
      base::length(interest_cell_type_branch_iGRN) == 0 ||
      base::is.null(base::names(interest_cell_type_branch_iGRN))) {
    base::stop("'interest_cell_type_branch_iGRN' must be a non-empty named list (output of get_branch_iGRN()).")
  }
  if (!base::is.numeric(interest_cell_type_tau) ||
      base::is.null(base::names(interest_cell_type_tau))) {
    base::stop("'interest_cell_type_tau' must be a named numeric vector (output of get_cell_type_tau()).")
  }
  if (!base::is.list(interest_cell_type_group) ||
      base::is.null(base::names(interest_cell_type_group))) {
    base::stop("'interest_cell_type_group' must be a named list (output of cell_grouping()).")
  }
  cell_types = base::names(interest_cell_type_branch_iGRN)
  if (!base::setequal(base::names(interest_cell_type_tau), cell_types) ||
      !base::setequal(base::names(interest_cell_type_group), cell_types)) {
    base::stop("Cell type names of interest_cell_type_branch_iGRN, interest_cell_type_tau and interest_cell_type_group must match.")
  }
  for (cell_type in cell_types) {
    group_ct = interest_cell_type_group[[cell_type]]
    if (!base::all(c("Branches_TFE_T", "Branches_TGE_T") %in% base::names(group_ct))) {
      base::stop("group element '", cell_type, "' must contain Branches_TFE_T and Branches_TGE_T.")
    }
    n_branch = base::length(interest_cell_type_branch_iGRN[[cell_type]])
    if (base::length(group_ct[["Branches_TFE_T"]]) != n_branch ||
        base::length(group_ct[["Branches_TGE_T"]]) != n_branch) {
      base::stop("Cell type '", cell_type, "': branch counts of iGRN, TFE_T and TGE_T must match.")
    }
  }
  if (!base::is.numeric(top_n_for_max_theta_in_each_regulon) ||
      base::length(top_n_for_max_theta_in_each_regulon) != 1 ||
      top_n_for_max_theta_in_each_regulon < 1 ||
      top_n_for_max_theta_in_each_regulon != base::round(top_n_for_max_theta_in_each_regulon)) {
    base::stop("'top_n_for_max_theta_in_each_regulon' must be a single positive integer.")
  }
  if (!base::is.numeric(retention_number_of_tops) ||
      base::length(retention_number_of_tops) != 1 ||
      retention_number_of_tops < 1 ||
      retention_number_of_tops != base::round(retention_number_of_tops)) {
    base::stop("'retention_number_of_tops' must be a single positive integer.")
  }
  rate_params = c(max_theta_retention_rate_in_tops,
                  max_theta_retention_rate_in_all,
                  top_ratio_in_refine)
  if (base::any(!base::is.numeric(rate_params)) ||
      base::any(rate_params < 0) || base::any(rate_params > 1)) {
    base::stop("'max_theta_retention_rate_in_tops', 'max_theta_retention_rate_in_all' and 'top_ratio_in_refine' must be in [0, 1].")
  }
  if (!base::is.logical(refine) || base::length(refine) != 1 ||
      !base::is.logical(keep_all_TF) || base::length(keep_all_TF) != 1) {
    base::stop("'refine' and 'keep_all_TF' must be single TRUE/FALSE values.")
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
  message("Run: Selecting Training Sets. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "4 Build Prediction Model"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Define the per-cell-type worker.
  select_training_set_process = function(branch_GRN,
                                         branch_tau,
                                         top_n_for_max_theta_in_each_regulon = 1,
                                         max_theta_retention_rate_in_tops = 1,
                                         max_theta_retention_rate_in_all = 1,
                                         retention_number_of_tops = 1000,
                                         refine = TRUE,
                                         keep_all_TF = TRUE,
                                         top_ratio_in_refine = 1,
                                         TFE_T,
                                         TGE_T) {
    n_branch = base::length(branch_GRN)

    # Build all TF-TG pairs per branch, attach theta_i, and drop pairs below tau.
    message("Building TF-TG pairs per branch and attaching theta_i.")
    all_combinations = base::list()
    for (n in base::seq_len(n_branch)) {
      key = base::paste0("branch_", n, "_GRN")
      comb = base::expand.grid(
        TG = base::rownames(branch_GRN[[n]]),
        TF = base::colnames(branch_GRN[[n]]),
        stringsAsFactors = FALSE
      )
      if (base::nrow(comb) > 0) {
        comb$theta_i = branch_GRN[[n]][base::as.matrix(comb[, c("TG", "TF")])]
      } else {
        comb$theta_i = base::numeric(0)
      }
      all_combinations[[key]] = base::subset(comb, theta_i >= branch_tau)
    }

    # Merge branches and deduplicate into the global candidate set.
    all_combinations[["GRN"]] = base::do.call(base::rbind, all_combinations)
    all_combinations[["GRN"]] =
      all_combinations[["GRN"]][!base::duplicated(all_combinations[["GRN"]]), ]

    # Multi-layer thresholding of the global candidate pairs.
    stratified_slice = function(data,
                                top_n_for_max_theta_in_each_regulon = 1,
                                max_theta_retention_rate_in_tops = 1,
                                max_theta_retention_rate_in_all = 1,
                                retention_number_of_tops = 1000,
                                refine = TRUE,
                                keep_all_TF = TRUE,
                                top_ratio_in_refine = 1) {
      # Empty candidate set: return the same three-column frame as the normal path.
      if (base::nrow(data) == 0) {
        message("No candidate pair at or above tau; returning an empty training set.")
        return(base::data.frame(
          TF = base::character(0), TG = base::character(0), theta_i = base::numeric(0),
          row.names = base::character(0)
        ))
      }

      # Clip the filtering parameters to valid ranges.
      message("Clipping filter parameters to valid ranges.")
      top_n_for_max_theta_in_each_regulon = base::ifelse(
        top_n_for_max_theta_in_each_regulon < 1.5,
        1,
        base::round(top_n_for_max_theta_in_each_regulon)
      )
      max_theta_retention_rate_in_tops =
        base::pmin(base::pmax(max_theta_retention_rate_in_tops, 0), 1)
      max_theta_retention_rate_in_all =
        base::pmin(base::pmax(max_theta_retention_rate_in_all, 0), 1)

      # Top-N targets per TF.
      message("Selecting the top ", top_n_for_max_theta_in_each_regulon,
              " interaction(s) by theta_i for each TF.")
      TFregulon_top_n = base::list()
      for (TF in base::unique(data$TF)) {
        TFregulon = data[data$TF == TF, ]
        TFregulon_top_n[[TF]] = dplyr::slice_max(
          .data = TFregulon, order_by = theta_i,
          n = top_n_for_max_theta_in_each_regulon, with_ties = FALSE
        )
      }
      TFregulon_top_n = data.table::rbindlist(TFregulon_top_n)

      # Top-N regulators per TG.
      message("Selecting the top ", top_n_for_max_theta_in_each_regulon,
              " interaction(s) by theta_i for each TG.")
      TGregulon_top_n = base::list()
      for (TG in base::unique(data$TG)) {
        TGregulon = data[data$TG == TG, ]
        TGregulon_top_n[[TG]] = dplyr::slice_max(
          .data = TGregulon, order_by = theta_i,
          n = top_n_for_max_theta_in_each_regulon, with_ties = FALSE
        )
      }
      TGregulon_top_n = data.table::rbindlist(TGregulon_top_n)

      # Initial training set: union of the per-TF and per-TG top-N pairs.
      message("Building the initial training set as the union of per-TF and per-TG top-N pairs.")
      train_set = base::unique(base::rbind(TFregulon_top_n, TGregulon_top_n))

      # Mark the three filter criteria (applied together only after optional refinement).
      message("Marking the three filter criteria on the initial training set.")
      train_set$max_theta_retention_rate_in_tops = base::ifelse(
        train_set$theta_i >= stats::quantile(
          train_set$theta_i, probs = 1 - max_theta_retention_rate_in_tops, na.rm = TRUE
        ), TRUE, FALSE
      )
      train_set$max_theta_retention_rate_in_all = base::ifelse(
        train_set$theta_i >= stats::quantile(
          data$theta_i, probs = 1 - max_theta_retention_rate_in_all, na.rm = TRUE
        ), TRUE, FALSE
      )
      train_set = dplyr::mutate(
        train_set,
        retention_number_of_tops = dplyr::row_number(dplyr::desc(theta_i)) <= retention_number_of_tops
      )
      message("Initial training set: ", base::nrow(train_set), " interaction(s), ",
              base::length(base::unique(train_set$TF)), " TF(s), ",
              base::length(base::unique(train_set$TG)), " TG(s).")

      # Optional mutual-top refinement followed by a theta_i quantile cutoff.
      if (refine == TRUE) {
        message("Running mutual-top refinement.")
        train_set0 = train_set
        train_set = train_set[
          (train_set$TF %in% base::unique(TGregulon_top_n$TF)) &
            (train_set$TG %in% base::unique(TFregulon_top_n$TG)),
        ]
        if (base::nrow(train_set) > 0) {
          message("Retaining the top ", top_ratio_in_refine * 100,
                  "% of interactions by theta_i after refinement.")
          train_set = train_set[train_set$theta_i >= stats::quantile(
            train_set$theta_i, probs = 1 - top_ratio_in_refine, na.rm = TRUE
          ), ]
        } else {
          message("No interaction survived mutual-top refinement before TF coverage rescue.")
        }
        message("After refinement: ", base::nrow(train_set), " interaction(s), ",
                base::length(base::unique(train_set$TF)), " TF(s), ",
                base::length(base::unique(train_set$TG)), " TG(s).")

        # Add back the strongest pre-refinement edge of each TF lost in refinement.
        if (keep_all_TF == TRUE) {
          TFs_filter = base::setdiff(base::unique(train_set0$TF), base::unique(train_set$TF))
          if (base::length(TFs_filter) > 0) {
            message("Adding back ", base::length(TFs_filter),
                    " TF(s) missing after refinement: ",
                    base::paste(base::sort(TFs_filter), collapse = ", "), ".")
            for (TF_filter in TFs_filter) {
              TF_filter_train_set = train_set0[train_set0$TF == TF_filter, , drop = FALSE]
              TF_filter_train_set =
                TF_filter_train_set[base::which.max(TF_filter_train_set$theta_i), , drop = FALSE]
              train_set = base::rbind(train_set, TF_filter_train_set)
            }
            message("After TF coverage rescue: ", base::nrow(train_set), " interaction(s), ",
                    base::length(base::unique(train_set$TF)), " TF(s), ",
                    base::length(base::unique(train_set$TG)), " TG(s).")
          } else {
            message("No TF lost in refinement needs to be added back.")
          }
        }
      }

      # Apply the three marked criteria together, then format the output.
      message("Applying the final triple filters and formatting the training set.")
      keep = base::rowSums(train_set[, c(
        "max_theta_retention_rate_in_tops",
        "max_theta_retention_rate_in_all",
        "retention_number_of_tops"
      )]) == 3
      train_set = train_set[keep, ]
      train_set = base::as.data.frame(train_set[, c("TF", "TG", "theta_i")])
      train_set = train_set[!base::duplicated(train_set), ]
      base::rownames(train_set) = base::paste0(train_set$TF, "_to_", train_set$TG)
      return(train_set)
    }

    slices_GRN = stratified_slice(
      data = all_combinations[["GRN"]],
      top_n_for_max_theta_in_each_regulon = top_n_for_max_theta_in_each_regulon,
      max_theta_retention_rate_in_tops = max_theta_retention_rate_in_tops,
      max_theta_retention_rate_in_all = max_theta_retention_rate_in_all,
      retention_number_of_tops = retention_number_of_tops,
      refine = refine,
      keep_all_TF = keep_all_TF,
      top_ratio_in_refine = top_ratio_in_refine
    )

    # Map the global set to each branch and append interleaved TFE_k / TGE_k group features.
    message("Appending group-resolved expression features (TFE and TGE) for each branch.")
    slices_branch_GRN = base::list()
    for (n in base::seq_len(n_branch)) {
      br = branch_GRN[[n]]
      tfe = TFE_T[[n]] # Rows are TFs, columns are cell groups.
      tge = TGE_T[[n]] # Rows are TGs, columns are cell groups.
      if (base::ncol(tfe) != base::ncol(tge)) {
        base::stop("Branch ", n, ": TFE_T and TGE_T have different numbers of cell groups.")
      }

      # Keep pairs whose TF and TG are both present in this branch's subnetwork.
      present = slices_GRN$TF %in% base::colnames(br) &
        slices_GRN$TG %in% base::rownames(br)
      sub = slices_GRN[present, , drop = FALSE]

      if (base::nrow(sub) > 0) {
        K = base::ncol(tfe)
        tfe_values = tfe[sub$TF, base::seq_len(K), drop = FALSE]
        tge_values = tge[sub$TG, base::seq_len(K), drop = FALSE]
        features = base::data.frame(row.names = base::rownames(sub))
        for (kk in base::seq_len(K)) {
          features[[base::paste0("TFE_", kk)]] = tfe_values[, kk]
          features[[base::paste0("TGE_", kk)]] = tge_values[, kk]
        }
        sub = base::cbind(sub, features)
      }
      slices_branch_GRN[[base::paste0("branch", n)]] = sub
    }
    slices_branch_GRN[["all"]] = slices_GRN
    return(slices_branch_GRN)
  }

  ### Select the training set of each cell type in parallel.
  interest_cell_type_branch_training_set = parallel::mclapply(
    X = cell_types,
    FUN = function(cell_type) {
      message("Selecting the training set for cell type ", cell_type, ".")
      select_training_set_process(
        branch_GRN = interest_cell_type_branch_iGRN[[cell_type]],
        branch_tau = interest_cell_type_tau[cell_type],
        top_n_for_max_theta_in_each_regulon = top_n_for_max_theta_in_each_regulon,
        max_theta_retention_rate_in_tops = max_theta_retention_rate_in_tops,
        max_theta_retention_rate_in_all = max_theta_retention_rate_in_all,
        retention_number_of_tops = retention_number_of_tops,
        refine = refine,
        keep_all_TF = keep_all_TF,
        top_ratio_in_refine = top_ratio_in_refine,
        TFE_T = interest_cell_type_group[[cell_type]][["Branches_TFE_T"]],
        TGE_T = interest_cell_type_group[[cell_type]][["Branches_TGE_T"]]
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_branch_training_set) = cell_types

  # Surface child-process errors instead of returning silent 'try-error' objects.
  failed = base::vapply(
    interest_cell_type_branch_training_set,
    function(x) methods::is(x, "try-error"),
    base::logical(1)
  )
  if (base::any(failed)) {
    base::stop("Training set selection failed for cell type(s): ",
               base::paste(base::names(failed)[failed], collapse = ", "), ".")
  }

  ### End.
  base::names(interest_cell_type_branch_training_set) = base::names(interest_cell_type_branch_iGRN)
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Selecting Training Sets. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_branch_training_set)
}
