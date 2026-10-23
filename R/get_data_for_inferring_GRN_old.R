#' @title Prepare Data for GRN Inference
#'
#' @description Builds inference-ready datasets by combining branch-level iGRN
#' TF-TG pairs with group expression features from cell grouping.
#' Partition-based parallel processing is supported, and both branch-level full
#' feature sets and a global de-duplicated pair set are returned.
#'
#' @param interest_cell_type_branch_iGRN List. Output of
#'   \code{\link{get_branch_iGRN}}, containing branch-specific iGRN matrices
#'   with target genes as rows and TFs as columns.
#' @param interest_cell_type_group List. Output of \code{\link{cell_grouping}},
#'   containing the aggregated TFE/TGE features for each branch.
#' @param n_part Positive integer. Number of partitions for parallel feature
#' binding.
#' Defaults to 1.
#' @param ncores Positive integer. Number of cores for parallel processing via
#'   \code{\link[parallel]{mclapply}}.
#' Defaults to 1; set to 1 on Windows.
#'
#' @return Nested list organized by cell type, each entry containing:
#' \itemize{
#'   \item \code{branch1}, \code{branch2}, ...: branch data frames with one row
#'     per TF-TG pair, including \code{TG}, \code{TF}, \code{theta_i} and the
#'     group features \code{TFE_K} and \code{TGE_K} for each group K.
#'   \item \code{GRN}: global de-duplicated pair set across branches, containing
#'     only \code{TG}, \code{TF} and \code{theta_i}.
#' }
#' Pairs with \code{theta_i = 0} are removed from both branch and global sets.
#'
#' @details For each cell type, the function proceeds as follows:
#' \enumerate{
#'   \item Generates every TG-TF combination from each branch iGRN matrix and
#'     extracts the regulatory strength \code{theta_i}.
#'   \item Binds the branch pairs and removes duplicates to form the global set.
#'   \item Splits each branch pair set into \code{n_part} partitions, and appends
#'     the TFE/TGE group features, in parallel when requested.
#'   \item Recombines the partitions into the complete branch datasets.
#'   \item Removes pairs with \code{theta_i = 0} from branch and global sets.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and switches to the \code{5 Infer GRN} directory,
#'     and restores the original working directory on exit.
#'   \item Partition-based parallelism speeds up feature binding for large pair
#'     sets; choose \code{n_part} according to dataset size and memory.
#'   \item Pairs with zero regulatory strength are excluded automatically.
#'   \item The \code{GRN} element gives a global, de-duplicated view without
#'     group features.
#'   \item Parallel processing relies on \code{\link[parallel]{mclapply}};
#'     set \code{ncores = 1} on Windows.
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
#' interest_cell_type_branch_data_for_GRN_infer = get_data_for_inferring_GRN_old(
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
get_data_for_inferring_GRN_old = function(
    interest_cell_type_branch_iGRN,
    interest_cell_type_group,
    n_part = 1,
    ncores = 1
) {
  ### Start.
  t_start = base::Sys.time()
  message("Run: Getting data for inferring GRN for each branch. Time: ", t_start, ".")
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
  interest_cell_type_branch_data_for_GRN_infer = base::list()
  for (cell_type in base::names(interest_cell_type_branch_iGRN)) {
    message("Process cell type ", cell_type, ".")
    interest_cell_type_branch_data_for_GRN_infer[[cell_type]] = base::list()
    branch_GRN = interest_cell_type_branch_iGRN[[cell_type]]
    n_branch = base::length(branch_GRN)

    # Generate every TG-TF combination per branch and extract theta_i.
    message("Generating all TG-TF combinations and their regulatory strength theta_i.")
    all_combinations = base::list()
    for (n in 1:n_branch) {
      branch_name = base::paste0("branch_", n, "_GRN")
      all_combinations[[branch_name]] = base::expand.grid(
        TG = base::rownames(branch_GRN[[n]]),
        TF = base::colnames(branch_GRN[[n]]),
        stringsAsFactors = FALSE
      )
      all_combinations[[branch_name]]$theta_i = base::apply(
        all_combinations[[branch_name]], 1, function(x) branch_GRN[[n]][x[1], x[2]]
      )
      base::rownames(all_combinations[[branch_name]]) = base::paste0(
        all_combinations[[branch_name]]$TF, "_to_", all_combinations[[branch_name]]$TG
      )
    }

    # Bind branch pairs and de-duplicate into the global set.
    all_combinations[["GRN"]] = base::do.call(base::rbind, all_combinations)
    all_combinations[["GRN"]] =
      all_combinations[["GRN"]][!base::duplicated(all_combinations[["GRN"]]), ]
    base::rownames(all_combinations[["GRN"]]) = base::paste0(
      all_combinations[["GRN"]]$TF, "_to_", all_combinations[["GRN"]]$TG
    )

    # Append the TFE/TGE group features for each branch.
    message("Adding cell group states (TFE, TGE).")
    for (n in 1:n_branch) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      GRN = all_combinations[[base::paste0("branch_", n, "_GRN")]]
      N_re_all = base::nrow(GRN)
      message("Branch ", n, " has ", N_re_all, " candidate pairs before filtering.")

      # Partition the pairs.
      message("Partitioning data.")
      each = base::floor(N_re_all / n_part)
      GRN_list = base::list()
      for (s in 1:n_part) {
        if (s != n_part) {
          GRN_list[[s]] = GRN[(each * (s - 1) + 1):(each * s), ]
        } else {
          GRN_list[[s]] = GRN[(each * (s - 1) + 1):N_re_all, ]
        }
      }

      # Add group features to one partition.
      process_GRN_slice = function(s,
                                   n,
                                   GRN_list,
                                   TFs = base::colnames(branch_GRN[[n]]),
                                   TGs = base::rownames(branch_GRN[[n]]),
                                   TFE_T = interest_cell_type_group[[cell_type]][["Branches_TFE_T"]],
                                   TGE_T = interest_cell_type_group[[cell_type]][["Branches_TGE_T"]]) {
        current_branch = GRN_list[[s]]
        for (TF_TG in base::rownames(current_branch)) {
          if (current_branch[TF_TG, "TF"] %in% TFs && current_branch[TF_TG, "TG"] %in% TGs) {
            # Fill TFE and TGE for every group.
            for (K in base::seq_len(base::ncol(TFE_T[[n]]))) {
              current_branch[TF_TG, c(base::paste0("TFE_", K), base::paste0("TGE_", K))] =
                c(TFE_T[[n]][current_branch[TF_TG, "TF"], K],
                  TGE_T[[n]][current_branch[TF_TG, "TG"], K])
            }
          } else {
            # Drop pairs absent from the current branch's gene sets.
            current_branch = current_branch[base::rownames(current_branch) != TF_TG, ]
          }
        }
        return(current_branch)
      }

      # Bind features across partitions, in parallel when requested.
      message("Adding cell group states (TFE, TGE) in parallel.")
      data_for_prediction_list = parallel::mclapply(
        X = base::seq_len(n_part),
        FUN = function(s) {
          process_GRN_slice(
            s, n, GRN_list = GRN_list,
            TFs = base::colnames(branch_GRN[[n]]),
            TGs = base::rownames(branch_GRN[[n]]),
            TFE_T = interest_cell_type_group[[cell_type]][["Branches_TFE_T"]],
            TGE_T = interest_cell_type_group[[cell_type]][["Branches_TGE_T"]]
          )
        },
        mc.cores = ncores
      )
      df = base::do.call(base::rbind, data_for_prediction_list)

      # Remove zero-strength pairs.
      message("Removing candidate pairs with theta_i equal to zero.")
      df = df[df$theta_i != 0, , drop = FALSE]
      interest_cell_type_branch_data_for_GRN_infer[[cell_type]][[base::paste0("branch", n)]] = df
      message("Branch ", n, " has ", base::nrow(df), " candidate pairs after filtering.")
    }

    # Filter the global set.
    df_all = all_combinations[["GRN"]]
    df_all = df_all[df_all$theta_i != 0, , drop = FALSE]
    interest_cell_type_branch_data_for_GRN_infer[[cell_type]][["GRN"]] = df_all
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Getting data for inferring GRN for each branch. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_branch_data_for_GRN_infer)
}
