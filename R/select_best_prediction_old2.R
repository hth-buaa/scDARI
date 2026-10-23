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
#'   \code{\link{filter_iGRN_by_ChIP_Atlas_priors}}. Cell-type-specific iGRN matrices, used
#'   to retrieve the observed regulatory strength \code{theta_i}.
#' @param Basic_Info List returned by \code{\link{get_expression_matrix}}. Basic information containing the Seurat object
#'   \code{scRNAseq_for_GRN} from which gene-level expression features are
#'   computed.
#' @param interest_cell_type_tau Named numeric vector. Output of \code{\link{get_cell_type_tau}}. Cell-type-specific
#'   binarization thresholds, used to shift normalized positive strengths above
#'   \code{tau}.
#' @param topn_max Positive integer. Maximum number of models averaged per
#' pair.
#' Defaults to 25.
#' @param topn_min Positive integer. Minimum number of models averaged per
#' pair.
#' Defaults to 5.
#' @param topn_ratio Numeric between 0 and 1. Number of models averaged as a
#' proportion of available trained models.
#' Defaults to 0.1.
#' @param fix_factor Numeric in (0, 1\]. Factor applied to the weights of
#' boundary strengths, with 1 applying no down-weighting; must be strictly
#' greater than 0.
#' Defaults to 1.0.
#' @param min_scale Numeric in (0, 1\]. Multiplier on the smallest positive
#' value when setting the min-max lower bound, applied to strength and loss
#' normalization.
#' For strengths, \code{min_scale = 1} maps the weakest positive edge of a TF
#' to \code{tau}, the strongest to just below 1, while absent pairs reset to 0;
#' values below 1 keep the weakest edge slightly above \code{tau}.
#' Defaults to 0.99.
#' @param each Positive integer. Maximum number of pairs per parallel chunk.
#' Defaults to 1000000.
#' @param ncores Positive integer. Number of cores for parallel processing via
#'   \code{\link[parallel]{mclapply}}.
#' Defaults to 1; set to 1 on Windows.
#'
#' @return The input list with a \code{pred_summary} frame per cell type.
#' Its columns are:
#' \itemize{
#'   \item \code{TF}, \code{TG}: the pair identifiers.
#'   \item \code{theta_p_branchN}, \code{loss_branchN}: branch-level averaged
#'     strength and loss.
#'   \item \code{FG.theta_p}, \code{FG.Loss}: cell-count weighted strength and
#'     loss across branches.
#'   \item \code{FG.theta_p_norm}, \code{FG.theta_p_norm.R}: within-TF
#'     normalized strength and its rank.
#'   \item \code{FG.Loss_norm}, \code{FG.Loss_norm.R}: within-TF normalized
#'     loss and its rank.
#'   \item \code{FG.theta_i}, \code{FG.theta_i.R}: observed strength and its
#'     normalized rank, with rank 0 for \code{theta_i = 0}.
#'   \item \code{F.*} and \code{G.*}: TF and TG expression features and their
#'     normalized ranks, including non-zero ratio, non-zero average and overall
#'     average.
#' }
#'
#' @details For each cell type, the function proceeds as follows:
#' \enumerate{
#'   \item Computes shared gene expression features once from the normalized
#'     scRNA-seq data.
#'   \item Splits the pairs of each branch into chunks.
#'   \item For each pair, selects \code{topn} lowest-loss models, where
#'     \code{topn} is clamped to
#'     \code{round(n_train * topn_ratio)} within \code{topn_min} and
#'     \code{topn_max}.
#'   \item Assigns triangular weights from loss, from 1 for the smallest loss to
#'     0 for the largest, or equal weights when all selected losses tie;
#'     weights of boundary strengths below 0.0101 or above 0.9899 are
#'     multiplied by \code{fix_factor}.
#'   \item Takes the weighted mean of the selected strengths and losses.
#'   \item Merges chunks, suffixes columns by branch, and removes pairs with
#'     strength below 0.011 or infinite loss.
#'   \item Integrates branches by cell-count weights.
#'   \item Normalizes strengths and losses within each TF: the lower bound is
#'     the smallest positive value times \code{min_scale}, the upper bound is
#'     the largest value plus 0.0001; positive strengths are shifted above
#'     \code{tau}, while absent pairs reset to 0, and losses are scaled without
#'     the shift.
#'   \item Computes normalized ranks, fills the observed \code{theta_i} and its
#'     rank, and attaches the gene features.
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and switches to the \code{5 Infer GRN}
#'     directory, and restores the original working directory on exit.
#'   \item \code{fix_factor} must be in (0, 1]; set it below 1 to discount
#'     boundary predictions.
#'   \item Boundary strengths near 0.01 or 0.99 reflect the optimization bounds
#'     and are treated as less reliable.
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
#' interest_cell_type_pGRN = select_best_prediction_old2(
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
select_best_prediction_old2 = function(
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
  message("The current working directory has been switched to: ", base::getwd(), ".")
  # Restore the working directory on normal exit or error.
  base::on.exit({
    base::setwd(original_dir)
    message("The working directory has been restored to: ", base::getwd(), ".")
  }, add = TRUE)

  ### Shared gene expression features, computed once for all cell types.
  message("Computing shared gene-level expression features.")
  scRNAseq = Seurat::NormalizeData(Basic_Info$scRNAseq_for_GRN, verbose = FALSE)
  rna_data = SeuratObject::GetAssayData(scRNAseq, layer = "data")
  nonzero_n = Matrix::rowSums(rna_data > 0)
  gene_features = base::data.frame(
    nonzero_freq = nonzero_n / base::ncol(rna_data),
    nonzero_mean = base::ifelse(nonzero_n > 0, Matrix::rowSums(rna_data) / nonzero_n, 0),
    overall_mean = Matrix::rowMeans(rna_data)
  )
  gene_features[base::paste0(base::colnames(gene_features), "-norm_rank")] =
    base::lapply(gene_features, function(x) {
      r = base::rank(x, ties.method = "min", na.last = "keep")
      r / base::max(r, na.rm = TRUE)
    })

  ### Process.
  for (cell_type in base::names(interest_cell_type_pGRN)) {
    branch_index = base::seq_along(
      interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]]
    )
    tau = base::as.numeric(interest_cell_type_tau[cell_type])

    for (n in branch_index) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      branch_result = interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]]
      re_pred = base::rownames(branch_result[["theta_p"]])
      N_pred_all = base::length(re_pred)
      message("Branch ", n, " has ", N_pred_all, " candidate pairs.")

      # Split pairs into chunks.
      n_part = base::ceiling(N_pred_all / each)
      re_pred_list = base::list()
      for (s in base::seq_len(n_part)) {
        if (s != n_part) {
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):(each * s)]
        } else {
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):N_pred_all]
        }
      }

      # Per-chunk inputs.
      pGRN_list = base::list()
      for (s in base::seq_len(n_part)) {
        chunk_rows = re_pred_list[[s]]
        pGRN_list[[s]] = base::list(
          theta_p = branch_result[["theta_p"]][chunk_rows, , drop = FALSE],
          loss = branch_result[["loss"]][chunk_rows, , drop = FALSE]
        )
      }

      # Number of models averaged, clamped to [topn_min, topn_max] and n_train.
      n_train = base::ncol(branch_result[["theta_p"]])
      topn = base::min(base::max(base::round(n_train * topn_ratio), topn_min),
                       topn_max, n_train)

      # Combine the topn lowest-loss models for every pair in one chunk.
      process_data_slice = function(pGRN_info, topn, fix_factor) {
        theta_frame = pGRN_info[["theta_p"]]
        loss_frame = pGRN_info[["loss"]]
        for (re in base::rownames(theta_frame)) {
          theta_vec = base::unlist(theta_frame[re, , drop = FALSE])
          loss_vec = base::unlist(loss_frame[re, , drop = FALSE])
          best_result = base::data.frame(base::matrix(
            NA, nrow = 1, ncol = 2,
            dimnames = base::list(re, c("theta_p", "loss"))
          ))

          loss_topn = base::sort(loss_vec)[base::seq_len(topn)]
          pattern = base::names(loss_topn)
          theta_topn = theta_vec[pattern]

          # Triangular weights from loss; equal weights when all losses tie.
          loss_range = base::diff(base::range(loss_topn))
          if (loss_range == 0) {
            w = base::rep(1, topn)
          } else {
            w = (base::max(loss_topn) - loss_topn) / loss_range
          }
          # Down-weight boundary strengths.
          boundary = theta_topn < 0.0101 | theta_topn > 0.9899
          w_fix = w
          w_fix[boundary] = w_fix[boundary] * fix_factor

          best_result[re, "theta_p"] = stats::weighted.mean(theta_topn, w_fix)
          best_result[re, "loss"] = stats::weighted.mean(loss_topn, w_fix)
          pGRN_info[[re]] = base::list(best_result = best_result)
        }

        # Merge pair results into a chunk-level frame.
        pair_results = pGRN_info[3:base::length(pGRN_info)]
        pGRN_info[["best_result"]] = base::as.data.frame(data.table::rbindlist(base::lapply(
          pair_results, function(x) x[["best_result"]]
        )))
        base::rownames(pGRN_info[["best_result"]]) = base::names(pair_results)
        # Remove boundary-invalid strengths and infinite losses.
        pGRN_info[["best_result"]] = pGRN_info[["best_result"]][
          pGRN_info[["best_result"]]$theta_p >= 0.011 &
            !base::is.infinite(pGRN_info[["best_result"]]$loss), , drop = FALSE]
        return(pGRN_info[["best_result"]])
      }

      prediction_list = parallel::mclapply(
        X = base::seq_len(n_part),
        FUN = function(s) {
          process_data_slice(pGRN_list[[s]], topn = topn, fix_factor = fix_factor)
        },
        mc.cores = ncores
      )
      base::gc()

      # Bind chunks and suffix the columns with the branch index.
      message("Summarizing best_result of this branch.")
      branch_best = base::do.call(base::rbind, prediction_list)
      base::colnames(branch_best) = base::paste0(base::colnames(branch_best), "_branch", n)
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]] =
        branch_best
    }

    # Build the integrated cell-type summary.
    message("Integrating branch results.")
    branch_blocks = interest_cell_type_pGRN[[cell_type]][base::paste0("branch", branch_index)]
    rows = base::sort(base::unique(base::Reduce(base::c, base::lapply(
      branch_blocks, function(x) base::rownames(x[["best_result"]])
    ))))
    branch_cols = base::Reduce(base::c, base::lapply(
      branch_blocks, function(x) base::colnames(x[["best_result"]])
    ))
    cols = base::c(
      "TF", "TG", branch_cols,
      "FG.theta_p", "FG.Loss",
      "FG.theta_p_norm", "FG.theta_p_norm.R", "FG.Loss_norm", "FG.Loss_norm.R",
      "FG.theta_i", "FG.theta_i.R",
      "F.NZER", "F.NZER.R", "F.NZEA", "F.NZEA.R", "F.EA", "F.EA.R",
      "G.NZER", "G.NZER.R", "G.NZEA", "G.NZEA.R", "G.EA", "G.EA.R"
    )
    pred_summary = base::data.frame(base::matrix(
      NA, nrow = base::length(rows), ncol = base::length(cols),
      dimnames = base::list(rows, cols)
    ))
    pred_summary[, c("TF", "TG")] = base::do.call(
      base::rbind, base::strsplit(base::rownames(pred_summary), "_to_")
    )
    for (n in branch_index) {
      branch_best = branch_blocks[[n]][["best_result"]]
      pred_summary[base::rownames(branch_best), base::colnames(branch_best)] = branch_best
    }

    # Integrate strength and loss by branch cell counts.
    n_branch = base::length(branch_index)
    weights = base::unlist(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]])
    if (n_branch > 1) {
      pred_summary$FG.theta_p = base::apply(
        pred_summary[, base::grep("^theta_p_branch([0-9]{1,2})$", base::names(pred_summary),
                                  value = TRUE)],
        1, function(x) {
          valid = base::which(!base::is.na(x))
          if (base::length(valid) == 0) 0
          else base::sum(x[valid] * weights[valid]) / base::sum(weights[valid])
        }
      )
      pred_summary$FG.Loss = base::apply(
        pred_summary[, base::grep("^loss_branch([0-9]{1,2})$", base::names(pred_summary),
                                  value = TRUE)],
        1, function(x) {
          valid = base::which(!base::is.na(x))
          if (base::length(valid) == 0) Inf
          else base::sum(x[valid] * weights[valid]) / base::sum(weights[valid])
        }
      )
    } else if (n_branch == 1) {
      pred_summary$FG.theta_p = pred_summary$theta_p_branch1
      pred_summary$FG.Loss = pred_summary$loss_branch1
    } else {
      base::stop("The number of branches is less than 1.")
    }

    # Within-TF min-max normalization of predicted strengths, tau-shifted.
    message("Normalizing predicted strengths within each TF.")
    pGRN = base::as.data.frame.matrix(stats::xtabs(FG.theta_p ~ TG + TF, data = pred_summary))
    pGRN_norm = pGRN
    for (tf in base::colnames(pGRN_norm)) {
      tf_tgs = pGRN_norm[, tf, drop = FALSE]
      base::colnames(tf_tgs)[1] = "theta_p"
      if (base::max(tf_tgs$theta_p) == 0) {
        tf_tgs$theta_p_norm = 0
      } else {
        max_before = base::max(tf_tgs$theta_p) + 0.0001
        min_before = base::min(tf_tgs$theta_p[tf_tgs$theta_p > 0]) * min_scale
        tf_tgs$theta_p_norm = tau + (1 - tau) *
          (tf_tgs$theta_p - min_before) / (max_before - min_before)
        tf_tgs$theta_p_norm[tf_tgs$theta_p_norm < 0] = 0
      }
      pGRN_norm[, tf] = tf_tgs[base::rownames(pGRN_norm), "theta_p_norm"]
    }
    pred_summary$FG.theta_p_norm = base::mapply(
      function(tf, tg) pGRN_norm[tg, tf],
      pred_summary$TF, pred_summary$TG, USE.NAMES = FALSE
    )
    positive_theta_p = base::sort(pred_summary$FG.theta_p_norm)
    positive_theta_p = positive_theta_p[positive_theta_p > 0]
    if (base::length(positive_theta_p) > 0) {
      rank_df = base::unique(base::data.frame(
        value = positive_theta_p,
        norm_rank = base::rank(positive_theta_p, ties.method = "min") /
          base::length(positive_theta_p)))
      pred_summary$FG.theta_p_norm.R =
        rank_df$norm_rank[base::match(pred_summary$FG.theta_p_norm, rank_df$value)]
    } else {
      pred_summary$FG.theta_p_norm.R = base::rep(NA_real_, base::nrow(pred_summary))
    }
    pred_summary$FG.theta_p_norm.R[pred_summary$FG.theta_p_norm <= 0] = 0

    # Within-TF min-max normalization of losses, without the tau shift.
    message("Normalizing losses within each TF.")
    pLoss = base::as.data.frame.matrix(stats::xtabs(FG.Loss ~ TG + TF, data = pred_summary))
    pLoss_norm = pLoss
    for (tf in base::colnames(pLoss_norm)) {
      tf_tgs = pLoss_norm[, tf, drop = FALSE]
      base::colnames(tf_tgs)[1] = "Loss"
      if (base::max(tf_tgs$Loss) == 0) {
        tf_tgs$Loss_norm = 0
      } else {
        max_before = base::max(tf_tgs$Loss) + 0.0001
        min_before = base::min(tf_tgs$Loss[tf_tgs$Loss > 0]) * min_scale
        tf_tgs$Loss_norm = (tf_tgs$Loss - min_before) / (max_before - min_before)
        tf_tgs$Loss_norm[tf_tgs$Loss_norm < 0] = 0
      }
      pLoss_norm[, tf] = tf_tgs[base::rownames(pLoss_norm), "Loss_norm"]
    }
    pred_summary$FG.Loss_norm = base::mapply(
      function(tf, tg) pLoss_norm[tg, tf],
      pred_summary$TF, pred_summary$TG, USE.NAMES = FALSE
    )
    positive_loss = base::sort(pred_summary$FG.Loss_norm)
    positive_loss = positive_loss[positive_loss > 0]
    if (base::length(positive_loss) > 0) {
      rank_df = base::unique(base::data.frame(
        value = positive_loss,
        norm_rank = base::rank(positive_loss, ties.method = "min") /
          base::length(positive_loss)
      ))
      pred_summary$FG.Loss_norm.R =
        rank_df$norm_rank[base::match(pred_summary$FG.Loss_norm, rank_df$value)]
    } else {
      pred_summary$FG.Loss_norm.R = base::rep(NA_real_, base::nrow(pred_summary))
    }

    # Observed theta_i from the iGRN matrix.
    message("Attaching observed theta_i and its rank.")
    iGRN = interest_cell_type_iGRN[[cell_type]]
    pred_summary$FG.theta_i = base::mapply(
      function(tf, tg) iGRN[tg, tf],
      pred_summary$TF, pred_summary$TG, USE.NAMES = FALSE
    )
    positive_theta_i = base::sort(base::unlist(iGRN))
    positive_theta_i = positive_theta_i[positive_theta_i > 0]
    if (base::length(positive_theta_i) > 0) {
      rank_df = base::unique(base::data.frame(
        value = positive_theta_i,
        norm_rank = base::rank(positive_theta_i, ties.method = "min") /
          base::length(positive_theta_i)))
      pred_summary$FG.theta_i.R =
        rank_df$norm_rank[base::match(pred_summary$FG.theta_i, rank_df$value)]
    } else {
      pred_summary$FG.theta_i.R = base::rep(NA_real_, base::nrow(pred_summary))
    }
    # Pairs with theta_i = 0 keep rank 0.
    pred_summary$FG.theta_i.R[pred_summary$FG.theta_i == 0] = 0

    # Attach the precomputed TF and TG gene features.
    message("Attaching gene-level expression features.")
    attach_feature = function(gene_col, feature, rank = FALSE) {
      col_name = base::ifelse(rank, base::paste0(feature, "-norm_rank"), feature)
      gene_features[pred_summary[[gene_col]], col_name]
    }
    pred_summary$F.NZER = attach_feature("TF", "nonzero_freq")
    pred_summary$F.NZER.R = attach_feature("TF", "nonzero_freq", rank = TRUE)
    pred_summary$G.NZER = attach_feature("TG", "nonzero_freq")
    pred_summary$G.NZER.R = attach_feature("TG", "nonzero_freq", rank = TRUE)
    pred_summary$F.NZEA = attach_feature("TF", "nonzero_mean")
    pred_summary$F.NZEA.R = attach_feature("TF", "nonzero_mean", rank = TRUE)
    pred_summary$G.NZEA = attach_feature("TG", "nonzero_mean")
    pred_summary$G.NZEA.R = attach_feature("TG", "nonzero_mean", rank = TRUE)
    pred_summary$F.EA = attach_feature("TF", "overall_mean")
    pred_summary$F.EA.R = attach_feature("TF", "overall_mean", rank = TRUE)
    pred_summary$G.EA = attach_feature("TG", "overall_mean")
    pred_summary$G.EA.R = attach_feature("TG", "overall_mean", rank = TRUE)

    interest_cell_type_pGRN[[cell_type]][["pred_summary"]] = pred_summary
  }

  ### End.
  t_end = base::Sys.time()
  message("Finish: Selecting best predictions. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start), ".")
  return(interest_cell_type_pGRN)
}
