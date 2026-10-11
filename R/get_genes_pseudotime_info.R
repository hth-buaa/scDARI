#' @title Extract Pseudotemporal Gene Expression Profiles along Pseudotime
#'
#' @description For each cell type and pseudotime branch,
#' this function calculates transcription factor (TF) and target gene (TG) expression importance scores from scRNA-seq data,
#' filters genes and cells, orders the retained cells by pseudotime,
#' and extracts pseudotime-series matrices of TF expression (TFE) and TG expression (TGE),
#' ready for downstream Gene Regulatory Network (GRN) inference.
#'
#' @param interest_cell_type_Branches List. Output from \code{\link{order_pseudotime_and_divide_branches}};
#' each cell type holds one ordered table of cells per branch.
#' @param interest_cell_type_data List. Output from \code{\link{get_interest_cell_type_data}}.
#' @param alpha_TF Numeric in \[0, 1\], default \code{0.63}. Quantile cutoff for TF filtering.
#' TFs whose importance score (TFIS) is strictly above this quantile are retained.
#' @param alpha_TG Numeric in \[0, 1\], default \code{0.69}. Quantile cutoff for TG filtering.
#' TGs whose importance score (TGIS) is strictly above this quantile are retained.
#' @param alpha_cell Numeric in \[0, 1\], default 1 (no filtering).
#' Maximum allowed zero-expression rate (computed on the scRNA-seq counts of the retained genes);
#' cells above this value are dropped.
#' @param w_TFNER Numeric in \[-0.2, 1\], default \code{0.69}. Weight of the max-normalized TF non-zero expression ratio (NER) in TFIS.
#' @param w_TFNERR Numeric in \[-0.2, 1\], default \code{1.00}. Weight of the rank-normalized TF non-zero expression ratio (NERR) in TFIS.
#' @param w_TFNEA Numeric in \[-0.2, 1\], default \code{0.67}. Weight of the max-normalized TF mean expression in non-zero cells (NEA) in TFIS.
#' @param w_TFNEAR Numeric in \[-0.2, 1\], default \code{0.00}. Weight of the rank-normalized TF non-zero expression average (NEAR) in TFIS.
#' @param w_TGNER Numeric in \[-0.2, 1\], default \code{0.30}. Weight of the max-normalized TG non-zero expression ratio (NER) in TGIS.
#' @param w_TGNERR Numeric in \[-0.2, 1\], default \code{0.88}. Weight of the rank-normalized TG non-zero expression ratio (NERR) in TGIS.
#' @param w_TGNEA Numeric in \[-0.2, 1\], default \code{0.23}. Weight of the max-normalized TG mean expression in non-zero cells (NEA) in TGIS.
#' @param w_TGNEAR Numeric in \[-0.2, 1\], default \code{0.12}. Weight of the rank-normalized TG non-zero expression average (NEAR) in TGIS.
#' @param ncores Positive integer. Number of cores for \code{\link[parallel]{mclapply}}.
#' Defaults to \code{1}. Fork-based parallelism is supported only on Unix-like systems;
#' it is forced to \code{1} on Windows.
#' @param sparsity_cut Nested list or \code{NULL} (default). Optional post-filtering by detection (non-zero) frequency.
#' Top-level names must match cell types;
#' each element is a data frame with rows \code{TFE}, \code{TGE} and columns \code{branch1}, \code{branch2}, ....
#' Values must be in \[0, 1\] and are the minimum required non-zero frequency per branch.
#' An invalid format is reset to \code{NULL} with a message.
#'
#' @details
#' For each cell type and branch, the workflow proceeds as follows:
#' \enumerate{
#'   \item Extract the scRNA-seq count matrix for cells assigned to the current branch.
#'   \item Compute four features separately for TFs and TGs: the max-normalized non-zero expression
#'         ratio (NER) and non-zero expression average (NEA), together with their rank-normalized
#'         versions (NERR and NEAR; 0 = lowest, 1 = highest).
#'   \item Compute the weighted importance scores TFIS and TGIS from the four expression features.
#'   \item Filter genes by score quantiles (\code{alpha_TF}, \code{alpha_TG}) and cells by the
#'         scRNA-seq zero-expression rate (\code{alpha_cell}).
#'   \item Map cells to integer pseudotime indices (\code{k~}) and order columns in ascending pseudotime.
#'   \item Normalize expression with \code{\link[Seurat]{NormalizeData}} and extract the TFE and TGE matrices.
#'   \item Optionally apply \code{sparsity_cut} to re-filter TFE/TGE by detection frequency and
#'         synchronize the retained TF/TG lists.
#'   \item Return the structured nested list.
#' }
#'
#' @note
#' \itemize{
#'   \item Alpha thresholds and feature weights directly control the number of retained genes and cells.
#'   \item Gene filtering uses a strict greater-than comparison against a quantile: \code{alpha_* = 1}
#'         retains no gene, and even \code{alpha_* = 0} drops the single lowest-scoring gene.
#'   \item Parallel computing via \code{\link[parallel]{mclapply}} works only on Unix-like systems;
#'         \code{ncores} is forced to 1 on Windows.
#'   \item The function creates the \code{2.2 Data Processing - Cell Grouping} directory and temporarily
#'         switches the working directory, restoring the original path on exit.
#'   \item This function uses scRNA-seq only and does not produce a chromatin activity (TGA) matrix.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_Branches = base::readRDS(base::paste0(
#'   "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'   "/interest_cell_type_Branches.rds"
#' ))
#' interest_cell_type_data = base::readRDS(base::paste0(
#'   "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'   "/interest_cell_type_data.rds"
#' ))
#' alpha_TF = 0.7; alpha_TG = 0.7; alpha_cell = 0.9
#' w_TFNER = 1; w_TFNERR = 0.6; w_TFNEA = 0.8; w_TFNEAR = 0.4
#' w_TGNER = 0.3; w_TGNERR = 0.1; w_TGNEA = 0.3; w_TGNEAR = 0.1
#' ncores = parallel::detectCores() - 1 # Linux
#' # ncores = 1                       # Windows
#'
#' # Example 1: without sparsity_cut
#' interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
#'   interest_cell_type_Branches = interest_cell_type_Branches,
#'   interest_cell_type_data = interest_cell_type_data,
#'   alpha_TF = alpha_TF, alpha_TG = alpha_TG, alpha_cell = alpha_cell,
#'   w_TFNER = w_TFNER, w_TFNERR = w_TFNERR, w_TFNEA = w_TFNEA, w_TFNEAR = w_TFNEAR,
#'   w_TGNER = w_TGNER, w_TGNERR = w_TGNERR, w_TGNEA = w_TGNEA, w_TGNEAR = w_TGNEAR,
#'   ncores = ncores
#' )
#'
#' # Example 2: with sparsity_cut (rows are TFE and TGE only)
#' sparsity_cut = list(
#'   cell_name1 = data.frame(
#'     branch1 = c(TFE = 0.00, TGE = 0.03),
#'     branch2 = c(TFE = 0.00, TGE = 0.03),
#'     branch3 = c(TFE = 0.00, TGE = 0.03),
#'     row.names = c("TFE", "TGE")
#'   ),
#'   cell_name2 = data.frame(
#'     branch1 = c(TFE = 0.00, TGE = 0.03),
#'     branch2 = c(TFE = 0.00, TGE = 0.03),
#'     branch3 = c(TFE = 0.00, TGE = 0.03),
#'     row.names = c("TFE", "TGE")
#'   )
#' )
#' interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
#'   interest_cell_type_Branches = interest_cell_type_Branches,
#'   interest_cell_type_data = interest_cell_type_data,
#'   alpha_TF = alpha_TF, alpha_TG = alpha_TG, alpha_cell = alpha_cell,
#'   w_TFNER = w_TFNER, w_TFNERR = w_TFNERR, w_TFNEA = w_TFNEA, w_TFNEAR = w_TFNEAR,
#'   w_TGNER = w_TGNER, w_TGNERR = w_TGNERR, w_TGNEA = w_TGNEA, w_TGNEAR = w_TGNEAR,
#'   ncores = ncores,
#'   sparsity_cut = sparsity_cut
#' )
#' # Rely solely on sparsity_cut (disable quantile/cell filtering as far as possible)
#' interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
#'   interest_cell_type_Branches = interest_cell_type_Branches,
#'   interest_cell_type_data = interest_cell_type_data,
#'   alpha_TF = 0, alpha_TG = 0, alpha_cell = 1,
#'   ncores = ncores,
#'   sparsity_cut = sparsity_cut
#' )
#'
#' # Check the minimum detection (non-zero) frequency across TFE/TGE matrices
#' min_nonzero_freq = 1.0
#' for (ct in base::names(interest_cell_type_genes_pseudotime_info)) {
#'   for (slot in c("Branches_TFE", "Branches_TGE")) {
#'     for (df in interest_cell_type_genes_pseudotime_info[[ct]][[slot]]) {
#'       if (base::nrow(df) > 0 && base::ncol(df) > 0) {
#'         min_nonzero_freq = base::min(min_nonzero_freq, base::rowSums(df != 0) / base::ncol(df))
#'       }
#'     }
#'   }
#' }
#' base::cat("min_nonzero_freq: ", min_nonzero_freq, "\n")
#' base::saveRDS(
#'   interest_cell_type_genes_pseudotime_info,
#'   file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#' }
get_genes_pseudotime_info = function(
    interest_cell_type_Branches,
    interest_cell_type_data,
    alpha_TF = 0.63, alpha_TG = 0.69, alpha_cell = 1.00,
    w_TFNER = 0.69, w_TFNERR = 1.00, w_TFNEA = 0.67, w_TFNEAR = 0.00,
    w_TGNER = 0.30, w_TGNERR = 0.88, w_TGNEA = 0.23, w_TGNEAR = 0.12,
    ncores = 1,
    sparsity_cut = NULL
) {
  ### Validate parameters.
  # Alpha thresholds must lie within [0, 1].
  if (alpha_TF < 0 || alpha_TF > 1 || alpha_TG < 0 || alpha_TG > 1 || alpha_cell < 0 || alpha_cell > 1) {
    base::stop("alpha_TF, alpha_TG, alpha_cell must be between 0 and 1.")
  }

  # The two inputs must describe the same set of cell types.
  cell_types = base::names(interest_cell_type_Branches)
  missing_in_data = base::setdiff(cell_types, base::names(interest_cell_type_data))
  if (base::length(missing_in_data) > 0) {
    base::stop("interest_cell_type_data is missing cell type(s): ",
               base::paste(missing_in_data, collapse = ", "), ".")
  }
  required_elements = c("interest_scRNAseq_for_GRN", "interest_gene_for_GRN",
                        "interest_TGs", "interest_TFs")
  for (cell_type in cell_types) {
    if (!base::all(required_elements %in% base::names(interest_cell_type_data[[cell_type]]))) {
      base::stop("interest_cell_type_data[['", cell_type,
                 "']] is malformed; use the output of get_interest_cell_type_data().")
    }
  }

  # Validate sparsity_cut; only the expression matrices TFE and TGE exist.
  if (!base::is.null(sparsity_cut)) {
    message("Validating sparsity_cut.")
    valid = TRUE
    if (!base::is.list(sparsity_cut) || base::is.null(base::names(sparsity_cut))) {
      valid = FALSE
    } else if (!base::setequal(base::names(sparsity_cut), cell_types)) {
      valid = FALSE
    } else {
      for (cell_type in cell_types) {
        x = sparsity_cut[[cell_type]]
        expected_rows = c("TFE", "TGE")
        expected_cols = base::paste0("branch", base::seq_along(interest_cell_type_Branches[[cell_type]]))
        x_mat = base::tryCatch(base::as.matrix(x), error = function(e) NULL)
        if (!base::is.data.frame(x) ||
            !base::setequal(base::rownames(x), expected_rows) ||
            !base::setequal(base::colnames(x), expected_cols) ||
            base::is.null(x_mat) || !base::is.numeric(x_mat) ||
            base::anyNA(x_mat) || base::any(x_mat < 0) || base::any(x_mat > 1)) {
          valid = FALSE
          break
        }
      }
    }
    if (!valid) {
      message("Invalid sparsity_cut format; resetting to NULL.")
      sparsity_cut = NULL
    }
  }

  # Fork-based parallelism is unavailable on Windows.
  if (.Platform$OS.type == "windows") {
    ncores = 1
    message("Parallel computing is not supported on Windows; ncores forced to 1.")
  }

  ### Start.
  t_start = base::Sys.time()
  message("Run: Obtaining Gene Expression/Activity Vectors along Pseudotime. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.2 Data Processing - Cell Grouping"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("Working directory set to: ", base::getwd(), ".")

  ### Define the per-cell-type processing function.
  get_genes_pseudotime_info_process = function(
    Branches,
    interest_scRNAseq_for_GRN,
    interest_gene_for_GRN, interest_TGs, interest_TFs,
    alpha_TF, w_TFNER, w_TFNERR, w_TFNEA, w_TFNEAR,
    alpha_TG, w_TGNER, w_TGNERR, w_TGNEA, w_TGNEAR,
    alpha_cell
  ) {
    # Per-branch storage.
    Branches_RNA = base::list()
    Branches_gene_for_GRN = base::list()
    Branches_TGs = base::list()
    Branches_TFs = base::list()
    Branches_n_cell = base::list()

    # Feature helper: max-normalized non-zero ratio/mean, their ranks, and the weighted importance score.
    calculate_features = function(exp_mat,  # gene expression matrix (genes x cells)
                                  genes,    # gene name vector (TF or TG)
                                  weights,  # 4 feature weights (NER, NERR, NEA, NEAR)
                                  prefix) { # column name prefix ("TF" or "TG")
      # Non-zero expression ratio (fraction of expressing cells).
      ner = base::rowMeans(exp_mat != 0)

      # Non-zero expression average (mean over expressing cells only).
      exp_na = exp_mat; exp_na[exp_na == 0] = NA
      nea = base::rowMeans(exp_na, na.rm = TRUE); rm(exp_na)

      # Fill NA with 0
      ner[base::is.na(ner)] = 0
      nea[base::is.na(nea)] = 0

      # Max-normalize both raw features to [0, 1].
      max_ner = base::max(ner)
      if (max_ner > 0) {
        ner = ner / max_ner
      } else {
        stop("Cell type '", cell_type, "': NER is zero for every gene, normalization failed.")
      }

      max_nea = base::max(nea)
      if (max_nea > 0) {
        nea = nea / max_nea
      } else {
        stop("Cell type '", cell_type, "': NEA is zero for every gene, normalization failed.")
      }

      # Rank-normalize (0 = lowest, 1 = highest).
      n = length(genes)
      nerr = base::rank(ner) / n
      near = base::rank(nea) / n

      # Weighted importance score.
      feat_mat = base::cbind(ner, nerr, nea, near)
      is_score = base::as.vector(feat_mat %*% weights)

      # Output feature data frame
      df = data.frame(
        feat_mat,
        IS = is_score,
        row.names = genes,
        check.names = FALSE
      )
      base::colnames(df) = paste0(prefix, c("NER", "NERR", "NEA", "NEAR", "IS"))

      return(df)
    }

    for (i in base::seq_along(Branches)) {
      message("[", cell_type, "] Branch ", i, ": extracting RNA counts.")
      # Restrict to the cells assigned to this branch.
      Branches_RNA[[i]] = interest_scRNAseq_for_GRN[, base::rownames(Branches[[i]])]
      rna_assay = SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "counts")

      message("[", cell_type, "] Branch ", i, ": computing TF/TG importance scores.")
      TF_features = calculate_features(
        exp_mat = base::as.matrix(rna_assay[interest_TFs, ]),
        genes = interest_TFs,
        weights = c(w_TFNER, w_TFNERR, w_TFNEA, w_TFNEAR),
        prefix = "TF"
      )
      TG_features = calculate_features(
        exp_mat = base::as.matrix(rna_assay[interest_TGs, ]),
        genes = interest_TGs,
        weights = c(w_TGNER, w_TGNERR, w_TGNEA, w_TGNEAR),
        prefix = "TG"
      )

      message("[", cell_type, "] Branch ", i, ": filtering genes by score quantile.")
      TFs_filter = base::rownames(TF_features)[TF_features$TFIS > stats::quantile(TF_features$TFIS, alpha_TF)]
      TGs_filter = base::rownames(TG_features)[TG_features$TGIS > stats::quantile(TG_features$TGIS, alpha_TG)]
      genes_filter = base::sort(base::unique(c(TFs_filter, TGs_filter)))
      if (base::length(genes_filter) == 0) {
        base::stop("[", cell_type, "] Branch ", i,
                   ": no gene passed the alpha_TF/alpha_TG quantile filters.")
      }
      Branches_RNA[[i]] = subset(Branches_RNA[[i]], features = genes_filter)

      message("[", cell_type, "] Branch ", i, ": filtering cells by zero-expression rate <= ", alpha_cell, ".")
      rna_counts = SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "counts")
      zero_freq_RNA = Matrix::colMeans(rna_counts == 0)
      cells_filter = base::names(zero_freq_RNA[zero_freq_RNA <= alpha_cell])
      Branches_RNA[[i]] = subset(Branches_RNA[[i]], cells = cells_filter)

      message("[", cell_type, "] Branch ", i, ": harmonizing TF/TG gene sets.")
      Branches_TFs[[i]] = TFs_filter[TFs_filter %in% base::rownames(Branches_RNA[[i]])]
      Branches_TGs[[i]] = TGs_filter[TGs_filter %in% base::rownames(Branches_RNA[[i]])]
      Branches_gene_for_GRN[[i]] = base::sort(base::union(Branches_TGs[[i]], Branches_TFs[[i]]))
      Branches_RNA[[i]] = Branches_RNA[[i]][Branches_gene_for_GRN[[i]], ]
      Branches_n_cell[[i]] = base::ncol(Branches_RNA[[i]])
    }
    rm(calculate_features)

    ### Map cells to integer pseudotime indices.
    message("[", cell_type, "] Mapping cells to pseudotime indices.")
    for (i in base::seq_along(Branches)) {
      message("Building cell-to-pseudotime map for branch ", i, ".")
      Branches[[i]]["cell"] = base::rownames(Branches[[i]])
      Branches[[i]]["k~"] = 1:base::nrow(Branches[[i]])
      name_mapping = stats::setNames(Branches[[i]]$`k~`, Branches[[i]]$cell)

      message("Renaming cells by pseudotime index for branch ", i, ".")
      base::colnames(Branches_RNA[[i]]) = name_mapping[base::colnames(Branches_RNA[[i]])]
    }

    ### Normalize expression (build the log-normalized 'data' layer).
    message("[", cell_type, "] Normalizing expression matrices.")
    for (i in base::seq_along(Branches)) {
      message("Normalizing scRNA-seq for branch ", i, ".")
      Branches_RNA[[i]] = Seurat::NormalizeData(Branches_RNA[[i]])
    }

    ### Extract ordered TFE / TGE pseudotime matrices.
    message("[", cell_type, "] Extracting TFE/TGE pseudotime matrices.")
    Branches_TFE = base::list()
    Branches_TGE = base::list()
    for (i in base::seq_along(Branches)) {
      norm_data = SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "data")

      # TF expression (TFE), ordered by ascending pseudotime index.
      Branches_TFE[[i]] = base::as.data.frame(norm_data)[Branches_TFs[[i]], ]
      Branches_TFE[[i]] = Branches_TFE[[i]][, base::order(base::as.integer(base::colnames(Branches_TFE[[i]])))]

      # TG expression (TGE), ordered by ascending pseudotime index.
      Branches_TGE[[i]] = base::as.data.frame(norm_data)[Branches_TGs[[i]], ]
      Branches_TGE[[i]] = Branches_TGE[[i]][, base::order(base::as.integer(base::colnames(Branches_TGE[[i]])))]
    }

    return(base::list("Branches_TFE" = Branches_TFE,
                      "Branches_TGE" = Branches_TGE,
                      "Branches_n_cell" = Branches_n_cell,
                      "Branches" = Branches,
                      "Branches_TFs" = Branches_TFs,
                      "Branches_TGs" = Branches_TGs))
  }

  ### Run all cell types (optionally in parallel).
  interest_cell_type_genes_pseudotime_info = parallel::mclapply(
    X = base::names(interest_cell_type_Branches),
    FUN = function(cell_type) {
      get_genes_pseudotime_info_process(
        Branches = interest_cell_type_Branches[[cell_type]],
        interest_scRNAseq_for_GRN = interest_cell_type_data[[cell_type]][["interest_scRNAseq_for_GRN"]],
        interest_gene_for_GRN = interest_cell_type_data[[cell_type]][["interest_gene_for_GRN"]],
        interest_TFs = interest_cell_type_data[[cell_type]][["interest_TFs"]],
        interest_TGs = interest_cell_type_data[[cell_type]][["interest_TGs"]],
        alpha_TF = alpha_TF,
        w_TFNER = w_TFNER, w_TFNERR = w_TFNERR, w_TFNEA = w_TFNEA, w_TFNEAR = w_TFNEAR,
        alpha_TG = alpha_TG,
        w_TGNER = w_TGNER, w_TGNERR = w_TGNERR, w_TGNEA = w_TGNEA, w_TGNEAR = w_TGNEAR,
        alpha_cell = alpha_cell
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_genes_pseudotime_info) = base::names(interest_cell_type_Branches)

  # Surface failures from parallel workers (mclapply stores them as 'try-error' otherwise).
  if (!base::is.null(sparsity_cut)) {
    message("Applying sparsity cutoff post-filtering.")
    for (cell_type in base::names(interest_cell_type_Branches)) {

      # Keep genes whose detection (non-zero) frequency meets the per-branch cutoff.
      for (feat in base::rownames(sparsity_cut[[cell_type]])) {
        slot = base::paste0("Branches_", feat)
        interest_cell_type_genes_pseudotime_info[[cell_type]][[slot]] = base::lapply(
          base::seq_along(interest_cell_type_genes_pseudotime_info[[cell_type]][[slot]]),
          function(n) {
            df = interest_cell_type_genes_pseudotime_info[[cell_type]][[slot]][[n]]
            cutoff = sparsity_cut[[cell_type]][feat, base::paste0("branch", n)]
            df[base::rowMeans(df != 0) >= cutoff, , drop = FALSE]
          }
        )
      }

      # Synchronize the TF/TG name lists with the filtered matrices.
      for (n in base::seq_along(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]])) {
        interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFs"]][[n]] =
          base::rownames(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFE"]][[n]])
        interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGs"]][[n]] =
          base::rownames(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGE"]][[n]])
      }
    }
  }

  ### End.
  base::setwd(original_dir)
  message("Working directory restored to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Obtaining Gene Expression/Activity Vectors along Pseudotime. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_genes_pseudotime_info)
}
