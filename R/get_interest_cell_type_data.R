#' @title Prepare Cell-Type-Specific scRNA-seq Data for GRN Construction
#'
#' @description For each cell type of interest, subset the GRN-ready scRNA-seq object to that cell type,
#' remove all-zero genes and cells, optionally filter genes and cells by their zero-expression rate,
#' and harmonize the transcription factor (TF) and target gene (TG) sets.
#' Cell types are processed independently (optionally in parallel),
#' and the results are returned as a named list ready for downstream Gene Regulatory Network (GRN) inference.
#'
#' @param Basic_Info List returned by \code{\link{get_expression_matrix}},
#' containing \code{scRNAseq_for_GRN}, \code{gene_for_GRN}, \code{TFs}, and \code{TGs}.
#' @param interest_cell_type Character vector of cell type labels to process.
#' Labels must exist in the \code{seurat_annotations} metadata column of \code{Basic_Info$scRNAseq_for_GRN};
#' duplicates are removed automatically.
#' @param alpha_gene Numeric in \[0, 1\].
#' Maximum allowed fraction of cells in which a gene has zero expression;
#' genes with a zero-expression rate above this value are dropped.
#' Defaults to \code{1} (rate filtering disabled; only all-zero genes are removed).
#' @param alpha_cell Numeric in \[0, 1\].
#' Maximum allowed fraction of retained genes in which a cell has zero expression;
#' cells with a zero-expression rate above this value are dropped.
#' Defaults to \code{1} (rate filtering disabled; only all-zero cells are removed).
#' @param ncores Single positive integer, number of parallel cores passed to \code{mc.cores} in \code{\link[parallel]{mclapply}}.
#' Defaults to \code{1} (sequential).
#' Values greater than 1 rely on forking and are effective on Linux/macOS;
#' on Windows the call silently falls back to sequential execution.
#'
#' @returns A named list with one element per cell type in \code{interest_cell_type}
#' (names are the cell type labels). Each element is a list of 4:
#' \itemize{
#'   \item \code{interest_scRNAseq_for_GRN}: Seurat object subsetted to the cell type and the retained genes and cells
#'   \item \code{interest_gene_for_GRN}: Sorted, de-duplicated union of retained TFs and TGs
#'   \item \code{interest_TFs}: Sorted TF symbols retained for the cell type
#'   \item \code{interest_TGs}: Sorted TG symbols retained for the cell type
#' }
#'
#' @details
#' The function executes the following steps:
#' \enumerate{
#'   \item Validates inputs: every requested cell type must be present in \code{seurat_annotations},
#'     \code{alpha_gene} and \code{alpha_cell} must lie in \[0, 1\], and cell type labels are de-duplicated;
#'     then creates the stage directory and temporarily switches the working directory (restored on exit).
#'   \item For each cell type:
#'   \itemize{
#'     \item Subsets the Seurat object to that cell type and removes all-zero genes and all-zero cells.
#'     \item Keeps genes whose zero-cell fraction is at most \code{alpha_gene},
#'       then keeps cells whose zero-gene fraction is at most \code{alpha_cell}.
#'     \item Intersects the retained genes with \code{gene_for_GRN} and rebuilds the TF/TG gene sets.
#'   }
#'   \item Processes all cell types with \code{\link[parallel]{mclapply}} and names the result by cell type;
#'     if any worker fails, execution stops and the offending cell type(s) are reported.
#'   \item Returns the nested named list (no files are written by this function).
#' }
#'
#' @note
#' \itemize{
#'   \item With the defaults \code{alpha_gene = alpha_cell = 1}, rate filtering is disabled
#'     and only genes/cells that are zero across the board are removed.
#'   \item Zero-expression rates are computed directly on the sparse \code{counts} matrix
#'     without coercing it to a dense matrix.
#'   \item The Seurat object must contain a \code{counts} layer and a \code{seurat_annotations} metadata column.
#'   \item The function returns an R list and does not save files itself;
#'     persist the result with \code{\link[base]{saveRDS}} if needed.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' expression_matrix_info = base::readRDS(
#'   "./1.3 Data Preprocessing - Filter Expression Matrix/expression_matrix_info.rds"
#' )
#'
#' interest_cell_type_data = get_interest_cell_type_data(
#'   Basic_Info = expression_matrix_info,
#'   interest_cell_type = c("CD14+ Mono", "Naive CD4 T"),
#'   alpha_gene = 0.95,
#'   alpha_cell = 0.95,
#'   ncores = 4
#' )
#'
#' base::saveRDS(
#'   interest_cell_type_data,
#'   file = paste0("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'                 "/interest_cell_type_data.rds")
#' )
#' }
get_interest_cell_type_data = function(
    Basic_Info,
    interest_cell_type,
    alpha_gene = 1,
    alpha_cell = 1,
    ncores = 1
) {
  ### Validate inputs before any processing.
  if (base::is.null(Basic_Info$scRNAseq_for_GRN)) {
    base::stop("Basic_Info must contain 'scRNAseq_for_GRN' (see get_expression_matrix()).")
  }
  interest_cell_type = base::unique(interest_cell_type)
  missing_cell_types = base::setdiff(
    interest_cell_type,
    Basic_Info$scRNAseq_for_GRN$seurat_annotations
  )
  if (base::length(missing_cell_types) > 0) {
    base::stop("Cell type(s) not found in seurat_annotations: ",
               base::paste(missing_cell_types, collapse = ", "), ".")
  }
  if (alpha_gene < 0 || alpha_gene > 1 || alpha_cell < 0 || alpha_cell > 1) {
    base::stop("alpha_gene and alpha_cell must be between 0 and 1.")
  }
  if (!(base::is.numeric(ncores) && base::length(ncores) == 1 && ncores >= 1)) {
    base::stop("ncores must be a single positive number.")
  }

  ### Start.
  t_start = base::Sys.time()
  message("Run: Extracting Data for Cell Types of Interest. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Process one cell type.
  get_interest_cell_type_data_process = function(cell_type,
                                                 alpha_gene = 1,
                                                 alpha_cell = 1,
                                                 scRNAseq_for_GRN,
                                                 gene_for_GRN,
                                                 TFs,
                                                 TGs) {
    message("The cell label is ", cell_type, ".")

    # Subset to the cell type, then drop all-zero genes and all-zero cells.
    message("Subsetting cells and removing all-zero genes and cells.")
    interest_scRNAseq_for_GRN = base::subset(scRNAseq_for_GRN, subset = seurat_annotations == cell_type)
    counts_mat = SeuratObject::GetAssayData(interest_scRNAseq_for_GRN, layer = "counts")
    interest_scRNAseq_for_GRN = interest_scRNAseq_for_GRN[Matrix::rowSums(counts_mat) > 0, ]
    counts_mat = SeuratObject::GetAssayData(interest_scRNAseq_for_GRN, layer = "counts")
    interest_scRNAseq_for_GRN = interest_scRNAseq_for_GRN[, Matrix::colSums(counts_mat) > 0]
    if (base::nrow(interest_scRNAseq_for_GRN) == 0 || base::ncol(interest_scRNAseq_for_GRN) == 0) {
      base::stop("No gene or cell left for cell type '", cell_type, "' after removing all-zero entries.")
    }

    # Keep genes whose zero-cell fraction is at most alpha_gene.
    message("Keeping genes with zero-expression rate <= ", alpha_gene, ".")
    counts_mat = SeuratObject::GetAssayData(interest_scRNAseq_for_GRN, layer = "counts")
    gene_zero_frequency = 1 - Matrix::rowSums(counts_mat != 0) / base::ncol(counts_mat)
    keep_genes = base::names(gene_zero_frequency[gene_zero_frequency <= alpha_gene])
    if (base::length(keep_genes) == 0) {
      base::stop("No gene left for cell type '", cell_type,
                 "' after gene filtering; relax alpha_gene.")
    }
    interest_scRNAseq_for_GRN = base::subset(interest_scRNAseq_for_GRN, features = keep_genes)

    # Keep cells whose zero-gene fraction is at most alpha_cell.
    message("Keeping cells with zero-expression rate <= ", alpha_cell, ".")
    counts_mat = SeuratObject::GetAssayData(interest_scRNAseq_for_GRN, layer = "counts")
    cell_zero_frequency = 1 - Matrix::colSums(counts_mat != 0) / base::nrow(counts_mat)
    keep_cells = base::names(cell_zero_frequency[cell_zero_frequency <= alpha_cell])
    if (base::length(keep_cells) == 0) {
      base::stop("No cell left for cell type '", cell_type,
                 "' after cell filtering; relax alpha_cell.")
    }
    interest_scRNAseq_for_GRN = base::subset(interest_scRNAseq_for_GRN, cells = keep_cells)

    # Harmonize the retained genes with gene_for_GRN and rebuild the TF/TG sets.
    message("Finalizing TF/TG gene sets.")
    retained_genes = base::intersect(gene_for_GRN, base::rownames(interest_scRNAseq_for_GRN))
    interest_TFs = base::sort(base::intersect(TFs, retained_genes))
    interest_TGs = base::sort(base::intersect(TGs, retained_genes))
    interest_gene_for_GRN = base::sort(base::union(interest_TFs, interest_TGs))
    interest_scRNAseq_for_GRN =
      interest_scRNAseq_for_GRN[base::sort(base::intersect(base::rownames(interest_scRNAseq_for_GRN), interest_gene_for_GRN)), ]
    message("Cell type '", cell_type, "': retained ",
            base::ncol(interest_scRNAseq_for_GRN), " cells and ",
            base::nrow(interest_scRNAseq_for_GRN), " genes.")

    # End.
    return(base::list("interest_scRNAseq_for_GRN" = interest_scRNAseq_for_GRN,
                      "interest_gene_for_GRN" = interest_gene_for_GRN,
                      "interest_TFs" = interest_TFs,
                      "interest_TGs" = interest_TGs))
  }

  ### Iterate over cell types (optionally in parallel).
  interest_cell_type_data = parallel::mclapply(
    X = interest_cell_type,
    FUN = function(cell_type) {
      get_interest_cell_type_data_process(
        alpha_gene = alpha_gene,
        alpha_cell = alpha_cell,
        cell_type = cell_type,
        scRNAseq_for_GRN = Basic_Info[["scRNAseq_for_GRN"]],
        gene_for_GRN = Basic_Info[["gene_for_GRN"]],
        TFs = Basic_Info[["TFs"]],
        TGs = Basic_Info[["TGs"]]
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_data) = interest_cell_type

  # Surface failures from parallel workers (mclapply stores them as 'try-error' otherwise).
  failed = interest_cell_type[base::vapply(
    interest_cell_type_data,
    function(x) base::inherits(x, "try-error"),
    logical(1)
  )]
  if (base::length(failed) > 0) {
    base::stop("Processing failed for cell type(s): ", base::paste(failed, collapse = ", "), ".")
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd())
  t_end = base::Sys.time()
  message("Finish: Extracting Data for Cell Types of Interest. Time: ", t_end)
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_data)
}
