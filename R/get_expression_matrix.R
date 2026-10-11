#' @title Extract the scRNA-seq Expression Matrix for GRN Construction
#'
#' @description Combine transcription factors (TFs) and target genes (TGs) into a unified gene set,
#' subset the scRNA-seq Seurat object to those genes,
#' and return the subsetted expression object together with the harmonized gene lists for downstream Gene Regulatory Network (GRN) construction.
#'
#' @param TFs Character vector of transcription factor gene symbols,
#' typically the output of \code{\link{get_TFs_from_priors}}.
#' @param TGs Character vector of target gene symbols,
#' typically the output of \code{\link{filter_TGs}}
#' @param scRNAseq Seurat object of single-cell RNA-seq data.
#' Genes absent from this object are dropped from the subsetted matrix.
#'
#' @returns
#' A list with 4 elements:
#' \itemize{
#'   \item \code{scRNAseq_for_GRN}: Seurat object of scRNA-seq data, subsetted to the union of TGs and TFs
#'   \item \code{gene_for_GRN}: Sorted character vector of the de-duplicated union of TGs and TFs
#'   \item \code{TGs}: Sorted, de-duplicated target gene symbols
#'   \item \code{TFs}: Sorted, de-duplicated transcription factor symbols
#' }
#'
#' @details
#' The function executes the following steps:
#' \enumerate{
#'   \item Creates the dedicated output directory \code{1.3 Data Preprocessing - Filter Expression Matrix} and temporarily switches the working directory (the original path is restored after completion)
#'   \item Combines TGs and TFs into a single sorted, de-duplicated gene set (\code{gene_for_GRN})
#'   \item Subsets the scRNA-seq Seurat object to \code{gene_for_GRN}
#'   \item Returns a list containing the subsetted expression object and the harmonized gene lists
#' }
#'
#' @note
#' \itemize{
#'   \item The returned \code{gene_for_GRN}, \code{TGs}, and \code{TFs} are derived from the inputs only;
#'     they are not intersected with the features present in \code{scRNAseq},
#'     whereas the subsetted \code{scRNAseq_for_GRN} retains only features available in the object.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' TFs = readRDS("./1.1 Data Preprocessing - Get TFs From Database/TFs.rds")
#' TGs = readRDS("./1.2 Data Preprocessing - Filter TGs/TFs.rds")
#' scRNAseq = base::readRDS("./scRNAseq.rds")
#'
#' expression_matrix_info = get_expression_matrix(
#'   TFs = TFs,
#'   TGs = TGs,
#'   scRNAseq = scRNAseq
#' )
#' base::saveRDS(
#'   expression_matrix_info,
#'   file = "./1.3 Data Preprocessing - Filter Expression Matrix/expression_matrix_info.rds"
#' )
#' }
get_expression_matrix = function(
    TFs,
    TGs,
    scRNAseq
) {
  ### Start.
  t_start = base::Sys.time()
  message("Run: Filtering Counts Matrices for Gene Expression. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "1.3 Data Preprocessing - Filter Expression Matrix"
  if (!base::dir.exists(new_folder)) {base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Extract the gene expression matrix.
  # Combine TGs and TFs into one sorted and de-duplicated gene set.
  gene_for_GRN = base::sort(base::unique(base::union(TGs, TFs)))

  # Subset the scRNA-seq object to the genes used for GRN construction.
  message("Extracting gene expression matrix.")
  scRNAseq_for_GRN = base::subset(scRNAseq, features = gene_for_GRN)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Filtering Counts Matrices for Gene Expression. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(base::list(
    "scRNAseq_for_GRN" = scRNAseq_for_GRN,
    "gene_for_GRN" = gene_for_GRN,
    "TGs" = base::sort(base::unique(TGs)),
    "TFs" = base::sort(base::unique(TFs))
  ))
}
