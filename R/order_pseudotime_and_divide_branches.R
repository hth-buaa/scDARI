#' @title Pseudotime Ordering and Interactive Branch Partitioning for scRNA-seq Data
#'
#' @description Builds developmental trajectories and assigns pseudotime for each cell type with Monocle2 (DDRTree).
#' It supports both automatic and manual DDRTree parameter selection,
#' exports trajectory visualizations and the CellDataSet (CDS),
#' then divides cells into branches through interactive input based on the generated plots.
#'
#' @param interest_cell_type_data List. Output of \code{\link{get_interest_cell_type_data}},
#' containing the filtered scRNA-seq object (\code{interest_scRNAseq_for_GRN}) and gene set (\code{interest_gene_for_GRN}) for each cell type.
#' @param auto_param_selection Logical. Whether to let DDRTree select its parameters automatically.
#' Defaults to \code{TRUE}.
#' When \code{FALSE}, the manual parameters below are forwarded to \code{\link[monocle]{reduceDimension}}.
#' @param maxIter Integer. Maximum DDRTree iterations; used only when \code{auto_param_selection = FALSE}.
#' Defaults to \code{20}. See \code{maxIter} in \code{\link[DDRTree]{DDRTree}}.
#' @param sigma Numeric. DDRTree bandwidth; used only when \code{auto_param_selection = FALSE}.
#' Defaults to \code{0.001}. See \code{sigma} in \code{\link[DDRTree]{DDRTree}}.
#' @param lambda Numeric or \code{NULL}. DDRTree regularization; used only when \code{auto_param_selection = FALSE}.
#' Defaults to \code{NULL}; in manual mode an explicit value may be required by DDRTree.
#' See \code{lambda} in \code{\link[DDRTree]{DDRTree}}.
#' @param ncenter Integer or \code{NULL}. Number of DDRTree centroids; used only when \code{auto_param_selection = FALSE}.
#' Defaults to \code{NULL}; in manual mode an explicit value (e.g. \code{3}) is recommended.
#' See \code{ncenter} in \code{\link[DDRTree]{DDRTree}}.
#' @param param.gamma Numeric. DDRTree gamma parameter; used only when \code{auto_param_selection = FALSE}.
#' Defaults to \code{10}. See \code{param.gamma} in \code{\link[DDRTree]{DDRTree}}.
#' @param tol Numeric. DDRTree convergence tolerance; used only when \code{auto_param_selection = FALSE}.
#' Defaults to \code{0.001}. See \code{tol} in \code{\link[DDRTree]{DDRTree}}.
#'
#' @return
#' A nested list named by cell type.
#' Each element is a list of branch data frames, in the order the branches were entered.
#' Every branch data frame has cells as rows (sorted by ascending pseudotime) and two columns:
#' \itemize{
#'   \item \code{Pseudotime}: pseudotime value assigned by Monocle.
#'   \item \code{State}: cell state (segment) assigned by the trajectory.
#' }
#'
#' @details
#' For each cell type, the function executes the following workflow:
#' \enumerate{
#'   \item Converts the scRNA-seq counts into a Monocle \code{CellDataSet} (negative-binomial family).
#'   \item Estimates size factors and dispersions for expression normalization.
#'   \item Sets the ordering genes to \code{interest_gene_for_GRN}.
#'   \item Reduces dimensions and builds the trajectory with DDRTree (automatic or manual parameters).
#'   \item Orders cells along the trajectory to assign pseudotime and states.
#'   \item Saves the CDS and trajectory plots (colored by State and Pseudotime) in the cell-type subdirectory.
#'   \item Interactively reads the number of branches and the ordered state IDs of each branch.
#'   \item Extracts the cells of every requested branch and sorts them by ascending pseudotime.
#' }
#'
#' @section Interactive Input:
#' The function pauses once per cell type. In an interactive R session:
#' \itemize{
#'   \item Enter the number of branches (a positive integer) based on the trajectory plots.
#'   \item For each branch, enter the space-separated integer state IDs that compose it (e.g. \code{1 2 3});
#'     states may be shared across branches (e.g. a common root state).
#' }
#' Inspect the images in the corresponding cell-type subdirectory before answering.
#' Entering an empty branch or an unknown state ID stops execution with the list of valid states.
#'
#' @note
#' \itemize{
#'   \item This is an interactive function relying on \code{\link[base]{readline}};
#'     it cannot be run non-interactively (e.g. in a batch \code{Rscript} job) and therefore processes cell types sequentially.
#'   \item Depends on \code{monocle}, \code{DDRTree} and \code{VGAM}.
#'     \code{igraph} version 2.0.3 and \code{dplyr} version 1.1.4 is recommended, as later versions may be incompatible with Monocle.
#'   \item Creates the directory \code{2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment} with one subdirectory per cell type;
#'     the CDS and plots are saved in the subdirectory and the working directory is restored on exit.
#'   \item If the default trajectory topology is unsatisfactory or the cell count is small,
#'     set \code{auto_param_selection = FALSE} and tune the DDRTree parameters (notably \code{ncenter}).
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform pseudotime analysis and branch partitioning
#' library(DDRTree)
#' packageVersion("igraph")
#' remove.packages("igraph")
#' install_version("igraph", version = "2.0.3")
#'
#' # Automatic DDRTree parameter selection (interactive prompts follow).
#' interest_cell_type_data = base::readRDS(paste0(
#'   "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'   "/interest_cell_type_data.rds"
#' ))
#' interest_cell_type_Branches = order_pseudotime_and_divide_branches(interest_cell_type_data)
#' # ...
#' # Please enter the number of branches based on the pseudotime 2D visualization:
#' #   6
#' # ...
#' # Please enter the integer vector of cell state order for branch 1 (space-separated):
#' #   1 11
#' # Please enter the integer vector of cell state order for branch 2 (space-separated):
#' #   1 2 10
#' # Please enter the integer vector of cell state order for branch 3 (space-separated):
#' #   1 2 3 9
#' # Please enter the integer vector of cell state order for branch 4 (space-separated):
#' #   1 2 3 4 8
#' # Please enter the integer vector of cell state order for branch 5 (space-separated):
#' #   1 2 3 4 5 6
#' # Please enter the integer vector of cell state order for branch 6 (space-separated):
#' #   1 2 3 4 5 7
#' # ...
#'
#' # Manual DDRTree parameter selection (interactive prompts follow), e.g. ncenter = 4:
#' interest_cell_type_Branches = order_pseudotime_and_divide_branches(
#'    interest_cell_type_data = interest_cell_type_data,
#'    auto_param_selection = FALSE,
#'    ncenter = 4
#' )
#' # ...
#' # Please enter the number of branches based on the pseudotime 2D visualization:
#' #    1
#' # Please enter the integer vector of cell state order for branch 1 (space-separated):
#' #    1
#' # ...
#'
#' base::saveRDS(
#'   interest_cell_type_Branches,
#'   file = paste0(
#'     "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'     "/interest_cell_type_Branches.rds"
#'   )
#' )
#' }
order_pseudotime_and_divide_branches = function(
    interest_cell_type_data,
    auto_param_selection = TRUE,
    maxIter = 20,
    sigma = 0.001,
    lambda = NULL,
    ncenter = NULL,
    param.gamma = 10,
    tol = 0.001
) {
  ### Validate parameters.
  if (!base::is.list(interest_cell_type_data) || base::is.null(base::names(interest_cell_type_data))) {
    base::stop("interest_cell_type_data must be a named list from get_interest_cell_type_data().")
  }
  if (!base::is.logical(auto_param_selection) || base::length(auto_param_selection) != 1 ||
      base::is.na(auto_param_selection)) {
    base::stop("auto_param_selection must be TRUE or FALSE.")
  }
  if (!auto_param_selection) {
    if (!(base::is.numeric(maxIter) && base::length(maxIter) == 1 && maxIter >= 1 &&
          maxIter == base::round(maxIter))) {
      base::stop("maxIter must be a positive integer when auto_param_selection = FALSE.")
    }
    if (!(base::is.numeric(sigma) && base::length(sigma) == 1 && sigma > 0)) {
      base::stop("sigma must be a single positive number when auto_param_selection = FALSE.")
    }
    if (!(base::is.numeric(tol) && base::length(tol) == 1 && tol > 0)) {
      base::stop("tol must be a single positive number when auto_param_selection = FALSE.")
    }
    if (!(base::is.numeric(param.gamma) && base::length(param.gamma) == 1 && param.gamma > 0)) {
      base::stop("param.gamma must be a single positive number when auto_param_selection = FALSE.")
    }
    if (!(base::is.null(lambda) || (base::is.numeric(lambda) && base::length(lambda) == 1))) {
      base::stop("lambda must be NULL or a single number when auto_param_selection = FALSE.")
    }
    if (!(base::is.null(ncenter) ||
          (base::is.numeric(ncenter) && base::length(ncenter) == 1 && ncenter >= 1 &&
           ncenter == base::round(ncenter)))) {
      base::stop("ncenter must be NULL or a positive integer when auto_param_selection = FALSE.")
    }
  }

  ### Start.
  t_start = base::Sys.time()
  message("Run: Performing pseudotime ordering and branch partitioning for cell types of interest. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Process one cell type.
  order_pseudotime_and_divide_branches_process = function (interest_scRNAseq_for_GRN, interest_gene_for_GRN) {
    # Build a Monocle CellDataSet from the sparse raw counts.
    message("Building the Monocle CellDataSet from counts.")
    joined = SeuratObject::JoinLayers(interest_scRNAseq_for_GRN)
    counts_mat = SeuratObject::GetAssayData(joined, layer = "counts")
    interest_cds = monocle::newCellDataSet(
      cellData = counts_mat,
      phenoData = methods::new("AnnotatedDataFrame", data = joined@meta.data),
      featureData = methods::new(
        "AnnotatedDataFrame",
        data = base::data.frame(gene_short_name = base::rownames(counts_mat),
                                row.names = base::rownames(counts_mat))
      ),
      lowerDetectionLimit = 0.1,
      expressionFamily = VGAM::negbinomial.size()
    )

    # Estimate size factors, dispersions, and set the ordering genes.
    message("Estimating size factors and dispersions.")
    interest_cds = BiocGenerics::estimateSizeFactors(interest_cds)
    message("It is recommended to use dplyr R package version 1.1.4, as later versions are not compatible with monocle.")
    interest_cds = BiocGenerics::estimateDispersions(interest_cds)

    message("Setting the ordering genes for trajectory inference.")
    interest_cds = monocle::setOrderingFilter(interest_cds, interest_gene_for_GRN)

    # Reduce dimensions with DDRTree (automatic or manual parameters).
    # library(DDRTree)
    message("Reducing dimensions with DDRTree.")
    if (auto_param_selection) {
      interest_cds = monocle::reduceDimension(
        interest_cds, reduction_method = "DDRTree", max_components = 2
      )
    } else {
      interest_cds = monocle::reduceDimension(
        cds = interest_cds, reduction_method = "DDRTree", max_components = 2,
        auto_param_selection = auto_param_selection,
        maxIter = maxIter, sigma = sigma, lambda = lambda,
        ncenter = ncenter, param.gamma = param.gamma, tol = tol
      )
    }

    # Order cells along the trajectory.
    message("Ordering cells along the trajectory.")
    message("It is recommended to use igraph R package version 2.0.3, as later versions are not compatible with monocle.")
    interest_cds = monocle::orderCells(interest_cds)

    # Save the CDS in the current cell-type subdirectory.
    message("Saving the CellDataSet as interest_cds.rds.")
    base::saveRDS(interest_cds, file = "interest_cds.rds")

    # Generating 2D visualizations by pseudotime and cell state saved in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment.
    message("Generating 2D visualizations by pseudotime and cell states saved in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment.")
    n_states = base::length(base::unique(Biobase::pData(interest_cds)$State))
    ggplot2::ggsave(
      "interest_AllCellState_in_trajectory.png", width = 5, height = 5,
      plot = monocle::plot_cell_trajectory(interest_cds, color_by = "State")
    )
    ggplot2::ggsave(
      "interest_EveryCellState_in_trajectory.png",
      width = 3 * base::ceiling(n_states / base::ceiling(base::sqrt(n_states))),
      height = 3 * base::ceiling(base::sqrt(n_states)),
      plot = monocle::plot_cell_trajectory(interest_cds, color_by = "State") +
        ggplot2::facet_wrap(~State, nrow = base::ceiling(base::sqrt(n_states)), axes = "all")
    )
    ggplot2::ggsave(
      "interest_Pseudotime_in_trajectory.png", width = 5, height = 5,
      plot = monocle::plot_cell_trajectory(interest_cds, color_by = "Pseudotime")
    )

    # Interactively read the number of branches.
    available_states = base::sort(base::unique(Biobase::pData(interest_cds)$State))
    message("Valid states in this cell type: ", base::paste(available_states, collapse = ", "), ".")
    message("Interactively input the number of branches (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).")
    n_branch = base::as.integer(base::readline(
      prompt = "Please enter the number of branches based on the trajectory plots (a positive integer, e.g., 3): "
    ))
    if (base::is.na(n_branch) || n_branch < 1) {
      stop("Invalid input: the number of branches must be a positive integer.")
    }

    # Interactively read the ordered state IDs of each branch.
    message("Interactively input the integer vector of cell state order for each branch (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).")
    branches = base::vector("list", n_branch)
    for (i in base::seq_len(n_branch)) {
      states_i = base::as.integer(base::strsplit(base::readline(
        prompt = base::sprintf("Please enter the space-separated state IDs for branch %d (e.g. 1 2 3): ", i)
      ), " ")[[1]])
      if (base::length(states_i) == 0 || base::any(base::is.na(states_i))) {
        base::stop("Invalid input: branch ", i, " requires a non-empty integer state vector.")
      }
      unknown_states = base::setdiff(states_i, available_states)
      if (base::length(unknown_states) > 0) {
        base::stop("Unknown state ID(s) ", base::paste(unknown_states, collapse = ", "),
                   " for branch ", i, ". Valid states are: ",
                   base::paste(available_states, collapse = ", "), ".")
      }
      branches[[i]] = states_i
    }

    # Extract each branch (cells of the requested states) and sort by ascending pseudotime.
    message("Extracting and ordering cells for each branch.")
    Branches = base::vector("list", n_branch)
    pdata = Biobase::pData(interest_cds)
    for (i in base::seq_along(branches)) {
      Branches[[i]] = pdata[pdata$State %in% branches[[i]], c("Pseudotime", "State")]
      Branches[[i]] = Branches[[i]][base::order(Branches[[i]][, 1]), ]
    }

    return(Branches)
  }

  ### Run each cell type inside its own subdirectory.
  interest_cell_type_Branches = stats::setNames(
    base::vector("list", base::length(interest_cell_type_data)),
    base::names(interest_cell_type_data)
  )
  for (name in base::names(interest_cell_type_Branches)) {
    # Start.
    message("Pseudotime ordering and branch partitioning for cell type: ", name, ".")
    original_dir1 = base::getwd()
    new_folder1 = name
    if (!base::dir.exists(new_folder1)) {
      base::dir.create(new_folder1, recursive = TRUE)
      message("Folder created: ", new_folder1, ".")
    } else {message("Folder already exists: ", new_folder, ".")}
    base::setwd(new_folder1)
    message("The current working directory has been switched to: ", base::getwd(), ".")

    # Do.
    interest_cell_type_Branches[[name]] = order_pseudotime_and_divide_branches_process(
      interest_scRNAseq_for_GRN = interest_cell_type_data[[name]]$interest_scRNAseq_for_GRN,
      interest_gene_for_GRN = interest_cell_type_data[[name]]$interest_gene_for_GRN
    )

    # End.
    base::setwd(original_dir1)
    message("The current working directory has been switched to: ", base::getwd(), ".")
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Performing pseudotime ordering and branch partitioning for cell types of interest. Time: ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_Branches)
}
