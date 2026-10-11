#' @title Compute Cell-Type-Specific GRN Binarization Thresholds
#'
#' @description Computes a cell-type-specific binarization threshold (\eqn{\tau}) for each initial
#' gene regulatory network. The threshold is derived from the minimum non-zero interaction strength,
#' with a fixed lower bound to avoid over-sensitivity.
#'
#' @param interest_cell_type_iGRN Named list. Output of
#' \code{\link{calculate_iGRN}} or
#' \code{\link{filter_iGRN_by_ChIP_Atlas_priors}}, with each element a target gene (rows) by TF
#' (columns) association matrix or data frame.
#'
#' @return A named numeric vector of thresholds (\eqn{\tau}) with cell types as names.
#' A cell type whose network contains no non-zero edge (e.g. a 0 x 0 matrix) receives \code{tau = 1},
#' so that no edge can pass binarization; a message is printed for such cases.
#'
#' @details
#' The threshold is calculated independently for each cell type as
#' \deqn{\tau = \max(\min(\theta > 0) - 0.005,\ 0.005)}
#' where \eqn{\theta} denotes all non-zero, non-missing interaction strengths in the network.
#'
#' The computation proceeds as follows:
#' \enumerate{
#'   \item Collects all values of the matrix and keeps the non-zero, non-missing ones.
#'   \item Subtracts a small offset of 0.005 from the minimum non-zero value.
#'   \item Enforces a lower bound of 0.005 to prevent an excessively small threshold.
#' }
#'
#' The resulting threshold is typically used to convert continuous interaction strengths into binary
#' edge calls (presence / absence) for downstream discrete modeling.
#'
#' @note
#' \itemize{
#'   \item Thresholds are cell-type-specific, adapting to the strength distribution of each network.
#'   \item The 0.005 lower bound prevents an overly sensitive threshold when the minimum non-zero
#'     value is extremely small.
#'   \item The function creates (if missing) and enters \code{3 get iGRN/}, restoring the original
#'     working directory on exit.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
#' interest_cell_type_tau = get_cell_type_tau(interest_cell_type_iGRN)
#' base::saveRDS(interest_cell_type_tau, file = "./3 get iGRN/interest_cell_type_tau.rds")
#' }
get_cell_type_tau = function(interest_cell_type_iGRN) {
  ### Validate inputs.
  if (!base::is.list(interest_cell_type_iGRN) || base::length(interest_cell_type_iGRN) == 0) {
    base::stop("'interest_cell_type_iGRN' must be a non-empty named list of iGRN matrices.")
  }
  if (base::is.null(base::names(interest_cell_type_iGRN))) {
    base::stop("'interest_cell_type_iGRN' must be named by cell type.")
  }
  for (cell_type in base::names(interest_cell_type_iGRN)) {
    if (!base::is.matrix(interest_cell_type_iGRN[[cell_type]]) &&
        !base::is.data.frame(interest_cell_type_iGRN[[cell_type]])) {
      base::stop(base::sprintf("Element '%s' must be a matrix or data frame.", cell_type))
    }
  }

  ### Start.
  t_start = base::Sys.time()
  message("Run: Computing the binarization threshold tau for each cell type. Time: ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Compute tau for each cell type.
  tau_offset = 0.005 # Subtracted from the minimum non-zero strength.
  tau_floor  = 0.005 # Enforced lower bound.
  tau_empty  = 1     # Threshold for networks with no non-zero edge (no edge can pass).

  interest_cell_type_tau = base::numeric(base::length(interest_cell_type_iGRN))
  base::names(interest_cell_type_tau) = base::names(interest_cell_type_iGRN)

  for (cell_type in base::names(interest_cell_type_iGRN)) {
    # unlist works identically for matrix and data.frame inputs.
    values = base::unlist(interest_cell_type_iGRN[[cell_type]], use.names = FALSE)
    positive_values = values[!base::is.na(values) & values > 0]
    if (base::length(positive_values) == 0) {
      interest_cell_type_tau[cell_type] = tau_empty
      message("[", cell_type, "] No non-zero edge; tau set to ", tau_empty, ".")
    } else {
      interest_cell_type_tau[cell_type] = base::max(
        base::min(positive_values) - tau_offset,
        tau_floor
      )
    }
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Computing the binarization threshold tau for each cell type. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start))
  return(interest_cell_type_tau)
}
