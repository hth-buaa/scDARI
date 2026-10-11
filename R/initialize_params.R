#' @title Initialize Constant Parameters for the scDARI Gene Regulatory Network Model
#'
#' @description Builds the constant parameter template of the scDARI gene regulatory
#' network model, shared by every TF-TG interaction.
#' Each trainable scalar receives an initial value with an upper and lower bound
#' for subsequent optimization.
#' No input is required.
#'
#' @return A data frame with one row per trainable parameter,
#' and three columns \code{initialization}, \code{upper}, \code{lower}.
#' The 11 row names are \code{r2}, \code{alpha1}, \code{alpha2}, \code{beta1},
#' \code{beta2}, \code{s}, \code{u1i}, \code{u2i}, \code{u1j}, \code{u2j},
#' and \code{v_tilde}, grouped as follows:
#' \itemize{
#'   \item \code{r2}: growth steepness of the effect efficiency above the threshold.
#'   \item \code{alpha1}, \code{alpha2}: conversion factors of the net-growth
#'     and feedback increments.
#'   \item \code{beta1}, \code{beta2}: correction terms of the two update equations.
#'   \item \code{s}: shape parameter of the feedback effect amplitude.
#'   \item \code{u1i}, \code{uij}: apparent dissociation constants of the TF
#'     and target Hill functions.
#'   \item \code{u2i}, \code{u2j}: Hill coefficients of the TF and target
#'     Hill functions.
#'   \item \code{v_tilde}: environmental demand saturation.
#' }
#'
#' @details The initial values and bounds are:
#' \itemize{
#'   \item \code{r2}: initialization 1, bounds \[0.370506, 1.113383\].
#'   \item \code{alpha1}, \code{alpha2}: initialization 1, bounds \[0.1, 50\].
#'   \item \code{beta1}, \code{beta2}: initialization 0, bounds \[-5, 5\].
#'   \item \code{s}: initialization 1.01, bounds \[1.00001, 4.9\].
#'   \item \code{u1i}, \code{u1i}: initialization 1, bounds \[0.1, 777600\].
#'   \item \code{u2j}, \code{u2j}: initialization 1, bounds \[0.1, 4.9\].
#'   \item \code{v_tilde}: initialization 0.5, bounds \[0.1, 0.9\].
#' }
#'
#' @note
#' \itemize{
#'   \item The function creates and enters \code{4 Build Prediction Model},
#'     restoring the original working directory on exit.
#'   \item The template is constant across interactions; per-interaction fitting
#'     starts from these initial values.
#'   \item The bounds are preset from biological interpretability and numerical
#'     stability, and generally do not need manual adjustment.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' initial_params = initialize_params()
#' base::saveRDS(
#'   initial_params,
#'   file = "./4 Build Prediction Model/initial_params.rds"
#' )
#' }
initialize_params = function() {
  ### Start.
  t_start = base::Sys.time()
  message("Run: Initializing model parameters. Time: ", t_start, ".")

  # Create and enter the working folder; the original directory is restored on exit.
  original_dir = base::getwd()
  new_folder = "4 Build Prediction Model"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder created: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  base::on.exit(base::setwd(original_dir), add = TRUE)
  message("Working directory set to: ", base::getwd(), ".")

  ### Build the constant parameter template shared by every TF-TG interaction.
  # Each row is ordered as c(initialization, upper, lower).
  param_matrix = base::rbind(
    # Effect efficiency parameter.
    r2      = c(1,   1.113383, 0.370506),

    # Effect conversion parameters.
    alpha1  = c(1,   50,       0.1),
    alpha2  = c(1,   50,       0.1),

    # Linear correction parameters.
    beta1   = c(0,   5,        -5),
    beta2   = c(0,   5,        -5),

    # Feedback shape parameter.
    s       = c(1.01, 4.9,     1.00001),

    # Hill apparent dissociation constants (TF, target).
    u1i     = c(1,   777600,   0.1),
    u1j     = c(1,   777600,   0.1),

    # Hill coefficients (TF, target).
    u2i     = c(1,   4.9,      0.1),
    u2j     = c(1,   4.9,      0.1),

    # Environmental demand saturation.
    v_tilde = c(0.5, 0.9,      0.1)
  )
  initial_params = base::data.frame(
    initialization = param_matrix[, 1],
    upper          = param_matrix[, 2],
    lower          = param_matrix[, 3],
    row.names      = base::rownames(param_matrix)
  )

  ### End.
  base::setwd(original_dir)
  message("Working directory restored to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Initializing model parameters. Time: ", t_end, ".")
  message("Running time: ", base::format(t_end - t_start), ".")
  return(initial_params)
}
