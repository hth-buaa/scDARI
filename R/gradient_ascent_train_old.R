#' @title Gradient Ascent for Local Parameter Refinement
#'
#' @description Performs gradient ascent to maximize the scDARI training objective.
#' It combines numerical gradients, a Brent line search for the learning rate,
#' bound clipping and early stopping to refine parameters after a global search.
#'
#' @param iterations Positive integer. Maximum number of gradient ascent
#' iterations.
#' @param alpha_lower Tiny positive numeric. Lower bound of the learning-rate
#' search range.
#' Defaults to 1e-5.
#' @param alpha_upper Tiny positive numeric. Upper bound of the learning-rate
#' search range.
#' Defaults to 1e-3.
#' @param alpha_guess Numeric within the search range. Initial guess for the
#' learning rate.
#' Defaults to 1e-4.
#' @param init_theta Numeric vector. Starting parameters, typically the solution
#' of a global search such as the genetic algorithm.
#' @param h_max Tiny positive numeric. Maximum finite-difference step size.
#' Defaults to 0.01.
#' @param h_min_percent Tiny positive numeric. Finite-difference step as a
#' fraction of the parameter range.
#' Defaults to 0.01.
#' @param ncores Positive integer. Number of cores for parallel gradient
#' evaluation via \code{\link[parallel]{mclapply}}.
#' Defaults to 1; set to 1 on Windows.
#' @param params_lower Numeric vector. Lower parameter bounds.
#' @param params_upper Numeric vector. Upper parameter bounds.
#' @param monitor Logical. Whether to report iteration progress with messages.
#' Defaults to FALSE.
#' @param params_names Character vector. Names corresponding to the parameter
#' elements.
#' @param train_data Data frame. Training data for one TF-TG pair.
#' @param tau Positive numeric. Regulatory binarization threshold.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths.
#' @param params_init Numeric vector. Initial parameter values.
#'
#' @return List with two elements:
#' \itemize{
#'   \item \code{parameters}: parameter vector achieving the reported objective.
#'   \item \code{objective_value}: maximum objective value achieved.
#' }
#' The search stops when an iteration does not improve the objective,
#' and the returned parameters remain consistent with the reported objective.
#'
#' @details Each iteration follows the workflow below:
#' \enumerate{
#'   \item Computes the numerical gradient, and skips zero-gradient coordinates.
#'   \item Selects the learning rate with a Brent line search via
#'     \code{\link[stats]{optim}}.
#'   \item Builds the candidate point and clips it to the parameter bounds.
#'   \item Commits the candidate only when the objective increases;
#'     otherwise the previous point is kept and the search stops.
#' }
#'
#' @note
#' \itemize{
#'   \item Gradients use forward finite differences, which are cheaper but less
#'     accurate than central differences.
#'   \item This function is intended for local refinement after a global search,
#'     not for training from a random initialization.
#'   \item Parallel gradients rely on \code{\link[parallel]{mclapply}};
#'     set \code{ncores = 1} on Windows.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' grad_result = gradient_ascent_train_old(
#'   iterations = 30,
#'   init_theta = as.numeric(ga_result@solution[1, ]),
#'   params_lower = param_lower,
#'   params_upper = param_upper,
#'   params_names = param_names,
#'   train_data = train_data,
#'   tau = tau,
#'   Tslot_K = Tslot_K,
#'   params_init = param_init,
#'   ncores = 1
#' )
#' }
gradient_ascent_train_old = function(
    iterations,
    alpha_lower = 1e-05,
    alpha_upper = 1e-03,
    alpha_guess = 1e-04,
    init_theta,
    h_max = 0.01,
    h_min_percent = 0.01,
    ncores = 1,
    params_lower,
    params_upper,
    monitor = FALSE,
    params_names,
    train_data,
    tau,
    Tslot_K,
    params_init
) {
  theta = base::as.numeric(init_theta)
  f0 = scDARI::f_train_cal_old(
    params = theta, params_names = params_names, train_data = train_data,
    tau = tau, Tslot_K = Tslot_K, params_init = params_init,
    params_upper = params_upper, params_lower = params_lower
  )

  for (i in base::seq_len(iterations)) {
    # Numerical gradient; zero-gradient coordinates are not updated.
    grad = scDARI::f_train_cal_gr_old(
      params = theta, h_max = h_max, h_min_percent = h_min_percent,
      ncores = ncores, params_names = params_names, train_data = train_data,
      tau = tau, Tslot_K = Tslot_K, params_upper = params_upper,
      params_lower = params_lower, params_init = params_init
    )
    grad_need_cal = grad != 0

    # Brent line search on the learning rate, minimizing the negated objective.
    best_alpha = stats::optim(par = alpha_guess, fn = function(al) {
      theta_new = theta
      theta_new[grad_need_cal] = base::pmin(
        base::pmax(theta[grad_need_cal] + al * grad[grad_need_cal],
                   params_lower[grad_need_cal]),
        params_upper[grad_need_cal])
      -scDARI::f_train_cal_old(
        params = theta_new, params_names = params_names, train_data = train_data,
        tau = tau, Tslot_K = Tslot_K, params_init = params_init,
        params_upper = params_upper, params_lower = params_lower)
    }, lower = alpha_lower, upper = alpha_upper, method = "Brent")

    # Candidate point at the selected learning rate.
    theta_new = theta
    theta_new[grad_need_cal] = base::pmin(
      base::pmax(theta[grad_need_cal] + best_alpha$par * grad[grad_need_cal],
                 params_lower[grad_need_cal]),
      params_upper[grad_need_cal])
    current_f = -best_alpha$value

    if (f0 >= current_f) {
      # No improvement: keep the previous point consistent with f0.
      if (monitor) {
        message("Grad | iter = ", i, " | Best = ", base::round(f0, 6),
                " | learning_rate = ", base::round(best_alpha$par, 6),
                " | ncores = ", ncores)
      }
      return(base::list(parameters = theta, objective_value = f0))
    }

    # Commit the improving step.
    if (monitor) {
      message("Grad | iter = ", i, " | Best = ", base::round(current_f, 6),
              " | learning_rate = ", base::round(best_alpha$par, 6),
              " | ncores = ", ncores)
    }
    theta = theta_new
    f0 = current_f
  }

  base::list(parameters = theta, objective_value = f0)
}

