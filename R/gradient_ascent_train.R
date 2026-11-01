#' @title Gradient Ascent for Local Parameter Refinement
#'
#' @description Maximizes the scDARI training objective with gradient ascent,
#' used to refine a candidate already found by a global search. Each iteration
#' computes a numerical gradient, selects the learning rate with a Brent line
#' search, clips the candidate to the parameter box and accepts it only when
#' the objective increases; the first non-improving iteration ends the search.
#'
#' @param iterations Positive integer. Hard upper bound on the number of
#' gradient-ascent iterations; training may stop earlier on a non-improving
#' step.
#' @param alpha_lower Tiny positive numeric. Lower end of the learning-rate
#' search interval. Defaults to 1e-5. Paired with \code{alpha_upper} and
#' enclosing \code{alpha_guess}.
#' @param alpha_upper Tiny positive numeric. Upper end of the learning-rate
#' search interval. Defaults to 1e-3. Paired with \code{alpha_lower}.
#' @param alpha_guess Numeric within the search interval. Starting point of
#' the Brent search and the value returned if no better rate exists. Defaults
#' to 1e-4. Closely related to the two interval bounds.
#' @param init_theta Numeric vector of length 11. Starting parameters,
#' typically the best solution of the genetic algorithm; the initial objective
#' is evaluated at this point.
#' @param h_max Tiny positive numeric. Upper cap on the finite-difference
#' step, forwarded to \code{\link{f_train_cal_gr}}. Defaults to 0.01.
#' @param h_min_percent Tiny positive numeric. Finite-difference step as a
#' fraction of each parameter's range, forwarded to
#' \code{\link{f_train_cal_gr}}. Defaults to 0.01. Closely related to
#' \code{h_max} and the bound vectors.
#' @param ncores Positive integer. Cores used for the parallel gradient,
#' forwarded to \code{\link{f_train_cal_gr}}. Defaults to 1; use 1 on Windows.
#' @param params_lower Numeric vector. Lower parameter bounds in the order of
#' \code{params_names}; used both for the gradient ranges and for clipping
#' candidate points. Paired with \code{params_upper}.
#' @param params_upper Numeric vector. Upper parameter bounds in the order of
#' \code{params_names}; paired with \code{params_lower}.
#' @param monitor Logical. Whether to print the iteration, objective, learning
#' rate and core count after each step. Defaults to \code{FALSE}, i.e. silent.
#' @param params_names Character vector. Names corresponding to the parameter
#' elements; used when the context has to be built.
#' @param train_data Data frame or one-row matrix. Training data for one
#' regulon; used when the context has to be built.
#' @param tau Numeric strictly between 0 and 1. Regulatory binarization
#' threshold; used when the context has to be built.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths; used
#' when the context has to be built.
#' @param params_init Numeric vector. Initial parameter values; retained for
#' the currently disabled boundary barrier and not used in the result.
#' @param ctx Optional list. Fixed context from \code{\link{f_train_context}}.
#' Defaults to \code{NULL}, so the context is built once internally and the
#' function remains usable standalone; a prebuilt context is the fast path.
#' Closely related to the context-building arguments.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{parameters}: parameter vector achieving the reported
#'     objective.
#'   \item \code{objective_value}: maximum objective value achieved.
#' }
#' The two elements are always consistent: on a non-improving iteration the
#' previous point and its objective are returned, and after exhausting
#' \code{iterations} the last accepted point is returned.
#'
#' @details Every iteration follows the same sequence:
#' \enumerate{
#'   \item Computes the numerical gradient with \code{\link{f_train_cal_gr}}
#'     and marks zero-gradient coordinates as fixed.
#'   \item Runs a Brent line search with \code{\link[stats]{optim}}, minimizing
#'     the negated objective over the learning-rate interval.
#'   \item Builds the candidate at the selected rate and clips it to the
#'     parameter box.
#'   \item Accepts the candidate only when it strictly improves the objective;
#'     otherwise the previous point is kept and the function returns.
#' }
#'
#' @note
#' \itemize{
#'   \item Gradients use forward finite differences, cheaper but less accurate
#'     than central differences.
#'   \item This is a local refinement step; starting from a random point rather
#'     than a GA solution is not the intended use.
#'   \item Parallel gradients use \code{\link[parallel]{mclapply}}; keep
#'     \code{ncores = 1} on Windows.
#' }
#'
#' @seealso Internally uses \code{\link{f_train_context}},
#'   \code{\link{f_train_core}} and \code{\link{f_train_cal_gr}}. Called by
#'   \code{\link{model_train}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' grad_result = gradient_ascent_train(
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
gradient_ascent_train = function(
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
    params_init,
    ctx = NULL
) {
  if (base::is.null(ctx)) {
    ctx = scDARI::f_train_context(
      params_names = params_names, train_data = train_data,
      tau = tau, Tslot_K = Tslot_K
    )
  }

  theta = base::as.numeric(init_theta)
  f0 = scDARI::f_train_core(params = theta, ctx = ctx)

  for (i in base::seq_len(iterations)) {
    # Numerical gradient; zero-gradient coordinates are not updated.
    grad = scDARI::f_train_cal_gr(
      params = theta, h_max = h_max, h_min_percent = h_min_percent,
      ncores = ncores, params_names = params_names, train_data = train_data,
      tau = tau, Tslot_K = Tslot_K, params_upper = params_upper,
      params_lower = params_lower, params_init = params_init, ctx = ctx
    )
    grad_need_cal = grad != 0

    # Brent line search on the learning rate, minimizing the negated objective.
    best_alpha = stats::optim(par = alpha_guess, fn = function(al) {
      theta_new = theta
      theta_new[grad_need_cal] = base::pmin(
        base::pmax(theta[grad_need_cal] + al * grad[grad_need_cal],
                   params_lower[grad_need_cal]),
        params_upper[grad_need_cal])
      -scDARI::f_train_core(params = theta_new, ctx = ctx)
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
