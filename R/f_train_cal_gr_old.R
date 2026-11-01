#' @title Numerical Gradient of the scDARI Training Objective
#'
#' @description Computes the numerical gradient of the scDARI training objective
#' using forward finite differences.
#' The partial derivatives with respect to the model parameters can be evaluated
#' in parallel, providing the search direction for gradient-based optimization.
#'
#' @param params Numeric vector. Model parameters for one TF-TG pair,
#' with order matching \code{params_names}.
#' @param h_max Tiny positive numeric. Maximum finite-difference step size.
#' Defaults to 0.01.
#' @param h_min_percent Tiny positive numeric. Step size as a fraction of the
#' parameter range, used together with \code{h_max}.
#' Defaults to 0.01.
#' @param ncores Positive integer. Number of cores for parallel evaluation of the
#' partial derivatives.
#' Defaults to 1.
#' @param params_names Character vector. Names corresponding to the elements of
#' \code{params}.
#' @param train_data Data frame. Training data for one pair, containing the prior
#' regulatory strength \code{theta_i} and group expression columns prefixed with
#' \code{TFE_} and \code{TGE_}.
#' @param tau Positive numeric. Regulatory binarization threshold.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths.
#' @param params_upper Numeric vector. Upper parameter bounds.
#' @param params_lower Numeric vector. Lower parameter bounds.
#' @param params_init Numeric vector. Initial parameter values.
#'
#' @return Numeric vector of the same length as \code{params}, the partial
#' derivatives of the training objective approximated by forward differences.
#' NA and infinite components are replaced with 0.
#'
#' @details The step size and derivatives are formed as follows:
#' \enumerate{
#'   \item For each parameter, the step is the minimum of \code{h_max} and
#'     \code{h_min_percent} times the parameter range.
#'   \item The baseline objective is computed once and reused for every partial
#'     derivative.
#'   \item Each perturbed point is kept within the bounds; when the parameter is
#'     already at its upper bound, the partial derivative is set to 0.
#'   \item Partial derivatives are evaluated with \code{\link[parallel]{mclapply}},
#'     and non-finite components are replaced with 0.
#' }
#'
#' @note
#' \itemize{
#'   \item Internally calls \code{\link{f_train_cal_old}}.
#'   \item Forward differences are cheaper but less accurate than central
#'     differences.
#'   \item Parallel evaluation relies on \code{\link[parallel]{mclapply}};
#'     set \code{ncores = 1} on Windows.
#' }
#'
#' @seealso \code{\link{f_train_cal_old}}, \code{\link{gradient_ascent_train_old}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' grad = f_train_cal_gr_old(
#'   params = param_vec,
#'   h_max = 0.01,
#'   h_min_percent = 0.01,
#'   ncores = 1,
#'   params_names = param_names,
#'   train_data = train_data,
#'   tau = tau,
#'   Tslot_K = Tslot_K,
#'   params_upper = param_upper,
#'   params_lower = param_lower,
#'   params_init = param_init
#' )
#' }
f_train_cal_gr_old = function(
    params,
    h_max = 0.01,
    h_min_percent = 0.01,
    ncores = 1,
    params_names = params_names,
    train_data = train_data,
    tau = tau,
    Tslot_K = Tslot_K,
    params_upper = params_upper,
    params_lower = params_lower,
    params_init = params_init
) {
  vec = base::as.numeric(params)
  n = base::length(vec)

  # Baseline objective at the current point.
  f0 = scDARI::f_train_cal_old(
    params = vec, params_names = params_names, train_data = train_data,
    tau = tau, Tslot_K = Tslot_K, params_init = params_init,
    params_upper = params_upper, params_lower = params_lower
  )

  # Per-parameter step: minimum of h_max and h_min_percent of the range.
  h_vec = base::pmin(base::rep(h_max, n),
                     (params_upper - params_lower) * h_min_percent)

  grad = base::unlist(parallel::mclapply(base::seq_len(n), function(i) {
    par_perturb = vec
    # Keep the perturbed point within bounds; at the upper bound, return 0.
    perturbed = base::min(vec[i] + h_vec[i], params_upper[i])
    if (perturbed == vec[i]) return(0)
    par_perturb[i] = perturbed
    gr = (scDARI::f_train_cal_old(
      params = par_perturb, params_names = params_names, train_data = train_data,
      tau = tau, Tslot_K = Tslot_K, params_init = params_init,
      params_upper = params_upper, params_lower = params_lower
    ) - f0) / (perturbed - vec[i])
    return(gr)
  }, mc.cores = ncores, mc.silent = TRUE, mc.allow.recursive = FALSE))

  # Non-finite components carry no usable direction.
  grad[base::is.na(grad) | base::is.infinite(grad)] = 0
  return(grad)
}
