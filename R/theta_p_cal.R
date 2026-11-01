#' @title Optimize the Regulatory Strength for One TF-TG Pair
#'
#' @description Holds the trained dynamic parameters fixed and treats the
#' regulatory strength as the only decision variable, using one-dimensional
#' Brent optimization to find the value that minimizes the expression
#' prediction error. It is called after parameter training to turn a fitted
#' dynamical model into a predicted, continuous regulatory strength.
#'
#' @param data Data frame or one-row matrix. Data for one regulon, containing
#' the prior strength \code{theta_i} and the group columns prefixed with
#' \code{TFE_} and \code{TGE_}. Defaults to a \code{train_data} variable in
#' the calling environment. Closely related to \code{Tslot_K}, which gives
#' the interval lengths for the same groups.
#' @param best_params Numeric vector of length 11. Trained model parameters;
#' all dynamic parameters are read from it and held fixed. Closely related to
#' \code{params_names}.
#' @param params_names Character vector. Names corresponding to the elements
#' of \code{best_params}, used to name the vector before each parameter is
#' read.
#' @param tau Numeric strictly between 0 and 1. Regulatory binarization
#' threshold. It also defines the optimization starting point, which is
#' clamped from \code{tau} to lie strictly inside the search interval.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths
#' between consecutive groups; the lag-1 subset multiplies the residuals.
#'
#' @return A single numeric value, the optimal regulatory strength returned
#' by the Brent search and always bounded within \[0.01, 0.99\].
#'
#' @details Only the regulatory strength varies, and the calculation proceeds
#' in three stages:
#' \enumerate{
#'   \item Computes once every theta-independent quantity: the two series,
#'     Hill saturations \eqn{v_i}, \eqn{v_j}, clamped \eqn{v_j}, boundary
#'     factor \eqn{rho}, feedback strength \eqn{S} and interval lengths.
#'   \item Defines the objective as the equal-weight prediction loss used in
#'     \code{\link{loss_cal}}; the effect efficiency \eqn{R(\theta)} from
#'     \code{\link{R_cal}} is the only term that changes with \eqn{theta}.
#'   \item Runs Brent optimization with \code{\link[stats]{optim}} on
#'     \[0.01, 0.99\] and returns the minimizing strength.
#' }
#'
#' @note
#' \itemize{
#'   \item The 0.01-0.99 bounds avoid the extreme strengths at which the
#'     Hill and gating functions can become numerically unstable.
#'   \item The starting point is clamped to \eqn{0.01 + 10^{-6}, 0.99 - 10^{-6}}
#'     so it is never placed exactly on an endpoint.
#'   \item A non-feasible objective evaluation returns the large positive
#'     fallback \eqn{10^{12}}.
#' }
#'
#' @seealso Internally uses \code{\link{R_cal}}, \code{\link{Hill_cal}} and
#'   \code{\link{S_cal}}. Called by \code{\link{model_train}}. Related loss
#'   function: \code{\link{loss_cal}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' predicted_theta = theta_p_cal(
#'   data = train_data,
#'   best_params = best_params,
#'   params_names = params_names,
#'   tau = tau,
#'   Tslot_K = Tslot_K
#' )
#' }
theta_p_cal = function (
    data = train_data,
    best_params,
    params_names,
    tau,
    Tslot_K
) {
  # Large positive fallback for non-feasible evaluations.
  bad_loss = 1e12

  pa_best = base::as.numeric(best_params)
  base::names(pa_best) = params_names

  # Fixed dynamic parameters.
  r2 = base::as.numeric(pa_best["r2"])
  a1 = base::as.numeric(pa_best["alpha1"])
  a2 = base::as.numeric(pa_best["alpha2"])
  b1 = base::as.numeric(pa_best["beta1"])
  b2 = base::as.numeric(pa_best["beta2"])
  s = base::as.numeric(pa_best["s"])
  u1i = base::as.numeric(pa_best["u1i"])
  u1j = base::as.numeric(pa_best["u1j"])
  u2i = base::as.numeric(pa_best["u2i"])
  u2j = base::as.numeric(pa_best["u2j"])
  v_tilde = base::as.numeric(pa_best["v_tilde"])

  TFE = base::as.numeric(data[, base::grep("TFE_", base::colnames(data), value = TRUE)])
  TGE = base::as.numeric(data[, base::grep("TGE_", base::colnames(data), value = TRUE)])
  Length = base::length(TFE)
  K_1 = base::seq_len(Length - 1L)

  vi = scDARI::Hill_cal(x = TFE[K_1], Dissociation_Constant = u1i, Hill_Coefficient = u2i)
  vj = scDARI::Hill_cal(x = TGE[K_1], Dissociation_Constant = u1j, Hill_Coefficient = u2j)
  vij = vi - vj
  vj_safe = base::pmin(base::pmax(vj, 1e-12), 1 - 1e-12)

  # Theta-independent boundary factor and feedback strength.
  rho_vec = TFE[K_1] / (TFE[K_1] + 0.01)
  S_vec = scDARI::S_cal(s = s, v = vj_safe, v_tilde = v_tilde)
  Tslot = Tslot_K[K_1]

  # Keep the starting point strictly inside the optimization interval.
  par0 = base::min(base::max(tau, 0.01 + 1e-6), 0.99 - 1e-6)

  result = stats::optim(par = par0, fn = function(theta) {
    R = scDARI::R_cal(
      theta_TF_TG = theta,
      tau = tau, r1 = 200, r2 = r2, r3 = 1, r4 = 1000, r5 = 1
    )
    if (!base::is.finite(R)) return(bad_loss)

    # Vectorized residuals over every lag-1 step.
    error_f1 = TGE[K_1] + R * a1 * vij * Tslot + b1 - TGE[K_1 + 1L]
    error_f2 = TFE[K_1] + rho_vec * R * a2 * S_vec * Tslot + b2 - TFE[K_1 + 1L]
    MSE = base::mean(error_f1^2 + error_f2^2)

    if (!base::is.finite(MSE)) return(bad_loss)
    return(MSE)
  }, lower = 0.01, upper = 0.99, method = "Brent")

  return(result$par)
}
