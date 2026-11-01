#' @title Optimize the Regulatory Strength for a Single TF-TG Pair
#'
#' @description Fixes the trained dynamic parameters and treats the regulatory
#' strength as the only variable, using Brent one-dimensional optimization to
#' find the value that minimizes the expression prediction error.
#'
#' @param data Data frame. Data for one TF-TG pair, containing the prior
#' regulatory strength \code{theta_i} and group expression columns prefixed with
#' \code{TFE_} and \code{TGE_}.
#' @param best_params Numeric vector. Trained model parameters.
#' @param params_names Character vector. Names corresponding to the elements of
#' \code{best_params}.
#' @param tau Positive numeric. Regulatory binarization threshold.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths.
#'
#' @return Single numeric value. Optimal regulatory strength, bounded within
#' \[0.01, 0.99\].
#'
#' @details All dynamic parameters are held fixed and only the regulatory
#' strength is optimized:
#' \enumerate{
#'   \item The theta-independent quantities, including the series, Hill
#'     saturations, boundary factor and feedback strength, are computed once.
#'   \item The objective is the equal-weight prediction loss used in
#'     \code{\link{loss_cal_old}}, with the effect efficiency \eqn{R(\theta)} as the
#'     only theta-dependent term.
#'   \item Brent optimization via \code{\link[stats]{optim}} returns the
#'     regulatory strength minimizing prediction error on \[0.01, 0.99\].
#' }
#'
#' @note
#' \itemize{
#'   \item Used to infer the regulatory strength after the dynamic parameters are trained.
#'   \item The 0.01-0.99 bounds avoid extreme values that may cause numerical instability.
#'   \item Non-feasible evaluations use the large positive fallback \eqn{10^{12}}.
#'   \item Computes the effect efficiency \eqn{R} via \code{\link{R_cal}}.
#'   \item Computes the Hill saturations \eqn{v_i}, \eqn{v_j} via \code{\link{Hill_cal}} and the net rate \eqn{v_i-v_j},
#'     and clamps \eqn{v_j} away from the endpoints for \eqn{S}.
#'   \item Computes the derivative of the unimodal feedback effect magnitude \eqn{S} via \code{\link{S_cal}}.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' predicted_theta = theta_p_cal_old(
#'   data = train_data,
#'   best_params = best_params,
#'   params_names = params_names,
#'   tau = tau,
#'   Tslot_K = Tslot_K
#' )
#' }
theta_p_cal_old = function (
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

  vi = scDARI::Hill_cal(x = TFE[-Length], Dissociation_Constant = u1i, Hill_Coefficient = u2i)
  vj = scDARI::Hill_cal(x = TGE[-Length], Dissociation_Constant = u1j, Hill_Coefficient = u2j)
  vij = vi - vj
  # Clamp saturation away from the endpoints for S, where 1/v or 1/(1-v) diverges.
  vj_safe = base::pmin(base::pmax(vj, 1e-12), 1 - 1e-12)

  # Theta-independent boundary factor and feedback strength.
  rho_vec = TFE[-Length] / (TFE[-Length] + 0.01)
  S_vec = base::sapply(vj_safe, function(v) {
    scDARI::S_cal(s = s, v = v, v_tilde = v_tilde)
  })

  # Keep the starting point strictly inside the optimization interval.
  par0 = base::min(base::max(tau, 0.01 + 1e-6), 0.99 - 1e-6)

  result = stats::optim(par = par0, fn = function(theta) {
    R = scDARI::R_cal(
      theta_TF_TG = theta,
      tau = tau, r1 = 200, r2 = r2, r3 = 1, r4 = 1000, r5 = 1
    )
    if (!base::is.finite(R)) return(bad_loss)

    MSE = base::mean(base::sapply(2:Length, function(K) {
      K_1 = K - 1
      Tslot_K_K_1 = Tslot_K[K_1]
      error_f1 = TGE[K_1] + R * a1 * vij[K_1] * Tslot_K_K_1 + b1 - TGE[K]
      error_f2 = TFE[K_1] + rho_vec[K_1] * R * a2 * S_vec[K_1] * Tslot_K_K_1 + b2 - TFE[K]
      error_f1 ^ 2 + error_f2 ^ 2
    }))
    if (!base::is.finite(MSE)) return(bad_loss)
    return(MSE)
  }, lower = 0.01, upper = 0.99, method = "Brent")

  return(result$par)
}
