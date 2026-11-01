#' @title Calculate the Training Objective for the scDARI Gene Regulatory Network Model
#'
#' @description Computes the training objective for a single TF-TG pair.
#' It sums the equal-weight mean squared prediction error of the two discrete
#' update equations, a penalty on the correction terms, and a discrete-stability
#' penalty, then negates and scales the result for maximization-based optimizers
#' such as genetic algorithms and gradient ascent.
#' Non-feasible individuals return a finite fallback value.
#'
#' @param params Numeric vector. Candidate model parameters for one TF-TG pair,
#' with order matching \code{params_names}.
#' @param params_names Character vector. Names corresponding to the elements of
#' \code{params}.
#' @param train_data Data frame. Training data for one pair, containing the prior
#' regulatory strength \code{theta_i} and group expression columns prefixed with
#' \code{TFE_} and \code{TGE_}.
#' @param tau Positive numeric. Regulatory binarization threshold, used to compute
#' the effect efficiency \eqn{R}.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths between
#' consecutive groups.
#' @param params_init Numeric vector. Initial parameter values, with order matching
#' \code{params_names}; retained for the currently disabled boundary barrier.
#' @param params_upper Numeric vector. Upper parameter bounds, with order matching
#' \code{params_names}; retained for the currently disabled boundary barrier.
#' @param params_lower Numeric vector. Lower parameter bounds, with order matching
#' \code{params_names}; retained for the currently disabled boundary barrier.
#'
#' @return A single finite numeric value, the negated and scaled training objective;
#' higher values correspond to lower loss and better fit.
#' The objective combines:
#' \itemize{
#'   \item Equal-weight mean squared error of the two update equations,
#'     averaged over lag-1 steps as \eqn{\mathrm{mean}(error_{f_1}^2 + error_{f_2}^2)}.
#'   \item Correction-term penalty \eqn{\beta_1^2 + \beta_2^2}.
#'   \item Discrete-stability penalty, derived from the eigenvalues of the
#'     linearization and the upper-quartile pseudotime step.
#' }
#' Non-feasible individuals, with non-finite intermediate values,
#' return the finite fallback \eqn{-10^{12}}.
#'
#' @details For one TF-TG pair, the function proceeds as follows:
#' \enumerate{
#'   \item Names the parameter vector and extracts the individual parameters, and extracts the TFE/TGE group series.
#'   \item For each lag-1 step, forms the prediction errors of \eqn{f_1} and
#'     \eqn{f_2}, and averages their squared errors into the MSE.
#'   \item Adds the correction-term penalty \eqn{\beta_1^2+\beta_2^2}.
#'   \item Computes the equilibrium level, local gains, Jacobi trace and
#'     determinant, discriminant and eigenvalues, and the stability penalty
#'     against the upper-quartile pseudotime step.
#'   \item Returns \eqn{-100\,(\mathrm{MSE} + \beta_1^2+\beta_2^2 + \text{stability penalty})}.
#' }
#'
#' @note
#' \itemize{
#'   \item The value is negated because genetic algorithms maximize, whereas the underlying loss is minimized.
#'   \item The two residuals are equally weighted; the weights may be adjusted according to measurement precision.
#'   \item The boundary barrier is currently disabled and retained for testing;
#'     \code{params_init}, \code{params_upper} and \code{params_lower} do not affect the returned value while it is disabled.
#'   \item This function is called repeatedly during optimization and is computationally intensive.
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
#' obj_value = f_train_cal_old(
#'   params = param_vec,
#'   params_names = param_names,
#'   train_data = train_data,
#'   tau = tau,
#'   Tslot_K = Tslot_K,
#'   params_init = param_init,
#'   params_upper = param_upper,
#'   params_lower = param_lower
#' )
#' }
f_train_cal_old = function(
    params,
    params_names,
    train_data,
    tau,
    Tslot_K,
    params_init,
    params_upper,
    params_lower
) {
  # Finite fallback for any non-feasible individual, so GA never receives NA/NaN/Inf.
  bad_fitness = -1e12

  base::names(params) = params_names
  base::names(params_init) = params_names
  base::names(params_upper) = params_names
  base::names(params_lower) = params_names

  TFE = base::as.numeric(train_data[, base::grep("TFE_", base::colnames(train_data), value = TRUE)])
  TGE = base::as.numeric(train_data[, base::grep("TGE_", base::colnames(train_data), value = TRUE)])
  Length = base::length(TFE)

  r2 = base::as.numeric(params["r2"])
  r2_i = base::as.numeric(params_init["r2"])
  r2_u = base::as.numeric(params_upper["r2"])
  r2_l = base::as.numeric(params_lower["r2"])
  R = scDARI::R_cal(
    theta_TF_TG = base::as.numeric(train_data[, "theta_i"]),
    tau = tau, r1 = 200, r2 = r2, r3 = 1, r4 = 1000, r5 = 1
  )
  if (!base::is.finite(R)) return(bad_fitness)

  a1 = base::as.numeric(params["alpha1"])
  a2 = base::as.numeric(params["alpha2"])
  a_i = base::as.numeric(params_init["alpha1"])
  a_u = base::as.numeric(params_upper["alpha1"])
  a_l = base::as.numeric(params_lower["alpha1"])

  b1 = base::as.numeric(params["beta1"])
  b2 = base::as.numeric(params["beta2"])

  s = base::as.numeric(params["s"])
  s_i = base::as.numeric(params_init["s"])
  s_u = base::as.numeric(params_upper["s"])
  s_l = base::as.numeric(params_lower["s"])

  u1i = base::as.numeric(params["u1i"])
  u1j = base::as.numeric(params["u1j"])
  u1_i = base::as.numeric(params_init["u1i"])
  u1_u = base::as.numeric(params_upper["u1i"])
  u1_l = base::as.numeric(params_lower["u1i"])

  u2i = base::as.numeric(params["u2i"])
  u2j = base::as.numeric(params["u2j"])
  u2_i = base::as.numeric(params_init["u2i"])
  u2_u = base::as.numeric(params_upper["u2i"])
  u2_l = base::as.numeric(params_lower["u2i"])

  v_tilde = base::as.numeric(params["v_tilde"])
  v_tilde_i = base::as.numeric(params_init["v_tilde"])
  v_tilde_u = base::as.numeric(params_upper["v_tilde"])
  v_tilde_l = base::as.numeric(params_lower["v_tilde"])

  vi = scDARI::Hill_cal(x = TFE[-Length], Dissociation_Constant = u1i, Hill_Coefficient = u2i)
  vj = scDARI::Hill_cal(x = TGE[-Length], Dissociation_Constant = u1j, Hill_Coefficient = u2j)
  vij = vi - vj
  # Clamp saturation away from the endpoints for S, where 1/v or 1/(1-v) diverges.
  vj_safe = base::pmin(base::pmax(vj, 1e-12), 1 - 1e-12)

  MSE = base::mean(base::sapply(2:Length, function(K) {
    K_1 = K - 1
    Tslot_K_K_1 = Tslot_K[K_1]
    TGE_K_1 = TGE[K_1]
    TGE_K = TGE[K]
    vij_K_1 = vij[K_1]
    error_f1 = TGE_K_1 + R * a1 * vij_K_1 * Tslot_K_K_1 + b1 - TGE_K
    TFE_K_1 = TFE[K_1]
    TFE_K = TFE[K]
    rho_K_1 = TFE_K_1 / (TFE_K_1 + 0.01)
    S_K_1 = scDARI::S_cal(s = s, v = vj_safe[K_1], v_tilde = v_tilde)
    error_f2 = TFE_K_1 + rho_K_1 * R * a2 * S_K_1 * Tslot_K_K_1 + b2 - TFE_K
    error_f1 ^ 2 + error_f2 ^ 2
  }))
  if (!base::is.finite(MSE)) return(bad_fitness)

  # Equilibrium-existence penalty on the correction terms.
  b_L2 = b1^2 + b2^2

  # # Boundary barrier penalty function.
  # boundary_penalty <- function(p, p_min, p_max, p0) {
  #   (p0 - p_min)^2 / (p - p_min) + (p_max - p0)^2 / (p_max - p) + p_min - p_max
  # }
  # r2_boundary_penalty = 0.0001 * boundary_penalty(p = r2, p_min = r2_l, p_max = r2_u, p0 = r2_i)
  # a1_boundary_penalty = 0.0001 * boundary_penalty(p = a1, p_min = a_l, p_max = a_u, p0 = a_i)
  # a2_boundary_penalty = 0.0001 * boundary_penalty(p = a2, p_min = a_l, p_max = a_u, p0 = a_i)
  # s_boundary_penalty = 0.0001 * boundary_penalty(p = s, p_min = s_l, p_max = s_u, p0 = s_i)
  # u1i_boundary_penalty = 0.0001 * boundary_penalty(p = u1i, p_min = u1_l, p_max = u1_u, p0 = u1_i)
  # u1j_boundary_penalty = 0.0001 * boundary_penalty(p = u1j, p_min = u1_l, p_max = u1_u, p0 = u1_i)
  # u2i_boundary_penalty = 0.0001 * boundary_penalty(p = u2i, p_min = u2_l, p_max = u2_u, p0 = u2_i)
  # u2j_boundary_penalty = 0.0001 * boundary_penalty(p = u2j, p_min = u2_l, p_max = u2_u, p0 = u2_i)
  # v_tilde_boundary_penalty = 0.0001 * boundary_penalty(p = v_tilde, p_min = v_tilde_l, p_max = v_tilde_u, p0 = v_tilde_i)
  # sum_boundary_penalty = r2_boundary_penalty + a1_boundary_penalty + a2_boundary_penalty + s_boundary_penalty +
  #   u1i_boundary_penalty + u1j_boundary_penalty + u2i_boundary_penalty + u2j_boundary_penalty + v_tilde_boundary_penalty

  # Conversion coefficients of the continuous system.
  c1 = R * a1
  c2 = R * a2

  # Equilibrium TF level and local gains.
  x_ast = (v_tilde * u1i / (1 - v_tilde))^(1 / u2i)
  gain_a = u2i / u1i^(1 / u2i) * v_tilde^(1 - 1 / u2i) * (1 - v_tilde)^(1 + 1 / u2i)
  gain_b = u2j / u1j^(1 / u2j) * v_tilde^(1 - 1 / u2j) * (1 - v_tilde)^(1 + 1 / u2j)
  gain_c = s * (1 / v_tilde^2 + 1 / (1 - v_tilde)^2)
  rho_ast = x_ast / (x_ast + 0.01)
  if (!base::all(base::is.finite(c(x_ast, gain_a, gain_b, gain_c, rho_ast, c1, c2)))) {
    return(bad_fitness)
  }

  # Trace, determinant and discriminant of the linearization.
  tr_J = -c1 * gain_b
  det_J = c1 * c2 * rho_ast * gain_a * gain_b * gain_c
  Delta = tr_J^2 - 4 * det_J

  # Discrete-stability penalty; non-finite discriminant is treated as unstable.
  if (!base::is.finite(Delta)) {
    stability_penalty = 1e6
  } else {
    if (Delta >= 0) {
      # Real eigenvalues.
      sqrt_Delta = base::sqrt(Delta)
      lambda1_real = (tr_J + sqrt_Delta) / 2
      lambda1_imag = 0
      lambda2_real = (tr_J - sqrt_Delta) / 2
      lambda2_imag = 0
    } else {
      # Conjugate complex eigenvalues.
      sqrt_neg_Delta = base::sqrt(-Delta)
      lambda1_real = tr_J / 2
      lambda1_imag = sqrt_neg_Delta / 2
      lambda2_real = tr_J / 2
      lambda2_imag = -sqrt_neg_Delta / 2
    }
    denom1 = lambda1_real^2 + lambda1_imag^2
    denom2 = lambda2_real^2 + lambda2_imag^2
    if (denom1 == 0 || denom2 == 0) {
      stability_penalty = 1e6
    } else {
      h_max_p1 = -lambda1_real / denom1
      h_max_p2 = -lambda2_real / denom2
      h_max_p = 2 * base::min(h_max_p1, h_max_p2)
      # Use the upper-quartile step length so a few extreme long segments
      # do not dominate the stability penalty.
      Tslot_K_threshold = base::as.numeric(stats::quantile(Tslot_K, probs = 0.75))
      stability_penalty = (Tslot_K_threshold - h_max_p) /
        (1 + 10^(-100 * (Tslot_K_threshold - h_max_p)))
    }
  }

  objective = -100 * (MSE + b_L2 + stability_penalty)
  # -100 * (MSE + b_L2 + sum_boundary_penalty + stability_penalty) (kept for later testing).
  if (!base::is.finite(objective)) return(bad_fitness)
  return(objective)
}
