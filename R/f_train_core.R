#' @title Evaluate the Training Objective from a Prebuilt Context
#'
#' @description Hot inner routine that scores one candidate parameter vector for
#' a single TF-TG pair. The candidate vector is the only quantity that changes
#' between calls; all fixed data are read from a context built once by
#' \code{\link{f_train_context}}. The routine combines prediction error,
#' correction-term regularization and a discrete-stability penalty into one
#' negated, scaled objective for maximization-based optimizers.
#'
#' @param params Numeric vector of length 11. Candidate model parameters for one
#' regulon; its elements are read through the positional indices in
#' \code{ctx$pos}, so its ordering need not be canonical. Closely related to
#' \code{ctx} and to the \code{params_names} used when that context was built.
#' @param ctx List. Fixed context returned by \code{\link{f_train_context}},
#' required and with no default because this routine is the performance-critical
#' inner function and never builds constants itself. It supplies the lag-1
#' expression series, interval lengths, the 75th-percentile step, prior
#' strength, threshold and positional indices. Closely related to \code{params}.
#'
#' @return A single finite numeric value. On a feasible evaluation it is
#' \eqn{-100\,(\mathrm{MSE} + \beta_1^2 + \beta_2^2 +
#' \text{stability penalty})}, where \eqn{\mathrm{MSE}} is the equal-weight
#' mean of the two squared update-equation residuals over lag-1 steps; higher
#' values mean lower loss and a better fit. If any intermediate value is
#' non-finite, the finite fallback \eqn{-10^{12}} is returned.
#'
#' @details Given a candidate and its fixed context, the evaluation proceeds in
#' order:
#' \enumerate{
#'   \item Unpacks the 11 parameters through \code{ctx$pos} and computes the
#'     effect efficiency \eqn{R} from the prior strength with
#'     \code{\link{R_cal}}.
#'   \item Computes the TF and target Hill saturations \eqn{v_i}, \eqn{v_j}
#'     with \code{\link{Hill_cal}}, their difference \eqn{v_i-v_j}, and clamps
#'     \eqn{v_j} away from 0 and 1.
#'   \item Forms the vectorized residuals of the two update equations over
#'     every lag-1 step and averages their squared errors into the MSE.
#'   \item Computes the equilibrium level, local gains, Jacobian trace and
#'     determinant, discriminant and eigenvalues, then the stability penalty
#'     against the context's 75th-percentile pseudotime step.
#'   \item Returns the negated, scaled objective, or the \eqn{-10^{12}}
#'     fallback when any required quantity is non-finite.
#' }
#'
#' @note
#' \itemize{
#'   \item The objective is negated because genetic algorithms and gradient
#'     ascent maximize, whereas the underlying loss is minimized.
#'   \item The two residuals are equally weighted; the \eqn{v_j} clamp and the
#'     \code{rho = TFE/(TFE + 0.01)} boundary factor keep the arithmetic inside
#'     \[0, 1\].
#'   \item A non-finite discriminant or a zero eigenvalue denominator raises the
#'     stability penalty to \eqn{10^6}.
#'   \item Always pass a context from \code{\link{f_train_context}}; calling
#'     this function without one is not supported.
#' }
#'
#' @seealso Internally uses \code{\link{R_cal}}, \code{\link{Hill_cal}} and
#'   \code{\link{S_cal}}. Called by \code{\link{f_train_cal}},
#'   \code{\link{f_train_cal_gr}} and \code{\link{gradient_ascent_train}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ctx = f_train_context(params_names, train_data, tau, Tslot_K)
#' objective = f_train_core(params = candidate_params, ctx = ctx)
#' }
f_train_core = function(params, ctx) {
  bad_fitness = -1e12

  p = base::as.numeric(params)
  pos = ctx$pos
  r2 = p[pos[1L]]
  a1 = p[pos[2L]]; a2 = p[pos[3L]]
  b1 = p[pos[4L]]; b2 = p[pos[5L]]
  s = p[pos[6L]]
  u1i = p[pos[7L]]; u1j = p[pos[8L]]
  u2i = p[pos[9L]]; u2j = p[pos[10L]]
  v_tilde = p[pos[11L]]

  TFE_cur = ctx$TFE_cur; TFE_next = ctx$TFE_next
  TGE_cur = ctx$TGE_cur; TGE_next = ctx$TGE_next

  R = scDARI::R_cal(
    theta_TF_TG = ctx$theta,
    tau = ctx$tau, r1 = 200, r2 = r2, r3 = 1, r4 = 1000, r5 = 1
  )
  if (!base::is.finite(R)) return(bad_fitness)

  vi = scDARI::Hill_cal(x = TFE_cur, Dissociation_Constant = u1i, Hill_Coefficient = u2i)
  vj = scDARI::Hill_cal(x = TGE_cur, Dissociation_Constant = u1j, Hill_Coefficient = u2j)
  vij = vi - vj
  vj_safe = base::pmin(base::pmax(vj, 1e-12), 1 - 1e-12)

  # Vectorized residuals over every lag-1 step.
  Tslot = ctx$Tslot
  error_f1 = TGE_cur + R * a1 * vij * Tslot + b1 - TGE_next
  rho = TFE_cur / (TFE_cur + 0.01)
  S_vec = scDARI::S_cal(s = s, v = vj_safe, v_tilde = v_tilde)
  error_f2 = TFE_cur + rho * R * a2 * S_vec * Tslot + b2 - TFE_next
  MSE = base::mean(error_f1^2 + error_f2^2)
  if (!base::is.finite(MSE)) return(bad_fitness)

  b_L2 = b1^2 + b2^2

  c1 = R * a1
  c2 = R * a2

  x_ast = (v_tilde * u1i / (1 - v_tilde))^(1 / u2i)
  gain_a = u2i / u1i^(1 / u2i) * v_tilde^(1 - 1 / u2i) * (1 - v_tilde)^(1 + 1 / u2i)
  gain_b = u2j / u1j^(1 / u2j) * v_tilde^(1 - 1 / u2j) * (1 - v_tilde)^(1 + 1 / u2j)
  gain_c = s * (1 / v_tilde^2 + 1 / (1 - v_tilde)^2)
  rho_ast = x_ast / (x_ast + 0.01)
  if (!base::all(base::is.finite(c(x_ast, gain_a, gain_b, gain_c, rho_ast, c1, c2)))) {
    return(bad_fitness)
  }

  tr_J = -c1 * gain_b
  det_J = c1 * c2 * rho_ast * gain_a * gain_b * gain_c
  Delta = tr_J^2 - 4 * det_J

  if (!base::is.finite(Delta)) {
    stability_penalty = 1e6
  } else {
    if (Delta >= 0) {
      sqrt_Delta = base::sqrt(Delta)
      lambda1_real = (tr_J + sqrt_Delta) / 2
      lambda1_imag = 0
      lambda2_real = (tr_J - sqrt_Delta) / 2
      lambda2_imag = 0
    } else {
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
      Tslot_K_threshold = ctx$q75
      stability_penalty = (Tslot_K_threshold - h_max_p) /
        (1 + 10^(-100 * (Tslot_K_threshold - h_max_p)))
    }
  }

  objective = -100 * (MSE + b_L2 + stability_penalty)
  if (!base::is.finite(objective)) return(bad_fitness)
  return(objective)
}
