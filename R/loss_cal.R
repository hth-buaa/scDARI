#' @title Calculate the Pure Prediction Loss for One TF-TG Pair
#'
#' @description Computes the equal-weight mean squared prediction error of the
#' two update equations for one regulon using fitted parameters. It contains
#' only prediction error, with no correction-term or stability penalty and no
#' negation, so it is the metric used to monitor fit and to choose the best
#' epoch.
#'
#' @param data Data frame or one-row matrix. Data for one regulon, containing
#' the prior strength \code{theta_i} and the group columns prefixed with
#' \code{TFE_} and \code{TGE_}. Defaults to a \code{train_data} variable in
#' the calling environment, so explicit passing is preferred in saved scripts.
#' Closely related to \code{Tslot_K}, which gives the interval lengths for
#' the same groups.
#' @param params Numeric vector of length 11. Fitted model parameters.
#' Defaults to \code{result_gradient$parameters} from the calling
#' environment; elements are read by name after \code{params_names} is
#' attached. Closely related to \code{params_names}.
#' @param params_names Character vector. Names corresponding to the elements
#' of \code{params}, used to name the vector before each parameter is read.
#' @param tau Numeric strictly between 0 and 1. Regulatory binarization
#' threshold, used to compute the effect efficiency \eqn{R}.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths
#' between consecutive groups; the lag-1 subset multiplies the residuals.
#'
#' @return A single finite numeric value, the equal-weight mean squared
#' prediction error \eqn{\mathrm{mean}(error_{f_1}^2 + error_{f_2}^2)} over
#' lag-1 steps; lower is better. A non-feasible evaluation, with a non-finite
#' \eqn{R} or MSE, returns the large positive fallback \eqn{10^{12}}.
#'
#' @details The loss is computed in one pass:
#' \enumerate{
#'   \item Locates the \code{TFE_}/\code{TGE_} columns by name and extracts
#'     the two group series.
#'   \item Computes \eqn{R}, the Hill saturations \eqn{v_i}, \eqn{v_j}, the
#'     clamped \eqn{v_j} and the feedback strength \eqn{S}.
#'   \item Forms the vectorized residuals of the two update equations over
#'     every lag-1 step and returns the mean of their squared errors.
#' }
#'
#' @note
#' \itemize{
#'   \item Differs from \code{\link{f_train_cal}}: this function excludes the
#'     correction and stability terms and is neither negated nor scaled.
#'   \item \code{data} must contain a \code{theta_i} column and
#'     \code{TFE_}/\code{TGE_} columns.
#'   \item The same function works on training, validation or test data.
#' }
#'
#' @seealso Internally uses \code{\link{R_cal}}, \code{\link{Hill_cal}} and
#'   \code{\link{S_cal}}. Called by \code{\link{model_train}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' train_loss = loss_cal(
#'   data = train_data,
#'   params = result_gradient$parameters,
#'   params_names = params_names,
#'   tau = tau,
#'   Tslot_K = Tslot_K
#' )
#' }
loss_cal = function (
    data = train_data,
    params = result_gradient$parameters,
    params_names,
    tau,
    Tslot_K
) {
  # Large positive fallback for non-feasible evaluations.
  bad_loss = 1e12
  base::names(params) = params_names

  TFE = base::as.numeric(data[, base::grep("TFE_", base::colnames(data), value = TRUE)])
  TGE = base::as.numeric(data[, base::grep("TGE_", base::colnames(data), value = TRUE)])
  Length = base::length(TFE)
  K_1 = base::seq_len(Length - 1L)

  r2 = base::as.numeric(params["r2"])
  R = scDARI::R_cal(
    theta_TF_TG = base::as.numeric(data[, "theta_i"]),
    tau = tau, r1 = 200, r2 = r2, r3 = 1, r4 = 1000, r5 = 1
  )
  if (!base::is.finite(R)) return(bad_loss)

  a1 = base::as.numeric(params["alpha1"])
  a2 = base::as.numeric(params["alpha2"])
  b1 = base::as.numeric(params["beta1"])
  b2 = base::as.numeric(params["beta2"])
  s = base::as.numeric(params["s"])
  u1i = base::as.numeric(params["u1i"])
  u1j = base::as.numeric(params["u1j"])
  u2i = base::as.numeric(params["u2i"])
  u2j = base::as.numeric(params["u2j"])
  v_tilde = base::as.numeric(params["v_tilde"])

  vi = scDARI::Hill_cal(x = TFE[K_1], Dissociation_Constant = u1i, Hill_Coefficient = u2i)
  vj = scDARI::Hill_cal(x = TGE[K_1], Dissociation_Constant = u1j, Hill_Coefficient = u2j)
  vij = vi - vj
  vj_safe = base::pmin(base::pmax(vj, 1e-12), 1 - 1e-12)

  # Vectorized residuals over every lag-1 step.
  Tslot = Tslot_K[K_1]
  error_f1 = TGE[K_1] + R * a1 * vij * Tslot + b1 - TGE[K_1 + 1L]
  rho = TFE[K_1] / (TFE[K_1] + 0.01)
  S_vec = scDARI::S_cal(s = s, v = vj_safe, v_tilde = v_tilde)
  error_f2 = TFE[K_1] + rho * R * a2 * S_vec * Tslot + b2 - TFE[K_1 + 1L]

  MSE = base::mean(error_f1^2 + error_f2^2)
  if (!base::is.finite(MSE)) return(bad_loss)
  return(MSE)
}
