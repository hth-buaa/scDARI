#' @title Calculate the Pure Prediction Loss for a Single TF-TG Pair
#'
#' @description Computes the equal-weight mean squared prediction error for one
#' TF-TG pair using fitted parameters.
#' It contains only the prediction error of the two update equations, without
#' regularization or stability terms, and is used for model evaluation.
#'
#' @param data Data frame. Data for one TF-TG pair, containing the prior
#' regulatory strength \code{theta_i} and group expression columns prefixed with
#' \code{TFE_} and \code{TGE_}.
#' @param params Numeric vector. Fitted model parameter vector.
#' @param params_names Character vector. Names corresponding to the elements of
#' \code{params}.
#' @param tau Positive numeric. Regulatory binarization threshold.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths.
#'
#' @return Single finite numeric value. Equal-weight mean squared prediction
#' error, averaged over lag-1 steps as
#' \eqn{\mathrm{mean}(error_{f_1}^2 + error_{f_2}^2)};
#' lower values indicate better prediction.
#' Non-feasible evaluations return the large positive fallback \eqn{10^{12}}.
#'
#' @details The pure prediction loss is computed as follows:
#' \enumerate{
#'   \item Extracts the TFE and TGE group series from \code{data}.
#'   \item For each lag-1 step, forms the prediction errors of \eqn{f_1} and
#'     \eqn{f_2}, and returns the mean of their squared errors.
#' }
#'
#' @note
#' \itemize{
#'   \item Differs from \code{\link{f_train_cal_old}}:
#'     this function excludes the correction and stability terms,
#'     and is neither negated nor scaled.
#'   \item \code{data} must contain columns prefixed with \code{TFE_},
#'     \code{TGE_} and a column named \code{theta_i}.
#'   \item It can be applied to training, validation or test datasets.
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
#' train_loss = loss_cal_old(
#'   data = train_data,
#'   params = result_gradient$parameters,
#'   params_names = params_names,
#'   tau = tau,
#'   Tslot_K = Tslot_K
#' )
#' }
loss_cal_old = function (
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

  if (!base::is.finite(MSE)) return(bad_loss)
  return(MSE)
}

