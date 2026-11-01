#' @title Calculate the Training Objective for One TF-TG Pair
#'
#' @description User-facing objective function for a single TF-TG pair, used as
#' the fitness function of the genetic algorithm and as the scoring function
#' inside gradient ascent. It ensures a fixed context is available, then
#' delegates the numerical evaluation to the hot routine
#' \code{\link{f_train_core}}. The objective combines equal-weight prediction
#' error, correction-term regularization and a discrete-stability penalty,
#' negated and scaled for maximization.
#'
#' @param params Numeric vector of length 11. Candidate model parameters for one
#' regulon, in the order given by \code{params_names}. Passed unchanged to
#' \code{\link{f_train_core}}.
#' @param params_names Character vector. Names corresponding to the elements of
#' \code{params}; used only when the context has to be built.
#' @param train_data Data frame or one-row matrix. Training data for one
#' regulon, containing the prior strength \code{theta_i} and the alternating
#' \code{TFE_k}/\code{TGE_k} group columns; used only when the context has to
#' be built.
#' @param tau Numeric strictly between 0 and 1. Regulatory binarization
#' threshold, used to compute the effect efficiency \eqn{R}; used only when the
#' context has to be built.
#' @param Tslot_K Positive numeric vector. Pseudotime interval lengths between
#' consecutive groups; used only when the context has to be built.
#' @param params_init Numeric vector. Initial parameter values in the order of
#' \code{params_names}. Retained for the currently disabled boundary barrier
#' and therefore not used in the returned value.
#' @param params_upper Numeric vector. Upper parameter bounds in the order of
#' \code{params_names}; retained for the disabled boundary barrier.
#' @param params_lower Numeric vector. Lower parameter bounds in the order of
#' \code{params_names}; retained for the disabled boundary barrier.
#' @param ctx Optional list. Fixed context from \code{\link{f_train_context}}.
#' Defaults to \code{NULL}, meaning the context is built once internally from
#' \code{params_names}, \code{train_data}, \code{tau} and \code{Tslot_K}, so
#' the function remains usable on its own. Supplying a prebuilt context avoids
#' rebuilding it on every call and is the fast path used by
#' \code{\link{model_train}}. Closely related to the four context-building
#' arguments above.
#'
#' @return A single finite numeric value, the negated and scaled training
#' objective; higher is better. It equals
#' \eqn{-100\,(\mathrm{MSE} + \beta_1^2 + \beta_2^2 +
#' \text{stability penalty})}, with \eqn{\mathrm{MSE}} the equal-weight mean of
#' the two squared update-equation residuals over lag-1 steps. A non-feasible
#' candidate, with non-finite intermediate values, returns the finite fallback
#' \eqn{-10^{12}}.
#'
#' @details The wrapper performs only two steps:
#' \enumerate{
#'   \item Builds the context with \code{\link{f_train_context}} when
#'     \code{ctx} is \code{NULL}; otherwise it reuses the supplied context.
#'   \item Calls \code{\link{f_train_core}} on \code{params} and returns its
#'     value unchanged.
#' }
#' The detailed construction of the MSE, the stability penalty and the
#' objective is documented in \code{\link{f_train_core}}.
#'
#' @note
#' \itemize{
#'   \item The objective is negated because the optimizers maximize while the
#'     underlying loss is minimized.
#'   \item The boundary barrier is currently disabled, so \code{params_init},
#'     \code{params_upper} and \code{params_lower} do not affect the result;
#'     they are kept for future testing.
#'   \item This function is called many thousands of times during training;
#'     pass a prebuilt \code{ctx} to avoid repeated context construction.
#' }
#'
#' @seealso Internally uses \code{\link{f_train_context}} and
#'   \code{\link{f_train_core}}. Called by \code{\link{model_train}} as the
#'   genetic-algorithm fitness. Related gradient helper:
#'   \code{\link{f_train_cal_gr}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj_value = f_train_cal(
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
f_train_cal = function(
    params,
    params_names,
    train_data,
    tau,
    Tslot_K,
    params_init,
    params_upper,
    params_lower,
    ctx = NULL
) {
  if (base::is.null(ctx)) {
    ctx = scDARI::f_train_context(
      params_names = params_names, train_data = train_data,
      tau = tau, Tslot_K = Tslot_K
    )
  }
  return(scDARI::f_train_core(params = params, ctx = ctx))
}
