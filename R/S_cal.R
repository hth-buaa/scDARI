#' @title Calculate the Feedback Strength S in the scDARI Gene Regulatory Network Model
#'
#' @description Computes \eqn{S_j}, the derivative of the unimodal feedback effect
#' magnitude \eqn{\widetilde S_j} with respect to the target degradation saturation
#' \eqn{v_j}. \eqn{S_j} is the feedback strength used in the transcription-factor update
#' equation \eqn{f_2}: it is positive below the demand saturation, negative above it,
#' and zero at the demand point.
#'
#' @param s Positive numeric. Shape parameter \eqn{s_j} of the feedback effect magnitude,
#' with \eqn{s_j \ge 1} for endpoint-bounded feedback.
#' @param v Numeric strictly between 0 and 1. Current target degradation saturation
#' \eqn{v_j = U_j(y)}, i.e. the Hill output of the target gene.
#' @param v_tilde Numeric strictly between 0 and 1. Environmental demand saturation
#' \eqn{\widetilde v}, at which the effect magnitude reaches its peak.
#'
#' @return Numeric. Feedback strength \eqn{S_j = \mathrm d\widetilde S_j/\mathrm d v_j}.
#' Positive when \eqn{v < v\_tilde} (the transcription factor is upregulated),
#' negative when \eqn{v > v\_tilde} (the transcription factor is downregulated),
#' and zero when \eqn{v = v\_tilde}.
#'
#' @note An equivalent closed form is:
#'
#'   Numerator = s * (v * vm)^(s - 1) * (v_tilde - v)
#'     * ((1 - v) * (1 - v_tilde) + v_tilde * v)
#'
#'   Denominator = (v_tilde * vm_tilde)^(s + 1)
#'     * exp(s * (2 * v_tilde - 1) * (v_tilde - v) / (v_tilde * vm_tilde))
#'
#'   S = Numerator / Denominator,
#'
#' where vm = 1 - v and vm_tilde = 1 - v_tilde.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' S_val = S_cal(
#'   s = 1.5,
#'   v = U_j_current,
#'   v_tilde = 0.5
#' )
#' }
S_cal = function(s, v, v_tilde) {
  # Complements of the current and demand saturations.
  vm = 1 - v
  vm_tilde = 1 - v_tilde

  # Effect magnitude S~ = h1(v; v_tilde)^s.
  numerator = v * vm
  denominator = v_tilde * vm_tilde
  exponent = exp(2 - v / v_tilde - vm / vm_tilde)
  S_tilde = (numerator / denominator * exponent)^s

  # Logarithmic derivative g = d ln h1 / dv.
  g = 1 / v - 1 / vm - 1 / v_tilde + 1 / vm_tilde

  # Feedback strength S = d S~ / dv = s * S~ * g.
  S = s * S_tilde * g
  return(S)
}
