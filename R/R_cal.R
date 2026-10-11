#' @title Calculate the Effect Efficiency R in the scDARI Gene Regulatory Network Model
#'
#' @description Computes the effect efficiency \eqn{R_{i,j}}, a smooth gating function
#' of the prior regulatory strength \eqn{\theta_{i,j}}.
#' It combines a small linear baseline with a double-sigmoid term,
#' and normalizes to exactly 1 at \eqn{\theta = 1},
#' so that strengths below the threshold \eqn{\tau} are strongly attenuated,
#' while strengths above \eqn{\tau} are retained and increase monotonically.
#' \eqn{R_{i,j}} multiplies both gene-expression update equations
#' and only rescales their time scale.
#'
#' @param theta_TF_TG Numeric vector. Prior regulatory strength
#' \eqn{\theta_{i,j}} in \[0, 1\].
#' @param tau Numeric between 0 and 1. Cell-type-specific binarization threshold;
#' it scales the linear baseline and is the center of the sharp gate.
#' @param r1 Positive numeric. Controls the below-threshold linear baseline;
#' larger values give a lower baseline, typically \eqn{r_1 > 100}.
#' Default is 200.
#' @param r2 Positive numeric. Controls the steepness of the rise above the threshold;
#' it is the only \eqn{r} parameter trained per interaction.
#' @param r3 Positive numeric. Position of the first sigmoid,
#' typically \eqn{r_3 \le 1}.
#' Default is 1.
#' @param r4 Positive numeric. Controls the steepness of the gate at \eqn{\tau},
#' typically \eqn{r_4 > 100}.
#' Default is 1000.
#' @param r5 Positive numeric. Amplitude scaling of the sigmoid term.
#' Default is 1.
#'
#' @return Numeric vector of the same length as \code{theta_TF_TG}.
#' Effect efficiency \eqn{R_{i,j}}: near 0 for \eqn{\theta} below \eqn{\tau},
#' rising sharply at \eqn{\tau} and monotonically above it,
#' with \eqn{R_{i,j}(1) = 1}.
#'
#' @details The calculation has two steps:
#' \enumerate{
#'   \item Builds the unnormalized value as the sum of the linear baseline
#'     \eqn{\theta/(r_1\tau)} and a double-sigmoid term;
#'     the first sigmoid, centered at \eqn{r_3}, shapes the rise above the threshold,
#'     and the second, centered at \eqn{\tau} with steepness \eqn{r_4 r_2},
#'     provides the sharp gate.
#'   \item Divides by \eqn{r_6}, the function value at \eqn{\theta = 1},
#'     which enforces \eqn{R_{i,j}(1) = 1}.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' R_val = R_cal(
#'   theta_TF_TG = c(0.2, 0.5, 0.8),
#'   tau = 0.3,
#'   r1 = 200,
#'   r2 = 1,
#'   r3 = 1,
#'   r4 = 1000,
#'   r5 = 1
#' )
#' }
R_cal = function(theta_TF_TG, tau, r1 = 200, r2, r3 = 1, r4 = 1000, r5 = 1) {
  # Unnormalized effect efficiency: linear baseline plus a double-sigmoid term.
  R_unnormalize = theta_TF_TG / (r1 * tau) + r5 / (
    (1 + 10^(r2 * (r3 - theta_TF_TG))) * (1 + 10^(r4 * r2 * (tau - theta_TF_TG)))
  )

  # Normalization constant: value at theta = 1, so R(1) = 1.
  r6 = 1 / (r1 * tau) + r5 / (
    (1 + 10^(r2 * (r3 - 1))) * (1 + 10^(r4 * r2 * (tau - 1)))
  )

  R = R_unnormalize / r6
  return(R)
}
