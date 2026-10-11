#' @title Calculate the Hill Function in the scDARI Gene Regulatory Network Model
#'
#' @description Computes the Hill function
#' \eqn{H_{K,n}(x) = x^n/(K + x^n)},
#' a saturating sigmoid used to describe the generation and degradation rates.
#' It returns the fractional saturation corresponding to the molecule level \eqn{x}.
#'
#' @param x Non-negative numeric vector. Concentration or expression level
#' of the regulatory molecule.
#' @param Dissociation_Constant Positive numeric. Apparent dissociation constant
#' \eqn{K} (\eqn{u_1}): the level at which half saturation occurs;
#' smaller values indicate a lower threshold.
#' @param Hill_Coefficient Positive numeric. Hill coefficient \eqn{n} (\eqn{u_2}):
#' values greater than 1 indicate positive cooperativity,
#' equal to 1 indicate no cooperativity,
#' and between 0 and 1 indicate negative cooperativity.
#'
#' @return Numeric vector in \[0, 1), the fractional saturation:
#' 0 at \eqn{x = 0}, increasing strictly, and approaching 1 as \eqn{x} grows.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' y = Hill_cal(
#'   x = expression_level,
#'   Dissociation_Constant = u1i,
#'   Hill_Coefficient = u2i
#' )
#' }
Hill_cal = function(x, Dissociation_Constant, Hill_Coefficient) {
  x = base::as.numeric(x)
  Dissociation_Constant = base::as.numeric(Dissociation_Constant)
  Hill_Coefficient = base::as.numeric(Hill_Coefficient)

  # Fractional saturation H_{K,n}(x) = x^n / (K + x^n).
  x_power = x^Hill_Coefficient
  y = x_power / (Dissociation_Constant + x_power)
  return(y)
}
