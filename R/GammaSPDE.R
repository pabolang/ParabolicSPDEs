#' Calculate Gamma 
#'
#' Approximation of the constant Gamma which results by non-negligible covariances of squared increments in the SPDE model.
#' @param n natrual number for the approximation of the infinite sum. Default is 10000.
#' @keywords Approximation of Gamma.
#' @references Bibinger, M. and Trabs, M. (2017) Volatility estimation for stochastic PDEs using high-frequency observations
#' @export
#' @seealso [ParabolicSPDEs::estimateParametersSPDE], [ParabolicSPDEs::kappa.test].
#' @examples
#' Gamma <- GammaSPDE()

GammaSPDE <- function(n=10000){
  I <- function(r){2*sqrt(r+1)-sqrt(r+2)-sqrt(r)}
  erg <- sum(I(0:n)^2)
  return(1/(pi)*(erg+2))
}
