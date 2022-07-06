#' Simulates multiple parabolic stochastic partial differential equations in one space dimension using Monte-Carlo repetitions
#'
#' Simulate Monte-Carlo samples of a SPDE model in one space dimension on a discrete grid with \code{N} temporal and \code{M} spatial points, where the grid points are equidistant within the unit square. Each repetition uses the "replacement"-method by default, see [ParabolicSPDEs::simulateSPDEmodel].
#' @param repetitions A natural number giving the number of the Monte-Carlo repetitions.
#' @param ... for further arguments see [ParabolicSPDEs::simulateSPDEmodel].
#' @keywords Monte-Calro samples SPDE.
#' @references Bibinger, M. and Bossert, P. (2022) Efficient parameter estimation for parabolic SPDEs based on a log-linear model for realized volatilities,
#' Hildebrand, F. (2020) On generating fully discrete samples of the stochastic heat equation on an interval.
#' @export
#' @seealso [ParabolicSPDEs::simulateSPDEmodel].
#' @return A list of numeric \code{NxM} matrices.
#' @examples
#' repetitions = 10
#' theta0 = 0
#' theta1 = 1
#' theta2 = 1
#' sigma = 0.5
#' numberSpatialPoints = 10
#' numberTemporalPoints = 1000
#' l <- MCSPDESamples(repetitions,theta0,theta1,theta2,sigma,numberSpatialPoints,numberTemporalPoints)


MCSPDESamples <- function(repetitions, theta0,theta1,theta2,
                                    sigma,numberSpatialPoints,numberTemporalPoints,L=10,xi=function(x){0*x},cutoff = 10000,method = "replacement"){

  MC_list <- lapply(1:repetitions, function(i){
    simulateSPDEmodel( theta0,theta1,theta2,
                        sigma,numberSpatialPoints,numberTemporalPoints,L,xi,cutoff,method)
  })
  return(MC_list)
}


