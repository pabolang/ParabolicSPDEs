#' Test for the curvature patameter in a parabolic SPDE model
#'
#' Aymptotic oracle and non-oracle tests for the curvature parameter in SPDE models. Non-oracle curvature tests are based on a log-linear model, see references Bibinger, M. and Bossert, P. (2022).
#' @param data a numerical \code{NxM} matrix. \code{N} denotes the temporal and \code{M} the spatial resolution of the grid.
#' @param curvature hypothesized curvature of the random field. Default is 0 (no curvature). Alternatively you can enter both hypothesized parameters \code{theta1} and \code{theta2}.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Default is \code{two.sided}.
#' @param conf.level confidence level for the returned confidence interval.
#' @param method a string indicating if a non-orcale or an oracle test is performed. Must either be \code{non-oracle} or \code{oracle}. When \code{orcale} is used, you must specify the normalized voltility by either enter \code{sigma0_squared} or the parameters \code{sigma} and \code{theta2}. Default is \code{non-oracle}.
#' @param spatialDelta a real number greater than zero and less than 1/2 for selecting only the data points which are delta away from the Dirichlet boundary condition. The default is 0.05.
#' @param ... further arguments depending on the chosen estimation method. See \code{estiamtionmethod} or the examples below.
#' @keywords Curvature Tests for SPDEs.
#' @references Bibinger, M. and Bossert, P. (2022) Efficient parameter estimation for parabolic SPDEs based on a log-linear model for realized volatilities.
#' @export
#' @seealso [ParabolicSPDEs::simulateSPDEmodel], [ParabolicSPDEs::estimateParametersSPDE].
#' @examples
#' # Data
#' theta0 = 0
#' theta1 = 0
#' theta2 = 1
#' sigma = 0.5
#' numberSpatialPoints = 100
#' numberTemporalPoints = 10000
#' spde <- simulateSPDEmodel(theta0,theta1,theta2,sigma,numberSpatialPoints,numberTemporalPoints)
#'
#'
#' kappa_test(spde,method = "non-oracle",curvature = -1,alternative = "g")






kappa_test <- function(data, curvature = 0,alternative = "two.sided",conf.level = 0.95,method = "non-oracle",spatialDelta = 0.05, sigma = NA, theta2 = NA, sigma0_squared = NA, theta1 = NA){
  dat <- data
  N <- dim(dat)[1] - 1
  M <- dim(dat)[2] - 1
  assumptionsTrue <- sqrt(N) >= M
  if(assumptionsTrue){
    DNAME <- deparse(substitute(data))
    if(alternative == "t"){
      alternative = "two.sided"
    }
    if(alternative == "l"){
      alternative = "less"
    }
    if(alternative == "g"){
      alternative = "greater"
    }
    
    if(!is.na(theta1) & !is.na(theta2)){
      curvature = theta1/theta2
    }
    
    if(method == "oracle"){
      c1 <- is.na(theta2) || is.na(sigma)
      try(if(c1 && is.na(sigma0_squared)) stop("Either provide theta2 and sigma or sigma0_squared!"))
      if(is.na(theta2)){
        sigma0_squared <- sigma0_squared
      } else {
        sigma0_squared <- sigma^2/sqrt(theta2)
      }
    }
    
    
    Gamma <- GammaSPDE()
    
    y <- seq(0,1,1/M)
    yWithoutBounds <- y[which(round(y,4) >= round(spatialDelta,4))[1] : (which( round(y,4) > round((1 - spatialDelta),4))[1]-1 )]
    m <- length(yWithoutBounds)
    
    if(method == "non-oracle"){
      est1 <- estimateParametersSPDE(data,estimationMethod = "both")
      est <- est1[2]
      a <- (1-2*spatialDelta)*sum(yWithoutBounds^2)/m
      b <- ((1-2*spatialDelta)*mean(yWithoutBounds))^2/(1-2*spatialDelta)
      u <- sqrt(N*m)/sqrt(Gamma*pi*(1-2*spatialDelta)/(a-b))
      val <- (est-curvature)*u
      if(alternative == "two.sided" ){
        if(val >= 0){
          PVal <- unname(pnorm(-val) + 1-pnorm(val))
        } else {
          PVal <- unname(pnorm(val)+1-pnorm(-val))
        }
      } else {
        if(alternative == "less" ){
          PVal <- pnorm(val)
        } else {
          if(alternative == "greater" ){
            PVal <- 1-pnorm(val)
          } else {
            print("Invalid alternative!")
          }
        }
      }
      
      v <- 1/u
      q <- pnorm(1-conf.level/2)
      ci <- c(est-q*v,est+q*v)
      attr(ci, "conf.level") <- conf.level
      
      names(est) <- names(curvature) <- "Curvature"
      names(val) <- "Lambda"
      sigma0_squared <- est1[1]
      names(sigma0_squared) <- "normalized volatility estimate"
      structure(list(statistic = val,
                     data.name = DNAME,
                     conf.int = ci,
                     parameter = sigma0_squared,
                     p.value = PVal,
                     estimate = est,
                     null.value = curvature,
                     alternative = alternative,
                     method = "Asymptotic Kappa Non-Oracle Test"),class = "htest")
      
      
    } else {
      if(method == "oracle"){
        est <- estimateParametersSPDE(data,sigma0_squared = sigma0_squared,estimationMethod = "OracleKappa")
        u <- sqrt(N*m)/sqrt(Gamma*pi/(sum(yWithoutBounds^2)/m))
        val <- (est-curvature)*u
        if(alternative == "two.sided" ){
          if(val >= 0){
            PVal <- unname(pnorm(-val) + 1-pnorm(val))
          } else {
            PVal <- unname(pnorm(val)+1-pnorm(-val))
          }
        } else {
          if(alternative == "less" ){
            PVal <- pnorm(val)
          } else {
            if(alternative == "greater" ){
              PVal <- 1-pnorm(val)
            } else {
              print("Invalid alternative!")
            }
          }
        }
        
        v <- 1/u
        q <- pnorm(1-conf.level/2)
        ci <- c(est-q*v,est+q*v)
        attr(ci, "conf.level") <- conf.level
        
        names(est) <- names(curvature) <- "Curvature"
        names(val) <- "Upsilon"
        names(sigma0_squared) <- "normalized volatility"
        structure(list(statistic = val,
                       parameter = sigma0_squared,
                       data.name = DNAME,
                       conf.int = ci,
                       p.value = PVal,
                       estimate = est,
                       null.value = curvature,
                       alternative = alternative,
                       method = "Asymptotic Kappa Oracle Test"),class = "htest")
        
      } else {
        print("No valid method. Check spelling or documentation.")
      }
    }
  } else{
    try(if(assumtionsTrue) stop("Assumptions are violated! Make sure N>=M^2"))
    print("Estimators by Trabs to be done!")
  }
}
