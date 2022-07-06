#' Estimate the parameters of a parabolic SPDE model
#'
#' Oracle and non-oralce estiamtions for SPDE models. All estimators are consistent and estimates the natrual parameters of the SPDE model. Non-Oracle estimations for parabolic SPDEs are based on a log-linear model, see references Bibinger, M. and Bossert, P. (2022).
#' @param data either a numerical \code{NxM} matrix or a list containing numerical \code{NxM} matrices. \code{N} denotes the temporal and \code{M} the spatial resolution of the grid.
#' @param estimationMethod a string indicating the parameter/parameters to be estimated. If only sigma is unknown choose "OracleSigma" and provide \code{theta1,theta2} respectively.
#' If \code{sigma and theta2} are known choose "OracleKappa" and provide the known parameters \code{sigma, theta2} or directly \code{sigma0_squared}.
#' If none of the parameters is known, choose "both". No additional information needs to be provided.
#' @param spatialDelta a real number greater than zero and less than 1/2 for selecting only the data points which are delta away from the Dirichlet boundary condition. Default is 0.05.
#' @param ... further arguments depending on the chosen estimation method. See \code{estiamtionmethod} or the examples below.
#' @keywords Parameter Estimation for SPDEs.
#' @references Bibinger, M. and Bossert, P. (2022) Efficient parameter estimation for parabolic SPDEs based on a log-linear model for realized volatilities,
#' Bibinger, M. and Trabs, M. (2017) Volatility estimation for stochastic PDEs using high-frequency observations,
#' Hildebrandt, F. and Trabs M. (2019) Parameter estimation for SPDEs based on discrete observations in time and space.
#' @export
#' @return a numeric vector. For method "OracleSigma" the returned value denotes the estimation of sigma^2. For method "OracleKappa" the returned value denotes the estimated quotient theta1/theta2.
#' For the method "both" two estimations are returned per data set, namely kappa = theta1/theta2 and sigma0_squared = sigma^2/(theta2)^(1/2).
#' See references for details on estimation methods.
#' @seealso [ParabolicSPDEs::simulateSPDEmodel], [ParabolicSPDEs::MCSPDESamples].
#' @examples
#' theta0 = 0
#' theta1 = 1
#' theta2 = 1
#' sigma = 0.5
#' numberSpatialPoints = 10
#' numberTemporalPoints = 1000
#' reputations = 10
#'
#' # Data
#' spde <- simulateSPDEmodel(theta0,theta1,theta2,sigma,numberSpatialPoints,numberTemporalPoints)
#' spde_list <- MCSPDESamples(reputations,theta0,theta1,theta2,sigma,numberSpatialPoints,numberTemporalPoints)
#'
#' # Estimation
#' estimateParametersSPDE(spde,estimationMethod = "OracleSigma",theta1=1,theta2=1)
#' estimateParametersSPDE(spde,estimationMethod = "OracleKappa",sigma=0.5,theta2=1)
#' estimateParametersSPDE(spde,estimationMethod = "both")
#'
#' estimateParametersSPDE(spde_list,estimationMethod = "OracleSigma",theta1=1,theta2=1)
#' estimateParametersSPDE(spde_list,estimationMethod = "OracleKappa",sigma=0.5,theta2=1)
#' estimateParametersSPDE(spde_list,estimationMethod = "both")










estimateParametersSPDE <- function(data, estimationMethod, spatialDelta = 0.05,theta1=NA,theta2=NA,sigma=NA,kappa=NA,sigma0_squared=NA){
  RV <- function(yPoint,y,dat){
    yIndex <- which(round(y,4) == round(yPoint,4))
    datPath <- dat[,yIndex]
    n <- length(datPath)
    sum <- (datPath[2]-datPath[1])^2
    for (i in 2:n) {
      sum <- sum + (datPath[i]-datPath[i-1])^2
    }
    return(sum)
  }
  estimator_sigmaOneSpatialPoint <- function(dat,yPoint,theta1,theta2){
    spatialPoints <- dim(dat)[2]-1
    y <- seq(0,1,1/spatialPoints)
    n <- dim(dat)[1]
    rv <- RV(yPoint = yPoint,y = y,dat = dat)
    scale <- sqrt(theta2*pi)*exp(yPoint*theta1/theta2)/sqrt(n)

    return(rv*scale)
  }
  # Check Assumptions by Bibinger, Bossert, Trabs
  if(is.list(data)){
    dat <- data[[1]]
    N <- dim(dat)[1] - 1
    M <- dim(dat)[2] - 1
    assumtionsTrue <- sqrt(N) >= M
  } else {
    dat <- data
    N <- dim(dat)[1] - 1
    M <- dim(dat)[2] - 1
    assumtionsTrue <- sqrt(N) >= M
  }

  if(assumtionsTrue){
    if(is.list(data)){
      require(parallel)
      require(pbapply)
      require(pbmcapply)
      numCores <- detectCores() - 1

      if(estimationMethod == "OracleSigma"){
        c1 <- is.na(theta1) && is.na(theta2)
        try(if(c1) stop("Provide theta1 and theta2!"))

        res <- pbmclapply(data, function(dat){
          M <- dim(dat)[2]-1
          spatialPoints <- M
          y <- seq(0,1,1/spatialPoints)
          yWithoutBounds <- y[which(round(y,4) >= round(spatialDelta,4))[1] : (which( round(y,4) > round((1 - spatialDelta),4))[1]-1 )]
          m <- length(yWithoutBounds)
          l <- lapply(yWithoutBounds, function(yy){
            estimator_sigmaOneSpatialPoint(dat,yy,theta1,theta2)
          })
          sum(unlist(l))/m
        },mc.cores = numCores)
        return(unlist(res))
      }
      if(estimationMethod == "OracleKappa"){
        c1 <- is.na(theta2) || is.na(sigma)
        try(if(c1 && is.na(sigma0_squared)) stop("Either provide theta2 and sigma or sigma0_squared!"))

        res <- lapply(data,function(dat){
          spatialPoints <- dim(dat)[2] - 1
          temporalPoints <- dim(dat)[1] - 1
          if(is.na(theta2)){
            sigma0_Squared <- sigma0_squared
          } else {
            sigma0_Squared <- sigma^2/sqrt(theta2)
          }
          try (if(sigma0_Squared<= 0) stop("theta2 and sigma must be greater 0!"))
          timeRep <- temporalPoints
          xrep <- spatialPoints
          time <- 1
          delta <- spatialDelta
          try (if(delta<0 || delta >= 0.5) stop("delta needs to be greater than 0 and less than 0.5!"))
          y <- seq(0,1,1/xrep)
          yWithoutBounds <- y[which(round(y,4) >= round(spatialDelta,4))[1] : (which( round(y,4) > round((time - spatialDelta),4))[1]-1 )]
          t <- seq(0,time, 1/timeRep)
          n <- dim(dat)[1]
          m <- length(yWithoutBounds)


          RV_list <- lapply(yWithoutBounds,function(yy){
            RV(yy,y,dat = dat)
          })
          RV_data <- unlist(RV_list)


          numeratorPart1_list <- lapply(1:m, function(k){
            log(RV_data[k]/sqrt(n))*yWithoutBounds[k]
          })
          numeratorPart1 <- sum(unlist(numeratorPart1_list))

          numeratorPart2_list <- lapply(yWithoutBounds, function(yy){
            yy*log(sigma0_Squared/(sqrt(pi)))
          })

          numeratorPart2 <- sum(unlist(numeratorPart2_list))
          numerator <- -numeratorPart1 + numeratorPart2
          denominator <- sum(yWithoutBounds^2)
          numerator / denominator
        })
        return(unlist(res))
      }
      if(estimationMethod == "both"){
        res <- pbmclapply(data,function(dat){
          delta <- spatialDelta
          try (if(delta<0 || delta >= 0.5) stop("delta needs to be greater than 0 and less than 0.5!"))
          spatialPoints <- dim(dat)[2] - 1
          temporalPoints <- dim(dat)[1] - 1
          y <- seq(0,1,1/spatialPoints)
          yWithoutBounds <- y[which(round(y,4) >= round(spatialDelta,4))[1] : (which( round(y,4) > round((1 - spatialDelta),4))[1]-1 )]
          t <- seq(0,1, 1/temporalPoints)
          n <- dim(dat)[1]
          m <- length(yWithoutBounds)

          RV_list <- lapply(yWithoutBounds,function(yy){
            log(RV(yy,y,dat)/sqrt(n))
          })
          RV_data <- unlist(RV_list)

          model <- lm(RV_data~yWithoutBounds)
          sigma_hat_squared <- exp(coef(model)[1])*sqrt(pi)
          names(sigma_hat_squared) <- "sigma^2_0"
          kappa_hat <- -coef(model)[2]
          names(kappa_hat) <- "kappa"

          return(c(sigma_hat_squared,kappa_hat))
        },mc.cores = numCores)
        return(unlist(res))
      }
    }
    else {
      if(estimationMethod == "OracleSigma"){
        c1 <- is.na(theta1) && is.na(theta2)
        try(if(c1) stop("Provide theta1 and theta2!"))

        dat <- data
        spatialPoints <- dim(dat)[2]-1
        y <- seq(0,1,1/spatialPoints)
        yWithoutBounds <- y[which(round(y,4) >= round(spatialDelta,4))[1] : (which( round(y,4) > round((1 - spatialDelta),4))[1]-1 )]
        m <- length(yWithoutBounds)
        l <- lapply(yWithoutBounds, function(yy){
          estimator_sigmaOneSpatialPoint(dat,yy,theta1,theta2)
        })
        return(sum(unlist(l))/m)
      }
      if(estimationMethod == "OracleKappa"){
        c1 <- is.na(theta2) || is.na(sigma)
        try(if(c1 && is.na(sigma0_squared)) stop("Either provide theta2 and sigma or sigma0_squared!"))
        dat <- data
        spatialPoints <- dim(dat)[2] - 1
        temporalPoints <- dim(dat)[1] - 1
        if(is.na(theta2)){
          sigma0_Squared <- sigma0_squared
        } else {
          sigma0_Squared <- sigma^2/sqrt(theta2)
        }
        try (if(sigma0_Squared<= 0) stop("theta2 and sigma must be greater 0!"))
        timeRep <- temporalPoints
        xrep <- spatialPoints
        time <- 1
        delta <- spatialDelta
        try (if(delta<0 || delta >= 0.5) stop("delta needs to be greater than 0 and less than 0.5!"))
        y <- seq(0,1,1/xrep)
        yWithoutBounds <- y[which(round(y,4) >= round(spatialDelta,4))[1] : (which( round(y,4) > round((time - spatialDelta),4))[1]-1 )]
        t <- seq(0,time, 1/timeRep)
        n <- dim(dat)[1]
        m <- length(yWithoutBounds)


        RV_list <- lapply(yWithoutBounds,function(yy){
          RV(yy,y,dat = dat)
        })
        RV_data <- unlist(RV_list)


        numeratorPart1_list <- lapply(1:m, function(k){
          log(RV_data[k]/sqrt(n))*yWithoutBounds[k]
        })
        numeratorPart1 <- sum(unlist(numeratorPart1_list))

        numeratorPart2_list <- lapply(yWithoutBounds, function(yy){
          yy*log(sigma0_Squared/(sqrt(pi)))
        })

        numeratorPart2 <- sum(unlist(numeratorPart2_list))
        numerator <- -numeratorPart1 + numeratorPart2
        denominator <- sum(yWithoutBounds^2)
        return(numerator / denominator)
      }
      if(estimationMethod == "both"){
        delta <- spatialDelta
        dat <- data
        try (if(delta<0 || delta >= 0.5) stop("delta needs to be greater than 0 and less than 0.5!"))
        spatialPoints <- dim(dat)[2] - 1
        temporalPoints <- dim(dat)[1] - 1
        y <- seq(0,1,1/spatialPoints)
        yWithoutBounds <- y[which(round(y,4) >= round(spatialDelta,4))[1] : (which( round(y,4) > round((1 - spatialDelta),4))[1]-1 )]
        t <- seq(0,1, 1/temporalPoints)
        n <- dim(dat)[1]
        m <- length(yWithoutBounds)

        RV_list <- lapply(yWithoutBounds,function(yy){
          log(RV(yy,y,dat)/sqrt(n))
        })
        RV_data <- unlist(RV_list)

        model <- lm(RV_data~yWithoutBounds)
        sigma_hat_squared <- exp(coef(model)[1])*sqrt(pi)
        names(sigma_hat_squared) <- "sigma^2_0"
        kappa_hat <- -coef(model)[2]
        names(kappa_hat) <- "kappa"

        return(c(sigma_hat_squared,kappa_hat))
      }
    }
  }
  else{
    print("Estimators by Trabs to be done!")
  }

}


