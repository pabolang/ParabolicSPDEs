#' Simulate a parabolic stochastic partial differential equations in one space dimension
#'
#' Simulate a sample of a SPDE model in one space dimension on a discrete grid with \code{N} temporal and \code{M} spatial points, where the grid points are equidistant within the unit square. The initial condition is set to be zero. For simulating SPDE samples with a general initial condition use the build-in cut-off method. Note that this method results in a systematic bias and dramatically increases computational costs. The SPDE model is using the Dirichlet boundary condition.
#' @param theta0 a real number which controls the drift of the solution field.
#' @param theta1 a real number which controls the curvature of the solution field.
#' @param theta2 a real number greater than zero which reduces the noise level of param sigma and the curvature effect of param theta1
#' @param sigma a real number greater than zero which controls the overall noise level of the solution field
#' @param numberOfSpatialPoints number of equidistant spatial points M in \code{[0,1]}.
#' @param numberOfTemporalPoints number of equidistant temporal points N in \code{[0,1]}.
#' @param L a natural number indicating the replacement bound LM dependent on multiples of \code{M}. Default is \code{L=10}.
#' @param xi initial condition. Default is \code{xi = 0}. For general initial condition choose method = "cutoff".
#' @param cutoff a natural number for the cut-off frequency of the Fourrier series. Only used when method = "cutoff".
#' @param method either "replacement" (Default) or "cutoff". Note, that the replacement method only allows the initial condition to be zero. For general initial conditions choose "cutoff"-method.
#' @keywords Sample of one SPDE
#' @references Bibinger, M. and Bossert, P. (2022) Efficient parameter estimation for parabolic SPDEs based on a log-linear model for realized volatilities,
#' Hildebrand, F. (2020) On generating fully discrete samples of the stochastic heat equation on an interval.
#' @export
#' @seealso [ParabolicSPDEs::MCSPDESamples], [ParabolicSPDEs::plotSPDE],[ParabolicSPDEs::estimateParametersSPDE].
#' @return A numeric \code{NxM} matrix.
#' @examples
#' theta0 = 0
#' theta1 = 1
#' theta2 = 1
#' sigma = 0.5
#' numberSpatialPoints = 10
#' numberTemporalPoints = 1000
#' spde <- simulateSPDEmodel(theta0,theta1,theta2,sigma,numberSpatialPoints,numberTemporalPoints)





simulateSPDEmodel <- function(theta0,theta1,theta2,sigma,numberSpatialPoints,numberTemporalPoints,L=10,xi=function(x){0*x},cutoff = 10000,method = "replacement"){
  require(parallel)
  require(pbapply)
  require(pbmcapply)
  numCores <- detectCores() - 1

  try (if(theta2 <= 0) stop("theta2 must be greater 0!"))
  try (if(sigma <= 0) stop("sigma must be greater 0!"))


  if(method == "replacement"){
    M <- numberSpatialPoints
    y <- seq(0,1,1/M)
    N <- numberTemporalPoints
    Delta <- 1/N
    kappa <- theta1/theta2
    LM <- L*M
    Gamma <- theta1^2/(4*theta2^2)-theta0/theta2
    Gamma0 <- sqrt(abs(Gamma))



    ## Step 1 ##
    Im <- function(m){
      ImPlus <- function(m){
        l <- lapply(0:LM, function(l){
          m+2*l*M
        })
        return(unlist(l))
      }
      ImMinus <- function(m){
        l <- lapply(0:LM, function(l){
          2*M-m+2*l*M
        })
        return(unlist(l))
      }
      erg <- c(ImPlus(m),ImMinus(m))
      return(subset(erg,0<erg & erg < LM))
    }

    el <- function(l,y){
      return(sqrt(2)*sin(pi*l*y)*exp(-kappa*y/2))
    }

    lambdal <- function(l){
      return(pi^2*theta2*l^2+theta1^2/(4*theta2)-theta0)
    }


    # in Diss x_k
    ul <- function(l){
      erg <- c(0,rep(NA,N))
      random <- rnorm(N)
      for (i in 1:N) {
        erg[i+1] <- exp(-lambdal(l)*Delta)*erg[i]+sigma*sqrt((1-exp(-2*lambdal(l)*Delta))/(2*lambdal(l)))*random[i]
      }
      return(erg)
    }

    ## Step 2 ##
    sm <- function(m){

      bm <- function(m){
        erg <- matrix(sqrt(2)*sin(pi*m*y),ncol=1)
        return(erg)
      }

      rho <- function(x,y){
        if(x<=y){
          if(Gamma<0){
            erg <- sin(Gamma0*(1-y))*sin(Gamma0*x)/(Gamma0*sin(Gamma0))
          } else{
            if (Gamma == 0){
              erg <- x*(1-y)
            } else {
              erg <- sinh(Gamma0*(1-y))*sinh(Gamma0*x)/(Gamma0*sinh(Gamma0))
            }
          }

        } else {
          x1 <- x
          y1 <- y
          x <- y1
          y <- x1
          if(Gamma<0){
            erg <- sin(Gamma0*(1-y))*sin(Gamma0*x)/(Gamma0*sin(Gamma0))
          } else{
            if (Gamma == 0){
              erg <- x*(1-y)
            } else {
              erg <- sinh(Gamma0*(1-y))*sinh(Gamma0*x)/(Gamma0*sinh(Gamma0))
            }
          }
        }
        return(erg*sigma^2/(2*theta2))
      }

      Sigma <- function(y){
        l1 <- lapply(y, function(y1){
          l2 <- lapply(y, function(y2){
            rho(y1,y2)
          })
          unlist(l2)
        })
        erg <- do.call(rbind,l1)
        return(matrix(erg,ncol=(M+1),nrow=(M+1),byrow=T))
      }

      rest <- function(m){
        I <- Im(m)
        l <- lapply(I, function(l){
          sigma^2/(2*lambdal(l))
        })
        erg <- sum(unlist(l))
        return(erg)
      }

      erg <- 1/M^2*( t(bm(m)) %*% Sigma(y) %*% bm(m) ) - rest(m)
      return(erg)

    }

    Rm <- function(m){
      sm <- sqrt(sm(m))
      erg <- c(0,rnorm(N,0,sm))
      return(erg)
    }


    ## Step 3 ##
    Um <- function(m){
      I <- Im(m)
      l <- lapply(I, function(l){
        ul(l)
      })
      ul_wholeTime <- do.call(rbind,l)
      erg <- colSums(ul_wholeTime)

      erg2 <- Rm(m)
      return(erg+erg2)
    }

    # zeit in Zeil, spalten sind alle m
    l <- pbmclapply(1:(M-1),function(m){
      Um(m)
    },mc.cores = 7)
    Um_data <- do.call(rbind,l)

    l <- pbmclapply(1:(M-1),function(m){
      el(m,y)
    },mc.cores = 2)
    em_data <- do.call(rbind,l)



    ## Output ##
    l1 <- pbmclapply(1:(M+1), function(y1){
      l2 <- lapply(1:(N+1), function(t1){
        l3 <- lapply(1:(M-1), function(m){
          Um_data[m,t1]*em_data[m,y1]
        })
        sum(unlist(l3))
      })
      unlist(l2)
    },mc.cores = 7)
    spde <- do.call(cbind,l1)
    return(spde)
  } else {
    if(method == "cutoff"){

      try (if(theta2 <= 0) stop("theta2 must be greater 0!"))
      try (if(sigma <= 0) stop("sigma must be greater 0!"))

      xrep <- numberSpatialPoints
      time <- 1
      timeRep <- numberTemporalPoints
      numberOfSum <- cutoff

      if(theta1 == 0){
        e_k <- function(k,y,theta1,theta2){ sqrt(2) * sin(pi * k * y)}
      } else {
        e_k <- function(k,y,theta1,theta2){ sqrt(2) * sin(pi * k * y) * exp(- y * theta1 / (2*theta2))}
      }

      lambda_k <- function(k,theta0,theta1,theta2){-theta0 + theta1^2 / (4*theta2) + theta2*pi^2*k^2}

      y <- seq(0,1,1/xrep)
      t <- seq(0,time,1/timeRep)

      isXiZero <- integrate(function(y){abs(xi(y))},0,1)$value == 0

      OUProcess <- function(k,t,theta0,theta1,theta2,sigma,xi,isXiZero = F){


        e_k <- function(k,y,theta1,theta2){ sqrt(2) * sin(pi * k * y) * exp(- y * theta1 / (2*theta2))}
        lambda_k <- function(k,theta0,theta1,theta2){-theta0 + theta1^2 / (4*theta2) + theta2*pi^2*k^2}


        if(isXiZero){
          x0 <- function(k){0*k}
        } else {
          x0 <- function(k){integrate(function(y){exp(y*theta1/theta2)*e_k(k,y,theta1,theta2)*xi(y)},0,1,subdivisions=2000,stop.on.error=F)$value}
        }

        normal <- rnorm((length(t) - 1))
        Delta <- t[2]
        x <- rep(NA,length(t))
        x[1] <- x0(k)


        for(i in 2:(length(t))){
          x[i] <- x[i-1]*exp(-lambda_k(k,theta0,theta1,theta2)*Delta) + sigma*sqrt((1-exp(-2*lambda_k(k,theta0,theta1,theta2)*Delta))/(2*lambda_k(k,theta0,theta1,theta2)))*normal[i-1]
        }
        return(x)


      }


      # Vertikal zeit; horizontal k
      x_kMatrix <- do.call(rbind,pbmclapply(1:numberOfSum,function(k){
        OUProcess(k,t,theta0,theta1,theta2,sigma,xi,isXiZero)
      },mc.cores = numCores)
      )

      e_kMatrix <- do.call(rbind, pbmclapply(1:numberOfSum,function(k){
        l <- lapply(y, function(yy){
          e_k(k,yy,theta1,theta2)
        })
        unlist(l)
      }))

      # Berechnung der Reihe mit Abbruch
      l <- pbmclapply(1:length(y),function(yy){

        lTime <- lapply(2:length(t), function(tt){

          lSum <- lapply(1:numberOfSum, function(kk){
            e_kMatrix[kk,yy]*x_kMatrix[kk,tt]
          })

          sum(unlist(lSum))
        })

        unlist(lTime)
      },mc.cores = numCores)

      z <- do.call(cbind,l)
      x <- seq(0,1,1/xrep)
      x0 <- xi(x)
      z <- rbind(x0,z)
      return(z)

    }
  }
}

