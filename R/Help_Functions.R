library(mvtnorm)
library(matlib)

GammaCalc <- function(n){
  I <- function(r){2*sqrt(r+1)-sqrt(r+2)-sqrt(r)}
  erg <- sum(I(0:n)^2)
  return(1/(pi)*(erg+2))
}


sumOfSquared_Scaling <- function(y){
  a <- sum(y^2)
  return(1/a)
}


doubleSumOfSquared_Scaling <- function(y){
  l2 <- lapply(y, function(y1){
    l1 <- lapply(y, function(y2){
      if(y1 != y2){}
      (y1-y2)^2
    })
    unlist(l1)
  })
  return(1/sum(unlist(l2)))
}

rgb2 <- function(r,g,b){
  return(rgb(r/255,g/255,b/255))
}

variance_kappa_numerator <- function(yWithoutBounds){
  y <- yWithoutBounds
  m <- length(y)
  list1 <- lapply(1:m, function(j){
    list2 <- lapply(1:m, function(l){
      y[l]-y[j]
    })
    sum(unlist(list2))^2
  })
  return(sum(unlist(list1)))
}

# realized volatility
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


asympVariance_etaMLE_sampling <- function(kappa,sigmaSquared,spatialDelta,kappaRange,sigmaRange){
  mu <- c(0,0)
  Gamma <- GammaCalc(100000)
  delta <- spatialDelta
  sigma11 <- function(sigmaSquared,delta){4*Gamma*pi*sigmaSquared^2*(1-delta+delta^2)/(1-2*delta)^2}
  sigma12 <- function(sigmaSquared,delta){6*Gamma*pi*sigmaSquared/(1-2*delta)^2}
  sigma22 <- function(delta){12*Gamma*pi/(1-2*delta)^2}
  
  sigma <- matrix(c(sigma11(sigmaSquared,delta),sigma12(sigmaSquared,delta),sigma12(sigmaSquared,delta),sigma22(delta)),byrow = T,ncol = 2)
  
  xMatrix <- as.matrix(expand.grid(sigmaRange, kappaRange))
  
  return(dmvnorm(x=xMatrix,mean=mu,sigma = sigma))
}

asympVariance_etaMLE <- function(kappa,sigmaSquared,spatialDelta){
  mu <- c(0,0)
  Gamma <- GammaCalc(100000)
  delta <- spatialDelta
  sigma11 <- function(sigmaSquared,delta){4*Gamma*pi*sigmaSquared^2*(1-delta+delta^2)/(1-2*delta)^2}
  sigma12 <- function(sigmaSquared,delta){6*Gamma*pi*sigmaSquared/(1-2*delta)^2}
  sigma22 <- function(delta){12*Gamma*pi/((1-2*delta)^2)}
  
  sigma <- matrix(c(sigma11(sigmaSquared,delta),sigma12(sigmaSquared,delta),sigma12(sigmaSquared,delta),sigma22(delta)),byrow = T,ncol = 2)
  
  
  
  return(sigma)
}


asympVariance_etaLSE_sampling <- function(kappa,sigmaSquared,spatialDelta,kappaRange,sigmaRange){
  delta <- spatialDelta
  c <- GammaCalc(10000)*pi*sigmaSquared^2
  f11 <- function(kappa,y){exp(-4*kappa*y)}
  f12 <- function(kappa,y){y*exp(-4*kappa*y)}
  f22 <- function(kappa,y){y^2*exp(-4*kappa*y)}
  g11 <- function(kappa,y){exp(-2*kappa*y)}
  g12 <- function(kappa,y){y*exp(-2*kappa*y)}
  g22 <- function(kappa,y){y^2*exp(-2*kappa*y)}
  
  u11 <- function(kappa,sigmaSquared,delta){integrate(function(y){f11(kappa,y)},delta,1-delta)$value}
  u12 <- function(kappa,sigmaSquared,delta){-sigmaSquared*integrate(function(y){f12(kappa,y)},delta,1-delta)$value}
  u22 <- function(kappa,sigmaSquared,delta){sigmaSquared^2*integrate(function(y){f22(kappa,y)},delta,1-delta)$value}
  
  v11 <- function(kappa,sigmaSquared,delta){integrate(function(y){g11(kappa,y)},delta,1-delta)$value}
  v12 <- function(kappa,sigmaSquared,delta){-sigmaSquared*integrate(function(y){g12(kappa,y)},delta,1-delta)$value}
  v22 <- function(kappa,sigmaSquared,delta){sigmaSquared^2*integrate(function(y){g22(kappa,y)},delta,1-delta)$value}
  
  U <- matrix(c(u11(kappa,sigmaSquared,delta),u12(kappa,sigmaSquared,delta),u12(kappa,sigmaSquared,delta),u22(kappa,sigmaSquared,delta)),ncol = 2,byrow = T)
  V <- matrix(c(v11(kappa,sigmaSquared,delta),v12(kappa,sigmaSquared,delta),v12(kappa,sigmaSquared,delta),v22(kappa,sigmaSquared,delta)),ncol = 2,byrow = T)
  erg <- inv(V) %*% U %*% inv(V)
  sigma <- c*erg
  print(sigma)
  xMatrix <- as.matrix(expand.grid(sigmaRange, kappaRange))
  
  return(dmvnorm(x=xMatrix,mean=c(0,0),sigma = sigma))
}

asympVariance_etaLSE <- function(kappa,sigmaSquared,spatialDelta){
  delta <- spatialDelta
  const <- GammaCalc(10000)*pi*sigmaSquared^2
  f11 <- function(kappa,y){exp(-4*kappa*y)}
  f12 <- function(kappa,y){y*exp(-4*kappa*y)}
  f22 <- function(kappa,y){y^2*exp(-4*kappa*y)}
  g11 <- function(kappa,y){exp(-2*kappa*y)}
  g12 <- function(kappa,y){y*exp(-2*kappa*y)}
  g22 <- function(kappa,y){y^2*exp(-2*kappa*y)}
  
  u11 <- function(kappa,sigmaSquared,delta){integrate(function(y){f11(kappa,y)},delta,1-delta)$value}
  u12 <- function(kappa,sigmaSquared,delta){-sigmaSquared*integrate(function(y){f12(kappa,y)},delta,1-delta)$value}
  u22 <- function(kappa,sigmaSquared,delta){sigmaSquared^2*integrate(function(y){f22(kappa,y)},delta,1-delta)$value}
  
  v11 <- function(kappa,sigmaSquared,delta){integrate(function(y){g11(kappa,y)},delta,1-delta)$value}
  v12 <- function(kappa,sigmaSquared,delta){-sigmaSquared*integrate(function(y){g12(kappa,y)},delta,1-delta)$value}
  v22 <- function(kappa,sigmaSquared,delta){sigmaSquared^2*integrate(function(y){g22(kappa,y)},delta,1-delta)$value}
  
  U <- matrix(c(u11(kappa,sigmaSquared,delta),u12(kappa,sigmaSquared,delta),u12(kappa,sigmaSquared,delta),u22(kappa,sigmaSquared,delta)),ncol = 2,byrow = T)
  V <- matrix(c(v11(kappa,sigmaSquared,delta),v12(kappa,sigmaSquared,delta),v12(kappa,sigmaSquared,delta),v22(kappa,sigmaSquared,delta)),ncol = 2,byrow = T)
  erg <- inv(V) %*% U %*% inv(V)
  sigma <- const*erg
  
  
  return(sigma)
}




convertPythonColor <- function(string){
  string <- gsub("rgb","rgb2",string)
  string <- gsub("'","",string)
  string <- gsub("\\[","c(",string)
  string <- gsub("\\]",")",string)
  return(print(string))
}


showColors_discrete <- function(colorPalette){
  n <- length(colorPalette)
  dat <- data.frame(x=seq(0,1,1/2),y=rep(1:n,each=length(seq(0,1,1/2))),group=as.factor(rep(1:n,each=length(seq(0,1,1/2)))))
  ggplot(dat,aes(x=x,y=y,group=group,color=group))+scale_color_manual(values = colorPalette)+geom_line(size=10)+theme_minimal()+labs(x="",y="")
}