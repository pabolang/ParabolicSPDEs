# ParabolicSPDEs
Simulate and plot parabolic SPDE models as introduced among others by [Bibinger, M. and Bossert, P. (2022) Efficient parameter estimation for parabolic SPDEs based on a log-linear model for realized volatilities](https://arxiv.org/abs/2207.00357). This R-package also provides estimation methods for estimating the natural parameters in this SPDE model.

## Installation

```r
require("devtools")
devtools::install_github('pabolang/ParabolicSPDEs')
```

---


## Usage
We consider the following parabolic stochastic partial differential equation 

$$ \text{d}X_t(y)=\bigg(\vartheta_2\frac{\partial^2 X_t(y)}{\partial y^2}+\vartheta_1\frac{\partial X_t(y)}{\partial y}+\vartheta_0X_t(y)\bigg)\text{d}t+\sigma\text{d}B_t(y),$$

where $\vartheta_0\in\mathbb{R}, \vartheta_1\in\mathbb{R},\vartheta_2>0,\sigma>0$ 
and a cylindrical Brownian motion $B_t(y)$ within a Sobolev space 
on the unit interval $\[0,1\]$. 
Further, we consider a Dirchilet boundary condition and an initial condition $\xi(y)=0$, 
where $(t,y)\in \[0,1\]^2$.


By using this package, we can simply simulate a SPDE model on a equidistant discrete $N\times M$ grid, 
where $N$ denotes the temporal and $M$ the spatial resolution, using the function `simulateSPDEmodel`:
```r
library(ParabolicSPDEs)
theta0 = 0
theta1 = 1
theta2 = 1
sigma = 0.5
numberSpatialPoints = 100
numberTemporalPoints = 10000

(spde <- simulateSPDEmodel(theta0,theta1,theta2,sigma,numberSpatialPoints,numberTemporalPoints))
```
The function `simulateSPDEmodel` returns a $N\times M$ matrix which we can be plotted using the following function:
```r
plotSPDE(spde)
```

<img width="709" alt="spde_plot" src="https://user-images.githubusercontent.com/78961989/177564860-d90f651e-763a-41b3-aa10-b3162fb571a2.png">


For creating multiple SPDE samples, use the function `MCSPDESamples`. 

This package also includes the function `estimateParametersSPDE` for estimating the parameters of a SPDE model. So far, only parametric estimators under the assumption $N\geq \sqrt{M}$ has been implemented. For more details on the estimators and the assumptions, see the references below. The function uses  the argument `estimationMethod` which includes an oracle estimation for the parameter $\sigma$ 
and the natural parameter $\kappa=\vartheta_1/\vartheta_2$, 
as well as an estimation for the parameter $(\sigma_0^2,\kappa)$, 
where $\sigma_0^2=\sigma^2/\sqrt{\vartheta_2}$.
```r
estimateParametersSPDE(spde, estimationMethod = "OracleSigma", theta1=1,theta2 = 1)
estimateParametersSPDE(spde, estimationMethod = "OracleKappa", sigma=0.5,theta2 = 1)
estimateParametersSPDE(spde, estimationMethod = "OracleKappa", sigma0_squared = 0.5^2/sqrt(1)))
estimateParametersSPDE(spde, estimationMethod = "both")
```
This function also supports a list of $N\times M$ matrices and returns the estimated parameters for each matrix respectively. 
Therefore, it is possible to create for example density plots for estimating the natural parameters of a SPDE:
```r
spde_list1 <- MCSPDESamples(reputations = 100,theta0 = 0,theta1 = 1,theta2 = 1,sigma = 0.5, numberSpatialPoints = 10, numberTemporalPoints = 1000)
spde_list2 <- MCSPDESamples(reputations = 100,theta0 = 0,theta1 = 1.1,theta2 = 1,sigma = 1, numberSpatialPoints = 10, numberTemporalPoints = 1000)
est1 <- estimateParametersSPDE(spde_list1, estimationMethod = "OracleKappa", theta2 = 1, sigma= 0.5)
est2 <- estimateParametersSPDE(spde_list2, estimationMethod = "OracleKappa", theta2 = 1, sigma = 1)

require(ggplot2)
dat <- data.frame(x=c(est1,est2),group=rep(1:2,each=length(est1)))
ggplot(dat,aes(x=x,fill=group,group=group,color=group))+
  geom_density(alpha=0.4)+
  theme_minimal()+
  labs(x="")+
  theme(legend.position = "none")
```

![dens](https://user-images.githubusercontent.com/78961989/177566297-9ae4c448-88fd-43ea-a3bb-21356528ae6a.png)

---

## To-Do's
- Implement estimators for violation of the assumption $N\geq \sqrt{M}$ according to Hildebrand and Trabs.
---

## References
[Bibinger, M. and Bossert, P. (2022) Efficient parameter estimation for parabolic SPDEs based on a log-linear model for realized volatilities.](https://arxiv.org/abs/2207.00357)

[Hildebrandt, F. and Trabs, M. (2019) Parameter estimation for SPDEs based on discrete observations in time and space](https://arxiv.org/abs/1910.01004)

[Hildebrandt, F. (2020) On generating fully discrete samples of the stochastic heat equation on an interval.](https://arxiv.org/abs/2001.03403)

[Bibinger, M. and Trabs, M. (2017) Volatility estimation for stochastic PDEs using high-frequency observations.](https://arxiv.org/abs/1710.03519)
