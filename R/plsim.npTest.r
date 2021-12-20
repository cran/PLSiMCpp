
#' @name plsim.npTest
#' 
#' @title Testing nonparametric component
#' 
#' @description Study the hypothesis test:
#' \deqn{H_0:\eta(u) = \theta_0+\theta_1u \ \mbox{ versus }\quad H_1:\ \eta(u)\ne \theta_0 + 
#' \theta_1u \ \mbox{for \ some \ } u }
#' where \eqn{\theta_0} and \eqn{\theta_1} are unknown constant parameters.
#' 
#' 
#' @param fit the result of function \link{plsim.est} or \link{plsim.vs.soft}.
#' @param significant_level significant level for the test, default: 0.05. 
#' 
#' @return
#' \item{T2}{test statistic.}
#' \item{pvalue}{p-value.}
#' \item{\eqn{delta}}{ \eqn{delta=1} means to reject the hypothesis \eqn{H_0}, while \eqn{\delta=0} means to accept the hypothesis \eqn{H_0}. }
#' 
#' @export
#' 
#' @examples
#' 
#' n = 50
#' sigma = 0.1
#'
#' alpha = matrix(1,2,1)
#' alpha = alpha/norm(alpha,"2")
#'
#' beta = matrix(4,1,1)
#'
#' x = matrix(1,n,1)
#'
#' z = matrix(runif(n*2),n,2)
#'
#' y = 4*((z%*%alpha-1/sqrt(2))^2) + x%*%beta + sigma*matrix(rnorm(n),n,1)
#'
#' # Obtain parameters in PLSiM using Profile Least Squares Estimator
#' fit_plsimest = plsim.est(x, z, y)
#' 
#' res_npTest_plsimest = plsim.npTest(fit_plsimest)
#' 
#' # Obtain parameters in PLSiM using Penalized Profile Least Squares Estimator
#' # with lambda set as 0.01
#' fit_plsim = plsim.vs.soft(x,z,y,lambda = 0.01)
#' 
#' res_npTest_plsim = plsim.npTest(fit_plsim)
#' 
#' @references
#'
#' H. Liang, X. Liu, R. Li, C. L. Tsai. \emph{Estimation and testing for partially linear single-index models}.
#' Annals of statistics, 2010, 38(6): 3811.
#' 

plsim.npTest=function(fit,significant_level=0.05)
{
  data = fit$data

  x = data$x
  y = data$y
  z = data$z

  h = data$h
  lambda = data$lambda
  zetaini = data$zetaini
  MaxStep = data$MaxStep
  flag = data$SiMflag
  useSCAD = data$useSCAD

  zeta = fit$zeta
  mse = fit$mse

  n = nrow(z)
  if(is.null(x)) dx= 0 else dx = ncol(x)
  dz = ncol(z)

  #under null hypothesis
  alpha = zeta[1:dz]
  beta = zeta[(dz+1):(dz+dx)]

  y_H0 = y-x%*%matrix(beta)
  x_H0 = cbind(z%*%matrix(alpha),matrix(1,n,1))
  theta = solve(t(x_H0)%*%x_H0)%*%t(x_H0)%*%y_H0
  RSS_H0 = sum((y-(x_H0%*%theta+x%*%matrix(beta)))^2)


  #under alternative hypothesis
  RSS_H1 = mse*n

  #For epa kernel
  rK = 0.2308
  dfn = rK*2*0.15/h


  T2 = rK/2*n*(RSS_H0-RSS_H1)/RSS_H1

  pvalue = pchisq(T2,df=dfn,lower.tail=F)
  
  if(pvalue<significant_level)
  {
    h = 1
  }
  else
  {
    h = 0
  }  

  return(list(T2=T2,pvalue=pvalue,delta=h))
}
