
#' @name plsim.MAVE
#' @aliases plsim.MAVE
#' @aliases plsim.MAVE.formula
#' @aliases plsim.MAVE.default
#' 
#' @title Minimum Average Variance Estimation
#'
#'
#' @description MAVE (Minimum Average Variance Estimation), proposed by Xia \emph{et al.} (2006)
#' to estimate parameters in PLSiM
#' \deqn{Y=\eta(Z^T\alpha)+X^T\beta+\epsilon.}
#' 
#' @usage plsim.MAVE(\dots)
#' 
#' \method{plsim.MAVE}{formula}(formula, data, \dots)
#' 
#' \method{plsim.MAVE}{default}(xdat=NULL, zdat, ydat, h=NULL, zeta_i=NULL, maxStep=100,
#' tol=1e-8, iniMethods="MAVE_ini", ParmaSelMethod="SimpleValidation", TestRatio=0.1, 
#' K = 3, seed=0, verbose=TRUE, \dots)
#'
#' @param formula a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment containing the variables in the model.
#' @param xdat input matrix (linear covariates). The model reduces to a single index model when \code{x} is NULL.
#' @param zdat input matrix (nonlinear covariates). \code{z} should not be NULL.
#' @param ydat input vector (response variable).
#' @param h a numerical value or a vector for bandwidth. If \code{h} is NULL, a default 
#' vector c(0.01,0.02,0.05,0.1,0.5) will be given. 
#' \link{plsim.bw} is employed to select the optimal bandwidth when \code{h} is a vector or NULL.
#' @param zeta_i initial coefficients, optional (default: NULL). It could be obtained by the function \code{\link{plsim.ini}}.
#' \code{zeta_i[1:ncol(z)]} is the initial coefficient vector \eqn{\boldmath{\alpha}_0},
#' and \code{zeta_i[(ncol(z)+1):(ncol(z)+ncol(x))]} is the initial coefficient vector \eqn{\boldmath{\beta}_0}.
#' @param maxStep the maximum iterations, default: 100.
#' @param tol convergence tolerance, default: 1e-8.
#' @param iniMethods string, optional (default: "SimpleValidation"). 
#' @param ParmaSelMethod the parameter for the function \link{plsim.bw}.
#' @param TestRatio the parameter for the function \link{plsim.bw}.
#' @param K the parameter for the function \link{plsim.bw}.
#' @param seed int, default: 0.
#' @param verbose bool, default: TRUE. Enable verbose output.
#' @param \dots additional arguments.
#'
#' @return
#' \item{eta}{estimated non-parametric part \eqn{\hat{\eta}(Z^T\boldmath{\hat{\alpha} })}.}
#' \item{zeta}{estimated coefficients. \code{zeta[1:ncol(z)]} is \eqn{\hat{\alpha}}, 
#' and \code{zeta[(ncol(z)+1):(ncol(z)+ncol(x))]} is \eqn{\hat{\beta}}.}
#' \item{data}{data information including \code{x}, \code{z}, \code{y}, bandwidth \code{h}, 
#' initial coefficients \code{zetaini} and iteration step \code{MaxStep}.}
#' \item{y_hat}{ \code{y}'s estimates.}
#' \item{mse}{mean squares erros between \code{y} and \code{y_hat}.}
#' \item{variance}{variance of \code{y_hat}.}
#' \item{r_square}{multiple correlation coefficient.}
#' \item{Z_alpha}{ \eqn{Z^T\boldmath{\hat{\alpha}}}.}
#'
#' @export
#'
#' @examples
#'
#' # EXAMPLE 1 (INTERFACE=FORMULA)
#' # To estimate parameters in partially linear single-index model using MAVE.
#' 
#' n = 30
#' sigma = 0.1
#'
#' alpha = matrix(1,2,1)
#' alpha = alpha/norm(alpha,"2")
#'
#' beta = matrix(4,1,1)
#'
#' x = matrix(1,n,1)
#' z = matrix(runif(n*2),n,2)
#' y = 4*((z%*%alpha-1/sqrt(2))^2) + x%*%beta + sigma*matrix(rnorm(n),n,1)
#'
#' fit = plsim.MAVE(y~x|z, h=0.1)
#' 
#' # EXAMPLE 2 (INTERFACE=DATA FRAME)
#' # To estimate parameters in partially linear single-index model using MAVE.
#' 
#' n = 30
#' sigma = 0.1
#'
#' alpha = matrix(1,2,1)
#' alpha = alpha/norm(alpha,"2")
#' beta = matrix(4,1,1)
#' 
#' x = rep(1,n)
#' z1 = runif(n)
#' z2 = runif(n) 
#' X = data.frame(x)
#' Z = data.frame(z1,z2)
#' 
#' x = data.matrix(X)
#' z = data.matrix(Z)
#' y = 4*((z%*%alpha-1/sqrt(2))^2) + x%*%beta + sigma*matrix(rnorm(n),n,1)
#' 
#' fit = plsim.MAVE(xdat=X, zdat=Z, ydat=y, h=0.1)
#'
#' @references
#'
#' Y. Xia, W. Härdle. \emph{Semi-parametric estimation of partially linear single-index models}.
#' Journal of Multivariate Analysis, 2006, 97(5): 1162-1184.
#'

plsim.MAVE = function(...)
{
  args = list(...)
  if (is(args[[1]],"formula"))
    UseMethod("plsim.MAVE",args[[1]])
  else
    UseMethod("plsim.MAVE")
}

plsim.MAVE.formula = function(formula,data,...)
{
  mf = match.call(expand.dots = FALSE)   
  m = match(c("formula","data"),
            names(mf), nomatch = 0) 
  mf = mf[c(1,m)]
  
  mf.xf = mf
  
  mf[[1]] = as.name("model.frame")
  mf.xf[[1]] = as.name("model.frame")
  
  chromoly = deal_formula(mf[["formula"]])
  
  if (length(chromoly) != 3)
    stop("invoked with improper formula, please see plsim.MAVE documentation for proper use")
  
  bronze = lapply(chromoly, paste, collapse = " + ")
  
  mf.xf[["formula"]] = as.formula(paste(" ~ ", bronze[[2]]),
                                  env = environment(formula))
  
  mf[["formula"]] = as.formula(paste(bronze[[1]]," ~ ", bronze[[3]]),
                               env = environment(formula))
  
  formula.all = terms(as.formula(paste(" ~ ",bronze[[1]]," + ",bronze[[2]], " + ",bronze[[3]]),
                                 env = environment(formula)))
  
  orig.class = if (missing(data))
    sapply(eval(attr(formula.all, "variables"), environment(formula.all)),class)
  else sapply(eval(attr(formula.all, "variables"), data, environment(formula.all)),class)
  
  arguments.mfx = chromoly[[2]]
  arguments.mf = c(chromoly[[1]],chromoly[[3]])
  
  
  mf[["formula"]] = terms(mf[["formula"]])
  mf.xf[["formula"]] = terms(mf.xf[["formula"]])
  
  if(all(orig.class == "ts")){
    arguments = (as.list(attr(formula.all, "variables"))[-1])
    attr(mf[["formula"]], "predvars") = bquote(.(as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments)))))[,.(match(arguments.mf,arguments)),drop = FALSE])
    attr(mf.xf[["formula"]], "predvars") = bquote(.(as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments)))))[,.(match(arguments.mfx,arguments)),drop = FALSE])
  }else if(any(orig.class == "ts")){
    arguments = (as.list(attr(formula.all, "variables"))[-1])
    arguments.normal = arguments[which(orig.class != "ts")]
    arguments.timeseries = arguments[which(orig.class == "ts")]
    
    ix <- sort(c(which(orig.class == "ts"),which(orig.class != "ts")),index.return = TRUE)$ix
    attr(mf[["formula"]], "predvars") = bquote((.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])[,.(match(arguments.mf,arguments)),drop = FALSE])
    attr(mf.xf[["formula"]], "predvars") = bquote((.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])[,.(match(arguments.mfx,arguments)),drop = FALSE])
  }
  
  
  mf = tryCatch({
    eval(mf,parent.frame())
  },error = function(e){
    NULL
  })
  
  mf.xf = tryCatch({
    eval(mf.xf,parent.frame())
  },error = function(e){
    NULL
  })
  
  if(is.null(mf)){
    cat( blue$bold("\n Z (")
         %+% black$bold("z")
         %+% blue$bold(") should not be NULL.\n")
         %+% blue$bold(" If Z is null, please utilize linear models, such as ")
         %+% black$bold("lm() ")
         %+% blue$bold("function. \n\n")
    )
    return(NULL)
  }
  else{
    ydat = model.response(mf)
  }
  
  xdat = mf.xf
  zdat = mf[, chromoly[[3]], drop = FALSE]
  
  ydat = data.matrix(ydat)
  if(!is.null(xdat) & is.null(dim(xdat[,1]))){
    xdat = data.matrix(xdat)
  }
  else if(!is.null(dim(xdat[,1]))){
    xdat = xdat[,1]
  }
  
  if(is.null(dim(zdat[,1]))){
    zdat = data.matrix(zdat)
  }
  else{
    zdat = zdat[,1]
  }
  
  
  fit = plsim.MAVE(xdat = xdat, zdat = zdat, ydat = ydat, ...)
  return(fit)
}

plsim.MAVE.default = function(xdat=NULL,zdat,ydat,h=NULL,zeta_i=NULL,maxStep=100,tol=1e-8,iniMethods="MAVE_ini",
                        ParmaSelMethod="SimpleValidation",TestRatio=0.1,K = 3,seed=0,verbose=TRUE,...)
{
  
  data = list(x=xdat,y=ydat,z=zdat)
  if ( is.null( .assertion_for_variables(data)) ) return(NULL)
  
  x = data$x
  y = data$y
  z = data$z
  
  if(is.data.frame(x))
    x = data.matrix(x)
  
  if(is.data.frame(z))
    z = data.matrix(z)
  
  if(is.data.frame(y))
    y = data.matrix(y)
  
  tmp = x
  x = z
  z = tmp
  
  
  n = nrow(x)
  dx = ncol(x)
  
  if( is.null(z) ) dz = 0 else dz = ncol(z)
  
  d_xz = dx + dz
  
  if( is.null(zeta_i) )
  {
    #cat( blue$bold("\n Adopt the ")
    #     %+%black$bold("MAVE_ini")
    #     %+%blue$bold(" method to initialize the coeffients. \n\n")
    #)
    zeta_i = plsim.ini(z,x,y,verbose=verbose)
  }
  else
  {
    if( length(zeta_i) < d_xz)
    {
      cat( blue$bold("\n The dimension of the coeffients (")
           %+%black$bold("zeta_i")
           %+%blue$bold(") is not right. \n\n")
      )
      
      return(NULL)
    }
  }
  
  
  if( is.null(h) )
  {
    #hList = c(0.01,0.02,0.05,0.1,0.5)
    
    res_MAVE_simple = plsim.bw(z,x,y,
                                    TargetMethod="MAVE",zeta_i=zeta_i,
                                    ParmaSelMethod=ParmaSelMethod,TestRatio=TestRatio,
                                    K=K,seed=seed,verbose=verbose)
    h = res_MAVE_simple$bandwidthBest
  }
  
  if( length(h) > 1)
  {
    res_MAVE_simple = plsim.bw(z,x,y,bandwidthList=h,
                                    TargetMethod="MAVE",zeta_i=zeta_i,
                                    ParmaSelMethod=ParmaSelMethod,TestRatio=TestRatio,
                                    K=K,seed=seed,verbose=verbose)
    h = res_MAVE_simple$bandwidthBest    
  }
  
  
  
  theta = zeta_i[1:dx]
  
  if(!is.null(z)) beta = zeta_i[(dx+1):d_xz]
  
  
  h2 = 2*h*h
  eyepq = diag(c(matrix(1,d_xz,1)))/n^2
  
  step1 = 0
  delta1 = 2
  
  while( (step1<maxStep) & (delta1>tol) )
  {
    step1 = step1 + 1
    
    theta_old = theta
    
    if(!is.null(z)) beta_old = beta
    
    tmp = .reshapeMatrix(x%*%theta,n)
    d = (tmp - t(tmp))^2
    ker = exp(-d/h2)
    ker = ker/.reshapeMatrix(matrix(rowSums(ker)),n)
    
    md = matrix(0,n,dx*dx)
    mc = matrix(0,n,dx)
    me = matrix(0,n,dx)
    
    D22 = 0
    C2 = 0
    mcz = matrix(0,n,dz)
    mdd = matrix(0,n,dx*dz)
    
    for(i in 1:n)
    {
      tmp = x - t(.reshapeMatrix(matrix(x[i,]),n))
      tmp1 = .reshapeMatrix(matrix(ker[,i]),dx)*tmp
      md[i,] = matrix(c(t(t(tmp1)%*%tmp)),1)
      mc[i,] = colSums(tmp1)
      me[i,] = t(y)%*%tmp1
      
      if( !is.null(z) )
      {
        z1 = .reshapeMatrix(matrix(ker[,i]),dz)*z
        D22 = D22 + t(z1)%*%z
        mcz[i,] = colSums(z1)
        C2 = C2 + t(y)%*%z1
        
        mdd[i,] = matrix(c(t(tmp)%*%z1),1,dx*dz)
      }
      
    }
    
    
    step2 = 0
    delta2 = 2
    
    
    
    while( (step2 < maxStep) & (delta2>tol) )
    {
      
      theta_old2 = theta
      if(!is.null(z)) beta_old2 = beta
      
      
      step2 = step2 + 1
      
      if(!is.null(z)) ye = y - z%*%beta else ye = y
      
      
      
      tmp = .reshapeMatrix(matrix(x%*%theta),n)
      d1 = tmp - t(tmp)
      ker1 = d1*ker
      s2 = matrix(rowSums(d1*ker1))
      s1 = matrix(rowSums(ker1))
      s = matrix(rowSums(ker))
      
      d = s2*s - s1*s1 + 1/n^2
      a = (ker%*%ye*s2 - ker1%*%ye*s1)/d
      b = (ker1%*%ye*s - ker%*%ye*s1)/d
      D = matrix(c(t(t(b*b)%*%md)),dx,dx)
      C = t(b)%*%me - t(b*a)%*%mc
      
      if(!is.null(z))
      {
        Cz = C2 - t(a)%*%mcz
        D12 = matrix(c(t(t(b)%*%mdd)),dx,dz)
        
        D = rbind(cbind(D,D12),cbind(t(D12),D22))
        C = cbind(C,Cz)
      }
      
      theta = matrix(solve(D+eyepq)%*%t(C))
      
      
      if(!is.null(z)) beta = theta[(dx+1):d_xz]
      
      theta = theta[1:dx]
      
      if(!is.null(theta))
      {
        if(theta[1]==0){
          si = theta[2]
        }
        else{
          si = theta[1]
        }
        theta = matrix(sign(theta[1])*theta/as.vector(sqrt(t(theta)%*%theta)))
      }
      
      
      if(!is.null(z))
      {
        delta2 = max( abs(cbind(t(matrix(theta)),t(matrix(beta))) - cbind(t(matrix(theta_old2)),t(matrix(beta_old2)))) )
      }
      else
      {
        delta2 = max( abs(theta - theta_old2) )
      }
      
      
      
      
    }
    
    if(!is.null(z))
    {
      delta1 = max( abs(cbind(t(matrix(theta)),t(matrix(beta))) - cbind(t(matrix(theta_old)),t(matrix(beta_old)))) )
    }
    else
    {
      delta1 = max( abs(theta - theta_old) )
    }
    
    
    
  }
  
  
  
  y_hat = a
  
  if(!is.null(z))
  {
    zeta = cbind(t(matrix(theta)),t(matrix(beta)))
    y_hat = y_hat + z%*%matrix(beta)
    SiMflag = 0
  }
  else
  {
    zeta = theta
    SiMflag = 1
  }
  
  #result = list()
  #result$eta = a
  #result$zeta = zeta
  #return(result)
  
  
  data = list(x=z,y=y,z=x,h=h,
              zetaini=zeta_i,MaxStep=maxStep,SiMflag = as.logical(SiMflag))
  
  
  fit = list(zeta=zeta,data=data,eta=a,
             Z_alpha=x%*%matrix(theta),
             variance=var(y_hat),r_square=.r_square(y,y_hat),
             mse=sum((y-y_hat)^2)/nrow(y),
             y_hat=y_hat)
  
  class(fit) = "pls"
  
  return(fit)
}

#
# MAVE=function(x,z,y,h,zeta_i=NULL,maxStep=100,tol=1e-8,iniMethods="MAVE")
# {
#
#   if( is.null(x) )
#   {
#     cat("x could be NULL. And it will be devoloped in the next release.")
#     return(NULL)
#   }
#
#   data = list(x=x,y=y,z=z)
#   if ( is.null( .assertion_for_variables(data)) ) return(NULL)
#
#
#   tmp = x
#   x = z
#   z = tmp
#
#
#   n = nrow(x)
#   dx = ncol(x)
#   dz = ncol(z)
#
#   d_xz = dx + dz
#
#   if( is.null(zeta_i) )
#   {
#     data_t = list(x=z,z=x,y=y)
#     class(data_t) = iniMethods
#
#     cat( blue$bold("\n Adopt the ")
#          %+%black$bold("MAVE_ini")
#          %+%blue$bold(" method to initialize the coeffients. \n\n")
#     )
#
#     zeta_i = plsim.ini(data_t)
#   }
#   else
#   {
#     if( length(zeta_i) < d_xz)
#     {
#       cat( blue$bold("\n The dimension of the coeffients (")
#            %+%black$bold("zeta_i")
#            %+%blue$bold(") is not right. \n\n")
#       )
#
#       return(NULL)
#     }
#   }
#
#
#
#   theta = zeta_i[1:dx]
#   beta = zeta_i[(dx+1):(dx+dz)]
#
#
#
#
#
#   h2 = 2*h*h
#   eyepq = diag(c(matrix(1,d_xz,1)))/n^2
#
#   step1 = 0
#   delta1 = 2
#
#   while( (step1<maxStep) & (delta1>tol) )
#   {
#     step1 = step1 + 1
#
#     theta_old = theta
#     beta_old = beta
#
#     tmp = .reshapeMatrix(x%*%theta,n)
#     d = (tmp - t(tmp))^2
#     ker = exp(-d/h2)
#     ker = ker/.reshapeMatrix(matrix(rowSums(ker)),n)
#
#     md = matrix(0,n,dx*dx)
#     mc = matrix(0,n,dx)
#     me = matrix(0,n,dx)
#
#     D22 = 0
#     C2 = 0
#     mcz = matrix(0,n,dz)
#     mdd = matrix(0,n,dx*dz)
#
#     for(i in 1:n)
#     {
#       tmp = x - t(.reshapeMatrix(matrix(x[i,]),n))
#       tmp1 = .reshapeMatrix(matrix(ker[,i]),dx)*tmp
#       md[i,] = matrix(c(t(t(tmp1)%*%tmp)),1)
#       mc[i,] = colSums(tmp1)
#       me[i,] = t(y)%*%tmp1
#
#       z1 = .reshapeMatrix(matrix(ker[,i]),dz)*z
#       D22 = D22 + t(z1)%*%z
#       mcz[i,] = colSums(z1)
#       C2 = C2 + t(y)%*%z1
#
#       mdd[i,] = matrix(c(t(tmp)%*%z1),1,dx*dz)
#
#     }
#
#
#     step2 = 0
#     delta2 = 2
#
#     while( (step2 < maxStep) & (delta2>tol) )
#     {
#
#
#       theta_old2 = theta
#       beta_old2 = beta
#
#       step2 = step2 + 1
#
#       ye = y - z%*%beta
#       tmp = .reshapeMatrix(matrix(x%*%theta),n)
#       d1 = tmp - t(tmp)
#       ker1 = d1*ker
#       s2 = matrix(rowSums(d1*ker1))
#       s1 = matrix(rowSums(ker1))
#       s = matrix(rowSums(ker))
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
#       d = s2*s - s1*s1 + 1/n^2
#       a = (ker%*%ye*s2 - ker1%*%ye*s1)/d
#       b = (ker1%*%ye*s - ker%*%ye*s1)/d
#       D = matrix(c(t(t(b*b)%*%md)),dx,dx)
#       C = t(b)%*%me - t(b*a)%*%mc
#
#       Cz = C2 - t(a)%*%mcz
#       D12 = matrix(c(t(t(b)%*%mdd)),dx,dz)
#
#       D = rbind(cbind(D,D12),cbind(t(D12),D22))
#       C = cbind(C,Cz)
#
#       theta = matrix(solve(D+eyepq)%*%t(C))
#       beta = theta[(dx+1):d_xz]
#       theta = theta[1:dx]
#       theta = matrix(sign(theta[1])*theta/sqrt(t(theta)%*%theta))
#
#       delta2 = max( abs(cbind(t(matrix(theta)),t(matrix(beta))) - cbind(t(matrix(theta_old2)),t(matrix(beta_old2)))) )
#
#
#     }
#
#     delta1 = max( abs(cbind(t(matrix(theta)),t(matrix(beta))) - cbind(t(matrix(theta_old)),t(matrix(beta_old)))) )
#
#
#   }
#
#
#
#
#
#   zeta = cbind(t(matrix(theta)),t(matrix(beta)))
#
#   result = list()
#   result$eta = a
#   result$zeta = zeta
#
#
#   return(result)
# }
