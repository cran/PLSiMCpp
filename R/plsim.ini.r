
#' @name plsim.ini
#' @aliases plsim.ini
#' @aliases plsim.ini.formula
#' @aliases plsim.ini.default
#' 
#' @title Initialize coefficients
#'
#' @description Xia \emph{et al.}'s MAVE method is used to obtain initialized 
#' coefficients \eqn{\alpha_0} and \eqn{\beta_0} for PLSiM 
#' \deqn{Y = \eta(Z^T\alpha) + X^T\beta + \epsilon}. 
#' 
#' @usage plsim.ini(\dots)
#' 
#' \method{plsim.ini}{formula}(formula, data, \dots)
#' 
#' \method{plsim.ini}{default}(xdat, zdat, ydat, Method="MAVE_ini", verbose = TRUE, \dots)
#' 
#' @param formula a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment containing the variables in the model.
#' @param xdat input matrix (linear covariates). The model reduces to a single index model when \code{x} is NULL.
#' @param zdat input matrix (nonlinear covariates). \code{z} should not be NULL.
#' @param ydat input vector (response variable).
#' @param Method string, optional (default="MAVE_ini").
#' @param verbose bool, default: TRUE. Enable verbose output.
#' @param \dots additional arguments.
#' 
#' @return
#' \item{zeta_i}{initial coefficients. \code{zeta_i[1:ncol(z)]} is the initial coefficient vector 
#' \eqn{\boldmath{\alpha}_0}, and \code{zeta_i[(ncol(z)+1):(ncol(z)+ncol(x))]} is the initial 
#' coefficient vector \eqn{\boldmath{\beta}_0}.}
#'
#' @export
#'
#' @examples
#' 
#' # EXAMPLE 1 (INTERFACE=FORMULA)
#' # To obtain initial values by using MAVE methods for partially
#' # linear single-index model.
#' 
#' n = 50
#' sigma = 0.1
#'
#' alpha = matrix(1,2,1)
#' alpha = alpha/norm(alpha,"2")
#' 
#' beta = matrix(4,1,1)
#' 
#' # Case1: Matrix Input
#' x = matrix(1,n,1)
#' z = matrix(runif(n*2),n,2)
#' y = 4*((z%*%alpha-1/sqrt(2))^2) + x%*%beta + sigma*matrix(rnorm(n),n,1)
#' 
#' zeta_i = plsim.ini(y~x|z)
#' 
#' # Case 2: Vector Input
#' x = rep(1,n)
#' z1 = runif(n)
#' z2 = runif(n) 
#' y = 4*((z%*%alpha-1/sqrt(2))^2) + x%*%beta + sigma*matrix(rnorm(n),n,1)
#' 
#' zeta_i = plsim.ini(y~x|z1+z2)
#' 
#' 
#' # EXAMPLE 2 (INTERFACE=DATA FRAME)
#' # To obtain initial values by using MAVE methods for partially
#' # linear single-index model.
#' 
#' n = 50
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
#' zeta_i = plsim.ini(xdat=X, zdat=Z, ydat=y)
#'
#' @references 
#' Y. Xia, W. Härdle. \emph{Semi-parametric estimation of partially linear single-index models}.
#' Journal of Multivariate Analysis, 2006, 97(5): 1162-1184.
#'
plsim.ini = function(...)
{
  args = list(...)
  if (is(args[[1]],"formula"))
    UseMethod("plsim.ini",args[[1]])
  else
    UseMethod("plsim.ini")
}

plsim.ini.formula = function(formula,data,...)
{
  mf = match.call(expand.dots = FALSE)   
  m = match(c("formula","data"), names(mf), nomatch = 0) 
  mf = mf[c(1,m)]
  
  mf.xf = mf
  
  mf[[1]] = as.name("model.frame")
  mf.xf[[1]] = as.name("model.frame")
  
  chromoly = deal_formula(mf[["formula"]])
  
  if (length(chromoly) != 3)
    stop("invoked with improper formula, please see plsim.ini documentation for proper use")
  
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
    
    ix = sort(c(which(orig.class == "ts"),which(orig.class != "ts")),index.return = TRUE)$ix
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
  
  zeta_i = plsim.ini(xdat = xdat, zdat = zdat, ydat = ydat, ...)
  
  return(zeta_i)
}


plsim.ini.default = function(xdat, zdat, ydat, 
                            Method="MAVE_ini", verbose = TRUE,...)
{
  if(verbose)
  {
    cat( blue$bold('\n Utilize the ')
         %+% black$bold('MAVE_ini')
         %+% blue$bold(' method to initialize coefficients.\n') )    
  }
  
  data = list(x=xdat, y=ydat, z=zdat)
  
  if ( is.null( .assertion_for_variables(data)) ) return(NULL)
  class(data) = Method
  
  x = data$z
  z = data$x
  y = data$y
  
  if(is.data.frame(x))
    x = data.matrix(x)
  
  if(is.data.frame(z))
    z = data.matrix(z)
  
  if(is.data.frame(y))
    y = data.matrix(y)
  
  n = nrow(x)
  dx = ncol(x)
  
  if( is.null(z) )
  {
    dz = 0
  }
  else
  {
    dz = ncol(z)
  }
  
  
  d_xz = dx + dz
  
  onep = matrix(1,dx,1)
  onen = matrix(1,n,1)
  
  B = diag(c(onep))
  nd = 1
  
  a = matrix(0,n,1)
  h2 = 2*n^(-2/(dx+4))
  
  if( !is.null(z) )
  {
    invzz = solve(t(z)%*%z)%*%t(z)
  }
  
  eyep1 = diag(c(matrix(1,dx+1,1)))/n^2
  
  
  for(iter in 1:dx)
  {
    
    ye = y - a
    
    if( !is.null(z) )
    {
      beta = invzz %*% ye
      ye = ye - z%*%beta
    }
    
    ab = matrix(1,dx,n)
    
    for(i in 1:n)
    {
      xi = x - t(.reshapeMatrix(matrix(x[i,]),n))
      kernel_xB = exp(-rowSums((xi%*%B)^2)/h2)
      onexi = cbind(onen,xi)
      xk = onexi*.reshapeMatrix(matrix(kernel_xB),dx+1)
      abi = solve(t(xk)%*%onexi+eyep1)%*%t(xk)%*%ye
      ab[,i] = abi[2:(dx+1)]
      a[i] = abi[1]
      
    }
    
    eigen_result = eigen(ab%*%t(ab))
    B0 = eigen_result$vectors
    D = eigen_result$values
    
    idx = order(D)
    B = B0
    B0 = B[,idx]
    
    if(dx == 1) B0 = matrix(B0)
    
    B[,1:nd] = B0[,dx:(dx-nd+1)]
    B = B[,1:max(1,dx-iter)]
    
  }
  
  alpha_i = B
  alpha_i = alpha_i/norm(alpha_i,"2")*sign(alpha_i[1])
  
  if( !is.null(z) ) beta_i = beta;
  
  if( !is.null(z) )
  {
    zeta_i = cbind(t(matrix(alpha_i)),t(matrix(beta_i)))
  }
  else
  {
    zeta_i = t(matrix(alpha_i))
  }
  
  return(zeta_i)
}
#
# plsim.ini.dbe=function(data)
# {
#   cat("Utilize the difference based method to initialize coeffiecnts.\n")
#
#   if ( is.null( .assertion_for_variables(data)) ) return(NULL)
#
#   x = data$x
#   z = data$z
#   y = data$y
#
#
#   n = nrow(z)
#   dz = ncol(z)
#
#   if( !is.null(x) )
#   {
#     dx = ncol(x)
#   }
#   else
#   {
#     dx = 0
#   }
#
#
#
#   yd = matrix(0,n,1)
#   zd = matrix(0,n,dz)
#
#   if( !is.null(x) )
#   {
#     xd = matrix(0,n,dx)
#   }
#
#   eyep1 = diag(c(matrix(1,dx+dz+1,1)))/n^2
#
#
#   for(k in 1:n)
#   {
#     zk = matrix(z[k,])
#     d = colSums((t(z) - zk%*%matrix(1,1,n))^2)
#     d[k] = max(d)
#     idx = which.min(d)
#
#     yd[k] = y[k] - y[idx]
#     zd[k,] = z[k,] - z[idx,]
#
#     if( !is.null(x) )
#     {
#       xd[k,] = x[k,] - x[idx,]
#     }
#
#   }
#
#   if( is.null(x) )
#   {
#     u = cbind(matrix(1,n,1),zd)
#   }
#   else
#   {
#     u = cbind(matrix(1,n,1),zd,xd)
#   }
#
#
#   b = solve(t(u)%*%u+eyep1)%*%t(u)%*%yd
#   alpha = matrix(b[2:(dz+1)])
#   alpha_i = alpha/norm(alpha,"2")*sign(alpha[1])
#
#   if( !is.null(x) )
#   {
#     beta_i = matrix(b[(dz+2):(dx+dz+1)])
#     zeta_i = cbind(t(matrix(alpha_i)),t(matrix(beta_i)))
#   }
#   else
#   {
#     zeta_i = t(matrix(alpha_i))
#   }
#
#
#   return(zeta_i)
# }
