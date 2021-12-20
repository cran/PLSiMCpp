#' @aliases summary.pls

summary.pls=function(object,...)
{
  data = object$data
  
  if(data$SiMflag == 0)
  {
    cat( blue$bold("\n Single Index Model"))
  }
  else
  {
    cat( blue$bold("\n Partial Linear Single-index Model "))
  }
  
  n = nrow(data$y)
  
  if(is.null(data$x))
  {
    dx = 0
  }
  else
  {
    dx = ncol(data$x)
  }
  
  dz = ncol(data$z)
  d = dx + dz
  
  

  cat( blue$bold("\n Regression Data: ")
       %+% blue$bold(as.character(n))
       %+% blue$bold(" training data points, in ")
       %+% blue$bold(as.character(d))
       %+% blue$bold(" variable(s)")
       )
  
  if( is.null(colnames(data$z)) )
  {
    colnames(data$z) = 1:dz
  }
  
  cat('\n         ')
  cat( blue$bold(paste(colnames(data$z), collapse= "\t")))
  cat( blue$bold("\n Alpha: "))  
  cat( blue$bold(paste( round(object$zeta[1:dz],4),collapse= "\t")) )
  
  if(!is.null(data$x))
  {
    if( is.null(colnames(data$x)) )
    {
      colnames(data$x) = 1:dx
    }
    
    cat('\n        ')
    cat( blue$bold(paste(colnames(data$x), collapse= "\t")))
    cat( blue$bold("\n Beta: ")) 
    cat( blue$bold(paste( round(object$zeta[(dz+1):(dz+dx)],4),collapse= "\t")) )
  }
  
  cat('\n\n')
  cat( blue$bold(' Mean Square Error: ')
       %+% blue$bold(round(object$mse,6))
     )
  cat('\n')
  cat( blue$bold(' R-squared: ')
       %+% blue$bold(round(object$r_square,4))
       )
  
}