

deal_formula = function(formula){
  tf <- as.character(formula)  
  tf <- tf[length(tf)]
  
  eval(parse(text=paste("c(",
                        ifelse(length(as.character(formula)) == 3,
                               'strsplit(as.character(formula)[2]," *[+] *"),',""),
                        'strsplit(strsplit(tf," *[|] *")[[1]]," *[+] *"))')))
}

.reshapeMatrix=function(x,n)
{
  return(matrix(rep.int(x,n),nrow(x),n))
}


.assertion_for_variables=function(data)
{
  if( is.null(data$y) )
  {
    cat( blue$bold("\n Y (")
         %+% black$bold("y")
         %+% blue$bold(") should not be NULL.\n\n") )
    
    return(NULL)
  }
  
  if( (!is.matrix(data$y))&(!is.data.frame(data$x)) )
  {
    cat( blue$bold("\n Y (")
         %+% black$bold("y")
         %+% blue$bold(") should be a matrix or dataframe.\n\n") )  
    return(NULL)
  }
  
  
  
  if( is.null(data$z) )
  {
    cat( blue$bold("\n Z (")
         %+% black$bold("z")
         %+% blue$bold(") should not be NULL.\n")
         %+% blue$bold(" If Z is null, please utilize linear models, such as ")
         %+% black$bold("lm() ")
         %+% blue$bold("function. \n\n")
    )
    
    return(NULL)
  }
  
  if((!is.matrix(data$z))&(!is.data.frame(data$z)))
  {
    cat( blue$bold("\n Z (")
         %+% black$bold("z")
         %+% blue$bold(") should be a matrix or dataframe.\n\n") )    
    return(NULL)
  }  
  
  
  if(!is.null(data$x))
  {
    if( (!is.matrix(data$x))&(!is.data.frame(data$x)) )
    {
      cat( blue$bold("\n X (")
           %+% black$bold("x")
           %+% blue$bold(") should be a matrix or dataframe.\n\n") )    
      return(NULL)
    }
  }
  
  
  if( !is.null(data$x) )
  {
    if(   length(data$y) !=  nrow(data$x) 
          #& nrow(data$x) > 1
    )
    {
      cat( blue$bold("\n The sample size of Y (")
           %+% black$bold("data$y")
           %+% blue$bold(") is not equal to that of X (")
           %+% black$bold("data$x")
           %+% blue$bold("). \n\n")
      )
      
      return(NULL)
    }
  }
  
  
  if( length(data$y) !=  nrow(data$z) )
  {
    cat( blue$bold("\n The sample size of Y (")
         %+% black$bold("data$y")
         %+% blue$bold(") is not equal to that of Z (")
         %+% black$bold("data$z")
         %+% blue$bold("). \n\n")
    )
    
    return(NULL)
  }
  
  
  return(TRUE)
  
}

.r_square=function(y,y_bar)
{
  y_ = mean(y)
  
  numerator = ( t(y - y_)%*%(y_bar - y_) )^2
  
  rs = numerator/(t(y - y_)%*%(y - y_)  *  t(y_bar - y_)%*%(y_bar - y_) )
  
  return(rs)
}
