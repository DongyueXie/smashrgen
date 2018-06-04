#' A funtion to reflect data
#' @param x: data vector
#' @param direct: left, right or both direction
#' @param len: the length to reflect
#' @return reflected vector
#' @export
reflect=function(x,direct='left',len){
  n=length(x)
  if(direct=='left'){
    x=c(rev(x[1:len]),x)
  }
  if(direct=='right'){
    x=c(x,rev(x[(n-len+1):n]))
  }
  if(direct=='both'){
    x=c(rev(x[1:len[1]]),x,rev(x[(n-len[2]+1):n]))
  }
  return(x)
}
