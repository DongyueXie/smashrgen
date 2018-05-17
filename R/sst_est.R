#' Function to estimate \sqrt(\sigma^2+s_t^2) together
#' @param x: normal data
#' @param method: rmad or smashu
#' @return estimate standard deviation \sqrt(\sigma^2+s_t^2)
#' @export


sst_est=function(x,method,filter.number=1,family="DaubExPhase"){
  if(method=='rmad'){
    sigma.est=rmad(x,filter.number=filter.number,family=family)
  }
  if(method=='smashu'){
    sigma.est=sqrt(smash.gaus(x,v.est=T))
  }
  return(sigma.est)
}
