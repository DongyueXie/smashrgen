#' Function to estimate \sqrt(\sigma^2+s_t^2) together
#' @param x: normal data
#' @param method: rmad or smashu
#' @return estimate standard deviation \sqrt(\sigma^2+s_t^2)
#' @export


sst_est=function(x,method,filter.number=1,family="DaubExPhase"){
  n = length(x)
  J = log2(n)
  if(method=='rmad'){
    x.w = wavethresh::wd(x, filter.number = filter.number,
                         family = family, type = "station")
    win.size = round(n/10)
    odd.boo = (win.size%%2 == 1)
    win.size = win.size + (1 - odd.boo)
    sigma.est = runmad(accessD(x.w, J - 1), win.size, endrule = "func")
    return(sigma.est)
  }
  if(method=='smashu'){
    sigma.est=sqrt(smash.gaus(x,v.est=T,joint=T)$var.res)
    return(sigma.est)
  }
}
