#' A function estimate \sigma^2 using method 4
#'
#' @param x: data
#' @param mu
#' @param st
#' @param method: 'smashc', 'moment', 'wls', 'ols'
#' @return estimated sd
#' @export

sigma_est=function(x,mu,st,method){
  if(method=='moment'){
    x.m=c(x[length(x)],x,x[1])
    st.m=c(st[length(x)],st,st[1])
    sg=c()
    for(i in 2:length(x)){
      sg=c(sg,
           ((x.m[i]-x.m[i+1])^2+(x.m[i]-x.m[i-1])^2-2*st.m[i]^2-st.m[i-1]^2-st.m[i+1]^2)/4)
    }
    return(sqrt(mean(ifelse(sg<0,1e-8,sg))))
  }else{
    zt=(x-mu)^2
    zt2=zt-st^2
    zt2.var=2/3*zt^4
    if(method=='wls'){
      wls.est=lm(zt2~1,weights = 1/zt2.var)$coefficients
      return(sqrt(ifelse(wls.est<0,1e-8,wls.est)))
    }
    if(method=='smashc'){
      smash.est=smash.gaus(zt2,sigma=sqrt(zt2.var))
      return(sqrt(mean(ifelse(smash.est<0,1e-8,smash.est))))
    }
    if(method=='ols'){
      return(sqrt(ifelse(mean(zt2)<0,1e-8,mean(zt2))))
    }
  }
}
