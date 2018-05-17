#' A function estimate \sigma^2 using method 4
#'
#' @param x: data
#' @param mu
#' @param st
#' @param method: 'eb', 'moment', 'wls', 'mle','huber'
#' @param var_est: method to estimate variance: 'rmad', 'smash', 'default'
#' @param k: parameter in huber m estimator
#' @return estimated sd
#' @export

sigma_est=function(x,mu,st,method,var_est='smash',k=NULL,family='DaubExPhase',filter.number=1){
  if(method=='moment'){
    x.m=c(x[length(x)],x,x[1])
    st.m=c(st[length(x)],st,st[1])
    sg=c()
    for(i in 2:length(x)){
      sg=c(sg,
           ((x.m[i]-x.m[i+1])^2+(x.m[i]-x.m[i-1])^2-2*st.m[i]^2-st.m[i-1]^2-st.m[i+1]^2)/4)
    }
    var.est=mean(sg)
    #return(sqrt(ifelse(moment.est<0,1e-8,moment.est)))
  }else{
    zt=(x-mu)^2
    zt2=zt-st^2
    if(var_est=='smash'){
      zt2.var=smash.gaus(zt2,v.est=T,family=family,filter.number = filter.number)
    }else if(var_est=='rmad'){
      zt2.var=(rmad(zt2,filter.number,family))^2
    }else{
      zt2.var=2/3*zt^2
    }
    if(method=='wls'){
      var.est=lm(zt2~1,weights = 1/zt2.var)$coefficients
      #return(sqrt(ifelse(wls.est<0,1e-8,wls.est)))
    }
    if(method=='huber'){
      if(is.null(k)){k=quantile(zt2,0.99)}
      var.est=hubers(zt2,k,s=sqrt(zt2.var))$mu
      #return(sqrt(mean(ifelse(smash.est<0,1e-8,smash.est))))
    }
    if(method=='mle'){
      var.est=mean(zt2)
      #return(sqrt(ifelse(mean(zt2)<0,1e-8,mean(zt2))))
    }
    if(method=='eb'){
      var.est=mean(ebnm_normal(zt2,sqrt(zt2.var))$posmean)
      #return(sqrt(mean(ifelse(ash.est<0,1e-8,ash.est))))
    }
  }
  return(sqrt(ifelse(var.est<0,1e-8,var.est)))
}
