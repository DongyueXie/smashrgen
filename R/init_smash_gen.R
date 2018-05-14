#' Initialize the input of smashgen algorithm

#' @param x: input data
#' @param ntri: Binomial number of trials
#' @param ashp: whether apply ash firstly
#' @param dist_family: family of ditribution: 'poisson','binomial'.
#' @return initilized m,s,y
#' @export


init_smash_gen=function(x,ntri,ashp,dist_family){
  n=length(x)
  if(dist_family=='poisson'){
    if(ashp){
      m0=ash(rep(0,n),1,lik=lik_pois(x))$result$PosteriorMean
    }else{
      m0=rep(mean(x),n)

    }
    s0=1/m0
    y0=log(m0)+(x-m0)/m0
    return(list(m0=m0,s0=s0,y0=y0))
  }
  if(dist_family=='binomial'){
    if(ashp){
      m0=ash(rep(0,n),1,lik=lik_binom(x,ntri))$result$PosteriorMean
    }else{
      m0=rep(mean(x/ntri),n)

    }
    s0=1/(ntri*m0*(1-m0))
    y0=log(m0/(1-m0))+(x/ntri-m0)/(m0*(1-m0))
    return(list(m0=m0,s0=s0,y0=y0))
  }
}


