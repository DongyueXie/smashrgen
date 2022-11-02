#' Update function for smashgen algo
#' @return mt,st,yt
#' @export

update_smash_gen=function(mu.hat,x,ntri,dist_family){
  #n=length(x)
  if(dist_family=='poisson'){
    mt=exp(mu.hat)
    return(list(mt=mt,st=1/mt,yt=log(mt)+(x-mt)/mt))
  }
  if(dist_family=='poi_binom'){
    #mt is the estimated p
    mt=exp(mu.hat)
    #this is p*ntri, for update purpose only
    mt.up=mt*ntri
    return(list(mt=mt,st=1/mt.up,yt=log(mt.up)+(x-mt.up)/mt.up))
  }
  if(dist_family=='binomial'){
    mt=exp(mu.hat)
    mt=mt/(1+mt)
    return(list(mt=mt,st=1/(ntri*mt*(1-mt)),yt=log(mt/(1-mt))+(x/ntri-mt)/(mt*(1-mt))))
  }
}
