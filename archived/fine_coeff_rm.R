#' A function that sets the finest level coefficience to zero
#'
#' @param: y0
#' @return wavelet
#' @export

fine_coeff_rm=function(y0,family,filter.number){
  wds=wd(y0,family = family,filter.number = filter.number)
  wtd=threshold(wds, levels = wds$nlevels-1,  policy="manual",value = Inf)
  return(wr(wtd))
}
