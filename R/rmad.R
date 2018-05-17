#' Running mad for estimating heteroscedastic variance
#' @param x: a vector of data
#' @param filter.number: as in smash.gaus
#' @param family: as n smash.gaus
#' @return a vector of estimated sd
#' @export

rmad=function(x,filter.number=1,family='DaubExPhase'){
  n = length(x)
  J = log2(n)
  x.w = wd(x, filter.number = filter.number, family = family, type = "station")
  win.size = round(n/10)
  odd.boo = (win.size%%2 == 1)
  win.size = win.size + (1 - odd.boo)
  sigma.est = runmad(accessD(x.w, J - 1), win.size, endrule = "func")
  return(sigma.est)
}
