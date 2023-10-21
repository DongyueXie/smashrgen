#'@title Apply ebps then run smash as a two-step procedure
#'@param x,s data vector and scaling factor. s can be a vector of the same length as x, or a scalar.
#'@param g_init a list of initial value of sigma2, and g_smooth. g_smooth is the initial prior g of the smoothing method. Can be NULL.
#'@param q_init a list of initial values of m, smooth. m is the posterior mean of mu, smooth the posterior mean of b. See the details below.
#'@param init_control See function ebps_init_control_default
#'@param general_control See function ebps_general_control_default
#'@param smooth_control See function ebps_smooth_control_default
#'@import smashr
#'@export
#'
#'

ebps_twostep = function(x,
                s = NULL,
                g_init = NULL,
                q_init = NULL,
                init_control = list(),
                general_control = list(),
                smooth_control = list(),
                homo = FALSE
){
  t_start = Sys.time()

  splitting_res = ebps(x,
                       s ,
                       g_init ,
                       q_init,
                       init_control,
                       general_control,
                       smooth_control)
  m = splitting_res$posterior$mean_log
  fit_smash = smash.gaus(m,homoskedastic = homo,v.est = T)
  t_end = Sys.time()
  return(list(posterior=list(mean_smooth = exp(fit_smash),
                             mean_log_smooth=fit_smash),
              log_likelihood = NULL,
              obj_trace = NULL,
              fitted_g = NULL,
              run_time = difftime(t_end,t_start,units='secs')))
}
