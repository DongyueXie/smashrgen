#'@title init ebps
#'@param g_init a list of sigma2, and g_smooth
#'@param q_init a list of m, smooth

ebps_init = function(x,s,
                     make_power_of_2,
                     g_init,
                     q_init,
                     m_init_method
){
  n = length(x)

  # m_init, sigma2_init, smooth_init
  if(!is.null(g_init$sigma2)){
    if(any(is.na(g_init$sigma2))){
      g_init$sigma2 = NULL
    }
  }

  if(is.null(q_init$m)){
    if(m_init_method == 'smash_poi'){
      q_init$m = smash.poiss(x/s,log=TRUE)
    }
    if(m_init_method=='vga'){
      if(is.null(g_init$sigma2)){
        if(is.null(q_init$smooth)){
          fit_init = ebpm_normal(x,s,g_init = list(mean=log(sum(x)/sum(s)),var=NULL),fix_g = c(T,F))
        }else{
          fit_init = ebpm_normal(x,s,g_init = list(mean=q_init$smooth,var=NULL),fix_g = c(T,F))
        }
        q_init$m = fit_init$posterior$mean_log
        g_init$sigma2 = fit_init$fitted_g$var
      }else{
        if(is.null(q_init$smooth)){
          fit_init = ebpm_normal(x,s,g_init = list(mean=log(sum(x)/sum(s)),var = g_init$sigma2),fix_g = c(T,T))
        }else{
          fit_init = ebpm_normal(x,s,g_init=list(mean = q_init$smooth,var = g_init$sigma2),fix_g = c(T,T))
        }
        q_init$m = fit_init$posterior$mean_log
      }
    }
  }

  if(is.null(g_init$sigma2)){
    if(is.null(smooth_init)){
      g_init$sigma2 = ebpm_normal(x,s,g_init = list(mean=log(sum(x)/sum(s)),var=NULL),fix_g=c(T,F))$fitted_g$var
    }else{
      g_init$sigma2 = var(m_init - smooth_init)
    }
  }


  J = log(n,2)
  n_orig = n
  if(ceiling(J)!=floor(J)){
    #stop('Length of x must be power of 2')
    # reflect
    if(make_power_of_2=='reflect'){
      x = reflect(x)
      s = reflect(s)$x
      q_init$m = reflect(q_init$m)$x
      if(is.numeric(q_init$smooth)){
        q_init$smooth = reflect(q_init$smooth)$x
      }

      idx = x$idx
      x = x$x
      n = length(x)
      J = log(n,2)
    }
    if(make_power_of_2=='extend'){

      x = extend(x)
      s = extend(s)$x
      q_init$m = extend(q_init$m)$x
      if(is.numeric(q_init$smooth)){
        q_init$smooth = extend(q_init$smooth)$x
      }

      idx = x$idx
      x = x$x
      n = length(x)
      J = log(n,2)
    }
  }else{
    idx = 1:n
  }


  return(list(q_init = q_init,
              g_init=g_init,
              idx=idx,
              x=x,
              s=s))
}



#
#   if(!is.numeric(m_init)|length(m_init)!=n){
#     if(m_init == 'smash_poi'){
#       m_init = smash.poiss(x,log=TRUE) - log(s)
#     }else if(m_init == 'logx'){
#       m_init = log(x/s)
#       if(min(x)==0){
#         idx0 = (x == 0)
#         if(ash_pm_init_for0){
#           x_pm = ash_pois(x,scale=s,link='identity')$result$PosteriorMean
#           m_init[idx0] = log(x_pm[idx0])
#         }else{
#           if(eps_for0 == 'estimate'){
#             eps_for0 = sum(round(x)==1)/sum(round(x)<=1)+0.1
#           }
#           m_init[idx0] = log((x[idx0]+eps_for0)/s[idx0])
#         }
#       }
#     }else if(m_init == 'vga'){
#       if(is.null(sigma2_init)){
#         if(is.null(smooth_init)){
#           fit_init = ebpm_normal(x,s,g_init = list(mean=log(sum(x)/sum(s)),var=NULL),fix_g = c(T,F))
#         }else{
#           fit_init = ebpm_normal(x,s,g_init = list(mean=smooth_init,var=NULL),fix_g = c(T,F))
#         }
#         m_init = fit_init$posterior$mean_log
#         sigma2_init = fit_init$fitted_g$var
#       }else{
#         if(is.null(smooth_init)){
#           fit_init = ebpm_normal(x,s,g_init = list(mean=log(sum(x)/sum(s)),var = sigma2_init),fix_g = c(T,T))
#         }else{
#           fit_init = ebpm_normal(x,s,g_init=list(mean = smooth_init,var = sigma2_init),fix_g = c(T,T))
#         }
#         m_init = fit_init$posterior$mean_log
#       }
#     }else{
#       stop('unknown init method of mu')
#     }
#   }
#
#   if(is.null(sigma2_init)){
#     #sigma2_init = var(m - ti.thresh(m,method='rmad'))
#     if(is.null(smooth_init)){
#       sigma2_init = ebpm_normal(x,s,g_init = list(mean=log(sum(x)/sum(s)),var=NULL),fix_g=c(T,F))$fitted_g$var
#       #sigma2_init = var(m - smash.gaus(m))
#     }else{
#       sigma2_init = var(m_init - smooth_init)
#     }
#
#   }
#sigma2 = sigma2_init
#m = m_init
#v = rep(sigma2/2,n)
