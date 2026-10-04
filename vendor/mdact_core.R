# Derived from author results_in_paper/MDACT_code/funcs.R (GPL-2).
# MDACT statistical formulas retained; F00 uses an equivalent analytic area.
mdact_null_cdf <- function(t, a, b, c) {
  # P(a*U + b*V + c*max(U,V)^2 <= t), independent U,V ~ Uniform(0,1).
  # Splitting at c*v^2+(a+b)*v=t gives the same conditional integral as F00.
  stopifnot(length(t) == 1L, is.finite(t), a > 0, b > 0, c > 0)
  if (t <= 0) return(0)
  if (t >= a + b + c) return(1)
  split <- 2 * t / (a + b + sqrt((a + b)^2 + 4 * c * t))
  plateau <- max(0, min(split, (t - c - a) / b))
  area <- plateau
  if (split > plateau) {
    d0 <- max(0, a^2 + 4 * c * (t - b * plateau))
    d1 <- max(0, a^2 + 4 * c * (t - b * split))
    # Factor the difference of powers to avoid subtracting close d^(3/2).
    average_root <- (2 / 3) * (d0 + sqrt(d0 * d1) + d1) / (sqrt(d0) + sqrt(d1))
    area <- area + (split - plateau) * (average_root - a) / (2 * c)
  }
  end <- min(1, 2 * t / (b + sqrt(b^2 + 4 * c * t)))
  if (end > split)
    area <- area + (end - split) / a *
      (t - b * (end + split) / 2 - c * (end^2 + end * split + split^2) / 3)
  min(1, max(0, area))
}

balancing_DACT_control_DR_adjust = function(p.M,p.Y,estws,significance_upper,control.method){
  sim.num <- length(p.M)
  
  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est
  
  pi.01.est.sd = pi.01.est/c
  pi.10.est.sd = pi.10.est/c
  pi.00.est.sd = pi.00.est/c
  
  #Ts <- pi.01.est*p.M + pi.10.est*p.Y + pi.00.est*pmax(p.M,p.Y)^2
  
  ss <- DR_DACT_thr_adjust(p.M, p.Y,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method)
  
  return(ss)
}

DR_DACT_thr_adjust <- function(p.M, p.Y,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method){
  sim.num <- length(p.M)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- max(1 - c, 0)
  
  Ts <- pi.01.est*p.M + pi.10.est*p.Y + pi.00.est*pmax(p.M,p.Y)^2
  
  F00 <- function(t) mdact_null_cdf(t, pi.01.est, pi.10.est, pi.00.est)
  
  F01 <- function(t){
    c0 <- pi.00.est * p.Y^2 + pi.01.est * p.Y
    c1 <- t- pi.10.est * p.Y
    c2 <- (t - pi.00.est * p.Y^2 - pi.10.est * p.Y)/pi.01.est
    m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.Y))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)
    
    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.Y)
  }
  
  F10 <- function(t){
    c0 <- pi.00.est * p.M^2 + pi.10.est * p.M
    c1 <- t- pi.01.est * p.M
    c2 <- (t - pi.00.est * p.M^2 - pi.01.est * p.M)/pi.10.est
    m2 <- (-pi.10.est+sqrt(pmax(0,pi.10.est^2 + 4*pi.00.est*(t-pi.01.est*p.M))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)
    
    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.M)
  }
  
  if(control.method == "FDR"){
    thr <- function(t) {
      c * (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) + 
                              pi.01.est/(pi.01.est+pi.11.est)*F01(t) + 
                              (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) - 
                                 pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t))) / max(mean(Ts < t), 1 / sim.num) - significance_upper
    }
    
    #t_s <- uniroot(thr, c(min(c(p.M, p.Y)), 1))$root
    sorted_index <- order(Ts) # Get the indices that sort 'Ts'
    lower_bound <- 1
    upper_bound <- length(Ts)
    
    while (upper_bound - lower_bound > 1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]
      
      x <- tryCatch(thr(mid_value), error = function(e) e)
      k <- 10
      while(inherits(x, "error") && k >= 0){
        x <- tryCatch(thr(round(mid_value,k)), error = function(e) e)
        k <- k-1
      }
      if(inherits(x, "error")) stop("MDACT integration failed after bounded retries: ", conditionMessage(x))
      if (x < 0) {
        lower_bound <- mid_index
      } else {
        upper_bound <- mid_index
      }
    }
    t_s <- Ts[sorted_index[lower_bound]]
    
    if(thr(t_s) < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }
  }
  if(control.method == "FWER"){
    thr <- function(t) {
      c * (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) + 
                              pi.01.est/(pi.01.est+pi.11.est)*F01(t) + 
                              (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) - 
                                 pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t)))* sim.num - significance_upper
    }
    
    
    sorted_index <- order(Ts) # Get the indices that sort 'Ts'
    lower_bound <- 1
    upper_bound <- length(Ts)
    
    while (upper_bound - lower_bound > 1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]
      
      x <- tryCatch(thr(mid_value), error = function(e) e)
      k <- 10
      while(inherits(x, "error") && k >= 0){
        x <- tryCatch(thr(round(mid_value,k)), error = function(e) e)
        k <- k-1
      }
      if(inherits(x, "error")) stop("MDACT integration failed after bounded retries: ", conditionMessage(x))
      if (x < 0) {
        lower_bound <- mid_index
      } else {
        upper_bound <- mid_index
      }
    }
    t_s <- Ts[sorted_index[lower_bound]]
    
    if(thr(t_s) < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }
  }
  
  if(control.method == "size"){
    thr <- function(t) {
      (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) + 
                          pi.01.est/(pi.01.est+pi.11.est)*F01(t) + 
                          (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) - 
                             pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t))) - significance_upper
    }
    x <- tryCatch(uniroot(thr, c(min(c(p.M, p.Y)), 1))$root, error = function(e) e)
    if(!inherits(x, "error")){
      t_s <- x
    }else{
      F_emp <- function(t) {
        F_emp <- (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) + 
                                     pi.01.est/(pi.01.est+pi.11.est)*F01(t) + 
                                     (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) - 
                                        pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t)))
        return(F_emp)
      }
      
      sorted_index <- order(Ts) # Get the indices that sort 'Ts'
      lower_bound <- 1
      upper_bound <- length(Ts)
      
      while (upper_bound - lower_bound > 1) {
        mid_index <- floor((lower_bound + upper_bound) / 2)
        mid_value <- Ts[sorted_index[mid_index]]
        
        if (F_emp(round(mid_value,6)) < significance_upper) {
          lower_bound <- mid_index
        } else {
          upper_bound <- mid_index
        }
      }
      t_s <- Ts[sorted_index[lower_bound]]
    }
    
    if(F_emp(round(t_s,6)) < significance_upper){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }
  }
}
