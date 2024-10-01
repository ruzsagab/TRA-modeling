setwd("/home/ruzsagab/Documents/Research/Behavioral_addictions")
plotDir <- "Graphs/"
par(mai=c(1.2,1.2,0.4,0.4))
if (!dir.exists(plotDir)) {dir.create(plotDir)}


### Data generating process for the observed consumption values as in the Theory of Rational Addiction ###
#   TRA consumption path: c_t = beta0 + beta1 * s_{t-1} + beta2 * s^2_{t-1} + sigma_u * u_t
#   Observation equation: o_t = c_t + sigma_w * w_t

TRA_betas <- function(cons_quadr, cons_extr, arg_extr) {
  beta2 <- cons_quadr
  beta1 <- (-2) * beta2 * arg_extr
  beta0 <- beta2 * arg_extr^2 + cons_extr
  return(c(beta0=beta0, beta1=beta1, beta2=beta2))
}

quadr.Fn <- function(x, beta0, beta1, beta2) {
  y <- beta0 + beta1 * x + beta2 * x^2
  return(y)
}

quadrTangent <- function(x, beta0, beta1, beta2) {
  y <- beta0 + beta1 * x + beta2 * x^2
  slope <- beta1 + beta2 * 2*x
  intercept <- y - x * slope
  return(c(intercept=intercept, slope=slope))
}

halfLife <- function(delta) {
  return(log(0.5)/log(1-delta))
}

TRA_diagnostics <- function(beta0, beta1, beta2, delta, cons_lim) {
  # Calculate the steady state: [s*, c*] such that s* = c* + (1-delta) * s* = beta0 + (1+beta1-delta) * s* + beta2 * s*^2
  if ((beta1-delta)^2 - 4*beta0*beta2 >= 0) {
    s_ss1 <- (-(beta1-delta) + sqrt((beta1-delta)^2 - 4*beta0*beta2)) / (2*beta2)
    s_ss2 <- (-(beta1-delta) - sqrt((beta1-delta)^2 - 4*beta0*beta2)) / (2*beta2)
    c_ss1 <- quadr.Fn(s_ss1, beta0, beta1, beta2)
    c_ss2 <- quadr.Fn(s_ss2, beta0, beta1, beta2)
    ss1 <- list(s=s_ss1, c=c_ss1, valid=(s_ss1>=0 && c_ss1>=0), stable=(beta1+2*beta2*s_ss1<delta), compl=(beta1+2*beta2*s_ss1>0))
    ss2 <- list(s=s_ss2, c=c_ss2, valid=(s_ss2>=0 && c_ss2>=0), stable=(beta1+2*beta2*s_ss2<delta), compl=(beta1+2*beta2*s_ss2>0))
    
    # Identify the stable [sss] and unstable [ssu] steady states 
    sss <- ssu <- NULL
    
    if (ss1$valid) {
      if (ss1$stable) {sss <- ss1[c("s","c","compl")]
      } else {ssu <- ss1[c("s","c","compl")]}
    }
    
    if (ss2$valid) {
      if (ss2$stable) {sss <- ss2[c("s","c","compl")]
      } else {ssu <- ss2[c("s","c","compl")]}
    } 
  } else {sss <- ssu <- ss1 <- ss2 <- NULL}
  
  # Identify the limiting steady state [ssl]
  s_lim <- cons_lim / delta
  c_lim <- cons_lim
  lim <- list(s=s_lim, c=c_lim, valid=(quadr.Fn(s_lim, beta0, beta1, beta2)>c_lim), stable=TRUE, compl=FALSE)
  if (lim$valid) {ssl <- lim} else {ssl <- NULL}
  
  # Calculate the point of abstinence: [s_a] such that for all s < s_a we have c(s) = beta0 + beta1 * s + beta2 * s^2 < 0
  if (beta2 <= 0 && beta1^2 - 4*beta0*beta2 >= 0) {
    s_a <- (-beta1 + sqrt(beta1^2 - 4*beta0*beta2)) / (2*beta2)
    root <- list(s=s_a, c=0, valid=(s_a>=0), compl=TRUE)
    
    if (root$valid) {abst <- root[c("s","c","compl")]
    } else {abst <- NULL}
  } else {abst <- root <- NULL}
  
  return(list(ss1=ss1, ss2=ss2, sss=sss, ssu=ssu, lim=lim, ssl=ssl, root=root, abst=abst))
}

### Function to plot the TRS consumption path
TRA_path <- function(cons_quadr, cons_extr, arg_extr, delta, simLab=NA, plotPDF=FALSE) {
  betas.Vect <- TRA_betas(cons_quadr, cons_extr, arg_extr)
  beta0 <- as.numeric(betas.Vect["beta0"])
  beta1 <- as.numeric(betas.Vect["beta1"])
  beta2 <- as.numeric(betas.Vect["beta2"])
  
  # Plot the consumption function
  diagn.List <- TRA_diagnostics(beta0, beta1, beta2, delta, cons_lim)
  sss <- diagn.List$sss
  ssu <- diagn.List$ssu
  ssl <- diagn.List$ssl
  abst <- diagn.List$abst
  
  s_dist <- max(ifelse(is.null(sss), 0, abs(sss$s - arg_extr)),
                ifelse(is.null(ssu), 0, abs(ssu$s - arg_extr)),
                ifelse(is.null(ssl), 0, abs(ssl$s - arg_extr)))
  z_path.Vect <- seq(from=0, to=arg_extr+s_dist, length.out=1000)
  c_path.Vect <- mapply(FUN=function(x) {min(cons_lim, max(0, quadr.Fn(x, beta0, beta1, beta2)))}, z_path.Vect)
  mainText <- paste0(ifelse(is.na(simLab), "", paste0("[",simLab,"] ")), "TRA consumption path")
  
  if (plotPDF) {pdf(file=paste0(plotDir,simLab,".pdf"), width=9, height=6, family="Helvetica")}
  plot(x=z_path.Vect, y=c_path.Vect, col="black", type="l", lwd=2, main=mainText, 
       xlim=c(0, max(z_path.Vect)), ylim=c(min(c_path.Vect), max(c_path.Vect)),
       xlab="consumption capital (S)", ylab="consumption (c)")
  abline(a=0, b=delta, col="red", lty=2)
  if (!is.null(sss)) {
    abline(v=sss$s, col="grey25", lty=2)
    text(x=sss$s, y=min(c_path.Vect), labels=paste0("sss=", round(sss$s,2)))
  }
  if (!is.null(ssu)) {
    abline(v=ssu$s, col="grey25", lty=2)
    text(x=ssu$s, y=min(c_path.Vect), labels=paste0("ssu=", round(ssu$s,2)))
  }
  if (!is.null(ssl)) {
    abline(v=ssl$s, col="grey25", lty=2)
    text(x=ssl$s, y=min(c_path.Vect), labels=paste0("ssl=", round(ssl$s,2)))
  }
  if (!is.null(abst$s)) {
    abline(v=abst$s, col="grey25", lty=2)
    text(x=abst$s, y=max(c_path.Vect), labels=paste0("abst=", round(abst$s,2)))
  }
  text(x=arg_extr, y=ifelse(beta2<0, 0.5*(min(c_path.Vect)+max(c_path.Vect)), max(c_path.Vect)),
       labels=paste0("\n\n", "beta2 = ", round(beta2, 4), ";   ", "beta1 = ", round(beta1, 4), "\n",
       "beta0 = ", round(beta0, 4), ";   ", "delta = ", round(delta, 4)))
  if (plotPDF) {dev.off()}
}

### Function to perform numeric simulations for the TRA model
TRA_sim <- function(cons_quadr, cons_extr, arg_extr, cons_lim, delta, s0, sigma_u, sigma_w, nObs,
                    seed=NULL, simLab=NA, plotPDF=FALSE) {
  betas.Vect <- TRA_betas(cons_quadr, cons_extr, arg_extr)
  beta0 <- as.numeric(betas.Vect["beta0"])
  beta1 <- as.numeric(betas.Vect["beta1"])
  beta2 <- as.numeric(betas.Vect["beta2"])
  
  # Deterministic consumption function
  cons.Fn <- function(s) {min(cons_lim, max(0, quadr.Fn(s, beta0, beta1, beta2)))}
  
  # Simulate the dynamics of (c, s)
  if (is.null(seed)) {set.seed(NULL); seed <- .Random.seed[626]}
  set.seed(seed)
  
  o.Vect <- c.Vect <- s.Vect <- c_det.Vect <- s_det.Vect <- numeric(nObs)
  u.Vect <- rnorm(nObs) # disturbance of consumption
  w.Vect <- rnorm(nObs) # observation noise
  
  c_det.Vect[1] <- cons.Fn(s0)
  s_det.Vect[1] <- c_det.Vect[1] + (1-delta) * s0
  c.Vect[1] <- max(0, cons.Fn(s0) + sigma_u * u.Vect[1])
  s.Vect[1] <- c.Vect[1] + (1-delta) * s0
  o.Vect[1] <- max(0, c.Vect[1] + sigma_w * w.Vect[1])
  
  for (i in 2:nObs) {
    c_det.Vect[i] <- cons.Fn(s_det.Vect[i-1])
    s_det.Vect[i] <- c_det.Vect[i] + (1-delta) * s_det.Vect[i-1]
    c.Vect[i] <- max(0, cons.Fn(s.Vect[i-1]) + sigma_u * u.Vect[i])
    s.Vect[i] <- c.Vect[i] + (1-delta) * s.Vect[i-1]
    o.Vect[i] <- max(0, c.Vect[i] + sigma_w * w.Vect[i])
  }
  
  # Define variables [z], [z_det] for the 1 period lagged consumption capital
  z.Vect <- c(s0, s.Vect[-nObs]) # z_t = s_{t-1}
  z_det.Vect <- c(s0, s_det.Vect[-nObs]) # z_det_t = s_det_{t-1} 
  
  # Determine the convergence and divergence points on the consumption path
  diagn.List <- TRA_diagnostics(beta0, beta1, beta2, delta, cons_lim)
  sss <- diagn.List$sss
  ssu <- diagn.List$ssu
  ssl <- diagn.List$ssl
  abst <- diagn.List$abst
  
  s_dist <- sigma_u/delta + max(ifelse(is.null(sss), 0, abs(sss$s - arg_extr)),
                                ifelse(is.null(ssu), 0, abs(ssu$s - arg_extr)),
                                ifelse(is.null(ssl), 0, abs(ssl$s - arg_extr)))
  z_path.Vect <- seq(from=0, to=arg_extr+s_dist, length.out=1000)
  c_path.Vect <- mapply(FUN=cons.Fn, z_path.Vect)
  
  # Determine the time index when (if ever) consumption has converged to a
  # stable [sss] or limiting [ssl] steady state or to the point of abstinence [abst]
  convergence <- "none"
  i_conv <- NA # Time index when consumption has converged
  nObs_nconv <- nObs # Number of observations before consumption has converged
  
  if (!is.null(sss)) {
    if (!is.null(ssu)) {sssReachable <- (s.Vect[nObs] - ssu$s) * (sss$s - ssu$s) >= 0
    } else {sssReachable <- TRUE}
    
    sssReached <- any(abs(s.Vect - sss$s) < sigma_u)
    if (sssReached) {
      convergence <- "sss"
      i_conv <- which(abs(s.Vect - sss$s) < sigma_u)[1]
      nObs_nconv <- i_conv-1
    }
  } else {sssReachable <- FALSE; sssReached <- FALSE}
  
  if (!is.null(ssl)) {
    if (!is.null(ssu)) {sslReachable <- (s.Vect[nObs] - ssu$s) * (ssl$s - ssu$s) >= 0
    } else {sslReachable <- TRUE}
    
    sslReached <- any(abs(s.Vect - ssl$s) < sigma_u)
    if (sslReached) {
      convergence <- "ssl"
      i_conv <- which(abs(s.Vect - ssl$s) < sigma_u)[1]
      nObs_nconv <- i_conv-1
    }
  } else {sslReachable <- FALSE; sslReached <- FALSE}
  
  if (!is.null(abst)) {
    if (!is.null(ssu)) {abstReachable <- (s.Vect[nObs] - ssu$s) * (abst$s - ssu$s) >= 0  
    } else {abstReachable <- (is.null(sss) && is.null(ssl))}
    
    abstReached <- any(s.Vect < abst$s)
    if (abstReached) {
      convergence <- "abst"
      i_conv <- which(s.Vect < abst$s)[1]
      nObs_nconv <- i_conv-1
    }
  } else {abstReachable <- FALSE; abstReached <- FALSE}
  
  convReached <- sssReached || sslReached || abstReached
  
  s_max <- max(s.Vect) + sigma_u/delta
  s_min <- min(s.Vect) - sigma_u/delta
  c_max <- max(c.Vect) + sigma_u
  c_min <- min(c.Vect) - sigma_u
  o_max <- max(max(o.Vect), max(c.Vect)) + sigma_u
  o_min <- min(min(o.Vect), min(c.Vect)) - sigma_u
  s_fin <- s_det.Vect[nObs]
  c_fin <- c_det.Vect[nObs]
  
  ssuPlot <- FALSE
  if (!is.null(ssu)) {ssuPlot <- ((ssu$s - s_fin) * (s0 - s_fin) >= 0)}

  s_extr <- ifelse(ssuPlot, ssu$s, ifelse(beta2!=0, (-beta1)/(2*beta2), 0))
  s_minn <- min(s_min, s_extr - sigma_u/delta)
  s_maxx <- max(s_max, s_extr + sigma_u/delta)
  c_segm.Vect <- mapply(FUN=cons.Fn, seq(from=s_minn, to=s_maxx, length.out=1000))
  c_minn <- min(c_min, min(c_segm.Vect) - sigma_u)
  c_maxx <- max(c_max, max(c_segm.Vect) + sigma_u)
  
  
    # Plot of (S, c)
  if (plotPDF) {pdf(file=paste0(plotDir,simLab,"_a.pdf"), width=9, height=6, family="Helvetica")}
  mainText <- paste0(ifelse(is.na(simLab), "", paste0("[",simLab,"] ")),
                     "Simulated consumption path", ", S0 = ", round(s0,2))
  yleg <- ifelse(c_fin > 0.5*(c_minn+c_maxx), 0.5*(c_minn+c_maxx), c_maxx)
  pch_dir <- ifelse(z_det.Vect[nObs]<s0, 60, 62)
  plot(x=z_path.Vect, y=c_path.Vect, xlim=c(s_minn, s_maxx), ylim=c(c_minn, c_maxx), col="orange", type="l", lty=1, lwd=2,
       xlab="consumption capital (S)", ylab="consumption (c)", main=mainText)
  points(x=z_det.Vect[1:nObs], y=c_det.Vect[1:nObs], col="orange", pch=pch_dir)
  points(x=z.Vect[1:min(nObs_nconv+1,nObs)], y=c.Vect[1:min(nObs_nconv+1,nObs)], col="blue4")
  if (convReached) {points(x=z.Vect[(i_conv+1):nObs], y=c.Vect[(i_conv+1):nObs], col="green4")}
  abline(a=0, b=delta, col="red", lty=2)
  legend(x=0.5*(s_minn+s_maxx), y=yleg,
         legend=c("stochastic (before conv.)","stochastic (converged)","deterministic path","steady state condition"),
         col=c("blue4","green4","orange","red"), lwd=c(NA,NA,2,1), pch=c(1,1,pch_dir,NA), lty=c(NA,NA,1,2),
         ncol=1, x.intersp=0.1, y.intersp=1.0, bty="n")
  if (plotPDF) {dev.off()}
  
  # Plot of c(t)
  if (plotPDF) {pdf(file=paste0(plotDir,simLab,"_c.pdf"), width=9, height=6, family="Helvetica")}
  mainText <- paste0(ifelse(is.na(simLab), "", paste0("[",simLab,"] ")),
                     "Simulated data for daily consumption", ", S0 = ", round(s0,2))
  yleg <- ifelse(c_fin > 0.5*(o_min+o_max), 0.5*(o_min+o_max), o_max)
  plot(x=1:nObs, y=o.Vect, col="lightblue", xlim=c(0, nObs), ylim=c(o_min, o_max), type="p",
       xlab="t", ylab="consumption", main=mainText)
  lines(x=1:nObs, y=o.Vect, col="lightblue")
  lines(x=1:nObs, y=c_det.Vect, col="orange", lwd=2)
  lines(x=1:min(nObs_nconv+1, nObs), y=c.Vect[1:min(nObs_nconv+1, nObs)], col="blue4")
  if (convReached) {lines(x=i_conv:nObs, y=c.Vect[i_conv:nObs], col="green4", lwd=1); abline(v=i_conv, col="grey25", lty=2)}
  if (sssReachable || sslReachable) {abline(h=ifelse(sssReachable, sss$c, ssl$c), lty=2, col="red")}
  legend(x=0.5*nObs, y=yleg,
         legend=c("observed","latent (before conv.)","latent (converged)","deterministic"),
         col=c("lightblue","blue4","green4","orange"), lty=c(1,1,1,1), lwd=c(1,1,1,2), pch=c(1,NA,NA,NA),
         ncol=1, x.intersp=0.1, y.intersp=1.0, bty="n")
  if (plotPDF) {dev.off()}
  
  # Plot of s(t)
  if (plotPDF) {pdf(file=paste0(plotDir,simLab,"_s.pdf"), width=9, height=6, family="Helvetica")}
  mainText <- paste0(ifelse(is.na(simLab), "", paste0("[",simLab,"] ")),
                     "Simulated data for consumption capital", ", S0 = ", round(s0,2))
  yleg <- ifelse(s_fin > 0.5*(s_min+s_max), 0.5*(s_min+s_max), s_max)
  plot(x=1:nObs, y=s_det.Vect[1:nObs], col="orange", xlim=c(0, nObs), ylim=c(s_min, s_max),
       type="l", lwd=2, xlab="t", ylab="consumption capital", main=mainText)
  lines(x=1:min(nObs_nconv+1, nObs), y=s.Vect[1:min(nObs_nconv+1, nObs)], col="blue4")
  if (convReached) {lines(x=i_conv:nObs, y=s.Vect[i_conv:nObs], col="green4"); abline(v=i_conv, col="grey25", lty=2)}
  if (sssReachable || sslReachable) {abline(h=ifelse(sssReachable, sss$s, ssl$s), lty=2, col="red")} 
  legend(x=0.5*nObs, y=yleg,
         legend=c("latent (before conv.)","latent (converged)","deterministic"),
         col=c("blue4","green4","orange"), lty=c(1,1,1), lwd=c(1,1,2), ncol=1, x.intersp=0.1, y.intersp=1.0, bty="n")
  if (plotPDF) {dev.off()}
  
  sim.List <- list(params=list(cons_quadr=cons_quadr, cons_extr=cons_extr, arg_extr=arg_extr,
                               delta=delta, s0=s0, sigma_u=sigma_u, sigma_w=sigma_w),
                   nObs=nObs, seed=seed, diagnostic=diagn.List, convergence=convergence, i_conv=i_conv,
                   noise=list(u=u.Vect, w=w.Vect),
                   data.full=list(c=c.Vect, s=s.Vect, o=o.Vect, c_det=c_det.Vect, s_det=s_det.Vect),
                   data.nconv=list(c=c.Vect[1:nObs_nconv], s=s.Vect[1:nObs_nconv], o=o.Vect[1:nObs_nconv]),
                   data.conv=NULL)
  
  if (convReached) {sim.List$data.conv <- list(c=c.Vect[(i_conv-1):nObs], s=s.Vect[(i_conv-1):nObs], o=o.Vect[(i_conv-1):nObs])} 
  return(sim.List)
}



### Input parameters for determining the TRA consumption function: c(s) = beta0 + beta1 * s + beta2 * s^2
if (!plotPDF) {
  cons_quadr <- 0.005 # quadratic coefficient [beta2] for the consumption fn.
  cons_extr <- 0.2 # extremum of the consumption fn. -> max or min (depending on cons_quadr being neg. or pos.)
  arg_extr <- 15 # maximizer / minimizer argument value of the consumption fn. -> arg_max or arg_min
  cons_lim <- 2 # upper limit on consumption
  delta <- 0.05 # rate of exponential decay
  cat("Half-life:", round(halfLife(delta), 2), "days") # half-life of consumption capital (in days)
  TRA_path(cons_quadr, cons_extr, arg_extr, delta)

  # Initial condition
  s0 <- 5 # initial value of consumption capital
  
  # Parameters for the disturbance / noise components
  sigma_u <- 0.05 # dispersion parameter for the disturbance of consumption
  sigma_w <- 0.075 # dispersion parameter for observation noise
  
  # Number of observations
  nObs <- 150
}

### Simulate data for various parameter settings
TRA_params.Matrix <- matrix(nrow=11, ncol=15)
colnames(TRA_params.Matrix) <- c("beta0","beta1","cons_quadr","cons_extr","arg_extr","cons_lim","delta",
                                 "s0","sigma_u","sigma_w","nObs","sss_s","ssu_s","ssl_s","abst_s")
rownames(TRA_params.Matrix) <- c("par.1",paste0("sim.1.",1:3), "par.2",paste0("sim.2.",1:2), "par.3",paste0("sim.3.",1:3))
TRA_params.Matrix["par.1", 3:11] <- c(-0.008, 1.25, 18, 2, 0.075, rep(NA, 4))
TRA_params.Matrix["sim.1.1", 3:11] <- c(-0.008, 1.25, 18, 2, 0.075,  9, 0.05, 0.075, 150)
TRA_params.Matrix["sim.1.2", 3:11] <- c(-0.008, 1.25, 18, 2, 0.075, 12, 0.05, 0.075, 150)
TRA_params.Matrix["sim.1.3", 3:11] <- c(-0.008, 1.25, 18, 2, 0.075, 20, 0.05, 0.075, 150)
TRA_params.Matrix["par.2", 3:11] <- c(-0.006, 0.90, 20, 2, 0.040, rep(NA, 4))
TRA_params.Matrix["sim.2.1", 3:11] <- c(-0.006, 0.90, 20, 2, 0.040, 13, 0.05, 0.075, 150)
TRA_params.Matrix["sim.2.2", 3:11] <- c(-0.006, 0.90, 20, 2, 0.040, 28, 0.05, 0.075, 150)
TRA_params.Matrix["par.3", 3:11] <- c(+0.005, 0.20, 15, 2, 0.050,  rep(NA, 4))
TRA_params.Matrix["sim.3.1", 3:11] <- c(+0.005, 0.20, 15, 2, 0.050,  0, 0.05, 0.075, 150)
TRA_params.Matrix["sim.3.2", 3:11] <- c(+0.005, 0.20, 15, 2, 0.050, 30, 0.05, 0.075, 150)
TRA_params.Matrix["sim.3.3", 3:11] <- c(+0.005, 0.20, 15, 2, 0.050, 33, 0.05, 0.075, 150)

plotPDF <- TRUE
for (r in 1:nrow(TRA_params.Matrix)) {
  for (parName in c("cons_quadr", "cons_extr", "arg_extr", "cons_lim", "delta", "s0", "sigma_u", "sigma_w", "nObs")) {
    assign(parName, as.numeric(TRA_params.Matrix[r, parName]))
  }
  
  betas.Vect <- TRA_betas(cons_quadr, cons_extr, arg_extr)
  beta0 <- as.numeric(betas.Vect["beta0"])
  beta1 <- as.numeric(betas.Vect["beta1"])
  beta2 <- as.numeric(betas.Vect["beta2"])
  
  diagn.List <- TRA_diagnostics(beta0, beta1, beta2, delta, cons_lim)
  sss <- diagn.List$sss
  ssu <- diagn.List$ssu
  ssl <- diagn.List$ssl
  abst <- diagn.List$abst

  TRA_params.Matrix[r, c("beta0","beta1")] <- c(beta0, beta1)
  for (x in c("sss","ssu","ssl","abst")) {
    if (!is.null(diagn.List[[x]])) {TRA_params.Matrix[r, paste0(x,"_s")] <- diagn.List[[x]]$s}
  }
  
  if (substr(rownames(TRA_params.Matrix)[r], start=1, stop=4) == "par.") {
    TRA_path(cons_quadr, cons_extr, arg_extr, delta,
             simLab=paste0("sim.", substr(rownames(TRA_params.Matrix)[r], start=5, stop=100), ".0"), plotPDF=plotPDF)
  } else {
    assign(paste0(rownames(TRA_params.Matrix)[r], ".List"),
           TRA_sim(cons_quadr, cons_extr, arg_extr, cons_lim, delta, s0, sigma_u, sigma_w, nObs,
                   seed=666+r, simLab=rownames(TRA_params.Matrix)[r], plotPDF=plotPDF))
  }
}  


# sim.List <- TRA_sim(cons_quadr=cons_quadr, cons_extr=cons_extr, arg_extr=arg_extr, delta=delta,
#                     s0=s0, sigma_u=sigma_u, sigma_w=sigma_w, nObs=nObs, seed=666)
# names(sim.List)
# sim.List$params
# sim.List$nObs
# sim.List$seed
# sim.List$diagnostic
# sim.List$convergence
# sim.List$i_conv
# sim.List$noise
# sim.List$data.full
# sim.List$data.nconv
# sim.List$data.conv
