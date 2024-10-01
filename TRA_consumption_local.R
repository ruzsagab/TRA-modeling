setwd("/home/ruzsagab/Documents/Research/Behavioral_addictions")
plotDir <- "Graphs/"
par(mai=c(1.2,1.2,0.4,0.4))
if (!dir.exists(plotDir)) {dir.create(plotDir)}


### Data generating process for the observed consumption values as in the Theory of Rational Addiction ###
#   TRA consumption path: c_t = beta0 + beta1 * s_{t-1} + beta2 * s^2_{t-1} + sigma_u * u_t
#   Observation equation: o_t = c_t + sigma_w * w_t

TRA_betas <- function(cons_quadr, sss_c, sss_slope, delta) {
  sss_s <- sss_c / delta
  beta2 <- cons_quadr
  beta1 <- sss_slope - (2 * cons_quadr * sss_s)
  beta0 <- sss_c - (beta1 * sss_s + beta2 * sss_s^2)
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

TRA_diagnostics <- function(cons_lim, cons_quadr, sss_c, sss_slope, delta) {
  if (sss_slope >= delta) {cat("Steady state stability condition (sss_slope < delta) is violated!"); return(NULL)}
  if (sss_c >= cons_lim) {cat("Steady state feasiblity condition (sss_c < cons_lim) is violated!"); return(NULL)}
  
  betas.Vect <- TRA_betas(cons_quadr, sss_c, sss_slope, delta)
  beta0 <- as.numeric(betas.Vect["beta0"])
  beta1 <- as.numeric(betas.Vect["beta1"])
  beta2 <- as.numeric(betas.Vect["beta2"])

  arg_extr <- ifelse(beta2!=0, ((-beta1) / (2*beta2)), NA)
  sss_s <- sss_c / delta
  ssu_s <- ifelse(beta2!=0, (-(beta1-delta) + sqrt((beta1-delta)^2 - 4*beta0*beta2)) / (2*beta2), NA)
  
  if (beta2!=0) {
    if (beta1^2 - 4*beta0*beta2 >= 0) {
      floor_root1 <- (-beta1 + sqrt(beta1^2 - 4*beta0*beta2)) / (2*beta2)
      floor_root2 <- (-beta1 - sqrt(beta1^2 - 4*beta0*beta2)) / (2*beta2)  
    } else {floor_root1 <- floor_root2 <- NA}
    
    if (beta1^2 - 4*(beta0-cons_lim)*beta2 >= 0) {
      ceiling_root1 <- (-beta1 + sqrt(beta1^2 - 4*(beta0-cons_lim)*beta2)) / (2*beta2)
      ceiling_root2 <- (-beta1 - sqrt(beta1^2 - 4*(beta0-cons_lim)*beta2)) / (2*beta2)  
    } else {ceiling_root1 <- ceiling_root2 <- NA}
  } else {
    floor_root1 <- floor_root2 <- (-beta0) / beta1
    ceiling_root1 <- ceiling_root2 <- (cons_lim - beta0) / beta1
  }
 
  if (beta2 > 0) {
    llim_s <- max(c(0, ifelse(sss_s>=arg_extr, arg_extr, NA), ifelse(sss_s>=floor_root1, floor_root1, NA), ceiling_root2), na.rm=TRUE)
    rlim_s <- min(c(ssu_s, ifelse(sss_s<=arg_extr, arg_extr, NA), ifelse(sss_s<=floor_root2, floor_root2, NA), ceiling_root1), na.rm=TRUE)
  } else if (beta2 < 0) {
    llim_s <- max(c(0, ssu_s, ifelse(sss_s>=arg_extr, arg_extr, NA), floor_root1, ifelse(sss_s>=ceiling_root2, ceiling_root2, NA)), na.rm=TRUE)
    rlim_s <- min(c(ifelse(sss_s<=arg_extr, arg_extr, NA), floor_root2, ifelse(sss_s<=ceiling_root1, ceiling_root1, NA)), na.rm=TRUE)
  } else {
    llim_s <- max(c(0, ifelse(sss_s>=floor_root1, floor_root1, NA), ifelse(sss_s>=ceiling_root1, ceiling_root1, NA)), na.rm=TRUE) 
    rlim_s <- min(c(ifelse(sss_s<=floor_root1, floor_root1, NA), ifelse(sss_s<=ceiling_root1, ceiling_root1, NA)), na.rm=TRUE)
  }
  
  return(list(sss_s=sss_s, ssu_s=ssu_s, arg_extr=arg_extr, floor_root1=floor_root1, floor_root2=floor_root2,
              ceiling_root1=ceiling_root1, ceiling_root2=ceiling_root2, llim_s=llim_s, rlim_s=rlim_s))
}

### Function to plot the TRS consumption path
TRA_path <- function(cons_lim, cons_quadr, sss_c, sss_slope, delta, simLab=NA, plotPDF=FALSE) {
  betas.Vect <- TRA_betas(cons_quadr, sss_c, sss_slope, delta)
  beta0 <- as.numeric(betas.Vect["beta0"])
  beta1 <- as.numeric(betas.Vect["beta1"])
  beta2 <- as.numeric(betas.Vect["beta2"])
  
  # Plot the consumption function
  diagn.List <- TRA_diagnostics(cons_lim, cons_quadr, sss_c, sss_slope, delta)
  arg_extr <- diagn.List$arg_extr
  sss_s <- diagn.List$sss_s
  llim_s <- diagn.List$llim_s
  rlim_s <- diagn.List$rlim_s
  llim_c <- quadr.Fn(llim_s, beta0, beta1, beta2)
  rlim_c <- quadr.Fn(rlim_s, beta0, beta1, beta2)
  dlim_s <- rlim_s - llim_s
  
  z_path.Vect <- seq(from=llim_s+0.1*dlim_s, to=rlim_s-0.1*dlim_s, length.out=160)
  z_lim_path.Vect <- seq(from=llim_s, to=rlim_s, length.out=200)
  c_path.Vect <- mapply(FUN=function(x) {min(cons_lim, max(0, quadr.Fn(x, beta0, beta1, beta2)))}, z_path.Vect)
  c_lim_path.Vect <- mapply(FUN=function(x) {min(cons_lim, max(0, quadr.Fn(x, beta0, beta1, beta2)))}, z_lim_path.Vect)
  mainText <- paste0(ifelse(is.na(simLab), "", paste0("[",simLab,"] ")), "TRA local consumption path")
  
  if (plotPDF) {pdf(file=paste0(plotDir,simLab,".pdf"), width=9, height=6, family="Helvetica")}
  plot(x=z_lim_path.Vect, y=c_lim_path.Vect, col="grey25", type="l",
       xlim=c(llim_s, rlim_s), ylim=c(min(llim_c,rlim_c), max(llim_c,rlim_c)),
       xlab="consumption capital (S)", ylab="consumption (c)", main=mainText)
  lines(x=z_path.Vect, y=c_path.Vect, col="black", lwd=3)
  abline(a=0, b=delta, col="red", lty=2)
  abline(v=c(llim_s, rlim_s), col="grey25", lty=2)
  points(x=sss_s, y=sss_c, col="red", cex=2.0, pch=3)
  text(x=0.5*(llim_s+rlim_s), y=ifelse(beta2<0, 0.5*(llim_c+rlim_c), max(llim_c,rlim_c)),
       labels=paste0("\n\n", "cons_quadr = ", round(cons_quadr, 4), ";   ", "sss_c = ", round(sss_c, 4), "\n",
                     "sss_slope = ", round(sss_slope, 4), ";   ", "delta = ", round(delta, 4)))
  if (plotPDF) {dev.off()}
}

### Function to perform numeric simulations for the TRA model
TRA_sim <- function(cons_lim, cons_quadr, sss_c, sss_slope, delta, s0=NA, sigma_u, sigma_w, nObs,
                    pathDir=ifelse(is.na(s0), NA, ifelse(s0<sss_c/delta, "right", "left")), seed=NULL, simLab=NA, plotPDF=FALSE) {
  betas.Vect <- TRA_betas(cons_quadr, sss_c, sss_slope, delta)
  beta0 <- as.numeric(betas.Vect["beta0"])
  beta1 <- as.numeric(betas.Vect["beta1"])
  beta2 <- as.numeric(betas.Vect["beta2"])

  # Diagnostics for the TRA consumption path
  diagn.List <- TRA_diagnostics(cons_lim, cons_quadr, sss_c, sss_slope, delta)
  arg_extr <- diagn.List$arg_extr
  sss_s <- diagn.List$sss_s
  llim_s <- diagn.List$llim_s
  rlim_s <- diagn.List$rlim_s
  llim_c <- quadr.Fn(llim_s, beta0, beta1, beta2)
  rlim_c <- quadr.Fn(rlim_s, beta0, beta1, beta2)
  dlim_s <- rlim_s - llim_s
  if (is.na(s0)) {
    if (is.na(pathDir)) {pathDir <- ifelse(rbinom(1,1,0.5), "right", "left")
    } else if (!pathDir %in% c("right", "left")) {pathDir <- ifelse(rbinom(1,1,0.5), "right", "left")}
    
    if (pathDir == "right") {s0 <- llim_s + 0.1*dlim_s
    } else {s0 <- rlim_s - 0.1*dlim_s}
  } else {s0 <- min(rlim_s - 0.1*dlim_s, max(llim_s + 0.1*dlim_s, s0))}
  
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
  z_path.Vect <- seq(from=llim_s, to=rlim_s, length.out=200)
  c_path.Vect <- mapply(FUN=cons.Fn, z_path.Vect)
  
  # Determine the time index when (if ever) consumption has converged to the stable steady state [sss]
  i_conv <- NA # Time index when consumption has converged
  nObs_nconv <- nObs # Number of observations before consumption has converged
  
  sssReached <- any(abs(s.Vect - sss_s) < sigma_u)
  if (sssReached) {
    i_conv <- which(abs(s.Vect - sss_s) < sigma_u)[1]
    nObs_nconv <- i_conv-1
  }
  
  s_max <- max(s.Vect) + sigma_u/delta
  s_min <- min(s.Vect) - sigma_u/delta
  c_max <- max(c.Vect) + sigma_u
  c_min <- min(c.Vect) - sigma_u
  o_max <- max(max(o.Vect), max(c.Vect)) + sigma_u
  o_min <- min(min(o.Vect), min(c.Vect)) - sigma_u
  
  # Plot of (S, c)
  if (plotPDF) {pdf(file=paste0(plotDir,simLab,"_a.pdf"), width=9, height=6, family="Helvetica")}
  mainText <- paste0(ifelse(is.na(simLab), "", paste0("[",simLab,"] ")),
                     "Simulated consumption path", ", S0 = ", round(s0,2))
  yleg <- ifelse(sss_c > 0.5*(c_min+c_max), 0.5*(c_min+c_max), c_max)
  pch_dir <- ifelse(z_det.Vect[nObs]<s0, 60, 62)
  plot(x=z_path.Vect, y=c_path.Vect, xlim=c(s_min, s_max), ylim=c(c_min, c_max), col="orange", type="l", lty=1, lwd=2,
       xlab="consumption capital (S)", ylab="consumption (c)", main=mainText)
  points(x=z_det.Vect[1:nObs], y=c_det.Vect[1:nObs], col="orange", pch=pch_dir)
  points(x=z.Vect[1:min(nObs_nconv+1,nObs)], y=c.Vect[1:min(nObs_nconv+1,nObs)], col="blue4")
  if (sssReached) {points(x=z.Vect[(i_conv+1):nObs], y=c.Vect[(i_conv+1):nObs], col="green4")}
  abline(a=0, b=delta, col="red", lty=2)
  legend(x=0.5*(s_min+s_max), y=yleg,
         legend=c("stochastic (before conv.)","stochastic (converged)","deterministic path","steady state condition"),
         col=c("blue4","green4","orange","red"), lwd=c(NA,NA,2,1), pch=c(1,1,pch_dir,NA), lty=c(NA,NA,1,2),
         ncol=1, x.intersp=0.1, y.intersp=1.0, bty="n")
  if (plotPDF) {dev.off()}
  
  # Plot of c(t)
  if (plotPDF) {pdf(file=paste0(plotDir,simLab,"_c.pdf"), width=9, height=6, family="Helvetica")}
  mainText <- paste0(ifelse(is.na(simLab), "", paste0("[",simLab,"] ")),
                     "Simulated data for daily consumption", ", S0 = ", round(s0,2))
  yleg <- ifelse(sss_c > 0.5*(o_min+o_max), 0.5*(o_min+o_max), o_max)
  plot(x=1:nObs, y=o.Vect, col="lightblue", xlim=c(0, nObs), ylim=c(o_min, o_max), type="p",
       xlab="t", ylab="consumption", main=mainText)
  lines(x=1:nObs, y=o.Vect, col="lightblue")
  lines(x=1:nObs, y=c_det.Vect, col="orange", lwd=2)
  lines(x=1:min(nObs_nconv+1, nObs), y=c.Vect[1:min(nObs_nconv+1, nObs)], col="blue4")
  if (sssReached) {lines(x=i_conv:nObs, y=c.Vect[i_conv:nObs], col="green4", lwd=1); abline(v=i_conv, col="grey25", lty=2)}
  abline(h=sss_c, lty=2, col="red")
  legend(x=0.5*nObs, y=yleg,
         legend=c("observed","latent (before conv.)","latent (converged)","deterministic"),
         col=c("lightblue","blue4","green4","orange"), lty=c(1,1,1,1), lwd=c(1,1,1,2), pch=c(1,NA,NA,NA),
         ncol=1, x.intersp=0.1, y.intersp=1.0, bty="n")
  if (plotPDF) {dev.off()}
  
  # Plot of s(t)
  if (plotPDF) {pdf(file=paste0(plotDir,simLab,"_s.pdf"), width=9, height=6, family="Helvetica")}
  mainText <- paste0(ifelse(is.na(simLab), "", paste0("[",simLab,"] ")),
                     "Simulated data for consumption capital", ", S0 = ", round(s0,2))
  yleg <- ifelse(sss_s > 0.5*(s_min+s_max), 0.5*(s_min+s_max), s_max)
  plot(x=1:nObs, y=s_det.Vect[1:nObs], col="orange", xlim=c(0, nObs), ylim=c(s_min, s_max),
       type="l", lwd=2, xlab="t", ylab="consumption capital", main=mainText)
  lines(x=1:min(nObs_nconv+1, nObs), y=s.Vect[1:min(nObs_nconv+1, nObs)], col="blue4")
  if (sssReached) {lines(x=i_conv:nObs, y=s.Vect[i_conv:nObs], col="green4"); abline(v=i_conv, col="grey25", lty=2)}
  abline(h=sss_s, lty=2, col="red")
  legend(x=0.5*nObs, y=yleg,
         legend=c("latent (before conv.)","latent (converged)","deterministic"),
         col=c("blue4","green4","orange"), lty=c(1,1,1), lwd=c(1,1,2), ncol=1, x.intersp=0.1, y.intersp=1.0, bty="n")
  if (plotPDF) {dev.off()}
  
  sim.List <- list(params=list(cons_lim=cons_lim, cons_quadr=cons_quadr, sss_c=sss_c, sss_slope=sss_slope,
                               delta=delta, s0=s0, sigma_u=sigma_u, sigma_w=sigma_w),
                   nObs=nObs, seed=seed, diagnostic=diagn.List, sssReached=sssReached, i_conv=i_conv,
                   noise=list(u=u.Vect, w=w.Vect),
                   data.full=list(c=c.Vect, s=s.Vect, o=o.Vect, c_det=c_det.Vect, s_det=s_det.Vect),
                   data.nconv=list(c=c.Vect[1:nObs_nconv], s=s.Vect[1:nObs_nconv], o=o.Vect[1:nObs_nconv]),
                   data.conv=NULL)
  
  if (sssReached) {sim.List$data.conv <- list(c=c.Vect[(i_conv-1):nObs], s=s.Vect[(i_conv-1):nObs], o=o.Vect[(i_conv-1):nObs])} 
  return(sim.List)
}

### Simulate data for various parameter settings
TRA_params.Matrix <- matrix(nrow=18, ncol=10)
colnames(TRA_params.Matrix) <- c("cons_quadr","sss_c","sss_slope","cons_lim","delta","pathDir","s0","sigma_u","sigma_w","nObs")
rownames(TRA_params.Matrix) <- c("par.1",paste0("sim.1.",1:2), "par.2",paste0("sim.2.",1:2), "par.3",paste0("sim.3.",1:2),
                                 "par.4",paste0("sim.4.",1:2), "par.5",paste0("sim.5.",1:2), "par.6",paste0("sim.6.",1:2))
TRA_params.Matrix["par.1", 1:5] <- c(-0.001, 0.5, -0.05, 1, 0.075)
TRA_params.Matrix["sim.1.1",] <- c(-0.001, 0.5, -0.05, 1, 0.075, 0, NA, 0.05, 0.075, 250)
TRA_params.Matrix["sim.1.2",] <- c(-0.001, 0.5, -0.05, 1, 0.075, 1, NA, 0.05, 0.075, 250)
TRA_params.Matrix["par.2", 1:5] <- c(0, 0.5, -0.05, 1, 0.075)
TRA_params.Matrix["sim.2.1",] <- c(0, 0.5, -0.05, 1, 0.075, 0, NA, 0.05, 0.075, 250)
TRA_params.Matrix["sim.2.2",] <- c(0, 0.5, -0.05, 1, 0.075, 1, NA, 0.05, 0.075, 250)
TRA_params.Matrix["par.3", 1:5] <- c(+0.001, 0.5, -0.05, 1, 0.075)
TRA_params.Matrix["sim.3.1",] <- c(+0.001, 0.5, -0.05, 1, 0.075, 0, NA, 0.05, 0.075, 250)
TRA_params.Matrix["sim.3.2",] <- c(+0.001, 0.5, -0.05, 1, 0.075, 1, NA, 0.05, 0.075, 250)
TRA_params.Matrix["par.4", 1:5] <- c(-0.001, 0.5, +0.05, 1, 0.075)
TRA_params.Matrix["sim.4.1",] <- c(-0.001, 0.5, +0.05, 1, 0.075, 0, NA, 0.05, 0.075, 250)
TRA_params.Matrix["sim.4.2",] <- c(-0.001, 0.5, +0.05, 1, 0.075, 1, NA, 0.05, 0.075, 250)
TRA_params.Matrix["par.5", 1:5] <- c(0, 0.5, +0.05, 1, 0.075)
TRA_params.Matrix["sim.5.1",] <- c(0, 0.5, +0.05, 1, 0.075, 0, NA, 0.05, 0.075, 250)
TRA_params.Matrix["sim.5.2",] <- c(0, 0.5, +0.05, 1, 0.075, 1, NA, 0.05, 0.075, 250)
TRA_params.Matrix["par.6", 1:5] <- c(+0.001, 0.5, +0.05, 1, 0.075)
TRA_params.Matrix["sim.6.1",] <- c(+0.001, 0.5, +0.05, 1, 0.075, 0, NA, 0.05, 0.075, 250)
TRA_params.Matrix["sim.6.2",] <- c(+0.001, 0.5, +0.05, 1, 0.075, 1, NA, 0.05, 0.075, 250)


plotPDF <- TRUE
for (r in 1:nrow(TRA_params.Matrix)) {
  for (parName in c("cons_quadr","sss_c","sss_slope","cons_lim","delta","sigma_u","sigma_w","nObs")) {
    assign(parName, as.numeric(TRA_params.Matrix[r, parName]))
    pathDir <- ifelse(as.logical(TRA_params.Matrix[r, "pathDir"]), "right", "left")
  }
  
  diagn.List <- TRA_diagnostics(cons_lim, cons_quadr, sss_c, sss_slope, delta)
  llim_s <- diagn.List$llim_s
  rlim_s <- diagn.List$rlim_s
  dlim_s <- rlim_s - llim_s
  s0 <- ifelse(pathDir=="right", llim_s + 0.1*dlim_s, rlim_s - 0.1*dlim_s)
  TRA_params.Matrix[r, "s0"] <- s0

  if (substr(rownames(TRA_params.Matrix)[r], start=1, stop=4) == "par.") {
    TRA_path(cons_lim, cons_quadr, sss_c, sss_slope, delta,
             simLab=paste0("sim.", substr(rownames(TRA_params.Matrix)[r], start=5, stop=100), ".0"), plotPDF=plotPDF)
  } else {
    assign(paste0(rownames(TRA_params.Matrix)[r], ".List"),
           TRA_sim(cons_lim, cons_quadr, sss_c, sss_slope, delta, s0, sigma_u, sigma_w, nObs, pathDir,
                   seed=666+r, simLab=rownames(TRA_params.Matrix)[r], plotPDF=plotPDF))
  }
}  


### Test the TRA_diagnostics() function
performTest <- FALSE
if (performTest) {
  params.Grid <- expand.grid(cons_quadr=c(-0.001,-0.005,0,0.001,0.005),
                               sss_c=c(0.2,0.5,0.8), sss_slope=c(-0.05,0.05), cons_lim=1, delta=0.075)
  
  for (r in 1:nrow(params.Grid)) {
    for (parName in c("cons_lim","cons_quadr","sss_c","sss_slope","delta")) {assign(parName, params.Grid[r, parName])}
    
    diagn.List <- TRA_diagnostics(cons_lim, cons_quadr, sss_c, sss_slope, delta)
    arg_extr <- diagn.List$arg_extr
    sss_s <- diagn.List$sss_s; ssu_s <- diagn.List$ssu_s
    llim_s <- diagn.List$llim_s; rlim_s <- diagn.List$rlim_s
    floor_root1 <- diagn.List$floor_root1; floor_root2 <- diagn.List$floor_root2
    ceiling_root1 <- diagn.List$ceiling_root1; ceiling_root2 <- diagn.List$ceiling_root2
    
    betas.Vect <- TRA_betas(cons_quadr, sss_c, sss_slope, delta)
    beta0 <- as.numeric(betas.Vect["beta0"])
    beta1 <- as.numeric(betas.Vect["beta1"])
    beta2 <- as.numeric(betas.Vect["beta2"])
    s_min <- min(c(0, ssu_s, floor_root1, floor_root2, ceiling_root1, ceiling_root2), na.rm=TRUE) 
    s_max <- max(c(0, ssu_s, floor_root1, floor_root2, ceiling_root1, ceiling_root2), na.rm=TRUE)
    
    x <- seq(from=s_min, to=s_max, length.out=100)
    y <- mapply(FUN=function(x) {quadr.Fn(x, beta0, beta1, beta2)}, x)
    x_lim <- seq(from=llim_s, to=rlim_s, length.out=100)
    y_lim <- mapply(FUN=function(x) {quadr.Fn(x, beta0, beta1, beta2)}, x_lim)
    
    plot(x, y, type="l", xlab="S", ylab="c", main=paste0("[",r,"]"))
    lines(x_lim, y_lim, col="green", lwd=3)
    lines(x_lim, y_lim)
    points(x=sss_s, y=sss_c, col="red", cex=2.0, pch=3)
    abline(v=c(llim_s, rlim_s), lwd=3, lty=2, col="green")
    abline(v=0, col="black")
    abline(v=arg_extr, col="purple")
    if (!is.na(floor_root1)) abline(v=c(floor_root1, floor_root2), col="grey50")
    if (!is.na(ceiling_root1)) abline(v=c(ceiling_root1, ceiling_root2), col="blue")
    if (!is.na(ssu_s)) abline(v=ssu_s, col="red")
    abline(h=0, col="grey50")
    abline(h=cons_lim, col="blue")
    abline(a=0, b=delta, col="red", lty=2)
  }
}
rSelected.Vect <- c(6,8,9,21,23,24)
