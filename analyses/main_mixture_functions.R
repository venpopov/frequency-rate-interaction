# Paul Bays' function for transforming K to SD (in radians)
k2sd <- function (K) {
  S <- matrix(0,1,length(K))
  for (j in 1:length(K)) {
    if (K[j]==0) S[j] = Inf
    if (is.infinite(K[j])) S[j] = 0    
    if (K[j] >= 0 & !is.infinite(K[j])) {
      S[j] = sqrt(-2*log(besselI(K[j],1)/besselI(K[j],0)));    
    }
  }
  S
}

sd2k <- function(S) {
  R = exp(-S**2/2)
  K = 1/(R**3-4*R**2+3*R)
  K[R < 0.85] = -0.4 + 1.39 * R[R < 0.85] + 0.43/(1-R[R <0.85])
  K[R < 0.53] = 2 * R[R < 0.53] + R[R < 0.53]**3 + (5*R[R < 0.53]**5)/6
  K
}


dst <- function(x, mu, sigma, df) { dt((x - mu)/sigma, df = df)/sigma }

# function for mixture model likelihood
mixtureLL <- function(dat, empirical=FALSE) {
  # transform x from degrees to radians
  rad = dat[c('y','V1','V2','V4','V5')] * pi / 180
  center=0
  if (empirical) {
    center = -1*pi/180
    rad[c('V1','V2','V4','V5')] = rad[c('V1','V2','V4','V5')]-1*pi/180
  }
  function(p_correct=2.5, p_other=1, sigma=9) {
    # trasnform laten prob
    p_c = exp(p_correct)/(exp(p_correct)+exp(p_other)+exp(0))
    p_o = exp(p_other)/(exp(p_correct)+exp(p_other)+exp(0))
    p_g = exp(0)/(exp(p_correct)+exp(p_other)+exp(0))
    # transform the normal sd into radians kappa for circular vonmises concentration parameter
    rad_sigma = sigma * pi /180
    kappa = (1/rad_sigma) ** 2
    l_norm <- brms::dvon_mises(rad$y, mu=center, kappa=kappa)
    l_norm1 <- brms::dvon_mises(rad$y, mu=rad$V1, kappa=kappa)
    l_norm2 <- brms::dvon_mises(rad$y, mu=rad$V2, kappa=kappa)
    l_norm4 <- brms::dvon_mises(rad$y, mu=rad$V4, kappa=kappa)
    l_norm5 <- brms::dvon_mises(rad$y, mu=rad$V5, kappa=kappa)
    l_unif <- dunif(rad$y, min=-pi, max=pi)
    likelihood <- p_c*l_norm + p_o/4*(l_norm1+l_norm2+l_norm4+l_norm5) + p_g*l_unif
    -sum(log(likelihood))
  }
}



# function for mixture model likelihood
mixtureLL_st <- function(dat, empirical=FALSE) {
  # transform x from degrees to radians
  rad = dat[c('y','V1','V2','V4','V5')] * pi / 180
  center = 0
  if (empirical) {
    center = -1*pi/180
    rad[c('V1','V2','V4','V5')] = rad[c('V1','V2','V4','V5')]-1*pi/180
  }
  function(p_correct=2.5, p_other=1, sigma=9, df=1) {
    # trasnform laten prob
    p_c = exp(p_correct)/(exp(p_correct)+exp(p_other)+exp(0))
    p_o = exp(p_other)/(exp(p_correct)+exp(p_other)+exp(0))
    p_g = exp(0)/(exp(p_correct)+exp(p_other)+exp(0))
    # transform the normal sd into radians kappa for circular vonmises concentration parameter
    rad_sigma = sigma * pi /180
    l_norm <- dst(rad$y, mu=center, sigma=rad_sigma, df=df)
    l_norm1 <- dst(rad$y, mu=rad$V1, sigma=rad_sigma, df=df)
    l_norm2 <- dst(rad$y, mu=rad$V2, sigma=rad_sigma, df=df)
    l_norm4 <- dst(rad$y, mu=rad$V4, sigma=rad_sigma, df=df)
    l_norm5 <- dst(rad$y, mu=rad$V5, sigma=rad_sigma, df=df)
    l_unif <- dunif(rad$y, min=-pi, max=pi)
    likelihood <- p_c*l_norm + p_o/4*(l_norm1+l_norm2+l_norm4+l_norm5) + p_g*l_unif
    if(sum(is.na(log(likelihood)))>0) {print(c(p_correct, p_other,sigma, df))}
    -sum(log(likelihood))
  }
}

# function for mixture model likelihood
mixtureLL_st4 <- function(dat, empirical=FALSE) {
  # transform x from degrees to radians
  rad = dat[c('y','V1','V2','V4','V5')] * pi / 180
  center = 0
  if (empirical) {
    center = -1*pi/180
    rad[c('V1','V2','V4','V5')] = rad[c('V1','V2','V4','V5')]-1*pi/180
  }
  function(p_correct=2.5, p_other_near=0.5, p_other_far=0.5, sigma=9, df=1) {
    # trasnform laten prob
    p_c = exp(p_correct)/(exp(p_correct)+exp(p_other_near)+exp(p_other_far)+exp(0))
    p_on = exp(p_other_near)/(exp(p_correct)+exp(p_other_near)+exp(p_other_far)+exp(0))
    p_of = exp(p_other_far)/(exp(p_correct)+exp(p_other_near)+exp(p_other_far)+exp(0))
    p_g = exp(0)/(exp(p_correct)+exp(p_other_near)+exp(p_other_far)+exp(0))
    # transform the normal sd into radians kappa for circular vonmises concentration parameter
    rad_sigma = sigma * pi /180
    l_norm <- dst(rad$y, mu=center, sigma=rad_sigma, df=df)
    l_norm1 <- dst(rad$y, mu=rad$V1, sigma=rad_sigma, df=df)
    l_norm2 <- dst(rad$y, mu=rad$V2, sigma=rad_sigma, df=df)
    l_norm4 <- dst(rad$y, mu=rad$V4, sigma=rad_sigma, df=df)
    l_norm5 <- dst(rad$y, mu=rad$V5, sigma=rad_sigma, df=df)
    l_unif <- dunif(rad$y, min=-pi, max=pi)
    likelihood <- p_c*l_norm + p_on/2*(l_norm2+l_norm4) + p_of/2*(l_norm1+l_norm5) + p_g*l_unif
    -sum(log(likelihood))
  }
}

# function for mixture model likelihood
mixtureLL4_st_vp <- function(dat, empirical=F) {
  # transform x from degrees to radians
  rad = dat[c('y','V1','V2','V4','V5')] * pi / 180
  function(p_correct=2.5, p_other=1, sigma_correct=9, sigma_wrong=9, df=1) {
    # trasnform laten prob
    p_c = exp(p_correct)/(exp(p_correct)+exp(p_other)+exp(0))
    p_o = exp(p_other)/(exp(p_correct)+exp(p_other)+exp(0))
    p_g = exp(0)/(exp(p_correct)+exp(p_other)+exp(0))
    # transform the normal sd into radians kappa for circular vonmises concentration parameter
    rad_sigma = sigma_correct * pi /180
    rad_sigma_wrong = sigma_wrong * pi /180
    center = ifelse(empirical, median(dat$angle_diff)*pi/180, 0)
    l_norm <- dst(rad$y, mu=center, sigma=rad_sigma, df=df)
    l_norm1 <- dst(rad$y, mu=rad$V1, sigma=rad_sigma_wrong, df=df)
    l_norm2 <- dst(rad$y, mu=rad$V2, sigma=rad_sigma_wrong, df=df)
    l_norm4 <- dst(rad$y, mu=rad$V4, sigma=rad_sigma_wrong, df=df)
    l_norm5 <- dst(rad$y, mu=rad$V5, sigma=rad_sigma_wrong, df=df)
    l_unif <- dunif(rad$y, min=-pi, max=pi)
    likelihood <- p_c*l_norm + p_o/4*(l_norm1+l_norm2+l_norm4+l_norm5) + p_g*l_unif
    -sum(log(likelihood))
  }
}


# function for mixture model likelihood
mixtureLL4 <- function(dat, empirical=F) {
  # transform x from degrees to radians
  rad = dat[c('y','V1','V2','V4','V5')] * pi / 180
  function(p_correct=2.5, p_other=1, sigma_correct=9, sigma_wrong=9) {
    # trasnform laten prob
    p_c = exp(p_correct)/(exp(p_correct)+exp(p_other)+exp(0))
    p_o = exp(p_other)/(exp(p_correct)+exp(p_other)+exp(0))
    p_g = exp(0)/(exp(p_correct)+exp(p_other)+exp(0))
    # transform the normal sd into radians kappa for circular vonmises concentration parameter
    rad_sigma = sigma_correct * pi /180
    kappa = (1/rad_sigma) ** 2
    rad_sigma_wrong = sigma_wrong * pi /180
    kappa_wrong = (1/rad_sigma_wrong) ** 2
    center = ifelse(empirical, median(dat$angle_diff)*pi/180, 0)
    l_norm <- brms::dvon_mises(rad$y, mu=center, kappa=kappa)
    l_norm1 <- brms::dvon_mises(rad$y, mu=rad$V1, kappa=kappa_wrong)
    l_norm2 <- brms::dvon_mises(rad$y, mu=rad$V2, kappa=kappa_wrong)
    l_norm4 <- brms::dvon_mises(rad$y, mu=rad$V4, kappa=kappa_wrong)
    l_norm5 <- brms::dvon_mises(rad$y, mu=rad$V5, kappa=kappa_wrong)
    l_unif <- dunif(rad$y, min=-pi, max=pi)
    likelihood <- p_c*l_norm + p_o/4*(l_norm1+l_norm2+l_norm4+l_norm5) + p_g*l_unif
    -sum(log(likelihood))
  }
}

mixtureLL4_ <- function(dat) {
  # transform x from degrees to radians
  rad = dat[c('y','V1','V2','V4','V5')] * pi / 180
  function(p_correct=0.7, p_other=0,  sigma_correct=9, sigma_wrong=9) {
    # constrain probabilities to be between 0 and 1
    if(any(p_correct+p_other >= 1) | any(p_correct+p_other < 0)) {
      return(10e6)
    }
    # transform the normal sd into radians kappa for circular vonmises concentration parameter
    rad_sigma = sigma_correct * pi /180
    kappa = (1/rad_sigma) ** 2
    rad_sigma_wrong = sigma_wrong * pi /180
    kappa_wrong = (1/rad_sigma_wrong) ** 2
    l_norm <- brms::dvon_mises(rad$y, mu=0, kappa=kappa)
    l_norm1 <- brms::dvon_mises(rad$y, mu=rad$V1, kappa=kappa_wrong)
    l_norm2 <- brms::dvon_mises(rad$y, mu=rad$V2, kappa=kappa_wrong)
    l_norm4 <- brms::dvon_mises(rad$y, mu=rad$V4, kappa=kappa_wrong)
    l_norm5 <- brms::dvon_mises(rad$y, mu=rad$V5, kappa=kappa_wrong)
    l_unif <- dunif(rad$y, min=-pi, max=pi)
    likelihood <- p_correct*l_norm + p_other/4*(l_norm1+l_norm2+l_norm4+l_norm5) + (1-p_correct-p_other)*l_unif
    -sum(log(likelihood))
  }
}

# function for mixture model likelihood
mixtureLL_fixed_centers <- function(dat) {
  # transform x from degrees to radians
  rad = dat[c('y','V1','V2','V4','V5')] * pi / 180
  function(p_correct=0.7, p_other=0,  sigma=9) {
    # constrain probabilities to be between 0 and 1
    if(any(p_correct+p_other >= 1) | any(p_correct+p_other < 0)) {
      return(10e6)
    }
    # transform the normal sd into radians kappa for circular vonmises concentration parameter
    rad_sigma = sigma * pi /180
    kappa = (1/rad_sigma) ** 2
    l_norm <- brms::dvon_mises(rad$y, mu=0, kappa=kappa)
    l_norm1 <- brms::dvon_mises(rad$y, mu=mean(rad$V1), kappa=kappa)
    l_norm2 <- brms::dvon_mises(rad$y, mu=mean(rad$V2), kappa=kappa)
    l_norm4 <- brms::dvon_mises(rad$y, mu=mean(rad$V4), kappa=kappa)
    l_norm5 <- brms::dvon_mises(rad$y, mu=mean(rad$V5), kappa=kappa)
    l_unif <- dunif(rad$y, min=-pi, max=pi)
    likelihood <- p_correct*l_norm + p_other/4*(l_norm1+l_norm2+l_norm4+l_norm5) + (1-p_correct-p_other)*l_unif
    -sum(log(likelihood))
  }
}

# function for mixture model likelihood
mixtureLL_noguess <- function(dat) {
  # transform x from degrees to radians
  rad = dat[c('y','V1','V2','V4','V5')] * pi / 180
  function(p_correct=0.7, sigma=9) {
    # constrain probabilities to be between 0 and 1
    if(any(p_correct >= 1) | any(p_correct < 0)) {
      return(10e6)
    }
    p_other = 1-p_correct
    # transform the normal sd into radians kappa for circular vonmises concentration parameter
    rad_sigma = sigma * pi /180
    kappa = (1/rad_sigma) ** 2
    l_norm <- brms::dvon_mises(rad$y, mu=0, kappa=kappa)
    l_norm1 <- brms::dvon_mises(rad$y, mu=rad$V1, kappa=kappa)
    l_norm2 <- brms::dvon_mises(rad$y, mu=rad$V2, kappa=kappa)
    l_norm4 <- brms::dvon_mises(rad$y, mu=rad$V4, kappa=kappa)
    l_norm5 <- brms::dvon_mises(rad$y, mu=rad$V5, kappa=kappa)
    l_unif <- dunif(rad$y, min=-pi, max=pi)
    likelihood <- p_correct*l_norm + p_other/4*(l_norm1+l_norm2+l_norm4+l_norm5) + (1-p_correct-p_other)*l_unif
    -sum(log(likelihood))
  }
}

# fit and return parameter estimates as a mixture of only one normal and one uniform (e.g. only accuracy or guessing)
fit_mixture2_noguess <- function(dat, init_values=list(p_correct=0.7, sigma=9)) {
  require(stats4)
  LL_resp <- mixtureLL_noguess(dat) 
  fit <- mle(LL_resp, start = init_values, method = "L-BFGS-B", lower = c(0, 2), upper = c(1, 30), nobs=nrow(dat))
  coef <- data.frame(t(fit@coef))
  coef$negll <- summary(fit)@m2logL
  coef$p_other <- 1-coef$p_correct
  coef$AIC <- AIC(fit)
  coef$BIC <- BIC(fit)
  return(round(coef,3))
}

# fit and return parameter estimates as a mixture of only one normal and one uniform (e.g. only accuracy or guessing)
fit_mixture2 <- function(dat) {
  require(stats4)
  LL_resp <- mixtureLL(dat) 
  fit <- mle(LL_resp, start = list(p_correct=0.89, sigma=9), method = "L-BFGS-B", lower = c(0, 2), upper = c(1, 30), nobs=nrow(dat))
  coef <- data.frame(t(fit@coef))
  coef$negll <- summary(fit)@m2logL
  coef$p_guess <- 1-coef$p_correct
  coef$AIC <- AIC(fit)
  coef$BIC <- BIC(fit)
  return(round(coef,3))
}

# fit and return parameter estimates as a mixture of 6 distributions - 1 for correct, 4 for failed bindings, 1 unifor for guessing
fit_mixture3 <- function(dat, init_values=list(p_correct=2.5, p_other=1, sigma=9), empirical=FALSE) {
  require(stats4)
  LL_resp <- mixtureLL(dat, empirical=empirical) 
  # debug(LL_resp)
  fit <- mle(LL_resp, start = init_values, method='L-BFGS-B', lower = c(-Inf,-Inf, 2), upper = c(Inf,Inf, 30), nobs=nrow(dat))
  coef <- data.frame(t(fit@coef))
  p_correct <- exp(coef$p_correct)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  p_other <- exp(coef$p_other)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  coef$p_correct <- p_correct
  coef$p_other <- p_other
  coef$p_guess <- 1-coef$p_correct-coef$p_other
  coef$negll <- summary(fit)@m2logL
  coef$AIC <- AIC(fit)
  coef$BIC <- BIC(fit)
  return(round(coef,3))
}

# fit and return parameter estimates as a mixture of 6 distributions - 1 for correct, 4 for failed bindings, 1 unifor for guessing
fit_mixture3_st <- function(dat, init_values=list(p_correct=2.5, p_other=1, sigma=9, df=1), empirical=FALSE) {
  require(stats4)
  LL_resp <- mixtureLL_st(dat, empirical=empirical) 
  # debug(LL_resp)
  fit <- mle(LL_resp, start = init_values, method='L-BFGS-B', lower = c(-100,-100, 2,0.001), upper = c(100,100, 30, 900), nobs=nrow(dat))
  coef <- data.frame(t(fit@coef))
  p_correct <- exp(coef$p_correct)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  p_other <- exp(coef$p_other)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  coef$p_correct <- p_correct
  coef$p_other <- p_other
  coef$p_guess <- 1-coef$p_correct-coef$p_other
  coef$negll <- summary(fit)@m2logL
  coef$AIC <- AIC(fit)
  coef$BIC <- BIC(fit)
  return(round(coef,3))
}

# fit and return parameter estimates as a mixture of 6 distributions - 1 for correct, 4 for failed bindings, 1 unifor for guessing
fit_mixture4_st <- function(dat, init_values=list(p_correct=2.5, p_other_near=0.5, p_other_far=0.5, sigma=9, df=1), empirical=FALSE) {
  require(stats4)
  LL_resp <- mixtureLL_st4(dat, empirical=empirical) 
  # debug(LL_resp)
  fit <- mle(LL_resp, start = init_values, method='L-BFGS-B', lower = c(-Inf,-Inf,-Inf, 2), upper = c(Inf,Inf,Inf, 30), nobs=nrow(dat))
  coef <- data.frame(t(fit@coef))
  p_correct <- exp(coef$p_correct)/(exp(coef$p_correct)+exp(coef$p_other_near)+exp(coef$p_other_far)+exp(0))
  p_other_near <- exp(coef$p_other_near)/(exp(coef$p_correct)+exp(coef$p_other_near)+exp(coef$p_other_far)+exp(0))
  p_other_far <- exp(coef$p_other_far)/(exp(coef$p_correct)+exp(coef$p_other_near)+exp(coef$p_other_far)+exp(0))
  coef$p_correct <- p_correct
  coef$p_other_near <- p_other_near
  coef$p_other_far <- p_other_far
  coef$p_guess <- 1-coef$p_correct-coef$p_other_near-coef$p_other_far
  coef$negll <- summary(fit)@m2logL
  coef$AIC <- AIC(fit)
  coef$BIC <- BIC(fit)
  return(round(coef,3))
}

# fit and return parameter estimates as a mixture of 6 distributions - 1 for correct, 4 for failed bindings, 1 unifor for guessing
fit_mixture4_st_vp <- function(dat, init_values=list(p_correct=2.5, p_other=1, sigma_correct=9, sigma_wrong=9, df=1), empirical=F) {
  require(stats4)
  LL_resp <- mixtureLL4_st_vp(dat, empirical=empirical) 
  # debug(LL_resp)
  fit <- mle(LL_resp, start = init_values, method = "L-BFGS-B", lower = c(-Inf,-Inf, 2, 2), upper = c(Inf,Inf, 30, 30), nobs=nrow(dat))
  coef <- data.frame(t(fit@coef))
  p_correct <- exp(coef$p_correct)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  p_other <- exp(coef$p_other)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  coef$p_correct <- p_correct
  coef$p_other <- p_other
  coef$p_guess <- 1-coef$p_correct-coef$p_other
  coef$negll <- summary(fit)@m2logL
  coef$AIC <- AIC(fit)
  coef$BIC <- BIC(fit)
  return(round(coef,3))
}

# fit and return parameter estimates as a mixture of 6 distributions - 1 for correct, 4 for failed bindings, 1 unifor for guessing
fit_mixture3_fixed_centers <- function(dat) {
  require(stats4)
  LL_resp <- mixtureLL_fixed_centers(dat) 
  # debug(LL_resp)
  fit <- mle(LL_resp, start = list(p_correct=0.7, p_other=0.1, sigma=9), method = "L-BFGS-B", lower = c(0,0, 2), upper = c(1, 1, 30), nobs=nrow(dat))
  coef <- data.frame(t(fit@coef))
  coef$p_guess <- 1-coef$p_correct-coef$p_other
  coef$negll <- summary(fit)@m2logL
  coef$AIC <- AIC(fit)
  coef$BIC <- BIC(fit)
  return(round(coef,3))
}

# fit and return parameter estimates as a mixture of 6 distributions - 1 for correct, 4 for failed bindings, 1 unifor for guessing
fit_mixture4 <- function(dat, init_values=list(p_correct=2.5, p_other=1, sigma_correct=9, sigma_wrong=9), empirical=F) {
  require(stats4)
  LL_resp <- mixtureLL4(dat, empirical=empirical) 
  # debug(LL_resp)
  fit <- mle(LL_resp, start = init_values, method = "L-BFGS-B", lower = c(-Inf,-Inf, 2, 2), upper = c(Inf,Inf, 30, 30), nobs=nrow(dat))
  coef <- data.frame(t(fit@coef))
  p_correct <- exp(coef$p_correct)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  p_other <- exp(coef$p_other)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  coef$p_correct <- p_correct
  coef$p_other <- p_other
  coef$p_guess <- 1-coef$p_correct-coef$p_other
  coef$negll <- summary(fit)@m2logL
  coef$AIC <- AIC(fit)
  coef$BIC <- BIC(fit)
  return(round(coef,3))
}
