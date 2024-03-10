pathogen_model <- function(t,y,parms){
  
  # read in initial states and their names
  States<-array(y, dim=dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  
  
  # unify the time unit of parameter inputs
  if(parms$time.step =='month'){
    period=12
    length.step=30.44 #days
  }else if(parms$time.step =='week'){
    period=52.1775
    length.step=7 #days
  }
  
  # waning rate of maternal immunity (by time step)
  omega = 1/(parms$DurationMatImmunityDays/length.step) # wanning immunity from M --> Xss
  
  # aging rate (by time step)
  mu = 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(WidthAgeClassMonth*4.345)} 
  
  # rate of recovery of first infection
  gamma1 = 1/(parms$dur.days1/length.step) 
  # rate of recovery of second infection
  gamma2 = 1/(parms$dur.days2/length.step)
  # rate of recovery of third infection
  gamma3 = 1/(parms$dur.days3/length.step)  
  gamma4 = gamma3  
  
   
  
  
  #gamma3 stands for rate of recovery from subsequent infection
  
  # Relative risk of infection (2nd)
  sigma1=parms$sigma1 
  # Relative risk of infection (3rd)
  sigma2=parms$sigma2
  # Relative risk of infection (4th+)
  sigma3=parms$sigma3
  
   
  # Relative infectiousness (2nd)
  rho1=parms$rho1 
  # Relative infectiousness (3rd+)
  rho2=parms$rho2
   
  
  #Pull out the states for the model as vectors
  M <-  States[,'M']  
  
  Xss  <-  States[,'Xss'] 
  Xi1s <-  States[,'Xi1s']  
  Xs1s <-  States[,'Xs1s']
  Xi2s <-  States[,'Xi2s']
  Xs2s <-  States[,'Xs2s']
  Xi3s <-  States[,'Xi3s']
  Xs3s <-  States[,'Xs3s']
  Xi4s <-  States[,'Xi4s']
  
  
  N_ages <- length(M) # the number of age groups
  
  ## parameter related to force of infection ################
  # per capita transmission probability
  baseline.txn.rate = parms$baseline.txn.rate # <<---------- estimated parameter 1
  
  # transmission probability per unit time
  b <- baseline.txn.rate / (parms$dur.days1/length.step) 
   
  q=parms$q # q depends on transmission type 
  # (whether depends on population density or not)
  contact = parms$contact # c2 is the contact matrix
  # transmission probability per unit time in each age group
  beta <- (b/100)/(sum(yinit.matrix)^(1-q))*contact 

  
  # 100 is a scaling factor for the contact matrix we choose
  # (see Ginny's paper and Matlab code for details)
  # this does not matter because most likely
  # you will need to estimate baseline.txn.rate
  
  Amp=parms$Amp # seasonal amplitude
  phi=parms$phi # seasonal phase shift
  #seasonality
  seasonal.txn <- (1+Amp*cos(2*pi*(t-phi*period)/period))
   
  
  # seasonal transmission probability
  beta_a_i <- seasonal.txn * beta 
  
  
  infectiousN <- (Xi1s+rho1*Xi2s+rho2*Xi3s+rho2*Xi4s) / sum(States)
  
   
  # for frequency dependent transmission
  
  lambda <- infectiousN %*% beta_a_i # force of transmission
  lambda <- as.vector(lambda) # vectorize force of transmission
  
  
  # create a matrix to record the changing variables
  dy <- matrix(NA, nrow=N_ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)
  
  period.birth.rate <- 
    log(parms$PerCapitaBirthsYear[t,]+1)/period
  # get period birth rate from annual birth rate
  
  
  #um is death rate
  um=parms$um
  #mu represents aging to the next class
  Aging.Prop <- c(0,mu[1:(N_ages-1)])
  
  #um is death rate
  
  dy[,'M'] <- period.birth.rate*sum(States) - 
    (omega+(mu+um))*M +
    Aging.Prop*c(0,M[1:(N_ages-1)]) 
  
  dy[,'Xss'] <- omega*M -
    (lambda) * Xss -
    (mu + um) * Xss + 
    Aging.Prop*c(0,Xss[1:(N_ages-1)])  
  
  dy[,'Xi1s'] <-   lambda*Xss - 
    (gamma1 + mu + um) * Xi1s + 
    Aging.Prop*c(0,Xi1s[1:(N_ages-1)]) 
  
  dy[,'Xs1s'] <- gamma1*Xi1s - 
    sigma1*lambda*Xs1s - 
    (mu+um)*Xs1s + 
    Aging.Prop*c(0,Xs1s[1:(N_ages-1)])  
  
  dy[,'Xi2s'] <- sigma1*lambda*Xs1s - 
    gamma2*Xi2s -
    (mu + um)*Xi2s + 
    Aging.Prop*c(0,Xi2s[1:(N_ages-1)]) 
  
  dy[,'Xs2s'] <- gamma2*Xi2s - 
    sigma2*lambda*Xs2s -
    (mu+um)*Xs2s + 
    Aging.Prop*c(0,Xs2s[1:(N_ages-1)]) 
  
  dy[,'Xi3s'] <- sigma2*lambda*Xs2s -
    (gamma3 + mu+um )*Xi3s +  
    Aging.Prop*c(0,Xi3s[1:(N_ages-1)]) 
  
  dy[,'Xs3s'] <- gamma3*Xi3s +  
    gamma4*Xi4s -
    sigma3*lambda*Xs3s -
    (mu + um)*Xs3s + 
    Aging.Prop*c(0,Xs3s[1:(N_ages-1)])  
  
  dy[,'Xi4s'] <- sigma3*lambda*Xs3s - 
    gamma4*Xi4s - 
    (mu + um)*Xi4s + 
    Aging.Prop*c(0,Xi4s[1:(N_ages-1)])  
  
  
  
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}






