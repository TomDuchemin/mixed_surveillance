# simulation procedure used for the paper
library(dplyr)
library(glmmTMB)

# simulate data without outbreaks
simulateData<-function(N,T,u_i,paramBeta,paramPhi){
  # First step: simulation of covariates
  ## parameters of the covariates distribution
  m1<- runif(n = N,min=30,max=50)/100 ; p2 <- runif(n = T)
  ## covariates
  X1 <- unlist(lapply(m1,function(x) rnorm(n=T,mean=x,sd=1)))
  X2 <- rep(unlist(lapply(p2,function(x) rbinom(n=1,prob=x,size=1))),N)
  ## trend
  t <- rep(seq(1,T),N)
  scaled_t <- (t-mean(t))/sd(t)
  ## random effect
  site <- as.factor(unlist(lapply(1:N,function(x) rep(x,T))))
  u_i<-u_i[site]
  ## Fourier coefficients for seasonality
  cos1<- cos(2*pi*t/52) ; cos2<- cos(2*2*pi*t/52)
  sin1<- sin(2*pi*t/52) ; sin2<- sin(2*2*pi*t/52)
  
  ## Data Set
  data <- data.frame(X1, X2, t, cos1, cos2, sin1, sin2, site, scaled_t, u_i,
                     obs = factor(1:length(X1)))
  
  ## Simulation of the outcome
  data_log<-paramBeta[1]+paramBeta[2]*data$X1+paramBeta[3]*data$X2+paramBeta[4]*data$t+paramBeta[5]*data$cos1+paramBeta[6]*data$cos2+paramBeta[7]*data$sin1+paramBeta[8]*data$sin2+data$u_i
  data$mean<-exp(data_log)
  ## we define the theta with an overdispersion Phi, to mimic the simulation of Noufaily, Statis. Med. (2012)
  theta <- mean(data$mean)/(paramPhi-1)
  data$outcome <- unlist(lapply(1:(N*T),function(x) rnbinom(1,mu=data$mean[x],size=theta)))

  return(data)
}


# add outbreaks to baseline data
adding_outbreak<-function(data, k_baseline=8, k_current=8,beta=paramBeta,
                          overdisp_phi=paramPhi){
  
  site<-data$site[!duplicated(data$site)]
  
  # covariates  to indicate where there are outbreaks
  data$outbreak<-FALSE
  # compute mean and sd of each observation
  mean_obs <- exp(paramBeta[1]+paramBeta[2]*data$X1+paramBeta[3]*data$X2+paramBeta[4]*data$t+paramBeta[5]*data$cos1+paramBeta[6]*data$cos2+paramBeta[7]*data$sin1+paramBeta[8]*data$sin2+data$u_i)+1
  theta <- mean(mean_obs)/(overdisp_phi-1)
  sd_obs <- sqrt(mean_obs+mean_obs^2/theta)
  # looping to add outbreaks

  for(i in 1:length(site)){
    # row numbers of the site we want to add outbreaks to
    index<-c(min(which(data$site==site[i])),max(which(data$site==site[i])))
    
    # week when the outbreaks occur
      t_baseline <- sample(index[1]:(index[2]-(49+1)),4)
      t_baseline <- t_baseline[order(t_baseline)]
      
      for(i_baseline in 1:4){
        # size of the outbreak
        size_baseline<-0
        while(size_baseline<2){
          size_baseline <- rpois(1,k_baseline*sd_obs[t_baseline[i_baseline]]) #size of the outbreak
        }
        
        # duration of outbreaks
        outbk<-rlnorm(size_baseline,meanlog=0,sdlog=0.5)
        h <- hist(outbk,breaks=seq(0,ceiling(max(outbk)),1),plot=FALSE)
        case <- h$count
        #add outbreak to the data
        data$outcome[t_baseline[i_baseline]:(t_baseline[i_baseline]+length(case)-1)]<-case+data$outcome[t_baseline[i_baseline]:(t_baseline[i_baseline]+length(case)-1)]
        
      }
    
    # Include the outbreak in current day
    t_current <- sample((index[2]-(49+1)):(index[2]-3),1)
    size_current<-0
    while(size_current<2){
      size_current <- rpois(1,k_baseline*sd_obs[t_current]) #size of the outbreak
    }
    outbk<-rlnorm(size_current,meanlog=0,sdlog=0.5)
    h <- hist(outbk,breaks=seq(0,ceiling(max(outbk)),1),plot=FALSE)
    case <- h$count
    duration_current<-length(case)
    duration_current<-ifelse(duration_current+t_current-1>index[2],index[2]-t_current+1,duration_current)
    t_end<-t_current+duration_current-1
    data$outcome[t_current:t_end]<-data$outcome[t_current:t_end]+case[1:duration_current]
    data$outbreak[t_current:t_end]<-TRUE # note when there is an outbreak
  }
  return(data)
}


## A simulation example: scenario 7 
set.seed(10)
### set parameters values
paramBeta<-c(1,-0.5,0.5,0.0075,0.02,0.02,0.02,0.02)
paramPhi<-1.5
N<-50
T<-6*52
### simulate random effects
u_i<-rnorm(50,0,0.5)
### simulate data without outbreaks
data<-simulateData(N,T,u_i,paramBeta,paramPhi)
### add outbreaks
data<-adding_outbreak(data)

summary(data)
summary(data$outcome[data$outbreak])
summary(data$outcome[!data$outbreak])
