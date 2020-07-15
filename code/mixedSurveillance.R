library(glmmTMB)

mixedSurveillance<-function(formula,data,week=max(data$t),alpha,s=1.5){
  
  ###
  # data must be ordered in site and time
  # time covariate must be named "t" and must represent weeks
  # the count outcome must be called "outcome"
  # the site factor used in the random effect must be called "site"
  # formula must not include seasonality and must include the random effect for the site
  ###
  
  # First step: include seasonality as in Noufaily, Statis. Med. (2012)
  dataRun<-data[data$t<=week,]
  ## attribute a number of week of year to each week of the data set
  factorSeas<-data.frame(t=min(dataRun$t):max(dataRun$t),
                         week=rep(1:52,ceiling((max(dataRun$t)-min(dataRun$t)+1)/52))[1:(max(dataRun$t)-min(dataRun$t)+1)])
  t0<-tail(factorSeas$week,1) #t0 is the week of interest                       
  ## attribute a factor to each week of the data set according to Noufaily, Statis. Med. (2012)
  factorSeas2<-data.frame(week=c(t0:52,1:(t0-1)),
                          factorSeas=as.factor(c(rep(1,4),unlist(lapply(seq(2,10),function(x) rep(x,5))),rep(1,3))))
  ## merge both dataset to attribute a factor to each week of the dataset
  factorSeas<-merge(factorSeas,factorSeas2,by="week")
  factorSeas<-factorSeas[order(factorSeas$t),]
  ## add the factor in the original dataset
  dataRun<-merge(dataRun,factorSeas,by="t")
  dataRun<-dataRun[order(dataRun$t),]
  levels(dataRun$factorSeas)<-paste("factorSeas",1:10,sep="")
  
  # Second step: drop the most recent observations from the model evaluation
  ## the last 26 weeks are dropped
  dataRun$keep<-unlist(lapply(dataRun$t,function(x) ifelse(x>max(dataRun$t)-26,0,1)))
  
  # Third step: running the model
  formula<-paste(formula,"factorSeas",sep="+")
  ## first run of the model
  m.glmm<-glmmTMB(as.formula(formula),data=dataRun[dataRun$keep==1,],family=nbinom2())
  ## computing the weights
  theta<-summary(m.glmm)$sigma
  pred<-predict(m.glmm,type="response")
  dataRun$W<-1
  dataRun[dataRun$keep==1,]$W<-weightsGlmm(dataRun[dataRun$keep==1,]$outcome,pred,theta,s=s)
  ## second run of the model (with weights)
  m.glmm<-glmmTMB(as.formula(formula),data=dataRun[dataRun$keep==1,],family=nbinom2(),weights=W)      
  
  # Fourth step: computing the upper-bound
  theta<-summary(m.glmm)$sigma
  pred<-predict(m.glmm,newdata=dataRun[dataRun$t==week,],type="response")
  dataRes <- data.frame(outcome=dataRun[dataRun$t==week,]$outcome,
                        pred=pred,
                        upp_pred=unlist(lapply(pred,function(x) qnbinom((1+alpha)/2 , mu=x,size=theta))),
                        t=rep(week,length(pred)),
                        site=dataRun[dataRun$t==week,]$site,
                        theta=rep(theta,length(pred)))
  
  return(dataRes)
  
}



# weights function
weightsGlmm<-function(Y,mu_pred,thetas,s,type="pearson"){
  # weights based on pearson residuals (baseline scenario)
  if(type=="pearson"){
    sd <- sqrt(mu_pred + mu_pred^2/thetas)
    r<-(Y-mu_pred)/sd
    weights<-ifelse(r<s,1,1/r^2)
    return(length(mu_pred)/sum(weights)*weights)}
  # weights based on ccdf
  if(type=="ccdf"){
      quant_pred<-unlist(lapply(1:length(mu_pred), function(x) 1-pnbinom(Y[x]-1,mu=mu_pred[x],size=thetas)))
      weights<-ifelse(quant_pred>s,1,quant_pred/s)
      return(length(mu_pred)/sum(weights)*weights)}
}
