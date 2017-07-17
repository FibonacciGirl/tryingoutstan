subject.data<-read.csv('generated_data/extracted-data-niw.csv')
subject.data$cond<-NULL



colnames(subject.data)<-c('j','t','p','rt')

subject.data$p <- sapply(subject.data$p, function(x){ if(x==1){return(0)};if(x==3){return(1)}})
subject.data<-subject.data[-which(subject.data$rt>2000),]


Tr = max(subject.data$t)
J  = length(unique(subject.data$j))
P  = length(unique(subject.data$p))
N  = length(subject.data$rt)
rt = subject.data$rt
jj = subject.data$j
tt = subject.data$t
pp = subject.data$p



alpha= 27
beta= 1

source('scripts/DBDA2E-utilities.R')
gammaShRaFromModeSD(25,5)


mu0 = 65
tau0= .8
mu1 = 1
tau1=1
mu2 = 10
tau2=33
mu3 = 27
tau3=26
mu4 = 27
tau4=26
mu5 = 3
tau5=1
mu6 = 1
tau6=1
mu7 = 0
tau7= Tr


alpha0 = 26
beta0= .25
alpha1 = 6
beta1 = 48
alpha2=6
beta2=48
alpha3=6
beta3=48
alpha4=6
beta4=48
alpha5=6
beta5=48
alpha6=6
beta6=48
alpha7=6
beta7=.5

a = 8
b = 1.5

p1=.5
p2=.5
p=c(p1,p2)
nu=J


nchains = 1

model.data<- list(rep(list(
  Tr = Tr,
  J=J,
  P=P,
  N=N,
  rt = rt,
  jj=jj,
  tt=tt,
  pp=pp,
  p=p,
  alpha,
  beta,
  a = a,
  b = b,
  mu0 = mu0,
  tau0 = tau0,
  mu1 = mu1,
  tau1 = tau1,
  mu2 = mu2,
  tau2 = tau2,
  mu3 = mu3,
  tau3 = tau3,
  mu4 = mu4,
  tau4 = tau4,
  mu5 = mu5,
  tau5 = tau5,
  mu6 = mu6,
  tau6 = tau6,
  mu7 = mu7,
  tau7 = tau7,

  
  alpha0 = alpha0,
  beta0 = beta0,
  alpha1 = alpha1,
  beta1 = beta1,
  alpha2 = alpha2,
  beta2 = beta2,
  alpha3 = alpha3,
  beta3 = beta3,
  alpha4 = alpha4,
  beta4 = beta4,
  alpha5 = alpha5,
  beta5=beta5,
  alpha6 = alpha6,
  beta6=beta6,
  alpha7 = alpha7,
  beta7=beta7
  
  
), nchains))

library(rstan)




rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

subject.fit.test <- stan(file = 'rt_comparison_constantRLR-logisticRLR-Vectorize.stan', data = model.data, iter = 5000,warmup =500 , 
                  chains = nchains, verbose = T)

write_json(subject.fit.test, 'rtu_intercept')

plot(subject.fit.test,par = c('rtu_intercept'))

plot(subject.fit.test,par = c('log_q_z1[11]','log_q_z2[11]'))
summary(subject.fit.test,par = c('pi0'))

remove(subject.fit.test)

pr<-matrix(nrow=J,ncol=2)
softmax2 <- function(x){ exp(x) / sum(exp(x))}
for(j in 1:J){
  for(t in 1:Tr){
    fit<-as.matrix(subject.fit.test, par = c(paste0('log_q_z1[',j,']'),paste0('log_q_z2[',j,']')))
  }
  probs<-matrix(nrow=4500,ncol=2)
  
  for(i in 1:4500){
    probs[i,]<-softmax2(fit[i,])
  }
  
  pr[j,1]<-mean(probs[,1])
  pr[j,2]<-mean(probs[,2])
}



pr



library(tidyr)
library(dplyr)

source('scripts/visualize_models.R')
plot.data.fits(subject.data = subject,fit=fit, model= c('power.logistic'))

subject<-subset(subject.data, subject == 2)
fit= list(power.logistic = list(intercept =.41 , base = .35 , rate = .81, upper.proportion = 1.03, lower.proportion = 33 , rate.1 = .435, split = 18,sigma = 40 ))



subject.data<- subject.data %>% spread('p','rt')
colnames(subject.data) <- c('subject', 't', 'unpredictable','predictable')


  
learners <- which(round(pr[,2],2) > .2)
p<-list()
plot.subject<-list()
i=255
for(l in learners){
  i=i+1
  plot.subject<-subset(subject.data, subject == l)

  p[[i]]<-plot.data.fits(plot.subject,fit=NULL)
}

library(gridExtra)
plot(arrangeGrob(grobs = p))
