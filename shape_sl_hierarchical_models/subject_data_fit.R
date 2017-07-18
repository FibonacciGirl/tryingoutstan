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





source('scripts/DBDA2E-utilities.R')
source('scripts/visualize_priors multi.R')
gammaShRaFromModeSD(25,5)

prior.1<-list('intercept' = list(fn = 'gamma', par1 =900, par2 =100 ),
              'proportion' = list(fn =  'gamma', par1 =1, par2 = .2),
              'base' = list(fn = 'unif', par1 =0, par2=1),
              'rate' = list(fn = 'exp', par1= 6,par2 =0),
              'sigma' = list(fn = 'gamma', par1=25, par2=5),
              'jump.proportion' = list(fn = 'exp', par1=3, par2=1),
              'split' = list(fn = 'unif', par1=0, par2=72) )

write_json(prior.1,'prior.1.json')
grid.arrange(grobs = visualize.priors.multi(list('prior.1.json')) )


alpha= 26.96
beta= 1.04

mu0 = 82.99
tau0= 0.09
mu1 = 1
tau1=1
mu2 = 6
tau2= 6
mu3 = 26.96
tau3=25.98
mu4 = 26.96
tau4= 25.98
mu5 = 6
tau5=6
mu6 = 1
tau6= 1
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
p_prior=c(p1,p2)
nu=J


nchains = 1

model.data.1<- list(rep(list(
  Tr = Tr,
  J=J,
  P=P,
  N=N,
  rt = rt,
  jj=jj,
  tt=tt,
  pp=pp,
  p_prior=p_prior,
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
model.data.1[[1]]$tau2


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

subject.fit.test <- stan(file = 'rt_comparison_constantRLR-logisticRLR-Vectorize.stan', data = model.data.1, iter = 10000,warmup =1000 , 
                  chains = nchains, verbose = T)


plot(subject.fit.test,par = c('rlr_logistic_intercept[5]'))
plot(subject.fit.test,par = c('pi0[201,1]','pi0[201,2]'))

plot(subject.fit.test,par = c('log_q_z1[11]','log_q_z2[11]'))
summary(subject.fit.test,par = c('pi0'))



pr<-matrix(nrow=J,ncol=2)
softmax2 <- function(x){ exp(x) / sum(exp(x))}
for(j in 1:J){
  for(t in 1:Tr){
    fit<-as.matrix(subject.fit.test, par = c(paste0('log_q_z1[',j,']'),paste0('log_q_z2[',j,']')))
  }
  probs<-matrix(nrow=900,ncol=2)
  
  for(i in 1:900){
    probs[i,]<-softmax2(fit[i,])
  }
  
  pr[j,1]<-mean(probs[,1])
  pr[j,2]<-mean(probs[,2])
}

pi<-matrix(nrow=J,ncol=2)

for(j in 1:J){
  fit<-as.matrix(subject.fit.test, par = c(paste0('pi0[',j,',1]') , paste0('pi0[',j,',2]')))
  

  pi[j,1]<-mean(fit[,1])
  pi[j,2]<-mean(fit[,2])
}




hist(pi[,1])
hist(pi[,2])

library(tidyr)
library(dplyr)

source('scripts/visualize_models.R')
plot.data.fits(subject.data = subject,fit=fit, model= c('power.logistic'))

subject<-subset(subject.data, subject == 2)
fit= list(power.logistic = list(intercept =.41 , base = .35 , rate = .81, upper.proportion = 1.03, lower.proportion = 33 , rate.1 = .435, split = 18,sigma = 40 ))



subject.data<- subject.data %>% spread('p','rt')
colnames(subject.data) <- c('subject', 't', 'unpredictable','predictable')





learners  
learners <- which(pi[,1] > .8)

length(learners)
p<-list()
plot.subject<-list()
i=0


while(i<9){
  i=i+1

  plot.subject<-subset(subject.data, subject == learners[i])
  p[[i]]<-plot.data.fits(subject.data = plot.subject)
}
plot.data.fits(subject.data = plot.subject)



library(gridExtra)
plot(arrangeGrob(grobs = p))
