library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



source('make_data.R')

gammaShRaFromModeSD(20,10)

Tr = 72
J= 30
alpha= 5
beta= .25
mu0 = 700
tau0= 100
mu1 = 26
tau1=51
mu2 = 26
tau2=51
mu3 = 26
tau3=51
mu4 = 26
tau4=51
mu5 = 26
tau5=51
mu6 = 26
tau6=51
mu7 = 26
tau7=51


alpha0 = 101
beta0= 1
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
beta7=48

a = 101
b = 10

p1=.5
p2=.5
p=c(p1,p2)


priors<- list(proportion = list( fn = 'unif',par1= .8, par2= 1.2),
              upper.proportion = list( fn = 'unif',par1= .5, par2= 2),
              lower.proportion = list( fn = 'unif',par1= .3, par2=.8 ),
              intercept = list(fn = 'gamma',par1 = 700, par2 = 100),
              base = list (fn = 'unif',par1 = .1 , par2 =.3 ),
              rate = list (fn = 'unif',par1 = .1 , par2 =.2 ),
              base.1 = list (fn = 'unif',par1 = .99 , par2 =1 ),
              rate.1 = list (fn = 'unif',par1 = 0.01 , par2 =0.1 ),
              jump = list (fn = 'unif',par1 = 0 , par2 =.75 ),
              split = list(fn = 'pois', par1 = 36 , par2 = 0),
              sigma = list (fn = 'gamma', par1 = 20, par2 =10 )
)



params<-list()


rtu<-array(dim=c(J,Tr))
rtp<-array(dim = c(J,Tr))
rt<-c()
jj<-c()
tt<-c()
pp<-c()
model<-list()

for(j in 1:J){
params[[j]]<-generate.parameter.values(priors)
select<-rbinom(1,1,.8)
if(select == 1){
  model[[j]]= c('power.power')
}
if(select == 0){
  model[[j]]= c('power.constant')
}
rt_data<-create.recovery.data(params[[j]],t=1:72,models=model[[j]])
rtu[j,]<-rt_data$unpredictable
rtp[j,]<-rt_data$predictable
rt<- c(rt,cbind(rtu[j,],rtp[j,]))
jj<-c(jj,rep(j, 2*Tr))
tt<-c(tt,rep(1:Tr,2))
pp<-c(pp,rep(0,Tr),rep(1,Tr))
}

rt<-c(rt)
jj<-as.vector(jj)
tt<-c(tt)
pp<-as.vector(pp)

N<-length(rt)
P=1
tt
dim(model.data.1[[1]]$rtu)

model.data.1<- list(list(
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
  beta6=beta6

  
  )
)

fit.vec.7 <- stan(file = 'rt_comparison_constantRLR-powerRLR-Vectorize-sparse-student-t.stan', data = model.data.1, iter = 10000,warmup =1000 , 
            chains = 1, verbose = T)

pairs(fit.vec.9)

plot(fit.vec.7,par= c('log_q_z1'))
plot(fit.vec.7,par= c('pi0'))
plot(fit.vec.7,par= c('sigma'))
model
mcmc.fit.vec.7 <- As.mcmc.list(fit.vec.7)


i=1
plot(mcmc.fit.vec.7[,(i+1),drop=F])
summary(mcmc.fit.vec.7[,(i+1),drop=F])

model[[19]]
params[[19]]$power.logistic
params[[i]]$power.constant


model
pr
source('visualize_models.R')
i=10
subject.data<-data.frame(t=1:72,predictable=rtp[i,],unpredictable=rtu[i,])
plot.data.fits(subject.data,fit=params[[i]],model=model[[i]])

params[[i]]

round(true)

pr<-matrix(nrow=J,ncol=2)
softmax2 <- function(x){ exp(x) / sum(exp(x))}
for(j in 1:J){
  for(t in 1:Tr){
  fit<-as.matrix(fit.vec.7, par = c(paste0('log_q_z1[',j,']'),paste0('log_q_z2[',j,']')))
  }
  probs<-matrix(nrow=9000,ncol=2)

  for(i in 1:9000){
    probs[i,]<-softmax2(fit[i,])
  }
  
  pr[j,1]<-mean(probs[,1])
  pr[j,2]<-mean(probs[,2])
}
pr[,2]


z.prob<-as.data.frame(pr)
colnames(z.prob)<- c('one','zero')

select<-c()
for(j in 1:J){
  if(model[[j]]=='power.power'){this=1}
  else{this=0}
  select<-c(select,this)
  
}
select
true<-cbind(select,pr[,2] )
round(true,2)

true.val
a<-as.matrix(z.prob[,1],nrow=J,ncol=1)
b<-as.matrix(z.prob[,2],nrow=J,ncol=1)
a
z.prob
data.plot<-data.frame(prob=rbind(a,b))
data.plot
data.plot$val<-rbind(as.matrix(rep(0,20),nrow=20,ncol=1),as.matrix(rep(1,20),nrow=20,ncol=1))

ggplot(data.plot)+
  geom_col(aes(x=val,y=prob),width=.5)+
  theme_bw()



plot.data.fits(params[[1]]$power.power,model=='power.constant')
params
subset(create.recovery.data(params,t=1:72),model== 'power.constant')



model.data.3chain<- list(list(
  Tr = Tr,
  J=J,
  rtu = rtu,
  rtp=rtp,
  p1=p1,
  p2=p2,
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
  beta6=beta6
),list(
  Tr = Tr,
  J=J,
  rtu = rtu,
  rtp=rtp,
  p1=p1,
  p2=p2,
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
  beta6=beta6
),list(
  Tr = Tr,
  J=J,
  rtu = rtu,
  rtp=rtp,
  p1=p1,
  p2=p2,
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
  beta6=beta6
)
)



