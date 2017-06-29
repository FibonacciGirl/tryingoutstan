library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())







Tr = 72
J= 80
alpha= 27
beta= .25
mu0 = 1
alpha0 = 1
alpha1 = 1
tau0= 100
beta0= 1
beta1 = 1
a = 101
b = 10
a0 = 27
b0 = .5
a1 = 2
b1 = 32
a2 = 2
b2=32

gammaShRaFromModeSD(50,10)

priors<- list(proportion = list( fn = 'unif',par1= .5, par2= 2),
              upper.proportion = list( fn = 'unif',par1= .5, par2= 2),
              lower.proportion = list( fn = 'unif',par1= .3, par2=.8 ),
              intercept = list(fn = 'gamma',par1 = 700, par2 = 100),
              base = list (fn = 'beta',par1 = 1 , par2 =1 ),
              rate = list (fn = 'beta',par1 = 1 , par2 =1 ),
              base.1 = list (fn = 'unif',par1 = 0 , par2 =1 ),
              rate.1 = list (fn = 'unif',par1 = .5 , par2 =1 ),
              jump = list (fn = 'unif',par1 = 0 , par2 =.75 ),
              split = list(fn = 'pois', par1 = 36 , par2 = 0),
              sigma = list (fn = 'gamma', par1 = 50, par2 =10 )
)



temp<-c()
for(i in 1:(J/4)){
params<-generate.parameter.values(priors)
temp<-c(temp,create.recovery.data(params,t=1:72)$unpredictable)

}
temp
params
rtu<-matrix(nrow=J,ncol=Tr)
for(j in 1:J){
  rtu[j,]<-temp[1:72]
  temp<-tail(temp, -72)
}
rtu
model.data<- list(list(
  Tr = Tr,
  J=J,
  rtu = rtu,
  alpha=alpha,
  beta=beta,
  mu0 = mu0,
  beta0 = beta0,
  beta1 = beta1,
  tau0=tau0,
  alpha0=alpha0,
  alpha1 = alpha1,
  a = a,
  b = b,
  a0 = a0,
  b0 = b0,
  a1 = a1,
  b1 = b1,
  a2 = a2,
  b2 = b2

),
list(
  Tr = Tr,
  J=J,
  rtu = rtu,
  alpha=alpha,
  beta=beta,
  mu0 = mu0,
  beta0 = beta0,
  beta1 = beta1,
  tau0=tau0,
  alpha0=alpha0,
  alpha1 = alpha1,
  a = a,
  b = b,
  a0 = a0,
  b0 = b0,
  a1 = a1,
  b1 = b1,
  a2 = a2,
  b2 = b2
  
),
list(
  Tr = Tr,
  J=J,
  rtu = rtu,
  alpha=alpha,
  beta=beta,
  mu0 = mu0,
  beta0 = beta0,
  beta1 = beta1,
  tau0=tau0,
  alpha0=alpha0,
  alpha1 = alpha1,
  a = a,
  b = b,
  a0 = a0,
  b0 = b0,
  a1 = a1,
  b1 = b1,
  a2 = a2,
  b2 = b2
  
))



fit.2 <- stan(file = 'rt_unpredictable_reparam.stan', data = model.data, iter = 10000,warmup =1000 , 
            chains = 1, verbose = T)


mcmc.fit.2 <- As.mcmc.list(fit.2)
plot(mcmc.fit.2)

plot.data.fits(subset(create.recovery.data(params,t=1:72),model=='power.constant'))
params
subset(create.recovery.data(params,t=1:72),model== 'power.constant')
