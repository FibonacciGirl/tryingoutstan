source('scripts/visualize_models.R')
source('scripts/visualize_priors multi.R')
source('scripts/make_data.R')
source('scripts/fit_conversion.R')
library(gridExtra)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



# enter values for prior distributions on types of parameters 

#means
prior.m<-list('intercept' = list(fn = 'gamma', par1 =900, par2 =100 ), 
            'base' = list(fn = 'unif', par1 =0, par2=1),
            'rate' = list(fn = 'gamma', par1= .3,par2 =.3),
            'proportion' = list(fn =  'gamma', par1 =1, par2 = .1),
            'jump' = list(fn = 'beta', par1=.5, par2=30),
            'split' = list(fn = 'unif', par1=0, par2=72),
            'sigma' = list(fn = 'gamma', par1=25, par2=10))

#variances
prior.v<-list('intercept' = list(fn = 'gamma', par1 =100, par2 =10 ), 
            'base' = list(fn = 'beta', par1 =.2, par2=7),
            'rate' = list(fn = 'beta', par1= .2,par2 =7),
            'proportion' = list(fn =  'beta', par1 =.1, par2 = 7),
            'jump' = list(fn = 'beta', par1=.1, par2=7),
            'split' = list(fn = 'gamma', par1=10, par2=10),
            'sigma' = list(fn = 'gamma', par1=10, par2=10))



write_json(prior.m, 'prior.m.json')
write_json(prior.v, 'prior.v.json')

##plot prior distributions 
prior.plots('prior.m.json')
prior.plots('prior.v.json')





##draw fake data from the priors
model = c('power.power') #specify model
t = 1:72 
n.subjects = 20

fake.data<-generate.fake.data(prior.m,model,t,n.subjects)

params<-fake.data$params
subject.data<-fake.data$fake.data


#view fake data plots
which.subject= 10

graph.fake.data(subject.data,params,which.subject)

##run stan on fake data
nchains = 20

fake.model.fit<-get.stan(fake.data= subject.data, prior.mean = prior.m, prior.var = prior.v, model=model, nchains)


plot(fake.model.fit,par=c('rtu_intercept[2]'))

## graph fit for subject

graph.fit.data(fake.model.fit,subject.data,which.subject,model)

##summary of subject fit
fit.summary<-summary( fake.model.fit , pars = c(paste0('sigma[',which.subject,']'),
                                     paste0('rtu_intercept[',which.subject,']'),
                                     paste0('rtu_base[',which.subject,']'),
                                     paste0('rtu_rate[',which.subject,']'),
                                     paste0('rlr_constant_intercept[',which.subject,']'),
                                     paste0('rlr_power_intercept[',which.subject,']'),
                                     paste0('rlr_power_base[',which.subject,']'),
                                     paste0('rlr_power_rate[',which.subject,']')))


##mcmc stuff

fit.mcmc<-As.mcmc.list(fake.model.fit, pars = c('sigma'))

plot(fit.mcmc)


##mean pi for all subjects
J=2*n.subjects
pi<-matrix(nrow=J,ncol=2)

for(j in 1:J){
  fit<-as.matrix(fake.model.fit, par = c(paste0('pi0[',j,',1]') , paste0('pi0[',j,',2]')))
  
  
  pi[j,1]<-mean(fit[,1])
  pi[j,2]<-mean(fit[,2])
}
pi[27,]


pi
