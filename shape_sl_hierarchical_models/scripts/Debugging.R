source('scripts/visualize_models.R')
source('scripts/visualize_priors multi.R')
source('scripts/make_data.R')
library(gridExtra)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



# enter values for prior distributions on types of parameters 

prior.m<-list('intercept' = list(fn = 'gamma', par1 =900, par2 =100 ), 
            'base' = list(fn = 'unif', par1 =0, par2=1),
            'rate' = list(fn = 'gamma', par1= .3,par2 =.3),
            'proportion' = list(fn =  'gamma', par1 =1, par2 = .1),
            'jump' = list(fn = 'beta', par1=.5, par2=30),
            'split' = list(fn = 'unif', par1=0, par2=72),
            'sigma' = list(fn = 'gamma', par1=25, par2=10))

prior.v<-list('intercept' = list(fn = 'gamma', par1 =100, par2 =10 ), 
            'base' = list(fn = 'beta', par1 =.2, par2=7),
            'rate' = list(fn = 'beta', par1= .2,par2 =7),
            'proportion' = list(fn =  'beta', par1 =.1, par2 = 7),
            'jump' = list(fn = 'beta', par1=.1, par2=7),
            'split' = list(fn = 'gamma', par1=10, par2=10),
            'sigma' = list(fn = 'gamma', par1=10, par2=10))




##plot prior distributions 
write_json(prior.m, 'prior.m.json')
write_json(prior.v, 'prior.v.json')
p<-visualize.priors.multi( list('prior.v.json'))
p<-visualize.priors.multi( list('prior.m.json'))
plot(arrangeGrob(grobs =p))



##draw fake data from the priors
model = c('power.logistic') #specify model
t = 1:72 
n.subjects = 20


fake.data.l<-fake.subject.data(prior.m,model,t,n.subjects)
fake.data.c<-fake.subject.data(prior.m,c('power.constant'),t,n.subjects)



#view fake data plots
params.l<-fake.data.l$params
params.c<-fake.data.c$params


fake.subject.l<-subset(fake.data.l$fake.data,model ==model)
fake.subject.c<-subset(fake.data.c$fake.data,model =='power.constant')
fake.subject.c$subject<-fake.subject.c$subject+20

fake.subject.c

which.subject= 3

pl<-plot.data.fits(subject.data=subset(fake.subject.l,subject == which.subject), fit = params.l[[which.subject]], model = model )
pc<-plot.data.fits(subject.data=subset(fake.subject.c,subject == which.subject), fit = params.c[[which.subject]], model = c('power.constant') )
pl
##convert fake data to model data format
nchains = 2


fake.model.fit<-get.stan(fake.data= rbind(fake.subject.l,fake.subject.c),prior.mean=prior.m,prior.var=prior.v,model=model, nchains)


plot(fake.model.fit, par = c('pi0[4,1]','pi0[4,2]'))

