source('scripts/visualize_models.R')
source('scripts/visualize_priors multi.R')
source('scripts/make_data.R')
library(gridExtra)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



# enter values for prior distributions on types of parameters 

prior<-list('intercept' = list(fn = 'gamma', par1 =900, par2 =100 ), 
            'base' = list(fn = 'unif', par1 =0, par2=1),
            'rate' = list(fn = 'gamma', par1= .3,par2 =.3),
            'proportion' = list(fn =  'gamma', par1 =1, par2 = .1),
            'jump' = list(fn = 'beta', par1=.5, par2=30),
            'split' = list(fn = 'unif', par1=0, par2=72),
            'sigma' = list(fn = 'gamma', par1=25, par2=10))


##plot prior distributions 
write_json(prior, 'prior.json')
p<-visualize.priors.multi(list('prior.json'))
plot(arrangeGrob(grobs =p))



##draw fake data from the priors
model = c('power.logistic') #specify model
t = 1:72 
n.subjects = 10


fake.data<-fake.subject.data(prior,model,t,n.subjects)


#view fake data plots
params<-subset(fake.data$params,model ==model)

fake.subject<-subset(fake.data$fake.data,model ==model)
fake.subject

which.subject= 1
params[[2]]

plot.data.fits(subject.data=subset(fake.subject,subject == which.subject), fit = params[[which.subject]], model = model )

##convert fake data to model data format
nchains = 1
fake.model.data<-get.plot.data(fake.data= fake.subject,prior.mean=prior,prior.var=prior,model=model, nchains)

fake.model.data[[1]]$tau2

##run stan
if(model == 'power.logistic'){
  f = 'rt_comparison_constantRLR_logisticRLR-Vectorize-sparse-student-t.stan'
}
if(model == 'power.power'){
  f = 'rt_comparison_constantRLR_powerRLR-Vectorize-sparse-student-t.stan'
}

model.data<-fake.model.data

subject.fit.test <- stan(file = f, data = fake.model.data, iter = 1000,warmup =100 , 
                         chains = nchains, verbose = T)

