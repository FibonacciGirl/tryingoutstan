source('scripts/visualize_models.R')
source('scripts/visualize_priors multi.R')
source('scripts/make_data.R')
source('scripts/Data.R')

subject.data<-read.csv('generated_data/extracted-data-niw.csv')
subject.data$cond<-NULL

colnames(subject.data)<-c('j','t','p','rt')
subject.data$p <- sapply(subject.data$p, function(x){ if(x==1){return(0)};if(x==3){return(1)}})
subject.data<-subject.data[-which(subject.data$rt>2000),]

prior.1<-list('intercept' = list(fn = 'gamma', par1 =800, par2 =200 ),
              'proportion' = list(fn =  'gamma', par1 =1, par2 = .2),
              'base' = list(fn = 'unif', par1 =0, par2=1),
              'rate' = list(fn = 'gamma', par1= .3,par2 =.1),
              'base.1' = list(fn = 'unif', par1 =0, par2=1),
              'rate.1' = list(fn = 'exp', par1= 3,par2 =.1),
              'sigma' = list(fn = 'gamma', par1=25, par2=5),
              'jump.proportion' = list(fn = 'exp', par1=20, par2=1),
              'split' = list(fn = 'unif', par1=0, par2=72) )

write_json(prior.1,'prior.1.json')

grid.arrange(grobs = visualize.priors.multi(list('prior.1.json')))

params<-generate.parameter.values(read_json('prior.1.json'))






plot.data.fits(fit = params, model = c('power.constant'))
plot.data.fits(fit = params, model = c('power.power'))
plot.data.fits(fit = params, model = c('power.logistic'))
plot.data.fits(fit = params, model = c('piecewise.power.constant'))
