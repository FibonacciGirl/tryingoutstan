source('scripts/visualize_models.R')
source('scripts/visualize_priors multi.R')
source('scripts/make_data.R')
source('scripts/Data.R')


library(geoR)

subject.data<-read.csv('generated_data/extracted-data-niw.csv')
subject.data$cond<-NULL

colnames(subject.data)<-c('j','t','p','rt')
subject.data$p <- sapply(subject.data$p, function(x){ if(x==1){return(0)};if(x==3){return(1)}})
subject.data<-subject.data[-which(subject.data$rt>2000),]

prior.1<-list('intercept' = list(fn = 'unif', par1 =800, par2 =801 ),
              'proportion' = list(fn =  'unif', par1 =.7, par2 =.71),
              'base' = list(fn = 'unif', par1 =.6, par2=.61),
              'rate' = list(fn = 'gamma', par1= .3,par2 =.35),
              'base.1' = list(fn = 'unif', par1 =.7, par2=.71),
              'rate.1' = list(fn = 'unif', par1= .5,par2 =.51),
              'sigma' = list(fn = 'gamma', par1=25, par2=5),
              'jump.proportion' = list(fn = 'unif', par1=0, par2=.3),
              'split' = list(fn = 'pois', par1=36, par2=72) )

write_json(prior.1,'prior.1.json')


#plot priors
prior.vals<-


grid.arrange(grobs = visualize.priors.multi(list('prior.1.json')))

pr<-convert.prior(read_json('prior.1.json'))

set.seed(15)
params<-generate.parameter.values(prior.1, model = c('power.constant','power.power','power.logistic'))
params$piecewise.power.constant=params$power.logistic

#plot model types
plot.data.fits(fit = params, model = c('power.constant'))
plot.data.fits(fit = params, model = c('power.power'))
plot.data.fits(fit = params, model = c('power.logistic'))
plot.data.fits(fit = params, model = c('piecewise.power.constant'))

#plot rlr
install.packages('RColorBrewer')
library(RColorBrewer)


t=1:72
plot.rlr('acclimation',params,t)
plot.rlr('power.constant',params,t)
plot.rlr('power.power',params,t)
plot.rlr('power.logistic',params,t)
plot.rlr('piecewise.power.constant',params,t)

scale_color_hue()

?ggplot

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


gg_color_hue(4)[1]

plot.rlr<-function(model,params,t){
  
  cols=gg_color_hue(2)

  
  t<-seq(from= min(t),to=max(t),by=.001)
  linetype= 'dashed'
  ylab = 'Relative Response Time'
  
  if(model == 'power.constant'){
    title = c('Constant Relative Learning Rate')
    par<-model.params(params, 'power.constant')
    rt<-mapply(rl.constant.model.predict,t,MoreArgs = list(par[4]) )
    ylim=c(0,1)
    col= cols[1]
  }
  if(model == 'power.power'){
    title = c('Power Law Relative Learning Rate')
    par<-model.params(params, 'power.power')
    rt<-mapply(rl.power.model.predict,t,MoreArgs = list(par[4], par[5],par[6]) )
    ylim=c(0,1)
    col= cols[1]
  }
  if(model == 'power.logistic'){
    title = c('Logistic Relative Learning Rate')
    par<-model.params(params, 'power.logistic')
    rt<-mapply(rl.logistic.model.predict,t,MoreArgs = list(par[4],par[5],par[6],par[7]) )
    ylim=c(0,1)
    col= cols[1]
  }
  if(model == 'piecewise.power.constant'){
    title = c('Constant Piecewise Relative Learning Rate')
    par<-model.params(params, 'piecewise.power.constant')
    rt<-mapply(rl.piecewise.constant.model.predict,t,MoreArgs = list(par[4],par[5],par[6]) )
  
    ylim=c(0,1)
    col= cols[1]
  }
  if(model == 'acclimation'){
    title = c('Power Law Acclimation')
    par<-model.params(params, 'power.constant')
    rt<-mapply(u.power.model.predict,t,MoreArgs = list(par[1],par[2],par[3]) )
    
    ylim=c(0,1000)
    col = cols[2]
    linetype= 'solid'
    ylab='Response Time'
  }
  
  
  plot.data<-data.frame(t=t,rt=rt)
  ?linetype
  p<-ggplot()+
    geom_line(data=plot.data,aes(x=t,y=rt),size=1,col=col,linetype= linetype)+
    ggtitle(title)+
    xlab ("Appearance Count")+
    ylab (ylab)+
    ylim(ylim)+
    theme_bw()
  

  return(p)
}

?scale_color_hue()
?geom_line

#plot subject examples
subject.data<- subject.data %>% spread('p','rt')
colnames(subject.data) <- c('subject', 't', 'unpredictable','predictable')


subject.index = c( 40,75, 19, 39, 157, 205)

plot.subject.examples<-function(subject.index){
  p<-list()
  i=0
  for(s in subject.index){
    i=i+1
    index.data<-subset(subject.data, subject == s)
    
    print(index.data)
    p[[i]]<-plot.data.fits(subject.data = index.data)
  }
  
  l<- ggplot()+
    grid.arrange(grobs = p)+
    theme_bw()
  
  return(l)
}


library(gridExtra)
plot.subject.examples(subject.index)



#plot group average

plot.group.average<-function(all.data){
  
  n=length(all.data)
  
  all.pred<-data.frame(t = all.data$t , rt = all.data$predictable, stimulus = rep('predictable',n))
  all.unpred<-data.frame(t = all.data$t, rt = all.data$unpredictable, stimulus = rep('unpredictable',n))
  
  all.plot.data<-rbind(all.pred,all.unpred)
  
  title<- 'Group Averages Response Rime'
  
  p<-ggplot(all.plot.data)+
    geom_smooth(aes(x=t,y=rt, col=stimulus), level = .99 )+
    ggtitle(title)+
    xlab ("Appearance Count")+
    ylab ('Response Time')+
    #ylim(0,2000)+
    theme_bw()
  
  return(p)
}
 ?geom_smooth
plot.group.average(subject.data)
