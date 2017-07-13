library(jsonlite)
library(ggplot2)
source('scripts/model_predict.R')
source('scripts/make_data.R')


plot.data.fits<-function(subject.data, fit=NULL,   model=  c('power.constant',
                                                             'power.logistic',
                                                             'power.power',
                                                             'piecewise.power.constant')){
    sub<- unique(subject.data$subject)
  
   title<-paste('Subject', sub)

  if(!is.null(fit) && any(!is.na(model))){
    
    line.fit <- expand.grid(t=1:72, model= model)
    
    
    line.fit$predictable <- mapply(function(m,t){
      if(m == 'power.constant'){
        params<-model.params(fit, 'power.constant')
        params<-head(params,-1)
        rtu<-mapply(u.power.model.predict,t, MoreArgs = list(params[1],params[2],params[3]))
        rl<-mapply(rl.constant.model.predict,t,MoreArgs= list(params[4]))
        rtp<-rl*rtu
        return(rtp)
      }
      if(m == 'power.logistic'){
        params<-model.params(fit, 'power.logistic')
        params<-head(params,-1)
        rtu<-mapply(u.power.model.predict,t,MoreArgs = list(params[1],params[2],params[3]))
        rl<-mapply(rl.logistic.model.predict,t,MoreArgs = list(params[4],params[5],params[6],params[7]))
        rtp<-rl*rtu
        return(rtp)
      }

      if(m == 'power.power'){
        params<-model.params(fit, 'power.power')
        params<-head(params,-1)
        rtu<-mapply(u.power.model.predict,t,MoreArgs = list(params[1],params[2],params[3]))
        rl<-mapply(rl.power.model.predict,t,MoreArgs = list(params[4],params[5],params[6]))
        rtp<-rl*rtu
        return(rtp)
      }
      if(m == 'piecewise.power.constant'){
        params<-model.params(fit, 'piecewise.power.constant')
        params<-head(params,-1)
      
        rtp<-numeric(length(t))
        t1<-head(t,params[6])
        t2<-tail(t,-params[6])
        
        
        rtu<-mapply(u.power.model.predict,t,MoreArgs = list(params[1],params[2],params[3]))
        rl1<-mapply(rl.constant.model.predict,t1,MoreArgs = list(params[4]))
        rl2<-mapply(rl.constant.model.predict,t2,MoreArgs = list(params[5]))
        rtp[t1]<-rl1*rtu[t1]
        rtp[t2]<-rl2*rtu[t2]
        return(rtp)
      }
    },line.fit$model, line.fit$t)
    
    line.fit$unpredictable <- mapply(function(m,t){
      if(m == 'power.constant'){
        params<-model.params(fit, 'power.constant')
        params<-head(params,-1)
        rtu<-mapply(u.power.model.predict,t, MoreArgs = list(params[1],params[2],params[3]))
        
        return(rtu)
      }

      if(m == 'power.logistic'){
        params<-model.params(fit, 'power.logistic')
        params<-head(params,-1)
        rtu<-mapply(u.power.model.predict,t,MoreArgs = list(params[1],params[2],params[3]))
        return(rtu)
      }

      if(m == 'power.power'){
        params<-model.params(fit, 'power.power')
        params<-head(params,-1)
        rtu<-mapply(u.power.model.predict,t,MoreArgs = list(params[1],params[2],params[3]))
        return(rtu)
      }

      if(m == 'piecewise.power.constant'){
        params<-model.params(fit, 'piecewise.power.constant')
        params<-head(params,-1)
        
        rtu<-mapply(u.power.model.predict,t,MoreArgs = list(params[1],params[2],params[3]))

        return(rtu)
      }

    },line.fit$model, line.fit$t)
    #line.fit$model<-sapply(line.fit$model,extract.model.name)
    
    # p<-list()
    # i<-0
    # for(m in model){
    #   i<-i+1
    # 
    #   plot.data<-subset(line.fit, model == m)
    # 
    #   print(subject.data)
    #   print(plot.data)
    # 
    # 
    #   p[[i]]<-ggplot()+
    #     geom_point(subject.data, aes(x=t,y=predictable))+
    #     geom_point(aes(x=t,y=unpredictable))+
    #     geom_line(plot.data, aes(x=t,y=predictable))+
    #     geom_line(plot.data, aes(x=t,y=unpredictable))+
    #     ggtitle(title=title, subtitle =  toString(m))+
    #     theme_bw()
    # }
    # 
    # p<-arrangeGrob(grobs = p)
    
    pred<-data.frame(t = line.fit$t , rt = line.fit$predictable, type = 'predictable',model = line.fit$model)
    unpred<-data.frame(t = line.fit$t, rt = line.fit$unpredictable, type = 'unpredictable',model = line.fit$model)
    subject.pred<-data.frame(t= subject.data$t, rt= subject.data$predictable, type = rep('predictable',length(subject.data$t)))
    subject.unpred<-data.frame(t= subject.data$t, rt = subject.data$unpredictable, type = rep('unpredictable',length(subject.data$t)))
    
    
    plot.data<-rbind(pred,unpred)
    sub.data<-rbind(subject.pred,subject.unpred)
    
    
    p<-ggplot()+
      geom_point(data=sub.data,aes(x=t, y =rt,col=type))+
      geom_line(data=plot.data,aes(x=t, y =rt,col=type,linetype=model ))+
      # ggtitle(title)+
      theme_bw()
    
    
    # p<-ggplot(subject.data, aes(x=t, y=predictable)) +
    #   geom_point()+
    #   geom_line(data=line.fit, aes(x=t, y=predictable, color=model), size=.5)+
    #   ggtitle(title)+
    #   theme_minimal()
    # up<-ggplot(subject.data, aes(x=t, y=unpredictable)) +
    #   geom_point()+
    #   geom_line(data=line.fit, aes(x=t, y=unpredictable, color=model), size=.5)+
    #   ggtitle(title)+
    #   theme_minimal()
  }
  # 
  else{
    subject.pred<-data.frame(t= subject.data$t, rt= subject.data$predictable, type = rep('predictable',length(subject.data$t)))
    subject.unpred<-data.frame(t= subject.data$t, rt = subject.data$unpredictable, type = rep('unpredictable',length(subject.data$t)))
    
    sub.data<-rbind(subject.pred,subject.unpred)
    
    
    p<-ggplot(sub.data) +
      geom_point(aes(x=t, y =rt,col=type))+
      ggtitle(title)+
      theme_minimal()
  }
  
  # p<-arrangeGrob(grobs = list(p,up))
  # return(plot(p))
  # 
  return(p)
}



plot.model.fit<-function(which.subject, model=  c('constant',
                                                  'power',
                                                  'exponential', 
                                                  #c('generalized.logistic'),
                                                  'constant.constant',
                                                  'constant.power',
                                                  'constant.exponential',
                                                  'power.constant',
                                                  'power.power',
                                                  'power.exponential',
                                                  'exponential.constant',
                                                  'exponential.power',
                                                  'exponential.exponential'))
{
  subject.data<-subset(model.data, subject == which.subject)
  
  fit.data<-read_json(paste('fit_data/Likelihood/subject',which.subject,'fit.likelihood.json',sep="."))
  p<-plot.data.fits(subject.data,fit.data)
  return(p)
}


all.data.fits.list<-function(which.subjects){
  likelihood.fits<-list()
  for(i in which.subjects){
    j<-which(which.subjects == i)
    
    likelihood.fits[[j]]<-paste('fit_data/Likelihood/subject',i,'fit.likelihood.json',sep='.')
  }
  return(likelihood.fits)
}


convert.recovery.data<-function(json){
  subject.data<-read_json(json)
  
  t<-c()
  model<-c()
  predictable<-c()
  unpredictable<-c()
  for(i in 1:length(subject.data)){
    t<-c(t,subject.data[[i]]$t)
    model<-c(model,subject.data[[i]]$model)
    predictable<-c(predictable,subject.data[[i]]$predictable)
    unpredictable<-c(unpredictable,subject.data[[i]]$unpredictable)
  }
  subject.data<-data.frame(t=t,model=model,unpredictable= unpredictable,predictable=predictable)
  return(subject.data)
}




lik.fit.function<-function(fit){
  
  lik<-c('power.constant'= fit$power.constant$Likelihood,
         'exponential.constant'=fit$exponential.constant$Likelihood,
         'power.logistic'= fit$power.logistic$Likelihood,
         'exponential.logistic'= fit$exponential.logistic$Likelihood,
         'power.power'= fit$power.power$Likelihood,
         'power.exponential'=fit$power.exponential$Likelihood,
         'exponential.power'=fit$exponential.power$Likelihood,
         'exponential.exponential'=fit$exponential.exponential$Likelihood,
         'piecewise.power.constant'=fit$piecewise.power.constant$Likelihood,
         'piecewise.exponential.constant'=fit$piecewise.exponential.constant$Likelihood
  )
  
  return(lik)
}


AICc<-function(likelihood, model){
  k<-numeric()
  
  if(model == 'power.constant'){
    k=5
  }
  if( model == 'exponential.constant'){
    k=5
  }
  if(model == 'power.power'){
    k=7
  }
  if(model =='power.exponential'){
    k=7
  }
  if(model == 'power.logistic'){
    k=8
  }
  if(model == 'exponential.logistic'){
    k=8
  }
  if(model == 'piecewise.exponential.constant'){
    k=7
  }
  if(model=='piecewise.power.constant'){
    k=7
  }
  if(model=='exponential.power'){
    k=7
  }
  if(model=='exponential.exponential'){
    k=7
  }
  
  n=72
  aic<-2*(k + likelihood)
  aicc<-AIC.correction(aic,n,k)
  return(aicc)  
}


AIC.correction<-function(aic,n,k){
  return(aic + (2*k*(k+1))/(n - k -1 ))
}


best.fit.AICc<-function(list.fits){
  n<-length(list.fits)
  best.fit<-array(dim=array(n,2))
  
  for(i in 1:n){
    
    fit<-read_json(list.fits[[i]])
    
    lik<-lik.fit.function(fit)
    
    model=c(  'power.constant',
              'exponential.constant',
              'power.logistic',
              'exponential.logistic',
              'power.power',
              'power.exponential',
              'exponential.power',
              'exponential.exponential',
              'piecewise.power.constant',
              'piecewise.exponential.constant')
    
    aicc<- mapply(AICc, lik, model)
    
    best.fit[i,1]<-names(which.min(aicc))
    best.fit[i,2]<-min(as.numeric(aicc))
  }
  best.fit<-data.frame(model = best.fit[,1], aic = best.fit[,2])
  
  return(best.fit)
}





##still need to fix to plot the figure
plot.group.model<-function(list.fits,model.data,model){
  p<-list()
  subjects<-which(best.fit.AICc(list.fits)$model == model)
  
  p[[1]]<-plot.group(model.data,subjects,model) 
  
  
  fits<-list()
  aicc<-c()
  i<-0
  for(s in subjects){
    i<-i+1
    fits[[i]]<-list.fits[[s]]
    
    lik<-lik.fit.function(fits[[i]])
    aicc<-mapply(AICc,lik,MoreArgs = list(model))
  }
  lik<-lik.fit.function(fits[[i]])
  aicc<-mapply(AIC.correction,lik,MoreArgs = list(model))
  
  
  
  subject.index<-list(s1<-which.min(maps$map),
                      s2<-which.min(maps$map[-s1]),
                      s3<-which.max(maps$map)[c(-s1,-s2)],
                      s4<-which.max(maps$map[c(-s1,-s2,-s3)]))
  
  
  i<-0
  for(s in subject.index){
    s.num<-as.numeric(s)
    
    i<-i+1
    
    
    sub<-maps$subject[s.num]
    
    
    s.data<-subset(model.data, subject == sub)
    
    
    fit<-read_json(fits[[s.num]])
    
    
    p[[i+1]]<-plot.data.fits(subject.data=s.data,fit,model)
    
  }
  
  multiplot(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]], layout = matrix(c(1,1,2,3,1,1,4,5),nrow=2,byrow=T))
  
}




plot.group<-function(all.data, fits, model){
  
  best.fit<-best.fit.AICc(fits)$model
  subjects<-which(best.fit == model)
  
  
  fit.data<-subset(all.data, subject %in% subjects)
  not.fit.data<-subset(all.data, !(subject %in% subjects))
  
  n.a<-length(all.data$t)
  n.m<-length(fit.data$t)
  n.n.m<-length(not.fit.data$t)
  
  all.data$group<-rep('all',n.a)
  fit.data$group<-rep('best fit subjects', n.m)
  not.fit.data$group<-rep('all with best fit subjects removed',n.n.m)
  
  
  all.pred<-data.frame(t = all.data$t , rt = all.data$predictable, type = rep('predictable',n.a),group = all.data$group)
  all.unpred<-data.frame(t = all.data$t, rt = all.data$unpredictable, type = rep('unpredictable',n.a),group = all.data$group)
  fit.pred<-data.frame(t= fit.data$t, rt= fit.data$predictable, type = rep('predictable',n.m), group = fit.data$group)
  fit.unpred<-data.frame(t= fit.data$t, rt = fit.data$unpredictable, type = rep('unpredictable',n.m), group = fit.data$group)
  not.fit.pred<-data.frame(t= not.fit.data$t, rt= not.fit.data$predictable, type = rep('predictable',n.n.m), group = not.fit.data$group)
  not.fit.unpred<-data.frame(t= not.fit.data$t, rt= not.fit.data$unpredictable, type = rep('unpredictable',n.n.m), group = not.fit.data$group)
  
  
  all.plot.data<-rbind(all.pred,all.unpred)
  fit.plot.data<-rbind(fit.pred,fit.unpred)
  not.fit.plot.data<-rbind(not.fit.pred,not.fit.unpred)
  
  plot.data<-rbind(all.plot.data,fit.plot.data,not.fit.plot.data)
  title<-extract.model.name(model)
  
  p<-ggplot(plot.data)+
    geom_smooth(aes(x=t,y=rt, col=group,linetype = type), se =F)+
    ggtitle(paste('All subjects vs. subjects best fit by', title,'model', sep=' '))+
    xlab ("Appearance Count")+
    ylab ('Response Time')+
    ylim(0,2000)+
    theme_bw()
  
  return(p)
}



exclude.subjects<-function(model.data){
  
  subjects<-c()
  n<-length(unique(model.data$subject))
  
  for(i in 1:n){
    which.subject = i
    subject.data<-subset(model.data, subject == which.subject)
    exclude<-0
    
    
    for(j in 1:6){
      time.data<-subset(subject.data$t, subject.data$t %in% c(((j-1)*(12)+1):(j*12)))
      
      if(length(time.data) < 4){
        exclude<-exclude+1
      }
    }
    if(exclude >= 1){
      model.data<-subset(model.data, !(subject %in% which.subject))
      subjects<-c(subjects,which.subject)
    }
    
  }
  #return(model.data)
  return(subjects)
}


