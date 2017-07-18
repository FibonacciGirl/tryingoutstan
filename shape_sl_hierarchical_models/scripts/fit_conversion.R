
prior.plots<-function(json){
  p<-visualize.priors.multi( list(json))
  
  plot<-plot(arrangeGrob(grobs = p))
  return(plot)
}


graph.fake.data<-function(subject.data,params, which.subject){

  plot.subject.data=subset(subject.data,subject == which.subject)
  model=unique(plot.subject.data$model)
  
  p<-plot.data.fits(subject.data=plot.subject.data, fit = params[[which.subject]], model = model )
  
  return(p)
}




extract.fit.parameters<-function(fit,which.subject,model){
  
  
  
  if(model == 'power.logistic'){
  fit.summary<-summary( fit , pars = c(paste0('sigma[',which.subject,']'),
                                       paste0('rtu_intercept[',which.subject,']'),
                          paste0('rtu_base[',which.subject,']'),
                          paste0('rtu_rate[',which.subject,']'),
                          paste0('rlr_constant_intercept[',which.subject,']'),
                          paste0('rlr_logistic_intercept[',which.subject,']'),
                          paste0('rlr_logistic_intercept2[',which.subject,']'),
                          paste0('rlr_logistic_rate[',which.subject,']'),
                          paste0('rlr_logistic_split[',which.subject,']')
                          ))
  }
  if(model == 'power.power'){
    fit.summary<-summary( fit , pars = c(paste0('sigma[',which.subject,']'),
                                         paste0('rtu_intercept[',which.subject,']'),
                                         paste0('rtu_base[',which.subject,']'),
                                         paste0('rtu_rate[',which.subject,']'),
                                         paste0('rlr_constant_intercept[',which.subject,']'),
                                         paste0('rlr_power_intercept[',which.subject,']'),
                                         paste0('rlr_power_base[',which.subject,']'),
                                         paste0('rlr_power_rate[',which.subject,']')))
  }
  
  model<-c(model, 'power.constant')
  
  parameter.values<- list()
  i=0
  
  for(m in model){
    i=i+1
    if(m == 'power.constant'){
      parameter.values[[i]] = list(intercept = fit.summary$summary[2],
                                            base= fit.summary$summary[3],
                                            rate= - fit.summary$summary[4],
                                            proportion = fit.summary$summary[5],
                                            sigma = fit.summary$summary[1])
    }
    if(m == 'piecewise.power.constant'){
     parameter.values[[i]] =  list(intercept = r.prior('intercept',priors),
                                             base= r.prior('base',priors), 
                                             rate= -r.prior('rate',priors),
                                             proportion = r.prior('proportion',priors),
                                             jump = r.prior('jump',priors),
                                             split = r.prior('split',priors),
                                             sigma=r.prior('sigma',priors))
    }
    if(m == 'power.power'){
      parameter.values[[i]] = list(intercept = fit.summary$summary[2],
                                            base= fit.summary$summary[3],
                                            rate= - fit.summary$summary[4],
                                            proportion = fit.summary$summary[6],
                                            base.1=  fit.summary$summary[7],
                                            rate.1= -fit.summary$summary[8],
                                            sigma = fit.summary$summary[1])
    }
    if(m == 'power.logistic'){
      parameter.values[[i]] =  list(intercept = fit.summary$summary[2],
                                             base= fit.summary$summary[3],
                                             rate= - fit.summary$summary[4],
                                             proportion = fit.summary$summary[6],
                                             jump = fit.summary$summary[7],
                                             rate.1= -fit.summary$summary[8],
                                             split = fit.summary$summary[9],
                                             sigma= fit.summary$summary[1])
      
    }
  }
  names(parameter.values) = model
  return(parameter.values)                        
  
}


graph.fit.data<-function(fit, subject.data, which.subject, model){
  
  explore<-extract.fit.parameters(fake.model.fit,which.subject,model)
  
  plot.data<-subset(subject.data, subject == which.subject)
  p<-plot.data.fits(subject.data=plot.data,fit = explore,model = c('power.constant',model))
  return(p)
}


