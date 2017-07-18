#function that excludes subjects based on criterion
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


