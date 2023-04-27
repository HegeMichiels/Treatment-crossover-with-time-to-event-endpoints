################################################
### Functions for the onestep estimator   #####
###        Hege Michiels                  #####
###############################################


##### Sandwich estimator ######

####Treatment policy #####

library(ivsacim)
library(pec)

#Returns Z-E(Z given T>t)
Z.term.t = function(Y.t,data){
  Z = data$Z
  pred = mean(Z*Y.t)/mean(Y.t)
  return(Z-pred)
}

#Function that returns the value of the estimating equation, for given value of beta
f.itt = function(p,dN, Y,hazards, data){
  return(mean(U.itt(p=p,dN=dN,Y=Y,hazards=hazards,data=data)))
}

#Returns Estimating equation for every patient integrated over all time points
U.itt = function(p,dN,Y,hazards,data){
  Z = data$Z
  integral = 0
  t.min = 0
  for(t in colnames(Y)){
    dN.t = dN[,t]
    Y.t = Y[,t]
    lambda.t = hazards[,t]
    Z.term.t = Z.term.t(Y.t=Y.t,data=data)
    integral = integral + U.itt.t(p=p,dN.t=dN.t,Y.t=Y.t,lambda.t=lambda.t,Z=Z,Z.term.t=Z.term.t,t=t,t.min=t.min)
    t.min = t
  }
  return(integral)
}

#Returns Estimating equation for every patient at time t
U.itt.t = function(p,dN.t,Y.t,lambda.t,Z,Z.term.t,t,t.min){
  beta.est = p
  return(Z.term.t*(dN.t-(Y.t*Z*beta.est+Y.t*lambda.t)*(as.numeric(t)-as.numeric(t.min))))
}

beta.estimate.itt  = function(data,lower,upper,hazards){
  #We first fit an ivsacim object to have dN and Y
  fit.ivsacim = ivsacim(time=data$time, event=data$status, instrument=data$Z, IV_valid = TRUE,treatment_init=data$Z, treatment_shift_time=data$shift ,n_sim = 1000)
  fit.ivsacim.by_prod = fit.ivsacim$by_prod
  stime = fit.ivsacim$stime
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  output = uniroot(f=f.itt,lower=lower,upper=upper,dN = dN, Y=Y,hazards = hazards,data=data)
  beta.est =output$root
  return(beta.est)
}

#Returns 
hazards.marginal = function(data){
  fit.ivsacim = ivsacim(time=data$time, event=data$status, instrument=data$Z, IV_valid = TRUE, treatment_init=data$Z, treatment_shift_time=data$shift ,n_sim = 1000)
  fit.ivsacim.by_prod = fit.ivsacim$by_prod
  stime = fit.ivsacim$stime
  
  #We fit an Aalen model to estimate the hazards
  fit.Aalen <- aalen(Surv(time, status) ~ 1, data = subset(data,Z==0),resample.iid=1)
  
  pred.Aalen = predict.timereg(fit.Aalen,resample.iid=1,newdata=data,times=stime)
  #gives survival predictions for survival models. Predictions given in matrix form with different subjects in different rows.
  survival.pred = pred.Aalen$S0 #gives survival predictions for survival models. Predictions given in matrix form with different subjects in different rows.
  colnames(survival.pred) = as.character(stime) #S(t given L)
  cum.hazards = apply(survival.pred, MARGIN = 2, FUN = function(x) -log(x))
  cum.hazards = t(apply(cum.hazards,MARGIN=1,FUN=cummax))
  hazards = cbind(cum.hazards[,1],t(apply(cum.hazards, 1, diff)))
  hazards = t(apply(hazards,MARGIN = 1,FUN = function(x) x/c(stime[1],diff(stime)))) #lambda(t given L)
  colnames(hazards) = as.character(stime)
  return(hazards)
}


#Function that returns the average of the derivative of the score function, for given value of beta
#p: value for beta
#fit.ivsacim.by_prod
f.itt.deriv = function(p,dN,Y,data){
  return(mean(U.itt.deriv(p=p,dN=dN,Y=Y,data=data)))
}

#Returns the derivative of the Estimating equation for every patient integrated over all time points
U.itt.deriv = function(p,dN,Y,data){
  integral = 0
  t.min = 0
  Z = data$Z
  for(t in colnames(Y)){
    dN.t = dN[,t]
    Y.t = Y[,t]
    Z.term.t = Z.term.t(Y.t=Y.t,data=data)
    integral = integral + U.itt.deriv.t(p=p,dN.t=dN.t,Y.t=Y.t,Z.term.t=Z.term.t,Z=Z,t=t,t.min=t.min,data=data)
    t.min = t
  }
  return(integral)
}

#Returns the derivative of the Estimating equation for every patient at time t
U.itt.deriv.t = function(p,dN.t,Y.t,Z.term.t,Z,t,t.min,data){
  beta.est = p
  diff = (as.numeric(t)-as.numeric(t.min))
  output = -Z.term.t*(Z*Y.t)*diff
  return(output)
}


p.itt.score = function(beta=0,fit.ivsacim.by_prod,data,stime,hazards){
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  
  U.vector = U.itt(p=beta,dN = dN, Y=Y,hazards = hazards,data=data)
  return(t.test(U.vector)$p.value)
}  



#Returns the SE for beta, using the sandwich estimator 
beta.itt.se.function = function(fit.ivsacim.by_prod,data,beta.est,hazards){
  dN = fit.ivsacim.by_prod$dN 
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  
  U.vector = U.itt(p=beta.est,dN = dN, Y=Y,hazards = hazards,data=data)
  var.U = var(U.vector)
  deriv.U.vector = U.itt.deriv(p=beta.est,dN = dN, Y=Y,data=data)
  var.beta = var.U/(dim(data)[1]*(mean(deriv.U.vector))^2)
  return(sqrt(var.beta))
}


f.p.value.itt.score = function(p,fit.ivsacim.by_prod,data,stime,hazards,p.value=0.05){
  beta = p
  return((p.itt.score(beta=beta,fit.ivsacim.by_prod=fit.ivsacim.by_prod,data=data,stime=stime,hazards=hazards)-p.value))
}

beta.bound.itt = function(lower,upper,fit.ivsacim.by_prod,data,stime,hazards,p.value=0.05){
  beta.bound = uniroot(f=f.p.value.itt.score,lower=lower,upper=upper,fit.ivsacim.by_prod=fit.ivsacim.by_prod,data=data,stime=stime,p.value=p.value,hazards=hazards)$root
  return(beta.bound)
}


####Ivsacim #####

add.Z.pred = function(data){
  fitZ.L <- glm(formula=Z~1, family="binomial", data=data)
  data$Z.pred <- predict(fitZ.L,newdata=data,type="response") #E(Z)
  return(data)
}

#Function that returns the value of the estimating equation, for given value of beta
#p: value for beta
#fit.ivsacim.by_prod
f.iv.score = function(p,dN,Y,D,D.int,data){
  return(mean(U.iv.score(p=p,dN=dN,D=D,D.int=D.int,Y=Y,data=data)))
}

#Returns Estimating equation for every patient integrated over all time points
U.iv.score = function(p,dN,D,D.int,Y,data){
  Z = data$Z
  Z.pred = data$Z.pred
  Z.c = Z-Z.pred
  integral = 0
  t.min = 0
  for(t in colnames(dN)){
    dN.t = dN[,t]
    D.t = D[,t]
    D.int.t = D.int[,t]
    Y.t = Y[,t]
    integral = integral + U.iv.score.t(p=p,dN.t=dN.t,D.t=D.t,D.int.t=D.int.t,Y.t=Y.t,Z.c=Z.c,t=t,t.min=t.min)
    t.min = t
  }
  return(integral)
}

#Returns Estimating equation for every patient at time t
U.iv.score.t = function(p,dN.t,D.t,D.int.t,Y.t,Z.c,t,t.min){
  beta.est = p
  return(Z.c*exp(beta.est*D.int.t)*(dN.t-(Y.t*D.t*beta.est)*(as.numeric(t)-as.numeric(t.min))))
}


beta.estimate.iv.score  = function(fit.ivsacim.by_prod,data,lower,upper,stime){
  data = add.Z.pred(data)
  
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime)
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  output = uniroot(f=f.iv.score,lower=lower,upper=upper,dN = dN,D=D, D.int=D.int,Y=Y,data=data)
  beta.est = output$root
  return(beta.est)
}


#Function that returns the average of the derivative of the score function, for given value of beta
#p: value for beta
#fit.ivsacim.by_prod
f.iv.deriv.score = function(p,dN,Y,D,D.int,data){
  return(mean(U.iv.deriv.score(p=p,dN=dN,D=D,D.int=D.int,Y=Y,data=data)))
}

#Returns the derivative of the score function for every patient integrated over all time points
U.iv.deriv.score = function(p,dN,D,D.int,Y,data){
  Z = data$Z
  Z.pred = data$Z.pred
  Z.c = Z-Z.pred
  integral = 0
  t.min = 0
  for(t in colnames(dN)){
    dN.t = dN[,t]
    D.t = D[,t]
    D.int.t = D.int[,t]
    Y.t = Y[,t]
    integral = integral + U.iv.deriv.score.t(p=p,dN.t=dN.t,D.t=D.t,D.int.t=D.int.t,Y.t=Y.t,Z.c=Z.c,t=t,t.min=t.min)
    t.min = t
  }
  return(integral)
}

#Returns Estimating equation for every patient at time t
U.iv.deriv.score.t = function(p,dN.t,D.t,D.int.t,Y.t,Z.c,t,t.min){
  beta.est = p
  diff = as.numeric(t)-as.numeric(t.min)
  return(Z.c*exp(beta.est*D.int.t)*Y.t*(D.int.t*(dN.t-beta.est*D.t*diff)-D.t*diff))
}

#Function that returns an estimate for beta
av.deriv.score.function.iv = function(fit.ivsacim.by_prod,data,beta.start,stime){
  data = add.Z.pred(data)
  
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  output = f.iv.deriv.score(p=beta.start,dN = dN, Y=Y,D=D,D.int=D.int,data=data)
  return(output)
}


#Returns the SE for beta, using the sandwich estimator 
beta.iv.se = function(fit.ivsacim.by_prod,data,beta.est){
  data = add.Z.pred(data)
  
  dN = fit.ivsacim.by_prod$dN 
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  
  U.vector = U.iv.score(p = beta.est,dN,D,D.int,Y,data)
  var.U = var(U.vector)
  deriv.U.vector = U.iv.deriv.score(p=beta.est,dN,D,D.int,Y,data)
  var.beta = var.U/(dim(data)[1]*(mean(deriv.U.vector))^2)
  return(sqrt(var.beta))
}

p.iv.score = function(beta=0,fit.ivsacim.by_prod,data,stime){
  data = add.Z.pred(data)
  
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  U.vector = U.iv.score(p=beta,dN = dN, Y=Y,D=D,D.int=D.int,data=data)
  return(t.test(U.vector)$p.value)
}


f.iv = function(p,fit.ivsacim.by_prod,data,stime,p.value=0.05){
  beta = p
  return((p.iv.score(beta=beta,fit.ivsacim.by_prod=fit.ivsacim.by_prod,data=data,stime=stime)-p.value))
}

beta.bound.iv = function(lower,upper,fit.ivsacim.by_prod,data,stime,p.value=0.05){
  beta.bound = uniroot(f=f.iv,lower=lower,upper=upper,fit.ivsacim.by_prod=fit.ivsacim.by_prod,data=data,stime=stime,p.value=p.value)$root
  return(beta.bound)
}

####Onestep estimator #####


#Returns Z-E(Z given T^0>t, L)
#Z.formula: Z ~ L 
#L.values = character vector with colnames of L values 
Z.term.strat = function(beta,Z,Y.t,D.int.t,Z.formula,L.values,data){
  if(is.null(L.values)){
    new.data = data.frame(cbind(Z = Z, Y.t = Y.t,D.int.t = D.int.t))
  }
  else{
  new.data = data.frame(cbind(Z = Z, Y.t = Y.t,D.int.t = D.int.t,data %>% select(L.values)))
  }
  fit.Z = glm(Z.formula,family=binomial,weights=Y.t*exp(beta*D.int.t),data=new.data)
  pred = predict(fit.Z,newdata = new.data,type="response")
  return(Z-pred)
}

#Function that returns the value of the estimating equation, for given value of beta
f = function(p,dN, Y,D,D.int,Z.formula,L.values,data){
  return(mean(U(p=p,dN=dN,D=D,D.int=D.int,Y=Y,Z.formula=Z.formula,L.values=L.values,data=data)))
}

#Returns Estimating equation for every patient integrated over all time points
U = function(p,dN,D,D.int,Y,Z.formula,L.values,data){
  integral = 0
  t.min = 0
  for(t in colnames(Y)){
    dN.t = dN[,t]
    D.t = D[,t]
    D.int.t = D.int[,t]
    Y.t = Y[,t]
    Z.term.t = Z.term.strat(beta=p,Z=data$Z,Y.t=Y.t,D.int.t=D.int.t,Z.formula=Z.formula,L.values=L.values,data=data)
    integral = integral + U.t(p=p,dN.t=dN.t,D.t=D.t,D.int.t=D.int.t,Y.t=Y.t,Z.term.t=Z.term.t,t=t,t.min=t.min)
    t.min = t
  }
  return(integral)
}

#Returns Estimating equation for every patient at time t
U.t = function(p,dN.t,D.t,D.int.t,Y.t,Z.term.t,t,t.min){
  beta.est = p
  return(Z.term.t*exp(beta.est*D.int.t)*(dN.t-(Y.t*D.t*beta.est)*(as.numeric(t)-as.numeric(t.min))))
}

av.function = function(fit.ivsacim.by_prod,data,beta.start,stime,Z.formula,L.values){
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  output = f(p=beta.start,dN = dN,Y=Y,D=D,D.int=D.int,Z.formula,L.values,data=data)
  return(output)
}



logit = function(x) log(x/(1-x))


#Function that returns the average of the derivative of the score function, for given value of beta
#p: value for beta
#fit.ivsacim.by_prod
f.deriv = function(p,dN, Y,D,D.int,Z.formula,L.values,data){
  return(mean(U.deriv(p=p,dN=dN,D=D,D.int=D.int,Y=Y,Z.formula=Z.formula,L.values=L.values,data=data)))
}

#Returns the derivative of the Estimating equation for every patient integrated over all time points
U.deriv = function(p,dN,D,D.int,Y,Z.formula,L.values,data){
  integral = 0
  t.min = 0
  Z = data$Z
  for(t in colnames(Y)){
    dN.t = dN[,t]
    D.t = D[,t]
    D.int.t = D.int[,t]
    Y.t = Y[,t]
    Z.term.t = Z.term.strat(beta=p,Z=data$Z,Y.t=Y.t,D.int.t=D.int.t,Z.formula=Z.formula,L.values=L.values,data=data)
    integral = integral + U.deriv.t(p=p,dN.t=dN.t,D.t=D.t,D.int.t=D.int.t,Y.t=Y.t,Z.term.t=Z.term.t,Z=Z,t=t,t.min=t.min,Z.formula=Z.formula,L.values=L.values,data=data)
    t.min = t
  }
  return(integral)
}

#Returns the derivative of the Estimating equation for every patient at time t
U.deriv.t = function(p,dN.t,D.t,D.int.t,Y.t,Z.term.t,Z,t,t.min,Z.formula,L.values,data){
  beta.est = p
  exp.term = exp(beta.est*D.int.t)
  deriv.Z.term.t = deriv.Z.term.t.strat(p=p,dN.t=dN.t,D.t=D.t,D.int.t=D.int.t,Y.t=Y.t,Z.term.t=Z.term.t,Z=Z,t=t,t.min=t.min,Z.formula=Z.formula,L.values=L.values,data=data)
  diff = (as.numeric(t)-as.numeric(t.min))
  output = -deriv.Z.term.t*exp.term*(dN.t-(Y.t*D.t*beta.est)*diff)+Z.term.t*exp.term*(D.int.t*(dN.t-Y.t*D.t*beta.est*diff)-Y.t*D.t*diff)
  return(output)
}

#Derivative of E(Z given T^0>=t,L) to beta
deriv.Z.term.t.strat = function(p,dN.t,D.t,D.int.t,Y.t,Z.term.t,Z,t,t.min,Z.formula,L.values,data){
  beta.est = p
  exp.term = exp(beta.est*D.int.t)
  pred.Z = Z-Z.term.t # E(Z given T^0>=t,L) = expit(gamma_t(beta)'L)
  
  deriv = deriv.gamma(p=p,dN.t=dN.t,D.t=D.t,D.int.t=D.int.t,Y.t=Y.t,Z.term.t=Z.term.t,Z=Z,t=t,t.min=t.min,Z.formula=Z.formula,L.values=L.values,data=data)
  L.matrix = cbind(1,data %>% select(L.values)) #Moet aangepast worden indien er interacties zitten in het model voor Z!!
  output = c((pred.Z/(1+exp(logit(pred.Z))))*t(deriv)%*%t(L.matrix))
  return(output)
}

#Returns dgamma_t(beta)/dbeta
deriv.gamma = function(p,dN.t,D.t,D.int.t,Y.t,Z.term.t,Z,t,t.min,Z.formula,L.values,data){
  beta.est = p
  exp.term = exp(beta.est*D.int.t)
  pred.Z = Z-Z.term.t # E(Z given T^0>=t,L) = expit(gamma_t(beta)'L)
  
  #part A
  part.1 =  Y.t*exp.term*D.int.t*Z.term.t
  if(is.null(L.values)){
    part.A = apply(rbind(part.1), MARGIN=1, FUN=mean)
  }else{
    data.new = data %>% select(L.values) #Moet aangepast worden indien er interacties zitten in het model voor Z!!
    matrix.2 = apply(data.new,MARGIN = 2,FUN = function(x) as.numeric(x)*part.1)
    part.A = apply(rbind(part.1,t(matrix.2)), MARGIN=1, FUN=mean)
  }
  
  #Part B
  part.1  = Y.t*exp.term*pred.Z/(1+exp(logit(pred.Z)))
  L.matrix = if(is.null(L.values)){cbind(rep(1,times=dim(data)[1]))}else{cbind(1,data.new)} #every row: c(1,STRATA1,STRATA2)
  part.B = t(as.matrix(L.matrix))%*%as.matrix(L.matrix*part.1)/dim(data)[1]
  #dgamma_t(beta)/dbeta
  deriv = solve(part.B)%*%part.A
  return(deriv)
}



#Function that returns an estimate for beta
av.deriv.function = function(fit.ivsacim.by_prod,data,beta.start,stime,Z.formula,L.values){
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  output = f.deriv(p=beta.start,dN = dN, Y=Y,D=D,D.int=D.int,Z.formula=Z.formula,L.values=L.values,data=data)
  return(output)
}


p.score = function(beta=0,fit.ivsacim.by_prod,data,stime,Z.formula,L.values){
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  U.vector = U(p=beta,dN = dN, Y=Y,D=D,D.int=D.int,data=data,Z.formula=Z.formula,L.values=L.values)
  return(t.test(U.vector)$p.value)
}

#Returns the SE for beta, using the sandwich estimator 
beta.se.function = function(fit.ivsacim.by_prod,data,beta.est,Z.formula,L.values){
  dN = fit.ivsacim.by_prod$dN 
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  U.vector = U(p=beta.est,dN = dN, Y=Y,D=D,D.int=D.int,data=data,Z.formula=Z.formula,L.values=L.values)
  var.U = var(U.vector)
  deriv.U.vector = U.deriv(p=beta.est,dN = dN, Y=Y,D=D,D.int=D.int,data=data,Z.formula=Z.formula,L.values=L.values)
  var.beta = var.U/(dim(data)[1]*(mean(deriv.U.vector))^2)
  return(sqrt(var.beta))
}


f.score = function(p,fit.ivsacim.by_prod,data,stime,Z.formula,L.values,p.value=0.05){
  beta = p
  return((p.score(beta=beta,fit.ivsacim.by_prod=fit.ivsacim.by_prod,data=data,stime=stime,Z.formula=Z.formula,L.values=L.values)-p.value))
}

beta.bound = function(lower,upper,fit.ivsacim.by_prod,data,stime,Z.formula,L.values,p.value=0.05){
  beta.bound = uniroot(f=f.score,lower=lower,upper=upper,fit.ivsacim.by_prod=fit.ivsacim.by_prod,data=data,stime=stime,p.value=p.value,Z.formula=Z.formula,L.values=L.values)$root
  return(beta.bound)
}


#####Onestep with hazard######


#Function that returns the value of the estimating equation, for given value of beta
f.hazard = function(p,dN, Y,D,D.int,Z.formula,L.values,hazards,data){
  return(mean(U.hazard(p=p,dN=dN,D=D,D.int=D.int,Y=Y,Z.formula=Z.formula,L.values=L.values,hazards=hazards,data=data)))
}

#Returns Estimating equation for every patient integrated over all time points
U.hazard = function(p,dN,D,D.int,Y,Z.formula,L.values,hazards,data){
  integral = 0
  t.min = 0
  for(t in colnames(Y)){
    dN.t = dN[,t]
    D.t = D[,t]
    D.int.t = D.int[,t]
    Y.t = Y[,t]
    lambda.t = hazards[,t]
    Z.term.t = Z.term.strat(beta=p,Z=data$Z,Y.t=Y.t,D.int.t=D.int.t,Z.formula=Z.formula,L.values=L.values,data=data)
    #Z.term.t = Z.term.t(Y.t,data)
    integral = integral + U.hazard.t(p=p,dN.t=dN.t,D.t=D.t,D.int.t=D.int.t,Y.t=Y.t,lambda.t=lambda.t,Z.term.t=Z.term.t,t=t,t.min=t.min)
    t.min = t
  }
  return(integral)
}

#Returns Estimating equation for every patient at time t
U.hazard.t = function(p,dN.t,D.t,D.int.t,Y.t,lambda.t,Z.term.t,t,t.min){
  beta.est = p
  return(Z.term.t*exp(beta.est*D.int.t)*(dN.t-(Y.t*D.t*beta.est+Y.t*lambda.t)*(as.numeric(t)-as.numeric(t.min))))
}

av.function.hazard = function(fit.ivsacim.by_prod,data,beta.start,stime,Z.formula,L.values,hazards){
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  output = f.hazard(p=beta.start,dN = dN,Y=Y,D=D,D.int=D.int,Z.formula=Z.formula,L.values=L.values,hazards=hazards,data=data)
  return(output)
}


hazards.conditional = function(data,T.formula){
  fit.ivsacim = ivsacim(time=data$time, event=data$status, instrument=data$Z, IV_valid = TRUE, treatment_init=data$Z, treatment_shift_time=data$shift ,n_sim = 1000)
  fit.ivsacim.by_prod = fit.ivsacim$by_prod
  stime = fit.ivsacim$stime
  
  #We fit an Aalen model to estimate the hazards
  fit.Aalen <- aalen(T.formula, data = subset(data,Z==0),resample.iid=1)
  
  pred.Aalen = predict.timereg(fit.Aalen,resample.iid=1,newdata=data,times=stime)
  #gives survival predictions for survival models. Predictions given in matrix form with different subjects in different rows.
  survival.pred = pred.Aalen$S0 #gives survival predictions for survival models. Predictions given in matrix form with different subjects in different rows.
  colnames(survival.pred) = as.character(stime) #S(t given L)
  cum.hazards = apply(survival.pred, MARGIN = 2, FUN = function(x) -log(x))
  cum.hazards = t(apply(cum.hazards,MARGIN=1,FUN=cummax))
  hazards = cbind(cum.hazards[,1],t(apply(cum.hazards, 1, diff)))
  hazards = t(apply(hazards,MARGIN = 1,FUN = function(x) x/c(stime[1],diff(stime)))) #lambda(t given L)
  colnames(hazards) = as.character(stime)
  return(hazards)
}



#Function that returns the average of the derivative of the score function, for given value of beta
#p: value for beta
#fit.ivsacim.by_prod
f.deriv.hazard = function(p,dN, Y,D,D.int,Z.formula,L.values,hazards,data){
  return(mean(U.deriv.hazard(p=p,dN=dN,D=D,D.int=D.int,Y=Y,Z.formula=Z.formula,L.values=L.values,hazards=hazards,data=data)))
}

#Returns the derivative of the Estimating equation for every patient integrated over all time points
U.deriv.hazard = function(p,dN,D,D.int,Y,Z.formula,L.values,hazards,data){
  integral = 0
  t.min = 0
  Z = data$Z
  for(t in colnames(Y)){
    dN.t = dN[,t]
    D.t = D[,t]
    D.int.t = D.int[,t]
    Y.t = Y[,t]
    lambda.t = hazards[,t]
    Z.term.t = Z.term.strat(beta=p,Z=data$Z,Y.t=Y.t,D.int.t=D.int.t,Z.formula=Z.formula,L.values=L.values,data=data)
    integral = integral + U.deriv.hazard.t(p=p,dN.t=dN.t,D.t=D.t,D.int.t=D.int.t,Y.t=Y.t,lambda.t=lambda.t,Z.term.t=Z.term.t,Z=Z,t=t,t.min=t.min,Z.formula=Z.formula,L.values=L.values,data=data)
    t.min = t
  }
  return(integral)
}

#Returns the derivative of the Estimating equation for every patient at time t
U.deriv.hazard.t = function(p,dN.t,D.t,D.int.t,Y.t,lambda.t,Z.term.t,Z,t,t.min,Z.formula,L.values,data){
  beta.est = p
  exp.term = exp(beta.est*D.int.t)
  deriv.Z.term.t = deriv.Z.term.t.strat(p=p,dN.t=dN.t,D.t=D.t,D.int.t=D.int.t,Y.t=Y.t,Z.term.t=Z.term.t,Z=Z,t=t,t.min=t.min,Z.formula=Z.formula,L.values=L.values,data=data)
  diff = (as.numeric(t)-as.numeric(t.min))
  output = -deriv.Z.term.t*exp.term*(dN.t-(Y.t*D.t*beta.est+Y.t*lambda.t)*diff)+Z.term.t*exp.term*(D.int.t*(dN.t-(Y.t*D.t*beta.est+Y.t*lambda.t)*diff)-Y.t*D.t*diff)
  return(output)
}


#Function that returns an estimate for beta
av.deriv.function.hazard = function(fit.ivsacim.by_prod,data,beta.start,stime,Z.formula,L.values,hazards){
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  output = f.deriv.hazard(p=beta.start,dN = dN, Y=Y,D=D,D.int=D.int,Z.formula=Z.formula,L.values=L.values,hazards=hazards,data=data)
  return(output)
}

p.score.hazard = function(beta=0,fit.ivsacim.by_prod,data,stime,Z.formula,L.values,hazards){
  dN = fit.ivsacim.by_prod$dN
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  U.vector = U.hazard(p=beta,dN = dN, Y=Y,D=D,D.int=D.int,data=data,Z.formula=Z.formula,L.values=L.values,hazards=hazards)
  return(t.test(U.vector)$p.value)
}

#Returns the SE for beta, using the sandwich estimator 
beta.se.function.hazard = function(fit.ivsacim.by_prod,data,beta.est,Z.formula,L.values,hazards){
  dN = fit.ivsacim.by_prod$dN 
  colnames(dN) = as.character(stime) 
  Y = fit.ivsacim.by_prod$risk_t
  colnames(Y) = as.character(stime) 
  D = fit.ivsacim.by_prod$D_status
  colnames(D) = as.character(stime)
  D.int = t(apply(D,MARGIN = 1,FUN = function(x) x*c(stime[1],diff(stime))))
  D.int = t(apply(D.int, MARGIN = 1, FUN = cumsum )) #int 0^t- D(s)
  
  U.vector = U.hazard(p=beta.est,dN = dN, Y=Y,D=D,D.int=D.int,data=data,Z.formula=Z.formula,L.values=L.values,hazards=hazards)
  var.U = var(U.vector)
  deriv.U.vector = U.deriv(p=beta.est,dN = dN, Y=Y,D=D,D.int=D.int,data=data,Z.formula=Z.formula,L.values=L.values)
  var.beta = var.U/(dim(data)[1]*(mean(deriv.U.vector))^2)
  return(sqrt(var.beta))
}

f.score.hazard = function(p,fit.ivsacim.by_prod,data,stime,Z.formula,L.values,p.value=0.05,hazards=hazards){
  beta = p
  return((p.score.hazard(beta=beta,fit.ivsacim.by_prod=fit.ivsacim.by_prod,data=data,stime=stime,Z.formula=Z.formula,L.values=L.values,hazards=hazards)-p.value))
}

beta.bound.hazard= function(lower,upper,fit.ivsacim.by_prod,data,stime,Z.formula,L.values,p.value=0.05,hazards=hazards){
  beta.bound = uniroot(f=f.score.hazard,lower=lower,upper=upper,fit.ivsacim.by_prod=fit.ivsacim.by_prod,data=data,stime=stime,p.value=p.value,Z.formula=Z.formula,L.values=L.values,hazards = hazards)$root
  return(beta.bound)
}


###Hypothetical curve ######

#Z.formula

hyp.surv.aalen = function(beta,times,data,T.formula){
  #Model for survival times
  fit.Aalen <- aalen(T.formula, data = subset(data,Z==0),resample.iid=1)
  pred.Aalen = predict.timereg(fit.Aalen,resample.iid=1,newdata=data,times=times)
  survival.pred = pred.Aalen$S0 #gives survival predictions for survival models. Predictions given in matrix form with different subjects in different rows.
  #Multiply every column by exp(-beta*t)
  survival.pred = survival.pred %*% diag(exp(-beta*times)) #In every column: exp(-beta*t)*P(T>t given Z=0,L)
  colnames(survival.pred) = as.character(times)
  #Take the avarage
  output = apply(survival.pred, MARGIN = 2, FUN = mean) 
  return(output)
}


hyp.surv.cox = function(beta,times,data,T.formula){
  fit.Cox <- coxph(T.formula,data = subset(data,Z==0),x=TRUE,y=TRUE)
  survival.pred = predictSurvProb(fit.Cox,newdata=data,times=times)
  colnames(survival.pred) = as.character(stime) #S(t given L, Z=0)
  #Multiply every column by exp(-beta*t)
  survival.pred = survival.pred %*% diag(exp(-beta*times)) #In every column: exp(-beta*t)*P(T>t given Z=0,L)
  colnames(survival.pred) = as.character(times)
  #Take the avarage
  output = apply(survival.pred, MARGIN = 2, FUN = mean) 
  return(output)
}  


surv.times = function(data,times){
  fit= survfit(Surv(time,status) ~ 1, data = data) #Kaplan Meier survival curve
  surv = summary(fit,time=times, extend = TRUE)$surv #Survival probabilities 
  return(surv)
}

library(reshape2)

plot.hyp.surv= function(overview){
  data.plot = overview
    plot.data <- melt(data.plot, measure.vars=c("y.control","y.exp","hypo.surv"), id.vars = c("time"), variable.name = "arm.y")
    
    plot.data$arm = ifelse(plot.data$arm.y=="y.control",yes="Observed control arm",no=plot.data$arm)
    plot.data$arm = ifelse(plot.data$arm.y=="y.exp",yes="Observed treatment arm",no=plot.data$arm)
    plot.data$arm = ifelse(plot.data$arm.y%in% c("hypo.surv.aalen","hypo.surv.cox","hypo.surv"),yes="Always control (without crossover)",no=plot.data$arm)
    plot.data$Observed = ifelse(plot.data$arm %in% c("Observed overall","Observed control arm","Observed treatment arm"), yes="Observed", no="Predicted")
    
    plot.data$arm <- factor(plot.data$arm, levels=c("Observed treatment arm", "Observed control arm", "Always control (without crossover)"))
    p = ggplot(data=plot.data,aes(x=time,y=value, color=arm,linetype=Observed))+ geom_line(size=1) #+geom_point(position=pd)
    p = p +ggtitle("Observed and hypothetical survival curves") 
    p = p + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.box="vertical", legend.margin=margin())
    p = p + xlab("Follow-up time") + ylab("Overall survival")
   p = p + labs(color='Arm') 
  return(p)
}

