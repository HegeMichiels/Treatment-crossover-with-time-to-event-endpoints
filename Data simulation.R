
#Data simulation ###

n = 500
Z = rbinom(n, 1, 0.5) #randomized arm: 1 control, 0 experimental treatment
gender <- rbinom(n,size=1,prob=0.5)
age <- rnorm(n,mean=50,sd = 7)


lambda=0.002
rho=1
C=60

#Crossover times 
V = runif(n=n)
Tshift = ifelse(Z==1, yes = (- log(V) / (3*lambda*exp(Z*0.35 + gender*0.4+age*0.01)))^(1 / rho), no=0) #only crossover for control patients 
hist(Tshift)
summary(Tshift)
Tshift = ifelse(Tshift<C,yes=Tshift,no=0) #if crossover happens after the last day of the trial, its removed
summary(Tshift)

#Survival times
V = runif(n=n)
Tsurv =  (- log(V) / (0.7*lambda * exp(Z*0.3 + gender*0.3+age*0.03)))^(1 / rho)
hist(Tsurv)
summary(Tsurv)
Tsurv = ifelse(Tshift ==0, yes= Tsurv,no = ifelse(Tshift<Tsurv, yes = Tsurv*runif(n,min=1,max=1.15),no = Tshift) )
hist(Tsurv)
summary(Tsurv)
Tshift = ifelse(Tshift<Tsurv,yes=Tshift,no=0)

#Admnistrative censoring
time <- pmin(Tsurv, C)
status <- as.numeric(Tsurv <= C)
data <-  data.frame(id=1:n, gender=gender,age=age,Z=Z,time=time, status=status,shift=Tshift)

hist(data$time)

setwd("~/OneDrive - UGent/Crossover/R-code")
write.table(data, file = "data.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
