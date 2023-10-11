library(lamW)
library(deSolve)
$
N <- 100000
I_max <- 4172
#E_0 <- 7
I_0 <- 83 #412    
R_0 <- 69 #345    #### look up in the simulation study design document Table X
gamma_SIR <- (0.05/5.1)/(0.05 + 1/5.1) ### possibly translating the gamma from SEIR model into gamma in SIR model using SEIR model's revocery rate 0.05 and latent period 5.1 days
beta_obj <- -(gamma_SIR/(1-I_max/N)) * lambertWm1(-((1-I_max/N)/(1-I_0/N))*exp(-1))

##objective simulation
parameters <- c(beta_para=beta_obj, gamma_para=gamma_SIR)
states <- c(S=100000-I_0-R_0, I=I_0, R=R_0)

SIR <- function(t, states, parameters) {
  with(
    as.list(c(states, parameters)),
    {
      dS <- -beta_para*S*I/N
      dI <- beta_para*S*I/N - gamma_para*I
      dR <- gamma_para*I
      
      list(c(dS, dI, dR))
    }
  )
}

times <- 1:730

out <- ode(y = states, times = times, func = SIR, parms = parameters)
out[,3]

plot(out[,2],type="l",ylim=c(0,100000),lwd=2,col="blue",xlab="Days forward (from May 15)",ylab="Number of people",cex.lab=1.4)
par(new=T)
plot(out[,3],type="l",ylim=c(0,100000),lwd=2,col="red",yaxt="n",xaxt="n",ylab="",xlab="")
par(new=T)
plot(out[,4],type="l",ylim=c(0,100000),lwd=2,col="green",yaxt="n",xaxt="n",ylab="",xlab="")
abline(h=I_max,lty=2,col="red")
legend(0,85000,lty=c(1,1,1,2),col=c("blue","red","green","red"), legend=c("Susceptible", "Infectious", "Removed", "I_max"))


###evaluation: data are taken from the simulation study design document Table X, 
###Training-5, where modeler assumes that 25% of the cases were not reported.
###The E compartment in SEIR model is added to the SIR model's S compartment
S_current <- 99666 + 81
#E_current <- 81
I_current <- 138
R_current <- 115

beta_t <- c()
S_t <- NULL
#E_t <- NULL
I_t <- NULL
R_t <- NULL
for (ti in 1:730) {
  controller <- out[ti,3]/N - I_current/N - (out[ti,2]/N - (S_current)/N) + beta_obj*out[ti,3]*out[ti,2]/(I_current*(S_current))  ####
  #controller_SEIR <- controller*0.05/gamma_SIR
  
  beta_t <- c(beta_t, controller)
  S_t <- c(S_t, S_current)
  #E_t <- c(E_t, E_current)
  I_t <- c(I_t, I_current)
  R_t <- c(R_t, R_current)
  
  parameters <- c(beta_para=unname(controller), gamma_para=gamma_SIR)
  states <- c(S=unname(S_current), I=unname(I_current), R=unname(R_current))
  
  
  times <- 1:2
  out1 <- ode(y = states, times = times, func = SIR, parms = parameters)
  
  S_current <- out1[2,2]
  #E_current <- out1[2,3]
  I_current <- out1[2,3]
  R_current <- out1[2,4]
}

dataset_baseline <- data.frame(out,beta_t,S_t,I_t,R_t)
write.csv(dataset_baseline, "/Users/gracejia/Documents/A-UW/covid19 NSF Project/baseline.csv", row.names = F)

plot(beta_t, type="l", xlab="Days forward (from May 15)",ylab="",cex.lab=1.4)
mtext("Recommended transmission rate\n(i.e. the controller)",2.1, line=2,cex=1.4)

plot(S_t,type="l",ylim=c(0,100000),lwd=2,col="blue",xlab="Days forward (from May 15)",ylab="Number of people",cex.lab=1.4)
par(new=T)
plot(I_t,type="l",ylim=c(0,100000),lwd=2,col="red",yaxt="n",xaxt="n",ylab="",xlab="")
par(new=T)
plot(R_t,type="l",ylim=c(0,100000),lwd=2,col="green",yaxt="n",xaxt="n",ylab="",xlab="")
abline(h=I_max,lty=2,col="red")
legend(0,85000,lty=c(1,1,1,1,2),col=c("blue","orange","red","green","red"), legend=c("Susceptible", "Exposed","Infectious", "Removed", "I_max"))

plot(out[,3],type="l",ylim=c(0,4200),lwd=2,col="red",ylab="Number of infectious people",xlab="Days forward",cex.lab=1.4,lty=2)
par(new=T)
plot(I_t,type="l",ylim=c(0,4200),lwd=2,col="red",yaxt="n",xaxt="n",ylab="",xlab="")
abline(h=I_max,lty=2,col="red",lwd=1)


beta_t_baseline <- beta_t
I_t_baseline <- I_t

###Evaluation with data issues
S_current <- 99666 + 81
#E_current <- 81
I_current <- 138
R_current <- 115

beta_t <- NULL
S_t <- NULL
#E_t <- NULL
I_t <- NULL
R_t <- NULL
for (ti in 1:730) {
  false_I <- I_current*0.15/0.75    ####
  false_S_E <- 100000 - false_I - R_current
  
  controller <- out[ti,3]/N - false_I/N - (out[ti,2]/N - false_S_E/N) + beta_obj*out[ti,3]*out[ti,2]/(false_I*false_S_E)
  
  if (controller > 0.196) controller <- 0.196
  
  beta_t <- c(beta_t, controller)
  S_t <- c(S_t, S_current)
  I_t <- c(I_t, I_current)
  R_t <- c(R_t, R_current)
  
  parameters <- c(beta_para=unname(controller), gamma_para=gamma_SIR)
  states <- c(S=unname(S_current), I=unname(I_current), R=unname(R_current))
  

  times <- 1:2
  out1 <- ode(y = states, times = times, func = SIR, parms = parameters)
  
  S_current <- out1[2,2]
  I_current <- out1[2,3]
  R_current <- out1[2,4]
}

dataset_scenario <- data.frame(out,beta_t,S_t,I_t,R_t)
write.csv(dataset_scenario, "/Users/gracejia/Documents/A-UW/covid19 NSF Project/scenario2_025.csv", row.names = F)
plot(beta_t, type="l", xlab="Days forward (from May 15)",ylab="",cex.lab=1.4)
par(new=T)
plot(beta_t_baseline, type="l",lty=2, xlab="",ylab="",cex.lab=1.4,ylim=c(0.057,0.086))
mtext("Allowable transmission rate",2.1, line=2.5,cex=1.4)

plot(S_t,type="l",ylim=c(0,100000),lwd=2,col="blue",xlab="Days forward (from May 15)",ylab="Number of people",cex.lab=1.4)
par(new=T)
plot(I_t,type="l",ylim=c(0,100000),lwd=2,col="red",yaxt="n",xaxt="n",ylab="",xlab="")
par(new=T)
plot(R_t,type="l",ylim=c(0,100000),lwd=2,col="green",yaxt="n",xaxt="n",ylab="",xlab="")
abline(h=I_max,lty=2,col="red")
legend(0,85000,lty=c(1,1,1,1,2),col=c("blue","orange","red","green","red"), legend=c("Susceptible", "Exposed","Infectious", "Removed", "I_max"))

#plot(I_t_baseline,type="l",ylim=c(0,4200),lwd=2,col="red",ylab="Number of infectious people",xlab="Days forward",cex.lab=1.4,lty=2)
plot(out[,3],type="l",ylim=c(0,4200),lwd=2,col="red",ylab="Number of infectious people",xlab="Days forward",cex.lab=1.4,lty=2)
par(new=T)
plot(I_t,type="l",ylim=c(0,4200),lwd=2,col="red",yaxt="n",xaxt="n",ylab="",xlab="")
par(new=T)
plot(I_t*0.15/0.25,type="l",ylim=c(0,4200),lwd=2,col="black",yaxt="n",xaxt="n",ylab="",xlab="")
abline(h=I_max,lty=2,col="red",lwd=1)

