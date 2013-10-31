
#==================== Extract the simulation results: Second derivative ====================#

#========== Note: 1<-our code with Gaussian Kernel; 2<- our code with Epan Kernel; ==========#
#==========       3<-locpoly (Gaussian Kernel); 4<-lpridge (Epan Kernel); =========#


load("SimuResultsDeriv2.RData")
time1 <- out$time1; time2 <- out$time2; time3 <- out$time3; time4 <- out$time4
timeBin <- out$timeBin; timeSmooth <- out$timeSmooth
sum(time1); sum(time2); sum(time3); sum(time4); sum(timeBin+timeSmooth)
[1] 0.153
[1] 0.018
[1] 0.037
[1] 0.051
[1] 0.227
# our code with Epan Kernel is really fast #

results1 <- out$results1; results2 <- out$results2
results3 <- out$results3; results4 <- out$results4
resultBin <- out$resultBin

#========== First compare our algorithm with lpridge (both with Epan Kernel) ==========#
D2dnorm <- function(x,mean,sd){
  temp1 = dnorm(x,mean,sd)
  temp2 = ((x-mean)/(sd^2))^2-1/(sd^2)
  out = temp1*temp2
  return(out)
}

Nsimu = 100
x_out = seq(0,10,by=0.1)
true = -sin(x_out)+30*D2dnorm(x_out,mean = 5,sd = 2/3)
ISE2 = ISE4 = ISEbin = rep(0,Nsimu)

for(i in 1:Nsimu){
  ISE2[i] <- sum((results2[[i]]-true)^2)
  ISE4[i] <- sum((results4[[i]]-true)^2)
  ISEbin[i] <- sum((resultBin[[i]]-true)^2)
}

mean(ISE2); mean(ISE4); mean(ISEbin)
[1] 11947.27
[1] 11947.27
[1] 14864.69
#===== Exactly the same result! Our code is slower! =====#



#========== Second compare our algorithm with locpoly (both with Gaussian Kernel) ==========#
ISE1 = ISE3 = rep(0,Nsimu)

for(i in 1:Nsimu){
  x_out <- results3[[i]]$x
  vtrue <- -sin(x_out)+30*D2dnorm(x_out,mean = 5,sd = 2/3)
  ISE1[i] <- sum((results1[[i]]-vtrue)^2)
  ISE3[i] <- sum((results3[[i]]$y-vtrue)^2)
}

mean(ISE1); mean(ISE3)
[1] 10865.27
[1] 10860.29
#===== Our code performs slightly worse and also slower (with the Gaussian Kernel) =====#
