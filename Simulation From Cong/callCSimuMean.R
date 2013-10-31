
#=============== Simulation Study to examine the R version which calls C: Mean Estimation ===============#
library(KernSmooth)
library(lpridge)

dyn.load('pace.so')
source("pace.R")
load("DataS.RData")

Nsimu = 100

#=============== First, run our code with binning ===============#

resultBin = list()
timeBin = timeSmooth = rep(0,Nsimu)
x_out = seq(0,10,by=0.1)

for(i in 1:Nsimu){
  x <- unlist(DataS[[i]]$x)
  y <- unlist(DataS[[i]]$y)
  timeBin[i] <- system.time(newData<-binning(x,y,num_bins=100))[3]
  x_in <- newData$newx
  y_in <- newData$newy
  w_in <- newData$count
  timeSmooth[i] <- system.time(resultBin[[i]]<-lwls(bandwidth=1,kernel=0,x_in,y_in,w_in,
                                                    x_out,power=1,cv_mode=0)$output)[3]
}
cat("Finished the first step.\n")


#========== Second, run our code without binning and compare with locpoly, lpridge ==========#

time1 = time2 = time3 = time4 = rep(0,Nsimu)
results1 = results2 = results3 = results4 = list()
# Note: 1<-our code with Gaussian Kernel; 2<- our code with Epan Kernel; 3<-locpoly; 4<-lpridge #
x_out = seq(0,10,by=0.1)

for(i in 1:Nsimu){
  print(i)
  x <- unlist(DataS[[i]]$x)
  y <- unlist(DataS[[i]]$y)
  w <- rep(1/length(x),length(x))
  time2[i] <- system.time(results2[[i]]<-lwls(bandwidth=1,kernel=0,x,y,w,x_out,power=1,cv_mode=0)$output)[3]
  time4[i] <- system.time(results4[[i]]<-lpridge(x=x,y=y,bandwidth=1,x.out=x_out,ridge=0)$est)[3]
  
  time3[i] <- system.time(results3[[i]]<-locpoly(x=x,y=y,bandwidth=1,gridsize=101))[3]
  time1[i] <- system.time(results1[[i]]<-lwls(bandwidth=1,kernel=2,x,y,w,results3[[i]]$x,power=1,cv_mode=0)$output)[3]
}

out <- list(timeBin=timeBin,timeSmooth=timeSmooth,resultBin=resultBin,time1=time1,time2=time2,time3=time3,
            time4=time4,results1=results1,results2=results2,results3=results3,results4=results4)

save(out,file="SimuResultsMean.RData") 
# correspond to "ridge=0" for lpridge, "f=0.1" for lowess and same output points for locpoly and our procedure #




