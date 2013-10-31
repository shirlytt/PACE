
#==================== Generate Simulation Data ====================#

#===== setting = 1 means "dense"; setting = 2 means "dense with missing"; =====%
#===== setting = 3 means "dense but irregular"; setting = 4 means "sparse" =====%
#===== shape=1 means the relatively smooth case; shape=2 means the case with a sharp peak =====%

SimuData <- function(n,setting,shape){
  
  x <- y <- list()
  scores <- matrix(0,n,10)
  
  for(i in 1:n){
    if(setting == 1){
      x[[i]] <- seq(0,10,length.out=51)
      mu <- x[[i]]+sin(x[[i]])
      k <- 1:CK
      const <- 9*((2/3)^(k-1))
      stemp <- rnorm(CK) # standard normal r.v. #
      score <- stemp*sqrt(const)
      scores[i,] <- score
      gtemp <- xeig(x[[i]],10,CK) # 10*51 matrix #
      gvalue <- as.vector(score%*%gtemp)
      merr <- rnorm(51) # measurement error #
      y[[i]] <- mu+gvalue+merr
    }
    
    if(setting == 2){
      x[[i]] <- seq(0,10,length.out=51)
      mu <- x[[i]]+sin(x[[i]])
      k <- 1:CK
      const <- 9*((2/3)^(k-1))
      stemp <- rnorm(CK) # standard normal r.v. #
      score <- stemp*sqrt(const)
      scores[i,] <- score
      gtemp <- xeig(x[[i]],10,CK) # 10*51 matrix #
      gvalue <- as.vector(score%*%gtemp)
      merr <- rnorm(51) # measurement error #
      y[[i]] <- mu+gvalue+merr
      missN <- min(rbinom(1,size=51,prob=0.1),16) # number of missing measurements #
      Ind <- sample(51,missN) # index of the non-missing measurements #
      # easy to use in data frame form
      y[[i]][sort(Ind)] = NA
    }
    
    if(setting == 3){
      num <- 34+sample(16,1) # number of measurements #
      tempT <- sort(10*runif(num))
      x[[i]] <- unique(tempT)
      num <- length(x[[i]])
      mu <- x[[i]]+sin(x[[i]])
      k <- 1:CK
      const <- 9*((2/3)^(k-1))
      stemp <- rnorm(CK) # standard normal r.v. #
      score <- stemp*sqrt(const)
      scores[i,] <- score
      gtemp <- xeig(x[[i]],10,CK) # 10*num matrix #
      gvalue <- as.vector(score%*%gtemp)
      merr <- rnorm(num) # measurement error #
      y[[i]] <- mu+gvalue+merr
    }
    
    if(setting == 4){
      num <- sample(4,1) # number of measurements #
      tempT <- sort(10*runif(num))
      x[[i]] <- unique(tempT)
      num <- length(x[[i]])
      mu <- x[[i]]+sin(x[[i]])
      k <- 1:CK
      const <- 9*((2/3)^(k-1))
      stemp <- rnorm(CK) # standard normal r.v. #
      score <- stemp*sqrt(const)
      scores[i,] <- score
      gtemp <- xeig(x[[i]],10,CK) # 10*num matrix #
      gvalue <- as.vector(score%*%gtemp)
      merr <- rnorm(num) # measurement error #
      y[[i]] <- mu+gvalue+merr
    }
  }
  
  if(shape == 2){
    for(i in 1:n){
      fvalue <- 30*dnorm(x[[i]],mean = 5,sd = 2/3)
      y[[i]] <- y[[i]]+fvalue
    }
  }
  
  data <- list(x=x,y=y,scores=scores)
  return(data)
}



#==================== Generate Simulation Data Set ====================#

source("xeig.R")

n = 200
Nsimu = 100
CK = 2
DataS = list()

# for(i in 1:Nsimu){
# DataS[[i]] = SimuData(n,1,1)
# }
# save(DataS,file="DataS1s.RData")
# for(i in 1:Nsimu){
# DataS[[i]] = SimuData(n,1,2)
# }
# save(DataS,file="DataS1ns.RData")
# for(i in 1:Nsimu){
# DataS[[i]] = SimuData(n,2,1)
# }
# save(DataS,file="DataS2s.RData")
# for(i in 1:Nsimu){
# DataS[[i]] = SimuData(n,2,2)
# }
# save(DataS,file="DataS2ns.RData")

for(i in 1:Nsimu){
DataS[[i]] = SimuData(n,3,1)
}
save(DataS,file="DataS3s.RData")

for(i in 1:Nsimu){
DataS[[i]] = SimuData(n,3,2)
}
save(DataS,file="DataS3ns.RData")

for(i in 1:Nsimu){
DataS[[i]] = SimuData(n,4,1)
}
save(DataS,file="DataS4s.RData")

for(i in 1:Nsimu){
DataS[[i]] = SimuData(n,4,2)
}
save(DataS,file="DataS4ns.RData")
