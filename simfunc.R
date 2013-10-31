xeig <- function(t,T,numPhi = 2){
  
  if(is.na(numPhi)) numPhi <- 2
  id <- which(t < 0 | t > T)
  if(length(id) > 0){
    warning("t must be in [0,T]. Invalid t's are removed!\n")
    t <- t[-id]
  } 
  
  if(numPhi == 1) phi <- -sqrt(2/T)*cos(2*pi*t/T)
  if(numPhi == 2){
    phi <- matrix(0,2,length(t))
    phi[1,] <- -sqrt(2/T)*cos(2*pi*t/T)
    phi[2,] <- sqrt(2/T)*sin(2*pi*t/T)
  }
  if(numPhi > 2){
    phi <- matrix(0,numPhi,length(t))
    id <- 1:numPhi
    oddID <- id%%2
    oddFactor <- 1:sum(oddID)
    evenID <- as.numeric(oddID == 0)
    evenFactor <- 1:sum(evenID)
    phiOdd <- matrix(0,sum(oddID),length(t))
    phiEven <- matrix(0,sum(evenID),length(t))
    
    for(i in 1:sum(oddID)){
      phiOdd[i,] <- -sqrt(2/T)*cos(2*oddFactor[i]*pi*t/T)
    }
    for(i in 1:sum(evenID)){
      phiEven[i,] <- sqrt(2/T)*sin(2*evenFactor[i]*pi*t/T)
    }
    phi[which(oddID == 1),] = phiOdd
    phi[which(evenID == 1),] = phiEven
  }
  return(phi)
}

#==================== Generate Simulation Data ====================#

#===== setting = 1 means "dense"; setting = 2 means "dense with missing"; =====%
#===== setting = 3 means "dense but irregular"; setting = 4 means "sparse" =====%
#===== shape=1 means the relatively smooth case; shape=2 means the case with a sharp peak =====%

SimuData <- function(n,setting,shape,error,CK){
  
  x <- y <- list()
  scores <- matrix(0,n,CK)
  
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
      merr <- rnorm(51,0,error) # measurement error #
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
      merr <- rnorm(51,0,error) # measurement error #
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
      merr <- rnorm(num,0,error) # measurement error #
      y[[i]] <- mu+gvalue+merr
    }
    
    if(setting == 4){
      num <- sample(4,1) + 1# number of measurements #
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
      merr <- rnorm(num,0,error) # measurement error #
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

