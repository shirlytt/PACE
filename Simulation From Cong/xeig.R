
#=============== Create eigenfunctions in the domain of [0,T] ===============#

#===== function phi = xeig(t,T,numPhi) =====#

#===========
#  Input:
#===========
# input t      : 1*m vector of time points for each of the eigenfunctions
# input T      : a real value number where T > 0, the data is to be 
#                sampled in the domain of [0,T] 
# input numPhi : positive integer for the number of eigenfunctions to be returned
#                default is set to be 2


#==========
#  Output:
#========== 
# output phi   : numPhi * m matrix of eigenfunctions with phi[k,] being the kth
#                eigenfunction

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





