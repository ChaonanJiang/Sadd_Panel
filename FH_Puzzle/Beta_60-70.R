library("R.matlab")
library("spdep") 
library("splm") 

FHD <- readMat("C:/Users/jiangc/Dropbox/Sadd_Panel_Latex/Papers/applications/Scaillet_Suggestions/FHD Data Debarsy Ertur/FHD.mat")
Wn <-  readMat("C:/Users/jiangc/Dropbox/Sadd_Panel_Latex/Papers/applications/Scaillet_Suggestions/FHD Data Debarsy Ertur/matrices.mat")

Yn <- FHD$data[1:264,3]  #### Investment rate
Xn <- FHD$data[1:264,5]  #### Saving rate


Wdia1 <- Wn$Wdia1
W <- Wn$Wdia1
Wn.a <- as.matrix(W)
Mn.a <- Wn.a
n <- 24
T <- 11
m <- n*(T-1)
MCsize <- 10000
lw <- mat2listw(Wn.a)
dt <- data.frame(id = rep(c(1:n),T), time =kronecker(1960:1970,rep(1,24)),Yn,Xn)

Y.tilde.nt<-matrix(Yn,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn,nrow=n,ncol=T)),T),n,T)
X.tilde.nt<-matrix(Xn,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Xn,nrow=n,ncol=T)),T),n,T)

sarar <- spml(formula = Yn~Xn, data = dt , listw = mat2listw(W), model = "within", spatial.error= "b",lag = T, LeeYu = T, Hess = F)
summary(sarar)
beta <- sarar$coefficients[3]
rho <- sarar$coefficients[2]
sig2 <- sarar$sigma2
lambda <- sarar$coefficients[1]
beta0 <- 0

###### log likelihood function 

log.lik<-function(beta,lambda,rho,sig2) { 
  
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta) 
  
  ell_nt<- -n*(T-1)/2*log(2*pi*sig2)+(T-1)*(log(det(Sn))+log(det(Rn)) ) - 0.5/sig2* sum(diag((t(V.tilde.nt) %*% (V.tilde.nt)))) 
  return(ell_nt) 
  
} 

###### The first derivative of log likelihood  

der1.log.lik <- function(beta,lambda,rho,sig2) { 
  
  Sn <- diag(n)-lambda*as.matrix(W) 
  Rn <- diag(n)-rho*Mn.a 
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  
  deriv1 <- c((1/sig2)*sum(diag(crossprod((Rn%*%X.tilde.nt),V.tilde.nt))),sum(diag(t(Yw)%*%V.tilde.nt))/sig2-(T-1)*sum(diag(Gn)),
              sum(diag(crossprod((Hn%*%V.tilde.nt),V.tilde.nt)))/sig2-(T-1)*sum(diag(Hn)),-n*(T-1)/(2*sig2)+1/(2*sig2*sig2)* sum(diag((t(V.tilde.nt) %*% (V.tilde.nt)))) )
  
  return(deriv1) 
  
} 

###### Score function S_n,t (containing the whole cross-sections) 

Score.nt <- function(t,beta,lambda,rho,sig2){ 
  
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  
  Scor.nt <- c(n*(T-1)*crossprod(Rn%*%X.tilde.nt[,t],V.tilde.nt[,t])/sig2, n*(T-1)*crossprod(Yw[,t],V.tilde.nt[,t])/sig2-n*(T-1)*(T-1)/T*sum(diag(Gn)),
               n*(T-1)*crossprod(Hn%*%V.tilde.nt[,t],V.tilde.nt[,t])/sig2-n*(T-1)*(T-1)/T*sum(diag(Hn)),n*(T-1)/(2*sig2*sig2)*(crossprod(V.tilde.nt[,t])-n*(T-1)*sig2/T))
  
  return(Scor.nt) 
} 


###### Score function S_i,t (for each individual i at time t) 

Score.it <- function(i,t,beta,lambda,rho,sig2) { 
  
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  
  Scor.it <- c(n*(T-1)*(Rn%*%X.tilde.nt)[i,t]*V.tilde.nt[i,t]/sig2, n*(T-1)*Yw[i,t]*V.tilde.nt[i,t]/sig2 - n*((T-1)*(T-1))/T*Gn[i,i],
               n*(T-1)*(Hn%*%V.tilde.nt)[i,t]*V.tilde.nt[i,t]/sig2 - n*((T-1)*(T-1))/T*Hn[i,i], n*(T-1)/(2*sig2*sig2)*(V.tilde.nt[i,t]*V.tilde.nt[i,t]-(T-1)*sig2/T))
  
  
  return(Scor.it) 
  
}


###### sum_t S_i,t (sum over time dimension) 

Sum.T.Scor <- function(i,beta,lambda,rho,sig2){ 
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  
  Score.iT <- c(n*(T-1)*((Rn%*%X.tilde.nt)%*%t(V.tilde.nt))[i,i]/sig2, n*(T-1)*(Yw%*%t(V.tilde.nt))[i,i]/sig2 - n*((T-1)*(T-1))*Gn[i,i],
                n*(T-1)*((Hn%*%V.tilde.nt)%*%t(V.tilde.nt))[i,i]/sig2 - n*((T-1)*(T-1))*Hn[i,i], n*(T-1)/(2*sig2*sig2)*((V.tilde.nt%*%t(V.tilde.nt))[i,i]-(T-1)*sig2))
  
  
  return(Score.iT) 
  
  
} 

###### The first derivative of S_i,t Negative

der1.Score.it <- function(i,t,beta,lambda,rho,sig2){ 
  
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn 
  
  a <- n*(T-1)*(Rn%*%X.tilde.nt)[i,t]*(Rn%*%X.tilde.nt)[i,t]/sig2  ###[1,1]
  b <- n*(T-1)*Yw[i,t]*(Rn%*%X.tilde.nt)[i,t]/sig2                 ###[2,1]
  c <- n*(T-1)*((Hn%*%V.tilde.nt)[i,t]*(Rn%*%X.tilde.nt)[i,t]+V.tilde.nt[i,t]*(Mn.a%*%X.tilde.nt)[i,t])/sig2  ###[3,1]
  d <- n*(T-1)*(Rn%*%X.tilde.nt)[i,t]*V.tilde.nt[i,t]/(sig2*sig2)   ###[4,1]
  e <- n*(T-1)*Yw[i,t]*Yw[i,t]/sig2 + n*((T-1)*(T-1))/T*G2n[i,i] ###[2,2]
  f <- n*(T-1)*(Hn%*%V.tilde.nt)[i,t]*Yw[i,t]/sig2+n*(T-1)*(Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t]*V.tilde.nt[i,t]/sig2 ###[2,3]
  g <- n*(T-1)*Yw[i,t]*V.tilde.nt[i,t]/(sig2*sig2)                 ###[2,4]
  h <- n*(T-1)*(Hn%*%V.tilde.nt)[i,t]*(Hn%*%V.tilde.nt)[i,t]/sig2 + n*((T-1)*(T-1))/T*H2n[i,i]   ###[3,3]
  k <- n*(T-1)*(Hn%*%V.tilde.nt)[i,t]*V.tilde.nt[i,t]/(sig2*sig2)  ###[3,4]
  j <- n*(T-1)*(V.tilde.nt[i,t]*V.tilde.nt[i,t])/(sig2*sig2*sig2)-n*(T-1)*(T-1)/(2*T*sig2*sig2)  ###[4,4]
  
  
  deriv.Scor.it <- matrix(c(a,b,c,d,b,e,f,g,c,f,h,k,d,g,k,j),4,4)  
  return(deriv.Scor.it) 
  
} 


######sum_t the first derivative of S_it 

Sum.T.Der1 <- function(i,beta,lambda,rho,sig2){ 
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  
  a <- -n*(T-1)*((Rn%*%X.tilde.nt)%*%t(Rn%*%X.tilde.nt))[i,i]/sig2  ###[1,1]
  b <- -n*(T-1)*(Yw%*%t(Rn%*%X.tilde.nt))[i,i]/sig2                 ###[2,1]
  c <- -n*(T-1)*(((Hn%*%V.tilde.nt)%*%t(Rn%*%X.tilde.nt))[i,i]+(V.tilde.nt%*%t(Mn.a%*%X.tilde.nt))[i,i])/sig2  ###[3,1]
  d <- -n*(T-1)*((Rn%*%X.tilde.nt)%*%t(V.tilde.nt))[i,i]/(sig2*sig2)   ###[4,1]
  e <- -n*(T-1)*(Yw%*%t(Yw))[i,i]/sig2 - n*((T-1)*(T-1))*G2n[i,i] ###[2,2]
  f <- -n*(T-1)*((Hn%*%V.tilde.nt)%*%t(Yw))[i,i]/sig2-n*(T-1)*((Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(V.tilde.nt))[i,i]/sig2 ###[2,3]
  g <- -n*(T-1)*(Yw%*%t(V.tilde.nt))[i,i]/(sig2*sig2)                 ###[2,4]
  h <- -n*(T-1)*((Hn%*%V.tilde.nt)%*%t(Hn%*%V.tilde.nt))[i,i]/sig2 - n*((T-1)*(T-1))*H2n[i,i]   ###[3,3]
  k <- -n*(T-1)*((Hn%*%V.tilde.nt)%*%t(V.tilde.nt))[i,i]/(sig2*sig2)  ###[3,4]
  j <- -n*(T-1)*((V.tilde.nt%*%t(V.tilde.nt))[i,i])/(sig2*sig2*sig2)+n*(T-1)*(T-1)/(2*sig2*sig2)  ###[4,4]
  
  
  deriv.Scor.it <- matrix(c(a,b,c,d,b,e,f,g,c,f,h,k,d,g,k,j),4,4)  
  return(deriv.Scor.it) 
  
  
} 

###### M matrix 


M.matrix <- function(i,beta,lambda,rho,sig2){
  
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn 
  
  G.dot <- Rn%*%Gn%*%solve(Rn)
  GX <- Rn%*%X.tilde.nt
  
  
  
  a <- n*(GX)[i,]%*%(t((GX)))[,i]/sig2  ###[1,1]
  b <- n*beta*(G.dot%*%GX)[i,]*(t(GX))[,i]/sig2                 ###[2,1]
  c <- 0                                                          ###[3,1]
  d <- 0                                                          ###[4,1]
  e <- n*(T-1)*(G.dot%*%t(G.dot)+G2n)[i,i]+n*beta*beta*(G.dot%*%GX)[i,]%*%(t(G.dot%*%GX))[,i]/sig2###[2,2]
  f <- n*(T-1)*(Hn%*%(t(G.dot)+G.dot))[i,i]                           ###[2,3]
  g <- n*(T-1)*G.dot[i,i]/(sig2)                 ###[2,4]
  h <- n*(T-1)*(Hn%*%t(Hn)+ H2n)[i,i]   ###[3,3]
  k <- n*(T-1)*(Hn)[i,i]/(sig2)  ###[3,4]
  l <- n*(T-1)/(2*sig2*sig2)  ###[4,4]
  
  
  M.1 <- matrix(c(a,b,0,0,b,e,f,g,0,f,h,k,0,g,k,l),4,4)  
  
  return(M.1) 
} 

M.matrix1 <- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  s.1 <- matrix(rep(0,16),4,4)
  
  for(l in 1:MCsize){
    V.tilde.nt <- matrix(rnorm(n*T,mean = 0, sd = sqrt((T-1)*sig2/T)),n,T)
    Y.tilde.nt <- Sn.inv%*%(Rn.inv%*%V.tilde.nt+X.tilde.nt*beta)
    Yw <- Rn%*%Wn.a%*%Y.tilde.nt
    a <- sum(diag((Rn%*%X.tilde.nt)%*%t(Rn%*%X.tilde.nt)))/sig2  ###[1,1]
    b <- sum(diag((Yw%*%t(Rn%*%X.tilde.nt))))/sig2                 ###[2,1]
    # c <- sum(diag((Hn%*%V.tilde.nt)%*%t(Rn%*%X.tilde.nt)+V.tilde.nt%*%t(Mn.a%*%X.tilde.nt)))/sig2  ###[3,1]
    c <-0
    # d <- sum(diag((Rn%*%X.tilde.nt)%*%t(V.tilde.nt)))/(sig2*sig2)   ###[4,1]
    d <- 0
    e <- sum(diag(Yw%*%t(Yw)))/sig2 +((T-1))*sum(diag(G2n)) ###[2,2]
    f <- sum(diag((Hn%*%V.tilde.nt)%*%t(Yw)))/sig2+sum(diag((Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(V.tilde.nt)))/sig2 ###[2,3]
    g <- sum(diag(Yw%*%t(V.tilde.nt)))/(sig2*sig2)                 ###[2,4]
    h <- sum(diag((Hn%*%V.tilde.nt)%*%t(Hn%*%V.tilde.nt)))/sig2 + ((T-1))*sum(diag(H2n))   ###[3,3]
    k <- sum(diag((Hn%*%V.tilde.nt)%*%t(V.tilde.nt)))/(sig2*sig2)  ###[3,4]
    j <- sum(diag(V.tilde.nt%*%t(V.tilde.nt)))/(sig2*sig2*sig2)-n*(T-1)/(2*sig2*sig2)  ###[4,4]
    s <- matrix(c(a,b,c,d,b,e,f,g,c,f,h,k,d,g,k,j),4,4)
    
    
    
    
    s.1 <- s+s.1 
  }
  s.1 <- s.1/MCsize
  return(s.1)
}

M <- M.matrix1(beta0,lambda,rho,sig2)

# s <- matrix(rep(0,4*4),4,4)
# 
# for(i in 1:n){
#   s <- s+ M.matrix(i,beta,lambda0,rho,sig2)
# }
# M1 <- s/n

###### IF function 

IF.iT <- function(i,beta,lambda,rho,sig2){ 
  IF <- (1/(T-1))*solve(M)%*%Sum.T.Scor(i,beta,lambda,rho,sig2)
  return(IF) 
  
} 


###### The second derivative of S_it 

#####1st 

der2.Score.it.1 <- function(i,t,beta,lambda,rho,sig2){ 
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn 
  
  a <- 2*n*(T-1)*((Mn.a%*%X.tilde.nt)[i,t]*((Rn%*%X.tilde.nt)[i,t]))/sig2   ####[1,3]
  b <- n*(T-1)*((Rn%*%X.tilde.nt)[i,t])*((Rn%*%X.tilde.nt)[i,t])/(sig2*sig2)    ####[1,4]
  c <- n*(T-1)*((Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t]*(Rn%*%X.tilde.nt)[i,t]+(Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*(Mn.a%*%X.tilde.nt)[i,t])/sig2   ####[2,3]
  d <- n*(T-1)*(Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*(Rn%*%X.tilde.nt)[i,t]/(sig2*sig2)   ####[2,4]
  e <- 2*n*(T-1)*(Hn%*%V.tilde.nt)[i,t]*(Mn.a%*%X.tilde.nt)[i,t]/sig2
  f <- n*(T-1)*((Hn%*%V.tilde.nt)[i,t]*(Rn%*%X.tilde.nt)[i,t]+V.tilde.nt[i,t]*(Mn.a%*%X.tilde.nt)[i,t])/(sig2*sig2)
  g <- 2*n*(T-1)*V.tilde.nt[i,t]*((Rn%*%X.tilde.nt)[i,t])/(sig2*sig2*sig2)
  
  Der2 <- matrix(c(0,0,a,b,0,0,c,d,a,c,e,f,b,d,f,g),4,4)
  return(Der2) 
  
} 

#### expectation
der2.Score.e1 <- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  s.1 <- matrix(rep(0,16),4,4)
  
  for(l in 1:MCsize){
    V.tilde.nt <- matrix(rnorm(n*T,mean = 0, sd = sqrt((T-1)*sig2/T)),n,T)
    Y.tilde.nt <- Sn.inv%*%(Rn.inv%*%V.tilde.nt+X.tilde.nt*beta)
    Yw <- Rn%*%Wn.a%*%Y.tilde.nt
    
    a <- 2*sum(diag((Mn.a%*%X.tilde.nt)%*%t((Rn%*%X.tilde.nt))))/sig2   ####[1,3]
    b <- sum(diag((Rn%*%X.tilde.nt)%*%t((Rn%*%X.tilde.nt))))/(sig2*sig2)    ####[1,4]
    c <- sum(diag((Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(Rn%*%X.tilde.nt)+(Rn%*%Wn.a%*%Y.tilde.nt)%*%t(Mn.a%*%X.tilde.nt)))/sig2   ####[2,3]
    d <- sum(diag((Rn%*%Wn.a%*%Y.tilde.nt)%*%t(Rn%*%X.tilde.nt)))/(sig2*sig2)   ####[2,4]
    e <- 2*sum(diag((Hn%*%V.tilde.nt)%*%t(Mn.a%*%X.tilde.nt)))/sig2
    f <- sum(diag((Hn%*%V.tilde.nt)%*%t(Rn%*%X.tilde.nt)+V.tilde.nt%*%t(Mn.a%*%X.tilde.nt)))/(sig2*sig2)
    g <- 2*sum(diag(V.tilde.nt%*%t((Rn%*%X.tilde.nt))))/(sig2*sig2*sig2)
    
    s <- matrix(c(0,0,a,b,0,0,c,d,a,c,e,f,b,d,f,g),4,4)
    
    
    
    s.1 <- s+s.1 
  }
  s.1 <- s.1/MCsize
  return(s.1)
}



der2.Score.m1 <- function(i,beta,lambda,rho,sig2){ 
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn 
  
  a <- 2*n*(((Mn.a%*%X.tilde.nt)%*%(t(Rn%*%X.tilde.nt)))[i,i])/sig2   ####[1,3]
  b <- n*(((Rn%*%X.tilde.nt)%*%(t(Rn%*%X.tilde.nt)))[i,i])/(sig2*sig2)    ####[1,4]
  c <- n*beta*((((Mn.a%*%Gn%*%X.tilde.nt)%*%t(Rn%*%X.tilde.nt))[i,i])+(((Rn%*%Gn%*%X.tilde.nt)%*%t(Mn.a%*%X.tilde.nt))[i,i]))/sig2   ####[2,3]
  d <- n*beta*(((Rn%*%Gn%*%X.tilde.nt)%*%t(Rn%*%X.tilde.nt))[i,i])/(sig2*sig2)   ####[2,4]
  e <- 0
  f <- 0
  g <- 0
  
  Der2 <- matrix(c(0,0,a,b,0,0,c,d,a,c,e,f,b,d,f,g),4,4)
  return(Der2) 
  
} 

# der2.Score1 <- der2.Score.e1(beta,lambda0,rho,sig2)

s <- matrix(rep(0,4*4),4,4)

for(i in 1:n){
  s <- s+ der2.Score.m1(i,beta0,lambda,rho,sig2)
}
der2.Score1 <- s/n
#####2nd 
der2.Score.it.2 <- function(i,beta,lambda,rho,sig2){ 
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  
  
  a <- n*(T-1)*((Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t]*(Rn%*%X.tilde.nt)[i,t]+(Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*(Mn.a%*%X.tilde.nt)[i,t])/sig2   ####[1,3]
  b <- n*(T-1)*(Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*(Rn%*%X.tilde.nt)[i,t]/(sig2*sig2)   ####[1,4]
  c <- -2*n*(T-1)*(T-1)*G3n[i,i]/T   ####[2,2]
  d <- 2*n*(T-1)*((Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t]*(Rn%*%Wn.a%*%Y.tilde.nt)[i,t])/sig2    ####[2,3]
  e <- n*(T-1)*((Rn%*%Mn.a%*%Y.tilde.nt)[i,t]*((Rn%*%Wn.a%*%Y.tilde.nt)[i,t]))/(sig2*sig2)          #####[2,4]
  f <- 2*n*(T-1)*((Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t])*((Hn%*%V.tilde.nt)[i,t])/sig2           #####[3,3]
  g <- n*(T-1)*((Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*((Hn%*%V.tilde.nt)[i,t])+(Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t]*V.tilde.nt[i,t])/(sig2*sig2)  ####[3,4]
  h <- n*(T-1)*2*(Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*V.tilde.nt[i,t]/(sig2*sig2*sig2)
  Der2 <- matrix(c(0,0,a,b,0,c,d,e,a,d,f,g,b,e,g,h),4,4)
  return(Der2) 
  
}

#### expectation 

der2.Score.e2 <- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  s.1 <- matrix(rep(0,16),4,4)
  
  for(l in 1:MCsize){
    V.tilde.nt <- matrix(rnorm(n*T,mean = 0, sd = sqrt((T-1)*sig2/T)),n,T)
    Y.tilde.nt <- Sn.inv%*%(Rn.inv%*%V.tilde.nt+X.tilde.nt*beta)
    Yw <- Rn%*%Wn.a%*%Y.tilde.nt
    
    a <- sum(diag((Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(Rn%*%X.tilde.nt)+(Rn%*%Wn.a%*%Y.tilde.nt)%*%t(Mn.a%*%X.tilde.nt)))/sig2   ####[1,3]
    b <- sum(diag(Rn%*%Wn.a%*%Y.tilde.nt%*%t(Rn%*%X.tilde.nt)))/(sig2*sig2)   ####[1,4]
    c <- -2*(T-1)*sum(diag(G3n))   ####[2,2]
    d <- 2*sum(diag((Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(Rn%*%Wn.a%*%Y.tilde.nt)))/sig2    ####[2,3]
    e <- sum(diag((Rn%*%Mn.a%*%Y.tilde.nt)%*%(t(Rn%*%Wn.a%*%Y.tilde.nt))))/(sig2*sig2)          #####[2,4]
    f <- 2*sum(diag((Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(Hn%*%V.tilde.nt)))/sig2           #####[3,3]
    g <- sum(diag((Rn%*%Wn.a%*%Y.tilde.nt)%*%t((Hn%*%V.tilde.nt))+(Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(V.tilde.nt)))/(sig2*sig2)  ####[3,4]
    h <- 2*sum(diag((Rn%*%Wn.a%*%Y.tilde.nt)%*%t(V.tilde.nt)))/(sig2*sig2*sig2)
    s <- matrix(c(0,0,a,b,0,c,d,e,a,d,f,g,b,e,g,h),4,4)
    
    
    
    
    s.1 <- s+s.1 
  }
  s.1 <- s.1/MCsize
  return(s.1)
}
der2.Score.m2 <- function(i,beta,lambda,rho,sig2){ 
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  G.dot <- Rn%*%Gn%*%solve(Rn)
  X.dot <- Rn%*%X.tilde.nt
  
  a <- n*beta*((((Mn.a%*%Gn%*%X.tilde.nt)%*%t(Rn%*%X.tilde.nt))[i,i])+(((Rn%*%Gn%*%X.tilde.nt)%*%t(Mn.a%*%X.tilde.nt))[i,i]))/sig2   ####[1,3]
  b <- n*beta*(((Rn%*%Gn%*%X.tilde.nt)%*%t(Rn%*%X.tilde.nt))[i,i])/(sig2*sig2)   ####[1,4]
  c <- -2*n*(T-1)*G3n[i,i]   ####[2,2]
  d <- 2*n*((T-1)*sig2*((Hn%*%G.dot)%*%t(G.dot))[i,i]+beta*beta*(Hn%*%G.dot%*%X.dot%*%t(G.dot%*%X.dot))[i,i])/sig2    ####[2,3]
  e <- n*((T-1)*sig2*((Rn%*%Mn.a%*%solve(Sn)%*%solve(Rn)%*%t(G.dot)))[i,i]+beta*beta*(((Rn%*%Mn.a%*%solve(Sn)%*%X.tilde.nt)%*%t(G.dot%*%X.dot))[i,i]))/(sig2*sig2)          #####[2,4]
  f <- 2*n*(T-1)*((Hn%*%G.dot%*%t(Hn))[i,i])           #####[3,3]
  g <- n*(T-1)*((G.dot%*%t(Hn))[i,i]+(Hn%*%G.dot)[i,i])/(sig2)  ####[3,4]
  h <- n*(T-1)*2*(G.dot[i,i])/(sig2*sig2)
  Der2 <- matrix(c(0,0,a,b,0,c,d,e,a,d,f,g,b,e,g,h),4,4)
  return(Der2) 
  
} 

# der2.Score2 <- der2.Score.e2(beta,lambda0,rho,sig2)

s <- matrix(rep(0,4*4),4,4)

for(i in 1:n){
  s <- s+ der2.Score.m2(i,beta0,lambda,rho,sig2)
}
der2.Score2 <- s/n
#####3rd 

der2.Score.it.3 <- function(i,t,beta,lambda,rho,sig2){ 
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  H3n <- H2n%*%Hn
  
  a <- 2*n*(T-1)*((Mn.a%*%X.tilde.nt)[i,t]*((Rn%*%X.tilde.nt)[i,t]))/sig2   ####[1,1]
  b <- n*(T-1)*((Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t]*(Rn%*%X.tilde.nt)[i,t]+(Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*(Mn.a%*%X.tilde.nt)[i,t])/sig2   ####[1,2]
  c <- 2*n*(T-1)*(Hn%*%V.tilde.nt)[i,t]*(Mn.a%*%X.tilde.nt)[i,t]/sig2     ####[1,3]
  d <- n*(T-1)*((Hn%*%V.tilde.nt)[i,t]*(Rn%*%X.tilde.nt)[i,t]+V.tilde.nt[i,t]*(Mn.a%*%X.tilde.nt)[i,t])/(sig2*sig2)   ####[1,4]
  e <- 2*n*(T-1)*((Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t]*(Rn%*%Wn.a%*%Y.tilde.nt)[i,t])/sig2    ####[2,2]
  f <- 2*n*(T-1)*((Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t])*((Hn%*%V.tilde.nt)[i,t])/sig2           #####[2,3]
  g <- n*(T-1)*((Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*((Hn%*%V.tilde.nt)[i,t])+(Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t]*V.tilde.nt[i,t])/(sig2*sig2)  ####[2,4]
  h <- -n*(T-1)*2*(T-1)*H3n[i,i]/T     ####[3,3]
  k <- n*(T-1)*((Hn%*%V.tilde.nt)[i,t]*(Hn%*%V.tilde.nt)[i,t])/(sig2*sig2)    #####[3,4]
  l <- n*(T-1)*2*((Hn%*%V.tilde.nt)[i,t]*V.tilde.nt[i,t])/(sig2*sig2*sig2)
  
  Der2 <- matrix(c(a,b,c,d,b,e,f,g,c,f,h,k,d,g,k,l),4,4)
  return(Der2) 
  
} 

#### expectation 
der2.Score.e3 <- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  H3n <- H2n%*%Hn
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  s.1 <- matrix(rep(0,16),4,4)
  
  for(l in 1:MCsize){
    V.tilde.nt <- matrix(rnorm(n*T,mean = 0, sd = sqrt((T-1)*sig2/T)),n,T)
    Y.tilde.nt <- Sn.inv%*%(Rn.inv%*%V.tilde.nt+X.tilde.nt*beta)
    Yw <- Rn%*%Wn.a%*%Y.tilde.nt
    
    a <- 2*sum(diag((Mn.a%*%X.tilde.nt)%*%t(Rn%*%X.tilde.nt)))/sig2   ####[1,1]
    b <- sum(diag((Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(Rn%*%X.tilde.nt)+(Rn%*%Wn.a%*%Y.tilde.nt)%*%t(Mn.a%*%X.tilde.nt)))/sig2   ####[1,2]
    c <- 2*sum(diag(Hn%*%V.tilde.nt%*%t(Mn.a%*%X.tilde.nt)))/sig2     ####[1,3]
    d <- sum(diag((Hn%*%V.tilde.nt)%*%t(Rn%*%X.tilde.nt)+V.tilde.nt%*%t(Mn.a%*%X.tilde.nt)))/(sig2*sig2)   ####[1,4]
    e <- 2*sum(diag((Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(Rn%*%Wn.a%*%Y.tilde.nt)))/sig2    ####[2,2]
    f <- 2*sum(diag((Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(Hn%*%V.tilde.nt)))/sig2           #####[2,3]
    g <- sum(diag((Rn%*%Wn.a%*%Y.tilde.nt%*%t(Hn%*%V.tilde.nt))+(Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(V.tilde.nt)))/(sig2*sig2)  ####[2,4]
    h <- -2*(T-1)*sum(diag(H3n))   ####[3,3]
    k <- sum(diag((Hn%*%V.tilde.nt)%*%t(Hn%*%V.tilde.nt)))/(sig2*sig2)    #####[3,4]
    l <- 2*sum(diag((Hn%*%V.tilde.nt)%*%t(V.tilde.nt)))/(sig2*sig2*sig2)
    
    s <- matrix(c(a,b,c,d,b,e,f,g,c,f,h,k,d,g,k,l),4,4)
    
    
    
    
    s.1 <- s+s.1 
  }
  s.1 <- s.1/MCsize
  return(s.1)
}
der2.Score.m3 <- function(i,beta,lambda,rho,sig2){ 
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  H3n <- H2n%*%Hn
  G.dot <- Rn%*%Gn%*%solve(Rn)
  X.dot <- Rn%*%X.tilde.nt
  
  a <- 2*n*((Mn.a%*%X.tilde.nt)%*%(t(Rn%*%X.tilde.nt)))[i,i]/sig2   ####[1,1]
  b <- n*beta*((((Mn.a%*%Gn%*%X.tilde.nt)%*%t(Rn%*%X.tilde.nt))[i,i])+(((Rn%*%Gn%*%X.tilde.nt)%*%t(Mn.a%*%X.tilde.nt))[i,i]))/sig2   ####[1,2]
  c <- 0    ####[1,3]
  d <- 0  ####[1,4]
  e <- 2*n*((T-1)*sig2*((Hn%*%G.dot)%*%t(G.dot))[i,i]+beta*beta*(Hn%*%G.dot%*%X.dot%*%t(G.dot%*%X.dot))[i,i])/sig2    ####[2,2]
  f <- 2*n*(T-1)*((Hn%*%G.dot%*%t(Hn))[i,i])           #####[2,3]
  g <- n*(T-1)*((G.dot%*%t(Hn))[i,i]+(Hn%*%G.dot)[i,i])/(sig2)  ####[2,4]
  h <- -n*(T-1)*2*H3n[i,i]     ####[3,3]
  k <- n*(T-1)*((Hn%*%t(Hn))[i,i])/(sig2)    #####[3,4]
  l <- n*(T-1)*2*(Hn[i,i])/(sig2*sig2)
  
  Der2 <- matrix(c(a,b,c,d,b,e,f,g,c,f,h,k,d,g,k,l),4,4)
  return(Der2) 
  
} 

# der2.Score3 <- der2.Score.e3(beta,lambda0,rho,sig2)

s <- matrix(rep(0,4*4),4,4)

for(i in 1:n){
  s <- s+ der2.Score.m3(i,beta0,lambda,rho,sig2)
}
der2.Score3 <- s/n
#####4th 

der2.Score.it.4 <- function(i,t,beta,lambda,rho,sig2){ 
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  H3n <- H2n%*%Hn
  
  a <- n*(T-1)*((Rn%*%X.tilde.nt)[i,t])*((Rn%*%X.tilde.nt)[i,t])/(sig2*sig2)   ####[1,1]
  b <- n*(T-1)*(Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*(Rn%*%X.tilde.nt)[i,t]/(sig2*sig2)     ####[1,2]
  c <- 2*n*(T-1)*((Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t])*((Hn%*%V.tilde.nt)[i,t])/sig2     ####[1,3]
  d <- 2*n*(T-1)*V.tilde.nt[i,t]*((Rn%*%X.tilde.nt)[i,t])/(sig2*sig2*sig2)   ####[1,4]
  e <- n*(T-1)*((Rn%*%Mn.a%*%Y.tilde.nt)[i,t]*((Rn%*%Wn.a%*%Y.tilde.nt)[i,t]))/(sig2*sig2)    ####[2,2]
  f <- n*(T-1)*((Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*((Hn%*%V.tilde.nt)[i,t])+(Mn.a%*%Wn.a%*%Y.tilde.nt)[i,t]*V.tilde.nt[i,t])/(sig2*sig2)           #####[2,3]
  g <- n*(T-1)*2*(Rn%*%Wn.a%*%Y.tilde.nt)[i,t]*V.tilde.nt[i,t]/(sig2*sig2*sig2)  ####[2,4]
  h <- n*(T-1)*((Hn%*%V.tilde.nt)[i,t]*(Hn%*%V.tilde.nt)[i,t])/(sig2*sig2)      ####[3,3]
  k <- n*(T-1)*2*((Hn%*%V.tilde.nt)[i,t]*V.tilde.nt[i,t])/(sig2*sig2*sig2)    #####[3,4]
  l <- -n*(T-1)*(T-1)/(T*sig2*sig2*sig2)+3*n*(T-1)*((V.tilde.nt)[i,t]*V.tilde.nt[i,t])/(sig2*sig2*sig2*sig2)   ####[4,4]
  
  Der2 <- matrix(c(a,b,c,d,b,e,f,g,c,f,h,k,d,g,k,l),4,4)
  return(Der2) 
  
} 

#### expectation 
der2.Score.e4 <- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  H3n <- H2n%*%Hn
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  s.1 <- matrix(rep(0,16),4,4)
  
  for(l in 1:MCsize){
    V.tilde.nt <- matrix(rnorm(n*T,mean = 0, sd = sqrt((T-1)*sig2/T)),n,T)
    Y.tilde.nt <- Sn.inv%*%(Rn.inv%*%V.tilde.nt+X.tilde.nt*beta)
    Yw <- Rn%*%Wn.a%*%Y.tilde.nt
    
    a <- sum(diag((Rn%*%X.tilde.nt)%*%t(Rn%*%X.tilde.nt)))/(sig2*sig2)   ####[1,1]
    b <- sum(diag(Rn%*%Wn.a%*%Y.tilde.nt%*%t(Rn%*%X.tilde.nt)))/(sig2*sig2)     ####[1,2]
    c <- 2*sum(diag((Mn.a%*%Wn.a%*%Y.tilde.nt)%*%t(Hn%*%V.tilde.nt)))/sig2     ####[1,3]
    d <- 2*sum(diag(V.tilde.nt%*%t(Rn%*%X.tilde.nt)))/(sig2*sig2*sig2)   ####[1,4]
    e <- sum(diag(Rn%*%Mn.a%*%Y.tilde.nt%*%t(Rn%*%Wn.a%*%Y.tilde.nt)))/(sig2*sig2)    ####[2,2]
    f <- sum(diag((Rn%*%Wn.a%*%Y.tilde.nt%*%t(Hn%*%V.tilde.nt))+(Mn.a%*%Wn.a%*%Y.tilde.nt%*%t(V.tilde.nt))))/(sig2*sig2)           #####[2,3]
    g <- 2*sum(diag(Rn%*%Wn.a%*%Y.tilde.nt%*%t(V.tilde.nt)))/(sig2*sig2*sig2)  ####[2,4]
    h <- sum(diag((Hn%*%V.tilde.nt)%*%t(Hn%*%V.tilde.nt)))/(sig2*sig2)      ####[3,3]
    k <- 2*sum(diag((Hn%*%V.tilde.nt)%*%t(V.tilde.nt)))/(sig2*sig2*sig2)    #####[3,4]
    l <- -n*(T-1)/(sig2*sig2*sig2)+3*sum(diag((V.tilde.nt)%*%t(V.tilde.nt)))/(sig2*sig2*sig2*sig2)   ####[4,4]
    
    s <- matrix(c(a,b,c,d,b,e,f,g,c,f,h,k,d,g,k,l),4,4)
    
    
    
    s.1 <- s+s.1 
  }
  s.1 <- s.1/MCsize
  return(s.1)
}
der2.Score.m4 <- function(i,beta,lambda,rho,sig2){ 
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt-X.tilde.nt*beta)
  Yw <- Rn%*%Wn.a%*%Y.tilde.nt 
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  H3n <- H2n%*%Hn
  G.dot <- Rn%*%Gn%*%solve(Rn)
  X.dot <- Rn%*%X.tilde.nt
  
  a <- n*(((Rn%*%X.tilde.nt)%*%t(Rn%*%X.tilde.nt))[i,i])/(sig2*sig2)   ####[1,1]
  b <- n*beta*((G.dot%*%X.dot)%*%(t(Rn%*%X.tilde.nt)))[i,i]/(sig2*sig2)     ####[1,2]
  c <- 2*n*(T-1)*((Hn%*%G.dot%*%t(Hn))[i,i])    ####[1,3]
  d <- 0  ####[1,4]
  e <- n*((T-1)*sig2*((Rn%*%Mn.a%*%solve(Sn)%*%solve(Rn)%*%t(G.dot)))[i,i]+beta*beta*(((Rn%*%Mn.a%*%solve(Sn)%*%X.tilde.nt)%*%t(G.dot%*%X.dot))[i,i]))/(sig2*sig2)   ####[2,2]
  f <- n*(T-1)*((G.dot%*%t(Hn))[i,i]+(Hn%*%G.dot)[i,i])/(sig2)         #####[2,3]
  g <- n*(T-1)*2*(G.dot[i,i])/(sig2*sig2)  ####[2,4]
  h <-  n*(T-1)*((Hn%*%t(Hn))[i,i])/(sig2)       ####[3,3]
  k <- n*(T-1)*2*(Hn[i,i])/(sig2*sig2)    #####[3,4]
  l <- 2*n*(T-1)/(sig2*sig2*sig2)   ####[4,4]
  
  Der2 <- matrix(c(a,b,c,d,b,e,f,g,c,f,h,k,d,g,k,l),4,4)
  return(Der2) 
  
} 

# der2.Score4 <- der2.Score.e4(beta,lambda0,rho,sig2)

s <- matrix(rep(0,4*4),4,4)

for(i in 1:n){
  s <- s+ der2.Score.m4(i,beta0,lambda,rho,sig2)
}
der2.Score4 <- s/n
###### Gamma function 

Gamma.ij <- function(i,j,beta,lambda,rho,sig2){ 
  
  gamma.1 <- t(IF.iT(j,beta,lambda,rho,sig2))%*%der2.Score1%*%IF.iT(i,beta,lambda,rho,sig2) 
  gamma.2 <- t(IF.iT(j,beta,lambda,rho,sig2))%*%der2.Score2%*%IF.iT(i,beta,lambda,rho,sig2)
  gamma.3 <- t(IF.iT(j,beta,lambda,rho,sig2))%*%der2.Score3%*%IF.iT(i,beta,lambda,rho,sig2)
  gamma.4 <- t(IF.iT(j,beta,lambda,rho,sig2))%*%der2.Score4%*%IF.iT(i,beta,lambda,rho,sig2)
  gamma <- c(gamma.1,gamma.2,gamma.3,gamma.4)
  return(gamma) 
  
}  

###### The second term of Von Mises expansion 

Phi.ij <- function(i,j,beta,lambda,rho,sig2){ 
  phi <- (IF.iT(i,beta,lambda,rho,sig2)+IF.iT(j,beta,lambda,rho,sig2)+solve(M)%*%Gamma.ij(i,j,beta,lambda,rho,sig2)
          +(1/(T-1))*solve(M)%*%(Sum.T.Der1(j,beta,lambda,rho,sig2)%*%IF.iT(i,beta,lambda,rho,sig2)+Sum.T.Der1(i,beta,lambda,rho,sig2)%*%IF.iT(j,beta,lambda,rho,sig2)) )
  return(phi) 
  
} 

###### g g^2 g^3 g^4


G.powers<- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  H3n <- H2n%*%Hn
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  G.dot <- Rn%*%Gn%*%Rn.inv
  X.dot <- Rn%*%X.tilde.nt
  
  Psi.mat <- matrix(rep(0,(4*MCsize)),MCsize,4)
  for(j in 1:MCsize){
    V.tilde <- matrix(rnorm(n*T,mean = 0, sd = sqrt((T-1)*sig2/T)),n,T)
    Y.tilde <- Sn.inv%*%(Rn.inv%*%V.tilde+X.tilde.nt*beta)
    Yw.1 <- Rn%*%Wn.a%*%Y.tilde
    # a <- c(sum(diag((Rn%*%X.tilde.nt)%*%t(V.tilde)))/sig2, sum(diag(Yw.1%*%t(V.tilde)))/sig2 - (T-1)*sum(diag(Gn)),
    #        sum(diag((Hn%*%V.tilde)%*%t(V.tilde)))/sig2 - (T-1)*sum(diag(Hn)), sum(diag((V.tilde%*%t(V.tilde))))/(2*sig2*sig2)-n*(T-1)/(2*sig2))
    # Psi.mat[j,1] <- t(solve(M)%*%a)%*%c(0,1,0,0)
    g1 <- vector(length = length(n))
    g2 <- vector(length = length(n))
    g3 <- vector(length = length(n))
    g4 <- vector(length = length(n))
    for(i in 1:n){
      
      b <- c(n*((Rn%*%X.tilde.nt)%*%t(V.tilde))[i,i]/sig2, n*(Yw.1%*%t(V.tilde))[i,i]/sig2 - n*((T-1))*Gn[i,i],
             n*((Hn%*%V.tilde)%*%t(V.tilde))[i,i]/sig2 - n*((T-1))*Hn[i,i], n/(2*sig2*sig2)*((V.tilde%*%t(V.tilde))[i,i]-(T-1)*sig2))
      g1[i] <- t(solve(M)%*%b)%*%c(1,0,0,0)
      g2[i] <- g1[i]*g1[i]
      g3[i] <- g2[i]*(t(solve(M)%*%b)%*%c(1,0,0,0))
      g4[i] <- g3[i]*(t(solve(M)%*%b)%*%c(1,0,0,0))
      
    }
    Psi.mat[j,1] <- sum(g1)
    Psi.mat[j,2] <- sum(g2)
    Psi.mat[j,3] <- sum(g3)
    Psi.mat[j,4] <- sum(g4)
  }
  Psi.pwr<- colMeans(Psi.mat)
  Psi.pwr[1] <- Psi.pwr[1]/(2*n)
  Psi.pwr[2] <- Psi.pwr[2]/(4*n)
  Psi.pwr[3] <- Psi.pwr[3]/(8*n)
  Psi.pwr[4] <- Psi.pwr[4]/(16*n)
  return(Psi.pwr)
}

g.powers <- G.powers(beta0,lambda,rho,sig2)

###### gamma
MCsize1 <- 200
Gam.Con1<- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  H3n <- H2n%*%Hn
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  G.dot <- Rn%*%Gn%*%Rn.inv
  X.dot <- Rn%*%X.tilde.nt
  M.inv <- solve(M)
  
  gam.mat <- matrix(rep(0,(4*MCsize1)),MCsize1,4) 
  
  for(k in 1:MCsize1){
    gamm1 <- matrix(rep(0,n*n),n,n)
    gamm2 <- matrix(rep(0,n*n),n,n)
    gamm3 <- matrix(rep(0,n*n),n,n)
    gamm4 <- matrix(rep(0,n*n),n,n)
    
    V.tilde <- matrix(rnorm(n*T,mean = 0, sd = sqrt((T-1)*sig2/T)),n,T)
    Y.tilde <- Sn.inv%*%(Rn.inv%*%V.tilde+X.tilde.nt*beta)
    Yw.1 <- Rn%*%Wn.a%*%Y.tilde
    
    for(i in 1:n){
      for(j in 1:n){
        
        a <- c(n*(((Rn%*%X.tilde.nt)%*%t(V.tilde))[i,i])/sig2, n*(Yw.1%*%t(V.tilde))[i,i]/sig2 - n*(T-1)*Gn[i,i],
               n*((Hn%*%V.tilde)%*%t(V.tilde))[i,i]/sig2 - n*(T-1)*Hn[i,i], n*((V.tilde%*%t(V.tilde))[i,i]-(T-1)*sig2)/(2*sig2*sig2))
        IF.i <- M.inv%*%a
        b <- c(n*(((Rn%*%X.tilde.nt)%*%t(V.tilde))[j,j])/sig2, n*(Yw.1%*%t(V.tilde))[j,j]/sig2 - n*(T-1)*Gn[j,j],
               n*((Hn%*%V.tilde)%*%t(V.tilde))[j,j]/sig2 - n*(T-1)*Hn[j,j], n*((V.tilde%*%t(V.tilde))[j,j]-(T-1)*sig2)/(2*sig2*sig2))
        IF.j <- M.inv%*%b
        # c <- c(n*(((Rn%*%X.tilde.nt)%*%t(V.tilde))[l,l])/sig2, n*(Yw.1%*%t(V.tilde))[l,l]/sig2 - n*(T-1)*Gn[l,l],
        #          n*((Hn%*%V.tilde)%*%t(V.tilde))[l,l]/sig2 - n*(T-1)*Hn[l,l], n*((V.tilde%*%t(V.tilde))[l,l]-(T-1)*sig2)/(2*sig2*sig2))
        # IF.l <- M.inv%*%c
        
        gamma.1 <- t(IF.j)%*%der2.Score1%*%IF.i
        gamma.2 <- t(IF.j)%*%der2.Score2%*%IF.i
        gamma.3 <- t(IF.j)%*%der2.Score3%*%IF.i
        gamma.4 <- t(IF.j)%*%der2.Score4%*%IF.i
        Gamma.ij <- c(gamma.1,gamma.2,gamma.3,gamma.4)
        
        a.1 <- -n*((X.dot)%*%t(X.dot))[i,i]/sig2  ###[1,1]
        b.1 <- -n*(Yw.1%*%t(X.dot))[i,i]/sig2                 ###[2,1]
        c.1 <- -n*(((Hn%*%V.tilde)%*%t(X.dot))[i,i]+(V.tilde%*%t(Mn.a%*%X.tilde.nt))[i,i])/sig2  ###[3,1]
        d.1 <- -n*((X.dot)%*%t(V.tilde))[i,i]/(sig2*sig2)   ###[4,1]
        e.1 <- -n*(Yw.1%*%t(Yw.1))[i,i]/sig2 - n*((T-1))*G2n[i,i] ###[2,2]
        f.1 <- -n*((Hn%*%V.tilde)%*%t(Yw.1))[i,i]/sig2-n*((Mn.a%*%Wn.a%*%Y.tilde)%*%t(V.tilde))[i,i]/sig2 ###[2,3]
        g.1 <- -n*(Yw.1%*%t(V.tilde))[i,i]/(sig2*sig2)                 ###[2,4]
        h.1 <- -n*((Hn%*%V.tilde)%*%t(Hn%*%V.tilde))[i,i]/sig2 - n*((T-1))*H2n[i,i]   ###[3,3]
        k.1 <- -n*((Hn%*%V.tilde)%*%t(V.tilde))[i,i]/(sig2*sig2)  ###[3,4]
        j.1 <- -n*((V.tilde%*%t(V.tilde))[i,i])/(sig2*sig2*sig2)+n*(T-1)/(2*sig2*sig2)  ###[4,4]
        
        ### already divided by (T-1)
        Sum.Der1.i <- matrix(c(a.1,b.1,c.1,d.1,b.1,e.1,f.1,g.1,c.1,f.1,h.1,k.1,d.1,g.1,k.1,j.1),4,4)  
        
        a.2 <- -n*((X.dot)%*%t(X.dot))[j,j]/sig2  ###[1,1]
        b.2 <- -n*(Yw.1%*%t(X.dot))[j,j]/sig2                 ###[2,1]
        c.2 <- -n*(((Hn%*%V.tilde)%*%t(X.dot))[j,j]+(V.tilde%*%t(Mn.a%*%X.tilde.nt))[j,j])/sig2 ###[3,1]
        d.2 <- -n*((X.dot)%*%t(V.tilde))[j,j]/(sig2*sig2)   ###[4,1]
        e.2 <- -n*(Yw.1%*%t(Yw.1))[j,j]/sig2 - n*((T-1))*G2n[j,j] ###[2,2]
        f.2 <- -n*((Hn%*%V.tilde)%*%t(Yw.1))[j,j]/sig2-n*((Mn.a%*%Wn.a%*%Y.tilde)%*%t(V.tilde))[j,j]/sig2 ###[2,3]
        g.2 <- -n*(Yw.1%*%t(V.tilde))[j,j]/(sig2*sig2)                 ###[2,4]
        h.2 <- -n*((Hn%*%V.tilde)%*%t(Hn%*%V.tilde))[j,j]/sig2 - n*((T-1))*H2n[j,j]   ###[3,3]
        k.2 <- -n*((Hn%*%V.tilde)%*%t(V.tilde))[j,j]/(sig2*sig2)  ###[3,4]
        j.2 <- -n*((V.tilde%*%t(V.tilde))[j,j])/(sig2*sig2*sig2)+n*(T-1)/(2*sig2*sig2)  ###[4,4]
        
        Sum.Der1.j <- matrix(c(a.2,b.2,c.2,d.2,b.2,e.2,f.2,g.2,c.2,f.2,h.2,k.2,d.2,g.2,k.2,j.2),4,4)  
        
        # ###### i,l
        # gamma.1i <- t(IF.l)%*%der2.Score1%*%IF.i
        # gamma.2i <- t(IF.l)%*%der2.Score2%*%IF.i
        # gamma.3i <- t(IF.l)%*%der2.Score3%*%IF.i
        # gamma.4i <- t(IF.l)%*%der2.Score4%*%IF.i
        # Gamma.il <- c(gamma.1i,gamma.2i,gamma.3i,gamma.4i)
        #  
        # 
        # a.3 <- -n*((X.dot)%*%t(X.dot))[l,l]/sig2  ###[1,1]
        # b.3 <- -n*(Yw.1%*%t(X.dot))[l,l]/sig2                 ###[2,1]
        # c.3 <- -n*((Hn%*%V.tilde)%*%t(X.dot))[l,l]/sig2  ###[3,1]
        # d.3 <- -n*((X.dot)%*%t(V.tilde))[l,l]/(sig2*sig2)   ###[4,1]
        # e.3 <- -n*(Yw.1%*%t(Yw.1))[l,l]/sig2 - n*((T-1))*G2n[l,l] ###[2,2]
        # f.3 <- -n*((Hn%*%V.tilde)%*%t(Yw.1))[l,l]/sig2-n*((Mn.a%*%Wn.a%*%Y.tilde)%*%t(V.tilde))[l,l]/sig2 ###[2,3]
        # g.3 <- -n*(Yw.1%*%t(V.tilde))[l,l]/(sig2*sig2)                 ###[2,4]
        # h.3 <- -n*((Hn%*%V.tilde)%*%t(Hn%*%V.tilde))[l,l]/sig2 - n*((T-1))*H2n[l,l]   ###[3,3]
        # k.3 <- -n*((Hn%*%V.tilde)%*%t(V.tilde))[l,l]/(sig2*sig2)  ###[3,4]
        # j.3 <- -n*((V.tilde%*%t(V.tilde))[l,l])/(sig2*sig2*sig2)+n*(T-1)/(2*sig2*sig2)  ###[4,4]
        # 
        # Sum.Der1.l <- matrix(c(a.3,b.3,c.3,d.3,b.3,e.3,f.3,g.3,c.3,f.3,h.3,k.3,d.3,g.3,k.3,j.3),4,4)
        # 
        # gam.il <- (IF.i+IF.l+M.inv.i%*%Gamma.il +(1/(T-1))*M.inv.i%*%(Sum.Der1.l%*%IF.i+Sum.Der1.i%*%IF.l) )[2]
        # 
        # ###### j,l
        # gamma.1j <- t(IF.l)%*%der2.Score.m1(i,beta,lambda,rho,sig2)%*%IF.j
        # gamma.2j <- t(IF.l)%*%der2.Score.m2(i,beta,lambda,rho,sig2)%*%IF.j
        # gamma.3j <- t(IF.l)%*%der2.Score.m3(i,beta,lambda,rho,sig2)%*%IF.j
        # gamma.4j <- t(IF.l)%*%der2.Score.m4(i,beta,lambda,rho,sig2)%*%IF.j
        # Gamma.jl <- c(gamma.1j,gamma.2j,gamma.3j,gamma.4j)
        # gam.jl <- (IF.j+IF.l+M.inv.j%*%Gamma.jl +(1/(T-1))*M.inv.j%*%(Sum.Der1.l%*%IF.j+Sum.Der1.j%*%IF.l) )[2]
        gamm1[i,j] <- (IF.i+IF.j+M.inv%*%Gamma.ij+M.inv%*%(Sum.Der1.j%*%IF.i+Sum.Der1.i%*%IF.j) )[1]
        gamm2[i,j] <- gamm1[i,j]*gamm1[i,j]
        gamm3[i,j] <- gamm1[i,j]*IF.i[1]*IF.j[1]
        gamm4[i,j] <- gamm3[i,j]*IF.i[1]
      }
    }
    gam.mat[k,1] <- sum(gamm1)
    gam.mat[k,2] <- sum(gamm2)
    gam.mat[k,3] <- sum(gamm3)
    gam.mat[k,4] <- sum(gamm4)
    # gam.mat[k,5] <- IF.i[2]*IF.j[2]*gam.il*gam.jl
    
    
  }
  
  
  
  gam.pwr<- colMeans(gam.mat)
  gam.pwr[1] <- gam.pwr[1]/(2*n*n)
  gam.pwr[2] <- gam.pwr[2]/(4*n*n)
  gam.pwr[3] <- gam.pwr[3]/(8*n*n)
  gam.pwr[4] <- gam.pwr[4]/(16*n*n)
  # gam.pwr[5] <- gam.pwr[5]/(16)
  
  return(gam.pwr)
  
}

gam.con1 <- Gam.Con1(beta0,lambda,rho,sig2)

Gam.Con2<- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  G2n <- Gn%*%Gn 
  H2n <- Hn%*%Hn
  G3n <- G2n%*%Gn
  H3n <- H2n%*%Hn
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  G.dot <- Rn%*%Gn%*%Rn.inv
  X.dot <- Rn%*%X.tilde.nt
  M.inv <- solve(M)
  
  gam.vec <- vector(length = 5)
  
  for(k in 1:5){
    gamm1 <- 0
    V.tilde <- matrix(rnorm(n*T,mean = 0, sd = sqrt((T-1)*sig2/T)),n,T)
    Y.tilde <- Sn.inv%*%(Rn.inv%*%V.tilde+X.tilde.nt*beta)
    Yw.1 <- Rn%*%Wn.a%*%Y.tilde
    
    for(i in 1:n){
      for(j in 1:n){
        for(l in 1:n){
          
          
          a <- c(n*(((Rn%*%X.tilde.nt)%*%t(V.tilde))[i,i])/sig2, n*(Yw.1%*%t(V.tilde))[i,i]/sig2 - n*(T-1)*Gn[i,i],
                 n*((Hn%*%V.tilde)%*%t(V.tilde))[i,i]/sig2 - n*(T-1)*Hn[i,i], n*((V.tilde%*%t(V.tilde))[i,i]-(T-1)*sig2)/(2*sig2*sig2))
          IF.i <- M.inv%*%a
          b <- c(n*(((Rn%*%X.tilde.nt)%*%t(V.tilde))[j,j])/sig2, n*(Yw.1%*%t(V.tilde))[j,j]/sig2 - n*(T-1)*Gn[j,j],
                 n*((Hn%*%V.tilde)%*%t(V.tilde))[j,j]/sig2 - n*(T-1)*Hn[j,j], n*((V.tilde%*%t(V.tilde))[j,j]-(T-1)*sig2)/(2*sig2*sig2))
          IF.j <- M.inv%*%b
          c <- c(n*(((Rn%*%X.tilde.nt)%*%t(V.tilde))[l,l])/sig2, n*(Yw.1%*%t(V.tilde))[l,l]/sig2 - n*(T-1)*Gn[l,l],
                 n*((Hn%*%V.tilde)%*%t(V.tilde))[l,l]/sig2 - n*(T-1)*Hn[l,l], n*((V.tilde%*%t(V.tilde))[l,l]-(T-1)*sig2)/(2*sig2*sig2))
          IF.l <- M.inv%*%c
          
          # gamma.1 <- t(IF.j)%*%der2.Score1%*%IF.i
          # gamma.2 <- t(IF.j)%*%der2.Score2%*%IF.i
          # gamma.3 <- t(IF.j)%*%der2.Score3%*%IF.i
          # gamma.4 <- t(IF.j)%*%der2.Score4%*%IF.i
          # Gamma.ij <- c(gamma.1,gamma.2,gamma.3,gamma.4)
          
          a.1 <- -n*((X.dot)%*%t(X.dot))[i,i]/sig2  ###[1,1]
          b.1 <- -n*(Yw.1%*%t(X.dot))[i,i]/sig2                 ###[2,1]
          c.1 <- -n*(((Hn%*%V.tilde)%*%t(X.dot))[i,i]+(V.tilde%*%t(Mn.a%*%X.tilde.nt))[i,i])/sig2  ###[3,1]
          d.1 <- -n*((X.dot)%*%t(V.tilde))[i,i]/(sig2*sig2)   ###[4,1]
          e.1 <- -n*(Yw.1%*%t(Yw.1))[i,i]/sig2 - n*((T-1))*G2n[i,i] ###[2,2]
          f.1 <- -n*((Hn%*%V.tilde)%*%t(Yw.1))[i,i]/sig2-n*((Mn.a%*%Wn.a%*%Y.tilde)%*%t(V.tilde))[i,i]/sig2 ###[2,3]
          g.1 <- -n*(Yw.1%*%t(V.tilde))[i,i]/(sig2*sig2)                 ###[2,4]
          h.1 <- -n*((Hn%*%V.tilde)%*%t(Hn%*%V.tilde))[i,i]/sig2 - n*((T-1))*H2n[i,i]   ###[3,3]
          k.1 <- -n*((Hn%*%V.tilde)%*%t(V.tilde))[i,i]/(sig2*sig2)  ###[3,4]
          j.1 <- -n*((V.tilde%*%t(V.tilde))[i,i])/(sig2*sig2*sig2)+n*(T-1)/(2*sig2*sig2)  ###[4,4]
          
          Sum.Der1.i <- matrix(c(a.1,b.1,c.1,d.1,b.1,e.1,f.1,g.1,c.1,f.1,h.1,k.1,d.1,g.1,k.1,j.1),4,4)  
          
          a.2 <- -n*((X.dot)%*%t(X.dot))[j,j]/sig2  ###[1,1]
          b.2 <- -n*(Yw.1%*%t(X.dot))[j,j]/sig2                 ###[2,1]
          c.2 <- -n*(((Hn%*%V.tilde)%*%t(X.dot))[j,j]+(V.tilde%*%t(Mn.a%*%X.tilde.nt))[j,j])/sig2  ###[3,1]
          d.2 <- -n*((X.dot)%*%t(V.tilde))[j,j]/(sig2*sig2)   ###[4,1]
          e.2 <- -n*(Yw.1%*%t(Yw.1))[j,j]/sig2 - n*((T-1))*G2n[j,j] ###[2,2]
          f.2 <- -n*((Hn%*%V.tilde)%*%t(Yw.1))[j,j]/sig2-n*((Mn.a%*%Wn.a%*%Y.tilde)%*%t(V.tilde))[j,j]/sig2 ###[2,3]
          g.2 <- -n*(Yw.1%*%t(V.tilde))[j,j]/(sig2*sig2)                 ###[2,4]
          h.2 <- -n*((Hn%*%V.tilde)%*%t(Hn%*%V.tilde))[j,j]/sig2 - n*((T-1))*H2n[j,j]   ###[3,3]
          k.2 <- -n*((Hn%*%V.tilde)%*%t(V.tilde))[j,j]/(sig2*sig2)  ###[3,4]
          j.2 <- -n*((V.tilde%*%t(V.tilde))[j,j])/(sig2*sig2*sig2)+n*(T-1)/(2*sig2*sig2)  ###[4,4]
          
          Sum.Der1.j <- matrix(c(a.2,b.2,c.2,d.2,b.2,e.2,f.2,g.2,c.2,f.2,h.2,k.2,d.2,g.2,k.2,j.2),4,4)  
          
          
          gamma.1i <- t(IF.l)%*%der2.Score1%*%IF.i
          gamma.2i <- t(IF.l)%*%der2.Score2%*%IF.i
          gamma.3i <- t(IF.l)%*%der2.Score3%*%IF.i
          gamma.4i <- t(IF.l)%*%der2.Score4%*%IF.i
          Gamma.il <- c(gamma.1i,gamma.2i,gamma.3i,gamma.4i)
          
          
          a.3 <- -n*((X.dot)%*%t(X.dot))[l,l]/sig2  ###[1,1]
          b.3 <- -n*(Yw.1%*%t(X.dot))[l,l]/sig2                 ###[2,1]
          c.3 <- -n*(((Hn%*%V.tilde)%*%t(X.dot))[l,l]+(V.tilde%*%t(Mn.a%*%X.tilde.nt))[l,l])/sig2  ###[3,1]
          d.3 <- -n*((X.dot)%*%t(V.tilde))[l,l]/(sig2*sig2)   ###[4,1]
          e.3 <- -n*(Yw.1%*%t(Yw.1))[l,l]/sig2 - n*((T-1))*G2n[l,l] ###[2,2]
          f.3 <- -n*((Hn%*%V.tilde)%*%t(Yw.1))[l,l]/sig2-n*((Mn.a%*%Wn.a%*%Y.tilde)%*%t(V.tilde))[l,l]/sig2 ###[2,3]
          g.3 <- -n*(Yw.1%*%t(V.tilde))[l,l]/(sig2*sig2)                 ###[2,4]
          h.3 <- -n*((Hn%*%V.tilde)%*%t(Hn%*%V.tilde))[l,l]/sig2 - n*((T-1))*H2n[l,l]   ###[3,3]
          k.3 <- -n*((Hn%*%V.tilde)%*%t(V.tilde))[l,l]/(sig2*sig2)  ###[3,4]
          j.3 <- -n*((V.tilde%*%t(V.tilde))[l,l])/(sig2*sig2*sig2)+n*(T-1)/(2*sig2*sig2)  ###[4,4]
          
          Sum.Der1.l <- matrix(c(a.3,b.3,c.3,d.3,b.3,e.3,f.3,g.3,c.3,f.3,h.3,k.3,d.3,g.3,k.3,j.3),4,4)
          
          gam.il <- (IF.i+IF.l+M.inv%*%Gamma.il +M.inv%*%(Sum.Der1.l%*%IF.i+Sum.Der1.i%*%IF.l) )[1]
          
          ###### j,l
          gamma.1j <- t(IF.l)%*%der2.Score.m1(i,beta,lambda,rho,sig2)%*%IF.j
          gamma.2j <- t(IF.l)%*%der2.Score.m2(i,beta,lambda,rho,sig2)%*%IF.j
          gamma.3j <- t(IF.l)%*%der2.Score.m3(i,beta,lambda,rho,sig2)%*%IF.j
          gamma.4j <- t(IF.l)%*%der2.Score.m4(i,beta,lambda,rho,sig2)%*%IF.j
          Gamma.jl <- c(gamma.1j,gamma.2j,gamma.3j,gamma.4j)
          gam.jl <- (IF.j+IF.l+M.inv%*%Gamma.jl +M.inv%*%(Sum.Der1.l%*%IF.j+Sum.Der1.j%*%IF.l) )[1]
          
          gamm1 <- gamm1+IF.i[1]*IF.j[1]*gam.il*gam.jl
          
        }
      }
    }
    gam.vec[k] <- gamm1
    
    
  }
  
  
  
  gam.pwr<- sum(gam.vec)/5
  gam.pwr <- gam.pwr/(16*n*n*n)
  
  
  return(gam.pwr)
  
}


Sigma2.nT <- (4/n)*g.powers[2]+(2/(n*(n-1)))*gam.con1[2] 

Sigma.nT <- sqrt(Sigma2.nT)
Sigma3.nT <- Sigma.nT*Sigma2.nT
Sigma4.nT <- Sigma2.nT*Sigma2.nT

###### The third order cumulant 

cum3.nt <-  (sqrt(g.powers[2]))^(-3)*(g.powers[3]+3*gam.con1[3]) 


###### The fourth order cumulant 

cum4.nt <- (sqrt(g.powers[2]))^(-4)*(g.powers[4] +12*gam.con1[4]+12*Gam.Con2(beta0,lambda,rho,sig2))-3 


###### cumulant generating function 

cgf <- function(u){
  c <- gam.con1[1]/n*u + 0.5* n*Sigma2.nT*u*u +(1/6)*(n^(1.5))*cum3.nt*(Sigma3.nT)*(u*u*u)+(1/24)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u*u*u)
  return(c)
}
seq.b <- seq(-99,99,1)
Cgf <- NULL
for(i in 1:length(seq.b)){
  Cgf[i] <- cgf(seq.b[i])
} 
plot(seq.b,Cgf)


###### The first derivative of cgf

Der1.cgf <- function(u){
  der1 <- gam.con1[1]/n + n*Sigma2.nT*u+(1/2)*(n^(1.5))*cum3.nt*(Sigma3.nT)*(u*u)+(1/6)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u*u)
  return(der1)
}
seq.b <- seq(-7,7,0.1)
der1.cgf <- NULL
for(i in 1:length(seq.b)){
  der1.cgf[i] <- Der1.cgf(seq.b[i])
}

plot(seq.b,der1.cgf)

###### The second derivative of cgf

Der2.cgf <- function(u){
  der2 <-  n*Sigma2.nT+(n^(1.5))*cum3.nt*(Sigma3.nT)*(u)+(1/2)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u)
  return(der2)
}
seq.b <- seq(-0.99,0.99,0.01)
der2.cgf <- NULL
for(i in 1:length(seq.b)){
  der2.cgf[i] <- Der2.cgf(seq.b[i])
}

plot(seq.b,der2.cgf)

Der2.cgf.Wang <- function(u,a){
  wn <- exp(-n*Sigma2.nT*a*a*u*u/2) 
  der2 <-  (n*Sigma2.nT+((n^(1.5))*cum3.nt*(Sigma3.nT)*(u)+(1/2)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u))*wn+
              2*((1/2)*(n^(1.5))*cum3.nt*(Sigma3.nT)*(u*u)+(1/6)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u*u))*(-n*Sigma2.nT*a*a*u)*wn+
              ((1/6)*(n^(1.5))*cum3.nt*(Sigma3.nT)*(u*u*u)+(1/24)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u*u*u))*wn*(n*Sigma2.nT*a*a-1)*n*Sigma2.nT*a*a)
  return(der2)
}
seq.b <- seq(-5,5,0.1)
der2.cgf.w <- NULL
for(i in 1:length(seq.b)){
  der2.cgf.w[i] <- Der2.cgf.Wang(seq.b[i],1.575)
}

plot(seq.b,der2.cgf.w)
###### find saddlepoints
Sad <- function(a){
  
  
  sad <- uniroot(function(u) Der1.cgf(u)-a, lower = -10000,upper = 10000)$root
  
  return(sad)
}



sad.grid<-seq(-0.99,0.99,0.01)
sad.vals<-NULL
for (i in 1:length(sad.grid)) {
  sad.vals[i]<-Sad(sad.grid[i])
}

plot(sad.grid,sad.vals)

###### pdf


p.nT <- function(a){
  p <- sqrt(n/(2*pi*Der2.cgf(Sad(a))))*exp(n*(cgf(Sad(a))-a*Sad(a)))
  return(p)
}

theta.grid<-seq(-0.9999,0.9999,by=0.0001)
p<-NULL
for (i in 1:length(theta.grid)) {
  p[i]<-p.nT(theta.grid[i])  
} 

c.int<-sum(p)*diff(theta.grid)[1]
p.std<-p/c.int

plot(theta.grid,p.std)

CDF.SAD1 <- function(b){
  theta.grid<-seq(-0.9999,0.9999,by=0.0001)
  p<-NULL
  for (i in 1:(10000*b+10000)) {
    p[i]<-p.nT(theta.grid[i])  
  } 
  
  c<-sum(p)*diff(theta.grid)[1]/c.int
  return(c)
}

1-CDF.SAD1(round(sarar$coefficients[3],4))   ####  sparse weights matrix:0.0455737, Wdia1: multiple solutions of sad


#### L-R 
CDF.SAD <- function(b){
  v <- Sad(b)
  c <- v*sqrt(Der2.cgf(v))
  r <- sign(v)*sqrt(2*n*(v*b-cgf(v)))
  p <- 1-pnorm(r)+dnorm(r)*(1/c-1/r)
  return(p)
}

CDF.SAD(sarar$coefficients[3])


#### Asymptotic varience

SIGMA <- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn.a
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn.a%*%solve(Rn)
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  G.dot <- Rn%*%Gn%*%Rn.inv
  X.dot <- Rn%*%X.tilde.nt
  
  a <- sum(diag(t(X.dot)%*%X.dot))/(m*sig2)
  b <- beta*(sum(diag(t(X.dot)%*%G.dot%*%X.dot)))/(m*sig2)
  c<- 1/n*(sum(diag(t(G.dot%*%X.dot)%*%(G.dot%*%X.dot)))*beta*beta/((T-1)*sig2)+sum(diag((G.dot+t(G.dot))%*%G.dot)))
  d <- sum(diag((t(Hn)+Hn)%*%G.dot))/n
  e <- sum(diag(G.dot))/(n*sig2)
  f <- sum(diag((t(Hn)+Hn)%*%Hn))/n
  g <- sum(diag(Hn))/(n*sig2)
  h <- 1/(2*sig2*sig2)
  Sig <-  matrix(c(a,b,0,0,b,c,d,e,0,d,f,g,0,e,g,h),4,4)
  return(Sig)
} 

Asy.sig2 <- solve(SIGMA(beta,lambda,rho,sig2))[1,1]/m
1-pnorm(sarar$coefficients[3],mean=0,sd=sqrt(Asy.sig2))   ### 0.04464198 

