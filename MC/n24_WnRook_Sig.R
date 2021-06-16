library("spdep") 
library("splm") 
library("graphics")


######### LOAD the MC results for MLE, z1 data from SAR_mle.rdata
#load("SAR_mle.rdata")

r1 <- 6 
c1 <-4 
n <- r1*c1   
T <- 2 
N <- n*T 

beta <- 0.0 
lambda1 <- 0.0
pho1 <- 0.0 
sig2 <- 1 
lambda0 <- 0.0
set.seed(1234) 

W <- cell2nb(r1,c1)    #rook type
Wl <- nb2listw(W)   #listw object 
Wn <- listw2dgCMatrix(Wl)  #sparse matrix
Wn.a <- as.matrix(Wn)
Mn.a <- Wn.a 
Cn0 <- runif(n,-0.1,0.1)    #fixed effects 
X <- rnorm(N)     #non stochastic time varying regressors 
V<- rnorm(N)       
Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 
Y.tilde.nt<-matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
X.tilde.nt <- matrix(X,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(X,nrow=n,ncol=T)),T),n,T)

r2 <- 50
n1 <- r2*r2   #n1 is 2500
N1 <- n1*T 
W1 <- cell2nb(r2,r2)    
Wl1 <- nb2listw(W1)   #listw object 
Wn1 <- listw2dgCMatrix(Wl1)  #sparse matrix 
Wn1.a <- as.matrix(Wn1)
Cn1 <- runif(n1, -0.1, 0.1)   #fixed effects 
X1 <- rnorm(N1)     #non stochastic time varying regressors 
V1 <- rnorm(N1)       
Yn11 <- solve((diag(N1)-kronecker(diag(T),lambda1*Wn1)))%*%(X1*beta + rep(Cn1,T) + V1) 
#M <- M.matrix(lambda1,Yn11)

r3 <- 35
n2 <- r3*r3   #n2 is 1225
N2 <- n2*T 
W2 <- cell2nb(r3,r3)    
Wl2 <- nb2listw(W2)   #listw object 
Wn2 <- listw2dgCMatrix(Wl2)  #sparse matrix 
Wn2.a <- as.matrix(Wn2)
Cn2 <- rnorm(n2)   #fixed effects 
X2 <- rnorm(N2)     #non stochastic time varying regressors 
V2 <- rnorm(N2)       
Yn12 <- solve((diag(N2)-kronecker(diag(T),lambda1*Wn2)))%*%(X2*beta + rep(Cn2,T) + V2) 


###### negative log-likelihood function 

mylog.lik.neg<-function(theta,Yn1) { 
  lambda <- theta[1]
  sig2 <- theta[2]
  Sn <- diag(n)-lambda*as.matrix(Wn) 

  Y.tilde.nt<-matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
  V.tilde.nt<-Sn %*% Y.tilde.nt 
  ell_nt<- n*(T-1)/2*log(2*pi*sig2)-(T-1)*log(det(Sn))+0.5* sum(diag((t(V.tilde.nt) %*% (V.tilde.nt))))/sig2
  return(ell_nt) 
  
} 


# seq.lam<-seq(-0.9,0.9,0.01) 
# lik2plot<-NULL 
# 
# for (l in 1:length((seq.lam))) { 
#   
#   lik2plot[l]<- -mylog.lik(seq.lam[l],Yn1) 
# } 

##### Plot of Log-Lik 

# plot(seq.lam,lik2plot) 
# abline(v=lambda1,col=2,lwd=2) 

###### NEW OPTIM FOR MLE #################### 
f<-nlminb(start = c(0,0.1), objective = mylog.lik.neg,lower = c(-0.99, 0.001), upper= c(0.99,Inf), Yn1=Yn1)
lambda.hat <- f$par[1]
sig2.hat <- f$par[2]
set.seed(1234)
MC.size <- 1000
lambda.hat <- numeric(length = MC.size)
sig2.hat <- numeric(length = MC.size)
for (i in 1:MC.size) {
  
  Cn0 <- runif(n,-0.1,0.1)    #fixed effects 
  X <- rnorm(N)     #non stochastic time varying regressors 
  V<- rnorm(N)       
  Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 
  Y.tilde.nt<-matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
  f<-nlminb(start = c(0,0.1), objective = mylog.lik.neg,lower = c(-0.99, 0.001), upper= c(0.99,Inf), Yn1=Yn1)
  lambda.hat[i] <- f$par[1]
  sig2.hat[i] <- f$par[2]
}
hist(lambda.hat,freq = F, 
     breaks = 80, xlim = c(-1,1),ylim = c(0,3),
     main = "",
     xlab="", col="gray")

###### transform Y_{nt} ###########
# Y.tilde.nt <- function(n,Yn1){
#   Y.tilde.t<-matrix(Yn1,nrow=n,ncol=T) 
#   Ybar<-rowSums(Y.tilde.t)/T
#   Y.tilde.n <- Y.tilde.t - matrix(rep(Ybar,T),n,T) 
#   return(Y.tilde.n)
# }



###### M matrix 
# 
# M.matrix <- function(lambda,Yn11){
#   Sn <- diag(n1)-lambda*Wn1.a
# 
#   Yw <- Wn1.a%*%Y.tilde.nt(n1,Yn11)
# 
#   Gn <- Wn1.a%*%solve(Sn)
#   G2n <- Gn%*%Gn
#   a <- (T-1)*sum(diag(G2n))
#   for(t in 1:T){
#     a <- a + t(Yw)[t,] %*% Yw[,t]/sig2
#     b <- sig2
#   }
#   b <- b-(T-1)/(T*sig2)
#   return(M)
# }

M.matrix1 <- function(lambda0,sig2.hat){
  Sn <- diag(n)-lambda0*Wn.a
  Gn <- Wn.a%*%solve(Sn)
  M <- array(dim = c(T,T,n)) 
  for (i in 1:n) {
    a <- (T-1)*t(Gn[i,])%*%((Gn[i,])+Gn[,i])
    b <- 1/sig2.hat*(T-1)*Gn[i,i]
    c <- (T-1)/(2*sig2.hat*sig2.hat)
    M[,,i] <- matrix(data = c(a,b,b,c),T,T)
  }
  return(M)
}

  
M <- M.matrix1(lambda1,sig2)

# M <- M.matrix(0,Yn1)

IF <- function(lambda,sig2){
  Sn <- diag(n)-lambda*Wn.a
  Gn <- Wn.a%*%solve(Sn)
  Yw <- Wn.a%*%Y.tilde.nt
  Vn <- Sn %*% Y.tilde.nt
  IF <- matrix(nrow = n, ncol = 2)
  for (i in 1:n) {
    
    psi.T1 <-  1/(sig2)*Yw[i,]%*%Vn[i,]-(T-1)*Gn[i,i]
    psi.T2 <- Vn[i,]%*%Vn[i,]/(2*sig2*sig2)-(T-1)/(2*sig2)
    psi.T <- c(psi.T1,psi.T2)
    IF[i,]<- solve(M[,,i])%*%psi.T
    
  }
  return(IF)
}



g.T <- IF(lambda1,sig2)[,1]/2  

gamma.ij <- function(i,j,lambda,sig2){
  Sn <- diag(n)-lambda*Wn.a
  Gn <- Wn.a%*%solve(Sn)
  Gn3 <- Gn%*%Gn%*%Gn
  Vn <- Sn %*% Y.tilde.nt
  
  IF.fct <- IF(lambda,sig2)
  
  a <- matrix(c(-2*Gn3[i,i], (T-1)*t(Gn[i,])%*%Gn[i,]/sig2,(T-1)*t(Gn[i,])%*%Gn[i,]/sig2, 2*(T-1)*Gn[i,i]/((sig2)^2)),2,2)
  b <- matrix(c((T-1)*t(Gn[i,])%*%Gn[i,]/sig2,2*(T-1)*Gn[i,i]/((sig2)^2),2*(T-1)*Gn[i,i]/((sig2)^2),2*(T-1)/((sig2)^3)),2,2)
  Ga <- c(IF.fct[j,]%*%a%*%IF.fct[i,],IF.fct[j,]%*%b%*%IF.fct[i,])
  
  psi.der.i <- matrix(c(-sum((Gn[i,]%*%Vn)^2)/sig2-(T-1)*Gn[i,]%*%Gn[,i], -Gn[i,]%*%Vn%*%Vn[i,]/(sig2*sig2),-Gn[i,]%*%Vn%*%Vn[i,]/(sig2*sig2), -Vn[i,]%*%Vn[i,]/(sig2^3)+(T-1)/(2*sig2*sig2)),2,2)
  psi.der.j <- matrix(c(-sum((Gn[j,]%*%Vn)^2)/sig2-(T-1)*Gn[j,]%*%Gn[,j], -Gn[j,]%*%Vn%*%Vn[j,]/(sig2*sig2),-Gn[j,]%*%Vn%*%Vn[j,]/(sig2*sig2), -Vn[j,]%*%Vn[j,]/(sig2^3)+(T-1)/(2*sig2*sig2)),2,2)
  phi <-  IF.fct[i,]+IF.fct[j,]+solve(M[,,i])%*%(Ga+  psi.der.i%*%IF.fct[j,]+psi.der.j%*%IF.fct[i,])
  return(phi[1]/2)
}

gamma.ij(1,2,lambda1,sig2)

cumulant <- function(MC.size,lambda1,sig2.hat) {
  U <- numeric(length = MC.size)
  g2 <- numeric(length = MC.size)
  g3 <- numeric(length = MC.size)
  g4 <- numeric(length = MC.size)
  
  gamma2 <- numeric(length = MC.size)
  g1g2gamma12 <- numeric(length = MC.size)
  g12g2gamma12 <- numeric(length = MC.size)
  g1g2gamma13gamma23 <- numeric(length = MC.size)
  for (k in 1:MC.size) {
    Cn0 <- runif(n,-1,1)#fixed effects 
    # X <- rnorm(N)     #non stochastic time varying regressors 
    V<- rnorm(N)       
    Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 
    Y.tilde.nt<-matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
    # f<-nlminb(start = c(0,0.1), objective = mylog.lik.neg,lower = c(-0.99, 0.001), upper= c(0.99,Inf), Yn1=Yn1)
    # lambda.hat <- f$par[1]
    # sig2.hat <- f$par[2]
    M <- M.matrix1(lambda1,sig2.hat)
    g.T <- IF(lambda1,sig2.hat)[,1]/2 
    g2[k] <- mean(g.T^2)
    g3[k] <- mean(g.T^3)
    g4[k] <- mean(g.T^4)
    
    S1 <- 0
    S2 <- 0
    S3 <- 0
    S4 <- 0
    S5 <- 0

    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        gamma <- gamma.ij(i,j,lambda1,sig2.hat) 
        S1 <- S1 + gamma
        S2 <- S2 + gamma^2
        S3 <- S3 + g.T[i]*g.T[j]*gamma
        S5 <- S5 + g.T[i]*g.T[j]*gamma*(g.T[i]+g.T[j])

        for (l in 1:n) {
          if((l != i) && (l!=j)){
            gamm.il <- gamma.ij(i,l,lambda1,sig2.hat) 
            gamm.jl <- gamma.ij(j,l,lambda1,sig2.hat) 
            S4 <- S4 + g.T[i]*g.T[j]*gamm.il*gamm.jl
           
          }
          
        }
      }
      
    }
    
    gamma2[k] <- S2*2/(n*(n-1))
    g1g2gamma12[k] <- S3*2/(n*(n-1))
    g1g2gamma13gamma23[k] <- S4*2/(n*(n-1)*(n-2)) 
    g12g2gamma12[k] <- S5/(n*(n-1))
    U[k] <- mean(g.T)*2 + 2*S1/(n*(n-1))
    print(k)
  }
  
  
  mu <- mean(U)
  sigma2.g <- mean(g2)
  sigma.g <- sqrt(sigma2.g)
  sigma2.nT <- 4*sigma2.g/n + 2/(n*(n-1))*mean(gamma2)
  sigma3.nT <- sqrt(sigma2.nT)*sigma2.nT
  cumu3 <- (mean(g3)+3*mean(g1g2gamma12))/(sigma2.g*sigma.g)
  cumu4 <- (mean(g4)+12*mean(g1g2gamma13gamma23)+12*mean(g12g2gamma12))/(sigma2.g*sigma2.g)-3
  cumulants <- c(mu,sigma2.g,sigma2.nT,cumu3,cumu4)
  names(cumulants) <- c("mu","sigma2.g","sigma2.nT","cumu3","cumu4")
  return(cumulants)
}

cumulant1 <- function(MC.size,lambda1) {
  U <- numeric(length = MC.size)
  g2 <- numeric(length = MC.size)
  g3 <- numeric(length = MC.size)
  g4 <- numeric(length = MC.size)
  
  gamma2 <- numeric(length = MC.size)
  g1g2gamma12 <- numeric(length = MC.size)
  g12g2gamma12 <- numeric(length = MC.size)
  g1g2gamma13gamma23 <- numeric(length = MC.size)
  for (k in 1:MC.size) {
    Cn0 <- runif(n,-0.1,0.1)   #fixed effects 
    # X <- rnorm(N)     #non stochastic time varying regressors 
    V<- rnorm(N)       
    Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 
    Y.tilde.nt<-matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
    f<-nlminb(start = c(0,0.1), objective = mylog.lik.neg,lower = c(-0.99, 0.001), upper= c(0.99,Inf), Yn1=Yn1)
    # lambda.hat <- f$par[1]
    sig2.hat <- f$par[2]
    M <- M.matrix1(lambda1,sig2.hat)
    g.T <- IF(lambda1,sig2.hat)[,1]/2 
    g2[k] <- mean(g.T^2)
    g3[k] <- mean(g.T^3)
    g4[k] <- mean(g.T^4)
    
    S1 <- 0
    S2 <- 0
    S3 <- 0
    S4 <- 0
    S5 <- 0
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        gamma <- gamma.ij(i,j,lambda1,sig2.hat) 
        S1 <- S1 + gamma
        S2 <- S2 + gamma^2
        S3 <- S3 + g.T[i]*g.T[j]*gamma
        S5 <- S5 + g.T[i]*g.T[j]*gamma*(g.T[i]+g.T[j])
        for (l in 1:n) {
          if((l != i) && (l!=j)){
            gamm.il <- gamma.ij(i,l,lambda1,sig2.hat) 
            gamm.jl <- gamma.ij(j,l,lambda1,sig2.hat) 
            S4 <- S4 + g.T[i]*g.T[j]*gamm.il*gamm.jl
          }
          
        }
      }
      
    }
    
    gamma2[k] <- S2
    g1g2gamma12[k] <- S3*2/(n*(n-1))
    g1g2gamma13gamma23[k] <- S4*2/(n*(n-1)*(n-2)) 
    g12g2gamma12[k] <- S5/(n*(n-1))
    U[k] <- mean(g.T)*2 + 2*S1/(n*(n-1))
    print(k)
  }
  
  
  mu <- mean(U)
  sigma2.g <- mean(g2)
  sigma.g <- sqrt(sigma2.g)
  sigma2.nT <- 4*sigma2.g/n + 4/(n^2*(n-1)^2)*mean(gamma2)
  sigma3.nT <- sqrt(sigma2.nT)*sigma2.nT
  cumu3 <- (mean(g3)+3*mean(g1g2gamma12))/(sigma2.g*sigma.g)
  cumu4 <- (mean(g4)+12*mean(g1g2gamma13gamma23)+12*mean(g12g2gamma12))/(sigma2.g*sigma2.g)-3
  cumulants <- c(mu,sigma2.g,sigma2.nT,cumu3,cumu4)
  names(cumulants) <- c("mu","sigma2.g","sigma2.nT","cumu3","cumu4")
  return(cumulants)
}
####### true sadd density cumulants ######
ptm <- proc.time()
set.seed(1234)
MC.size <- 50

U <- numeric(length = MC.size)
g2 <- numeric(length = MC.size)
g3 <- numeric(length = MC.size)
g4 <- numeric(length = MC.size)

gamma2 <- numeric(length = MC.size)
g1g2gamma12 <- numeric(length = MC.size)
g12g2gamma12 <- numeric(length = MC.size)
g1g2gamma13gamma23 <- numeric(length = MC.size)
for (k in 1:MC.size) {
  Cn0 <- runif(n,-1,1)#fixed effects 
  # X <- rnorm(N)     #non stochastic time varying regressors 
  V<- rnorm(N)       
  Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 
  Y.tilde.nt<-matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
  # f<-nlminb(start = c(0,0.1), objective = mylog.lik.neg,lower = c(-0.99, 0.001), upper= c(0.99,Inf), Yn1=Yn1)
  # lambda.hat <- f$par[1]
  # sig2.hat <- f$par[2]
  M <- M.matrix1(lambda1,sig2)
  g.T <- IF(lambda1,sig2)[,1]/2 
  g2[k] <- mean(g.T^2)
  g3[k] <- mean(g.T^3)
  g4[k] <- mean(g.T^4)
  
  S1 <- 0
  S2 <- 0
  S3 <- 0
  S4 <- 0
  S5 <- 0
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      gamma <- gamma.ij(i,j,lambda1,sig2) 
      S1 <- S1 + gamma
      S2 <- S2 + gamma^2
      S3 <- S3 + g.T[i]*g.T[j]*gamma
      S5 <- S5 + g.T[i]*g.T[j]*gamma*(g.T[i]+g.T[j])
      
      for (l in 1:n) {
        if((l != i) && (l!=j)){
          gamm.il <- gamma.ij(i,l,lambda1,sig2) 
          gamm.jl <- gamma.ij(j,l,lambda1,sig2) 
          S4 <- S4 + g.T[i]*g.T[j]*gamm.il*gamm.jl
          
        }
        
      }
    }
    
  }
  
  gamma2[k] <- S2*2/(n*(n-1))
  g1g2gamma12[k] <- S3*2/(n*(n-1))
  g1g2gamma13gamma23[k] <- S4*2/(n*(n-1)*(n-2)) 
  g12g2gamma12[k] <- S5/(n*(n-1))
  U[k] <- mean(g.T)*2 + 2*S1/(n*(n-1))
  print(k)
}


mu <- mean(U)
sigma2.g <- mean(g2)
sigma.g <- sqrt(sigma2.g)
sigma2.nT <- 4*sigma2.g/n + 2/(n*(n-1))*mean(gamma2)
sigma3.nT <- sqrt(sigma2.nT)*sigma2.nT
cumu3 <- (mean(g3)+3*mean(g1g2gamma12))/(sigma2.g*sigma.g)
cumu4 <- (mean(g4)+12*mean(g1g2gamma13gamma23)+12*mean(g12g2gamma12))/(sigma2.g*sigma2.g)-3

###### cumulant generating function ######
cgf <- function(u){
  c <-  mu*u + 0.5* n*sigma2.nT*u*u+(1/6)*(n^(1.5))*cumu3*(sigma3.nT)*(u*u*u) +(1/24)*(n*n)*cumu4*(sigma2.nT*sigma2.nT)*(u*u*u*u)
  return(c)
}



###### The first derivative of cgf #####

Der1.cgf <- function(u){
  der1 <- mu + n*sigma2.nT*u+(1/2)*(n^(1.5))*cumu3*(sigma3.nT)*(u*u)+(1/6)*(n*n)*cumu4*(sigma2.nT*sigma2.nT)*(u*u*u)
  return(der1)
}


###### The second derivative of cgf ######

Der2.cgf <- function(u){
  der2 <- n*sigma2.nT+(n^(1.5))*cumu3*(sigma3.nT)*(u)+(1/2)*(n*n)*cumu4*(sigma2.nT*sigma2.nT)*(u*u)
  return(der2)
}

Der2.cgf.Wang <- function(u,a){
  wn <- exp(-n*sigma2.nT*a*a*u*u/2) 
  der2 <-  (n*sigma2.nT+((n^(1.5))*cumu3*(sigma3.nT)*(u)+(1/2)*(n*n)*cumu4*(sigma2.nT*sigma2.nT)*(u*u))*wn+
              2*((1/2)*(n^(1.5))*cumu3*(sigma3.nT)*(u*u)+(1/6)*(n*n)*cumu4*(sigma2.nT*sigma2.nT)*(u*u*u))*(-n*sigma2.nT*a*a*u)*wn+
              ((1/6)*(n^(1.5))*cumu3*(sigma3.nT)*(u*u*u)+(1/24)*(n*n)*cumu4*(sigma2.nT*sigma2.nT)*(u*u*u*u))*wn*(n*sigma2.nT*a*a-1)*n*sigma2.nT*a*a)
  return(der2)
}
###### find saddlepoints 

Sad <- function(a){

  sad <- uniroot(function(u) Der1.cgf(u)-a, lower = -100,upper = 100,tol = 0.0001)$root
  
  return(sad)
}

##### plot of cgf

sad.grid<-seq(-0.99,0.99,0.01)
cgf.vals<-NULL
for (i in 1:length(sad.grid)) {
  cgf.vals[i]<-cgf(sad.grid[i])
}

plot(sad.grid,cgf.vals)
plot(sad.grid,Der1.cgf(sad.grid))
plot(sad.grid,Der2.cgf(sad.grid))

#### plot of saddlepoints
sad.vals<-NULL
for (i in 1:length(sad.grid)) {
  sad.vals[i]<-Sad(sad.grid[i])
  
}

plot(sad.grid,sad.vals)


cgf(sad.vals)
Der2.cgf(sad.vals)

###### pdf ###########

p.nT <- function(a){
  p <- sqrt(n/(2*pi*Der2.cgf(Sad(a))))*exp(n*(cgf(Sad(a))-a*Sad(a)))
  return(p)
}

theta.grid<-seq(-0.99,0.99,by=0.01)
p<-NULL
for (i in 1:length(theta.grid)) {
  p[i]<-p.nT(theta.grid[i])  
} 
 
c.int<-sum(p)*diff(theta.grid)[1]
p.std<-p/c.int

#plot(theta.grid,p.std)
while (is.na(p.std)) {
  cumulants <- cumulant(MC.size,lambda1,sig2)
  mu <- cumulants[1]
  sigma2.g <- cumulants[2]
  sigma.g <- sqrt(sigma2.g)
  sigma2.nT <- cumulants[3]
  sigma3.nT <- sqrt(sigma2.nT)*sigma2.nT
  cumu3 <- cumulants[4]
  cumu4 <- cumulants[5]
  p<-NULL
  for (i in 1:length(theta.grid)) {
    p[i]<-p.nT(theta.grid[i])  
  } 
  
  c.int<-sum(p)*diff(theta.grid)[1]
  p.std<-p/c.int
}

plot(theta.grid,p.std)

####### fbplot ############
#set.seed(1234)
fb.size <- 30
theta.grid<-seq(-0.99,0.99,by=0.01)
p.fb1 <- matrix(ncol = fb.size,nrow = length(theta.grid))
MC.size <- 50
lambda.hat <- numeric(length = fb.size)
sig2.hat <- numeric(length = fb.size)
for (s in 1:fb.size) {
  
  
  Cn0 <- runif(n,-1,1)   #fixed effects 
  X <- rnorm(N)     #non stochastic time varying regressors 
  V<- rnorm(N)       
  Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 
  Y.tilde.nt<-matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
  f<-nlminb(start = c(0,0.1), objective = mylog.lik.neg,lower = c(-0.99, 0.001), upper= c(0.99,Inf), Yn1=Yn1)
  lambda.hat[s] <- f$par[1]
  sig2.hat[s] <- f$par[2]
  
  
  U <- numeric(length = MC.size)
  g2 <- numeric(length = MC.size)
  g3 <- numeric(length = MC.size)
  g4 <- numeric(length = MC.size)
  
  gamma2 <- numeric(length = MC.size)
  g1g2gamma12 <- numeric(length = MC.size)
  g12g2gamma12 <- numeric(length = MC.size)
  g1g2gamma13gamma23 <- numeric(length = MC.size)
  for (k in 1:MC.size) {
    Cn0 <- runif(n,-1,1)#fixed effects 
    # X <- rnorm(N)     #non stochastic time varying regressors 
    V<- rnorm(N)       
    Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 
    Y.tilde.nt<-matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
    # f<-nlminb(start = c(0,0.1), objective = mylog.lik.neg,lower = c(-0.99, 0.001), upper= c(0.99,Inf), Yn1=Yn1)
    # lambda.hat <- f$par[1]
    # sig2.hat <- f$par[2]
    M <- M.matrix1(lambda1,sig2.hat[s])
    g.T <- IF(lambda1,sig2.hat[s])[,1]/2 
    g2[k] <- mean(g.T^2)
    g3[k] <- mean(g.T^3)
    g4[k] <- mean(g.T^4)
    
    S1 <- 0
    S2 <- 0
    S3 <- 0
    S4 <- 0
    S5 <- 0
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        gamma <- gamma.ij(i,j,lambda1,sig2.hat[s]) 
        S1 <- S1 + gamma
        S2 <- S2 + gamma^2
        S3 <- S3 + g.T[i]*g.T[j]*gamma
        S5 <- S5 + g.T[i]*g.T[j]*gamma*(g.T[i]+g.T[j])
        
        for (l in 1:n) {
          if((l != i) && (l!=j)){
            gamm.il <- gamma.ij(i,l,lambda1,sig2.hat[s]) 
            gamm.jl <- gamma.ij(j,l,lambda1,sig2.hat[s]) 
            S4 <- S4 + g.T[i]*g.T[j]*gamm.il*gamm.jl
            
          }
          
        }
      }
      
    }
    
    gamma2[k] <- S2*2/(n*(n-1))
    g1g2gamma12[k] <- S3*2/(n*(n-1))
    g1g2gamma13gamma23[k] <- S4*2/(n*(n-1)*(n-2)) 
    g12g2gamma12[k] <- S5/(n*(n-1))
    U[k] <- mean(g.T)*2 + 2*S1/(n*(n-1))
    print(k)
  }
  
  
  mu <- mean(U)
  sigma2.g <- mean(g2)
  sigma.g <- sqrt(sigma2.g)
  sigma2.nT <- 4*sigma2.g/n + 2/(n*(n-1))*mean(gamma2)
  sigma3.nT <- sqrt(sigma2.nT)*sigma2.nT
  cumu3 <- (mean(g3)+3*mean(g1g2gamma12))/(sigma2.g*sigma.g)
  cumu4 <- (mean(g4)+12*mean(g1g2gamma13gamma23)+12*mean(g12g2gamma12))/(sigma2.g*sigma2.g)-3
  
  
  p1<-NULL
  for (i in 1:length(theta.grid)) {
    p1[i]<-p.nT(theta.grid[i])  
  } 
  c.int<-sum(p1)*diff(theta.grid)[1]
  p.fb1[,s]<-p1/c.int
  print(s)
}

proc.time() - ptm
save.image(file='n24_WnRook_Sig1.RData')

#library("fda")
lamda.hat1 <- lambda.hat

lambda.hat.c <- cbind(lambda.hat,lambda.hat1)
p.fb.c <- cbind(p.fb.c,p.fb50)

p.fb.omit <- p.fb.c[,!is.na(p.fb.c[1,])]

my.hist <- hist(lambda.hat,freq = F, 
     breaks = 80, xlim = c(-1,1),ylim = c(0,2.5),
     main = "",
     xlab="", col="gray")
plot(my.hist$breaks,
     c(my.hist$density,0)
     ,type="s",xlim=c(-0.99,0.99),ylim = c(0,2.5),lwd=2,xlab=" ", ylab="Density", main=" ",col="gray52")
lines(theta.grid,p.std,col="blue",lwd=2,ylim = c(0,2.5))

par(new=TRUE)
fbplot(p.fb.omit,method = "MBD",ylim = c(0,2.5),xaxt = 'n',color=NA,outliercol=NA,barcol="orange3",ylab=" ",xlab=" " )
