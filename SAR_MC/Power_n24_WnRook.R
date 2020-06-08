library("spdep") 
library("splm") 
library("microbenchmark")
library("graphics")
library("doSNOW")
library("latex2exp")

r1 <- 6 
c1 <-4 
n <- r1*c1   
T <- 2 
N <- n*T 

beta <- 0.0 
lambda1 <- 0.0 
pho1 <- 0.0 
sig <- 1 

set.seed(1234) 

W <- cell2nb(r1,c1)    
Wl <- nb2listw(W)   #listw object 
Wn <- listw2dgCMatrix(Wl)  #sparse matrix
Wn.a <- as.matrix(Wn)
Cn0 <- rnorm(n,0,sig)   #fixed effects 
X <- rnorm(N,0,sig)     #non stochastic time varying regressors 
V<- rnorm(N,0,sig)       
Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 

r2 <- 50
n1 <- r2*r2   #n is 10000 
N1 <- n1*T 
W1 <- cell2nb(r2,r2)    
Wl1 <- nb2listw(W1)   #listw object 
Wn1 <- listw2dgCMatrix(Wl1)  #sparse matrix 
Wn1.a <- as.matrix(Wn1)
Cn1 <- rnorm(n1,0,sig)   #fixed effects 
X1 <- rnorm(N1,0,sig)     #non stochastic time varying regressors 
V1 <- rnorm(N1,0,sig)       
Yn11 <- solve((diag(N1)-kronecker(diag(T),lambda1*Wn1)))%*%(X1*beta + rep(Cn1,T) + V1) 
#M <- M.matrix(lambda1,Yn11)

r3 <- 35
n2 <- r3*r3   #n is 10000 
N2 <- n2*T 
W2 <- cell2nb(r3,r3)    
Wl2 <- nb2listw(W2)   #listw object 
Wn2 <- listw2dgCMatrix(Wl2)  #sparse matrix 
Wn2.a <- as.matrix(Wn2)
Cn2 <- rnorm(n2,0,sig)   #fixed effects 
X2 <- rnorm(N2,0,sig)     #non stochastic time varying regressors 
V2 <- rnorm(N2,0,sig)       
Yn12 <- solve((diag(N2)-kronecker(diag(T),lambda1*Wn2)))%*%(X2*beta + rep(Cn2,T) + V2) 


###### log likelihood function 

mylog.lik<-function(lambda,Yn1) { 
  
  Sn <- diag(n)-lambda*as.matrix(Wn) 
  Y.tilde.t<-matrix(0,nrow=n,ncol=T) 
  index<-seq(1,n,by=1) 
  Y.tilde.t[,1]<-Yn1[index] 
  for (i in 2:T) 
  { 
    a<-(i-1)*n+1 
    b<- n*i  
    index<-seq(a,b,1) 
    Y.tilde.t[,i]<-Yn1[index]   
  } 
  
  Ybar<-apply(Y.tilde.t,1,mean) 
  Y.tilde.nt <- Y.tilde.t - matrix(rep(Ybar,T),n,T) 
  V.tilde.nt<-Sn %*% Y.tilde.nt 
  ell_nt<- (T-1)*log(det(Sn))  - 0.5* sum(diag((t(V.tilde.nt) %*% (V.tilde.nt))))/(sig*sig) 
  return(ell_nt) 
  
} 


seq.lam<-seq(-0.9,0.9,0.01) 
lik2plot<-NULL 

for (l in 1:length((seq.lam))) { 
  
  lik2plot[l]<- -mylog.lik(seq.lam[l],Yn1) 
} 

##### Plot of Log-Lik 

plot(seq.lam,lik2plot) 
abline(v=lambda1,col=2,lwd=2) 

###### NEW OPTIM FOR MLE #################### 
myMLE.lambda<-optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn1, maximum = T)$maximum 

#### check with the spml algm 
d1 <- data.frame(id = rep(c(1:n),T),time = c(rep(1,n),rep(2,n)), Yn1@x,X) 
feml1 <- spml(Yn1.x~X,  listw = Wl, data = d1, model = "within", spatial.error = "none",lag = T, LeeYu = T, Hess = F) 

R.MLE <-feml1$coefficients[1] 

check<-cbind(myMLE.lambda,R.MLE) 
check 








Y.tilde.nt <- function(n,Yn1){
  Y.tilde.t<-matrix(Yn1,nrow=n,ncol=T) 
  Ybar<-rowSums(Y.tilde.t)/T
  Y.tilde.n <- Y.tilde.t - matrix(rep(Ybar,T),n,T) 
  return(Y.tilde.n)
}

###### The first derivative of log likelihood  

der1.log.lik <- function(lambda,Yn1,n,Wn.a) { 
  
  Sn <- diag(n)-lambda*Wn.a 
  V.tilde.nt<-Sn%*%Y.tilde.nt(n,Yn1)
  Yw <- Wn.a%*%Y.tilde.nt(n,Yn1) 
  Gn <- Wn.a%*%solve(Sn)
  deriv1 <- sum(diag(t(Yw)%*%V.tilde.nt))/(sig*sig)-(T-1)*sum(diag(Gn)) 
  return(deriv1) 
  
} 



###### Score function S_n,t (containing the whole cross-sections) 

Score.nt <- function(t,lambda,Yn1,n,Wn.a){ 
  
  Sn <- diag(n)-lambda*Wn.a 
  V.tilde.nt<-Sn%*%Y.tilde.nt(n,Yn1)
  Yw <- Wn.a%*%Y.tilde.nt(n,Yn1) 
  Gn <- Wn.a%*%solve(Sn)

  Scor.nt <- n*(T-1)*t(Yw)[t,] %*% V.tilde.nt[,t]/(sig*sig)-n*((T-1)*(T-1))/T*sum(diag(Gn)) 
  return(Scor.nt) 
} 



###### Score function S_i,t (for each individual i at time t) 

Score.it <- function(i,t, lambda, Yn1,n,Wn.a) { 
  
  Sn <- diag(n)-lambda*Wn.a 
  
  V.tilde.nt<-Sn%*%Y.tilde.nt(n,Yn1)
  Yw <- Wn.a%*%Y.tilde.nt(n,Yn1) 
  Gn <- Wn.a%*%solve(Sn)

  Scor.it <- n*(T-1)*Yw[i,t]* V.tilde.nt[i,t]/(sig*sig) - ((T-1)*(T-1))/T*sum(diag(Gn)) 
  return(Scor.it) 
  
} 



###### The first derivative of S_i,t with respect to lambda 

der1.Score.it <- function(i,t, lambda, Yn1,n,Wn.a){ 
  
  Sn <- diag(n)-lambda*Wn.a
  V.tilde.nt<-Sn%*%Y.tilde.nt(n,Yn1)
  Yw <- Wn.a%*%Y.tilde.nt(n,Yn1) 
  Gn <- Wn.a%*%solve(Sn)
  G2n <- Gn%*%Gn  
  deriv.Scor.it <- -n*(T-1)*(Yw[i,t])*(Yw[i,t])/(sig*sig) - (T-1)*(T-1)/T*sum(diag(G2n)) 
  return(deriv.Scor.it) 
  
} 


###### sum_t S_i,t (sum over time dimension) 

Sum.T.Scor <- function(i,lambda,Yn1,n,Wn.a){ 
   
    Score.iT <- Score.it(i,1,lambda,Yn1,n,Wn.a)+Score.it(i,2,lambda,Yn1,n,Wn.a)   
    S <- sum(Score.iT)
    return(S)
  
  
} 


######sum_t the first derivative of S_it 
Sum.T.Der1 <- function(i,lambda,Yn1,n,Wn.a){ 
 
   Der1 <- der1.Score.it(i,1,lambda,Yn1,n,Wn.a) +der1.Score.it(i,2,lambda,Yn1,n,Wn.a)
   S <- sum(Der1)
   return(S)
 
  
} 


###### M matrix 


M.matrix <- function(lambda,Yn11){ 
  Sn <- diag(n1)-lambda*Wn1.a 
 
  Yw <- Wn1.a%*%Y.tilde.nt(n1,Yn11) 
  Gn <- Wn1.a%*%solve(Sn)
  G2n <- Gn%*%Gn
  M <- (T-1)*sum(diag(G2n)) 
  for(t in 1:T){ 
    M <- M + t(Yw)[t,] %*% Yw[,t]/(sig*sig) 
  } 
  
  return(M) 
} 



M <- M.matrix(lambda1,Yn11)


###### IF function 

IF.iT <- function(i,lambda,Yn1,n,Wn.a){ 
  IF <- (1/(T-1))*Sum.T.Scor(i,lambda,Yn1,n,Wn.a)/M
  return(IF) 
  
} 


###### The second derivative of S_it 

der2.Score.it <- function(lambda,n,Wn.a){ 
  Sn <- diag(n)-lambda*Wn.a 
  Gn <- Wn.a %*% solve(Sn) 
  Der2 <- -(T-1)*(T-1)/T*sum(diag(2*Gn%*%Gn%*%Gn)) 
  return(Der2) 
  
} 


###### Gamma function 

Gamma.ij <- function(i,j, lambda, Yn1,n,Wn.a){ 
  Sn <- diag(n)-lambda*Wn.a 
  Gn <- Wn.a %*% solve(Sn) 
  gamma <- IF.iT(j, lambda, Yn1,n,Wn.a)*(-(T-1)*sum(diag(2*Gn%*%Gn%*%Gn)))*IF.iT(i,lambda,Yn1,n,Wn.a) 
  return(gamma) 
  
}  

###### The second term of Von Mises expansion 

Phi.ij <- function(i,j,lambda,Yn1,n,Wn.a){ 
  phi <- IF.iT(i,lambda,Yn1,n,Wn.a)+IF.iT(j,lambda,Yn1,n,Wn.a)+Gamma.ij(i,j,lambda,Yn1,n,Wn.a)/M+((1/(T-1))*Sum.T.Der1(j,lambda,Yn1,n,Wn.a)*IF.iT(i,lambda,Yn1,n,Wn.a)+(1/(T-1))*Sum.T.Der1(i,lambda,Yn1,n,Wn.a)*IF.iT(j,lambda,Yn1,n,Wn.a))/M 
  return(phi) 
  
} 

###### psi function 

Psi.iT <- function(i,lambda,Yn1,n,Wn.a){ 
  
  psi <- (1/2) * IF.iT(i,lambda,Yn1,n,Wn.a) 
  return(psi) 
  
} 

###### gamma function 

gam.ij <- function(i,j,lambda,Yn1,n,Wn.a){ 
  gam <- 1/2* Phi.ij(i,j,lambda,Yn1,n,Wn.a) 
  return(gam) 
  
} 







###### Psi^2,Psi^3,Psi^4 
# Psi.powers <- function(lambda,Yn1,n,Wn.a) { 
#   Psi.mat <- matrix(rep(0,(3*n)),n,3)
#   
#   for(i in 1:n){ 
#     Psi.mat[i,1]<-(Sum.T.Scor(i,lambda,Yn1,n,Wn.a))*(Sum.T.Scor(i,lambda,Yn1,n,Wn.a))
#     Psi.mat[i,2]<- Psi.mat[i,1]*(Sum.T.Scor(i,lambda,Yn1,n,Wn.a))
#     Psi.mat[i,3]<- Psi.mat[i,2]*(Sum.T.Scor(i,lambda,Yn1,n,Wn.a))
#    
#   }
#   Psi.pwr<- colSums(Psi.mat)/n
#   a <- M*(T-1)
#   Psi.pwr[1] <- Psi.pwr[1]/(4*a*a)
#   Psi.pwr[2] <- Psi.pwr[2]/(8*a*a*a)
#   Psi.pwr[3] <- Psi.pwr[3]/(16*a*a*a*a)
#   return(Psi.pwr)}



G.powers<- function(lambda,Yn1,n,Wn.a){
  Sn <- diag(n)-lambda*Wn.a 
  V.tilde.nt<-Sn%*%Y.tilde.nt(n,Yn1)
  Yw <- Wn.a%*%Y.tilde.nt(n,Yn1) 
  Gn <- Wn.a%*%solve(Sn)
  
  Psi.mat <- matrix(rep(0,(4*n)),n,4)
  for(i in 1:n){
    ##Psi.mat[i,1]  <- Score.iT:  Sum of Score.it over T divided by (T-1)
    Psi.mat[i,1] <- n*(Yw[i,1]* V.tilde.nt[i,1]+ Yw[i,2]* V.tilde.nt[i,2]- (T-1)*Gn[i,i]) 
    
    Psi.mat[i,2]<-( Psi.mat[i,1])*( Psi.mat[i,1])
    Psi.mat[i,3]<- Psi.mat[i,2]*( Psi.mat[i,1])
    Psi.mat[i,4]<- Psi.mat[i,3]*( Psi.mat[i,1])
  }
  Psi.pwr<- colSums(Psi.mat)/n
  Psi.pwr[1] <- Psi.pwr[1]/(2*M)
  Psi.pwr[2] <- Psi.pwr[2]/(4*M*M)
  Psi.pwr[3] <- Psi.pwr[3]/(8*M*M*M)
  Psi.pwr[4] <- Psi.pwr[4]/(16*M*M*M*M)
  return(Psi.pwr)
}



g.powers <- G.powers(lambda1,Yn11,n1,Wn1.a)




#Psi.powers(lambda1,Yn11,n1,Wn1)->ccc


######  
# 
# gam <- function(lambda,Yn11,n,Wn){ 
#   
#   gam.mat <- matrix(rep(0,(3*n)),n,3)
#   for(i in 1:r2){ 
#     for(j in 1:r2){ 
#       gam.mat[(i-1)*r2+j,1] <- Phi.ij(i,j,lambda,Yn11,n,Wn)*Phi.ij(i,j,lambda,Yn11,n,Wn)
#       gam.mat[(i-1)*r2+j,2] <- Sum.T.Scor(i,lambda,Yn11,n,Wn) * Sum.T.Scor(j,lambda,Yn11,n,Wn) * Phi.ij(i,j,lambda,Yn11,n,Wn)
#       gam.mat[(i-1)*r2+j,3] <- (Sum.T.Scor(i,lambda,Yn11,n,Wn))*gam.mat[(i-1)*r2+j,2]
#        } 
#   } 
#   gam.mat<- colSums(gam.mat)/n
#   gam.mat[1] <- gam.mat[1]/4
#   gam.mat[2] <- gam.mat[2]/(8*(T-1)*(T-1)*M*M)
#   gam.mat[3] <- gam.mat[3]/(16*(T-1)*(T-1)*(T-1)*M*M*M)
#   return(gam.mat) 
# } 


Gam.Con1<- function(lambda,Yn1,n,Wn.a){
  Sn <- diag(n)-lambda*Wn.a 
  V.tilde.nt<-Sn%*%Y.tilde.nt(n,Yn1)
  Yw <- Wn.a%*%Y.tilde.nt(n,Yn1) 
  Gn <- Wn.a%*%solve(Sn)
  G2n <- Gn %*%Gn
  G3n <- G2n%*%Gn
  g <- (T-1)*sum(diag(2*G3n))
  
  Scor.iT <- vector(length = n)
  Der.Scor.iT <- vector(length = n)
  for(i in 1:n){
    Scor.iT[i] <- n*(Yw[i,1]* V.tilde.nt[i,1]+ Yw[i,2]* V.tilde.nt[i,2]- (T-1)*Gn[i,i])
    Der.Scor.iT[i] <- -n*((Yw[i,1])*(Yw[i,1])+(Yw[i,2])*(Yw[i,2]) - (T-1)*G2n[i,i])
  }
  gamm.ij <- matrix(0,n,n)
  
  for(i in 1:n){
    for(j in 1:n){
      gamm.ij[i,j]<-Scor.iT[i]+ Scor.iT[j]-Scor.iT[j]*n*(T-1)*2*G3n[i,i]* Scor.iT[i]/(M*M)+Der.Scor.iT[j]* Scor.iT[i]/M+Der.Scor.iT[i]* Scor.iT[j]/M
    }
  }
  
  gam.mat <- matrix(rep(0,(4*n*n)),n*n,4)    
  for(i in 1:n){
    for(j in 1:n){
      gam.mat[(i-1)*n+j,1] <- gamm.ij[i,j]
      gam.mat[(i-1)*n+j,2] <- gamm.ij[i,j]*gamm.ij[i,j]     ##/4M2
      gam.mat[(i-1)*n+j,3] <- Scor.iT[i]*Scor.iT[j]*gamm.ij[i,j]    ##/8M3
      gam.mat[(i-1)*n+j,4] <- gam.mat[(i-1)*n+j,2]*Scor.iT[i]     ##/16M4
    }
  }
  
  gam.pwr<- colSums(gam.mat)/(n*n)
  gam.pwr[1] <- gam.pwr[1]/(2*M)
  gam.pwr[2] <- gam.pwr[2]/(4*M*M)
  gam.pwr[3] <- gam.pwr[3]/(8*M*M*M)
  gam.pwr[4] <- gam.pwr[4]/(16*M*M*M*M)
  
  # gam2.s <- 0
  # for(i in 1:n){
  #   for(j in 1:n){
  #     for(k in 1:n){
  #       gam2.s <- gam2.s +   Scor.iT[i]*Scor.iT[j]*gamm.ij[i,k]*gamm.ij[j,k]
  #     }
  #   }
  # }
  # gam2 <- gam2.s/(16*n*n*n*M*M*M*M)
  # 
  return(gam.pwr)
  
}

gam.con1 <- Gam.Con1(lambda1,Yn11,n1,Wn1.a)


###### The expected value of psi_i psi_j gamma_ik gamma_jk 

Gam.Con2 <- function(lambda,Yn1,n,Wn.a){
  Sn <- diag(n)-lambda*Wn.a 
  V.tilde.nt<-Sn%*%Y.tilde.nt(n,Yn1)
  Yw <- Wn.a%*%Y.tilde.nt(n,Yn1) 
  Gn <- Wn.a%*%solve(Sn)
  G2n <- Gn %*%Gn
  G3n <- G2n%*%Gn
  g <- (T-1)*sum(diag(2*G3n))
  
  Scor.iT <- vector(length = n)
  Der.Scor.iT <- vector(length = n)
  for(i in 1:n){
    Scor.iT[i] <- n1*(Yw[i,1]* V.tilde.nt[i,1]+ Yw[i,2]* V.tilde.nt[i,2])- (T-1)*sum(diag(Gn))
    Der.Scor.iT[i] <- -n1*((Yw[i,1])*(Yw[i,1])+(Yw[i,2])*(Yw[i,2]) )- (T-1)*sum(diag(G2n))
  }
  gamm.ij <- matrix(0,n,n)
  
  for(i in 1:n){
    for(j in 1:n){
      gamm.ij[i,j]<-  Scor.iT[i]+ Scor.iT[j]-  Scor.iT[j]*n1*(T-1)*2*G3n[i,i]* Scor.iT[i]/(M*M)+Der.Scor.iT[j]* Scor.iT[i]/M+Der.Scor.iT[i]* Scor.iT[j]/M
    }
  }
  
  gam2.s <- 0
  for(i in 1:n){
    for(j in 1:n){
      for(k in 1:n){
        gam2.s <- gam2.s +   Scor.iT[i]*Scor.iT[j]*gamm.ij[i,k]*gamm.ij[j,k]
      }
    }
  }
  gam2 <- gam2.s/(16*n*n*n*M*M*M*M)
  
  return(gam2)
  
}


###### Total variance of U-stats 

Sigma2.nT <- (4/n)*g.powers[2]+(2/(n*(n-1)))*gam.con1[2] 

Sigma.nT <- sqrt(Sigma2.nT)
Sigma3.nT <- Sigma.nT*Sigma2.nT
Sigma4.nT <- Sigma2.nT*Sigma2.nT





###### The third order cumulant 

cum3.nt <-  (sqrt(g.powers[2]))^(-3)*(g.powers[3]+3*gam.con1[3]) 


###### The fourth order cumulant 

cum4.nt <- (sqrt(g.powers[2]))^(-4)*(g.powers[4]+12*gam.con1[4]+12*Gam.Con2(lambda1,Yn12,n2,Wn2.a))
 
  




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
  der1 <- gam.con1[1]/n + n*Sigma2.nT*u +(1/2)*(n^(1.5))*cum3.nt*(Sigma3.nT)*(u*u) +(1/6)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u*u)
  return(der1)
}

seq.b <- seq(-0.99,0.99,0.01)
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

Der3.cgf <- function(u){
  der3 <-  (n^(1.5))*cum3.nt*(Sigma3.nT)+(n*n)*cum4.nt*(Sigma4.nT)*(u)
  
  return(der3)
}

seq.b <- seq(-99,99,1)
der2.cgf <- NULL
for(i in 1:length(seq.b)){
  der2.cgf[i] <- Der2.cgf(seq.b[i])
}

plot(seq.b,der2.cgf)
###### find saddlepoints


# Sad <- function(a){
#   coefficients <- c((gam.con1[1]/n-a),(n*Sigma2.nT),((1/2)*(n^(1.5))*cum3.nt*(Sigma3.nT)),((1/6)*(n*n)*cum4.nt*(Sigma4.nT)))
#   sad <- polyroot(coefficients)
#   return(sad)
# }


Sad <- function(a){


  sad <- uniroot(function(u) Der1.cgf(u)-a, lower = -10000,upper = 10000)$root
  
  return(sad)
}

sad.grid<-seq(-0.99,0.99,0.01)
cgf.vals<-NULL
for (i in 1:length(sad.grid)) {
  cgf.vals[i]<-cgf(sad.grid[i])
}

plot(sad.grid,cgf.vals)

sad.grid<-seq(-0.99,0.99,0.01)
sad.vals<-NULL
for (i in 1:length(sad.grid)) {
  sad.vals[i]<-Sad(sad.grid[i])
}

plot(sad.grid,sad.vals)

h <- function(a){
  h <- cgf(Sad(a))-a*Sad(a)
  return(h)
}

###### pdf


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
plot(theta.grid,p.std)


###### Asymptotic variance
Gn <- Wn.a%*%solve(diag(n)-lambda1*Wn.a)
asy.sigma2 <- sum(diag((t(Gn)+Gn)%*%Gn))/(n*n*(T-1))



### LOAD the MC results for MLE

Z1.plot<-(z1-mean(z1))
my.hist<-hist(Z1.plot,freq = F, breaks = 80, xlim=c(-1,1),
              main="",
              xlab="", col="gray", ylim = c(0,2.5))

# theta.grid<-my.hist$breaks[2:14]
# my.hist.dens<-my.hist$density

plot(theta.grid,p.std,typ='l',col="blue",ylim = c(0,2.5),lwd=2)
lines(theta.grid,dnorm(theta.grid,0,sqrt(asy.sigma2)),col="red",lwd=2, lty=2)
# lines(theta.grid,my.hist.dens,col="blue")
#abline(v=quantile(Z1.plot,0.95),col="magenta",lwd=3)


####################### KERNEL DENSITY

ks.den<-density(Z1.plot,bw=0.15)

lines(ks.den$x,ks.den$y,typ="l",ylim=c(0,3))

theta.grid<-ks.den$x
p<-NULL
for (i in 1:length(theta.grid)) {
  p[i]<-p.nT(theta.grid[i])  
} 


c.int<-sum(p)*diff(theta.grid)[1]
p.std<-p/c.int

lines(theta.grid,p.std,typ='l',col="blue",ylim = c(0,2.5))
lines(theta.grid,dnorm(theta.grid,0,sqrt(asy.sigma2)),col="red")



####################### CDF of SAD, kER

CDF.SAD <- function(b){
  v <- Sad(b)
  c <- v*sqrt(Der2.cgf(v))
  r <- sign(v)*sqrt(2*n*(v*b-cgf(v)))
  p <- 1-pnorm(r)+dnorm(r)*(1/c-1/r)
  return(p)
}




CDF.SAD1 <- function(b){
  theta.grid<-seq(-0.999,0.999,by=0.001)
  p<-NULL
  for (i in 1:(1000*b+1000)) {
    p[i]<-p.nT(theta.grid[i])  
  } 
  
  c<-sum(p)*diff(theta.grid)[1]/c.int
  return(c)
  
}


Sur.Lap <- function(a){
  lap <- sqrt(2*pi/(-n*(1-2/Der2.cgf(a)+2*Der3.cgf(a)*(a-Der1.cgf(a))/(Der2.cgf(a)*Der2.cgf(a)*Der2.cgf(a)))))*p.nT(a)
  
  return(lap)
  
}




CDF.KER <- function(b){
  
  p<-NULL
  for (i in 1:length(ks.den$x)) {
    if(ks.den$x[i]<= b)
      p[i]<-ks.den$y[i] 
    else 
      break
      
  } 
  
  c<-sum(p)*diff(ks.den$x)[1]
  return(c)
  
}


seq.b <- seq(0.3,0.99,0.01)
lug <- vector(length = length(seq.b))
sad1 <- vector(length = length(seq.b))
lap <- vector(length = length(seq.b))

for(i in 1:length(seq.b)){
  lug[i] <- CDF.SAD(seq.b[i])
  sad1[i] <- 1-CDF.SAD1(seq.b[i])
  lap[i] <- Sur.Lap(seq.b[i])
}





plot(seq.b,lap,col="coral3", type="l", lty=5, lwd=2,xlab = "a",
     ylab = "P",main="Survival Function")
lines(seq.b,sad1,lwd=2,lty=4,col="olivedrab")
lines(seq.b,lug,col="blue",lwd=2)
legend(0.75,0.35,c("laplace","numerical intergration","LR"),lty=c(5,4,1),col=c("coral3","olivedrab","blue"),lwd = c(2,2,2))
abline(0.05,0)


####################### PP-Plot
CDF.EMP <- ecdf((z1-mean(z1)))


seq.b <- seq(-0.99,0.99,0.01)
cdf.emp <- vector(length = length(seq.b))
cdf.sad <- vector(length = length(seq.b))
cdf.ker <- vector(length = length(seq.b))

for(i in 1:length(seq.b)){
  cdf.emp[i] <- CDF.EMP(seq.b[i])
  cdf.sad[i] <- CDF.SAD1(seq.b[i])
  cdf.ker[i] <- CDF.KER(seq.b[i])
}


cdf.asy1 <- pnorm(seq.b,0,sqrt(asy.sigma2))

plot(seq.b,cdf.ker)
lines(seq.b,cdf.emp,col="purple")
plot(cdf.emp,cdf.asy1,col="red", type="l", lty=2, lwd=2,xlab = "P(Empirical)",
     ylab = "P(models)",main="")
lines(cdf.emp,cdf.sad,col="blue",lwd=2)
abline(0,1,type="l",lty=4,lwd=3)


########## Left tail

CDF.EMP <- ecdf((z1-mean(z1)))


seq.c <- seq(-0.78,-0.57,0.005)

cdf.asy <- pnorm(seq.c,0,sqrt(asy.sigma2))
rela.lug <- vector(length = length(seq.c))
rela.asy <- vector(length = length(seq.c))
rela.edg <- vector(length = length(seq.c))

for(i in 1:length(seq.c)){
  rela.lug[i] <- (1-CDF.SAD(seq.c[i]))/(CDF.EMP(seq.c[i]))-1
  rela.asy[i] <- (cdf.asy[i])/(CDF.EMP(seq.c[i]))-1
  rela.edg[i] <- (CDF.EDG(seq.c[i]))/(CDF.EMP(seq.c[i]))-1
}


smooth <- smooth.spline(seq.c,abs(rela.lug),df=16)
smooth1 <- smooth.spline(seq.c,abs(rela.edg),df=16)
plot(seq.c,abs(rela.asy),col="red", type="l", lty=2, lwd=2, ylim = c(0,1),xlab = "b",
     ylab = "Relative error",main="")
lines(smooth,col="blue",lwd=2)
lines(smooth1,col="orange",type="l", lty=3,  lwd=2)
# lines(seq.c,abs(rela.lug),col="blue", lwd=2)


####### Right tail

CDF.EMP <- ecdf((z1-mean(z1)))


seq.c <- seq(0.5,0.70,0.005)

cdf.asy <- pnorm(seq.c,0,sqrt(asy.sigma2))
rela.lug <- vector(length = length(seq.c))
rela.asy <- vector(length = length(seq.c))
rela.edg <- vector(length = length(seq.c))

for(i in 1:length(seq.c)){
  rela.lug[i] <- (CDF.SAD(seq.c[i]))/(1-CDF.KER(seq.c[i]))-1
  rela.asy[i] <- (1-cdf.asy[i])/(1-CDF.KER(seq.c[i]))-1
  rela.edg[i] <- (1-CDF.EDG(seq.c[i]))/(1-CDF.KER(seq.c[i]))-1

}


smooth <- smooth.spline(seq.c,abs(rela.lug),df=16)
smooth1 <- smooth.spline(seq.c,abs(rela.edg),df=16)
plot(smooth1,col="dark green",pch=18,  lwd=2, ylim = c(0,1),xlab = "b",
     ylab = "Relative error",main="")
lines(smooth,col="blue",lwd=2)
lines(seq.c,abs(rela.asy),col="red",type="l",lty=2,  lwd=2)
# lines(seq.c,abs(rela.lug),col="blue", lwd=2)



# ############# Edgeworth Expansion
# 
# 
# Gn <- Wn.a%*%(diag(n)-lambda1*Wn.a)
# Gn2 <- Gn%*%Gn
# Gn3 <- Gn2%*%Gn
# trGn <- sum(diag(Gn))
# trGn2.s <- sum(diag(Gn%*%Gn+crossprod(Gn))) 
# a  <- (T-1)*(trGn2.s-2/n*(trGn*trGn))
# f  <- function(u){
#   fvalue <- a^(-3/2)*(T-1)/3*(8*trGn*trGn*trGn/(n*n)-6*trGn*trGn2.s/n+sum(diag(Gn3+3*Gn2%*%t(Gn)))+
#                           (sum(diag(2*Gn3+3*crossprod(Gn,Gn2)))-3*trGn*(sum(diag(Gn2))+trGn2.s)
#                            +4*trGn*trGn*trGn/(n*n))*u*u)
#   return(fvalue)
#   }
# 
# 
# ######### CDF of lambda-lambda1
# CDF.EGE <- function(x){
#   p <- pnorm(sqrt(a)*x,0,1)+f(sqrt(a)*x)*dnorm(sqrt(a)*x,0,1)
#   return(p)
# }





########## Edgeworth Expansion

CDF.EDG <- function(x){
  p <- pnorm(x/Sigma.nT,0,1)+dnorm(x/Sigma.nT,0,1)*n^(-0.5)*cum3.nt/6*((x/Sigma.nT)*(x/Sigma.nT)-1)
  return(p)
}


cdf.edg <- vector(length=length(seq.b))
for(i in 1:length(seq.b)){
  cdf.edg[i] <- CDF.EGE(seq.b[i])

}
plot(seq.b,cdf.edg,col="dark green",pch=18,lwd=2,xlab = "",ylab = "CDF of Models",main = "")
lines(seq.b,cdf.emp,lwd=3,type = "l",lty=4)
lines(seq.b,cdf.sad,col="blue",lwd=2)
lines(seq.b,cdf.asy1,col= "red",lwd=2,type = "l",lty=2)

#################### Two side test

########### Quantile

OBJ.sad.975 <- function(b){
  o <- (CDF.SAD(b)-0.025)*(CDF.SAD(b)-0.025)
  return(o)
}

OBJ.sad.025 <- function(b){
  o <- (CDF.SAD(b)-0.975)*(CDF.SAD(b)-0.975)
  return(o)
}

OBJ.sad.95 <- function(b){
  o <- (CDF.SAD(b)-0.05)*(CDF.SAD(b)-0.05)
  return(o)
}

OBJ.sad.05 <- function(b){
  o <- (CDF.SAD(b)-0.95)*(CDF.SAD(b)-0.95)
  return(o)
}

OBJ.edg.025 <- function(b){
  o <- (CDF.EDG(b)-0.025)*(CDF.EDG(b)-0.025)
  return(o)
}

OBJ.edg.975 <- function(b){
  o <- (CDF.EDG(b)-0.975)*(CDF.EDG(b)-0.975)
  return(o)
}

OBJ.edg.05 <- function(b){
  o <- (CDF.EDG(b)-0.05)*(CDF.EDG(b)-0.05)
  return(o)
}

OBJ.edg.95 <- function(b){
  o <- (CDF.EDG(b)-0.95)*(CDF.EDG(b)-0.95)
  return(o)
}


######## Quantile of SAD at 5% 10% siginificant level

q.sad.975 <- optimize(OBJ.sad.975,c(0.2,1))$minimum
q.sad.025 <- optimize(OBJ.sad.025,c(-1,-0.2))$minimum
q.sad.95 <- optimize(OBJ.sad.95,c(0.2,1))$minimum
q.sad.05 <- optimize(OBJ.sad.05,c(-1,-0.2))$minimum
q.sad.10 <- optimize(function(b) (CDF.SAD(b)-0.9)*(CDF.SAD(b)-0.9), c(-1,-0.2))$minimum



######## Quantile of EDG at 5% 10%

q.edg.975 <- optimize(OBJ.edg.975,c(0,1))$minimum
q.edg.025 <- optimize(OBJ.edg.025,c(-1,0))$minimum
q.edg.95 <- optimize(OBJ.edg.95,c(0,1))$minimum
q.edg.05 <- optimize(OBJ.edg.05,c(-1,0))$minimum
q.edg.10 <- optimize(function(b) (CDF.EDG(b)-0.10)*(CDF.EDG(b)-0.10),c(-1,0))$minimum


######## Quantile of asyptotic normal 5% 10%
q.asy.975 <- qnorm(0.975,0,sqrt(asy.sigma2))
q.asy.025 <- qnorm(0.025,0,sqrt(asy.sigma2))
q.asy.95 <- qnorm(0.95,0,sqrt(asy.sigma2))
q.asy.05 <- qnorm(0.05,0,sqrt(asy.sigma2))
q.asy.10 <- qnorm(0.10,0,sqrt(asy.sigma2))

lambda2 <- 0.2
lambda <- 0.1
MC.size <- 25000
r1 <- 6
c1 <- 4

cl <- makeCluster(3)

registerDoSNOW(cl)

set.seed(123)

Powers <- foreach(i = 1:MC.size,.combine='+')%dopar%{
  library("spdep")
  library("splm")
  W <- cell2nb(r1,c1)   
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  Yn.1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda2*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  Yn.2 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  z <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn, maximum = T)$maximum 
  z.1 <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn.1, maximum = T)$maximum
  z.2 <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn.2, maximum = T)$maximum
  
 
    
    ####SAD LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if( z > q.sad.05) sad.5 <- 0
    else sad.5<-1
    
    if(z.1 > q.sad.05) sad1.5 <- 0
    else sad1.5 <-1
    
    if(z.2 > q.sad.05) sad2.5 <- 0
    else sad2.5 <-1
  
    if( z > q.sad.10) sad.10 <- 0
    else sad.10 <-1
  
    if(z.1 > q.sad.10) sad1.10 <- 0
    else sad1.10 <-1
  
    if(z.2 > q.sad.10) sad2.10 <- 0
    else sad2.10 <-1
  
    ####EDG LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
  if( z > q.edg.05) edg.5 <- 0
  else edg.5<-1
  
  if(z.1 > q.edg.05) edg1.5 <- 0
  else edg1.5 <-1
  
  if(z.2 > q.edg.05) edg2.5 <- 0
  else edg2.5 <-1
  
  if( z > q.edg.10) edg.10 <- 0
  else edg.10 <-1
  
  if(z.1 > q.edg.10) edg1.10 <- 0
  else edg1.10 <-1
  
  if(z.2 > q.edg.10) edg2.10 <- 0
  else edg2.10 <-1
    
    
    ####Asy LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
  if( z > q.asy.05) asy.5 <- 0
  else asy.5<-1
  
  if(z.1 > q.asy.05) asy1.5 <- 0
  else asy1.5 <-1
  
  if(z.2 > q.asy.05) asy2.5 <- 0
  else asy2.5 <-1
  
  if( z > q.asy.10) asy.10 <- 0
  else asy.10 <-1
  
  if(z.1 > q.asy.10) asy1.10 <- 0
  else asy1.10 <-1
  
  if(z.2 > q.asy.10) asy2.10 <- 0
  else asy2.10 <-1
  
  return(c(sad.5,sad.10,sad1.5,sad1.10,sad2.5,sad2.10,edg.5,edg.10,edg1.5,edg1.10,edg2.5,edg2.10,
               asy.5,asy.10,asy1.5,asy1.10,asy2.5,asy2.10))
  
  
}


Powers <- Powers/MC.size
names(Powers) <- c("sad.5","sad.10","sad1.5","sad1.10","sad2.5","sad2.10","edg.5","edg.10","edg1.5","edg1.10","edg2.5","edg2.10",
        "asy.5","asy.10","asy1.5","asy1.10","asy2.5","asy2.10")

stopCluster(cl)

 








seq.d <- seq(0.005,0.1,0.005)
q.sad.u <- vector(length = length(seq.d))
q.sad.l <- vector(length = length(seq.d))
q.edg.u <- vector(length = length(seq.d))
q.edg.l <- vector(length = length(seq.d))
q.asy.u <- vector(length = length(seq.d))
q.asy.l <- vector(length = length(seq.d))
for(i in 1:length(seq.d)){
  q.sad.u[i] <- optimize(function(b) (CDF.SAD(b)-seq.d[i]/2)*(CDF.SAD(b)-seq.d[i]/2), lower = 0,upper = 1)$minimum
  q.sad.l[i] <- optimize(function(b) (CDF.SAD(b)+seq.d[i]/2-1)*(CDF.SAD(b)+seq.d[i]/2-1), lower = -1,upper = 0)$minimum
  q.edg.u[i] <- optimize(function(b) (CDF.EDG(b)+seq.d[i]/2-1)*(CDF.EDG(b)+seq.d[i]/2-1), lower = 0,upper = 1)$minimum
  q.edg.l[i] <- optimize(function(b) (CDF.EDG(b)-seq.d[i]/2)*(CDF.EDG(b)-seq.d[i]/2), lower = -1,upper = 0)$minimum
  q.asy.u[i] <- qnorm((1-seq.d[i]/2),0,sqrt(asy.sigma2))
  q.asy.l[i] <- qnorm((seq.d[i]/2),0,sqrt(asy.sigma2))
  
}

lambda2 <- 0.2
lambda <- 0.1
MC.size <- 25000
r1 <- 6
c1 <- 4


cl <- makeCluster(3)

registerDoSNOW(cl)

set.seed(123)

Powers.sum <- foreach(i = 1:MC.size,.combine='+')%dopar%{
  library("spdep")
  library("splm")
  W <- cell2nb(r1,c1)   
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  Yn.1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda2*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  Yn.2 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  z <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn, maximum = T)$maximum #c(feml1$coefficients[1])
  z.1 <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn.1, maximum = T)$maximum
  z.2 <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn.2, maximum = T)$maximum
  
  x.sad <- vector(length = length(seq.d))
  x.sad.1 <- vector(length = length(seq.d))
  x.sad.2 <- vector(length = length(seq.d))
  x.edg <- vector(length = length(seq.d))
  x.edg.1 <- vector(length = length(seq.d))
  x.edg.2 <- vector(length = length(seq.d))
  x.asy <- vector(length = length(seq.d))
  x.asy.1 <- vector(length = length(seq.d))
  x.asy.2 <- vector(length = length(seq.d))
  
  
  for(j in 1:length(seq.d)){
    
  
  ####SAD LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
  if( (q.sad.u[j] > z)&(z > q.sad.l[j])) x.sad.1[j] <- 0
  else x.sad.1[j] <-1
  
  if((q.sad.u[j] > z.1)&(z.1 > q.sad.l[j])) x.sad.2[j] <- 0
  else x.sad.2[j] <-1
 
  if((q.sad.u[j] > z.2)&(z.2 > q.sad.l[j])) x.sad[j] <- 0
  else x.sad[j] <-1
 
  ####EDG LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
  if((q.edg.u[j] > z)&(z> q.edg.l[j]) ) x.edg.1[j] <- 0
  else x.edg.1[j] <-1
 
  if( (q.edg.u[j] > z.1)&(z.1> q.edg.l[j]) ) x.edg.2[j] <- 0
  else x.edg.2[j] <-1
  
  if( (q.edg.u[j] > z.2)&(z.2> q.edg.l[j]) ) x.edg[j] <- 0
  else x.edg[j] <-1

  
  ####Asy LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
  if((q.asy.u[j] > z)&(z> q.asy.l[j]) ) x.asy.1[j] <- 0
  else x.asy.1[j] <-1

  if((q.asy.u[j] > z.1)&(z.1> q.asy.l[j]) ) x.asy.2[j] <- 0
  else x.asy.2[j] <-1
  
  if((q.asy.u[j] > z.2)&(z.2> q.asy.l[j]) ) x.asy[j] <- 0
  else x.asy[j] <-1
  
  }
  return(cbind(x.sad,x.sad.1,x.sad.2,x.edg,x.edg.1,x.edg.2,x.asy,x.asy.1,x.asy.2))
  
  
}


Powers.sum <- Powers.sum/MC.size

stopCluster(cl)

plot(seq.d,Powers.sum[,8],type = "l",lty=2,lwd=2,col="red")
lines(seq.d,Powers.sum[,5],lwd=2,col="dark green")
 
Delta.sad.1 <- vector(length=length(seq.d))
Delta.sad.2 <- vector(length=length(seq.d))
Delta.edg.1 <- vector(length=length(seq.d))
Delta.edg.2 <- vector(length=length(seq.d))
Delta.asy.1 <- vector(length=length(seq.d))
Delta.asy.2 <- vector(length=length(seq.d))

for(i in 1:length(seq.d)){
  Delta.sad.1[i] <- qnorm(Powers.sum[i,2],0,1)-qnorm(Powers.sum[i,1],0,1)
  Delta.sad.2[i] <- qnorm(Powers.sum[i,3],0,1)-qnorm(Powers.sum[i,1],0,1)
  Delta.edg.1[i] <- qnorm(Powers.sum[i,5],0,1)-qnorm(Powers.sum[i,4],0,1)
  Delta.edg.2[i] <- qnorm(Powers.sum[i,6],0,1)-qnorm(Powers.sum[i,4],0,1)
  Delta.asy.1[i] <- qnorm(Powers.sum[i,8],0,1)-qnorm(Powers.sum[i,7],0,1)
  Delta.asy.2[i] <- qnorm(Powers.sum[i,9],0,1)-qnorm(Powers.sum[i,7],0,1)
}

plot(seq.d,Delta.edg.1,pch=18,col="dark green",ylim=c(0.37,0.45),lwd=2,xlab = "size",ylab = "Estimate of Delta",main = "")
lines(seq.d,Delta.asy.1,type = "l",lty=2,col="red",lwd=2)
lines(seq.d,Delta.sad.1,col="blue",lwd=2)  

plot(seq.d,Delta.edg.2,pch=18,col="dark green",lwd=2,ylim=c(0.3,0.75),xlab = "size",ylab = "Estimate of Delta",main = "")
lines(seq.d,Delta.asy.2,type = "l",lty=2,col="red",lwd=2)
lines(seq.d,Delta.sad.2,col="blue",lwd=2) 




############## one side test  LEFT

seq.d <- seq(0.01,0.1,0.005)
q.sad.u <- vector(length = length(seq.d))
q.sad.l <- vector(length = length(seq.d))
q.edg.u <- vector(length = length(seq.d))
q.edg.l <- vector(length = length(seq.d))
q.asy.u <- vector(length = length(seq.d))
q.asy.l <- vector(length = length(seq.d))
for(i in 1:length(seq.d)){
  q.sad.u[i] <- optimize(function(b) (CDF.SAD1(b)+seq.d[i]-1)*(CDF.SAD1(b)+seq.d[i]-1), lower = 0,upper = 1)$minimum
  q.sad.l[i] <- optimize(function(b) (CDF.SAD1(b)-seq.d[i])*(CDF.SAD1(b)-seq.d[i]), lower = -1,upper = 0)$minimum
  q.edg.u[i] <- optimize(function(b) (CDF.EDG(b)+seq.d[i]-1)*(CDF.EDG(b)+seq.d[i]-1), lower = 0,upper = 1)$minimum
  q.edg.l[i] <- optimize(function(b) (CDF.EDG(b)-seq.d[i])*(CDF.EDG(b)-seq.d[i]), lower = -1,upper = 0)$minimum
  q.asy.u[i] <- qnorm((1-seq.d[i]),0,sqrt(asy.sigma2))
  q.asy.l[i] <- qnorm((seq.d[i]),0,sqrt(asy.sigma2))
  
}

lambda2 <- 0.2
lambda <- 0.1
MC.size <- 25000
r1 <- 6
c1 <- 4


cl <- makeCluster(3)

registerDoSNOW(cl)

set.seed(123)

Powers.sum.l <- foreach(i = 1:MC.size,.combine='+')%dopar%{
  library("spdep")
  library("splm")
  W <- cell2nb(r1,c1)   
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn <- solve((diag(N)-kronecker(diag(T),lambda*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  Yn.1 <- solve((diag(N)-kronecker(diag(T),lambda2*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  Yn.2 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  z <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn, maximum = T)$maximum #c(feml1$coefficients[1])
  z.1 <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn.1, maximum = T)$maximum
  z.2 <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn.2, maximum = T)$maximum
  
  x.sad <- vector(length = length(seq.d))
  x.sad.1 <- vector(length = length(seq.d))
  x.sad.2 <- vector(length = length(seq.d))
  x.edg <- vector(length = length(seq.d))
  x.edg.1 <- vector(length = length(seq.d))
  x.edg.2 <- vector(length = length(seq.d))
  x.asy <- vector(length = length(seq.d))
  x.asy.1 <- vector(length = length(seq.d))
  x.asy.2 <- vector(length = length(seq.d))
  
  
  for(j in 1:length(seq.d)){
    
    
    ####SAD LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z > q.sad.l[j]) x.sad.1[j] <- 0
    else x.sad.1[j] <-1
    
    if(z.1 > q.sad.l[j]) x.sad.2[j] <- 0
    else x.sad.2[j] <-1
    
    if(z.2 > q.sad.l[j]) x.sad[j] <- 0
    else x.sad[j] <-1
    
    ####EDG LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z> q.edg.l[j] ) x.edg.1[j] <- 0
    else x.edg.1[j] <-1
    
    if( z.1> q.edg.l[j]) x.edg.2[j] <- 0
    else x.edg.2[j] <-1
    
    if( z.2> q.edg.l[j] ) x.edg[j] <- 0
    else x.edg[j] <-1
    
    
    ####Asy LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z> q.asy.l[j] ) x.asy.1[j] <- 0
    else x.asy.1[j] <-1
    
    if(z.1> q.asy.l[j] ) x.asy.2[j] <- 0
    else x.asy.2[j] <-1
    
    if(z.2> q.asy.l[j]) x.asy[j] <- 0
    else x.asy[j] <-1
    
  }
  return(cbind(x.sad,x.sad.1,x.sad.2,x.edg,x.edg.1,x.edg.2,x.asy,x.asy.1,x.asy.2))
  
  
}


Powers.sum.l <- Powers.sum.l/MC.size

stopCluster(cl)

plot(seq.d,Powers.sum.l[,8],type = "l",lty=2,lwd=2,col="red")
lines(seq.d,Powers.sum.l[,5],lwd=2,col="dark green")
lines(seq.d,Powers.sum.l[,2],col="blue",lwd=2)

Delta.sad.1l <- vector(length=length(seq.d))
Delta.sad.2l <- vector(length=length(seq.d))
Delta.edg.1l <- vector(length=length(seq.d))
Delta.edg.2l <- vector(length=length(seq.d))
Delta.asy.1l <- vector(length=length(seq.d))
Delta.asy.2l <- vector(length=length(seq.d))

for(i in 1:length(seq.d)){
  Delta.sad.1l[i] <- -qnorm(Powers.sum.l[i,2],0,1)+qnorm(Powers.sum.l[i,1],0,1)
  Delta.sad.2l[i] <- -qnorm(Powers.sum.l[i,3],0,1)+qnorm(Powers.sum.l[i,1],0,1)
  Delta.edg.1l[i] <- -qnorm(Powers.sum.l[i,5],0,1)+qnorm(Powers.sum.l[i,4],0,1)
  Delta.edg.2l[i] <- -qnorm(Powers.sum.l[i,6],0,1)+qnorm(Powers.sum.l[i,4],0,1)
  Delta.asy.1l[i] <- -qnorm(Powers.sum.l[i,8],0,1)+qnorm(Powers.sum.l[i,7],0,1)
  Delta.asy.2l[i] <- -qnorm(Powers.sum.l[i,9],0,1)+qnorm(Powers.sum.l[i,7],0,1)
}

plot(seq.d,Delta.edg.1l,pch=18,col="dark green",ylim=c(0.37,0.53),lwd=2,xlab = "size",ylab = "Estimate of Delta",main = "")
lines(seq.d,Delta.asy.1l,type = "l",lty=2,col="red",lwd=2)
lines(seq.d,Delta.sad.1l,col="blue",lwd=2)  

plot(seq.d,Delta.edg.2l,pch=18,col="dark green",lwd=2,ylim=c(0.7,0.82),xlab = "size",ylab = "Estimate of Delta",main = "")
lines(seq.d,Delta.asy.2l,type = "l",lty=2,col="red",lwd=2)
lines(seq.d,Delta.sad.2l,col="blue",lwd=2) 

############ Right side test


cl <- makeCluster(3)

registerDoSNOW(cl)

set.seed(123)

Powers.sum.r <- foreach(i = 1:MC.size,.combine='+')%dopar%{
  library("spdep")
  library("splm")
  W <- cell2nb(r1,c1)   
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n,0,sig)   #fixed effects
  X <- rnorm(N,0,sig) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N,0,sig)      
  Yn <- solve((diag(N)-kronecker(diag(T),lambda*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  Yn.1 <- solve((diag(N)-kronecker(diag(T),lambda2*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  Yn.2 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  z <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn, maximum = T)$maximum #c(feml1$coefficients[1])
  z.1 <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn.1, maximum = T)$maximum
  z.2 <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn.2, maximum = T)$maximum
  
  x.sad <- vector(length = length(seq.d))
  x.sad.1 <- vector(length = length(seq.d))
  x.sad.2 <- vector(length = length(seq.d))
  x.edg <- vector(length = length(seq.d))
  x.edg.1 <- vector(length = length(seq.d))
  x.edg.2 <- vector(length = length(seq.d))
  x.asy <- vector(length = length(seq.d))
  x.asy.1 <- vector(length = length(seq.d))
  x.asy.2 <- vector(length = length(seq.d))
  
  
  for(j in 1:length(seq.d)){
    
    
    ####SAD LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z < q.sad.u[j]) x.sad.1[j] <- 0
    else x.sad.1[j] <-1
    
    if(z.1 < q.sad.u[j]) x.sad.2[j] <- 0
    else x.sad.2[j] <-1
    
    if(z.2 < q.sad.u[j]) x.sad[j] <- 0
    else x.sad[j] <-1
    
    ####EDG LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z < q.edg.u[j] ) x.edg.1[j] <- 0
    else x.edg.1[j] <-1
    
    if( z.1 < q.edg.u[j]) x.edg.2[j] <- 0
    else x.edg.2[j] <-1
    
    if( z.2 < q.edg.u[j] ) x.edg[j] <- 0
    else x.edg[j] <-1
    
    
    ####Asy LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z < q.asy.u[j] ) x.asy.1[j] <- 0
    else x.asy.1[j] <-1
    
    if(z.1 < q.asy.u[j] ) x.asy.2[j] <- 0
    else x.asy.2[j] <-1
    
    if(z.2 < q.asy.u[j]) x.asy[j] <- 0
    else x.asy[j] <-1
    
  }
  return(cbind(x.sad,x.sad.1,x.sad.2,x.edg,x.edg.1,x.edg.2,x.asy,x.asy.1,x.asy.2))
  
  
}


Powers.sum.r <- Powers.sum.r/MC.size

stopCluster(cl)

####### estimated size 

smooth.sad.size <- smooth.spline(Powers.sum.r[,1],df=8)
smooth.edg.size <- smooth.spline(Powers.sum.r[,4],df=9)
smooth.asy.size <- smooth.spline(Powers.sum.r[,7],df=9)
plot(smooth.edg.size,ylim=c(0,0.18),pch=18,lwd=2,col="dark green", axes = F, frame.plot = TRUE,ylab=TeX('Empirical size $\\hat{\\alpha}$'),xlab=TeX('nominal size $\\alpha$'))
lines(smooth.asy.size,type = "l",lty=2,lwd=2,col="red")
lines(smooth.sad.size,col="blue",lwd=2)
lines(seq(1,19,1),seq.d,type="l",lty=4,lwd=2)
axis(1, at = seq(1,19,1), labels = seq.d)
axis(2,at = seq(0,0.18,0.02),labels= seq(0,0.18,0.02))
segments(9,0,9,Powers.sum.r[9,7],lty=3,col = "grey60")
segments(0,Powers.sum.r[9,7],9,Powers.sum.r[9,7],lty=3,col = "grey60")
segments(0,Powers.sum.r[9,4],9,Powers.sum.r[9,4],lty=3,col = "grey60")
segments(0,Powers.sum.r[9,1],9,Powers.sum.r[9,1],lty=3,col = "grey60")
points(9,Powers.sum.r[9,7],pch=16,cex=1)
text(9,Powers.sum.r[9,7],labels = "Asymptotic",cex=0.8,pos = 3)
points(9,Powers.sum.r[9,4],pch=16,cex=1)
text(9,Powers.sum.r[9,4],labels = "Edgeworth",cex=0.8,pos = 3)
points(9,Powers.sum.r[9,1],pch=16,cex=1)
text(9,Powers.sum.r[9,1],labels = "Saddlepoint",cex=0.8,pos = 3)

####### Empirical Power 
###lambda=0.1

smooth.sad.pow <- smooth.spline(Powers.sum.r[,2],df=8)
smooth.edg.pow <- smooth.spline(Powers.sum.r[,5],df=9)
smooth.asy.pow <- smooth.spline(Powers.sum.r[,8],df=9)
plot(smooth.edg.pow,pch=18,ylim=c(0,0.31),lwd=2,col="dark green", axes = F, frame.plot = TRUE,ylab=TeX('Empirical power $\\hat{\\beta}$'),xlab=TeX('nominal size $\\alpha$'))
lines(smooth.asy.pow,type = "l",lty=2,lwd=2,col="red")
lines(smooth.sad.pow,col="blue",lwd=2)
axis(1, at = seq(1,19,1), labels = seq.d)
axis(2,at = seq(0,0.30,0.05),labels= seq(0,0.30,0.05))

### lambda =0.2
smooth.sad.pow1 <- smooth.spline(Powers.sum.r[,3],df=8)
smooth.edg.pow1<- smooth.spline(Powers.sum.r[,6],df=9)
smooth.asy.pow1 <- smooth.spline(Powers.sum.r[,9],df=9)
plot(smooth.edg.pow1,ylim=c(0,0.46),pch=18,lwd=2,col="dark green", axes = F, frame.plot = TRUE,ylab=TeX('Empirical power $\\hat{\\beta}$'),xlab=TeX('nominal size $\\alpha$'))
lines(smooth.asy.pow1,type = "l",lty=2,lwd=2,col="red")
lines(smooth.sad.pow1,col="blue",lwd=2)
axis(1, at = seq(1,19,1), labels = seq.d)
axis(2,at = seq(0,0.46,0.05),labels= seq(0,0.46,0.05))

plot(seq.d,Powers.sum.r[,8],type = "l",lty=2,lwd=2,col="red")
lines(seq.d,Powers.sum.r[,5],lwd=2,col="dark green")
lines(seq.d,Powers.sum.r[,2],col="blue",lwd=2)

Delta.sad.1r <- vector(length=length(seq.d))
Delta.sad.2r <- vector(length=length(seq.d))
Delta.edg.1r <- vector(length=length(seq.d))
Delta.edg.2r <- vector(length=length(seq.d))
Delta.asy.1r <- vector(length=length(seq.d))
Delta.asy.2r <- vector(length=length(seq.d)) 

for(i in 1:length(seq.d)){
  Delta.sad.1r[i] <- qnorm(Powers.sum.r[i,2],0,1)-qnorm(Powers.sum.r[i,1],0,1)
  Delta.sad.2r[i] <- qnorm(Powers.sum.r[i,3],0,1)-qnorm(Powers.sum.r[i,1],0,1)
  Delta.edg.1r[i] <- qnorm(Powers.sum.r[i,5],0,1)-qnorm(Powers.sum.r[i,4],0,1)
  Delta.edg.2r[i] <- qnorm(Powers.sum.r[i,6],0,1)-qnorm(Powers.sum.r[i,4],0,1)
  Delta.asy.1r[i] <- qnorm(Powers.sum.r[i,8],0,1)-qnorm(Powers.sum.r[i,7],0,1)
  Delta.asy.2r[i] <- qnorm(Powers.sum.r[i,9],0,1)-qnorm(Powers.sum.r[i,7],0,1)
}

smooth.sad.1r <- smooth.spline(seq.d,Delta.sad.1r,df=6)
smooth.edg.1r <- smooth.spline(seq.d,Delta.edg.1r,df=7)
smooth.asy.1r <- smooth.spline(seq.d,Delta.asy.1r,df=8)
plot(smooth.edg.1r,pch=18,col="dark green",ylim=c(0.38,0.47),lwd=2,xlab = "nominal size",ylab = "Estimate of delta",main = "")
lines(smooth.asy.1r,type = "l",lty=2,col="red",lwd=2)
lines(smooth.sad.1r,col="blue",lwd=2)  
smooth.sad.2r <- smooth.spline(seq.d,Delta.sad.2r,df=8)
smooth.edg.2r <- smooth.spline(seq.d,Delta.edg.2r,df=8)
smooth.asy.2r <- smooth.spline(seq.d,Delta.asy.2r,df=8)
plot(smooth.edg.2r,pch=18,col="dark green",lwd=2,ylim=c(0.78,0.945),xlab = "nominal size",ylab = "Estimate of delta",main = "")
lines(smooth.asy.2r,type = "l",lty=2,col="red",lwd=2)
lines(smooth.sad.2r,col="blue",lwd=2)  


Roc.sad.1r <- vector(length=length(seq.d))
Roc.sad.2r <- vector(length=length(seq.d))
Roc.edg.1r <- vector(length=length(seq.d))
Roc.edg.2r <- vector(length=length(seq.d))
Roc.asy.1r <- vector(length=length(seq.d))
Roc.asy.2r <- vector(length=length(seq.d)) 

for(i in 1:length(seq.d)){
  Roc.sad.1r[i] <- pnorm(Delta.sad.1r[i]+qnorm(0.05,0,1),0,1)
  Roc.sad.2r[i] <- pnorm(Delta.sad.2r[i]+qnorm(0.1,0,1),0,1)
  Roc.edg.1r[i] <- pnorm(Delta.edg.1r[i]+qnorm(0.05,0,1),0,1)
  Roc.edg.2r[i] <- pnorm(Delta.edg.2r[i]+qnorm(0.1,0,1),0,1)
  Roc.asy.1r[i] <- pnorm(Delta.asy.1r[i]+qnorm(0.05,0,1),0,1)
  Roc.asy.2r[i] <- pnorm(Delta.asy.2r[i]+qnorm(0.1,0,1),0,1)
}

s.Roc.sad.1r <- smooth.spline(seq.d,Roc.sad.1r,df=7)
s.Roc.edg.1r <- smooth.spline(seq.d,Roc.edg.1r,df=8)
s.Roc.asy.1r <- smooth.spline(seq.d,Roc.asy.1r,df=8)
plot(s.Roc.edg.1r,pch=18,lwd=2,col="dark green",ylim=c(0.102,0.120),xlab = "nominal size",ylab = "R(nominal size)",main = "")
lines(s.Roc.asy.1r,type = "l",lty=2,lwd=2,col="red")
lines(s.Roc.sad.1r,col="blue",lwd=2)

s.Roc.sad.2r <- smooth.spline(seq.d,Roc.sad.2r,df=9)
s.Roc.edg.2r <- smooth.spline(seq.d,Roc.edg.2r,df=9)
s.Roc.asy.2r <- smooth.spline(seq.d,Roc.asy.2r,df=9)
plot(s.Roc.edg.2r,col="dark green",pch=18,ylim=c(0.308,0.371),lwd=2,xlab = "nominal size",ylab = "R(nominal size)",main = "")
lines(s.Roc.asy.2r,type = "l",lty=2,lwd=2,col="red")
lines(s.Roc.sad.2r,col="blue",lwd=2)


####### wrong part

#### fixed alpha

q.sad.a <- vector(length = length(seq.d))
q.edg.a <- vector(length = length(seq.d))
q.asy.a <- vector(length = length(seq.d))

for(i in 1:length(seq.d)){
  q.sad.a[i] <- optimize(function(b) (CDF.SAD1(b)+Powers.sum.r[i,1]-1)*(CDF.SAD1(b)+Powers.sum.r[i,1]-1), lower = 0,upper = 1)$minimum
  
  q.edg.a[i] <- optimize(function(b) (CDF.EDG(b)+Powers.sum.r[i,4]-1)*(CDF.EDG(b)+Powers.sum.r[i,4]-1), lower = 0,upper = 1)$minimum
  
  q.asy.a[i] <- qnorm((1-Powers.sum.r[i,7]),0,sqrt(asy.sigma2))
  
  
}


cl <- makeCluster(3)

registerDoSNOW(cl)

set.seed(123)

Powers.sum.a <- foreach(i = 1:MC.size,.combine='+')%dopar%{
  library("spdep")
  library("splm")
  W <- cell2nb(r1,c1)   
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n,0,sig)   #fixed effects
  X <- rnorm(N,0,sig) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N,0,sig)      
  Yn <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  Yn.1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda2*Wn)))%*%(X*beta + rep(Cn0,T) + V)
 
  z <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn, maximum = T)$maximum #c(feml1$coefficients[1])
  z.1 <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn.1, maximum = T)$maximum
 
  
  x.sad.1 <- vector(length = length(seq.d))
  x.sad.2 <- vector(length = length(seq.d))
  
  x.edg.1 <- vector(length = length(seq.d))
  x.edg.2 <- vector(length = length(seq.d))
  
  x.asy.1 <- vector(length = length(seq.d))
  x.asy.2 <- vector(length = length(seq.d))
  
  
  
  for(j in 1:length(seq.d)){
    
    
    ####SAD LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z < q.sad.a[j]) x.sad.1[j] <- 0
    else x.sad.1[j] <-1
    
    if(z.1 < q.sad.a[j]) x.sad.2[j] <- 0
    else x.sad.2[j] <-1
    
    
    ####EDG LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z < q.edg.a[j] ) x.edg.1[j] <- 0
    else x.edg.1[j] <-1
    
    if( z.1 < q.edg.a[j]) x.edg.2[j] <- 0
    else x.edg.2[j] <-1
    
  
    
    
    ####Asy LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z < q.asy.a[j] ) x.asy.1[j] <- 0
    else x.asy.1[j] <-1
    
    if(z.1 < q.asy.a[j] ) x.asy.2[j] <- 0
    else x.asy.2[j] <-1
    
   
    
  }
  return(cbind(x.sad.1,x.sad.2,x.edg.1,x.edg.2,x.asy.1,x.asy.2))
  
  
}


Powers.sum.a <- Powers.sum.a/MC.size

stopCluster(cl)

plot(seq.d,Powers.sum.a[,3],col="dark green",pch=18,lwd=2,ylim=c(0,0.40),xlab = "nominal size",ylab = "size adjusted power",main = "")
lines(seq.d,Powers.sum.a[,5],type = "l",lty=2,lwd=2,col="red")
lines(seq.d,Powers.sum.a[,1],col="blue",lwd=2)

plot(seq.d,Powers.sum.a[,4],col="dark green",pch=18,lwd=2,ylim=c(0,0.55),xlab = "nominal size",ylab = "size adjusted power",main = "")
lines(seq.d,Powers.sum.a[,6],type = "l",lty=2,lwd=2,col="red")
lines(seq.d,Powers.sum.a[,2],col="blue",lwd=2)

##### fixed hat alpha
cl <- makeCluster(3)

registerDoSNOW(cl)

set.seed(123)

Sur.hat <- foreach(i = 1:MC.size,.combine='+')%dopar%{
  library("spdep")
  library("splm")
  W <- cell2nb(r1,c1)   
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n,0,sig)   #fixed effects
  X <- rnorm(N,0,sig) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N,0,sig)      
  Yn <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  z <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn, maximum = T)$maximum #c(feml1$coefficients[1])

  quan <- seq(0.28,0.49,0.001)
  
  x.emp <- vector(length = length(quan))
 
  
  
  
  for(j in 1:length(quan)){
    
    
    ####SAD LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z < quan[j]) x.emp[j] <- 0
    else x.emp[j] <-1
    
    
    
    
  }
  return(x.emp)
  
  
}


Sur.hat <- Sur.hat/MC.size

stopCluster(cl)

alpha <- seq(0.01,0.1,0.01)
quan <- seq(0.28,0.49,0.001)
Quan.hat <- vector(length= length(alpha))
for(i in 1:length(alpha)){
  Quan.hat[i] <- quan[which.min(abs(Sur.hat-alpha[i]))]
}

cl <- makeCluster(3)

registerDoSNOW(cl)

set.seed(123)

Powers.wo <- foreach(i = 1:MC.size,.combine='+')%dopar%{
  library("spdep")
  library("splm")
  W <- cell2nb(r1,c1)   
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n,0,sig)   #fixed effects
  X <- rnorm(N,0,sig) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N,0,sig)      
  Yn <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  Yn.1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda2*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  
  z <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn, maximum = T)$maximum #c(feml1$coefficients[1])
  z.1 <-  optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn.1, maximum = T)$maximum
  
  
  x.sad.1 <- vector(length = length(alpha))
  x.sad.2 <- vector(length = length(alpha))
  
  x.edg.1 <- vector(length = length(alpha))
  x.edg.2 <- vector(length = length(alpha))
  
  x.asy.1 <- vector(length = length(alpha))
  x.asy.2 <- vector(length = length(alpha))
  
  
  
  for(j in 1:length(alpha)){
    
    
    ####SAD LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z < Quan.hat[j]) x.sad.1[j] <- 0
    else x.sad.1[j] <-1
    
    if(z.1 < Quan.hat[j]) x.sad.2[j] <- 0
    else x.sad.2[j] <-1
    
    
    ####EDG LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z < Quan.hat[j] ) x.edg.1[j] <- 0
    else x.edg.1[j] <-1
    
    if( z.1 < Quan.hat[j]) x.edg.2[j] <- 0
    else x.edg.2[j] <-1
    
    
    
    
    ####Asy LAMBDA=0.1, 5% 10%  LAMBDA=0.2, 5% 10%
    if(z <Quan.hat[j] ) x.asy.1[j] <- 0
    else x.asy.1[j] <-1
    
    if(z.1 < Quan.hat[j] ) x.asy.2[j] <- 0
    else x.asy.2[j] <-1
    
    
    
  }
  return(cbind(x.sad.1,x.sad.2,x.edg.1,x.edg.2,x.asy.1,x.asy.2))
  
  
}


Powers.wo <- Powers.wo/MC.size

stopCluster(cl)

plot(alpha,Powers.wo[,3],col="dark green",pch=18,lwd=3,xlab = "empirical size",ylab = "size adjusted power",main = "")
lines(alpha,Powers.wo[,5],type = "l",lty=2,lwd=3,col="red")
lines(alpha,Powers.wo[,1],col="blue",lwd=2)

plot(alpha,Powers.wo[,4],col="dark green",pch=18,lwd=3,xlab = "empirical size",ylab = "size adjusted power",main = "")
lines(alpha,Powers.wo[,6],type = "l",lty=2,lwd=3,col="red")
lines(alpha,Powers.wo[,2],col="blue",lwd=2)