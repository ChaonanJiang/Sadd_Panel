library("spdep") 
library("splm") 
library("graphics")



r1 <- 6 
c1 <-4 
n <- r1*c1   
T <- 2 
N <- n*T 

beta <- 0.0 
lambda1 <- 0.2 
pho1 <- 0.0 
sig <- 1 

set.seed(1234) 

W <- cell2nb(r1,c1)    
Wl <- nb2listw(W)   #listw object 
Wn <- listw2dgCMatrix(Wl)  #sparse matrix
Wn.a <- as.matrix(Wn)
Cn0 <- rnorm(n)   #fixed effects 
X <- rnorm(N)     #non stochastic time varying regressors 
V<- rnorm(N)       
Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 
Y.tilde.nt<-matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)

r2 <- 50
n1 <- r2*r2   #n1 is 2500
N1 <- n1*T 
W1 <- cell2nb(r2,r2)    
Wl1 <- nb2listw(W1)   #listw object 
Wn1 <- listw2dgCMatrix(Wl1)  #sparse matrix 
Wn1.a <- as.matrix(Wn1)
Cn1 <- rnorm(n1)   #fixed effects 
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


###### log likelihood function 

mylog.lik<-function(lambda,Yn1) { 
  
  Sn <- diag(n)-lambda*as.matrix(Wn) 

  Y.tilde.nt<-matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
  V.tilde.nt<-Sn %*% Y.tilde.nt 
  ell_nt<- (T-1)*log(det(Sn))  - 0.5* sum(diag((t(V.tilde.nt) %*% (V.tilde.nt)))) 
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
  deriv1 <- sum(diag(t(Yw)%*%V.tilde.nt))-(T-1)*sum(diag(Gn)) 
  return(deriv1) 
  
} 



###### Score function S_n,t (containing the whole cross-sections) 

Score.nt <- function(t,lambda,Yn1,n,Wn.a){ 
  
  Sn <- diag(n)-lambda*Wn.a 
  V.tilde.nt<-Sn%*%Y.tilde.nt(n,Yn1)
  Yw <- Wn.a%*%Y.tilde.nt(n,Yn1) 
  Gn <- Wn.a%*%solve(Sn)

  Scor.nt <- n*(T-1)*t(Yw)[t,] %*% V.tilde.nt[,t]-n*((T-1)*(T-1))/T*sum(diag(Gn)) 
  return(Scor.nt) 
} 



###### Score function S_i,t (for each individual i at time t) 

Score.it <- function(i,t, lambda, Yn1,n,Wn.a) { 
  
  Sn <- diag(n)-lambda*Wn.a 
  
  V.tilde.nt<-Sn%*%Y.tilde.nt(n,Yn1)
  Yw <- Wn.a%*%Y.tilde.nt(n,Yn1) 
  Gn <- Wn.a%*%solve(Sn)

  Scor.it <- n*(T-1)*Yw[i,t]* V.tilde.nt[i,t] - ((T-1)*(T-1))/T*sum(diag(Gn)) 
  return(Scor.it) 
  
} 



###### The first derivative of S_i,t with respect to lambda 

der1.Score.it <- function(i,t, lambda, Yn1,n,Wn.a){ 
  
  Sn <- diag(n)-lambda*Wn.a
  V.tilde.nt<-Sn%*%Y.tilde.nt(n,Yn1)
  Yw <- Wn.a%*%Y.tilde.nt(n,Yn1) 
  Gn <- Wn.a%*%solve(Sn)
  G2n <- Gn%*%Gn  
  deriv.Scor.it <- -n*(T-1)*(Yw[i,t])*(Yw[i,t]) - (T-1)*(T-1)/T*sum(diag(G2n)) 
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
ptm <- proc.time()

M.matrix <- function(lambda,Yn11){ 
  Sn <- diag(n1)-lambda*Wn1.a 
 
  Yw <- Wn1.a%*%Y.tilde.nt(n1,Yn11) 
  Gn <- Wn1.a%*%solve(Sn)
  G2n <- Gn%*%Gn
  M <- (T-1)*sum(diag(G2n)) 
  for(t in 1:T){ 
    M <- M + t(Yw)[t,] %*% Yw[,t] 
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
  phi <- IF.iT(i,lambda,Yn1,n,Wn.a)+IF.iT(j,lambda,Yn1,n,Wn.a)+Gamma.ij(i,j,lambda,Yn1,n,Wn.a)/M +
  ((1/(T-1))*Sum.T.Der1(j,lambda,Yn1,n,Wn.a)*IF.iT(i,lambda,Yn1,n,Wn.a)+(1/(T-1))*Sum.T.Der1(i,lambda,Yn1,n,Wn.a)*IF.iT(j,lambda,Yn1,n,Wn.a))/M 
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


###### E(g_i),E(g_i^2), E(g_i^3), E(g_i^4)



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

###### E(g_igamma_ij)
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
  c <-  gam.con1[1]/n*u + 0.5* n*Sigma2.nT*u*u+(1/6)*(n^(1.5))*cum3.nt*(Sigma3.nT)*(u*u*u) +(1/24)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u*u*u)
          return(c)
}
 
  


###### The first derivative of cgf

Der1.cgf <- function(u){
  der1 <- gam.con1[1]/n + n*Sigma2.nT*u+(1/2)*(n^(1.5))*cum3.nt*(Sigma3.nT)*(u*u)+(1/6)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u*u)
  return(der1)
}


###### The second derivative of cgf

Der2.cgf <- function(u){
  der2 <- n*Sigma2.nT+(n^(1.5))*cum3.nt*(Sigma3.nT)*(u)+(1/2)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u)
  return(der2)
}

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
proc.time() - ptm
###### Asymptotic variance
Gn <- Wn.a%*%solve(diag(n)-lambda1*Wn.a)
asy.sigma2 <- sum(diag((t(Gn)+Gn)%*%Gn))/(n*n*(T-1))

######  Bootstrap
ptm <- proc.time()
set.seed(1234)
N.C <- 100
boot.times <- 1000
lambda.hat <- vector(length = boot.times*N.C)
for (j in 1:N.C) {

  Cn0 <- rnorm(n)   #fixed effects 
  X <- rnorm(N)     #non stochastic time varying regressors 
        
  
  for (i in 1:boot.times) {
    V<- rnorm(N)      
    Yn1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
    lambda.hat[i+boot.times*(j-1)]<-optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn1, maximum = T)$maximum 
  }
  #lines(density(lambda.hat-mean(lambda.hat)),lwd=2,lty=3,col="red")
   
}

proc.time() - ptm







lambda.hat2 <- matrix(lambda.hat,boot.times,N.C)
grid <- seq(-0.99,0.99,0.01)
DEN <- NULL
den.1 <- matrix(rep(0,N.C*length(grid)),length(grid),N.C)
for (j in 1:N.C) {
  DEN <- density(lambda.hat2[,j]-mean(lambda.hat2[,j]))
  for (i in 1:length(grid)) {
  den.1[i,j] <- approx(DEN$x,DEN$y,xout=grid[i])$y
  if(is.na(den.1[i,j])){den.1[i,j]=0}
  
}
}



### LOAD the MC results for MLE


Z1.plot<-(z1-mean(z1))
my.hist<-hist(Z1.plot,freq = F, breaks = 80, xlim=c(-0.99,0.99),
              main=" ",
              xlab=" ", col="gray", ylim = c(0,2.5))
plot(my.hist$breaks,
     c(my.hist$density,0)
     ,type="s",xlim=c(-0.99,0.99),ylim = c(0,2.5),lwd=2,xlab=" ", ylab="Density", main=" ",col="gray52")
par(new=TRUE)
fbplot(den.1,method = "MBD",ylim = c(0,2.5),xaxt = 'n',color=NA,outliercol=NA,barcol="orange3",ylab=" ",xlab=" " )

# theta.grid<-my.hist$breaks[2:14]
# my.hist.dens<-my.hist$density
#lines(density(lambda.hat-mean(lambda.hat)),lwd=4,lty=3,col="purple")
lines(theta.grid,p.std,typ='l',col="blue",ylim = c(0,2.5),lwd=2)
lines(theta.grid,dnorm(theta.grid,0,sqrt(asy.sigma2)),col="red",lwd=2, lty=2)
# lines(theta.grid,my.hist.dens,col="blue")
#abline(v=quantile(Z1.plot,0.95),col="magenta",lwd=3)
Z1.t.plot<-(z1.t-mean(z1.t))
my.hist<-hist(Z1.t.plot,freq = F, breaks = 80, xlim=c(-0.99,0.99),
              main=" ",
              xlab=" ", col="gray", ylim = c(0,2.5))

Z1.t6.plot<-(z1.t6-mean(z1.t6))
my.hist<-hist(Z1.t6.plot,freq = F, breaks = 80, xlim=c(-0.99,0.99),
              main=" ",
              xlab=" ", col="gray", ylim = c(0,2.5))

####################### KERNEL DENSITY

ks.den<-density(Z1.plot,bw=0.13)

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



####################### CDF of SAD, kER, BOOTSTRAP

CDF.SAD <- function(b){
  v <- Sad(b)
  c <- v*sqrt(Der2.cgf(v))
  r <- sign(v)*sqrt(2*n*(v*b-cgf(v)))
  p <- 1-pnorm(r)+dnorm(r)*(1/c-1/r)
  return(p)
}




CDF.SAD1 <- function(b){
  theta.grid<-seq(-0.99,0.99,by=0.01)
  p<-NULL
  for (i in 1:(100*b+100)) {
    p[i]<-p.nT(theta.grid[i])  
  } 
  
  c<-sum(p)*diff(theta.grid)[1]/c.int
  return(c)
  
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



CDF.BOO <- ecdf(lambda.hat-mean(lambda.hat))

####################### PP-Plot
CDF.EMP <- ecdf((z1-mean(z1)))


seq.b <- seq(-0.99,0.99,0.01)
cdf.emp <- vector(length = length(seq.b))
cdf.sad <- vector(length = length(seq.b))
cdf.ker <- vector(length = length(seq.b))
cdf.boo1 <- vector(length = length(seq.b))

for(i in 1:length(seq.b)){
  cdf.emp[i] <- CDF.EMP(seq.b[i])
  cdf.sad[i] <- CDF.SAD1(seq.b[i])
  cdf.ker[i] <- CDF.KER(seq.b[i])
  cdf.boo1[i] <- CDF.BOO(seq.b[i])
}


cdf.asy1 <- pnorm(seq.b,0,sqrt(asy.sigma2))

plot(seq.b,cdf.ker)
lines(seq.b,cdf.emp,col="purple")
plot(cdf.emp,cdf.asy1,col="red", type="l", lty=2, lwd=2,xlab = "P(Empirical)",
     ylab = "P(models)",main="")
lines(cdf.emp,cdf.sad,col="blue",lwd=2)
#lines(cdf.emp,cdf.boo1,col="purple",lwd=4,lty=3)
abline(0,1,type="l",lty=4,lwd=3)


########## Left tail

CDF.EMP <- ecdf((z1-mean(z1)))


seq.c <- seq(-0.5,-0.25,0.005)

cdf.asy <- pnorm(seq.c,0,sqrt(asy.sigma2))
rela.lug <- vector(length = length(seq.c))
rela.asy <- vector(length = length(seq.c))
rela.edg <- vector(length = length(seq.c))
rela.boo <- vector(length = length(seq.c))

for(i in 1:length(seq.c)){
  rela.lug[i] <- (CDF.SAD1(seq.c[i]))/(CDF.EMP(seq.c[i]))-1
  rela.asy[i] <- (cdf.asy[i])/(CDF.EMP(seq.c[i]))-1
  rela.edg[i] <- (CDF.EDG(seq.c[i]))/(CDF.EMP(seq.c[i]))-1
  rela.boo[i] <- (CDF.BOO(seq.c[i]))/(CDF.EMP(seq.c[i]))-1
}


smooth <- smooth.spline(seq.c,abs(rela.lug),df=16)
smooth1 <- smooth.spline(seq.c,abs(rela.edg),df=16)
smooth2 <- smooth.spline(seq.c,abs(rela.boo),df=26)
plot(smooth1,col="dark green",pch=18,  lwd=2, ylim = c(0,1.1),xlab = "z",
     ylab = "Relative error",main="")
lines(smooth,col="blue",lwd=2)
lines(seq.c,abs(rela.asy),col="red",type="l",lty=2,  lwd=2)
#lines(smooth2,lwd=4,lty=3,col="purple")

####### Right tail

CDF.EMP <- ecdf((z1-mean(z1)))


seq.c <- seq(0.275,0.575,0.005)

cdf.asy <- pnorm(seq.c,0,sqrt(asy.sigma2))
rela.lug <- vector(length = length(seq.c))
rela.asy <- vector(length = length(seq.c))
rela.edg <- vector(length = length(seq.c))
rela.boo <- vector(length = length(seq.c))

for(i in 1:length(seq.c)){
  rela.lug[i] <- (1-CDF.SAD1(seq.c[i]))/(1-CDF.EMP(seq.c[i]))-1
  rela.asy[i] <- (1-cdf.asy[i])/(1-CDF.EMP(seq.c[i]))-1
  rela.edg[i] <- (1-CDF.EDG(seq.c[i]))/(1-CDF.EMP(seq.c[i]))-1
  rela.boo[i] <- (1-CDF.BOO(seq.c[i]))/(1-CDF.EMP(seq.c[i]))-1
}


smooth <- smooth.spline(seq.c,abs(rela.lug),df=16)
smooth1 <- smooth.spline(seq.c,abs(rela.edg),df=20)
smooth2 <- smooth.spline(seq.c,abs(rela.boo),df=20)
plot(smooth1,col="dark green",pch=18,  lwd=2, ylim = c(0,1.2),xlab = "b",
     ylab = "Relative error",main="")
lines(smooth,col="blue",lwd=2)
lines(seq.c,abs(rela.asy),col="red",type="l",lty=2,  lwd=2)
#lines(smooth2,lwd=4,lty=3,col="purple")



   

########## Edgeworth Expansion

CDF.EDG <- function(x){
  p <- pnorm(x/Sigma.nT,0,1)-dnorm(x/Sigma.nT,0,1)*n^(-0.5)*cum3.nt/6*((x/Sigma.nT)*(x/Sigma.nT)-1)
  return(p)
}


cdf.edg <- vector(length=length(seq.b))
for(i in 1:length(seq.b)){
  cdf.edg[i] <- CDF.EDG(seq.b[i])

}
plot(seq.b,cdf.edg,col="dark green",pch=18,lwd=2,xlab = "",ylab = "CDF of Models",main = "")
lines(seq.b,cdf.emp,lwd=3,type = "l",lty=4)
lines(seq.b,cdf.sad,col="blue",lwd=2)
lines(seq.b,cdf.asy1,col= "red",lwd=2,type = "l",lty=2)


############### Quantile

OBJ.sad.975 <- function(b){
  o <- (CDF.SAD(b)-0.025)*(CDF.SAD(b)-0.025)
  return(o)
}

OBJ.sad.025 <- function(b){
  o <- (CDF.SAD(b)-0.975)*(CDF.SAD(b)-0.975)
  return(o)
}

OBJ.edg <- function(b,p){
  o <- (1-CDF.EDG(b)-p)*(1-CDF.EDG(b)-p)
  return(o)
}




 
