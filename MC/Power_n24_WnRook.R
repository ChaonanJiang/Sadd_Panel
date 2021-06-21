library("spdep") 
library("splm") 
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
n1 <- r2*r2   #n1 is 2300
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
n2 <- r3*r3   #n2 is 1225 
N2 <- n2*T 
W2 <- cell2nb(r3,r3)    
Wl2 <- nb2listw(W2)   #listw object 
Wn2 <- listw2dgCMatrix(Wl2)  #sparse matrix 
Wn2.a <- as.matrix(Wn2)
Cn2 <- rnorm(n2,0,sig)   #fixed effects 
X2 <- rnorm(N2,0,sig)     #non stochastic time varying regressors 
V2 <- rnorm(N2,0,sig)       
Yn12 <- solve((diag(N2)-kronecker(diag(T),lambda1*Wn2)))%*%(X2*beta + rep(Cn2,T) + V2) 


###### log-likelihood function 

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


###### transform Y_{nt} ###########
Y.tilde.nt <- function(n,Yn1){
  Y.tilde.t<-matrix(Yn1,nrow=n,ncol=T) 
  Ybar<-rowSums(Y.tilde.t)/T
  Y.tilde.n <- Y.tilde.t - matrix(rep(Ybar,T),n,T) 
  return(Y.tilde.n)
}



###### M matrix 

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


###### The expected value of g_{i,T}, g^2_{i,T}, g^3_{i,T}, g^4_{i,T}

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

###### The expected value of gamma_{i,j,T}, gamma^2_{i,j,T},      ######
###### g_{i,T}g_{j,T}gamma_{i,j,T}, g^2_{i,T}g_{j,T}gamma_{i,j,T} ######

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


###### The expected value of g_{i,T}g_{j,T}gamma_{i,k,T}gamma_{j,k,T}

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

Sigma2.nT <- (4/n)*g.powers[2]+(2/(n*(n-1)))*gam.con1[2]   #sigma^2

Sigma.nT <- sqrt(Sigma2.nT)    #sigma
Sigma3.nT <- Sigma.nT*Sigma2.nT    #sigma^3
Sigma4.nT <- Sigma2.nT*Sigma2.nT   #sigma^4





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

Sad <- function(a){
  
  sad <- uniroot(function(u) Der1.cgf(u)-a, lower = -10000,upper = 10000)$root
  
  return(sad)
}

##### plot of cgf

sad.grid<-seq(-0.99,0.99,0.01)
cgf.vals<-NULL
for (i in 1:length(sad.grid)) {
  cgf.vals[i]<-cgf(sad.grid[i])
}

plot(sad.grid,cgf.vals)

#### plot of saddlepoints
sad.vals<-NULL
for (i in 1:length(sad.grid)) {
  sad.vals[i]<-Sad(sad.grid[i])
}

plot(sad.grid,sad.vals)


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


###### Asymptotic variance
Gn <- Wn.a%*%solve(diag(n)-lambda1*Wn.a)
asy.sigma2 <- sum(diag((t(Gn)+Gn)%*%Gn))/(n*n*(T-1))




####################### CDF of SAD ###############

CDF.SAD1 <- function(b){
  theta.grid<-seq(-0.99,0.99,by=0.01)
  p<-NULL
  for (i in 1:(100*b+100)) {
    p[i]<-p.nT(theta.grid[i])  
  } 
  
  c<-sum(p)*diff(theta.grid)[1]/c.int
  return(c)
  
}    

########## Edgeworth Expansion

CDF.EDG <- function(x){
  p <- pnorm(x/Sigma.nT,0,1)-dnorm(x/Sigma.nT,0,1)*n^(-0.5)*cum3.nt/6*((x/Sigma.nT)*(x/Sigma.nT)-1)
  return(p)
}


############## one side test #####################
######## get quantiles #######
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


############ Right side test ###############

lambda2 <- 0.2
lambda <- 0.1
MC.size <- 25000
r1 <- 6
c1 <- 4

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
  Yn <- solve((diag(N)-kronecker(diag(T),lambda*Wn)))%*%(X*beta + rep(Cn0,T) + V)       #lambda=0.1
  Yn.1 <- solve((diag(N)-kronecker(diag(T),lambda2*Wn)))%*%(X*beta + rep(Cn0,T) + V)    #lambda2=0.2
  Yn.2 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)    #lambda1=0.0
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
    
    
    ####Saddlepoint
    if(z < q.sad.u[j]) x.sad.1[j] <- 0
    else x.sad.1[j] <-1     #lambda=0.1
    
    if(z.1 < q.sad.u[j]) x.sad.2[j] <- 0
    else x.sad.2[j] <-1     #lambda=0.2
    
    if(z.2 < q.sad.u[j]) x.sad[j] <- 0
    else x.sad[j] <-1       #lambda=0.0
    
    ####Edgeworth
    if(z < q.edg.u[j] ) x.edg.1[j] <- 0
    else x.edg.1[j] <-1     #lambda=0.1
    
    if( z.1 < q.edg.u[j]) x.edg.2[j] <- 0
    else x.edg.2[j] <-1     #lambda=0.2
    
    if( z.2 < q.edg.u[j] ) x.edg[j] <- 0
    else x.edg[j] <-1       #lambda=0.0
    
    
    ####Asymptotic normal
    if(z < q.asy.u[j] ) x.asy.1[j] <- 0
    else x.asy.1[j] <-1     #lambda=0.1
    
    if(z.1 < q.asy.u[j] ) x.asy.2[j] <- 0
    else x.asy.2[j] <-1     #lambda=0.2
    
    if(z.2 < q.asy.u[j]) x.asy[j] <- 0
    else x.asy[j] <-1       #lambda=0.0
    
  }
  return(cbind(x.sad,x.sad.1,x.sad.2,x.edg,x.edg.1,x.edg.2,x.asy,x.asy.1,x.asy.2))
  
  
}


Powers.sum.r <- Powers.sum.r/MC.size

stopCluster(cl)

####### estimated size ################

#############################################################################
# Figure 5: Estimated size versus nominal size between 1% and 10% under     #
# saddlepoint (continuous line), Edgeworth (dotted line with diamonds) and  # 
# first-order asymptotic approximation (dotted line).                       #
# Sample size is n = 24 and Wn is Rook, while lambda_0 = 0:0.               #
#############################################################################

smooth.sad.size <- smooth.spline(Powers.sum.r[,1],df=8)  #saddlepoint
smooth.edg.size <- smooth.spline(Powers.sum.r[,4],df=9)  #edgeworth
smooth.asy.size <- smooth.spline(Powers.sum.r[,7],df=9)  #asymptotic normal
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

