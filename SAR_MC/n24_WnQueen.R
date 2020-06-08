library("spdep") 
library("splm") 

######### LOAD the MC results for MLE, z1.2 data from SAR_mle.rdata
load("SAR_mle.rdata")


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

W <- cell2nb(r1,c1,type = "queen",torus = F)   #queen type 
Wl <- nb2listw(W)   #listw object 
Wn <- listw2dgCMatrix(Wl)  #sparse matrix
Wn.a <- as.matrix(Wn)
Cn0 <- rnorm(n)   #fixed effects 
X <- rnorm(N)     #non stochastic time varying regressors 
V<- rnorm(N)       
Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 

r2 <- 50
n1 <- r2*r2   #n is 2500
N1 <- n1*T 
W1 <- cell2nb(r2,r2,type = "queen",torus = F)    #queen type
Wl1 <- nb2listw(W1)   #listw object 
Wn1 <- listw2dgCMatrix(Wl1)  #sparse matrix 
Wn1.a <- as.matrix(Wn1)
Cn1 <- rnorm(n1)   #fixed effects 
X1 <- rnorm(N1)     #non stochastic time varying regressors 
V1 <- rnorm(N1)       
Yn11 <- solve((diag(N1)-kronecker(diag(T),lambda1*Wn1)))%*%(X1*beta + rep(Cn1,T) + V1) 


r3 <- 35
n2 <- r3*r3   #n is 1225
N2 <- n2*T 
W2 <- cell2nb(r3,r3,type = "queen",torus = F)    #queen type
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


######### transform Y_{nt} ###########

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

ptm <- proc.time()

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

Sigma.nT <- sqrt(Sigma2.nT)     #sigma
Sigma3.nT <- Sigma.nT*Sigma2.nT     #sigma^3
Sigma4.nT <- Sigma2.nT*Sigma2.nT    #sigma^4





###### The third order cumulant 

cum3.nt <-  (sqrt(g.powers[2]))^(-3)*(g.powers[3]+3*gam.con1[3]) 


###### The fourth order cumulant 

cum4.nt <- (sqrt(g.powers[2]))^(-4)*(g.powers[4]+12*gam.con1[4]+12*Gam.Con2(lambda1,Yn12,n2,Wn2.a))


###### cumulant generating function 

cgf <- function(u){
  c <- gam.con1[1]/n*u + 0.5* n*Sigma2.nT*u*u+(1/6)*(n^(1.5))*cum3.nt*(Sigma3.nT)*(u*u*u) +(1/24)*(n*n)*cum4.nt*(Sigma4.nT)*(u*u*u*u)
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

seq.b <- seq(-2,2,0.01)
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

proc.time() - ptm         

plot(theta.grid,p.std)
###### Asymptotic variance
Gn <- Wn.a%*%solve(diag(n)-lambda1*Wn.a)
asy.sigma2 <- sum(diag((t(Gn)+Gn)%*%Gn))/(n*n*(T-1))

######  Bootstrap
ptm <- proc.time()
set.seed(1234)
N.C <- 100
boot.times <- 499    #B=49,499 oor 999
lambda.hat <- vector(length = boot.times*N.C)
for (j in 1:N.C) {
  
  Cn0 <- rnorm(n)   #fixed effects 
  X <- rnorm(N)     #non stochastic time varying regressors 
         
  for (i in 1:boot.times) {
    V<- rnorm(N)  
    Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V) 
    lambda.hat[i+boot.times*(j-1)]<-optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn1, maximum = T)$maximum 
  }
  
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

#################################################################################
# Figure 4: Density plot for saddlepoint (continuous line) vs asymptotic normal #
# (dotted line) probability approximation to the exact density for the MLE,     #
# for n=24, Wn=Queen.                                                                #
#################################################################################

Z11.plot<-(z1.2-mean(z1.2))/(1)
my.hist<-hist(Z11.plot,freq = F, 
              breaks = 80, xlim = c(-1,1),ylim = c(0,3),
              main = "",
              xlab="", col="gray")
lines(theta.grid,p.std,typ='l',col="blue",ylim = c(0,2),lwd=2)
lines(theta.grid,dnorm(theta.grid,0,sqrt(asy.sigma2)),
      col="red",lwd=2, lty=2)

#################################################################################
# Figure 7: SAR(1) model.                                                       #
# Density plots for saddlepoint (continuous line) vs the functional boxplot of  #
# the parametric bootstrap probability approximation to the exact density (as   #
# expressed by the histogram and obtained using MC with size 25000) for the MLE,#
# n=24, Wn=Queen, lambda_0=0.2.                                                 #
# Zoom on the right tail. In each plot, we display the functional central curve #
# (dotted line with crosses), the first and third functional quartile (two-dash #
# lines)                                                                        #
#################################################################################


source("fbplot.R")
plot(my.hist$breaks,
     c(my.hist$density,0)
     ,type="s",xlim=c(-0.99,0.99),ylim = c(0,1.8),lwd=2,xlab=" ", ylab="Density", main=" ",col="gray52")
lines(theta.grid,p.std,typ='l',col="blue",ylim = c(0,2),lwd=2)
par(new=TRUE)
fbplot(den.1,method = "MBD",ylim = c(0,1.8),xaxt = 'n',color=NA,outliercol=NA,barcol="orange3",ylab=" ",xlab=" " )



par(new=TRUE)
fbplot(den.1,method = "MBD",ylim = c(0,4),xaxt = 'n',color=NA,outliercol=NA,barcol="orange3",ylab=" ",xlab=" " )

####################### CDF of saddlepoint approximation ###############

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


###############################################################################
# Figure 3: PP-plot for saddlepoint (continuous line) vs asymptotic normal    #
# (dotted line) probability approximation for n=24 and Wn=Rook                #
###############################################################################

CDF.EMP <- ecdf((z1.2-mean(z1.2)))

seq.b <- seq(-0.99,0.99,0.01)
cdf.emp <- vector(length = length(seq.b))
cdf.sad <- vector(length = length(seq.b))


for(i in 1:length(seq.b)){
  cdf.emp[i] <- CDF.EMP(seq.b[i])
  cdf.sad[i] <- CDF.SAD1(seq.b[i])
}


cdf.asy1 <- pnorm(seq.b,0,sqrt(asy.sigma2))

plot(cdf.emp,cdf.asy1,col="red", type="l", lty=2, lwd=2,xlab = "P(Empirical)",
     ylab = "P(models)",main="")      #asymptotic normal
lines(cdf.emp,cdf.sad,col="blue",lwd=2)      #saddlepoint
abline(0,1,type="l",lty=4,lwd=3)






###############################################################################
# Figure 5: Relative error (in absolute value) for the approximate left tail  #
# probability, as obtained using the Gaussian asymptotic theory (dotted line),#
# the Edgeworth approximation (dotted line with diamonds) and saddlepoint     #
# approximation (continuous line) for the MLE. n=24 and Wn=queen              #
###############################################################################


CDF.EMP <- ecdf((z1.2-mean(z1.2))/(1))
cdf.asy <- pnorm(seq.c,0,sqrt(asy.sigma2))

seq.c <- seq(-0.9,-0.65,0.005)

rela.sad <- vector(length = length(seq.c))
rela.asy <- vector(length = length(seq.c))
rela.edg <- vector(length = length(seq.c))


for(i in 1:length(seq.c)){
  rela.sad[i] <- (CDF.SAD1(seq.c[i]))/(CDF.EMP(seq.c[i]))-1
  rela.asy[i] <- (cdf.asy[i])/(CDF.EMP(seq.c[i]))-1
  rela.edg[i] <- (CDF.EDG(seq.c[i]))/(CDF.EMP(seq.c[i]))-1
}


smooth <- smooth.spline(seq.c,abs(rela.sad),df=16)
smooth1 <- smooth.spline(seq.c,abs(rela.edg),df=16)
smooth2 <- smooth.spline(seq.c,abs(rela.boo),df=8)
plot(smooth1,col="dark green",pch=18,  lwd=2, ylim = c(0,1.1),xlab = "z",
     ylab = "Relative error",main="")     #edgeworth
lines(smooth,col="blue",lwd=2)         #saddlepint
lines(seq.c,abs(rela.asy),col="red",type="l",lty=2,  lwd=2)   #asymptotic normal

