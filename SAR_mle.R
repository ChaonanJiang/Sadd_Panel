###### Figure 2: SAR(1) model: QQ-plot vs normal of the MLE, for different sample sizes and different weight matrices ######

library("spdep")
library("splm")

################# n = 24 ################

r1 <- 6 
c1 <- 4
n <- 24 

T <- 2    ### time dimension 

##### initialize parameters 

beta <- 0.0
lambda1 <- 0.2
pho1 <- 0.0
sig <- 1.0
lambda2 <- 0.0
pho2 <- 0.0
N <- n*T
MC.size<-5000   #### Monte Carlo size

##### log-likelihood function 
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
  ell_nt<- (T-1)*log(det(Sn))  - 0.5* sum(diag((t(V.tilde.nt) %*% (V.tilde.nt))))
  return(ell_nt)
  
}



################# SETTING 1 (W = rook)

set.seed(123)

z1 <- matrix(rep(0,MC.size), MC.size,1)

### MC repeat 5000 times
for(i in 1:MC.size){
  W <- cell2nb(r1,c1)   # rook type
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  
  z1[i,]= optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn1, maximum = T)$maximum #c(feml1$coefficients[1])
  print(i)
}



################# SETTING 2 (W = queen)

set.seed(123)

z1.2 <- matrix(rep(0,MC.size), MC.size,1)

### MC repeat 5000 times
for(i in 1:MC.size){
  W <- cell2nb(r1,c1,type = "queen")   # queen type
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  
#   d1 <- data.frame(id = rep(c(1:n),T),time = c(rep(1,n),rep(2,n)), Yn1@x,X)
#   feml1 <- spml(Yn1.x~X,  listw = Wl, data = d1, model = "within", spatial.error = "none",lag = T, LeeYu = T, Hess = F)
#   summary(feml1)
  
  z1.2[i,]=optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn1, maximum = T)$maximum #c(feml1$coefficients[1])
  print(i)
}




################# SETTING 3 (W = queen w torus)

set.seed(123)

z1.3 <- matrix(rep(0,MC.size), MC.size,1)

# MC repeat 1000 times
for(i in 1:MC.size){
  W <- cell2nb(r1,c1,type = "queen", torus=T)   #queen with torus type
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  
#   d1 <- data.frame(id = rep(c(1:n),T),time = c(rep(1,n),rep(2,n)), Yn1@x,X)
#   feml1 <- spml(Yn1.x~X,  listw = Wl, data = d1, model = "within", spatial.error = "none",lag = T, LeeYu = T, Hess = F)
#   summary(feml1)
  
  z1.3[i,]= optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn1, maximum = T)$maximum #c(feml1$coefficients[1])
  print(i)
}


########## QQ-plots of MLE lambda hat for n=24 ########

par(mfrow=c(1,3))

qqnorm((z1-mean(z1))/mad(z1,constant=1.44),col="blue", main ="")
abline(0,1,col = 2)
qqnorm((z1.2-mean(z1.2))/mad(z1.2,constant=0.95),col="blue", main ="")
abline(0,1,col = 2)
qqnorm((z1.3-mean(z1.3))/mad(z1.3,constant=0.95),col="blue", main ="")
abline(0,1,col = 2)




################# n = 100 ###############

r1 <- 10
c1 <- 10
n <- r1*c1

T <-2   #### time dimension

### initialize parameters
beta <- 0.0
lambda1 <- 0.2
pho1 <- 0.0
sig <- 1.0
lambda2 <- 0.0
pho2 <- 0.0
N <- n*T
MC.size<-5000

################# SETTING 1 (W = rook)

set.seed(123)

z11 <- matrix(rep(0,MC.size), MC.size,1)

# MC repeat 5000 times
for(i in 1:MC.size){
  W <- cell2nb(r1,c1)   # rook type
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  
  z11[i,]= optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn1, maximum = T)$maximum #c(feml1$coefficients[1])
  print(i)
}



################# SETTING 2 (W = queen)

set.seed(123)

z11.2 <- matrix(rep(0,MC.size), MC.size,1)

# MC repeat 5000 times
for(i in 1:MC.size){
  W <- cell2nb(r1,c1,type = "queen")   # queen type
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  
  
  z11.2[i,]= optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn1, maximum = T)$maximum #c(feml1$coefficients[1])
  print(i)
}




################# SETTING 3 (W = queen w torus)

set.seed(123)

z11.3 <- matrix(rep(0,MC.size), MC.size,1)

# MC repeat 5000 times
for(i in 1:MC.size){
  W <- cell2nb(r1,c1,type = "queen", torus=T)    # queen with torus type
  Wl <- nb2listw(W)   #listw object
  Wn <- listw2dgCMatrix(Wl)  #sparse matrix
  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N) #matrix(1,nrow=N,ncol=1)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((Diagonal(N)-kronecker(Diagonal(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + V)
  

  z11.3[i,]= optimize(mylog.lik,interval=c(-0.99,0.99), Yn1=Yn1, maximum = T)$maximum #c(feml1$coefficients[1])
  print(i)
}


########## QQ-plots of MLE lambda hat for n=100 ########


par(mfrow=c(1,3))

qqnorm((z11-mean(z11))/mad(z11,constant=1.44),col="blue",main="")
abline(0,1,col = 2)
qqnorm((z11.2-mean(z11.2))/mad(z11.2,constant=0.95),col="blue",main="")
abline(0,1,col = 2)
qqnorm((z11.3-mean(z11.3))/mad(z11.3,constant=1.04),col="blue",main="")
abline(0,1,col = 2)

