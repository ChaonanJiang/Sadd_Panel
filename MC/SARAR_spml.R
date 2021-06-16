###### Figure 2: SARAR(1,1) model: QQ-plot vs normal of the MLE lambda hat, for different sample sizes and different weight matrices ######
library("spdep")
library("splm")


################## n = 24 #############


r1 <- 6
c1 <-4
n <- r1*c1  
T <- 5

#define theta1 and theta2
beta <- 1.0
lambda1 <- 0.2
pho1 <- 0.5
sig <- 1.0
lambda2 <- 0.5
pho2 <- 0.2
N <- n*T
MC.size<-5000




############## SETTING 1.1 (W is rook)
ptm <- proc.time()
set.seed(123)
W <- cell2nb(r1,c1)   
Wl <- nb2listw(W)   #listw object
Wn <- as.matrix(listw2dgCMatrix(Wl))  #sparse matrix

z1 <- matrix(rep(0,MC.size*4), MC.size,4)

# MC repeat 1000 times
for(i in 1:MC.size){

Cn0 <- rnorm(n)   #fixed effects
X <- rnorm(N)     #non stochastic time varying regressors
V<- rnorm(N)      
Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + solve((diag(N)-kronecker(diag(T),pho1*Wn)))%*%V)

d1 <- data.frame(id = rep(c(1:n),T), time =kronecker(1:T,rep(1,n)),Yn1,X)
feml1 <- spml(Yn1~X,  listw = Wl, data = d1, model = "within", spatial.error = "b",lag = T, LeeYu = T, Hess = F)
summary(feml1)



z1[i,]=c(feml1$coefficients,feml1$sigma2)
print(i)

}

proc.time() - ptm
############## SETTING 1.2 (W is queen)

set.seed(123)
W <- cell2nb(r1,c1,type = "queen",torus = F)   
Wl <- nb2listw(W)   #listw object
Wn <- as.matrix(listw2dgCMatrix(Wl))  #sparse matrix
z2 <- matrix(rep(0,MC.size*4), MC.size,4)
# MC repeat 1000 times
for(i in 1:MC.size){

  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + solve((diag(N)-kronecker(diag(T),pho1*Wn)))%*%V)
  
  d1 <- data.frame(id = rep(c(1:n),T), time =kronecker(1:T,rep(1,n)),Yn1,X)
  feml1 <- spml(Yn1~X,  listw = Wl, data = d1, model = "within", spatial.error = "b",lag = T, LeeYu = T, Hess = F)
  summary(feml1)
  
  
  
  z2[i,]=c(feml1$coefficients,feml1$sigma2)
  
  
}
 

############## SETTING 1.3 (W is queen with torus)

set.seed(123)
W <- cell2nb(r1,c1,type = "queen",torus = T)   
Wl <- nb2listw(W)   #listw object
Wn <- as.matrix(listw2dgCMatrix(Wl))  #sparse matrix
z3 <- matrix(rep(0,MC.size*4), MC.size,4)
# MC repeat 1000 times
for(i in 1:MC.size){

  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + solve((diag(N)-kronecker(diag(T),pho1*Wn)))%*%V)
  
  d1 <- data.frame(id = rep(c(1:n),T), time =kronecker(1:T,rep(1,n)),Yn1,X)
  feml1 <- spml(Yn1~X,  listw = Wl, data = d1, model = "within", spatial.error = "b",lag = T, LeeYu = T, Hess = F)
  summary(feml1)
  
  
  
  z3[i,]=c(feml1$coefficients,feml1$sigma2)
  
  
}


################## n = 100 #############



############## SETTING 2.1 (W is rook)
r2 <- 10
n<- r2*r2
N <- n*T

W <- cell2nb(r2,r2)   
Wl <- nb2listw(W)   #listw object
Wn <- as.matrix(listw2dgCMatrix(Wl)) #sparse matrix

z11 <- matrix(rep(0,MC.size*4), MC.size,4)
set.seed(123)
# MC repeat 1000 times
for(i in 1:MC.size){

  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + solve((diag(N)-kronecker(diag(T),pho1*Wn)))%*%V)
  
  d1 <- data.frame(id = rep(c(1:n),T), time =kronecker(1:T,rep(1,n)),Yn1,X)
  feml1 <- spml(Yn1~X,  listw = Wl, data = d1, model = "within", spatial.error = "b",lag = T, LeeYu = T, Hess = F)
  summary(feml1)
  
  
  
  z11[i,]=c(feml1$coefficients,feml1$sigma2)
  
  
}


############## SETTING 2.2 (W is queen)

set.seed(123)
W <- cell2nb(r2,r2,type = "queen",torus = F)   
Wl <- nb2listw(W)   #listw object
Wn <- as.matrix(listw2dgCMatrix(Wl)) #sparse matrix

z21 <- matrix(rep(0,MC.size*4), MC.size,4)
# MC repeat 1000 times
for(i in 1:MC.size){

  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + solve((diag(N)-kronecker(diag(T),pho1*Wn)))%*%V)
  
  d1 <- data.frame(id = rep(c(1:n),T), time =kronecker(1:T,rep(1,n)),Yn1,X)
  feml1 <- spml(Yn1~X,  listw = Wl, data = d1, model = "within", spatial.error = "b",lag = T, LeeYu = T, Hess = F)
  summary(feml1)
  
  
  
  z21[i,]=c(feml1$coefficients,feml1$sigma2)
  
  
}


############## SETTING 2.2 (W is queen with torus)

set.seed(123)

W <- cell2nb(r2,r2,type = "queen",torus = T)   
Wl <- nb2listw(W)   #listw object
Wn <- as.matrix(listw2dgCMatrix(Wl)) #sparse matrix

z31 <- matrix(rep(0,MC.size*4), MC.size,4)
# MC repeat 1000 times
for(i in 1:MC.size){

  Cn0 <- rnorm(n)   #fixed effects
  X <- rnorm(N)     #non stochastic time varying regressors
  V<- rnorm(N)      
  Yn1 <- solve((diag(N)-kronecker(diag(T),lambda1*Wn)))%*%(X*beta + rep(Cn0,T) + solve((diag(N)-kronecker(diag(T),pho1*Wn)))%*%V)
  
  d1 <- data.frame(id = rep(c(1:n),T), time =kronecker(1:T,rep(1,n)),Yn1,X)
  feml1 <- spml(Yn1~X,  listw = Wl, data = d1, model = "within", spatial.error = "b",lag = T, LeeYu = T, Hess = F)
  summary(feml1)
  
  
  
  z31[i,]=c(feml1$coefficients,feml1$sigma2)
  
  
}

############## PLOTS



par(mfrow=c(1,3))
qqnorm((z1[,1]-mean(z1[,1]))*sqrt(24),col="blue", main="",ylim=c(-4,4))
abline(0,1,col = 2)
qqnorm((z2[,1]-mean(z2[,1]))*sqrt(24),col="blue", main="",ylim=c(-6,5))
abline(0,1,col = 2)
qqnorm((z3[,1]-mean(z3[,1]))*sqrt(24),col="blue", main="",ylim=c(-6.3,5.2))
abline(0,1,col = 2)


par(mfrow=c(1,3))
qqnorm((z11[,1]-mean(z11[,1]))*sqrt(100),col="blue", main="",ylim = c(-4,4))
abline(0,1,col = 2)
qqnorm((z21[,1]-mean(z21[,1]))*sqrt(100),col="blue", main="",ylim = c(-6,5))
abline(0,1,col = 2)
qqnorm((z31[,1]-mean(z31[,1]))*sqrt(100),col="blue", main="",ylim = c(-6.3,5.2))
abline(0,1,col = 2)


save.image(file='spml.RData')