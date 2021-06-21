library("spdep") 
library("splm") 
library("graphics")
library("numDeriv")
library("R.matlab")

rm(list=ls())


##################################################

rho<-0.5

r1 <- 6 
c1 <- 4
n <- r1*c1   
T <- 4
N <- n*T 
m <- n*(T-1)
###### initialize parameters ######
beta <- 0
lambda <- 0.2 
lambda0 <- 0.0
rho <- 0.5 
sig <- 1 

###### set weight matrices ######

W <- cell2nb(r1,c1,type="queen")
#rook type
Wl <- nb2listw(W)      #listw object 
Wn <- listw2dgCMatrix(Wl)  #sparse matrix
Wn.a <- as.matrix(Wn)
Mn <- Wn.a

#Wn <-  readMat("matrices.mat")
#Wn.a <- as.matrix(Wn$W)
#Mn <- Wn.a


K.size <- 100


###### negative log-likelihood function ######

log.lik.neg <-function(theta) { 
  #sig <- theta[3]
  Sn <- diag(n)-theta[1]*Wn.a 
  Rn <- diag(n)-theta[2]*Mn
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt) 
  
  ell_nt<- n*(T-1)/2*log(2*pi*sig^2)-(T-1)*(log(det(Sn))+log(det(Rn)) ) +0.5/(sig^2)* sum(diag((t(V.tilde.nt) %*% (V.tilde.nt)))) 
  return(ell_nt) 
  
} 


  
######### generate Y, X for MC.size times ##########
set.seed(1234)
MC.size <- 5
Y <- array(dim = c(n,T,MC.size))  # transformed Ynt for MC.size times
# X <- array(dim = c(n,T,MC.size))  # only suitable for dim(beta) = 1, transformed Xnt for MC.size times
# V <- array(dim = c(n,T,MC.size))  # transformed Vnt for MC.size times
Cn0 <- matrix(nrow = n, ncol = MC.size)
Vnt <-matrix(nrow = N, ncol = MC.size)
for (i in 1:MC.size) {
  
  Cn0[,i] <- rnorm(n)   # fixed effects 
#   Xnt <- runif(N,-1,1)     # non stochastic time varying regressors
#   X[,,i] <- matrix(Xnt,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Xnt,nrow=n,ncol=T)),T),n,T) 
  Vnt[,i] <- rnorm(N, sd= sig)     # error term  
  # Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn)))%*%(rep(Cn0[,i],T) + solve((diag(N)-kronecker(diag(T),rho*Mn)))%*%Vnt[,i]) 
  # Y[,,i] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
  # 
}

   
###### c.g.f. K ########
K.psi.neg <- function(nu, lambda, theta2){  ## consider a simple case theta2 = rho
  Sn <- diag(n)-lambda*Wn.a
  Rn <- diag(n)-theta2*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  ln.T <- numeric(n)
  V <- array(dim = c(n,T,MC.size))
  Y <- array(dim = c(n,T,MC.size))
  
  for (i in 1:n) {
    exp.T <- numeric(MC.size)
    # Y <- array(dim = c(n,T,MC.size))
    # V <- array(dim = c(n,T,MC.size))
    for (j in 1:MC.size) {
      psi.t <- matrix(nrow = 2, ncol = T)
      Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0[,j],T) + solve((diag(N)-kronecker(diag(T),theta2*Mn)))%*%Vnt[,j])
      Y[,,j] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
      V[,,j] <- Rn%*%(Sn %*% Y[,,j])
      for (t in 1:T) {
        #psi.t[,t] <- c(1/(sig^2)*((Rn%*%Gn%*%solve(Rn)%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*((Hn%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Hn[i,i])
       
        psi.t[,t] <-  c(1/(sig^2)*Rn[i,]%*%(Gn%*%solve(Rn)%*%V[,t,j])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*Hn[i,]%*%V[,t,j]*V[i,t,j]-(T-1)/T*Hn[i,i])
      }
      exp.T[j] <- exp(t(nu)%*%rowSums(psi.t))
    }
    ln.T[i] <- log(sum(exp.T)/MC.size)
  }
  K.psi <- sum(ln.T)/n    #### didn't add minus sign
  return(K.psi)
}


K.psi.neg1 <- function(nu, lambda, theta2){  ## consider a simple case theta2 = rho
  Sn <- diag(n)-lambda*Wn.a
  Rn <- diag(n)-theta2*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  ln.T <- numeric(n)
  V <- array(dim = c(n,T,MC.size))
  Y <- array(dim = c(n,T,MC.size))
  VnTMC <- NULL
  for (j in 1:MC.size) {
    Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0[,j],T) + solve((diag(N)-kronecker(diag(T),theta2*Mn)))%*%Vnt[,j])
    Y[,,j] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
    V[,,j] <- Rn%*%(Sn %*% Y[,,j])
    VnTMC <- rbind(VnTMC,V[,,j]) 
  }
  
  for (i in 1:n) {
    
    # Y <- array(dim = c(n,T,MC.size))
    # V <- array(dim = c(n,T,MC.size))
   
     
      
        #psi.t[,t] <- c(1/(sig^2)*((Rn%*%Gn%*%solve(Rn)%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*((Hn%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Hn[i,i])
    A <- Rn[i,]%*%(Gn%*%solve(Rn))
    
    psi.T1 <-  1/(sig^2)*diag(kronecker(diag(MC.size),A)%*%VnTMC%*%V[i,,])-rep((T-1)*Gn[i,i],MC.size)
    psi.T2 <- 1/(sig^2)*diag(kronecker(diag(MC.size),t(Hn[i,]))%*%VnTMC%*%V[i,,])-rep((T-1)*Hn[i,i],MC.size)
    psi.T <- rbind(psi.T1,psi.T2)
    exp.T <- exp(t(nu)%*%psi.T)
    
    ln.T[i] <- log(sum(exp.T)/MC.size)
  }
  K.psi <- sum(ln.T)/n    #### didn't add minus sign
  return(K.psi)
}

###### inf sup K #######

K.sup <- function(lambda, theta2){
  nlm <- nlm(f= K.psi.neg1, p = c(0,0), lambda = lambda, theta2 = theta2, ndigit = 2, gradtol = 1e-1,stepmax = 0.1285,steptol=1e-1)
  # nlm <- optim(par = c(0,0),fn = K.psi.neg,  lambda = lambda, theta2 = theta2, method = "BFGS")
  #nlm <- nlminb(start = c(0,0), objective = K.psi.neg,lambda = lambda, theta2 = theta2,control = list(rel.tol=1e-6, step.min=1e-4,step.max= 0.5))
  k.sup <- -nlm$minimum
  nu.est <- nlm$estimate
  # k.sup <- -nlm$objective
  # nu.est <- nlm$par
  print(nu.est)
  return(k.sup)  # add minus sign
}


K.sup.grad <- function(lambda, theta2){
  nu.start <- c(0,0)
  K.g <- grad(K.psi.neg, x= nu.start, lambda = lambda, theta2 = theta2)   
  K.j <- hessian(K.psi.neg,x=nu.start, lambda = lambda, theta2 = theta2)
  nu <- nu.start-solve(K.j) %*% K.g
  while (norm(nu-nu.start) > 0.01) {
    nu.start <- nu
    K.g <- grad(K.psi.neg, x= nu.start, lambda = lambda, theta2 = theta2)   
    K.j <- hessian(K.psi.neg,x=nu.start, lambda = lambda, theta2 = theta2)
    nu <- nu.start-solve(K.j) %*% K.g
  }
  return(-K.psi.neg(nu,lambda,theta2))
  
}



K.inf <- function(lambda){
  #return(optimize(K.sup,interval=c(-0.99,0.99), lambda = lambda)$objective)
  #return(optim(fn= K.sup, par = 0, lambda = lambda, lower = -0.99, upper = 0.99, method = "L-BFGS-B")$value)
  f <- nlminb(start = 0, objective = K.sup, lambda = lambda, lower = -0.99, upper = 0.99,control = list(rel.tol=1e-2, step.min=0.01,step.max= 0.13,eval.max=10))
  k.inf <- f$objective
  theta2.est <- f$par
  # seq.theta2 <- seq(-0.9,0.9,0.1)
  # k.start <- K.sup(lambda, seq.theta2[1])
  # for (i in 2:length(seq.theta2)) {
  #   k.inf1 <- K.sup(lambda, seq.theta2[i])
  #   if(k.inf1 < k.start){
  #     k.start <- k.inf1
  #     theta21.est <- seq.theta2[i]
  #   }
  # 
  # 
  # }

  
  return(c(k.inf,theta2.est))
}

# seq.theta2 <- seq(-0.9,0.9,0.1)
# k.start <- K.sup(0.3, seq.theta2[1])
# for (i in 2:length(seq.theta2)) {
#   k.inf1 <- K.sup(0.3, seq.theta2[i])
#   if(k.inf1 < k.start){
#     k.start <- k.inf1
#     theta21.est <- seq.theta2[i]
#   }
# 
# 
# }

SIGMA <- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  G.dot <- Rn%*%Gn%*%Rn.inv
  X.dot <- Rn%*%X.tilde.nt
  
#  a <- sum(diag(t(X.dot)%*%X.dot))/(m*sig2)
#  b <- beta*(sum(diag(t(X.dot)%*%G.dot%*%X.dot)))/(m*sig2)
  c<- 1/n*(sum(diag(t(G.dot%*%X.dot)%*%(G.dot%*%X.dot)))*beta*beta/((T-1)*sig2)+sum(diag((G.dot+t(G.dot))%*%G.dot)))
  d <- sum(diag((t(Hn)+Hn)%*%G.dot))/n
#  e <- sum(diag(G.dot))/(n*sig2)
  f <- sum(diag((t(Hn)+Hn)%*%Hn))/n
#  g <- sum(diag(Hn))/(n*sig2)
#  h <- 1/(2*sig2*sig2)
 # Sig <-  matrix(c(a,b,0,0,b,c,d,e,0,d,f,g,0,e,g,h),4,4)
  Sig <-  matrix(c(c,d,d,f),2,2)
  return(Sig)
} 

ptm <- proc.time()


set.seed(1234)


lambda.hat <- numeric(K.size)
rho.hat <- numeric(K.size)
theta2.est <- numeric(K.size)
test.stat <- numeric(K.size)
test.wald <- numeric(K.size)
for (k in 1:K.size) {
  Cn0.h <- runif(n,-1,1)
   Xnt <- runif(N,-1,1)
  # Xnt <- rnorm(N)
  Ent.h <- rnorm(N, sd= sig)
  Yn <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0.h,T) + solve((diag(N)-kronecker(diag(T),rho*Mn)))%*%Ent.h) 
  Y.tilde.nt<-matrix(Yn,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn,nrow=n,ncol=T)),T),n,T)
  X.tilde.nt<-matrix(Xnt,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Xnt,nrow=n,ncol=T)),T),n,T)
  # dt <- data.frame(id = rep(c(1:n),T), time =kronecker(1:T,rep(1,n)),Yn,Xnt)
  # sarar <- spml(formula = Yn~Xnt, data = dt, listw = Wl, model = "within", spatial.error= "b",lag = T, LeeYu = T, Hess = F)
   # lambda.hat[k] <- sarar$coefficients[1]
   # rho.hat[k] <- sarar$coefficients[2]
  # test.wald[k] <- (lambda.hat-lambda0)^2/sarar$vcov[1,1]
  
 
  
  para <- nlminb(start = c(0,0.1), objective = log.lik.neg, lower = c(-0.99,-0.99), upper = c(0.99,0.99),control = list(rel.tol=1e-2,step.max= 0.02,eval.max=10))$par
  lambda.hat[k] <- para[1]
  rho.hat[k] <- para[2]

  lamb.cov <- solve(SIGMA(beta,lambda.hat[k],rho.hat[k],sig^2))[2,2]/m
  test.wald[k] <- (lambda.hat[k]-lambda0)^2/lamb.cov

  # nll <- log.lik.neg(c(-0.99,-0.99))
  # seq.a <- seq(-0.99,0.99,0.01)a
  # for (i in 1:length(seq.a)) {
  #   for (j in 1:length(seq.a)) {
  #     if(nll>log.lik.neg(c(seq.a[i],seq.a[j]))){nll <- log.lik.neg(c(seq.a[i],seq.a[j]))
  #     lambda.hat1 <- seq.a[i]
  #     rho.hat <- seq.a[j]}
  # 
  #   }
  # 
  # }
  #  print(lambda.hat1)

  k.inf <- K.inf(lambda.hat[k])
  h.hat <- k.inf[1]
  theta2.est[k] <- k.inf[-1]
  test.stat[k] <- 2*n*h.hat
  print(test.stat[k])
  print(k)
}
# seq.lam <- seq(-0.3,0.3,0.1)
# T.S <- numeric(length = length(seq.lam))
# for (i in 1:length(seq.lam)) {
#   T.S[i] <- K.inf(seq.lam[i])[1]
#   print(i)
#   
# }
# T.S*2*n
proc.time() - ptm

# quantile(test.stat,0.9)
# quantile(test.wald,0.9)
# qchisq(0.9,df=1)

t_s_09505<-quantile(test.stat,0.95)
t_wald09505<-quantile(test.wald,0.95)
qchisq(0.95,df=1)

pchisq(t_s_09505,1)

t_s_097505<-quantile(test.stat,0.975)
t_wald097505<-quantile(test.wald,0.975)
qchisq(0.975,df=1)

pchisq(t_s_097505,1)

#######################################################
#######################################################
#######################################################

rho <- 0.25 


r1 <- 6 
c1 <- 4
n <- r1*c1   
T <- 4
N <- n*T 
m <- n*(T-1)
###### initialize parameters ######
beta <- 0
lambda <- 0.2 
lambda0 <- 0.0
rho <- 0.25 
sig <- 1 

###### set weight matrices ######

W <- cell2nb(r1,c1,type="queen")
#rook type
Wl <- nb2listw(W)      #listw object 
Wn <- listw2dgCMatrix(Wl)  #sparse matrix
Wn.a <- as.matrix(Wn)
Mn <- Wn.a

#Wn <-  readMat("matrices.mat")
#Wn.a <- as.matrix(Wn$W)
#Mn <- Wn.a


###### negative log-likelihood function ######

log.lik.neg <-function(theta) { 
  #sig <- theta[3]
  Sn <- diag(n)-theta[1]*Wn.a 
  Rn <- diag(n)-theta[2]*Mn
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt) 
  
  ell_nt<- n*(T-1)/2*log(2*pi*sig^2)-(T-1)*(log(det(Sn))+log(det(Rn)) ) +0.5/(sig^2)* sum(diag((t(V.tilde.nt) %*% (V.tilde.nt)))) 
  return(ell_nt) 
  
} 



######### generate Y, X for MC.size times ##########
set.seed(1234)
MC.size <- 100
Y <- array(dim = c(n,T,MC.size))  # transformed Ynt for MC.size times
# X <- array(dim = c(n,T,MC.size))  # only suitable for dim(beta) = 1, transformed Xnt for MC.size times
# V <- array(dim = c(n,T,MC.size))  # transformed Vnt for MC.size times
Cn0 <- matrix(nrow = n, ncol = MC.size)
Vnt <-matrix(nrow = N, ncol = MC.size)
for (i in 1:MC.size) {
  
  Cn0[,i] <- rnorm(n)   # fixed effects 
  #   Xnt <- runif(N,-1,1)     # non stochastic time varying regressors
  #   X[,,i] <- matrix(Xnt,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Xnt,nrow=n,ncol=T)),T),n,T) 
  Vnt[,i] <- rnorm(N, sd= sig)     # error term  
  # Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn)))%*%(rep(Cn0[,i],T) + solve((diag(N)-kronecker(diag(T),rho*Mn)))%*%Vnt[,i]) 
  # Y[,,i] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
  # 
}


###### c.g.f. K ########
K.psi.neg <- function(nu, lambda, theta2){  ## consider a simple case theta2 = rho
  Sn <- diag(n)-lambda*Wn.a
  Rn <- diag(n)-theta2*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  ln.T <- numeric(n)
  V <- array(dim = c(n,T,MC.size))
  Y <- array(dim = c(n,T,MC.size))
  
  for (i in 1:n) {
    exp.T <- numeric(MC.size)
    # Y <- array(dim = c(n,T,MC.size))
    # V <- array(dim = c(n,T,MC.size))
    for (j in 1:MC.size) {
      psi.t <- matrix(nrow = 2, ncol = T)
      Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0[,j],T) + solve((diag(N)-kronecker(diag(T),theta2*Mn)))%*%Vnt[,j])
      Y[,,j] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
      V[,,j] <- Rn%*%(Sn %*% Y[,,j])
      for (t in 1:T) {
        #psi.t[,t] <- c(1/(sig^2)*((Rn%*%Gn%*%solve(Rn)%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*((Hn%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Hn[i,i])
        
        psi.t[,t] <-  c(1/(sig^2)*Rn[i,]%*%(Gn%*%solve(Rn)%*%V[,t,j])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*Hn[i,]%*%V[,t,j]*V[i,t,j]-(T-1)/T*Hn[i,i])
      }
      exp.T[j] <- exp(t(nu)%*%rowSums(psi.t))
    }
    ln.T[i] <- log(sum(exp.T)/MC.size)
  }
  K.psi <- sum(ln.T)/n    #### didn't add minus sign
  return(K.psi)
}


K.psi.neg1 <- function(nu, lambda, theta2){  ## consider a simple case theta2 = rho
  Sn <- diag(n)-lambda*Wn.a
  Rn <- diag(n)-theta2*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  ln.T <- numeric(n)
  V <- array(dim = c(n,T,MC.size))
  Y <- array(dim = c(n,T,MC.size))
  VnTMC <- NULL
  for (j in 1:MC.size) {
    Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0[,j],T) + solve((diag(N)-kronecker(diag(T),theta2*Mn)))%*%Vnt[,j])
    Y[,,j] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
    V[,,j] <- Rn%*%(Sn %*% Y[,,j])
    VnTMC <- rbind(VnTMC,V[,,j]) 
  }
  
  for (i in 1:n) {
    
    # Y <- array(dim = c(n,T,MC.size))
    # V <- array(dim = c(n,T,MC.size))
    
    
    
    #psi.t[,t] <- c(1/(sig^2)*((Rn%*%Gn%*%solve(Rn)%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*((Hn%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Hn[i,i])
    A <- Rn[i,]%*%(Gn%*%solve(Rn))
    
    psi.T1 <-  1/(sig^2)*diag(kronecker(diag(MC.size),A)%*%VnTMC%*%V[i,,])-rep((T-1)*Gn[i,i],MC.size)
    psi.T2 <- 1/(sig^2)*diag(kronecker(diag(MC.size),t(Hn[i,]))%*%VnTMC%*%V[i,,])-rep((T-1)*Hn[i,i],MC.size)
    psi.T <- rbind(psi.T1,psi.T2)
    exp.T <- exp(t(nu)%*%psi.T)
    
    ln.T[i] <- log(sum(exp.T)/MC.size)
  }
  K.psi <- sum(ln.T)/n    #### didn't add minus sign
  return(K.psi)
}

###### inf sup K #######

K.sup <- function(lambda, theta2){
  nlm <- nlm(f= K.psi.neg1, p = c(0,0), lambda = lambda, theta2 = theta2, ndigit = 2, gradtol = 1e-1,stepmax = 0.18,steptol=1e-1)
  # nlm <- optim(par = c(0,0),fn = K.psi.neg,  lambda = lambda, theta2 = theta2, method = "BFGS")
  #nlm <- nlminb(start = c(0,0), objective = K.psi.neg,lambda = lambda, theta2 = theta2,control = list(rel.tol=1e-6, step.min=1e-4,step.max= 0.5))
  k.sup <- -nlm$minimum
  nu.est <- nlm$estimate
  # k.sup <- -nlm$objective
  # nu.est <- nlm$par
  print(nu.est)
  return(k.sup)  # add minus sign
}


K.sup.grad <- function(lambda, theta2){
  nu.start <- c(0,0)
  K.g <- grad(K.psi.neg, x= nu.start, lambda = lambda, theta2 = theta2)   
  K.j <- hessian(K.psi.neg,x=nu.start, lambda = lambda, theta2 = theta2)
  nu <- nu.start-solve(K.j) %*% K.g
  while (norm(nu-nu.start) > 0.01) {
    nu.start <- nu
    K.g <- grad(K.psi.neg, x= nu.start, lambda = lambda, theta2 = theta2)   
    K.j <- hessian(K.psi.neg,x=nu.start, lambda = lambda, theta2 = theta2)
    nu <- nu.start-solve(K.j) %*% K.g
  }
  return(-K.psi.neg(nu,lambda,theta2))
  
}



K.inf <- function(lambda){
  #return(optimize(K.sup,interval=c(-0.99,0.99), lambda = lambda)$objective)
  #return(optim(fn= K.sup, par = 0, lambda = lambda, lower = -0.99, upper = 0.99, method = "L-BFGS-B")$value)
  f <- nlminb(start = 0, objective = K.sup, lambda = lambda, lower = -0.99, upper = 0.99,control = list(rel.tol=1e-2, step.min=0.01,step.max= 0.04,eval.max=10))
  k.inf <- f$objective
  theta2.est <- f$par
  # seq.theta2 <- seq(-0.9,0.9,0.1)
  # k.start <- K.sup(lambda, seq.theta2[1])
  # for (i in 2:length(seq.theta2)) {
  #   k.inf1 <- K.sup(lambda, seq.theta2[i])
  #   if(k.inf1 < k.start){
  #     k.start <- k.inf1
  #     theta21.est <- seq.theta2[i]
  #   }
  # 
  # 
  # }
  
  
  return(c(k.inf,theta2.est))
}

# seq.theta2 <- seq(-0.9,0.9,0.1)
# k.start <- K.sup(0.3, seq.theta2[1])
# for (i in 2:length(seq.theta2)) {
#   k.inf1 <- K.sup(0.3, seq.theta2[i])
#   if(k.inf1 < k.start){
#     k.start <- k.inf1
#     theta21.est <- seq.theta2[i]
#   }
# 
# 
# }

SIGMA <- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  G.dot <- Rn%*%Gn%*%Rn.inv
  X.dot <- Rn%*%X.tilde.nt
  
  #  a <- sum(diag(t(X.dot)%*%X.dot))/(m*sig2)
  #  b <- beta*(sum(diag(t(X.dot)%*%G.dot%*%X.dot)))/(m*sig2)
  c<- 1/n*(sum(diag(t(G.dot%*%X.dot)%*%(G.dot%*%X.dot)))*beta*beta/((T-1)*sig2)+sum(diag((G.dot+t(G.dot))%*%G.dot)))
  d <- sum(diag((t(Hn)+Hn)%*%G.dot))/n
  #  e <- sum(diag(G.dot))/(n*sig2)
  f <- sum(diag((t(Hn)+Hn)%*%Hn))/n
  #  g <- sum(diag(Hn))/(n*sig2)
  #  h <- 1/(2*sig2*sig2)
  # Sig <-  matrix(c(a,b,0,0,b,c,d,e,0,d,f,g,0,e,g,h),4,4)
  Sig <-  matrix(c(c,d,d,f),2,2)
  return(Sig)
} 

ptm <- proc.time()


set.seed(1234)


lambda.hat <- numeric(K.size)
rho.hat <- numeric(K.size)
theta2.est <- numeric(K.size)
test.stat <- numeric(K.size)
test.wald <- numeric(K.size)
for (k in 1:K.size) {
  Cn0.h <- runif(n,-1,1)
  Xnt <- runif(N,-1,1)
  # Xnt <- rnorm(N)
  Ent.h <- rnorm(N, sd= sig)
  Yn <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0.h,T) + solve((diag(N)-kronecker(diag(T),rho*Mn)))%*%Ent.h) 
  Y.tilde.nt<-matrix(Yn,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn,nrow=n,ncol=T)),T),n,T)
  X.tilde.nt<-matrix(Xnt,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Xnt,nrow=n,ncol=T)),T),n,T)
  # dt <- data.frame(id = rep(c(1:n),T), time =kronecker(1:T,rep(1,n)),Yn,Xnt)
  # sarar <- spml(formula = Yn~Xnt, data = dt, listw = Wl, model = "within", spatial.error= "b",lag = T, LeeYu = T, Hess = F)
  # lambda.hat[k] <- sarar$coefficients[1]
  # rho.hat[k] <- sarar$coefficients[2]
  # test.wald[k] <- (lambda.hat-lambda0)^2/sarar$vcov[1,1]
  
  
  
  para <- nlminb(start = c(0,0.1), objective = log.lik.neg, lower = c(-0.99,-0.99), upper = c(0.99,0.99),control = list(rel.tol=1e-2,step.max= 0.04,eval.max=10))$par
  lambda.hat[k] <- para[1]
  rho.hat[k] <- para[2]
  
  lamb.cov <- solve(SIGMA(beta,lambda.hat[k],rho.hat[k],sig^2))[2,2]/m
  test.wald[k] <- (lambda.hat[k]-lambda0)^2/lamb.cov
  
  # nll <- log.lik.neg(c(-0.99,-0.99))
  # seq.a <- seq(-0.99,0.99,0.01)a
  # for (i in 1:length(seq.a)) {
  #   for (j in 1:length(seq.a)) {
  #     if(nll>log.lik.neg(c(seq.a[i],seq.a[j]))){nll <- log.lik.neg(c(seq.a[i],seq.a[j]))
  #     lambda.hat1 <- seq.a[i]
  #     rho.hat <- seq.a[j]}
  # 
  #   }
  # 
  # }
  #  print(lambda.hat1)
  
  k.inf <- K.inf(lambda.hat[k])
  h.hat <- k.inf[1]
  theta2.est[k] <- k.inf[-1]
  test.stat[k] <- 2*n*h.hat
  print(test.stat[k])
  print(k)
}
# seq.lam <- seq(-0.3,0.3,0.1)
# T.S <- numeric(length = length(seq.lam))
# for (i in 1:length(seq.lam)) {
#   T.S[i] <- K.inf(seq.lam[i])[1]
#   print(i)
#   
# }
# T.S*2*n
proc.time() - ptm

# quantile(test.stat,0.9)
# quantile(test.wald,0.9)
# qchisq(0.9,df=1)


t_s_095025<-quantile(test.stat,0.95)
t_wald095025<-quantile(test.wald,0.95)
qchisq(0.95,df=1)

t_s_0975025<-quantile(test.stat,0.975)
t_wald0975025<-quantile(test.wald,0.975)
qchisq(0.975,df=1)


##################################################
##################################################
##################################################

rho<-0.75

r1 <- 6 
c1 <- 4
n <- r1*c1   
T <- 4
N <- n*T 
m <- n*(T-1)
###### initialize parameters ######
beta <- 0
lambda <- 0.2 
lambda0 <- 0.0
rho <- 0.75 
sig <- 1 

###### set weight matrices ######

W <- cell2nb(r1,c1,type="queen")
#rook type
Wl <- nb2listw(W)      #listw object 
Wn <- listw2dgCMatrix(Wl)  #sparse matrix
Wn.a <- as.matrix(Wn)
Mn <- Wn.a

#Wn <-  readMat("matrices.mat")
#Wn.a <- as.matrix(Wn$W)
#Mn <- Wn.a


###### negative log-likelihood function ######

log.lik.neg <-function(theta) { 
  #sig <- theta[3]
  Sn <- diag(n)-theta[1]*Wn.a 
  Rn <- diag(n)-theta[2]*Mn
  V.tilde.nt<-Rn%*%(Sn %*% Y.tilde.nt) 
  
  ell_nt<- n*(T-1)/2*log(2*pi*sig^2)-(T-1)*(log(det(Sn))+log(det(Rn)) ) +0.5/(sig^2)* sum(diag((t(V.tilde.nt) %*% (V.tilde.nt)))) 
  return(ell_nt) 
  
} 



######### generate Y, X for MC.size times ##########
set.seed(1234)
MC.size <- 100
Y <- array(dim = c(n,T,MC.size))  # transformed Ynt for MC.size times
# X <- array(dim = c(n,T,MC.size))  # only suitable for dim(beta) = 1, transformed Xnt for MC.size times
# V <- array(dim = c(n,T,MC.size))  # transformed Vnt for MC.size times
Cn0 <- matrix(nrow = n, ncol = MC.size)
Vnt <-matrix(nrow = N, ncol = MC.size)
for (i in 1:MC.size) {
  
  Cn0[,i] <- rnorm(n)   # fixed effects 
  #   Xnt <- runif(N,-1,1)     # non stochastic time varying regressors
  #   X[,,i] <- matrix(Xnt,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Xnt,nrow=n,ncol=T)),T),n,T) 
  Vnt[,i] <- rnorm(N, sd= sig)     # error term  
  # Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn)))%*%(rep(Cn0[,i],T) + solve((diag(N)-kronecker(diag(T),rho*Mn)))%*%Vnt[,i]) 
  # Y[,,i] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
  # 
}


###### c.g.f. K ########
K.psi.neg <- function(nu, lambda, theta2){  ## consider a simple case theta2 = rho
  Sn <- diag(n)-lambda*Wn.a
  Rn <- diag(n)-theta2*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  ln.T <- numeric(n)
  V <- array(dim = c(n,T,MC.size))
  Y <- array(dim = c(n,T,MC.size))
  
  for (i in 1:n) {
    exp.T <- numeric(MC.size)
    # Y <- array(dim = c(n,T,MC.size))
    # V <- array(dim = c(n,T,MC.size))
    for (j in 1:MC.size) {
      psi.t <- matrix(nrow = 2, ncol = T)
      Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0[,j],T) + solve((diag(N)-kronecker(diag(T),theta2*Mn)))%*%Vnt[,j])
      Y[,,j] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
      V[,,j] <- Rn%*%(Sn %*% Y[,,j])
      for (t in 1:T) {
        #psi.t[,t] <- c(1/(sig^2)*((Rn%*%Gn%*%solve(Rn)%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*((Hn%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Hn[i,i])
        
        psi.t[,t] <-  c(1/(sig^2)*Rn[i,]%*%(Gn%*%solve(Rn)%*%V[,t,j])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*Hn[i,]%*%V[,t,j]*V[i,t,j]-(T-1)/T*Hn[i,i])
      }
      exp.T[j] <- exp(t(nu)%*%rowSums(psi.t))
    }
    ln.T[i] <- log(sum(exp.T)/MC.size)
  }
  K.psi <- sum(ln.T)/n    #### didn't add minus sign
  return(K.psi)
}


K.psi.neg1 <- function(nu, lambda, theta2){  ## consider a simple case theta2 = rho
  Sn <- diag(n)-lambda*Wn.a
  Rn <- diag(n)-theta2*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  ln.T <- numeric(n)
  V <- array(dim = c(n,T,MC.size))
  Y <- array(dim = c(n,T,MC.size))
  VnTMC <- NULL
  for (j in 1:MC.size) {
    Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0[,j],T) + solve((diag(N)-kronecker(diag(T),theta2*Mn)))%*%Vnt[,j])
    Y[,,j] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
    V[,,j] <- Rn%*%(Sn %*% Y[,,j])
    VnTMC <- rbind(VnTMC,V[,,j]) 
  }
  
  for (i in 1:n) {
    
    # Y <- array(dim = c(n,T,MC.size))
    # V <- array(dim = c(n,T,MC.size))
    
    
    
    #psi.t[,t] <- c(1/(sig^2)*((Rn%*%Gn%*%solve(Rn)%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*((Hn%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Hn[i,i])
    A <- Rn[i,]%*%(Gn%*%solve(Rn))
    
    psi.T1 <-  1/(sig^2)*diag(kronecker(diag(MC.size),A)%*%VnTMC%*%V[i,,])-rep((T-1)*Gn[i,i],MC.size)
    psi.T2 <- 1/(sig^2)*diag(kronecker(diag(MC.size),t(Hn[i,]))%*%VnTMC%*%V[i,,])-rep((T-1)*Hn[i,i],MC.size)
    psi.T <- rbind(psi.T1,psi.T2)
    exp.T <- exp(t(nu)%*%psi.T)
    
    ln.T[i] <- log(sum(exp.T)/MC.size)
  }
  K.psi <- sum(ln.T)/n    #### didn't add minus sign
  return(K.psi)
}

###### inf sup K #######

K.sup <- function(lambda, theta2){
  nlm <- nlm(f= K.psi.neg1, p = c(0,0), lambda = lambda, theta2 = theta2, ndigit = 2, gradtol = 1e-1,stepmax = 0.06,steptol=1e-1)
  # nlm <- optim(par = c(0,0),fn = K.psi.neg,  lambda = lambda, theta2 = theta2, method = "BFGS")
  #nlm <- nlminb(start = c(0,0), objective = K.psi.neg,lambda = lambda, theta2 = theta2,control = list(rel.tol=1e-6, step.min=1e-4,step.max= 0.5))
  k.sup <- -nlm$minimum
  nu.est <- nlm$estimate
  # k.sup <- -nlm$objective
  # nu.est <- nlm$par
  print(nu.est)
  return(k.sup)  # add minus sign
}


K.sup.grad <- function(lambda, theta2){
  nu.start <- c(0,0)
  K.g <- grad(K.psi.neg, x= nu.start, lambda = lambda, theta2 = theta2)   
  K.j <- hessian(K.psi.neg,x=nu.start, lambda = lambda, theta2 = theta2)
  nu <- nu.start-solve(K.j) %*% K.g
  while (norm(nu-nu.start) > 0.01) {
    nu.start <- nu
    K.g <- grad(K.psi.neg, x= nu.start, lambda = lambda, theta2 = theta2)   
    K.j <- hessian(K.psi.neg,x=nu.start, lambda = lambda, theta2 = theta2)
    nu <- nu.start-solve(K.j) %*% K.g
  }
  return(-K.psi.neg(nu,lambda,theta2))
  
}



K.inf <- function(lambda){
  #return(optimize(K.sup,interval=c(-0.99,0.99), lambda = lambda)$objective)
  #return(optim(fn= K.sup, par = 0, lambda = lambda, lower = -0.99, upper = 0.99, method = "L-BFGS-B")$value)
  f <- nlminb(start = 0, objective = K.sup, lambda = lambda, lower = -0.99, upper = 0.99,control = list(rel.tol=1e-2, step.min=0.01,step.max= 0.04,eval.max=10))
  k.inf <- f$objective
  theta2.est <- f$par
  # seq.theta2 <- seq(-0.9,0.9,0.1)
  # k.start <- K.sup(lambda, seq.theta2[1])
  # for (i in 2:length(seq.theta2)) {
  #   k.inf1 <- K.sup(lambda, seq.theta2[i])
  #   if(k.inf1 < k.start){
  #     k.start <- k.inf1
  #     theta21.est <- seq.theta2[i]
  #   }
  # 
  # 
  # }
  
  
  return(c(k.inf,theta2.est))
}

# seq.theta2 <- seq(-0.9,0.9,0.1)
# k.start <- K.sup(0.3, seq.theta2[1])
# for (i in 2:length(seq.theta2)) {
#   k.inf1 <- K.sup(0.3, seq.theta2[i])
#   if(k.inf1 < k.start){
#     k.start <- k.inf1
#     theta21.est <- seq.theta2[i]
#   }
# 
# 
# }

SIGMA <- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  Sn.inv <- solve(Sn)
  Rn.inv <- solve(Rn)
  G.dot <- Rn%*%Gn%*%Rn.inv
  X.dot <- Rn%*%X.tilde.nt
  
  #  a <- sum(diag(t(X.dot)%*%X.dot))/(m*sig2)
  #  b <- beta*(sum(diag(t(X.dot)%*%G.dot%*%X.dot)))/(m*sig2)
  c<- 1/n*(sum(diag(t(G.dot%*%X.dot)%*%(G.dot%*%X.dot)))*beta*beta/((T-1)*sig2)+sum(diag((G.dot+t(G.dot))%*%G.dot)))
  d <- sum(diag((t(Hn)+Hn)%*%G.dot))/n
  #  e <- sum(diag(G.dot))/(n*sig2)
  f <- sum(diag((t(Hn)+Hn)%*%Hn))/n
  #  g <- sum(diag(Hn))/(n*sig2)
  #  h <- 1/(2*sig2*sig2)
  # Sig <-  matrix(c(a,b,0,0,b,c,d,e,0,d,f,g,0,e,g,h),4,4)
  Sig <-  matrix(c(c,d,d,f),2,2)
  return(Sig)
} 

ptm <- proc.time()


set.seed(1234)


lambda.hat <- numeric(K.size)
rho.hat <- numeric(K.size)
theta2.est <- numeric(K.size)
test.stat <- numeric(K.size)
test.wald <- numeric(K.size)
for (k in 1:K.size) {
  Cn0.h <- runif(n,-1,1)
  Xnt <- runif(N,-1,1)
  # Xnt <- rnorm(N)
  Ent.h <- rnorm(N, sd= sig)
  Yn <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0.h,T) + solve((diag(N)-kronecker(diag(T),rho*Mn)))%*%Ent.h) 
  Y.tilde.nt<-matrix(Yn,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn,nrow=n,ncol=T)),T),n,T)
  X.tilde.nt<-matrix(Xnt,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Xnt,nrow=n,ncol=T)),T),n,T)
  # dt <- data.frame(id = rep(c(1:n),T), time =kronecker(1:T,rep(1,n)),Yn,Xnt)
  # sarar <- spml(formula = Yn~Xnt, data = dt, listw = Wl, model = "within", spatial.error= "b",lag = T, LeeYu = T, Hess = F)
  # lambda.hat[k] <- sarar$coefficients[1]
  # rho.hat[k] <- sarar$coefficients[2]
  # test.wald[k] <- (lambda.hat-lambda0)^2/sarar$vcov[1,1]
  
  
  
  para <- nlminb(start = c(0,0.1), objective = log.lik.neg, lower = c(-0.99,-0.99), upper = c(0.99,0.99),control = list(rel.tol=1e-2,step.max= 0.04,eval.max=10))$par
  lambda.hat[k] <- para[1]
  rho.hat[k] <- para[2]
  
  lamb.cov <- solve(SIGMA(beta,lambda.hat[k],rho.hat[k],sig^2))[2,2]/m
  test.wald[k] <- (lambda.hat[k]-lambda0)^2/lamb.cov
  
  # nll <- log.lik.neg(c(-0.99,-0.99))
  # seq.a <- seq(-0.99,0.99,0.01)a
  # for (i in 1:length(seq.a)) {
  #   for (j in 1:length(seq.a)) {
  #     if(nll>log.lik.neg(c(seq.a[i],seq.a[j]))){nll <- log.lik.neg(c(seq.a[i],seq.a[j]))
  #     lambda.hat1 <- seq.a[i]
  #     rho.hat <- seq.a[j]}
  # 
  #   }
  # 
  # }
  #  print(lambda.hat1)
  
  k.inf <- K.inf(lambda.hat[k])
  h.hat <- k.inf[1]
  theta2.est[k] <- k.inf[-1]
  test.stat[k] <- 2*n*h.hat
  print(test.stat[k])
  print(k)
}
# seq.lam <- seq(-0.3,0.3,0.1)
# T.S <- numeric(length = length(seq.lam))
# for (i in 1:length(seq.lam)) {
#   T.S[i] <- K.inf(seq.lam[i])[1]
#   print(i)
#   
# }
# T.S*2*n
proc.time() - ptm

# quantile(test.stat,0.9)
# quantile(test.wald,0.9)
# qchisq(0.9,df=1)



t_s_095075<-quantile(test.stat,0.95)
t_wald095075<-quantile(test.wald,0.95)
qchisq(0.95,df=1)

t_s_0975075<-quantile(test.stat,0.975)
t_wald0975075<-quantile(test.wald,0.975)
qchisq(0.975,df=1)


########################################

c(t_s_095025,t_s_0975025)


c(t_wald095025,t_wald0975025)

c(t_s_09505,t_s_097505)
c(t_wald09505,t_wald097505)


c(t_s_095075,t_s_0975075)
c(t_wald095075,t_wald0975075)


