#####################################################################
############# Test for the null hypothese lambda = 0 ################
#####################################################################
library("spdep") 
library("splm") 
library("graphics")
library("numDeriv")
library("R.matlab")



#############load data and weight matrices #########

FHD <- readMat("FHD.mat")
Wn <-  readMat("matrices.mat")

############# load data ###############

index <- 1  ## time index: 1 -> 1960-1970, 2 -> 1971-1985, 3 -> 1986-2000

if(index == 1){
  Yn <- FHD$data[1:264,3]  #### Investment rate
  Xn <- FHD$data[1:264,5]  #### Saving rate
  timeline <- 1960:1970
}else if(index == 2){
  Yn <- FHD$data[265:624,3]  #### Investment rate
  Xn <- FHD$data[265:624,5]  #### Saving rate
  timeline <- 1971:1985
}else if(index == 3){
  Yn <- FHD$data[625:984,3]  #### Investment rate
  Xn <- FHD$data[625:984,5]  #### Saving rate
  timeline <- 1986:2000
}else{
  print("Wrong data")
}

W.index <- 1 ## weight matrix index: 1 -> 7 nearest neighbors, 2 -> inverse distance

if(W.index == 1){
  W <- Wn$W   # 7NN
}else if(W.index==2){
  W <- Wn$Wdia1  # inverse distance 
}else{
  print("Wrong weight matrix")
}



Wn.a <- as.matrix(W)
Mn <- Wn.a
N <- length(Yn)
n <- 24
T <- N/n
m <- n*(T-1)

lambda0 <- 0
rho0 <- 0

lw <- mat2listw(Wn.a)
dt <- data.frame(id = rep(c(1:n),T), time =kronecker(timeline,rep(1,24)),Yn,Xn)

######## transform Y_{nt} and X_{nt} ########

Y.tilde.nt<-matrix(Yn,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn,nrow=n,ncol=T)),T),n,T)
X.tilde.nt<-matrix(Xn,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Xn,nrow=n,ncol=T)),T),n,T)


#########################################################################
# Table 2: SARAR(1,1) model: Maximum likelihood estimates of Parameters # 
# beta, lambda, rho. Standard errors are between brackets.              #
#########################################################################

sarar <- spml(formula = Yn~Xn, data = dt , listw = mat2listw(W), model = "within", spatial.error= "b",lag = T, LeeYu = T, Hess = F)
summary(sarar)
beta.hat <- sarar$coefficients[3]
rho.hat <- sarar$coefficients[2]
sig2.hat <- sarar$sigma2
lambda.hat <- sarar$coefficients[1]


  
######### generate Y, X for MC.size times ##########
set.seed(1234)
MC.size <- 50

Cn0 <- matrix(nrow = n, ncol = MC.size)
for (i in 1:MC.size) {
  
  Cn0[,i] <- rnorm(n)   # fixed effects 
  
}

   
###### c.g.f. K ########
K.psi.neg <- function(nu, lambda, theta2){  ## consider a simple case theta2 = rho
  beta <- theta2[1] 
  rho <- theta2[2]
  sig2 <- theta2[3]
  Sn <- diag(n)-lambda*Wn.a
  Rn <- diag(n)-rho*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  ln.T <- numeric(n)
  V <- array(dim = c(n,T,MC.size))
  Y <- array(dim = c(n,T,MC.size))
  Vnt <- matrix(rnorm(N*MC.size,sd = sqrt(sig2)),nrow = N, ncol = MC.size)
  for (i in 1:n) {
    exp.T <- numeric(MC.size)
    # Y <- array(dim = c(n,T,MC.size))
    # V <- array(dim = c(n,T,MC.size))
    for (j in 1:MC.size) {
      psi.t <- matrix(nrow = 4, ncol = T)
      Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0[,j],T) + Xn*beta + solve((diag(N)-kronecker(diag(T),rho*Mn)))%*%Vnt[,j])
      Y[,,j] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
      V[,,j] <- Rn%*%(Sn %*% Y[,,j]-X.tilde.nt*beta) 
      for (t in 1:T) {
        
        #psi.t[,t] <- c(1/(sig^2)*((Rn%*%Gn%*%solve(Rn)%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*((Hn%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Hn[i,i])
       
        psi.t[,t] <-  c(1/sig2*Rn[i,]%*%X.tilde.nt[,t]*V[i,t,j],
          1/(sig2)*Rn[i,]%*%(Gn%*%X.tilde.nt[,t]*beta+Gn%*%solve(Rn)%*%V[,t,j])*V[i,t,j]-(T-1)/T*Gn[i,i],
          1/(sig2)*Hn[i,]%*%V[,t,j]*V[i,t,j]-(T-1)/T*Hn[i,i],
          1/(2*sig2*sig2)*(V[i,t,j]*V[i,t,j]-(T-1)/T*sig2))
      }
      exp.T[j] <- exp(t(nu)%*%rowSums(psi.t))
    }
    ln.T[i] <- log(mean(exp.T))
  }
  K.psi <- sum(ln.T)/n    #### didn't add minus sign
  return(K.psi)
}


K.psi.neg1 <- function(nu, lambda, theta2){  
  beta <- theta2[1] 
  rho <- theta2[2]
  sig2 <- theta2[3]
  Sn <- diag(n)-lambda*Wn.a
  Rn <- diag(n)-rho*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
  ln.T <- numeric(n)
  V <- array(dim = c(n,T,MC.size))
  Y <- array(dim = c(n,T,MC.size))
  Vnt <- matrix(rnorm(N*MC.size,sd = sqrt(sig2)),nrow = N, ncol = MC.size)
  VnTMC <- NULL
  X.rep <- NULL
  for (j in 1:MC.size) {
    Yn1 <- solve((diag(N)-kronecker(diag(T),lambda0*Wn.a)))%*%(rep(Cn0[,j],T) + Xn*beta + solve((diag(N)-kronecker(diag(T),rho*Mn)))%*%Vnt[,j])
    Y[,,j] <- matrix(Yn1,nrow=n,ncol=T)-matrix(rep(rowMeans(matrix(Yn1,nrow=n,ncol=T)),T),n,T)
    V[,,j] <- Rn%*%(Sn %*% Y[,,j]-X.tilde.nt*beta)
    VnTMC <- rbind(VnTMC,V[,,j]) 
    X.rep <- rbind(X.rep,X.tilde.nt)
  }
  
  for (i in 1:n) {
    
    # Y <- array(dim = c(n,T,MC.size))
    # V <- array(dim = c(n,T,MC.size))
   
     
      
        #psi.t[,t] <- c(1/(sig^2)*((Rn%*%Gn%*%solve(Rn)%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Gn[i,i],1/(sig^2)*((Hn%*%V[,t,j])[i])*V[i,t,j]-(T-1)/T*Hn[i,i])
    A <- Rn[i,]%*%(Gn%*%solve(Rn))

    B <- Rn[i,]%*%Gn*beta
    psi.T1 <- 1/(sig2)*diag(kronecker(diag(MC.size),t(Rn[i,]))%*%X.rep%*%V[i,,])
    psi.T2 <- 1/(sig2)*diag(kronecker(diag(MC.size),A)%*%VnTMC%*%V[i,,])+1/(sig2)*diag(kronecker(diag(MC.size),B)%*%X.rep%*%V[i,,])-rep((T-1)*Gn[i,i],MC.size)
    psi.T3 <- 1/(sig2)*diag(kronecker(diag(MC.size),t(Hn[i,]))%*%VnTMC%*%V[i,,])-rep((T-1)*Hn[i,i],MC.size)
    psi.T4 <- 1/(2*sig2*sig2)*diag(t(V[i,,])%*%V[i,,])-rep((T-1)/(2*sig2),MC.size)
    psi.T <- rbind(psi.T1,psi.T2,psi.T3,psi.T4)
    exp.T <- exp(t(nu)%*%psi.T)
    
    ln.T[i] <- log(sum(exp.T)/MC.size)
  }
  K.psi <- sum(ln.T)/n    #### didn't add minus sign
  return(K.psi)
}

###### inf sup K #######
set.seed(1234)
start.value <- c(0,0,0,0)
K.sup <- function(lambda, theta2){
   # nlm <- nlm(f= K.psi.neg1, p = runif(4,-1,1), lambda = lambda, theta2 = theta2, stepmax = 0.5,
   #            steptol = 0.01,iterlim = 10)
   #nlm <- optim(par = runif(4,-1,1),fn = K.psi.neg1,  lambda = lambda, theta2 = theta2, method = "BFGS")
    nlm <- nlminb(start = start.value, objective = K.psi.neg1,lambda = lambda, theta2 = theta2,
                  control = list(rel.tol= 1e-2,iter.max=50,eval.max = 50))
  
  # k.sup <- -nlm$minimum   #nlm
  # nu.est <- nlm$estimate
  
   # k.sup <- -nlm$value    #optim
   # nu.est <- nlm$par
  
  k.sup <- -nlm$objective
  nu.est <- nlm$par
  # print(nu.est)
  return(k.sup)  # add minus sign
}

# start.value1 <- matrix(data = runif(4*2,-0.1,0.1),ncol = 4)
# out <-  apply(start.value1, 1, nlminb, objective=K.psi.neg1, lambda = lambda.hat, theta2 = c(beta.hat,rho.hat,sig2.hat),
#               control = list(rel.tol= 1e-2))
# str(out)
# 
# out1 <-  apply(start.value, 1, optim, fn=K.psi.neg1, lambda = lambda.hat, theta2 = c(beta.hat,rho.hat,sig2.hat),
#                method = "BFGS")
# str(out1)
# 
# out2 <-  apply(start.value, 1, optim, fn=K.psi.neg1, lambda = lambda.hat, theta2 = c(beta.hat,rho.hat,sig2.hat),
#                )
# str(out2)

K.inf <- function(lambda){
  f <- nlminb(start = c(beta.hat,rho.hat,sig2.hat), objective = K.sup, lambda = lambda,
              lower = c(-9999,-0.99,0.001), upper = c(9999,0.99,9999),
              control = list(rel.tol= 1e-2,iter.max=50,eval.max = 50))
  k.inf <- f$objective
  theta2.est <- f$par
  # f <- optim(par = c(0,0,0.01), fn = K.sup, method = "L-BFGS-B", lower = c(-9999,-0.99,0.001), upper = c(9999,0.99,9999))
  # k.inf <- f$value
  # theta <- f$par
  
  return(c(k.inf,theta2.est))
}


SIGMA <- function(beta,lambda,rho,sig2){
  Sn <- diag(n)-lambda*Wn.a 
  Rn <- diag(n)-rho*Mn
  Gn <- Wn.a%*%solve(Sn)
  Hn <- Mn%*%solve(Rn)
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

lamb.cov <- (solve(SIGMA(beta.hat,lambda.hat,rho.hat,sig2.hat)))[2,2]/m
test.wald <- (lambda.hat-lambda0)^2/lamb.cov

ptm <- proc.time()

k.inf <- K.inf(lambda.hat)
h.hat <- k.inf[1]
theta2.est <- k.inf[-1]
test.stat <- 2*n*h.hat

proc.time() - ptm

#############################################################
# Table 3: SARAR(1,1) model: p-values of Saddlepoint (SADn) #
# and Wald (ASY) tests for several composite hypotheses.    #
#############################################################
p.sad <- pchisq(test.stat,3,lower.tail = F)
p.wald <- pchisq(test.wald,3,lower.tail = F)

