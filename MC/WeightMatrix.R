#################### Figure 1: Different types of weight matrices for n=24 and n=100 ##################

library("splm") 

#################### n is 24 ##################

r1 <- 6  ### number of rows in the grid
c1 <- 4  ### number of columns in the grid
##### Room type
nb3rt <- cell2nb(r1, c1)  
xyc <- attr(nb3rt, "region.id")
xy <- matrix(as.integer(unlist(strsplit(xyc, ":"))), ncol=2, byrow=TRUE)

##### Queen type
nb3q <- cell2nb(r1, c1,type = "queen")  
xycq <- attr(nb3q, "region.id")
xyq <- matrix(as.integer(unlist(strsplit(xycq, ":"))), ncol=2, byrow=TRUE)

##### Queen with torus type
nb3qt <- cell2nb(r1, c1,type = "queen",torus = T)   
xycqt <- attr(nb3qt, "region.id")
xyqt <- matrix(as.integer(unlist(strsplit(xycqt, ":"))), ncol=2, byrow=TRUE)

##### Plot Rook, Queen, Queen with torus weight matrices for n is 24
par(mfrow=c(1,3))
plot(nb3rt, xy)
plot(nb3q, xyq)
plot(nb3qt, xyqt)


################### n is 100 ####################

r2 <- 10
c2 <- 10

##### Rook type
nb10rt <- cell2nb(r2, c2)   
xyc1 <- attr(nb10rt, "region.id")
xy1 <- matrix(as.integer(unlist(strsplit(xyc1, ":"))), ncol=2, byrow=TRUE)

##### Queen type
nb10q <- cell2nb(r2, c2,type = "queen")  
xycq1 <- attr(nb10q, "region.id")
xyq1 <- matrix(as.integer(unlist(strsplit(xycq1, ":"))), ncol=2, byrow=TRUE)

##### Queen with torus type
nb10qt <- cell2nb(r2, c2,type = "queen",torus = T)   
xycqt1 <- attr(nb10qt, "region.id")
xyqt1 <- matrix(as.integer(unlist(strsplit(xycqt1, ":"))), ncol=2, byrow=TRUE)

##### Plot Rook, Queen, Queen with torus weight matrices for n is 100
par(mfrow=c(1,3))
plot(nb10rt, xy1)
plot(nb10q, xyq1)
plot(nb10qt, xyqt1)
