#################### n is 24

r1 <- 6
c1 <- 4
nb3rt <- cell2nb(r1, c1)
xyc <- attr(nb3rt, "region.id")
xy <- matrix(as.integer(unlist(strsplit(xyc, ":"))), ncol=2, byrow=TRUE)
nb3q <- cell2nb(r1, c1,type = "queen")
xycq <- attr(nb3q, "region.id")
xyq <- matrix(as.integer(unlist(strsplit(xycq, ":"))), ncol=2, byrow=TRUE)
nb3qt <- cell2nb(r1, c1,type = "queen",torus = T)
xycqt <- attr(nb3qt, "region.id")
xyqt <- matrix(as.integer(unlist(strsplit(xycqt, ":"))), ncol=2, byrow=TRUE)
par(mfrow=c(1,3))
plot(nb3rt, xy)
plot(nb3q, xyq)
plot(nb3qt, xyqt)


################### n is 100

r2 <- 10
c2 <- 10

nb10rt <- cell2nb(r2, c2)
xyc1 <- attr(nb10rt, "region.id")
xy1 <- matrix(as.integer(unlist(strsplit(xyc1, ":"))), ncol=2, byrow=TRUE)
nb10q <- cell2nb(r2, c2,type = "queen")
xycq1 <- attr(nb10q, "region.id")
xyq1 <- matrix(as.integer(unlist(strsplit(xycq1, ":"))), ncol=2, byrow=TRUE)
nb10qt <- cell2nb(r2, c2,type = "queen",torus = T)
xycqt1 <- attr(nb10qt, "region.id")
xyqt1 <- matrix(as.integer(unlist(strsplit(xycqt1, ":"))), ncol=2, byrow=TRUE)
par(mfrow=c(1,3))
plot(nb10rt, xy1)
plot(nb10q, xyq1)
plot(nb10qt, xyqt1)