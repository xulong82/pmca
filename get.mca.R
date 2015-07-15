
get.mca <- function(X,Y) {
  if (ncol(X) != ncol(Y) & nrow(X)==nrow(Y)) {X <- t(X); Y=t(Y)} 
  X.s <- t(scale(t(X),scale=F))
  Y.s <- t(scale(t(Y),scale=F))
  n1 <- nrow(X)
  n2 <- nrow(Y)
  N <- ncol(X)
  C <- X.s%*%t(Y.s)/N
  svd <- svd(C)
  u <- svd$u
  v <- svd$v
  sigma <- svd$d
  Px <- t(u)%*%X.s
  Zx <- X%*%t(Px)
  Zy <- Y%*%t(Px)
  Zx <- t(scale(t(Zx),scale=T,center=F))
  Zy <- t(scale(t(Zy),scale=T,center=F))
  results <- list(u=u,v=v,sigma=sigma,Px=Px,Zx=Zx,Zy=Zy)
  return(results)
}
