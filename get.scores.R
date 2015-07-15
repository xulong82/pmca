
get.scores <- function(Zx,Zy) {
  p <- nrow(Zx)
  q <- nrow(Zy)
  n <- ncol(Zx)-1
  score <- array(dim=c(q,n,p))
  for (i in 1:p) {
    for (j in 1:n) {
      zvec <- rep(Zx[i,j],q)
      score[,j,i] <- abs(zvec-Zy[,j])
    }
  }
  score
}