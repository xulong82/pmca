

match.patterns <- function(score,w,rowterms=rownames(Y)) {
  n <- length(w)
  if (length(dim(score))==2) {
    g <- vector("list",n)
    g[[1]] <- rowterms[which(score[,1]<=w[1])]
    for (j in 2:n) {
      g[[j]] <- intersect(g[[j-1]],rowterms[which(score[,j]<=w[j])])
    }
  } else {
    p <- dim(score)[3]
    g <- vector("list",p)
    for (i in 1:p) {
      g[[i]] <- vector("list",n)
      g[[i]][[1]] <- rowterms[which(score[,1,i]<=w[1])]
      for (j in 2:n) {
        g[[i]][[j]] <- intersect(g[[i]][[j-1]],rowterms[which(score[,j,i]<=w[j])])
      }
    }
  }
  g
}