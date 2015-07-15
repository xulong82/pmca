

iterative.proc <- function(scores,alpha=.05,w,tau=1,method="overall",by=3,rowterms=rownames(Y),plot=TRUE) {
  B <- dim(scores)[4]; p <- dim(scores)[3]; J <- dim(scores)[2]; q <- dim(scores)[1]
  w <- w[1:J]
  if (method=="overall") {
    fpr <- array(dim=c(J,B))
    done <- FALSE
    while (!done) {
      wopt <- w/tau
      for (b in 1:B) {
        g <- match.patterns(scores[,,1,b],w=wopt)
        fpr[,b] <- unlist(lapply(g,length))/q
      }
      if (mean(fpr[by,])<=alpha) {
        done=TRUE
      } else {
        tau=tau+.1
      }
    }
    if (plot) {
      tit <- paste0("Histogram of fpr for Component ",by)
      hist(fpr[by,],xlab=paste0("estimated FPR for B=",B),main=tit)
      abline(v=alpha,col="red",lwd=2)
    }
    FPR <- rowMeans(fpr)
  } else { #FPR for each rowterm of X
    fpr <- array(dim=c(p,J,B))
    FPR <- array(dim=c(p,J))
    done <- FALSE
    while (!done) {
      wopt <- w/tau
      for (b in 1:B) {
        g <- match.patterns(scores[,,,b],w=wopt)
        for (i in 1:p) {
          fpr[i,,b] <- unlist(lapply(g[[i]],length))/q
        }
      }
      for (i in 1:p) {
        FPR[i,] <- rowMeans(fpr[i,,])
      }
      if (all(FPR[,by]<=alpha)) {
        done=TRUE
      } else {
        tau <- tau + .1
      }
    }
    if (plot) {
      tit <- paste0("Histogram of FPR for Rowterm ",p," and Component ",by)
      hist(fpr[i,by,],xlab=paste0("estimated FPR for B=",B),main=tit)
      abline(v=alpha,col="red",lwd=2)
    }
  }
  result <- list(FPR=FPR,wopt=wopt,tau=tau)
  result
} 