

pmca <- function(X,Y,method="overall",alpha=.05,B=1000,by=NULL,plot=TRUE) {
  mca.real <- get.mca(X,Y)
  if (is.null(by)) {
    by <- which(mca.real$sigma^2<=.10)[1]
  }
  scores.real <- get.scores(mca.real$Zx,mca.real$Zy)
  scores.neg <- get.scores(-mca.real$Zx,mca.real$Zy)
  set.seed(B)
  method="each"
  scores.rand <- permutation.proc(X,Y,method=method,B=B)
  w <- apply(mca.real$Zx,2,sd)
  set.seed(B)
  it.result <- iterative.proc(scores.rand,alpha=alpha,w=w,method=method,by=by,plot=plot)
  g <- match.patterns(scores.real,w=it.result$wopt)  
  get.inter(g)
  h <- match.patterns(scores.neg,w=it.result$wopt)  
  fpr <- it.result$FPR
  if (!is.null(dim(fpr))) {
    rownames(fpr) <- rownames(X)
    colnames(fpr) <- paste0("component=",1:ncol(fpr))
  } else {
    names(fpr) <- paste0("component=",1:length(fpr))
  }
  names(it.result$wopt) <- paste0("component=",1:length(it.result$wopt))
  result <- list(mapped=g,antimapped=h,FPR=fpr,wopt=it.result$wopt,tau=it.result$tau)
}