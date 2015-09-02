scriptdir <- "/Volumes/FreeAgent GoFlex Drive/jax/meiosis/scripts/PMCA/"
setwd(scriptdir)
source("get.mca.R")
source("get.inter.R")
source("get.scores.R")
source("iterative.proc.R")
source("permutation.proc.R")
source("pmca.R")
source("match.patterns.R")

setwd("/Volumes/FreeAgent GoFlex Drive/jax/Xulong/gareth/")
load("dataForPMCA_chosenGenes.RData")

mca.real <- get.mca(X,Y) # get Zx, Zy, sigma, etc.
plot(mca.real$sigma) # plot of singular values of the covariance matrix

# by: by which component do you want the FPR <= alpha
# method: "overall" if you want an overall FPR, "each" if you want a FPR for each rowterm of X
# B: number of permutations
# alpha: want the FPR <= alpha
# plot: TRUE will plot one of the FPR distributions across all B permutations.
by=6; method="each"; B=1000; alpha=.05; plot=TRUE;

# get scores for the real data. Lower scores --> closer pattern match
scores.real <- get.scores(mca.real$Zx,mca.real$Zy)
scores.neg <- get.scores(-mca.real$Zx,mca.real$Zy) #for anti-associated lists
set.seed(B) # for reproducibility
scores.rand <- permutation.proc(X,Y,method=method,B=B) # get scores for all B permutations
w <- apply(mca.real$Zx,2,sd) # get starting window vector
set.seed(B) # for reproducibilty
# tau: controls the width of the window (w/tau) 
#     if you are not getting good results, try making tau smaller or larger (larger = more strict)
it.result <- iterative.proc(scores.rand,alpha=alpha,w=w,method=method,by=by,plot=plot,tau=1) 

# get wopt (optimal window vector)
# g = list of rownames(Y) that match patterns of rownames(X)
#   g[[i]][[j]]: is all the rownames(Y) that match the pattern of rownames(X)[i] when the component = j
#   example: If you want to know what cellTypes map to "APP5m1558.2014" (row 3) 
#         when you use 4 components, g[[3]][[4]]
g <- match.patterns(scores.real,w=it.result$wopt)  

## get.inter(g,J=by): will return the intersections across all rownames(X)
#   so inter[i,j] = number of rownames(Y) that mapped to both rownames(X)[i] & rownames(X)[j]
get.inter(g,4)

## same thing for h but h is the anti-associated rowterms(Y)
h <- match.patterns(scores.neg,w=it.result$wopt)  

## false positive rate for all components
fpr <- it.result$FPR
if (!is.null(dim(fpr))) {
  rownames(fpr) <- rownames(X)
  colnames(fpr) <- paste0("component=",1:ncol(fpr))
} else {
  names(fpr) <- paste0("component=",1:length(fpr))
}
names(it.result$wopt) <- paste0("component=",1:length(it.result$wopt))
result <- list(mapped=g,antimapped=h,FPR=fpr,wopt=it.result$wopt,tau=it.result$tau) # RESULTS

## use this to get interesting intersections (vary J)
inter = get.inter(g,5)
rownames(inter) <- colnames(inter) <- rownames(X)
inter <- inter[which(rowSums(inter)>0),]
inter <- inter[,which(colSums(inter)>0)]
inter

mapped.mat <- data.frame(sample=rownames(inter),cellType=rep("nn",nrow(inter)))
mapped.mat$cellType <- as.character(mapped.mat$cellType)
J <- 5 #component
for (i in 1:nrow(inter)) {
  d <- paste(g[[which(rownames(X)==rownames(inter)[i])]][[J]],collapse=",")  
  mapped.mat$cellType[i] <- d
}
mapped.mat
save(result,mapped.mat,file=paste0("PMCAresult_",Sys.Date(),".RData"))

