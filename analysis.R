rm(list = ls())
setwd("~/Dropbox/GitHub/pmca")
load("dataForPMCA_chosenGenes.RData")

cell = rownames(Y) # single cell RNA-seq
sample = rownames(X) # time serial AD bulk RNA-seq

geno = gsub(".*(WT|APP).*", "\\1", sample)
month = gsub(".*(2m|4m|5m|6m).*", "\\1", sample)
group = c("WT2m", "WT4m", "WT5m", "WT6m", "APP2m", "APP4m", "APP5m", "APP6m")
group = factor(paste0(geno, month), levels = group)

# X matrix as per sample
wt2m = colMeans(X[group == "WT2m", ])
X = sweep(X, 2, wt2m, "-") # WT:2m reference

# X matrix as per group
X = t(sapply(levels(group), function(x) colMeans(X[group == x, ])))

X = t(sapply(levels(group), function(x) { # median / mad: not stable??
  apply(X[group == x, ], 2, function(xx) median(xx) / (1.4826 * mad(xx)))
}))

X = sweep(X, 2, X[1, ], "-")[-1, ] # WT:2m reference

# MCA
Xp = X - rowMeans(X)
Yp = Y - rowMeans(Y)
# C = Xp %*% t(Yp) / ncol(X) # R's
C = Xp %*% t(Yp) / (ncol(X) - 1)

svd.c = svd(C)
u = svd.c$u # left singular vector
v = svd.c$v # right signular vector
d = svd.c$d # cross-covariance

plot(cumsum(d)/sum(d), type = "b", ylim = c(0, 1)) 
text(x = 1:48, y = cumsum(d)/sum(d), labels = 1:48, col = "red", adj = -1)

A = t(Xp) %*% u # projections of X onto u
B = t(Yp) %*% v # projections of Y onto v

# A and B are linear combinations of X and Y, and vice versa
# total covariance of A and B are the same as total covariance of X and Y
# cov(A, B) is diagnal, and is the largest in the first component

# Signular vectors (u, v) describe patterns of anomalies in each dataset that tend to occur simultaneously
# Coefficients (A, B) represent amplitudes of the respective patterns in each sample

Zx = X %*% A
Zx <- t(scale(t(Zx),scale = T, center = F))

Zy = Y %*% A # Robyn
Zy <- t(scale(t(Zy),scale = T, center = F))

# Zy = Y %*% B # MCA

source("pmca.R")
source("get.mca.R")
source("get.inter.R")
source("get.scores.R")
source("match.patterns.R")
source("iterative.proc.R")
source("permutation.proc.R")

# by: by which component do you want the FPR <= alpha
# method: "overall" if you want an overall FPR, "each" if you want a FPR for each rowterm of X
# B: number of permutations
# alpha: want the FPR <= alpha
# plot: TRUE will plot one of the FPR distributions across all B permutations.

by = 2; method = "each"; B = 1000; alpha = .05; plot = TRUE;

scores.pos <- get.scores(Zx, Zy) # distance
scores.neg <- get.scores(-Zx, Zy) # anti-associated
# score[i,k,j]: distance of i-th cell in Y and j-th sample in X in regard to k-th principle

set.seed(B)
scores.ran <- permutation.proc(X, Y, method = method, B = B) # get scores for permutations

w <- apply(Zx, 2, sd) # get starting window vector
# stringencies inverse-correlate with principle component variations

# get wopt (optimal window vector)
# tau: controls the width of the window (w/tau, larger tau = more strict) 
set.seed(B)
it.result <- iterative.proc(scores.ran, alpha, w, method = method, by = by, plot=plot, tau=0.3) 
it.result$tau

g <- match.patterns(scores.pos, w=it.result$wopt) # associated
g <- match.patterns(scores.neg, w=it.result$wopt) # anti-associated

# g[[i]][[j]]: rownames(Y) that match the pattern of rownames(X)[i] when the component = j
# example: all the cellTypes map to "APP5m1558.2014" (row 3) when you use 4 components, g[[3]][[4]]

# get.inter(g, J=by): returns the intersections across all rownames(X)
# inter[i,j] = number of rownames(Y) that mapped to both rownames(X)[i] & rownames(X)[j]
J = 1
inter = get.inter(g, J)
rownames(inter) <- colnames(inter) <- rownames(X)
inter <- inter[which(rowSums(inter) > 0), ]
inter <- inter[, which(colSums(inter) > 0)]
inter

mapped <- data.frame(sample = rownames(inter), cellType = rep("nn", nrow(inter)))
mapped$cellType = as.character(mapped$cellType)

for (i in 1:nrow(inter)) {
  d <- paste(g[[which(rownames(X)==rownames(inter)[i])]][[J]], collapse=",")  
  mapped$cellType[i] <- d
}

mapped
