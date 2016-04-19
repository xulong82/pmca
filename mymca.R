rm(list = ls()); setwd("~/Dropbox/GitHub/pmca")
load("./dataForPMCA_chosenGenes.RData")

(gene = colnames(X)) # genes
(cell = rownames(Y)) # single cell RNA-seq
(sample = rownames(X)) # time serial AD bulk RNA-seq

geno = gsub(".*(WT|APP).*", "\\1", sample)
month = gsub(".*(2m|4m|5m|6m).*", "\\1", sample)
group = c("WT2m", "WT4m", "WT5m", "WT6m", "APP2m", "APP4m", "APP5m", "APP6m")
group = factor(paste0(geno, month), levels = group)

wt2m = colMeans(X[group == "WT2m", ])
X = sweep(X, 2, wt2m, "-") # WT:2m reference

X = t(sapply(levels(group), function(x) colMeans(X[group == x, ])))
X = sweep(X, 2, X[1, ], "-")[-1, ] # WT:2m reference

Xp = X - rowMeans(X)
Yp = Y - rowMeans(Y)
C = Xp %*% t(Yp) / (ncol(X) - 1)

svd.c = svd(C) # decompose the covariance
u = svd.c$u # left singular vector
v = svd.c$v # right signular vector
d = svd.c$d # cross-covariance

A = t(Xp) %*% u 
B = t(Yp) %*% v

pdf("variance.pdf", width = 6, height = 4)

plot(cumsum(d^2)/sum(d^2), type = "b", ylim = c(0, 1), ylab = "Cumulative variance (%)") 
text(x = 1:8, y = cumsum(d^2)/sum(d^2), labels = 1:8, col = "red", adj = c(0, 1))

dev.off()

pdf("mode.pdf", width = 8, height = 5)

plot(scale(u[, 1]), type = "b", xlab = "bulk sample")
text(x = 1:8, y = scale(u[, 1]), labels = rownames(X), col = "red", adj = 0)
plot(scale(v[, 1]), type = "b", xlab = "cell")
text(x = 1:48, y = scale(v[, 1]), labels = rownames(Y), col = "red", adj = 1)

plot(scale(A[, 1]), type = "b", xlab = "gene", ylab = "1st A/B vectors", ylim = c(-8, 3))
lines(scale(B[, 1]), type = "b", col = "red")
lab = gene; gene[which(abs(scale(B[, 1])) < 1)] = ""
text(x = 1:108, y = scale(B[, 1]), labels = lab, col = "red", adj = 1)

dev.off()

Zx = Xp %*% A
Zx <- t(scale(t(Zx), scale = T, center = F))

Zy = Y %*% A
Zy <- t(scale(t(Zy), scale = T, center = F)) 

source("pmca.R")
source("get.mca.R")
source("get.inter.R")
source("get.scores.R")
source("match.patterns.R")
source("iterative.proc.R")
source("permutation.proc.R")

by = 2 # by: by which component do you want the FPR <= alpha
method = "each" # method: "overall" if you want an overall FPR, "each" if you want a FPR for each rowterm of X
B = 1e3 # B: number of permutations
alpha = .05 # alpha: want the FPR <= alpha

my.scores.pos <- lapply(1:by, function(x) { sapply(Zx[, x], function(y) abs(y - Zy[, x])) })

scores.pos <- get.scores(Zx, Zy) # distance
scores.neg <- get.scores(-Zx, Zy) # anti-associated

# score[i,k,j]: distance of i-th cell in Y and j-th sample in X in regard to k-th principle

set.seed(B)
scores.ran <- permutation.proc(X, Y, method = method, B = B) # get scores for permutations

summary(scores.ran[1, 1, 1, ])
quantile(scores.ran[1, 1, 1, ], 0.05)
scores.pos[1, 1, 1]

pdf("density.pdf", width = 7, height = 5)
par(mfrow = c(2, 3))
lapply(1:6, function(x) plot(density(scores.ran[1, x, 1, ]), main = paste0("cell 1 to sample 1 on PC", x)))
dev.off()

# shuffle the dynamics for each cell over all genes, compute score, and repeat

w <- apply(Zx, 2, sd) # get starting window vector
w = w^2
w = d

# stringencies inverse-correlate with principle component variations

# get wopt (optimal window vector)
# tau: controls the width of the window (w/tau, larger tau = more strict) 

set.seed(B)
it.result <- iterative.proc(scores.ran, alpha, w, method = method, by = by, tau=0.3) 
it.result$tau

g <- match.patterns(scores.pos, w=it.result$wopt) # associated
g2 <- match.patterns(scores.neg, w=it.result$wopt) # anti-associated

# g[[i]][[j]]: rownames(Y) that match the pattern of rownames(X)[i] when the component = j
# example: all the cellTypes map to "APP5m1558.2014" (row 3) when you use 4 components, g[[3]][[4]]

# get.inter(g, J=by): returns the intersections across all rownames(X)
# inter[i,j] = number of rownames(Y) that mapped to both rownames(X)[i] & rownames(X)[j]
by = 2
inter = get.inter(g, by)
rownames(inter) <- colnames(inter) <- rownames(X)
inter <- inter[which(rowSums(inter) > 0), ]
inter <- inter[, which(colSums(inter) > 0)]
inter

mapped <- data.frame(sample = rownames(inter), cellType = rep("nn", nrow(inter)))
mapped$cellType = as.character(mapped$cellType)

for (i in 1:nrow(inter)) {
  d <- paste(g[[which(rownames(X)==rownames(inter)[i])]][[by]], collapse=",")  
  mapped$cellType[i] <- d
}

mapped
