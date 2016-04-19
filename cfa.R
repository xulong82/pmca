rm(list = ls())
setwd("~/Dropbox/GitHub/pmca")
load("./dataForPMCA_chosenGenes.RData")

(gene = colnames(X)) # genes
(cell = rownames(Y)) # single cell RNA-seq
(sample = rownames(X)) # time serial AD bulk RNA-seq

geno = gsub(".*(WT|APP).*", "\\1", sample)
month = gsub(".*(2m|4m|5m|6m).*", "\\1", sample)
group = c("WT2m", "WT4m", "WT5m", "WT6m", "APP2m", "APP4m", "APP5m", "APP6m")
group = factor(paste0(geno, month), levels = group)

X = t(sapply(levels(group), function(x) colMeans(X[group == x, ])))

Xp = X - rowMeans(X)
Yp = Y - rowMeans(Y)
C = Xp %*% t(Yp) / (ncol(X) - 1)

svd.c = svd(C) # decompose the covariance
u = svd.c$u # left singular vector
v = svd.c$v # right signular vector
d = svd.c$d # cross-covariance

pdf("variance.pdf", width = 6, height = 4)
plot(cumsum(d^2)/sum(d^2), type = "b", ylim = c(0, 1), ylab = "Cumulative variance (%)") 
text(x = 1:8, y = cumsum(d^2)/sum(d^2), labels = 1:8, col = "red", adj = -1)
dev.off()

plot(scale(u[, 1]), type = "b")
plot(scale(v[, 1]), type = "b")
text(x = 1:48, y = scale(v[, 1]), labels = 1:48, col = "red", adj = 1)

cell[c(37, 38, 39, 45, 46)]

A = t(Xp) %*% u
B = t(Yp) %*% v

plot(scale(A[, 1]), type = "b")
lines(scale(B[, 1]), type = "b", col = "red")
abline(h = -1)

gene[which(abs(scale(B[, 1])) > 1)]

plot(scale(A[, 2]), type = "b")
lines(scale(B[, 2]), type = "b", col = "red")

cor(A[, 1], B[, 1])
cor(A[, 2], B[, 2])

my_cfa_permute <- function(X, Y) {
  Yp = Y[, sample(1:ncol(Y))]
  C = X %*% t(Yp) / (ncol(X) - 1)
  svd(C)$d
}

cfa_permute <- replicate(1e3, my_cfa_permute(Xp, Yp))
(cfa_permute_cutoff <- apply(cfa_permute, 1, function(x) quantile(x, 0.99)))

n.pair = min(which(cfa_permute_cutoff > d)) - 1

p1x = u[, 1]
names(p1x) = rownames(X)
(p1x = sort(abs(p1x), decreasing = T))
x1 = X[names(p1x), ] 

x2 = sapply(1:nrow(x1), function(x) p1x[1:x] %*% x1[1:x, ])
x2 = apply(x2, 2, function(x) cor(x, x2[, 8]))
plot(x2, type = "b", ylim = c(0, 1))
text(x = 1:8, y = x2, labels = rownames(x1), col = "red", adj = 0)
p1 = v[, 1]
names(p1) = cell
(p1 = sort(abs(p1), decreasing = T))
y1 = Y[names(p1), ] 

y2 = sapply(1:nrow(y1), function(x) p1[1:x] %*% y1[1:x, ])
y2 = apply(y2, 2, function(x) cor(x, y2[, 48]))
plot(y2, type = "b", ylim = c(0, 1))
text(x = 1:48, y = y2, labels = 1:48, col = "red", adj = 0)
