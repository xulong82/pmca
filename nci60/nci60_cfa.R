rm(list = ls())
setwd("~/Dropbox/GitHub/pmca/nci60")

rm(list = ls())
load("nci60.rdt")

X = nci60$rna
Y = nci60$protein

Xp = X - rowMeans(X)
Yp = Y - rowMeans(Y)
C = Xp %*% t(Yp) / (ncol(X) - 1)

svd.c = svd(C) # decompose the covariance
u = svd.c$u # left singular vector
v = svd.c$v # right signular vector
d = svd.c$d # cross-covariance

pdf("variance.pdf", width = 6, height = 4)
plot(cumsum(d^2)/sum(d^2), type = "b", ylim = c(0, 1), ylab = "Cumulative variance (%)") 
text(x = 1:48, y = cumsum(d^2)/sum(d^2), labels = 1:48, col = "red", adj = -1)
dev.off()

barplot(u[, 1])
barplot(v[, 1])

# column vectors of u correspond to the structures in X that explain the covariance
# column vectors of v correspond to the structures in Y that explain the covariance

A = t(Xp) %*% u # projections of X onto u
B = t(Yp) %*% v # projections of Y onto v

plot(scale(A[, 1]), type = "b")
lines(scale(B[, 1]), type = "b", col = "red")

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
cfa_permute_cutoff <- apply(cfa_permute, 1, function(x) quantile(x, 0.999))

n.pair = min(which(cfa_permute_cutoff > d)) - 1

p1 = u[, 1]
names(p1) = rownames(X)
(p1 = sort(abs(p1), decreasing = T))
x1 = X[names(p1), ] 

x2 = sapply(1:nrow(x1), function(x) p1[1:x] %*% x1[1:x, ])
x2 = apply(x2, 2, function(x) cor(x, x2[, ncol(x2)]))
plot(x2)
abline(v = 0.05 * length(x2))
