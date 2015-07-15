

get.inter <- function(g,J=by) {
  inter <- array(dim=c(length(g),length(g)))
  for (i in 1:length(g)) {
    for (j in 1:length(g)) {
      inter[i,j] <- length(intersect(g[[i]][[J]],g[[j]][[J]]))
    }
  }
  inter
}