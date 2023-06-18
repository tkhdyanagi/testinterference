exposure_a <- function(Z, A, i) {
  Z[i]
}

exposure_b <- function(Z, A, i) {
  as.integer(Z + A %*% Z > 0)[i]
}

exposure_c <- function(Z, A, i) {
  cbind(Z, as.integer(A %*% Z > 0))[i, ]
}

exposure_d <- function(Z, A, i) {
  cbind(Z, A %*% Z)[i, ]
}

constant <- function(Z, A, i) {
  rep(0, length(i))
}
