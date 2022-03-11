
pair_indices_all <- function(n, k) {
  n_total <- choose(n, k)
  
  res <- matrix(nrow = n_total, ncol = k)
  res[1, ] <- seq_len(n)[1:2]
  
  idx_pair <- 2L
  e <- 0
  h <- k
  a <- seq_len(k)
  while (a[1L] != n-k+1) {
    if (e < n-h) {
      h <- 1L
      e <- a[k]
      j <- 1L
    } else {
      e <- a[k-h]
      h <- h+1L
      j <- 1L:h
    }
    a[k-h+j] <- e+j
    res[idx_pair, ] <- x[a]
    idx_pair <- idx_pair+1L
  }
  res
}



