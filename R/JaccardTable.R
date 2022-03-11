
.JaccardTable <- function(
  c1, c2
) {
  j <- matrix(NA, nrow = nlevels(c1), ncol = nlevels(c2), dimnames = list(levels(c1), levels(c2)))
  for (idx_c1 in levels(c1)) {
    for (idx_c2 in levels(c2)) {
      t1 <- c1 == idx_c1
      t2 <- c2 == idx_c2
      j[idx_c1, idx_c2] <- sum(t1 & t2) / sum(t1 | t2)
    }
  }
  j
}
