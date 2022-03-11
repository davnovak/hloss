
ComputeDistances <- function(tax, leaves_only = FALSE) UseMethod('ComputeDistances', tax)

ComputeDistances.Taxonomy <- function(tax, leaves_only = FALSE) {
  nodes <- if (leaves_only) tax$node_paths[tax$is_leaf] else tax$node_paths
  n <- length(nodes)
  
  ## Compute pair-distances
  M <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:i) {
      if (i != j) {
        M[i, j] <- M[j, i] <- Path(tax, from = nodes[i], to = nodes[j], only_total = TRUE, weighted = TRUE)
      }
    }
  }
  M <- M / max(M, na.rm = TRUE)
  colnames(M) <- rownames(M) <- nodes
  
  ## Compute self-distances (local tree densities)
  D <- vector(mode = 'numeric', length = n)
  for (i in 1:n) {
    node <- nodes[i]
    d <- M[node, ]
    d <- d[!is.na(d)]
    D[i] <- -sum(log(d))
  }
  D <- D / max(D)
  diag(M) <- D
  tax$d <- M
  tax$changes_made <- FALSE
  invisible(tax)
}



