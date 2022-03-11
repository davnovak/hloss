
#' Compute Rand index
#'
#' @param pred numeric or character vector of predicted clusters/populations
#' @param true numeric or character vector of predicted clusters/populations
#'
#' @export
RI <- function(pred, true) {
  pred <- as.numeric(as.factor(pred))
  true <- as.numeric(as.factor(true))
  
  diff_pred <- abs(sapply(pred, function(x) x - pred))
  diff_pred[diff_pred > 1] <- 1
  
  diff_true <- abs(sapply(true, function(x) x - true))
  diff_true[diff_true > 1] <- 1
  
  1 - (sum(abs(diff_pred-diff_true))/2) / (choose(nrow(diff_pred), 2))
}

#' Compute Hierarchy-aware Rand index
#'
#' Computes the Hierarchy-aware Rand index (hRI) based on a vector of cluster assignment, a vector of ground-truth labels and a \code{Taxonomy} object.
#' 
#' Since computation of hierarchical RI is expensive, exact solution is computed for cluster/label vectors of up to `max_npt` points.
#' If the number of points is higher, an average value for multiple runs on random subsets of the vectors is computed.
#' Size of the subset is then determined by the parameter `approx_batch` and the number of runs is given by `approx_niter`.
#'
#' @param tax object of class \code{Taxonomy}
#' @param pred numeric or character vector of predicted clusters
#' @param true character vector of true node labels from \code{Taxonomy}
#' @param max_npt integer: maximum number of points for which to compute an exact solution. Default value is \code{1000}
#' @param approx_batch integer: size of random subset per run for approximate computation of the index value. Default value is \code{500}
#' @param approx_niter integer: number of random subsets (runs) to be taken for approximate computation of the index value. Default value is \code{5000}
#'
#' @export
hRI <- function(tax, pred, true, max_npt = 1000, approx_batch = 500, approx_niter = 5000) {
  if (is.null(tax$d)) {
    stop('Pair- and self-distances need to be computed first using ComputeDistanes')
  }
  if (tax$changes_made) {
    stop('Changes were made since last computation of pair- and self-distances; recompute them first')
  }
  n <- length(pred)
  if (length(true) != n) {
    stop('Vectors pred and true need to be the same length')
  }
  pred <- as.integer(as.factor(pred))
  true <- GetPathsFromSubstrings(tax, true)
  vals <- hRI_terms(tax, pred, true, max_npt, approx_batch, approx_niter)
  (vals[1] + vals[2]) / sum(vals) # (matches) / (matches + mismatches)
}

#' Compute adjusted Rand index
#'
#' @param pred numeric or character vector of predicted clusters/populations
#' @param true numeric or character vector of predicted clusters/populations
#'
#' @export
ARI <- function(pred, true) {
  if (all(as.numeric(as.factor(pred)) == as.numeric(as.factor(true)))) {
    return(1)
  }
  pred <- as.numeric(as.factor(pred))
  true <- as.numeric(as.factor(true))
  
  counts <- table(pred, true)
  if (nrow(counts) == 1 && ncol(counts) == 1) {
    return(1)
  }
  
  aa <- sum(choose(counts, 2))
  bb <- sum(choose(rowSums(counts), 2))-aa
  cc <- sum(choose(colSums(counts), 2))-aa
  dd <- choose(sum(counts), 2)-aa-bb-cc
  
  (aa-(aa+bb)*(aa+cc)/(aa+bb+cc+dd))/((aa+bb+aa+cc)/2-(aa+bb)*(aa+cc)/(aa+bb+cc+dd))
}

#' Compute Hierarchy-aware adjusted Rand index (hARI)
#'
#' Computes the adjusted version of the Hierarchy-aware Rand index (hRI).
#' The adjusted index is computed as `AdjustedIndex = (Index - ExpectedIndex)/(MaxIndex - ExpectedIndex)`.
#'
#' The `ExpectedIndex` is calculated as hARI for a simulated random clustering.
#' The `MaxIndex` is \code{1}.
#'
#' Since computation of hierarchical (A)RI is expensive, exact solution is computed for cluster/label vectors of up to `max_npt` points.
#' If the number of points is higher, an average value for multiple runs on random subsets of the vectors is computed.
#' Size of the subset is then determined by the parameter `approx_batch` and the number of runs is given by `approx_niter`.
#'
#' @param tax object of class \code{Taxonomy}
#' @param pred vector of predicted node labels from \code{Taxonomy}
#' @param true vector of true node labels from \code{Taxonomy}
#' @param seed integer: random seed for sampling random points to calculate `ExpectedIndex`
#' @param max_npt integer: maximum number of points for which to compute an exact solution. Default value is \code{1000}
#' @param approx_batch integer: size of random subset per run for approximate computation of the index value. Default value is \code{500}
#' @param approx_niter integer: number of random subsets (runs) to be taken for approximate computation of the index value. Default value is \code{5000}
#'
#' @export
hARI <- function(tax, pred, true, seed = 42, max_npt = 1000, approx_batch = 500, approx_niter = 5000) {
  if (is.null(tax$d)) {
    stop('Pair- and self-distances need to be computed first')
  }
  if (tax$changes_made) {
    stop('Changes were made since last computation of pair- and self-distances; recompute them first')
  }
  n <- length(pred)
  if (length(true) != n) {
    stop('Vectors pred and true need to be the same length')
  }
  
  pred <- as.integer(as.factor(pred))
  true <- GetPathsFromSubstrings(tax, true)
  
  if (all(as.integer(as.factor(pred)) == as.integer(as.factor(true)))) {
    return(1)
  }
  
  Index <- hRI(tax, pred, true, max_npt, approx_batch, approx_niter)
  
  set.seed(seed)
  pred <- sample(1:tax$n, size = length(true), replace = TRUE)
  ExpectedIndex <- hRI(tax, pred, true, max_npt, approx_batch, approx_niter)
  
  MaxIndex <- 1
  
  (Index - ExpectedIndex) / (MaxIndex - ExpectedIndex)
}

RI_terms <- function(pred, true) {
  
  pred <- as.integer(as.factor(pred))
  true <- as.integer(as.factor(true))
  
  if (length(unique(pred)) == 1 && length(unique(true)) == 1) {
    return(1)
  }
  
  counts <- table(pred, true)
  aa <- sum(choose(counts, 2))
  bb <- sum(choose(rowSums(counts), 2)) - aa
  cc <- sum(choose(colSums(counts), 2)) - aa
  dd <- choose(sum(counts), 2) - aa - bb - cc
  
  c(aa, dd, cc, bb) # flipped
}

hRI_terms <- function(tax, pred, true, max_npt = 1000, approx_batch = 500, approx_niter = 5000) {
  aa <- bb <- cc <- dd <- 0
  n <- length(pred)
  
  true <- as.integer(factor(true, levels = rownames(tax$d)))
  
  single_iter <- function(tax, pred, true, n) {
    pair_idcs <- t(combn(n, 2))
    
    pred_diff <- apply(matrix(pred[pair_idcs], ncol = 2), 1, anyDuplicated) == 0
    true_diff <- apply(matrix(true[pair_idcs], ncol = 2), 1, anyDuplicated) == 0
    
    idcs_SameSame <- pair_idcs[!pred_diff & !true_diff, ] # (a) Same-Same MATCH
    idcs_DiffDiff <- pair_idcs[ pred_diff &  true_diff, ] # (b) Diff-Diff MATCH
    idcs_SameDiff <- pair_idcs[!pred_diff &  true_diff, ] # (c) Same-Diff MISMATCH
    idcs_DiffSame <- pair_idcs[ pred_diff & !true_diff, ] # (d) Diff-Same MISMATCH
    
    aa <- sum(    tax$d[matrix(true[idcs_SameSame], ncol = 2)])
    bb <- sum(1 - tax$d[matrix(true[idcs_DiffDiff], ncol = 2)])
    cc <- sum(    tax$d[matrix(true[idcs_SameDiff], ncol = 2)])
    dd <- sum(1 - tax$d[matrix(true[idcs_DiffSame], ncol = 2)])
    
    c(aa, bb, cc, dd)
  }
  
  if (n <= max_npt) { # exact solution
    return(single_iter(tax, pred, true, n))
    
  } else { # approximate solution
    res <- vector(mode = 'list', length = approx_niter)
    mult <- n / approx_batch
    for (approx_iter in approx_niter) {
      set.seed(approx_iter)
      subset_idcs <- sample(1:n, size = approx_batch, replace = FALSE)
      
      res[[approx_iter]] <- single_iter(tax, pred[subset_idcs], true[subset_idcs], approx_batch) * mult
    }
    res <- do.call(rbind, res)
    return(colMeans(res))
    
  }
}

#' Compute Accuracy
#'
#' @param pred numeric or character vector of predicted populations
#' @param true numeric or character vector of predicted populations
#'
#' @export
OverallAcc <- function(pred, true) {
  mean(pred == true)
}

#' Compute Hierarchy-aware Accuracy
#'
#' @param tax object of class \code{Taxonomy} with computed distances
#' @param pred numeric or character vector of predicted populations
#' @param true numeric or character vector of predicted populations
#'
#' @export
hOverallAcc <- function(tax, pred, true) {
  pred <- as.integer(factor(pred, levels = rownames(tax$d)))
  true <- as.integer(factor(true, levels = rownames(tax$d)))
  
  ## Treat pair-distances as penalties and self-distances as rewards
  d <- 1 - tax$d
  diag(d) <- 1 - diag(d)
    
  unadj <- mean(d[cbind(pred, true)])
  upper_bound <- mean(d[cbind(true, true)])
  
  unadj / upper_bound
}

#' Compute Accuracy per each ground-truth label
#'
#' Computes accuracy of classification per each cell population contained in the ground-truth vector of labels.
#'
#' @param tax object of class \code{Taxonomy} with computed distances
#' @param pred character vector of predicted populations
#' @param true character vector of predicted populations
#'
#' @export
AccPerMatch <- function(tax, pred, true) {
  pops <- rownames(tax$d)
  n_pops <- length(pops)
  
  res <- rep(NA, n_pops)
  names(res) <- pops
  
  for (pop in pops) {
    idcs <- which(true == pop)
    res[pop] <- sum(pred[idcs] == pop) / length(idcs)
  }
  
  res
}

#' Compute Hierarchy-aware Accuracy per each ground-truth label
#'
#' Computes accuracy of classification per each cell population contained in the ground-truth vector of labels, using a hierarchy-aware scoring scheme.
#'
#' @param tax object of class \code{Taxonomy} with computed distances
#' @param pred character vector of predicted populations
#' @param true character vector of predicted populations
#'
#' @export
hAccPerMatch <- function(tax, pred, true, scale = FALSE) {
  pops <- rownames(tax$d)
  n_pops <- length(pops)
  
  d <- -tax$d
  diag(d) <- -diag(d)
  
  res <- rep(NA, n_pops)
  names(res) <- pops
  
  for (pop in pops) {
    idcs <- which(true == pop)
    
    if (scale) {
      lower_bound <- min(d[pop, ])
      upper_bound <- max(d[pop, ])
      unadj <- mean(d[cbind(pred[idcs], true[idcs])])
      res[pop] <- (unadj - lower_bound) / (upper_bound - lower_bound)
    } else {
      res[pop] <- mean(d[cbind(pred[idcs], true[idcs])])
    }
  }
  
  res
}

#' Compute F1 score per each ground-truth label
#'
#' Computes the F1 score (harmonic mean of Precision and Recall) for each cell ground-truth population matched to predicted labels.
#'
#' @param tax object of class \code{Taxonomy} with computed distances
#' @param pred character vector of predicted populations
#' @param true character vector of predicted populations
#'
#' @export
F1PerMatch <- function(tax, pred, true) {
  pops <- rownames(tax$d)
  n_pops <- length(pops)
  
  res <- rep(NA, n_pops)
  names(res) <- pops
  
  for (pop in pops) {
    idcs_pred <- pred == pop
    idcs_true <- true == pop
    if (sum(idcs_pred) == 0 || sum(idcs_true) == 0) {
      res[pop] <- 0
    } else {
      precision <- sum(idcs_pred & idcs_true) / sum(idcs_pred)
      recall <- sum(idcs_pred & idcs_true) / sum(idcs_true)
      f1 <- 2 * (precision * recall) / (precision + recall)
      res[pop] <- f1
    }
  }
  
  res
}

#' Compute Hierarchy-aware F1 score per each ground-truth label
#'
#' Computes the F1 score (harmonic mean of Precision and Recall) for each cell ground-truth population matched to predicted labels, using a hierarchy-aware scoring scheme.
#'
#' @param tax object of class \code{Taxonomy} with computed distances
#' @param pred character vector of predicted populations
#' @param true character vector of predicted populations
#'
#' @export
hF1PerMatch <- function(tax, pred, true) {
  pops <- rownames(tax$d)
  n_pops <- length(pops)
  
  res <- rep(NA, n_pops)
  names(res) <- pops
  
  for (pop in pops) {
    idcs_pred <- pred == pop
    idcs_true <- true == pop
    
    if (sum(idcs_pred) == 0 || sum(idcs_true) == 0) {
      res[pop] <- 0
    } else {
      tp <- which(idcs_pred & idcs_true)
      tp_plus_fp <- which(idcs_pred)
      tp_plus_fn <- which(idcs_true)
      
      numerator <- sum(tax$d[cbind(pred[tp], true[tp])])
      
      precision <- numerator / sum(tax$d[cbind(pred[tp_plus_fp], true[tp_plus_fp])])
      recall <- numerator / sum(tax$d[cbind(pred[tp_plus_fn], true[tp_plus_fn])])
      
      f1 <- 2 * (precision * recall) / (precision + recall)
      res[pop] <- f1
    }
  }
  
  res
}

