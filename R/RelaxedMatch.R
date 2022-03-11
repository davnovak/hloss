
.RelaxedMatch <- function(
  c1, c2, obj = 'f1', unassigned, reference_is_c2 = FALSE
) {
  ## Create integer mapping of each level of c1 and c2 with NAs for unassigned's

  if (reference_is_c2) {
    levels(c2)[levels(c2) %in% unassigned] <- NA
    c1[is.na(c2)] <- NA
  } else {
    levels(c1)[levels(c1) %in% unassigned] <- NA
    c2[is.na(c1)] <- NA
  }

  l_c1 <- levels(c1)
  l_c2 <- levels(c2)

  c1 <- as.integer(c1)
  c2 <- as.integer(c2)

  ## Store group names and number of groups

  n_c1 <- length(l_c1)
  n_c2 <- length(l_c2)

  ## Iteratively look for optimal match from c2 for each group in c1
  ## Compute Precision, Recall and F1 for each match that is made

  matches <- rep(NA, n_c1)
  names(matches) <- l_c1

  Scores <- matrix(NA, ncol = 3, nrow = n_c1, dimnames = list(l_c1, c('Precision', 'Recall', 'F1')))
  for (idx_c1 in seq_along(l_c1)) {
    stats <-
      t(sapply(seq_along(l_c2), function(idx_c2) {
        TruePositive <- sum(c1 == idx_c1 & c2 == idx_c2, na.rm = TRUE)
        Positive     <- sum(c2 == idx_c2, na.rm = TRUE)
        True         <- sum(c1 == idx_c1, na.rm = TRUE)

        Precision <- if (Positive == 0) 0 else TruePositive / Positive
        Recall    <- if (True == 0) 0 else TruePositive / True

        list(
          precision = Precision,
          recall    = Recall,
          f1        = if (Precision + Recall == 0) 0 else 2 * Precision * Recall / (Precision + Recall)
        )
      }))
    idx_match <- which.max(stats[, obj])[1]
    matches[idx_c1] <- if (stats[idx_match, obj] > 0) l_c2[idx_match] else NA
    Scores[idx_c1, ] <- unlist(stats[idx_match, ])
  }

  ## Get Precision, Recall, F1 and number of matched data points for each matched pair of groups

  Precision.PerMatch <- Scores[, 'Precision']
  Recall.PerMatch    <- Scores[, 'Recall']
  F1.PerMatch        <- Scores[, 'F1']
  NPoints.PerMatch   <- sapply(match(matches, l_c2), FUN = function(idx_c2) sum(c2 == idx_c2, na.rm = TRUE))

  ## Compute mean stats values across matches, unweighted and weighted by size of group in c1 (reference)

  Precision.UnweightedMean <- mean(Precision.PerMatch)
  Recall.UnweightedMean    <- mean(Recall.PerMatch)
  F1.UnweightedMean        <- mean(F1.PerMatch)

  NPoints.c1 <- rep(0, n_c1)
  names(NPoints.c1) <- names(matches)
  tt <- table(c1)
  NPoints.c1[names(tt)] <- tt
  weights <- NPoints.c1 / sum(NPoints.c1)

  Precision.WeightedMean <- sum(Precision.PerMatch * weights)
  Recall.WeightedMean    <- sum(Recall.PerMatch * weights)
  F1.WeightedMean        <- sum(F1.PerMatch * weights)

  ## Return matching scheme, stats and which objective was used for matching

  list(
    Matches                  = matches,
    Objective                = obj,
    Precision.PerMatch       = Precision.PerMatch,
    Recall.PerMatch          = Recall.PerMatch,
    F1.PerMatch              = F1.PerMatch,
    Precision.UnweightedMean = Precision.UnweightedMean,
    Recall.UnweightedMean    = Recall.UnweightedMean,
    F1.UnweightedMean        = F1.UnweightedMean,
    Precision.WeightedMean   = Precision.WeightedMean,
    Recall.WeightedMean      = Recall.WeightedMean,
    F1.WeightedMean          = F1.WeightedMean
  )
}
