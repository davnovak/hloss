
#' Match clusters across two clusterings and plot results
#'
#' Compares two labels of cluster assignment per data point (or a vector of ground-truth labels and a clustering vector) \code{c1} and \code{c2}, matching groups in each vector to each other while maximising the value of an evaluation metric \code{obj}.
#' The evaluation metric \code{obj} is either \code{f1} (default), \code{precision} or \code{recall}.
#'
#' Three approaches are used to solve the cluster-cluster (or label-cluster) matching problem. All of them seek to maximise the total value of \code{obj}.
#' Approach (i) gives 1-to-1 matches, whereby each group in \code{c1} is matched to a (different) group in \code{c2}. (In the special case where the number of groups in \code{c1} is equal to the number of groups in \code{c2}, this guarantees no unmatched groups.)
#' Approach (ii) uses a relaxed fixed-\code{c1} matching, whereby each group in \code{c1} is matched to the group in \code{c2} that maximises \code{obj} value of the match. This can result in 1-to-many matches.
#' Approach (iii) uses a relaxed fixed-\code{c2} matching, which mirrors approach (ii).
#'
#' If \code{c1} is in fact a vector of ground-truth labels (or manual annotation of each data point), there may be de-facto unlabelled data points in the original data.
#' \code{unassigned} is an optional vector of the labels given to data points which don't belong to an annotated population.
#' If specified, the unassigned groups in \code{c1} are left out of the evaluation: points that are unassigned are ignored in constructing the contingency tables for each match and groups in \code{c2} may not be matched to these unassigned points.
#'
#' In addition to evaluation results, a heatmap showing agreement between \code{c1} and \code{c2} and agreement between the different matching approaches is produced by default.
#'
#' @param c1 factor, numeric or character vector: assignment of each data point to a cluster or otherwise defined population
#' @param c2 factor, numeric or character vector: assignment of each data point to a cluster
#' @param obj string: evaluation metric used for matching groups in \code{c1} and \code{c2}; one of \code{f1} (default), \code{precision} and \code{recall}
#' @param title string: tile of Jaccard similarity heatmap plot (default value is '\code{Jaccard heatmap}')
#' @param unassigned optional string vector: names of levels of \code{c1} denoting unlabelled data points
#' @param generate_plot logical: whether a Jaccard heatmap-style plot should be generated (default value is \code{TRUE})
#' @param c1_name optional string: name of the \code{c1} vector to be used in text of the plot (default value is '\code{c1}')
#' @param c2_name optional string: name of the \code{c2} vector to be used in text of the plot (default value is '\code{c2}')
#' @param verbose logical: indicates whether to display progress messages (default value is \code{FALSE})
#'
#' @returns list of results for evaluation approach (i) \code{Results.Bijective}, approach (ii) \code{Results.FixedC1} and approach \code{Results.FixedC2}, as well as a Jaccard similarity heatmap diagram \code{Plot} (if produced)
#'
#' @export
Jacc <- function(
  c1, c2, obj = 'f1', title = 'Jaccard heatmap', unassigned = NULL, generate_plot = TRUE, c1_name = 'c1', c2_name = 'c2', verbose = FALSE
) {
  ## Checks

  stopifnot(length(c1) == length(c2))
  stopifnot(length(c1) > 0)

  if (!is.factor(c1))
    c1 <- as.factor(c1)

  if (!is.factor(c2))
    c2 <- as.factor(c2)

  stopifnot(is.character(obj))
  stopifnot(is.atomic(obj))
  stopifnot(length(obj) == 1)
  stopifnot(tolower(obj) %in% c('f1', 'precision', 'recall'))

  if (is.null(title))
    title <- ''
  stopifnot(is.character(title))
  stopifnot(is.atomic(title))

  if (!is.null(unassigned))
    stopifnot(all(unassigned %in% levels(c1)))

  stopifnot(is.logical(verbose))
  stopifnot(is.atomic(verbose))

  ## Construct Jaccard table

  n_steps <- if (generate_plot) 5 else 4

  if (verbose)
    message('(1/', n_steps, ') Constructing Jaccard similarity table')
  jt <- .JaccardTable(c1, c2)
  counts <- table(c1, c2)

  ## Evaluation: bijective matching and relaxed matching (in both directions)

  if (verbose)
    message(paste0('(2/', n_steps, ') Matching groups bijectively (1-to-1)'))
  eval_bijective <- .BijectiveMatch(c1, c2, obj, unassigned)

  if (verbose)
    message('(3/', n_steps, ') Matching groups for fixed c1 groups (1-to-many)')
  eval_fixedc1   <- .RelaxedMatch(c1, c2, obj, unassigned)

  if (verbose)
    message('(4/', n_steps, ') Matching groups for fixed c2 groups (1-to-many)')
  eval_fixedc2   <- .RelaxedMatch(c2, c1, obj, unassigned, reference_is_c2 = TRUE)

  p <- NULL

  if (generate_plot) {
    ## Create visual representation of matched-group overlaps

    if (verbose)
      message('(5/5) Producing Jaccard similarity heatmap')
    p <- .JaccardSimilarityHeatmap(jt, counts, eval_bijective, eval_fixedc1, eval_fixedc2, unassigned, title, c1_name = c1_name, c2_name = c2_name)
  }

  list(
    Results.Bijective = eval_bijective,
    Results.FixedC1   = eval_fixedc1,
    Results.FixedC2   = eval_fixedc2,
    Plot              = p
  )
}
