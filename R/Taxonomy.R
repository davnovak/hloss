
Node <- function(name, parent = c(), count = NA) {
  x          <- new.env(hash = TRUE)
  x$name     <- name
  x$parent   <- parent
  x$children <- c()
  x$count    <- count
  structure(x, class = 'Node')
}

print.Node <- function(x) {
  if (length(x$parent) == 0) {
    message('Root node ', x$name, '')
  } else {
    message('Node ', x$name, '')
    message('-> parent node:\n\t', x$parent)
  }
  message('-> child nodes:\n\t', paste(x$children, collapse = ',\n\t'))
  message('-> cell count:\n\t', x$count)
}

#' Instantiate a taxonomic tree
#'
#' Constructor for a \code{Taxonomy} object.
#' Call \code{AddNode} on a \code{Taxonomy} object to build it up.
#' Alternatively, use function \code{TaxonomyFromPaths} to construct the tree in one call.
#'
#' To change length of branches, use \code{SetDistanceToChild} or \code{SetDistancesFromParent}.
#'
#' Use \code{print} to print a read-out of paths in the tree and its properties.
#' Use \code{plot} to plot the tree structure.
#'
#' @param name string: name of taxonomy
#'
#' @export
Taxonomy <- function(name) {
  x <- new.env(hash = TRUE)
  x$n                 <- 0
  x$depth             <- 0
  x$node_depths       <- c()
  x$name              <- name
  x$node_names        <- c()
  x$node_paths        <- c()
  x$parent            <- c(0)
  x$is_leaf           <- c()
  x$distance          <- c() # to parent
  x$count             <- c()
  x$nodes             <- list()
  x$d                 <- NULL # distance matrix
  x$df                <- NULL # flipped-distance matrix
  x$dmax              <- NULL # max distance
  x$changes_made      <- TRUE
  x$use_cell_counts   <- FALSE
  structure(x, class = 'Taxonomy')
}

print.Taxonomy <- function(tax) {
  message('Taxonomy "', tax$name, '"')
  message('-> number of leaves: ', sum(tax$is_leaf))
  message('-> number of nodes: ', tax$n)
  message('-> depth: ', tax$depth)
  if (!is.na(tax$count[1])) {
    message('-> total cell count: ', tax$count[1])
  }
  message('• ', paste(sort(tax$node_paths), collapse = '\n• '))
}

plot.Taxonomy <- function(tax, show_weights = TRUE, cex_vertex = 0.9, cex_edge = 1.1) {
  names <- tax$node_names
  parent <- tax$parent
  parent[1] <- 1
  parent <- names[parent]
  
  g <- igraph::graph_from_data_frame(data.frame('parent' = parent, 'name' = tax$node_names))
  g <- g - igraph::edge(paste0(tax$node_names[1], '|', tax$node_names[1]))
  
  igraph::V(g)$color <- c('#758deb', rep(NA, times = length(tax$nodes) - 1))
  igraph::V(g)$label.color <- 'black'
  igraph::V(g)$label.cex <- cex_vertex
  
  igraph::V(g)$frame.color <- 'white'
  
  if (show_weights) {
    igraph::E(g)$color <- '#ebebeb'
    igraph::E(g)$label.cex <- cex_edge
  }
  
  w <- tax$distance[2:tax$n]
  w <- w - min(w) + 1
  w <- w * 4 - 0.5
  
  old_par <- par(no.readonly = TRUE)
  par(mar = c(0, 0, 0, 0))
  igraph::plot.igraph(
    x = g,
    edge.width = if (show_weights) w else 1,
    edge.label = if (show_weights) tax$distance[2:tax$n] else c()
  )
  dummy <- gc(verbose = FALSE)
  do.call(par, old_par)
}

GetPathFromName <- function(tax, name) {
  path <- tax$node_paths[tax$node_names == name]
  if (length(path) == 0) {
    stop('Could not find path for node name "', name, '" in taxonomy "', tax$name, '"')
  } else if (length(path) > 1) {
    stop('Node named "', name, '" in taxonomy "', tax$name, '" can be matched to multiple paths:\n\t', paste(path, collapse = ',\n\t'))
  }
  path
}

GetPathFromSubstring <- function(tax, s) {
  if (grepl('/', s, fixed = TRUE)) {
    exact_match <- tax$node_paths[tax$node_paths == s]
    if (length(exact_match) == 0) {
      path <- tax$node_paths[grep(s, tax$node_paths, fixed = TRUE)]
      path <- path[
        grepl(paste0(s, '$'), path)
      ]
      if (length(path) == 0) {
        stop('Could not find path for substring "', s, '" in taxonomy "', tax$name, '"')
      } else if (length(path) > 1) {
        stop('Node substring "', s, '" in taxonomy "', tax$name, '" can be matched to multiple paths:\n\t', paste(path, collapse = ',\n\t'))
      }
    } else {
     path <- exact_match 
    }
  } else {
    path <- GetPathFromName(tax, s)
  }
  path
}

GetPathsFromSubstrings <- function(tax, s) {
  s <- as.factor(s)
  levels(s) <- sapply(levels(s), function(x) GetPathFromSubstring(tax, x))
  s
}

CheckIfCountValid <- function(tax, name, count, idx_parent) {
  idx_siblings   <- which(tax$parent == idx_parent)
  count_parent   <- tax$count[idx_parent]
  if (!is.na(count_parent)) {
    if (count_parent < count) {
      stop('Child node ', name, ' may not have a higher cell count (', count, ') than its parent ', tax$node_paths[idx_parent], ' (', count_parent, ')')
    }
    if (length(idx_siblings) > 0) {
      counts_siblings <- tax$count[idx_siblings]
      if (any(!is.na(counts_siblings))) {
        if (sum(counts_siblings) + count > count_parent) {
          stop('Child nodes of the node ', tax$node_paths[idx_parent], ' may not have a higher summed cell count than the parent node itself')
        }
      }
    }
  }
}

#' Add node to a taxonomic tree
#'
#' Adds a node to a \code{Taxonomy} object if possible.
#'
#' To change length of a specific branch, use \code{SetDistanceToChild}.
#' To change length of branches leading from a parent node, use \code{SetDistancesFromParent}.
#'
#' @param tax object of class \code{Taxonomy}
#' @param name string: name of new node
#' @param parent string: name of parent node (if non-root)
#' @param distance numeric: length of branch leading from node to this new node (if non-root). Defaults to \code{1}
#' @param count integer: count of cells in dataset associated with the added node
#'
#' @export
AddNode <- function(tax, name, parent = c(), distance = 1, count = NA) UseMethod('AddNode', tax)

AddNode.Taxonomy <- function(tax, name, parent = c(), distance = 1, count = NA) {
  if (length(parent) == 0 && tax$n > 0) {
    stop('Taxonomy "', tax$name, '" has already got a root node "', tax$members[[1]]$name, '"')
  }
  if (!length(parent) %in% c(0, 1)) {
    stop('Zero or one parent must be specified for any node')
  }
  
  is_root <- length(parent) == 0
  
  tax$node_names <- c(tax$node_names, name)
  tax$nodes      <- c(tax$nodes, Node(name = name, parent = parent))
  tax$n          <- tax$n + 1
  tax$is_leaf    <- c(tax$is_leaf, TRUE)
  tax$count      <- c(tax$count, count)
  
  depth <- 1
  if (is_root) {
    tax$node_paths <- c(tax$node_paths, name)
    tax$distance   <- c(tax$distance, 0)
    tax$node_depths <- c(tax$node_depths, depth)
  } else {
    parent         <- GetPathFromSubstring(tax, parent)
    idx_parent     <- which(tax$node_paths == parent)
    
    CheckIfCountValid(tax, name, count, idx_parent)
    
    tax$parent     <- c(tax$parent, idx_parent)
    path <- paste0(tax$node_paths[idx_parent], '/', name)
    tax$node_paths <- c(tax$node_paths, path) 
    tax$distance   <- c(tax$distance, distance)
    tax$nodes[[idx_parent]]$children <- c(tax$nodes[[idx_parent]]$children, path)
    depth          <- tax$node_depths[idx_parent] + 1
    tax$node_depths<- c(tax$node_depths, depth)
    tax$is_leaf[idx_parent] <- FALSE
  }
  tax$depth <- max(tax$node_depths)
  
  dup <- unique(tax$node_paths[duplicated(tax$node_paths)])
  if (length(dup) > 0) {
    stop('Node already exists: ', dup)
  }
  
  tax$changes_made <- TRUE
  
  invisible(tax)
}

SetCount <- function(tax, name, count) UseMethod('SetCount', tax)

SetCount.Taxonomy <- function(tax, name, count) {
  path <- GetPathFromSubstring(tax, name)
  idx_node <- which(tax$node_paths == path)
  idx_parent <- tax$parent[idx_node]
  
  CheckIfCountValid(tax, path, count, idx_parent)
  
  tax$count[idx_node] <- count
  
  invisible(tax)
}

#' Set length of branches leading from a parent node
#'
#' @param tax object of class \code{Taxonomy}
#' @param name string: path (or uniquely matchable substring thereof) of a parent node
#' @param distance numeric: branch length
#'
#' @export
SetDistancesFromParent <- function(tax, name, distance) UseMethod('SetDistancesFromParent', tax)

SetDistancesFromParent.Taxonomy <- function(tax, name, distance) {
  path <- GetPathFromSubstring(tax, name)
  idx_node <- which(tax$node_paths == path)
  idx_children <- which(tax$parent == idx_node)
  tax$distance[idx_children] <- distance
  invisible(tax)
}

#' Set total distance across branches at particular depth of taxonomic tree
#'
#' Branches are weighted uniformly, with a total value for a certain level of depth.
#' For a particular \code{depth}, branches leading from \code{depth} to \code{depth+1} are changed.
#'
#' @param tax object of class \code{Taxonomy}
#' @param depth integer: desired depth (root is at 1)
#' @param value numeric: desired sum of weights of branches at the desired depth
#' @param round_digits numeric: to how many decimal places should the uniform weight be rounded (or \code{NULL} to not apply \code{round}). Defaults to 3
#'
#' @export
SetTotalWeightAtLevel <- function(tax, depth, value, round_digits = 3) UseMethod('SetTotalWeightAtLevel', tax)

SetTotalWeightAtLevel.Taxonomy <- function(tax, depth, value, round_digits = 3) {
  idcs <- which(tax$node_depths == depth + 1)
  n <- length(idcs)
  value <- value / n
  if (!is.null(round_digits)) {
    value <- round(value, round_digits)
  }
  tax$distance[idcs] <- value
  invisible(tax)
}


#' Set length of a branch
#'
#' @param tax object of class \code{Taxonomy}
#' @param name string: path (or uniquely matchable substring thereof) of a node
#' @param distance numeric: branch length
#'
#' @export
SetDistanceToChild <- function(tax, name, distance) UseMethod('SetDistanceToChild', tax)

SetDistanceToChild.Taxonomy <- function(tax, name, distance) {
  path <- GetPathFromSubstring(tax, name)
  tax$distance[tax$node_paths == path] <- distance
  invisible(tax)
}

#' Build a taxonomic tree from node paths
#'
#' Constructs a \code{Taxonomy} object with a hierarchy of labels.
#' Requires a vector of paths to each class, eg. \code{c'Lymphocyte', 'Lymphocyte/T Cell', 'Lymphocyte/B Cell', 'Lymphocyte/T Cell/CD4+', 'Lymphocyte/T Cell/CD8+'}.
#' Assumes default branch lengths (all \code{1}), which can then be changed using \code{SetDistanceToChild} or \code{SetDistancesFromParent}.
#'
#' Use \code{print} on a \code{Taxonomy} object to print a read-out of paths in the tree and its properties.
#' Use \code{plot} on a \code{Taxonomy} object to plot the tree structure.
#'
#' @param name string: name of taxonomy
#' @param paths string vector: path to each node to be added
#'
#' @export
TaxonomyFromPaths <- function(name, paths) {
  p   <- strsplit(paths, '/')
  p   <- p[order(sapply(p, length))]
  tax <- Taxonomy(name = name)
  AddNode(tax, name = p[[1]][1])
  for (idx in 2:length(p)) {
    l <- length(p[[idx]])
    AddNode(tax, name = p[[idx]][l], parent = paste0(p[[idx]][1:(l-1)], collapse = '/'))
  }
  tax
}

#' Build a calibrated taxonomic tree
#'
#' Constructs a \code{Taxonomy} object with a hierarchy of labels.
#' Requires a vector of paths to each class, eg. \code{c'Lymphocyte', 'Lymphocyte/T Cell', 'Lymphocyte/B Cell', 'Lymphocyte/T Cell/CD4+', 'Lymphocyte/T Cell/CD8+'}.
#'
#' Counts associated with each node can also be given, in which case weights of branches leading to each node are set as inversely proportional to the number of cells associated with it.
#' (This is proposed for large datasets where population sizes are considered representative and indicate how rare or common each population typically is.)
#'
#' Otherwise, we set weights such that the sum of weights of branches coming out of each parent node is \code{1} and the weights are distributed uniformly among the branches.
#'
#' Use \code{print} on a \code{Taxonomy} object to print a read-out of paths in the tree and its properties.
#' Use \code{plot} on a \code{Taxonomy} object to plot the tree structure.
#'
#' @param name string: name of taxonomy
#' @param paths string vector: path to each node to be added
#' @param counts integer vector or \code{NULL}: cell count associated with each node to be added. Defaults to \code{NULL}
#' @param round_digits integer: how many decimal places to round weight values to, or \code{NULL}. Defaults to \code{3}
#'
#' @export
CalibratedTaxonomy <- function(name, paths, counts = NULL, round_digits = 3) {
  p <- strsplit(paths, '/')
  o <- order(sapply(p, length))
  p <- p[o]
  use_cell_counts <- !is.null(counts)
  if (use_cell_counts) {
    counts <- counts[o]
  }
  tax <- Taxonomy(name = name)
  tax$use_cell_counts <- use_cell_counts
  AddNode(tax, name = p[[1]][1], count = if (use_cell_counts) counts[1] else NA)
  for (idx in 2:length(p)) {
    l <- length(p[[idx]])
    AddNode(tax, name = p[[idx]][l], parent = paste0(p[[idx]][1:(l-1)], collapse = '/'), count = if (use_cell_counts) counts[idx] else NA)
  }
  Calibrate(tax, across = 'subtree', from = 1, to = 1, f = NULL, use_cell_counts = use_cell_counts, round_digits = round_digits)
  tax
}

FormatPath <- function(path) {
  strsplit(path ,'/')[[1]]
}

GetFormattedPaths <- function(tax) {
  strsplit(tax$node_paths, '/')
}

#' Get step count & distance between two nodes in a taxonomic tree
#'
#' Finds the shortest path between two nodes in a \code{Taxonomy} object.
#'
#' @param tax \code{Taxonomy} object
#' @param from string: uniquely matchable substring identifying start node
#' @param to string: uniquely matchable substring identifying end node
#' @param only_total logical: whether to only return a single number (total distance), as opposed to rootward and leafward distances. Defaults to \code{FALSE}
#' @param weighted logical: whether to include weighted distances as well (or only the weighted distance, if \code{only_total} is \code{TRUE}. Defaults to \code{TRUE}
#' 
#'
#' @return If \code{only_total} is \code{TRUE} a single number (total distance, weighted if \code{weighted} is \code{TRUE}) is returned, otherwise a vector of \code{rootward} unweighted distance and \code{leafward} unweighted distance is returned (if \code{weighted} is \code{TRUE}, then weighted distances \code{rootward_w} and \code{leafward_w} are also returned).
#'
#' @export
Path <- function(tax, from, to, only_total = FALSE, weighted = TRUE, exp_decay = FALSE) {
  from <- FormatPath(GetPathFromSubstring(tax, from))
  to <- FormatPath(GetPathFromSubstring(tax, to))
  
  suppressWarnings(overlap <- sum(from == to))
  # mrca <- paste(from[1:overlap], collapse = '/') # most recent common ancestor
  
  n_from <- length(from)
  n_to <- length(to)
  
  rootward <- n_from - overlap
  leafward <- n_to - overlap
  
  if (weighted) {
    if (rootward == 0) {
      rootward_d <- c()
      
      # rootward_w <- 0
    } else {
      waypoints_rootward <- sapply((overlap+1):n_from, function(n) paste(from[1:n], collapse = '/'))
      idcs <- match(waypoints_rootward, tax$node_paths)
      
      rootward_d <- tax$distance[idcs]
      # w <- dexp(rootward_d, rate = mean(rootward_d))
      # rootward_w <- sum(rootward_d * w) * multiplier_rootward
      
      # rootward_w <- sum(tax$distance[idcs]) * multiplier_rootward
    }
    if (leafward == 0) {
      leafward_d <- c()
      
      # leafward_d <- 0
    } else {
      waypoints_leafward <- sapply((overlap+1):n_to, function(n) paste(to[1:n], collapse = '/'))
      idcs <- match(waypoints_leafward, tax$node_paths)
      
      leafward_d <- tax$distance[idcs]
      # w <- dexp(leafward_d, rate = mean(leafward_d))
      # leafward_w <- sum(leafward_d * w) * multiplier_leafward
      
      # leafward_w <- sum(tax$distance[idcs]) * multiplier_leafward
    }
  }
  
  if (only_total) {
    if (weighted) {
      
      distances <- c(rootward_d, leafward_d)
      weights <- if (exp_decay) dexp(distances, rate = mean(distances)) else 1
      sum(distances * weights)
      
      # rootward_w + leafward_w
    } else {
      rootward + leafward
    }
  } else {
    if (weighted) {
      
      distances <- c(rootward_d, leafward_d)
      weights <- if (exp_decay) dexp(distances, rate = mean(distances)) else 1
      rootward_w <- if (rootward == 0) 0 else distances[1:rootward] * weights[1:rootward]
      leafward_w <- if (leafward == 0) 0 else distances[(rootward+1):(rootward+leafward)] * weights[(rootward+1):(rootward+leafward)]
      
      c('rootward' = rootward, 'leafward' = leafward, 'rootward_w' = rootward_w, 'leafward_w' = leafward_w)
    } else {
      c('rootward' = rootward, 'leafward' = leafward)
    }
  }
}

#' Get step count & distance from node to root
#'
#' Computes number of steps and total distance from node to root in a \code{Taxonomy} object.
#'
#' @param tax \code{Taxonomy} object
#' @param node string: uniquely matchable substring identifying start node
#' @param single logical: whether to only return a single number (total distance), as opposed to rootward and leafward distances. Defaults to \code{FALSE}
#' @param weighted logical: whether to include weighted distance (or only the weighted distance, if \code{single} is \code{TRUE}. Defaults to \code{TRUE}
#' 
#' @return If \code{single} is \code{TRUE} a single number (total distance, weighted if \code{weighted} is \code{TRUE}) is returned, otherwise a vector of \code{unweighted} and \code{weighted} distance is returned.
#'
#' @export
PathToRoot <- function(tax, node, single = FALSE, weighted = TRUE) {
  node <- FormatPath(GetPathFromSubstring(tax, node))
  
  unweighted <- length(node) - 1
  
  if (single && !weighted) {
    return(unweighted)
  }
  
  if (unweighted > 0) {
    waypoints <- sapply((1:unweighted)+1, function(n) paste(node[1:n], collapse = '/'))
    weighted <- sum(tax$distance[match(waypoints, tax$node_paths)])
  } else {
    weigthed <- 0
  }
  
  if (single) {
    return(weighted)
  }
  
  c('unweighted' = unweighted, 'weighted' = weighted)
}

CalibrateBySubtreeRecursively <- function(tax, s, idx_parent, use_cell_counts, round_digits) {
  idcs_children <- match(tax$nodes[[idx_parent]]$children, tax$node_paths)
  n_children <- length(idcs_children)
  if (n_children > 0) { # stop condition
    if (use_cell_counts && all(!is.na(tax$count[idcs_children]))) {
      vals <- 1 / tax$count[idcs_children]
      vals <- vals / sum(vals) * s[tax$node_depths[idx_parent]]
      if (!is.null(round_digits)) {
        vals <- round(vals, round_digits)
      }
      tax$distance[idcs_children] <- vals
    } else {
      val <- s[tax$node_depths[idx_parent]] / n_children
      if (!is.null(round_digits)) {
        val <- round(val, round_digits)
      }
      tax$distance[idcs_children] <- val
    }
    for (idx_child in idcs_children) {
      CalibrateBySubtreeRecursively(tax, s, idx_child, use_cell_counts, round_digits)
    }
  }
  invisible(tax)
}

#' Calibrate \code{Taxonomy} branch weights automatically
#' 
#' Calibrates branch weights of a \code{Taxonomy} automatically based on the tree topology.
#' 
#' We either fix the sum of branch weights at each depth level or fix the sum of branch weights coming out of each parent node.
#' For the first option, set parameter \code{across} to \code{level}.
#' For the second option, set parameter \code{across} to \code{subtree}.
#' 
#' You can use different values for the sum of each depth level (or for subtrees at different depth levels).
#' To do this, use the parameters \code{from} and \code{to} to create a series of values through linear interpolation.
#' To transform the values using a function (for example \code{exp}), pass the function as parameter \code{f}.
#' 
#' If calibrating by subtree, branch weights can be distributed among child nodes uniformly, or as inverserly proportional to counts
#' of cells associated with each child node. This way, misidentifying a rare population 
#'
#' @param tax \code{Taxonomy} object
#' @param across character: level of calibration (either \code{level} or \code{subtree}). Defaults to \code{subtree}
#' @param from numeric: root-level summed weight. Defaults to \code{1}
#' @param to numeric: leaf-level summed weight. Defaults to \code{1}
#' @param f function: function to transform the sequence \code{[from, to]} of length \code{tax$depth} or \code{NULL} to keep the sequence. Defaults to \code{NULL}
#' @param use_cell_counts logical: if \code{across == 'subtree'}, whether to divide branch weights within child nodes based on known cell counts associated with them (inverse proportionality). Alternatively, weights are uniform. Defaults to \code{TRUE}
#' @param round_digits numeric: to how many decimal places should the uniform weight be rounded (or \code{NULL} to not apply \code{round}). Defaults to \code{3}
#'
#' @export
Calibrate <- function(tax, across = 'subtree', from = 1, to = 1, f = NULL, use_cell_counts = TRUE, round_digits = 3) UseMethod('Calibrate', tax)

Calibrate.Taxonomy <- function(tax, across = 'subtree', from = 1, to = 1, f = NULL, use_cell_counts = TRUE, round_digits = 3) {
  across <- match.arg(across, choices = c('level', 'subtree'))
  s <- seq(from = from, to = to, length.out = tax$depth)
  if (!is.null(f)) {
    s <- f(s)
  }
  if (across == 'level') {
    for (idx in seq_along(s)) {
      SetTotalWeightAtLevel(tax, depth = idx, value = s[idx], round_digits = round_digits)
    }
  } else if (across == 'subtree') {
    CalibrateBySubtreeRecursively(tax, s, idx_parent = 1, use_cell_counts, round_digits)
  }
  invisible(tax)
}

GetChildren <- function(tax, idx) {
  idcs <- which(tax$parent == idx)
  if (length(idcs) > 0) {
    c(unlist(sapply(idcs, function(x) GetChildren(tax, x))), idcs)
  } else {
    c()
  }
}

#' Generate a named vector of sampling probabilites
#' 
#' Uses a `Taxonomy` object with known counts of cells associated with nodes in the taxonomy and produces a vector of probabilites for sampling a cell by node of taxonomy.
#' This is to generate a random sampling as baseline to be compared with clustering.
#'
#' @param tax \code{Taxonomy} object
#'
#' @export
SamplingProbabilities <- function(tax) {
  if (!tax$use_cell_counts) {
    rep(1, tax$n)
  }
  probs <- tax$count
  for (idx in 1:length(probs)) {
    probs[idx] <- probs[idx] - sum(tax$count[which(tax$parent == idx)])
  }
  probs / sum(probs)
}

