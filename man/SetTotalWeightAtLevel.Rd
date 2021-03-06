% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Taxonomy.R
\name{SetTotalWeightAtLevel}
\alias{SetTotalWeightAtLevel}
\title{Set total distance across branches at particular depth of taxonomic tree}
\usage{
SetTotalWeightAtLevel(tax, depth, value, round_digits = 3)
}
\arguments{
\item{tax}{object of class \code{Taxonomy}}

\item{depth}{integer: desired depth (root is at 1)}

\item{value}{numeric: desired sum of weights of branches at the desired depth}

\item{round_digits}{numeric: to how many decimal places should the uniform weight be rounded (or \code{NULL} to not apply \code{round}). Defaults to 3}
}
\description{
Branches are weighted uniformly, with a total value for a certain level of depth.
For a particular \code{depth}, branches leading from \code{depth} to \code{depth+1} are changed.
}
