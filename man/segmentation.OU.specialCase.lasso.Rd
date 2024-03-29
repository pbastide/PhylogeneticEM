% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/M_step.R
\name{segmentation.OU.specialCase.lasso}
\alias{segmentation.OU.specialCase.lasso}
\title{Segmentation in the OU special case, using lasso regression}
\usage{
segmentation.OU.specialCase.lasso(
  phylo,
  nbr_of_shifts,
  D,
  Xp,
  penscale = rep(1, (nrow(phylo$edge) + 1)),
  ...
)
}
\arguments{
\item{phylo}{a phylogenetic tree}

\item{nbr_of_shifts}{Number of shifts on the phylogeny allowed}
}
\value{
List containing : beta_0 : the optimal value at the root
                          shifts : list containing the computed tau and delta
                          costs : vector of costs
}
\description{
\code{segmentation.OU.specialCase.lasso} performs the segmentation using a 
lasso regression to select for the edges where the shifts are added.
}
\details{
This function re-write the sum of costs to be minimized as a least squares 
 regression problem, and uses a lasso regression to solve it. It uses
 functions \code{incidence.matrix.full} to express the problem as a 
 linear model.
}
\keyword{internal}
