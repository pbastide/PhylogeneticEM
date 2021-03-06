% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shifts_manipulations.R
\name{update.incidence.matrix.full}
\alias{update.incidence.matrix.full}
\title{Update function for incidence matrix}
\usage{
update.incidence.matrix.full(daughtersParams, parent, ...)
}
\arguments{
\item{daughtersParams}{: rows of updated matrix corresponding to the
daughters of the current node.}
}
\value{
Vector of length ntaxa + Nnode, indicating to which lineages the 
branch above the current node belongs to.
}
\description{
\code{update.incidence.matrix.full} updates the matrix initialized in 
\code{init.incidence.matrix.full} for the computation of the incidence matrix
in \code{incidence.matrix.full}.
}
\details{
A node belongs to all the lineages of its daughters.
}
\keyword{internal}
