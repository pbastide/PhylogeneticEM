% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsimonyNumber.R
\name{transform_shifts_values}
\alias{transform_shifts_values}
\title{Transform the shift values}
\usage{
transform_shifts_values(shifts, from = 0, to, phylo)
}
\arguments{
\item{shifts}{the shifts on the original process}

\item{from}{alpha value of the original process. If equals 0, then the original process is 
taken to be a BM.}

\item{to}{alpha value of the destination process}

\item{phylo}{the phylogenetic tree (un-scaled)}
}
\description{
\code{transform_shifts_values} takes the shifts generating a given expectation structure
given an OU with alpha = from, and gives back the equivalent shifts values that produce the
same structure with an OU with alpha = to. If from or to is 0, then the process is supposed
to be a BM.
}
\keyword{internal}
