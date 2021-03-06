% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsimonyNumber.R
\name{extract.equivalent_shifts}
\alias{extract.equivalent_shifts}
\title{Extract the shifts values for one trait.}
\usage{
\method{extract}{equivalent_shifts}(x, trait = 1, what = c("shifts_values", "root_values"), ...)
}
\arguments{
\item{x}{an object of class \code{equivalent_shifts}, result of
function \code{\link{equivalent_shifts}}}

\item{trait}{the number of the trait to be extracted. Default to 1.}

\item{what}{one of "shifts_values" or "root_values".}

\item{...}{unused.}
}
\value{
A matrix with the values of the shifts (\code{what = "shifts_values"}) or
the root (\code{what = "root_values"}) for the trait for each equivalent
configuration. Each column is one configuration.
}
\description{
\code{extract.equivalent_shifts} takes an object of class
\code{equivalent_shifts}, result of function \code{\link{equivalent_shifts}},
and returns the shifts of root values for a given trait.
}
\seealso{
\code{\link{equivalent_shifts}}, \code{\link{plot.equivalent_shifts}},
\code{\link{equivalent_shifts_edges}}
}
