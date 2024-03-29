% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_EM.R
\name{compute_gauss_lasso}
\alias{compute_gauss_lasso}
\title{Do a lm on top of a lasso regression.}
\usage{
compute_gauss_lasso(
  Ypt,
  Xp,
  delta,
  root,
  projection = which(rowSums(delta) != 0)
)
}
\arguments{
\item{Ypt}{(transformed) data}

\item{Xp}{(transformed) matrix of regression}

\item{delta}{regression coefficients obtained with a lasso regression}

\item{root}{the position of the root (intercept) in delta}
}
\value{
Named list, with "E0.gauss" the intercept (value at the root);
"shifts.gauss" the list of shifts found on the branches; and "residuals" the
residuals of the regression
}
\description{
\code{compute_gauss_lasso} takes the variables selected by a lasso procedure, and
uses them to do a simple linear least square regression. Function used is
\code{lm} for non-transformed data (root = NULL), and \code{lm.fit} for
transformed data (root = an integer).
}
\details{
Depending on the value of root, the behavior is different. If root is null, then
we fit a linear regression with an intercept. If root is equal to an integer,
then the "intercept" column of the matrix Xp (that has possibly been trough a 
multiplication with a Cholesky decomposition of the variance) is included, rather
than the intercept.
}
\keyword{internal}
