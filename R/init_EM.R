# {initializations of the EM}
# Copyright (C) {2014} {SR, MM, PB}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
###############################################################################
## Here are the functions to initialize the EM.
## Dependencies : generic_functions.R
## : simulate.R
## : shifts_manipulations.R
###############################################################################

######################################################
## Defaults initializations
######################################################
##
# init.EM.default (process, ...)
# PARAMETERS:
#            @process (string) Random process to simulate. (See note above)
# RETURNS:
#            (list) list of initial parameters for the model
# DEPENDENCIES:
#            init.EM.default.BM, init.EM.default.OU
# PURPOSE:
#            Make a simple default initialization of the parameters of the model
# NOTES:
#            Useless raw, but can be re-used later as a first step (to define a structure)
# REVISIONS:
#            22/05/14 - Initial release
#            28/05/14 - Add OU
##
init.EM.default <- function(process){
  if (process == "BM"){
    return(init.EM.default.BM)
  } else if (process == "OU"){
    return(init.EM.default.OU)
  }
}

init.EM.default.BM <- function(phylo = NULL,
                               Y_data = matrix(NA, 1, length(phylo$tip.label)),
                               p = nrow(Y_data),
                               variance.init = diag(1, p, p),
                               random.init = FALSE,
                               value.root.init = rep(0, p),
                               exp.root.init = rep(1, p),
                               var.root.init = diag(1, p, p),
                               edges.init = NULL,
                               values.init = matrix(0, p, length(edges.init)),
                               relativeTimes.init = NULL,
                               nbr_of_shifts = length(edges.init),
                               subtree.list = NULL,
                               sBM_variance = FALSE, ...) {
  if (random.init) {
    value.root.init <- NA
    if (sBM_variance){
      var.root.init <- phylo$root.edge * variance.init
    }
  } else {
    exp.root.init <- NA
    var.root.init <- NA
  }
  # Always start with some shifts, in case of default initialization (if number of shifts different from 0)
  if (length(edges.init) < nbr_of_shifts){
    miss <- nbr_of_shifts - length(edges.init)
    auth_edges <- which(!(1:nrow(phylo$edge) %in% edges.init))
    new_edges <- sample(auth_edges, miss)
    edges.init <- c(edges.init, new_edges)
    # edges.init <- c(edges.init, sample_shifts_edges(phylo, miss, part.list = subtree.list))
  }
  # If not enought values, complete with 0s
  if (is.null(values.init) || is.vector(values.init)){
    n_shifts_provided <- length(values.init)
  } else {
    n_shifts_provided <- ncol(values.init)
  }
  if (n_shifts_provided < nbr_of_shifts){
    miss <- nbr_of_shifts - ncol(values.init)
    values.init <- cbind(values.init, matrix(0, ncol = miss, nrow = p))
  }
  params_init = list(variance = variance.init,
                     root.state = list(random = random.init,
                                       value.root = value.root.init,
                                       exp.root = exp.root.init,
                                       var.root = var.root.init),
                     shifts = list(edges = edges.init,
                                   values = values.init,
                                   relativeTimes = relativeTimes.init))
  params_init <- check_dimensions(p,
                                  params_init$root.state,
                                  params_init$shifts,
                                  params_init$variance)
  params_init$root.state <- test.root.state(params_init$root.state, "BM")
  params_init$variance <- as(params_init$variance, "dpoMatrix")
  return(params_init)
}

init.EM.default.OU <- function(phylo = NULL,
                               Y_data = matrix(NA, 1, length(phylo$tip.label)),
                               p = nrow(Y_data),
                               variance.init = diag(1, p, p),
                               random.init = TRUE,
                               stationary.root.init = TRUE,
                               value.root.init = rep(1, p),
                               exp.root.init = rep(1, p),
                               var.root.init = diag(1, p, p),
                               edges.init = NULL,
                               values.init = matrix(0, p, length(edges.init)),
                               relativeTimes.init = NULL,
                               selection.strength.init=1,
                               optimal.value.init = rep(0, p),
                               nbr_of_shifts = length(edges.init),
                               subtree.list = NULL, ...) {
  if (random.init) {
    value.root.init <- NA
    if (stationary.root.init) {
      exp.root.init <- optimal.value.init
      variance.init <- var.root.init * (2 * selection.strength.init)
    }
  } else {
    exp.root=NA
    var.root=NA
  }
  # Always start with some shifts, in case of default initialization (if number of shifts different from 0)
  if (length(edges.init) < nbr_of_shifts){
    miss <- nbr_of_shifts - length(edges.init)
    auth_edges <- which(!(1:nrow(phylo$edge) %in% edges.init))
    new_edges <- sample(auth_edges, miss)
    edges.init <- c(edges.init, new_edges)
    # edges.init <- c(edges.init, sample_shifts_edges(phylo, miss, part.list = subtree.list))
  }
  # If not enought values, complete with 0s
  if (is.null(values.init) || is.vector(values.init)){
    n_shifts_provided <- length(values.init)
  } else {
    n_shifts_provided <- ncol(values.init)
  }
  if (n_shifts_provided < nbr_of_shifts){
    miss <- nbr_of_shifts - ncol(values.init)
    values.init <- cbind(values.init, matrix(0, ncol = miss, nrow = p))
  }
  params_init <- list(variance = variance.init,
                      root.state = list(random = random.init,
                                        stationary.root = stationary.root.init,
                                        value.root = value.root.init,
                                        exp.root = exp.root.init,
                                        var.root = var.root.init),
                      shifts = list(edges = edges.init,
                                    values = values.init,
                                    relativeTimes = relativeTimes.init),
                      selection.strength = selection.strength.init,
                      optimal.value = optimal.value.init)
  params_init <- check_dimensions(p,
                                  params_init$root.state,
                                  params_init$shifts,
                                  params_init$variance,
                                  params_init$selection.strength,
                                  params_init$optimal.value)
  params_init$root.state <- test.root.state(params_init$root.state, "OU",
                                            variance = variance.init,
                                            selection.strength = params_init$selection.strength,
                                            optimal.value = optimal.value.init)
  params_init$variance <- as(params_init$variance, "dpoMatrix")
  return(params_init)
}

##
#' @title Create an object params_process
#'
#' @description
#' \code{params_process} creates or extracts a set of parameters of class
#' \code{params_process}.
#'
#' @param x an S3 object.
#' @param ... further arguments to be passed to the specific method.
#' 
#' @return An S3 object of class \code{params_process}. This is essentially a list containing the following entries:
#' \describe{
#' \item{process}{The model used. One of "BM" (for a full BM
#' model, univariate or multivariate); "OU" (for a full OU model, univariate or
#' multivariate); or "scOU" (for a "scalar OU" model).}
#' \item{p}{Dimension of the trait.}
#' \item{root.state}{List describing the state of the root, with:
#' \describe{
#'  \item{random}{random state (TRUE) or deterministic state (FALSE)}
#'  \item{value.root}{if deterministic, value of the character at the root}
#'  \item{exp.root}{if random, expectation of the character at the root}
#'  \item{var.root}{if random, variance of the character at the root (pxp matrix)}
#'  }}
#' \item{shifts}{List with position and values of the shifts:
#' \describe{
#'  \item{edges}{vector of the K id of edges where the shifts are}
#'  \item{values}{matrix p x K of values of the shifts on the edges
#'   (one column = one shift)}
#'  \item{relativeTimes}{vector of dimension K of relative time of the shift from the parent node of edges}
#'  }}
#' \item{variance}{Variance-covariance matrix size p x p.}
#' \item{selection.strength}{Matrix of selection strength size p x p (OU).}
#' \item{optimal.value}{Vector of p optimal values at the root (OU).}
#' }
#' 
#' @seealso \code{\link{params_process.character}},
#' \code{\link{params_process.PhyloEM}},
#' \code{\link{params_BM}}, \code{\link{params_OU}}
#' \code{\link{simul_process.params_process}}
#' 
#' @export
##
params_process <- function(x, ...) UseMethod("params_process")

##
#' @export
#' @method print params_process
##
print.params_process <- function(x, ...){
  if (x$process == "BM"){
    cat(paste0("\n", ncol(x$variance), " dimensional BM process with a "))
    if (x$root.state$random){
      cat("random root.\n")
      cat("Root expectations:\n")
      print(x$root.state$exp.root)
      cat("\n")
      cat("Root variance:\n")
      print(as.matrix(x$root.state$var.root))
      cat("\n")
    } else {
      cat("fixed root.\n")
      cat("\nRoot value:\n")
      print(x$root.state$value.root)
      cat("\n")
    }
    cat("Process variance:\n")
    print(as.matrix(x$variance))
    cat("\n")
  } else {
    if (x$process == "OU"){
      cat(paste0("\n", ncol(x$variance), " dimensional OU process with a "))
    } else if (x$process == "scOU"){
      cat(paste0("\n", ncol(x$variance), " dimensional scOU process with a "))
    }
    if (x$root.state$random){
      if (x$root.state$stationary.root){
        cat("random stationary root.\n\n") 
      } else {
        cat("random root.\n\n") 
      }
      cat("Root expectations:\n")
      print(x$root.state$exp.root)
      cat("\n")
      cat("Root variance:\n")
      print(as.matrix(x$root.state$var.root))
      cat("\n")
    } else {
      cat("fixed root.\n\n") 
      cat("\nRoot value:\n")
      print(x$root.state$value.root)
      cat("\n")
    }
    cat("Process variance:\n")
    print(as.matrix(x$variance))
    cat("\n")
    cat("Process selection strength:\n")
    print(as.matrix(x$selection.strength))
    cat("\n")
    cat("Process root optimal values:\n")
    print(x$optimal.value)
    cat("\n")
  }
  
  if (length(x$shifts$edges) > 0){
    cat(paste0("Shifts positions on branches: ", paste(x$shifts$edges, collapse = ", "), " \n"))
    cat("Shifts values:\n")
    colnames(x$shifts$values) <- x$shifts$edges
    print(x$shifts$values)
    cat("\n")
  } else {
    cat(paste0("This process has no shift.\n"))
  }
}

##
#' @title Create an object \code{params_process}
#'
#' @description
#' \code{params_process} creates a coherent object params_process from user 
#' provided values of the parameters.
#'
#' @param x one of "BM" or "OU"
#' @param ... specified parameters, see functions \code{\link{params_BM}} and 
#' \code{\link{params_OU}} for details.
#' 
#' @return an object of class \code{params_process}.
#' 
#' @seealso \code{\link{params_BM}}, \code{\link{params_OU}}
#' 
#' @export
#'
##
params_process.character <- function(x, ...){
  if (x == "BM"){
    res <- params_BM(...)
  } else if (x == "OU"){
    res <- params_OU(...)
  }
  return(res)
}

##
#' @title Create an object \code{params_process} for a BM
#'
#' @description
#' \code{params_BM} creates a coherent object params_process from user 
#' provided values of the parameters. Non specified parameters are set to 
#' default values.
#'
#' @param p the dimension (number of traits) of the parameters. Default to 1.
#' @param variance the variance (rate matrix) of the BM. Default to 
#' \code{diag(1, p, p)}.
#' @param random whether the root of the BM is random (TRUE) or fixed (FALSE).
#' Default to FALSE.
#' @param value.root if random=FALSE, the root value. Default to 0.
#' @param exp.root if random=TRUE, the root expectation. Default to 0.
#' @param var.root if random=TRUE, the root variance. Default to
#' \code{diag(1, p, p)}.
#' @param edges a vector of edges where the shifts occur. Default to NULL
#' (no shift).
#' @param values a matrix of shift values, with p lines and as many columns as
#' the number of shifts. Each column is the p values for one shift. Default to
#' \code{matrix(0, p, length(edges))}.
#' @param relativeTimes (unused) the relative position of the shift on the
#' branch, between 0 (beginning of the branch) and 1 (end of the branch). Default
#' to 0.
#' @param nbr_of_shifts the number of shifts to use (randomly drawn). Use only
#' if \code{edges} is not specified. In that case, a phylogenetic tree must be
#' provided (to allow a random sampling of its edges).
#' @param phylo a phylogenetic tree of class \code{\link[ape]{phylo}}. Needed only if
#' the shifts edges are not specified, or if sBM_variance=TRUE. Default to NULL.
#' If sBM_variance=TRUE, it must have a specified value for the root branch
#' length (slot root.edge).
#' @param sBM_variance if the root is random, does it depend on the length of the
#' root edge ? (For equivalent purposes with a rescaled OU). Default to FALSE. If
#' TRUE, a phylogenetic tree with root edge length must be provided.
#' @param ... unused.
#' 
#' @return an object of class \code{params_process}.
#' 
#' @seealso \code{\link{params_process}}, \code{\link{params_OU}}
#' 
#' @export
#'
##
params_BM <- function(p = 1,
                      variance = diag(1, p, p),
                      random = FALSE,
                      value.root = rep(0, p),
                      exp.root = rep(0, p),
                      var.root = diag(1, p, p),
                      edges = NULL,
                      values = matrix(0, p, length(edges)),
                      relativeTimes = NULL,
                      nbr_of_shifts = length(edges),
                      phylo = NULL,
                      sBM_variance = FALSE,
                      ...) {
  if (random) {
    value.root <- NA
    if (sBM_variance){
      var.root <- phylo$root.edge * variance
    }
  } else {
    exp.root <- NA
    var.root <- NA
  }
  # Always start with some shifts, in case of default initialization (if number of shifts different from 0)
  if (length(edges) < nbr_of_shifts){
    miss <- nbr_of_shifts - length(edges)
    auth_edges <- which(!(1:nrow(phylo$edge) %in% edges))
    new_edges <- sample(auth_edges, miss)
    edges <- c(edges, new_edges)
    # edges <- c(edges, sample_shifts_edges(phylo, miss, part.list = subtree.list))
  }
  # If not enought values, complete with 0s
  if (is.null(values) || is.vector(values)){
    n_shifts_provided <- length(values)
  } else {
    n_shifts_provided <- ncol(values)
  }
  if (n_shifts_provided < nbr_of_shifts){
    miss <- nbr_of_shifts - ncol(values)
    values <- cbind(values, matrix(0, ncol = miss, nrow = p))
  }
  params_init = list(variance = variance,
                     root.state = list(random = random,
                                       value.root = value.root,
                                       exp.root = exp.root,
                                       var.root = var.root),
                     shifts = list(edges = edges,
                                   values = values,
                                   relativeTimes = relativeTimes))
  params_init <- check_dimensions(p,
                                  params_init$root.state,
                                  params_init$shifts,
                                  params_init$variance)
  params_init$root.state <- test.root.state(params_init$root.state, "BM")
  params_init$variance <- as(params_init$variance, "dpoMatrix")
  params_init$process <- "BM"
  class(params_init) <- "params_process"
  return(params_init)
}

##
#' @title Create an object \code{params_process} for an OU
#'
#' @description
#' \code{params_OU} creates a coherent object params_process from user 
#' provided values of the parameters. Non specified parameters are set to 
#' default values.
#'
#' @param p the dimension (number of traits) of the parameters. Default to 1.
#' @param variance the variance (rate matrix) of the BM. Default to 
#' \code{diag(1, p, p)}.
#' @param selection.strength the selection strength matrix. Default to 
#' \code{diag(1, p, p)}.
#' @param optimal.value the vector of the optimal values at the root. Default
#' to \code{rep(0, p)}.
#' @param random whether the root of the OU is random (TRUE) or fixed (FALSE).
#' Default to TRUE.
#' @param stationary.root whether the root of the OU is stationary (TRUE) or not.
#' Default to TRUE.
#' @param value.root if random=FALSE, the root value. Default to 0.
#' @param exp.root if random=TRUE, the root expectation. Default to 0. If
#' stationary.root=TRUE, default to \code{optimal.value}.
#' @param var.root if random=TRUE, the root variance. Default to
#' \code{diag(1, p, p)}. If stationary.root=TRUE, default to
#' the stationary variance computed from \code{variance} and
#' \code{selection.strength}, see function
#' \code{\link{compute_stationary_variance}}.
#' @param edges a vector of edges where the shifts occur. Default to NULL
#' (no shift).
#' @param values a matrix of shift values, with p lines and as many columns as
#' the number of shifts. Each column is the p values for one shift. Default to
#' \code{matrix(0, p, length(edges))}.
#' @param relativeTimes (unused) the relative position of the shift on the
#' branch, between 0 (beginning of the branch) and 1 (end of the branch). Default
#' to 0.
#' @param nbr_of_shifts the number of shifts to use (randomly drawn). Use only
#' if \code{edges} is not specified. In that case, a phylogenetic tree must be
#' provided (to allow a random sampling of its edges).
#' @param phylo a phylogenetic tree of class \code{\link[ape]{phylo}}. Needed only if
#' the shifts edges are not specified.
#' @param ... unused.
#' 
#' @return an object of class \code{params_process}.
#' 
#' @seealso \code{\link{params_process}}, \code{\link{params_BM}}
#' 
#' @export
#'
##
params_OU <- function(p = 1,
                      variance = diag(1, p, p),
                      selection.strength = diag(1, p, p),
                      optimal.value = rep(0, p),
                      random = TRUE,
                      stationary.root = TRUE,
                      value.root = rep(0, p),
                      exp.root = rep(0, p),
                      var.root = diag(1, p, p),
                      edges = NULL,
                      values = matrix(0, p, length(edges)),
                      relativeTimes = NULL,
                      nbr_of_shifts = length(edges),
                      phylo = NULL,
                      ...) {
  if (random) {
    value.root <- NA
    if (stationary.root) {
      exp.root <- optimal.value
      var.root <- compute_stationary_variance(variance, selection.strength)
    }
  } else {
    exp.root=NA
    var.root=NA
  }
  # Always start with some shifts, in case of default initialization (if number of shifts different from 0)
  if (length(edges) < nbr_of_shifts){
    miss <- nbr_of_shifts - length(edges)
    auth_edges <- which(!(1:nrow(phylo$edge) %in% edges))
    new_edges <- sample(auth_edges, miss)
    edges <- c(edges, new_edges)
    # edges <- c(edges, sample_shifts_edges(phylo, miss, part.list = subtree.list))
  }
  # If not enought values, complete with 0s
  if (is.null(values) || is.vector(values)){
    n_shifts_provided <- length(values)
  } else {
    n_shifts_provided <- ncol(values)
  }
  if (n_shifts_provided < nbr_of_shifts){
    miss <- nbr_of_shifts - ncol(values)
    values <- cbind(values, matrix(0, ncol = miss, nrow = p))
  }
  params_init <- list(variance = variance,
                      root.state = list(random = random,
                                        stationary.root = stationary.root,
                                        value.root = value.root,
                                        exp.root = exp.root,
                                        var.root = var.root),
                      shifts = list(edges = edges,
                                    values = values,
                                    relativeTimes = relativeTimes),
                      selection.strength = selection.strength,
                      optimal.value = optimal.value)
  params_init <- check_dimensions(p,
                                  params_init$root.state,
                                  params_init$shifts,
                                  params_init$variance,
                                  params_init$selection.strength,
                                  params_init$optimal.value)
  params_init$root.state <- test.root.state(params_init$root.state, "OU",
                                            variance = variance,
                                            selection.strength = params_init$selection.strength,
                                            optimal.value = optimal.value)
  params_init$variance <- as(params_init$variance, "dpoMatrix")
  params_init$process <- "OU"
  class(params_init) <- "params_process"
  return(params_init)
}


#########################################################
## Lasso initializations
#########################################################

##
#' @title Do a lasso regression with the number of non-zero variables fixed.
#'
#' @description
#' \code{lasso_regression_K_fixed} does the following regression :
#' ||Yp-Xp.delta|| + lambda |delta|_1 using the function \code{glmnet::glmnet} of 
#' package \code{glmnet}, where delta is a vector representing the shifts 
#' occurring on the branches. It does a gauss lasso regression using function 
#' \code{lm} on top of it. This function is used in functions 
#' \code{init.EM.lasso}, \code{segmentation.OU.specialCase.lasso}, ...
#'
#' @details
#' lambda is chosen so that delta has the right number of non zero components.
#' If not possible, either temporarily raise the number of shifts and then select
#' only the shifts with the highest modulus, or if not possible, throw an error.
#'
#' @param Yp (transformed) data
#' @param Xp (transformed) matrix of regression
#' @param K number of non-zero components allowed
#' 
#' @return E0.gauss the intercept (value at the root)
#' @return shifts.gauss the list of shifts found on the branches
#' 
#' @keywords internal
#'
#06/10/14 - Initial release
##

lasso_regression_K_fixed.glmnet_multivariate <- function(Yp, Xp, K,
                                                         root = NULL,
                                                         penscale = rep(1, ncol(Xp)),
                                                         K_lag = 0) {
  # library(glmnet)
  ## Lag
  K_original <- K
  K <- K + K_lag
  ## Dim
  p <- nrow(Yp)
  ## Root is the intercept, should be excluded from varaiable selection
  # In that case, project Yp on the orthogonal of the root
  if (!is.null(root)){
    #     L <- Xp[ , root]
    #     norme_L <- drop(crossprod(L))
    #     Xp_noroot <- Xp[ , -root, drop = FALSE]
    #     Xp_orth <- Xp_noroot - (tcrossprod(L) %*% Xp_noroot) / norme_L
    #     Yp_orth <- Yp - crossprod(Yp, L) / (norme_L) * L
    #     intercept <- FALSE
    #     penscale <- penscale[-root]
    Xp_orth <- Xp
    Yp_orth <- Yp
    intercept <- FALSE
    penscale[root] <- 0
    # K_original <- K
    K <- K + 1 ## the root does not count in the number of non zero parameters.
  } else {
    Xp_orth <- Xp
    Yp_orth <- Yp
    intercept <- TRUE
    # K_original <- K
  }
  ## fit
  fit <- glmnet::glmnet(x = 0 + Xp_orth,
                        y = t(Yp_orth),
                        family = "mgaussian",
                        alpha = 1,
                        nlambda = 500,
                        intercept = intercept,
                        dfmax = K + 10,
                        penalty.factor = penscale)
  df <- fit$df
  ## Find the lambda that gives the right number of ruptures
  # Check that lambda goes far enought
  if (K > max(df)) {
    fit <- glmnet::glmnet(x = 0 + Xp_orth,
                          y = t(Yp_orth),
                          family = "mgaussian",
                          alpha = 1,
                          nlambda = 500,
                          intercept = intercept,
                          lambda.min.ratio = 0,
                          penalty.factor = penscale)
    df <- fit$df
  }
  if (K > max(df)) {
    stop("Lasso regression failed. There are too many variables.")
  }
  # ## If the right lambda does not exists, find it.
  # count <- 0
  # fit_tmp <- fit
  # while (!any(df == K) && count < 500) {
  #   count <- count + 1
  #   K_inf <- max(K - 2, min(df))
  #   while (!any(K_inf == df) && (K_inf >= 0)) {
  #     K_inf <- K_inf - 1
  #   }
  #   lambda_inf <- fit$lambda[tail(which(K_inf == df), n = 1)]
  #   K_sup <- min(K + 2, max(df))
  #   if (K > max(df)){
  #     fit <- fit_tmp
  #     df <- fit_tmp$df
  #     break
  #   }
  #   while (!any(K_sup == df) && (K_sup <= max(df))) {
  #     K_sup <- K_sup + 1
  #   }
  #   lambda_sup <- fit$lambda[head(which(K_sup == df), n = 1)]
  #   lambda <- seq(from = lambda_inf, to = lambda_sup, length.out = 100)
  #   fit_tmp <- fit
  #   fit <- glmnet::glmnet(x = 0 + Xp_orth,
  #                 y = t(Yp_orth),
  #                 family = "mgaussian",
  #                 alpha = 1,
  #                 lambda = lambda,
  #                 intercept = intercept,
  #                 penalty.factor = penscale)
  #   df <- fit$df
  # }
  # rm(fit_tmp)
  # ## If the right lambda does not exists, raise the number of shifts
  # K_2 <- K
  # while (!any(df == K_2) && K_2 <= min(dim(Xp))) {
  #   if (K_2 == K){
  #     warning("During lasso regression, could not find the right lambda for the number of shifts K. Temporarly raised it to do the lasso regression, and furnishing the K largest coefficients.")
  #   }
  #   K_2 <- K_2 + 1
  # }
  ## If the right lambda does not exists, raise the number of shifts
  K_2 <- K
  if (K_2 <= max(df)){ # If K < max, raise it till find a right one
    while (!any(df == K_2) && K_2 < min(dim(Xp))) {
      K_2 <- K_2 + 1
    }
  } else { # If K_original < max(df) but K > max(df), decrease K
    while (!any(df == K_2) && K_2 > K_original + 1) {
      K_2 <- K_2 - 1
    }
  }
  ## If could not find the right lambda, do a default initialization
  if (!any(df == K_2)) {
    stop("Lasso initialization failed : could not find a satisfying number of shifts.")
  }
  ## Select the row with the right number of coefficients
  index <- min(which(df == K_2))
  if (p == 1){
    delta <- fit$beta[, index]
    delta <- matrix(delta, nrow = length(delta))
  } else {
    delta <- sapply(fit$beta, function(z) z[, index])
  }
  # Check that the matrix is of full rank
  if (K_lag == 0){ # Only if K_lag = 0
    projection <- which(rowSums(delta) != 0)
    Xproj <- Xp[ , projection, drop = FALSE]
    if (dim(Xproj)[2] != qr(Xproj)$rank) {
      warning("The solution fund by lasso had non independent vectors. Had to modify this solution.")
      # Re-do a fit and try again.
      fit <- glmnet::glmnet(x = 0 + Xp_orth,
                            y = t(Yp_orth),
                            family = "mgaussian",
                            alpha = 1,
                            nlambda = 500,
                            intercept = intercept,
                            dfmax = K + 10,
                            penalty.factor = penscale)
      delta <- try(find_independent_regression_vectors.glmnet_multivariate(Xp, K, fit, root))
      if (inherits(delta, "try-error")) stop("The selected variables do not produce a full rank regression matrix !")
    }
  }
  # ## If we had to raise the number of shifts, go back to the initial number, taking the K largest shifts
  # if (is.null(root)){
  #   edges <- order(-rowSums(abs(delta)))[1:K]
  #   delta.bis <- matrix(0, dim(delta)[1], dim(delta)[2])
  #   delta.bis[edges, ] <- delta[edges, ]
  # } else {
  #   edges <- order(-rowSums(abs(delta[-root, , drop = F])))[1:K_original]
  #   delta.bis <- matrix(0, dim(delta)[1] - 1, dim(delta)[2])
  #   delta.bis[edges, ] <- delta[edges, ]
  # }
  # ## Gauss lasso
  # return(compute_gauss_lasso(t(Yp), Xp, delta.bis, root))
  ## Gauss lasso, going back to the right number of shifts
  if (is.null(root)){
    vals <- rowSums(abs(delta))
    edges <- unname(which(vals != 0))
  } else {
    vals <- rowSums(abs(delta[-root, , drop = F]))
    edges <- unname(which(vals != 0))
  }
  ## Handle case wher combinatoire is too expensive
  K_choose <- K_original
  ed_or <- order(-vals)
  n_ed <- length(edges)
  fixed_edges <- NULL
  while((choose(length(edges), K_choose) > 50000) && (K_choose > 0)){
    K_choose <- K_choose - 1
    fixed_edges <- ed_or[1:(K_original - K_choose)]
    edges <- ed_or[(K_original - K_choose + 1):n_ed]
  }
  posibilities <- combn(edges, K_choose, simplify = FALSE) # All possible combinaisons
  Ypt <- t(Yp) # Do the transpose only once
  fun <- function(posi){
    posi <- c(fixed_edges, posi)
    if (is.null(root)){
      delta.bis <- matrix(0, dim(delta)[1], dim(delta)[2])
      delta.bis[posi, ] <- delta[posi, ]
    } else {
      delta.bis <- matrix(0, dim(delta)[1] - 1, dim(delta)[2])
      delta.bis[posi, ] <- delta[posi, ]
    }
    return(suppressWarnings(compute_gauss_lasso(Ypt, Xp, delta.bis, root, projection = posi)))
  }
  res_try <- lapply(posibilities, fun)
  scores <- sapply(res_try, function(z) sum(z$residuals^2))
  ## Check that the matrix is of full rank
  order_scores <- order(scores)
  for (os in order_scores){
    res <- res_try[[os]]
    projection <- which(res$beta[, index] != 0)
    Xproj <- Xp[ , projection, drop = FALSE]
    if (dim(Xproj)[2] == qr(Xproj)$rank) return(res)
  }
  stop(paste0("At lasso regression, could not find K=", K_original, " independent edges (to produce a full rank matrix)"))
}

# lasso_regression_K_fixed.grplasso <- function(Yvec, Xkro, K,
#                                               root,
#                                               penscale = rep(1,
#                                                              length(unique(group))),
#                                               group = 1:ncol(Xkro),
#                                               p_dim,
#                                               K_lag = 0) {
#   ## Lag
#   K_original <- K
#   K <- K + K_lag
#   ## Root is the intercept, should be excluded from varaiable selection
#   Xp_orth <- Xkro
#   Yp_orth <- Yvec
#   index <- group
#   index[index == root] <- NA
#   K <- K + 1 ## the root does not count in the number of non zero parameters.
#   ## fit
#   lambda_max <- grplasso::lambdamax(Xp_orth, Yp_orth, model = grplasso::LinReg(),
#                                     index = group, center = FALSE, standardize = FALSE)
#   lambda_min <- 0.01 * lambda_max
#   lambda_grid <- exp(seq(log(lambda_max), log(lambda_min), length.out = 100))
#   co <- capture.output(
#     suppressWarnings(
#     fit <- grplasso::grplasso(x = Xp_orth,
#                               y = Yp_orth,
#                               index = group,
#                               model = grplasso::LinReg(),
#                               lambda = lambda_grid,
#                               center = FALSE, standardize = FALSE)
#     )
#   )
#   df <- apply(fit$coefficients, 2, function(z) length(unique(group[z != 0])))
#   ## Find the lambda that gives the right number of ruptures
#   # Check that lambda goes far enought
#   if (K_original + 1 > max(df)) {
#     lambda_min <- 0.0001 * lambda_max
#     lambda_grid <- exp(seq(log(lambda_max), log(lambda_min), length.out = 100))
#     co <- capture.output(
#       suppressWarnings(
#       fit <- grplasso::grplasso(x = Xp_orth,
#                                 y = Yp_orth,
#                                 index = group,
#                                 model = grplasso::LinReg(),
#                                 lambda = lambda_grid,
#                                 center = FALSE, standardize = FALSE)
#       )
#     )
#     df <- apply(fit$coefficients, 2, function(z) length(unique(group[z != 0])))
#   }
#   if (K_original + 1 > max(df)) {
#     stop("Lasso regression failed. There are too many variables.")
#   }
#   ## If the right lambda does not exists, raise the number of shifts
#   K_2 <- K
#   if (K_2 <= max(df)){ # If K < max, raise it till find a right one
#     while (!any(df == K_2) && K_2 < max(group)) {
#       K_2 <- K_2 + 1
#     }
#   } else { # If K_original < max(df) but K > max(df), decrease K
#     while (!any(df == K_2) && K_2 > K_original + 1) {
#       K_2 <- K_2 - 1
#     }
#   }
#   ## If could not find the right lambda, do a default initialization
#   if (!any(df == K_2)) {
#     stop("Lasso initialization failed : could not find a satisfying number of shifts.")
#   }
#   ## Select the row with the right number of coefficients
#   index <- min(which(df == K_2))
#   if (p_dim == 1){
#     delta <- fit$coefficients[, index]
#     delta <- matrix(delta, nrow = length(delta))
#   } else {
#     delta <- fit$coefficients[, index]
#     delta <- t(matrix(delta, nrow = p_dim))
#   }
#   ## Check that the matrix is of full rank
#   if (K_lag == 0){
#     projection <- which(fit$coefficients[, index] != 0)
#     Xproj <- Xkro[ , projection, drop = FALSE]
#     if (dim(Xproj)[2] != qr(Xproj)$rank) {
#       warning("The solution fund by lasso had non independent vectors. Had to modify this solution.")
#       # Try again.
#       delta <- try(find_independent_regression_vectors.grplasso(Xkro, K,
#                                                                 fit, root,
#                                                                 p_dim, group))
#       if (inherits(delta, "try-error")) stop("The selected variables do not produce a full rank regression matrix !")
#     }
#   }
#   vals <- rowSums(abs(delta[-root, , drop = F]))
#   edges <- unname(which(vals != 0))
#   ## Handle case wher combinatoire is too expensive
#   K_choose <- K_original
#   ed_or <- order(-vals)
#   n_ed <- length(edges)
#   fixed_edges <- NULL
#   while((choose(length(edges), K_choose) > 50000) && (K_choose > 0)){
#     K_choose <- K_choose - 1
#     fixed_edges <- ed_or[1:(K_original - K_choose)]
#     edges <- ed_or[(K_original - K_choose + 1):n_ed]
#   }
#   posibilities <- combn(edges, K_choose, simplify = FALSE) # All possible combinaisons
#   fun <- function(posi){
#     posi <- c(fixed_edges, posi)
#     if (is.null(root)){
#       delta.bis <- matrix(0, dim(delta)[1], dim(delta)[2])
#       delta.bis[posi, ] <- delta[posi, ]
#     } else {
#       delta.bis <- matrix(0, dim(delta)[1] - 1, dim(delta)[2])
#       delta.bis[posi, ] <- delta[posi, ]
#     }
#     return(suppressWarnings(compute_gauss_lasso.grplasso(Yvec, Xkro, delta.bis,
#                                                          root, group, p_dim)))
#   }
#   res_try <- lapply(posibilities, fun)
#   scores <- sapply(res_try, function(z) sum(z$residuals^2))
#   ## Check that the matrix is of full rank
#   order_scores <- order(scores)
#   for (os in order_scores){
#     res <- res_try[[os]]
#     projection <- which(res$beta[, index] != 0)
#     Xproj <- Xkro[ , projection, drop = FALSE]
#     if (dim(Xproj)[2] == qr(Xproj)$rank) return(res)
#   }
#   stop(paste0("At lasso regression, could not find K=", K_original, " independent edges (to produce a full rank matrix)"))
# }

lasso_regression_K_fixed.gglasso <- function(Yvec, Xkro, K,
                                             root = NULL,
                                             penscale = rep(1,
                                                            length(unique(group))),
                                             group = 1:ncol(Xkro),
                                             p_dim,
                                             K_lag = 0) {
  ## Lag
  K_original <- K
  K <- K + K_lag
  ## Root is the intercept, should be excluded from varaiable selection
  if (!is.null(root)){
    Xp_orth <- Xkro
    Yp_orth <- Yvec
    intercept <- FALSE
    penscale[root] <- 0
    K <- K + 1 ## the root does not count in the number of non zero parameters.
  } else {
    Xp_orth <- Xkro
    Yp_orth <- Yvec
    intercept <- TRUE
  }
  ## fit
  co <- capture.output(
    fit <- gglasso::gglasso(x = Xp_orth,
                            y = Yp_orth,
                            group = group,
                            loss = "ls",
                            nlambda = 500,
                            intercept = intercept,
                            dfmax = K + 10,
                            pf = penscale)
  )
  df <- apply(fit$beta, 2, function(z) length(unique(group[z != 0])))
  ## Find the lambda that gives the right number of ruptures
  # Check that lambda goes far enought
  if (K_original + 1 > max(df)) {
    co <- capture.output(
      fit <- gglasso::gglasso(x = Xp_orth,
                              y = Yp_orth,
                              group = group,
                              loss = "ls",
                              nlambda = 500,
                              intercept = intercept,
                              dfmax = K + 10,
                              pf = penscale,
                              lambda.factor = 0)
    )
    df <- apply(fit$beta, 2, function(z) length(unique(group[z != 0])))
  }
  if (K_original + 1 > max(df)) {
    stop("Lasso regression failed. There are too many variables.")
  }
  # ## If the right lambda does not exists, find it.
  # count <- 0
  # fit_tmp <- fit
  # while (!any(df == K) && count < 500) {
  #   count <- count + 1
  #   K_inf <- max(K - 2, min(df))
  #   while (!any(K_inf == df) && (K_inf >= 0)) {
  #     K_inf <- K_inf - 1
  #   }
  #   lambda_inf <- fit$lambda[tail(which(K_inf == df), n = 1)]
  #   K_sup <- min(K + 2, max(df))
  #   if (K > max(df)){
  #     fit <- fit_tmp
  #     df <- fit_tmp$df
  #     break
  #   }
  #   while (!any(K_sup == df) && (K_sup <= max(df))) {
  #     K_sup <- K_sup + 1
  #   }
  #   lambda_sup <- fit$lambda[head(which(K_sup == df), n = 1)]
  #   lambda <- seq(from = lambda_inf, to = lambda_sup, length.out = 100)
  #   fit_tmp <- fit
  #   co <- capture.output(
  #   fit <- gglasso::gglasso(x = Xp_orth,
  #                  y = Yp_orth,
  #                  group = group,
  #                  loss = "ls",
  #                  intercept = intercept,
  #                  pf = penscale,
  #                  lambda = lambda)
  #   )
  #   df_prev <- df
  #   df <- apply(fit$beta, 2, function(z) length(unique(group[z != 0])))
  #   if (identical(df, df_prev)) break
  # }
  # rm(fit_tmp)
  ## If the right lambda does not exists, raise the number of shifts
  K_2 <- K
  if (K_2 <= max(df)){ # If K < max, raise it till find a right one
    while (!any(df == K_2) && K_2 < max(group)) {
      K_2 <- K_2 + 1
    }
  } else { # If K_original < max(df) but K > max(df), decrease K
    while (!any(df == K_2) && K_2 > K_original + 1) {
      K_2 <- K_2 - 1
    }
  }
  ## If could not find the right lambda, do a default initialization
  if (!any(df == K_2)) {
    stop("Lasso initialization failed : could not find a satisfying number of shifts.")
  }
  ## Select the row with the right number of coefficients
  index <- min(which(df == K_2))
  if (p_dim == 1){
    delta <- fit$beta[, index]
    delta <- matrix(delta, nrow = length(delta))
  } else {
    delta <- fit$beta[, index]
    delta <- t(matrix(delta, nrow = p_dim))
  }
  ## Check that the matrix is of full rank
  if (K_lag == 0){
    projection <- which(fit$beta[, index] != 0)
    Xproj <- Xkro[ , projection, drop = FALSE]
    if (dim(Xproj)[2] != qr(Xproj)$rank) {
      warning("The solution fund by lasso had non independent vectors. Had to modify this solution.")
      # Re-do a fit and try again.
      co <- capture.output(
        fit <- gglasso::gglasso(x = Xp_orth,
                                y = Yp_orth,
                                group = group,
                                loss = "ls",
                                nlambda = 500,
                                intercept = intercept,
                                dfmax = K + 10,
                                pf = penscale)
      )
      delta <- try(find_independent_regression_vectors.gglasso(Xkro, K,
                                                               fit, root,
                                                               p_dim, group))
      if (inherits(delta, "try-error")) stop("The selected variables do not produce a full rank regression matrix !")
    }
  }
  if (is.null(root)){
    vals <- rowSums(abs(delta))
    edges <- unname(which(vals != 0))
  } else {
    vals <- rowSums(abs(delta[-root, , drop = F]))
    edges <- unname(which(vals != 0))
  }
  ## Handle case wher combinatoire is too expensive
  K_choose <- K_original
  ed_or <- order(-vals)
  n_ed <- length(edges)
  fixed_edges <- NULL
  while((choose(length(edges), K_choose) > 50000) && (K_choose > 0)){
    K_choose <- K_choose - 1
    fixed_edges <- ed_or[1:(K_original - K_choose)]
    edges <- ed_or[(K_original - K_choose + 1):n_ed]
  }
  posibilities <- combn(edges, K_choose, simplify = FALSE) # All possible combinaisons
  fun <- function(posi){
    posi <- c(fixed_edges, posi)
    if (is.null(root)){
      delta.bis <- matrix(0, dim(delta)[1], dim(delta)[2])
      delta.bis[posi, ] <- delta[posi, ]
    } else {
      delta.bis <- matrix(0, dim(delta)[1] - 1, dim(delta)[2])
      delta.bis[posi, ] <- delta[posi, ]
    }
    return(suppressWarnings(compute_gauss_lasso.gglasso(Yvec, Xkro, delta.bis,
                                                        root, group, p_dim)))
  }
  res_try <- lapply(posibilities, fun)
  scores <- sapply(res_try, function(z) sum(z$residuals^2))
  ## Check that the matrix is of full rank
  order_scores <- order(scores)
  for (os in order_scores){
    res <- res_try[[os]]
    projection <- which(res$beta[, index] != 0)
    Xproj <- Xkro[ , projection, drop = FALSE]
    if (dim(Xproj)[2] == qr(Xproj)$rank) return(res)
  }
  stop(paste0("At lasso regression, could not find K=", K_original, " independent edges (to produce a full rank matrix)"))
}


# lasso_regression_K_fixed.glmnet <- function (Yp, Xp, K, intercept.penalty = FALSE ) {
#   ## Penalty on the first coordinate = intercept : force first cooerdinate to be null
#   excl <- NULL
#   if (intercept.penalty) excl <- c(1)
#   ## fit
#   fit <- glmnet::glmnet(x = 0 + Xp, y = Yp, alpha = 1, exclude = excl)
#   ## Find the lambda that gives the right number of ruptures
#   # Check that lambda goes far enought
#   if (K > max(fit$df)) {
#     fit <- glmnet::glmnet(x = 0 + Xp, y = Yp, alpha = 1, lambda.min.ratio = 0, exclude = excl)
#   }
#   if (K > max(fit$df)) {
#     stop("Lasso regression failed. There are too many variables.")
#   }
#   ## If the right lambda does not exists, find it.
#   count <- 0
#   while (sum(fit$df == K) == 0 && count < 500) {
#     count <- count + 1
#     K_inf <- K-1
#     while ((sum(K_inf == fit$df) == 0) && (K_inf >= 0)) {
#       K_inf <- K_inf - 1
#     }
#     lambda_inf <- fit$lambda[tail(which(K_inf == fit$df), n=1)]
#     K_sup <- K + 1
#     while ((sum(K_sup == fit$df) == 0) && (K_sup <= max(fit$df))) {
#       K_sup <- K_sup + 1
#     }
#     lambda_sup <- fit$lambda[head(which(K_sup == fit$df), n=1)]
#     lambda <- seq(from = lambda_inf, to = lambda_sup, length.out = 100)
#     fit <- glmnet::glmnet(x = 0 + Xp, y = Yp, alpha = 1, lambda = lambda, exclude = excl)
#   }
#   ## If the right lambda does not exists, raise the number of shifts
#   K_2 <- K
#   while (sum(fit$df == K_2) == 0 && K_2 < 500) {
#     warning("During lasso regression, could not find the right lambda for the number of shifts K. Temporarly raised it to do the lasso regression, and furnishing the K largest coefficients.")
#     K_2 <- K_2 + 1
#   }
#   ## If could not find the right lambda, do a default initialization
#   if (sum(fit$df == K_2) == 0) {
#     stop("Lasso initialization fail : could not find a satisfying number of shifts.")
#   } else {
#     delta <- coef(fit, s = fit$lambda[min(which(fit$df == K_2))])
#     E0 <- delta[1]; # Intercept
#     delta <- delta[-1];
#     ## Gauss lasso
#     projection <- which(delta != 0)
#     Xproj <- 0 + Xp[, projection]
#     fit.gauss <- lm(Yp ~ Xproj)
#     delta.gauss <- rep(0, dim(Xp)[2])
#     E0.gauss <- coef(fit.gauss)[1]; names(E0.gauss) <- NULL
#     delta.gauss[projection] <- coef(fit.gauss)[-1]
#     # If lm fails to find some coeeficients, put them to 0
#     delta.gauss[is.na(delta.gauss)] <- 0.1
#     ## If we had to raise the number of shifts, go back to the initial number, taking the K largest shifts
#     edges <- order(-abs(delta.gauss))[1:K]
#     delta.gauss.final <- rep(0, length(delta.gauss))
#     delta.gauss.final[edges] <- delta.gauss[edges]
#     shifts.gauss <- shifts.vector_to_list(delta.gauss.final);
#     return(list(E0.gauss = E0.gauss, shifts.gauss = shifts.gauss))
#   }
# }

# lasso_regression_K_fixed.quadrupen <- function (Yp, Xp, K, root = NULL, penscale = rep(1, ncol(Xp))) {
#   ## Root is the intercept, should be excluded from varaiable selection
#   # In that case, project Yp on the orthogonal of the root
#   if (!is.null(root)){
#     L <- Xp[ , root]
#     norme_L <- drop(crossprod(L))
#     Xp_noroot <- Xp[ , -root, drop = FALSE]
#     Xp_orth <- Xp_noroot - (tcrossprod(L) %*% Xp_noroot) / norme_L
#     Yp_orth <- Yp - crossprod(Yp, L) / (norme_L) * L
#     intercept <- FALSE
#     penscale <- penscale[-root]
#   } else {
#     Xp_orth <- Xp
#     Yp_orth <- Yp
#     intercept <- TRUE
#   }
#   ## fit
#   fit <- elastic.net(x = 0 + Xp_orth, y = Yp_orth, lambda2 = 0, nlambda1 = 500, intercept = intercept, max.feat = K + 10, penscale = penscale)
#   df <- rowSums(fit@active.set)
#   ## Find the lambda that gives the right number of ruptures
#   # Check that lambda goes far enought
#   if (K > max(df)) {
#     fit <- elastic.net(x = 0 + Xp_orth, y = Yp_orth, lambda2 = 0, min.ratio = 10^(-10), intercept = intercept, penscale = penscale)
#     df <- rowSums(fit@active.set)
#   }
#   if (K > max(df)) {
#     stop("Lasso regression failed. There are too many variables.")
#   }
#   ## If the right lambda does not exists, find it.
#   count <- 0
#   while (!any(df == K) && count < 500) {
#     count <- count + 1
#     K_inf <- K - 1
#     while (!any(K_inf == df) && (K_inf >= 0)) {
#       K_inf <- K_inf - 1
#     }
#     lambda_inf <- fit@lambda1[tail(which(K_inf == df), n = 1)]
#     K_sup <- K + 1
#     while (!any(K_sup == df) && (K_sup <= max(df))) {
#       K_sup <- K_sup + 1
#     }
#     lambda_sup <- fit@lambda1[head(which(K_sup == df), n = 1)]
#     lambda <- seq(from = lambda_inf, to = lambda_sup, length.out = 100)
#     fit <- elastic.net(x = 0 + Xp_orth, y = Yp_orth, lambda1 = lambda, lambda2 = 0, intercept = intercept)
#     df <- rowSums(fit@active.set)
#   }
#   ## If the right lambda does not exists, raise the number of shifts
#   K_2 <- K
#   while (!any(df == K_2) && K_2 <= min(dim(Xp))) {
#     warning("During lasso regression, could not find the right lambda for the number of shifts K. Temporarly raised it to do the lasso regression, and furnishing the K largest coefficients.")
#     K_2 <- K_2 + 1
#   }
#   ## If could not find the right lambda, do a default initialization
#   if (!any(df == K_2)) {
#     stop("Lasso initialization failed : could not find a satisfying number of shifts.")
#   }
#   ## Select the row with the right number of coefficients
#   index <- min(which(df == K_2))
#   delta <- fit@coefficients[index,]
#   # If we put aside the root, replace it in the coefficients
#   if (!is.null(root)){
#     #deltabis <- unname(c(delta[1:(root - 1)], 0, tail(delta, n = max(0, length(delta) - root + 1))))
#     delta <- append(delta, 0, after = root - 1)
#   }
#   # Check that the matrix is of full rank
#   projection <- which(delta != 0)
#   Xproj <- Xp[ , projection, drop = FALSE]
#   if (dim(Xproj)[2] != qr(Xproj)$rank) {
#     warning("The solution fund by lasso had non independent vectors. Had to modify this solution.")
#     # Re-do a fit and try again.
#     fit <- elastic.net(x = 0 + Xp_orth, y = Yp_orth, lambda2 = 0, nlambda1 = 500, intercept = intercept, max.feat = K + 10)
#     delta <- try(find_independent_regression_vectors(Xp, K, fit, root))
#     if (inherits(delta, "try-error")) stop("The selected variables do not produce a full rank regression matrix !")
#   }
#   ## If we had to raise the number of shifts, go back to the initial number, taking the K largest shifts
#   edges <- order(-abs(delta))[1:K]
#   delta.bis <- rep(0, length(delta))
#   delta.bis[edges] <- delta[edges]
#   ## Gauss lasso
#   return(compute_gauss_lasso(t(Yp), Xp, delta.bis, root))
# }

##
#' @title Given a regularization path, find K selected independent variables.
#'
#' @description
#' \code{find_independent_regression_vectors} tries to find a situation where K variables
#' are selected, so that the selected columns of matrix Xp are independent.
#'
#' @details
#' To do that, if a set of selected is not independent, we go back to the previous selected 
#' variables, and forbid the moves that led to non-independence for the rest of the path.
#'
#' @param Xp (transformed) matrix of regression
#' @param K number of non-zero components allowed
#' @param root integer, position of the root column (intercept) excluded from the fit. null 
#' if no root column.
#' 
#' @return delta a vector of regression with K non-zero coefficients.
#' 
#' @keywords internal
#'
##
find_independent_regression_vectors.glmnet_multivariate <- function(Xp, K, fit, root){
  if (!is.list(fit$beta)){
    p <- 1
  } else {
    p <- length(fit$beta)
  }
  if (p == 1){
    deltas <- fit$beta
    deltas <- array(deltas, dim = c(1, dim(deltas)))
  } else {
    # library(plyr)
    deltas <- plyr::laply(fit$beta, function(z) as.matrix(z))
  }
  nsets <- dim(deltas)[3]
  #   if (!is.null(root)){
  #     deltas <- apply(deltas, 1, function(z) append(z, 0, after = root - 1))
  #   }
  projections <- t(apply(deltas, 3, function(z) return(colSums(z) != 0)))
  check_independence <- function(projection, Xp){
    Xproj <- Xp[ , projection, drop = FALSE]
    return(dim(Xproj)[2] == qr(Xproj)$rank)
  }
  for (i in 1:nsets){
    # If not independent : go back to the previous state.
    if (!check_independence(projections[i, ], Xp)){
      # Variables that were activated or inactivated
      changes <- xor(projections[i - 1, ], projections[i, ])
      # Activated variables : inactivate them for the futur
      new_vars <- changes & projections[i, ]
      projections[i:nsets, new_vars] <- FALSE
      # Inactivated variables : re-activate them for the futur
      del_vars <- changes & projections[i - 1, ]
      projections[i:nsets, del_vars] <- TRUE
    }
  }
  ## Find the right number of selected variables
  n_select <- rowSums(projections)
  right_ones <- n_select >= K
  if (!any(right_ones)){
    stop("Could not find K independent vectors in the regression path provided.")
  } else {
    right_one <- which(right_ones)[1]
    ## If too many, take the K largests.
    delta_ind <- t(matrix(deltas[, , right_one], nrow = dim(deltas)[1]))
    edges <- which(projections[right_one, ])
    delta.bis <- matrix(0, dim(delta_ind)[1], dim(delta_ind)[2])
    values <- delta_ind[edges, , drop = FALSE]
    values[rowSums(values) == 0, ] <- 1 # If projection selected edges not initially present
    delta.bis[edges, ] <- values
    return(delta.bis)
#     return(matrix(rep(0 + projections[right_one, ], 3),
#                   nrow = length(projections[right_one, ])))
  }
}

# find_independent_regression_vectors.grplasso <- function(Xkro, K, fit, root, p, group){
#   nsets <- ncol(fit$coefficients)
#   projections <- t(fit$coefficients != 0)
#   check_independence <- function(projection, Xkro){
#     Xproj <- Xkro[ , projection, drop = FALSE]
#     return(dim(Xproj)[2] == qr(Xproj)$rank)
#   }
#   for (i in 1:nsets){
#     # If not independent : go back to the previous state.
#     if (!check_independence(projections[i, ], Xkro)){
#       # Variables that were activated or inactivated
#       changes <- xor(projections[i - 1, ], projections[i, ])
#       # Activated variables : inactivate them for the futur
#       new_vars <- changes & projections[i, ]
#       projections[i:nsets, new_vars] <- FALSE
#       # Inactivated variables : re-activate them for the futur
#       del_vars <- changes & projections[i - 1, ]
#       projections[i:nsets, del_vars] <- TRUE
#     }
#   }
#   ## Find the right number of selected variables
#   n_select <- apply(projections, 1, function(z) length(unique(group[z])))
#   right_ones <- n_select >= K
#   if (!any(right_ones)){
#     stop("Could not find K independent vectors in the regression path provided.")
#   } else {
#     right_one <- which(right_ones)[1]
#     ## If too many, take the K largests.
#     delta_ind <- t(matrix(fit$coefficients[, right_one], nrow = p))
#     edges <- which(rowSums(delta_ind) != 0)
#     delta.bis <- matrix(0, dim(delta_ind)[1], dim(delta_ind)[2])
#     values <- delta_ind[edges, , drop = FALSE]
#     values[rowSums(values) == 0, ] <- 1 # If projection selected edges not initially present
#     delta.bis[edges, ] <- values
#     return(delta.bis)
#   }
# }

find_independent_regression_vectors.gglasso <- function(Xkro, K, fit, root, p, group){
  nsets <- ncol(fit$beta)
  #   if (!is.null(root)){
  #     deltas <- apply(deltas, 1, function(z) append(z, 0, after = root - 1))
  #   }
  projections <- t(fit$beta != 0)
  check_independence <- function(projection, Xkro){
    Xproj <- Xkro[ , projection, drop = FALSE]
    return(dim(Xproj)[2] == qr(Xproj)$rank)
  }
  for (i in 1:nsets){
    # If not independent : go back to the previous state.
    if (!check_independence(projections[i, ], Xkro)){
      # Variables that were activated or inactivated
      changes <- xor(projections[i - 1, ], projections[i, ])
      # Activated variables : inactivate them for the futur
      new_vars <- changes & projections[i, ]
      projections[i:nsets, new_vars] <- FALSE
      # Inactivated variables : re-activate them for the futur
      del_vars <- changes & projections[i - 1, ]
      projections[i:nsets, del_vars] <- TRUE
    }
  }
  ## Find the right number of selected variables
  # n_select <- rowSums(projections)
  n_select <- apply(projections, 1, function(z) length(unique(group[z])))
  right_ones <- n_select >= K
  if (!any(right_ones)){
    stop("Could not find K independent vectors in the regression path provided.")
  } else {
    right_one <- which(right_ones)[1]
    ## If too many, take the K largests.
    delta_ind <- t(matrix(fit$beta[, right_one], nrow = p))
    edges <- which(rowSums(delta_ind) != 0)
    delta.bis <- matrix(0, dim(delta_ind)[1], dim(delta_ind)[2])
    values <- delta_ind[edges, , drop = FALSE]
    values[rowSums(values) == 0, ] <- 1 # If projection selected edges not initially present
    delta.bis[edges, ] <- values
    return(delta.bis)
    #     return(matrix(rep(0 + projections[right_one, ], 3),
    #                   nrow = length(projections[right_one, ])))
  }
}

# find_independent_regression_vectors.quadrupen <- function(Xp, K, fit, root){
#   deltas <- fit@coefficients
#   nsets <- dim(deltas)[1]
#   if (!is.null(root)){
#     deltas <- apply(deltas, 1, function(z) append(z, 0, after = root - 1))
#   }
#   projections <- t(apply(deltas, 1, function(z) return(z != 0)))
#   check_independence <- function(projection, Xp){
#     Xproj <- Xp[ , projection, drop = FALSE]
#     return(dim(Xproj)[2] == qr(Xproj)$rank)
#   }
#   for (i in 1:nsets){
#     # If not independent : go back to the previous state.
#     if (!check_independence(projections[i, ], Xp)){
#       # Variables that were activated or inactivated
#       changes <- xor(projections[i - 1, ], projections[i, ])
#       # Activated variables : inactivate them for the futur
#       new_vars <- changes && projections[i, ]
#       projections[i:nsets, new_vars] <- 0
#       # Inactivated variables : re-activate them for the futur
#       del_vars <- changes && projections[i - 1, ]
#       projections[i:nsets, del_vars] <- 1
#     }
#   }
#   ## Find the right number of selected variables
#   n_select <- rowSums(projections)
#   right_ones <- n_select == K
#   if (!any(right_ones)){
#     stop("Could not find K independent vectors in the regression path provided.")
#   } else {
#     right_one <- which(right_ones)
#     return(projections[right_one, ])
#   }
# }

##
#' @title Do a lm on top of a lasso regression.
#'
#' @description
#' \code{compute_gauss_lasso} takes the variables selected by a lasso procedure, and
#' uses them to do a simple linear least square regression. Function used is
#' \code{lm} for non-transformed data (root = NULL), and \code{lm.fit} for
#' transformed data (root = an integer).
#'
#' @details
#' Depending on the value of root, the behavior is different. If root is null, then
#' we fit a linear regression with an intercept. If root is equal to an integer,
#' then the "intercept" column of the matrix Xp (that has possibly been trough a 
#' multiplication with a Cholesky decomposition of the variance) is included, rather
#' than the intercept.
#'
#' @param Yp (transformed) data
#' @param Xp (transformed) matrix of regression
#' @param delta regression coefficients obtained with a lasso regression
#' @param root the position of the root (intercept) in delta
#' 
#' @return Named list, with "E0.gauss" the intercept (value at the root);
#' "shifts.gauss" the list of shifts found on the branches; and "residuals" the
#' residuals of the regression
#' 
#' @keywords internal
#'
##
compute_gauss_lasso <- function (Ypt, Xp, delta, root,
                                 projection = which(rowSums(delta) != 0)) {
  # projection <- which(rowSums(delta) != 0)
  if (is.null(root)) { # If no one is excluded, "real" intercept
    Xproj <- 0 + Xp[, projection, drop = FALSE]
    fit.gauss <- lm(Ypt ~ Xproj)
    delta.gauss <- matrix(0, dim(Xp)[2], dim(Ypt)[2])
    coefs_gauss <- matrix(coef(fit.gauss), ncol = dim(Ypt)[2])
    E0.gauss <- coefs_gauss[1, ]; names(E0.gauss) <- NULL
    delta.gauss[projection, ] <- coefs_gauss[-1, ]
  } else { # take intercept (root) into consideration
    Xproj <- 0 + Xp[, c(root, projection), drop = FALSE]
    fit.gauss <- lm.fit(Xproj, Ypt)
    delta.gauss <- matrix(0, dim(Xp)[2] - 1, dim(Ypt)[2])
    coefs_gauss <- matrix(coef(fit.gauss), ncol = dim(Ypt)[2])
    E0.gauss <- coefs_gauss[1, ];# names(E0.gauss) <- NULL
    delta.gauss[projection, ] <- coefs_gauss[-1, ]
  }
  # If lm fails to find some coefficients
  if (anyNA(delta.gauss)) {
    warning("There were some NA in the lm fit for the Gauss Lasso. These were replaced with the values obtained from the Lasso. This is not the optimal solution.")
    delta.gauss[is.na(delta.gauss)] <- delta[is.na(delta.gauss)]
  }
  # Result
  # shifts.gauss <- shifts.matrix_to_list(t(delta.gauss));
  return(list(E0.gauss = E0.gauss, 
              delta.gauss = delta.gauss,
              residuals = residuals(fit.gauss)))
}

# compute_gauss_lasso.grplasso <- function (Yvec, Xkro, delta, root, group, p,
#                                           projection = which(as.vector(t(delta)) != 0)) {
#     Xproj <- 0 + Xkro[, c(which(group == root), projection), drop = FALSE]
#     fit.gauss <- lm(Yvec ~ Xproj - 1)
#     delta.gauss <- rep(0, dim(delta)[1] * dim(delta)[2])
#     coefs_gauss <- coef(fit.gauss)
#     E0.gauss <- coefs_gauss[1:p]; names(E0.gauss) <- NULL
#     delta.gauss[projection] <- coefs_gauss[-(1:p)]
#     delta.gauss <- t(matrix(delta.gauss, ncol = nrow(delta)))
#   # If lm fails to find some coefficients
#   if (anyNA(delta.gauss)) {
#     warning("There were some NA in the lm fit for the Gauss Lasso. These were replaced with the values obtained from the Lasso. This is not the optimal solution.")
#     delta.gauss[is.na(delta.gauss)] <- delta[is.na(delta.gauss)]
#   }
#   # Result
#   # shifts.gauss <- shifts.matrix_to_list(t(delta.gauss));
#   return(list(E0.gauss = E0.gauss, 
#               delta.gauss = delta.gauss,
#               residuals = residuals(fit.gauss)))
# }

compute_gauss_lasso.gglasso <- function (Yvec, Xkro, delta, root, group, p,
                                         projection = which(as.vector(t(delta)) != 0)) {
  # projection <- which(as.vector(t(delta)) != 0)
  if (is.null(root)) { # If no one is excluded, "real" intercept
    Xproj <- 0 + Xkro[, projection, drop = FALSE]
    fit.gauss <- lm(Yvec ~ Xproj)
    delta.gauss <- rep(0, dim(delta)[1] * dim(delta)[2])
    coefs_gauss <- coef(fit.gauss)
    E0.gauss <- coefs_gauss[1:p]; names(E0.gauss) <- NULL
    delta.gauss[projection] <- coefs_gauss[-(1:p)]
    delta.gauss <- t(matrix(delta.gauss, ncol = nrow(delta)))
  } else { # take intercept (root) into consideration
    Xproj <- 0 + Xkro[, c(which(group == root), projection), drop = FALSE]
    fit.gauss <- lm(Yvec ~ Xproj - 1)
    delta.gauss <- rep(0, dim(delta)[1] * dim(delta)[2])
    coefs_gauss <- coef(fit.gauss)
    E0.gauss <- coefs_gauss[1:p]; names(E0.gauss) <- NULL
    delta.gauss[projection] <- coefs_gauss[-(1:p)]
    delta.gauss <- t(matrix(delta.gauss, ncol = nrow(delta)))
  }
  # If lm fails to find some coefficients
  if (anyNA(delta.gauss)) {
    warning("There were some NA in the lm fit for the Gauss Lasso. These were replaced with the values obtained from the Lasso. This is not the optimal solution.")
    delta.gauss[is.na(delta.gauss)] <- delta[is.na(delta.gauss)]
  }
  # Result
  # shifts.gauss <- shifts.matrix_to_list(t(delta.gauss));
  return(list(E0.gauss = E0.gauss,
              delta.gauss = delta.gauss,
              residuals = residuals(fit.gauss)))
}

##
#' @title Initialization of the shifts using Lasso.
#'
#' @description
#' \code{init.EM.lasso} does the following regression :
#' ||Y_data-T.delta||_(Sigma_YY^(-1)) + lambda |delta|_1
#' using the function \code{glmnet::glmnet} of package \code{glmnet},
#' through function \code{lasso_regression_K_fixed}. 
#' T is the incidence matrix of the tree, and 
#' delta the vectorial representation of the shifts (see functions 
#' \code{incidence.matrix} and \code{shifts.list_to_vector} for further details).
#'
#' @details
#' A Cholesky decomposition of function Sigma_YY^(-1) is used.
#' lambda is chosen so that delta has the right number of non zero components.
#'
#' @param Y_data data at the tips.
#' @param times_shared (matrix) : times of shared ancestry, result of function 
#' \code{compute_times_ca}.
#' @param distances_phylo (matrix) : phylogenetic distance, result of function 
#' \code{compute_dist_phy}
#' @param nbr_of_shifts number of shifts used in the EM algorithm
#' 
#' @return params_init the list of initial parameters to be used, in the right
#'  format.
#'  
#' @keywords internal
#'
#18/06/14 - Initial release
#06/10/14 - Externalization of function lasso
##
init.EM.lasso <- function(phylo,
                          Y_data,
                          Y_data_imp = Y_data,
                          Y_data_vec_known = as.vector(Y_data),
                          process,
                          times_shared = compute_times_ca(phylo),
                          distances_phylo,
                          nbr_of_shifts,
                          K_lag_init = 0,
                          use_sigma = TRUE,
                          params_sigma = NULL,
                          variance.init = diag(1, p, p),
                          random.init = FALSE,
                          value.root.init = rep(0, p),
                          exp.root.init = rep(1, p),
                          var.root.init = diag(1, p, p),
                          edges.init = NULL,
                          values.init = matrix(0, p, length(edges.init)),
                          relativeTimes.init = NULL,
                          selection.strength.init = 1,
                          optimal.value.init = rep(0, p),
                          T_tree = incidence.matrix(phylo),
                          subtree.list = NULL,
                          miss = FALSE,
                          sBM_variance = FALSE,
                          stationary.root.init = FALSE,
                          impute_init_Rphylopars = FALSE,
                          masque_data,
                          independent = FALSE,
                          ...) {
  ntaxa <- length(phylo$tip.label)
  p <- nrow(Y_data)
  init.EM.default <- init.EM.default(process)
  ## If no shifts, hasty fix for initial value (if missing values)
  # TO DO : take variance matrix into account
  if (nbr_of_shifts == 0 && any(is.na(Y_data_imp)) && !impute_init_Rphylopars){
    E0 <- rowMeans(Y_data, na.rm = TRUE)
    params_init <- init.EM.default(Y_data = Y_data,
                                   value.root.init = E0, 
                                   exp.root.init = E0, 
                                   optimal.value.init = E0,
                                   edges.init = NULL, 
                                   values.init = NULL, 
                                   relativeTimes.init = NULL, 
                                   selection.strength.init = selection.strength.init, 
                                   random.init = random.init, 
                                   var.root.init = var.root.init,
                                   variance.init = variance.init,
                                   stationary.root.init = stationary.root.init,
                                   sBM_variance = sBM_variance,
                                   phylo = phylo, ...)
    return(params_init)
  }
  ## If missing data, impute them using Rphylopars
  if (impute_init_Rphylopars && any(is.na(Y_data_imp))){
    message("Imputing data for lasso initialization.")
    Y_data_imp <- impute.data.Rphylopars(phylo, Y_data, process, random.init)
  }
  ## Actualization of incidence matrix
  Tr <- T_tree
  if (independent){
    ac_tree <- lapply(selection.strength.init,
                      function(z) return(incidence_matrix_actualization_factors(tree = phylo, selection.strength = z, times_shared = times_shared)))
    Tr <- lapply(ac_tree, function(z) return(T_tree * z))
  } else {
    ac_tree <- incidence_matrix_actualization_factors(tree = phylo, 
                                                      selection.strength = selection.strength.init,
                                                      times_shared = times_shared)
    Tr <- T_tree * ac_tree
  }
  ## Choose the norm :
  if (use_sigma) {
    if (is.null(params_sigma)){
      # Initialize Sigma with default parameters
      params_sigma <- init.EM.default(Y_data = Y_data_imp,
                                      selection.strength.init = selection.strength.init, 
                                      random.init = random.init,
                                      stationary.root.init = stationary.root.init,
                                      var.root.init = var.root.init,
                                      variance.init = variance.init,
                                      value.root.init = value.root.init,
                                      exp.root.init = exp.root.init,
                                      edges.init = edges.init,
                                      values.init = values.init,
                                      relativeTimes.init = relativeTimes.init,
                                      optimal.value.init = optimal.value.init,
                                      nbr_of_shifts = nbr_of_shifts + K_lag_init,
                                      phylo = phylo,
                                      sBM_variance = sBM_variance)
    }
    if (!any(is.na(Y_data_imp)) && !independent){ # If there are no NA, do matrix computations
      # Choose process
      compute_tree_correlations_matrix  <- switch(process, 
                                                  BM = compute_tree_correlations_matrix.BM,
                                                  OU = compute_tree_correlations_matrix.scOU,
                                                  scOU = compute_tree_correlations_matrix.scOU)
      Fm <- compute_tree_correlations_matrix(times_shared = times_shared,
                                             distances_phylo = distances_phylo,
                                             params_old = params_sigma)
      Fm_YY <- extract.variance_covariance(Fm, what="YY",
                                           masque_data = c(rep(TRUE, ntaxa),
                                                           rep(FALSE, dim(Fm)[1] - ntaxa)))
      # Cholesky of tree-correlations
      Fm_chol <- t(chol(Fm_YY))
      Fm_chol_inv <- t(backsolve(t(Fm_chol), diag(ncol(Fm_chol)))) # Fm_YY_inv = t(Fm_chol_inv)%*%Fm_chol_inv
      # Cholesky of rate matrix
      R <- params_sigma$variance
      R_chol <- t(chol(R)) # R = R_chol %*% t(R_chol)
      R_chol_inv <- t(backsolve(t(R_chol), diag(ncol(R_chol)))) # R_YY_inv = t(R_chol_inv)%*%R_chol_inv
      # Transform Y_data and T
      Tr <- cbind(Tr, rep(1, dim(Tr)[1])) # Here we use hypothesis : beta_0 = mu if OU
      Tp <- Fm_chol_inv %*% Tr
      Yp <- R_chol_inv %*% Y_data_imp %*% t(Fm_chol_inv)
      fit <- try(lasso_regression_K_fixed.glmnet_multivariate(Yp = Yp, Xp = Tp,
                                          K = nbr_of_shifts,
                                          K_lag = K_lag_init,
                                          root = dim(Tr)[2]))
      chol_data <- TRUE
      } else { # If there are some NAs, use vectors.
        attr(params_sigma, "p_dim") <- p
        fun <- function(i){
          return(compute_mean_variance.simple(phylo = phylo,
                                              times_shared = times_shared,
                                              distances_phylo = distances_phylo,
                                              process = process,
                                              params_old = params_sigma,
                                              masque_data = masque_data))
        }
        if (independent){
          params_list <- split_params_independent(params_sigma)
          # Compute Sigma_YY^{-1/2} for each trait
          masque_data_matr <- matrix(masque_data,
                                     ncol = length(phylo$tip.label) + phylo$Nnode)
          moments <- vector(mode = "list", length = length(params_list))
          for (i in 1:length(params_list)){
            moments[[i]] <- compute_mean_variance.simple(phylo = phylo,
                                                         times_shared = times_shared,
                                                         distances_phylo = distances_phylo,
                                                         process = process,
                                                         params_old = params_list[[i]],
                                                         masque_data = masque_data_matr[i,])$Sigma_YY_chol_inv
          }
          Vp <- bdiag(moments)
          # Normalize data
          missbis <- as.vector(t(matrix(miss, nrow = p)))
          Yp <- Vp %*% as.vector(t(Y_data))[!missbis]
          # Regressor
          Tr <- lapply(Tr, function(z) return(cbind(z, rep(1, dim(z)[1]))))
          Xp <- bdiag(Tr)
          # Normalize predictor
          row_names_X <- 1:nrow(Xp)
          Xp <- Xp[!miss, ]
          Xp <- Vp %*% Xp
          # Reorder matrices
          corrdata <- as.vector(sapply(1:ntaxa,
                                       function(z) ((0:(p-1)) * ntaxa + z)))
          # corrdata <- corrdata[!miss]
          row_names_X <- row_names_X[!miss]
          corrdata <- as.vector(stats::na.exclude(match(corrdata, row_names_X)))
          corrreg <- as.vector(sapply(1:((nrow(phylo$edge) + 1)),
                                      function(z) ((0:(p-1)) * (nrow(phylo$edge) + 1) + z)))
          Xp <- Xp[corrdata, corrreg]
          # Data
          Ytemp <- rep(NA, ntaxa * p)
          Ytemp[!missbis] <- as.vector(Yp)
          corrdata <- as.vector(sapply(1:ntaxa,
                                       function(z) ((0:(p-1)) * ntaxa + z)))
          corrdata <- corrdata[!miss]
          Yp <- Ytemp[corrdata]
          # root
          root <- ncol(Tr[[1]])
        } else { # Case BM with missing values
          # Normalize data
          Vp <- compute_mean_variance.simple(phylo = phylo,
                                             times_shared = times_shared,
                                             distances_phylo = distances_phylo,
                                             process = process,
                                             params_old = params_sigma,
                                             masque_data = masque_data)$Sigma_YY_chol_inv
          Yp <- Vp %*% Y_data_vec_known
          # Regressor
          Tr <- cbind(Tr, rep(1, nrow(Tr)))
          Xp <- kronecker(Tr, diag(rep(1, p)))
          Xp <- Xp[masque_data[1:(p*ntaxa)], ]
          # Normalize predictor
          Xp <- Vp %*% Xp
          # Root
          root <- ncol(Tr)
        }
      # Fit
      group <- rep(1:root, each = p)
      fit <- try(lasso_regression_K_fixed.gglasso(Yvec = as.vector(Yp),
                                                  Xkro = as.matrix(Xp),
                                                  K = nbr_of_shifts,
                                                  K_lag = K_lag_init,
                                                  root = root,
                                                  group = group,
                                                  p_dim = p))
      chol_data <- FALSE
    }
  } else {
    # Return untransformed Y_data and T
    if (!any(is.na(Y_data_imp)) && !independent){
      Tr <- cbind(Tr, rep(1, dim(Tr)[1])) # Here we use hypothesis : beta_0 = mu if OU
      Yp <- Y_data_imp
      fit <- try(lasso_regression_K_fixed.glmnet_multivariate(Yp = Yp, Xp = Tp,
                                                              K = nbr_of_shifts,
                                                              K_lag = K_lag_init,
                                                              root = dim(Tr)[2]))
    } else {
      if (independent){
        Tr <- lapply(Tr, function(z) return(cbind(z, rep(1, dim(z)[1]))))
        Xp <- bdiag(Tr)
        # Reorder matrices
        corrdata <- as.vector(sapply(1:ntaxa,
                                     function(z) ((0:(p-1)) * ntaxa + z)))
        corrdata <- corrdata[!miss]
        corrreg <- as.vector(sapply(1:((nrow(phylo$edge) + 1)),
                                    function(z) ((0:(p-1)) * (nrow(phylo$edge) + 1) + z)))
        Xp <- Xp[corrdata, corrreg]
        # Data
        Ytemp <- rep(NA, ntaxa * p)
        missbis <- as.vector(t(matrix(miss, nrow = p)))
        Ytemp[!missbis] <- as.vector(t(Y_data))[!missbis]
        Yp <- Ytemp[corrdata]
        # root
        root <- ncol(Tr[[1]])
      } else { # Case BM with missing values
        Yp <- Y_data_vec_known
        # Regressor
        Xp <- kronecker(Tr, diag(rep(1, p)))
        Xp <- Xp[masque_data[1:(p*ntaxa)], ]
        # Root
        root <- ncol(Tr)
      }
      # Fit
      group <- rep(1:root, each = p)
      fit <- try(lasso_regression_K_fixed.gglasso(Yvec = as.vector(Yp),
                                                  Xkro = as.matrix(Xp),
                                                  K = nbr_of_shifts,
                                                  K_lag = K_lag_init,
                                                  group = group,
                                                  root = root,
                                                  p_dim = p))
    }
    chol_data <- FALSE
  }
  ## Fit
  if (inherits(fit, "try-error")) {
    warning("Lasso initialization fail : could not find a satisfying number of shifts. Proceeding to a default initialization.")
    return(init.EM.default(Y_data = Y_data,
                           value.root.init = value.root.init, 
                           exp.root.init = exp.root.init, 
                           optimal.value.init = optimal.value.init,
                           edges.init = edges.init, 
                           values.init = values.init, 
                           relativeTimes.init = relativeTimes.init, 
                           selection.strength.init = selection.strength.init, 
                           random.init = random.init, 
                           var.root.init = var.root.init,
                           variance.init = variance.init,
                           stationary.root.init = stationary.root.init,
                           nbr_of_shifts = nbr_of_shifts,
                           phylo = phylo,
                           sBM_variance = sBM_variance, ...))
  } else { 
    E0.gauss <- fit$E0.gauss
    delta.gauss <- t(fit$delta.gauss)
    if (chol_data){
      E0.gauss <- as.vector(R_chol %*% E0.gauss)
      delta.gauss <- R_chol %*% delta.gauss
    }
    shifts.gauss <- shifts.matrix_to_list(delta.gauss)
#     ## If OU, apply the correct factor to shifts
#     if (process == "OU" && !is.null(times_shared) && !is.null(shifts.gauss$values)){
#       parents <- phylo$edge[shifts.gauss$edges,1]
#       factors <- compute_actualization_factors(selection.strength = selection.strength.init, 
#                                                t_tree = t_tree, 
#                                                times_shared = times_shared, 
#                                                parents = parents)
#       shifts.gauss$values <- shifts.gauss$values/factors
#     }
    params_init <- init.EM.default(Y_data = Y_data,
                                   value.root.init = E0.gauss, 
                                   exp.root.init = E0.gauss, 
                                   optimal.value.init = E0.gauss,
                                   edges.init = shifts.gauss$edges, 
                                   values.init = shifts.gauss$values, 
                                   relativeTimes.init = shifts.gauss$relativeTimes, 
                                   selection.strength.init = selection.strength.init, 
                                   random.init = random.init, 
                                   var.root.init = var.root.init,
                                   variance.init = variance.init,
                                   stationary.root.init = stationary.root.init,
                                   sBM_variance = sBM_variance,
                                   phylo = phylo, ...)
    return(params_init)
  }
}

compute_actualization_factors <- function(selection.strength, 
                                                t_tree, 
                                                times_shared, 
                                                parents){
  factors <- exp(- selection.strength * (t_tree - diag(times_shared[parents, parents, drop = FALSE])))
  return(1-factors)
}

###############################################
## Robust initialization of alpha (and gamma)
###############################################

# init.alpha.BM <- function(method.init.alpha){
#   return(function(...) return(0))
# }
# 
# init.alpha.OU <- function(method.init.alpha){
#   return(switch(method.init.alpha, 
#                 default = init.alpha.default,
#                 estimation = init.alpha.estimation))
# }

init.alpha.gamma.BM <- function(method.init.alpha){
  fun <- function(random.root, init.var.root, ...){
    if (random.root){ 
      gamma_0 <- init.var.root
    } else {
      gamma_0 <- NULL
    }
    return(list(alpha_0 = NULL,
                gamma_0 = gamma_0))
  }
  return(fun)
}

init.alpha.gamma.OU <- function(method.init.alpha){
  return(switch(method.init.alpha, 
                default = init.alpha.gamma.default,
                estimation = init.alpha.gamma.estimation))
}

init.alpha.default <- function(init.selection.strength, known.selection.strength, alpha_known, ...){
  if (alpha_known) {
    return(matrix(known.selection.strength, 1, length(known.selection.strength)))
  } else {
    return(matrix(init.selection.strength, 1, length(init.selection.strength)))
  }
}


# init.alpha.estimation <- function(phylo, Y_data, nbr_of_shifts, distances_phylo, max_triplet_number, ...){
#   ## Initialize a vector with the group of each tip
#   tips_groups <- rep(0, length(phylo$tip.label))
#   names(tips_groups) <- phylo$tip.label
#   ## Initialize shifts by a lasso without sigma
#   if (nbr_of_shifts > 0) {
#     lasso <- init.EM.lasso(phylo=phylo, Y_data=Y_data, process="OU", nbr_of_shifts=nbr_of_shifts, use_sigma=FALSE)
#     ## Roeorder phylo and trace edges
#     phy <- reorder(phylo, order = "cladewise")
#     edges_shifts <- correspondanceEdges(edges=lasso$shifts$edges,from=phylo,to=phy)
#     ## Set groups of tips (one group = all the tips under a given shift)
#     Tr <- incidence.matrix(phy)
#     for (ed in order(edges_shifts)) { # Do it in order so that edges with higher numbers erase groups (edges closer from the tips)
#       ed_sh <- edges_shifts[ed]
#       tips_groups[phy$tip.label[Tr[,ed_sh]]] <- ed
#     }
#   } else {
#     edges_shifts <- NULL
#   }
#   ## For each group, take all the triplets of tips to estimate the covariance sigma_ij
#   cor_hat <- NULL # estimations from trilpets of pairs corelations
#   square_diff <- NULL # (Y_i-Y_j)^2
#   dists <- NULL # corresponding phylogenetic distances between pairs
#   hat_gam <- rep(0, length(edges_shifts)+1)
#   for (grp in 0:length(edges_shifts)) {
#     tips <- which(tips_groups==grp)
#     hat_gam[grp+1] <- var(Y_data[tips])
#     #     if (length(tips) > 2){
#     # #       ## Genrate all combinations
#     # #       Combinations <- combn(x=tips, m=3)
#     # #       Ncomb <- ncol(Combinations)
#     # #       ## If too many of them, sample from them
#     # #       if (Ncomb > max_triplet_number) {
#     # #         Combinations <- Combinations[,sample(Ncomb, max_triplet_number)]
#     # #       }
#     #       Ncomb <- choose(length(tips), 3)
#     #       ## If too many of them, sample from them (with replacement)
#     #       if (Ncomb > max_triplet_number) {
#     #         Combinations <- replicate(max_triplet_number, sample(tips, 3))
#     #       } else {
#     #         Combinations <- combn(x=tips, m=3)
#     #       }
#     #       temp_cor <- apply(X=Combinations, MARGIN=2, FUN=estimate_covariance_from_triplet, Y_data=Y_data, distances_phylo=distances_phylo)
#     #       temp_dist <- apply(X=Combinations, MARGIN=2, FUN=compute_dist_triplet, distances_phylo=distances_phylo)
#     #       temp_cor <- as.vector(temp_cor); temp_dist <- as.vector(temp_dist)
#     #       cor_hat <- c(cor_hat, temp_cor)
#     #       dists <- c(dists, temp_dist)
#     #     }
#     #   }
#     if (length(tips) > 1){
#       Z <- outer(Y_data[tips], Y_data[tips], function(x,y){x-y} )
#       square_diff <- c(square_diff, (Z[upper.tri(Z)])^2)
#       Z <- distances_phylo[tips,tips]
#       dists <- c(dists, Z[upper.tri(Z)])
#     }
#   }
#   #     ## Empirical mean of estimators for corelations estimated several times
#   #     un_dists <- unique(dists)
#   #     un_cor_hat <- rep(NA, length(un_dists))
#   #   }
#   #   for (dd in 1:length(un_dists)){
#   #     un_cor_hat[dd] <- mean(cor_hat[dists==un_dists[dd]])
#   #   }
#   #   pos_cor_hat <- cor_hat[cor_hat>=0]
#   #   dists_pos <- dists[cor_hat>=0]
#   ## Fit sigma_ij against gamma*exp(-alpha*d_ij)
#   #  fit <- nls(pos_cor_hat ~ gam*exp(-alpha*dists_pos), start=list(gam=mean(hat_gam), alpha=1))
#   #  fit <- nls(square_diff ~ gam*(1-exp(-alpha*dists)), start=list(gam=mean(hat_gam), alpha=1))
#   df <- data.frame(square_diff=square_diff, dists=dists)
#   fit.rob <- try(robustbase::nlrob(square_diff ~ gam*(1-exp(-alpha*dists)),
#                                    data = df,
#                                    start = list(gam = mean(hat_gam, na.rm=TRUE),
#                                                 alpha = 1)))
#   if (inherits(fit.rob, "try-error")) {
#     warning("Robust estimation of alpha failed")
#     return(init.alpha.default(...))
#   } else { 
#     return(unname(coef(fit.rob)["alpha"]))
#   }
# }

# compute_dist_triplet <- function(distances_phylo,v) {
#   phylo_dist <- c(distances_phylo[v[1],v[2]], distances_phylo[v[2],v[3]], distances_phylo[v[1],v[3]])
#   return(phylo_dist)
# }

# estimate_covariance_from_triplet <- function(Y_data, distances_phylo, v){
#   Y_i <- Y_data[v[1]]; Y_j <- Y_data[v[2]]; Y_k <- Y_data[v[3]];
#   hat_gamma <- var(Y_data[v])
#   cor <- hat_gamma - (1/9) * ((Y_i + Y_j + Y_k)^2 + (Y_i - Y_j)^2 + (Y_j - Y_k)^2 + (Y_i - Y_k)^2)
#   hat_sigmas <- c(Y_i*Y_j, Y_j*Y_k, Y_i*Y_k) + cor
#   return(hat_sigmas)
# }


init.alpha.gamma.default <- function(init.selection.strength, known.selection.strength,
                                     alpha_known, init.var.root, ...){
  if (!is.vector(init.var.root)){
    gamma_0 <- diag(init.var.root) 
  } else {
    gamma_0 <- init.var.root
  }
  gamma_0 <- matrix(gamma_0, 1, length(gamma_0))
  return(list(alpha_0 = init.alpha.default(init.selection.strength,
                                           known.selection.strength,
                                           alpha_known),
              gamma_0 = gamma_0))
}

##
#' @title Initialization the selection strength alpha using robust estimation
#'
#' @description
#' \code{init.alpha.estimation} fits (Y_i-Y_j)^2 ~ gamma^2(1-exp(-alpha*d_ij))
#' for all couples of tips (i,j) that have the same mean, i.e than are not
#' separated by a shift. Shifts are initialized thanks to a lasso
#' (function \code{init.EM.lasso}).
#'
#' @details
#' Function \code{robustbase::nlrob} is used for the robust fit.
#'
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param Y_data data at the tips.
#' @param nbr_of_shifts : number of shifts wanted
#' @param distances_phylo (matrix) : phylogenetic distance, result of function 
#' \code{compute_dist_phy}
#' @param nbr_of_shifts number of shifts used in the EM algorithm
#' 
#' @return params_init the list of initial parameters to be used, in the right
#'  format.
#'  
#' @keywords internal
#'
#10/07/14 - Initial release
##
init.alpha.gamma.estimation <- function(phylo, 
                                        Y_data, 
                                        nbr_of_shifts, 
                                        times_shared,
                                        distances_phylo, 
                                        T_tree,
                                        subtree.list,
                                        max_triplet_number, 
                                        alpha_known,
                                        method.init.alpha.estimation,
                                        tol, h_tree,
                                        miss,
                                        masque_data,
                                        independent, ...){
  ## Initialize a vector with the group of each tip
  tips_groups <- rep(0, length(phylo$tip.label))
  names(tips_groups) <- phylo$tip.label
  p <- nrow(Y_data)
  ## Initialize shifts by a lasso without sigma
  if (nbr_of_shifts > 0) {
    lasso <- init.EM.lasso(phylo = phylo,
                           Y_data = Y_data,
                           process = "OU",
                           nbr_of_shifts = nbr_of_shifts,
                           use_sigma = FALSE,
                           random.init = TRUE,
                           stationary.root.init = TRUE,
                           times_shared = times_shared,
                           distances_phylo = distances_phylo,
                           T_tree = T_tree,
                           subtree.list = subtree.list,
                           miss = miss,
                           impute_init_Rphylopars = FALSE,
                           masque_data = masque_data,
                           independent = independent,
                           selection.strength.init = rep(1, p))
    ## Roeorder phylo and trace edges
    phy <- reorder(phylo, order = "cladewise")
    edges_shifts <- correspondanceEdges(edges = lasso$shifts$edges,
                                        from = phylo, to = phy)
    ## Set groups of tips (one group = all the tips under a given shift)
    Tr <- incidence.matrix(phy)
    for (ed in order(edges_shifts)) { # Do it in order so that edges with higher numbers erase groups (edges closer from the tips)
      ed_sh <- edges_shifts[ed]
      tips_groups[phy$tip.label[Tr[, ed_sh]]] <- ed
    }
  } else {
    edges_shifts <- NULL
  }
  ## For each group, take all the triplets of tips to estimate the covariance sigma_ij
  cor_hat <- NULL # estimations from trilpets of pairs corelations
  square_diff <- vector("list", p) # (Y_i-Y_j)^2
  dists <- NULL # corresponding phylogenetic distances between pairs
  hat_gam <- matrix(NA, nrow = length(edges_shifts)+1, ncol = p)
  hat_gam_mad <- matrix(NA, nrow = length(edges_shifts)+1, ncol = p)
  for (grp in 0:length(edges_shifts)) {
    tips <- which(tips_groups==grp)
    if (length(tips) > 1){
      for (l in 1:p){
        hat_gam[grp+1, l] <- var(na.omit(Y_data[l, tips]))
        hat_gam_mad[grp+1, l] <- mad(Y_data[l, tips], na.rm = TRUE)^2
        Z <- outer(Y_data[l, tips], Y_data[l, tips],
                   function(x,y){x-y} )
        square_diff[[l]] <- c(square_diff[[l]], (Z[upper.tri(Z)])^2)
      }
      Z <- distances_phylo[tips,tips]
      dists <- c(dists, Z[upper.tri(Z)])
    }
  }
  ## Estimation of gamma
  gamma_0 <- matrix(NA, nrow = length(method.init.alpha.estimation) + 2, ncol = p)
  rownames(gamma_0) <- c("var", "mad", method.init.alpha.estimation)
  gamma_0["var", ] <- colMeans(hat_gam, na.rm = TRUE) # Simple variance
  gamma_0["mad", ] <- robustbase::colMedians(hat_gam_mad, na.rm = TRUE) # MAD
               
  ## Estimation of alpha
  # Supress couple "too far away"
  too_far <- (dists > h_tree)
  dists <- dists[!too_far]
  square_diff <- do.call(rbind, square_diff)
  square_diff <- square_diff[, !too_far, drop = FALSE]
  if (alpha_known) {
    return(list(alpha_0 = init.alpha.gamma.default(alpha_known, ...)$alpha_0,
                gamma_0 = gamma_0[c("var", "mad")]))
  } else {
    alpha_0 <- matrix(NA, nrow = length(method.init.alpha.estimation), ncol = p)
    rownames(alpha_0) <- method.init.alpha.estimation
    for (method in method.init.alpha.estimation){
      estimate.alpha  <- switch(method, 
                                regression = estimate.alpha.regression,
                                regression.MM = estimate.alpha.regression.MM,
                                median = estimate.alpha.median)
      
      for (l in 1:p){
        mask <- !is.na(square_diff[l, ])
        ag_0_try <- try(estimate.alpha(square_diff[l, mask],
                                       dists[mask],
                                       gamma_0["mad", l],
                                       tol, h_tree), silent = TRUE)
        
        if (inherits(ag_0_try, "try-error")) {
          message(paste0("Robust estimation of alpha by ", method, " failed."))
          alpha_0[method, l] <- NA # init.alpha.gamma.default(alpha_known, ...)$alpha_0
          gamma_0[method, l] <- NA
        } else {
          alpha_0[method, l] <- ag_0_try[["alpha_0"]]
          gamma_0[method, l] <- ag_0_try[["gamma_0"]]
        }
      }
    }
    return(list(alpha_0 = alpha_0, 
                gamma_0 = gamma_0))
  }
}

init.variance.BM.estimation <- function(phylo, 
                                        Y_data, 
                                        Y_data_imp,
                                        Y_data_vec_known,
                                        nbr_of_shifts, 
                                        times_shared,
                                        distances_phylo, 
                                        h_tree,
                                        random.root,
                                        T_tree,
                                        subtree.list,
                                        miss,
                                        impute_init_Rphylopars,
                                        masque_data,
                                        ...){
  ## Initialize a vector with the group of each tip
  tips_groups <- rep(0, length(phylo$tip.label))
  names(tips_groups) <- phylo$tip.label
  ## Initialize shifts by a lasso with default parameters
  if (nbr_of_shifts > 0) {
    lasso <- init.EM.lasso(phylo = phylo,
                           Y_data = Y_data,
                           Y_data_imp = Y_data_imp,
                           Y_data_vec_known = Y_data_vec_known,
                           process = "BM", 
                           times_shared = times_shared,
                           distances_phylo = distances_phylo, 
                           nbr_of_shifts = nbr_of_shifts, 
                           random.init = random.root,
                           use_sigma = TRUE,
                           T_tree = T_tree,
                           h_tree = h_tree,
                           subtree.list = subtree.list,
                           miss = miss,
                           impute_init_Rphylopars = impute_init_Rphylopars,
                           masque_data = masque_data,
                           ...)
    ## Roeorder phylo and trace edges
    phy <- reorder(phylo, order = "cladewise")
    edges_shifts <- correspondanceEdges(edges = lasso$shifts$edges,
                                        from = phylo,
                                        to = phy)
    ## Set groups of tips (one group = all the tips under a given shift)
    Tr <- incidence.matrix(phy)
    for (ed in order(edges_shifts)) { # Do it in order so that edges with higher numbers erase groups (edges closer from the tips)
      ed_sh <- edges_shifts[ed]
      tips_groups[phy$tip.label[Tr[,ed_sh]]] <- ed
    }
  } else {
    edges_shifts <- NULL
  }
  ## Recenter each group around its mean.
  centered_data <- NULL
  for (grp in 0:length(edges_shifts)) {
    tips <- which(tips_groups == grp)
    if (length(tips) > 1){
      centered_data <- cbind(centered_data,
                             Y_data[, tips, drop = F] - rowMeans(Y_data[, tips, drop = F],
                                                                 na.rm = TRUE))
    }
  }
  # centered_data <- centered_data[, colSums(is.na(centered_data)) < 1]
  R_0 <- try(suppressWarnings(robustbase::covMcd(t(centered_data),
                                                 nsamp = "deterministic")))
  # Robust did not fail
  if (!inherits(R_0, "try-error")) {
    Cov0 <- R_0$cov
    if (any(is.na(Cov0))) {
      warning("The initial estimation of the variance by covMcd gave some NAs. Replacing them by default value of 0.1.")
      Cov0[is.na(Cov0)] <- 0.1
    }
    Cov0 <- 1 / (h_tree + phylo$root.edge) * Cov0
    Cov0 <- try(suppressWarnings(nearPD(Cov0))) # Make sure the matrix is positive definite
    if (!inherits(Cov0, "try-error")) {
      return(Cov0$mat)
    }    
  }
  # Robust did fail
  warning("Robust intial estimation of the variance with covMcd failed or did not give a positive definite matrix. Doing a standard variance initialization with function cov.")
  Cov0 <- cov(t(centered_data), use = "na.or.complete")
  if (any(is.na(Cov0))) {
    warning("The initial estimation of the variance by covMcd gave some NAs. Replacing them by default value of 0.1.")
    Cov0[is.na(Cov0)] <- 0.1
  }
  Cov0 <- 1 / (h_tree + phylo$root.edge) * Cov0
  Cov0 <- try(suppressWarnings(nearPD(Cov0))) # Make sure the matrix is positive definite
  if (!inherits(Cov0, "try-error")) {
    return(Cov0$mat)
  }
  # Everything failed
  warning("Standard cov failed too. Returning default initialization with identity matrix for the covariance.")
  return(diag(rep(1, nrow(Y_data))))
}

## Regression on normalized half life to have the good tolerance.
estimate.alpha.regression <- function (square_diff, dists, gamma_0, tol, h_tree) {
  tol_t_half <- tol$normalized_half_life * h_tree
  df <- data.frame(square_diff = square_diff,
                   dists = dists)
  fit.rob <- robustbase::nlrob(square_diff ~ 2 * gam * (1 - exp(-log(2) / t_half * dists)),
                   data = df,
                   start = list(gam = gamma_0, t_half = log(2)),
                   tol = tol_t_half,
                   control = nls.control(tol = tol_t_half,
                                         maxiter = 100))
  gamma_0_M <- unname(coef(fit.rob)["gam"])
  if (gamma_0_M < 0){
    return(list(alpha_0 = NA,
                gamma_0 = NA))
  }
  return(list(alpha_0 = log(2) / unname(coef(fit.rob)["t_half"]),
              gamma_0 = gamma_0_M))
}

estimate.alpha.regression.MM <- function (square_diff, dists, gamma_0,
                                          tol, h_tree) {
  tol_t_half <- tol$normalized_half_life * h_tree
  df <- data.frame(square_diff = square_diff,
                   dists = dists)
  set.seed(18051220)
  low_bound = c(gam = unname(gamma_0)/5,
                t_half = 0.01 * h_tree)
  up_bound = c(gam = 5 * unname(gamma_0),
               t_half = 10 * h_tree)
  fit.rob <- robustbase::nlrob(square_diff ~ (2 * gam * (1 - exp(-log(2) / t_half * dists))),
                   data = df,
                   tol = tol_t_half,
                   lower = low_bound,
                   upper = up_bound,
                   method = "MM")
  gamma_0_M <- unname(coef(fit.rob)["gam"])
  if (gamma_0_M < 0){
    return(list(alpha_0 = NA,
                gamma_0 = NA))
  }
  return(list(alpha_0 = log(2) / unname(coef(fit.rob)["t_half"]),
              gamma_0 = gamma_0_M))
}

estimate.alpha.median <- function (square_diff, dists, gamma_0, ...) {
  delta <- 1 - square_diff / (2 * gamma_0)
  delta <- sapply(delta, function(x) max(0, x))
  rap <- -log(delta) / dists
  med <- median(rap)
  if (is.infinite(med)){
    rap[is.infinite(rap)] <- NA
    med <- max(rap, na.rm = TRUE)
  }
  return(list(alpha_0 = med,
              gamma_0 = gamma_0))
}


##
#' @title Initial imputation of missing data for lasso
#'
#' @description
#' \code{impute.data.Rphylopars} uses function \code{phylopars} from package \code{Rphylopars}
#' to impute missing data.
#' 
#' @details 
#' This function assume that there are no shifts on the tree. It is only a first approximation
#' for initialization purposes.
#'
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param Y_data data at the tips.
#' @param process the stochastic process
#' @param random.init whether root is random or fixed.
#' 
#' @return Y_data_imp the imputed data using Rphylopars
#' 
#' @keywords internal
##

impute.data.Rphylopars <- function(phylo, Y_data, process, random.init){
  if (!requireNamespace("Rphylopars", quietly = TRUE)) {
    stop("Rphylopars is needed when option 'impute_init_Rphylopars' is set to TRUE. Please install it.",
         call. = FALSE)
  } else {
    message("Using Rphylopars for initial data imputation.")
  }
  process_Rphylopars <- choose_process_Rphyopars(process, random.init)
  # library(Rphylopars)
  trait_data <- as.data.frame(t(Y_data))
  trait_data <- cbind(phylo$tip.label, trait_data)
  colnames(trait_data)[1] <- "species"
#   trait_data <- as.data.frame(t(Y_data))
#   trait_data[ , "species"] <- phylo$tip.label
  fit_phylopars <- Rphylopars::phylopars(trait_data,
                                         phylo,
                                         model = process_Rphylopars,
                                         pheno_error = FALSE,
                                         phylo_correlated = TRUE,
                                         pheno_correlated = FALSE,
                                         REML = TRUE,
                                         # optim_limit = 50,
                                         # BM_first = TRUE,
                                         usezscores = TRUE)
#   data_phylopars <- try(phylopars.predict(fit_phylopars, nodes = NULL))
#   if (inherits(data_phylopars, "try-error")) { # If fails, replace with mean of the trait
#     warning("The RPhyloPars imputation failed. Taking the mean of each trait for missing data for initialization.")
#     Y_data_imp <- Y_data
#     for (j in 1:(dim(Y_data_imp)[1])){
#       Y_data_imp[j, is.na(Y_data_imp[j, ])] <- mean(Y_data_imp[j, ], na.rm = TRUE)
#     }
#   } else {
#     Y_data_imp <- t(unname(as.matrix(data_phylopars$predicted)))
#   }
  Y_data_imp <- t(fit_phylopars$anc_recon[1:ncol(Y_data), ])
  return(Y_data_imp)
}

choose_process_Rphyopars <- function(process, random.init){
  if (process == "BM"){
    return(process)
  } else if (process == "OU"){
    stop("Rphylopars imputation only works for scalar OU. Could not do the Lasso initialization.")
  } else if (process == "scOU"){
    if (random.init){
      return("OUrandomRoot")
    } else {
      return("OUfixedRoot")
    }
  }
}