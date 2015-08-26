# {E step}
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
## Here are the functions to compute the E step.
## Dependencies : generic_functions.R
## : simulate.R
## : shifts_manipulations.R
###############################################################################

##
#' @title E step
#'
#' @description
#' \code{compute_E.simple} computes the E step in the simple case where the invert matrix Sigma_YY_inv is given
#'
#' @details
#' This function takes parameters sim, Sigma and Sigma_YY_inv from  
#' \code{compute_mean_variance.simple}. It uses functions 
#' \code{extract.variance_covariance}, \code{extract.covariance_parents}, and
#'  \code{extract.simulate} to extract the needed quantities from these objects.
#'
#' @param phylo Input tree.
#' @param Y_data : vector indicating the data at the tips
#' @param sim (list) : result of function \code{simulate}
#' @param Sigma : variance-covariance matrix, result of function \code{compute_variance_covariance}
#' @param Sigma_YY_inv : invert of the variance-covariance matrix of the data
#' 
#' @return conditional_law_X (list) : list of conditionnal statistics :
#'                   "expectation" : vector of length ntaxa+nNodes, with ntaxa
#' fisrt values set to Y_data (tips), and from ntaxa+1 of conditional expectation
#' of the nodes conditionned to the tips E[Z_j|Y]
#'                   "variances" : vector of length ntaxa+nNodes with ntaxa 0
#' (tips) and conditional variance of the nodes conditionned to the tips 
#' Var[Z_j|Y]
#'                  "covariances" : vector of length ntaxa+nNodes with ntaxa 0
#' (tips) and conditional covariance of the nodes and their parents conditionned 
#' to the tips Cov[Z_j,Z_pa(j)|Y], with NA for the root.
#'                   "optimal.values" : vector of length ntaxa+nNodes of optimal
#' values beta(t_j)
#' 
#' 21/05/14 - Initial release (non fonctionnal)
#' 22/05/14 - Minimal "working" release
#' 02/06/14 - Add case OU
#' 29/09/14 - Reshape to externalize computation of Sigma_YY_inv
##
compute_E.simple <- function (phylo, Y_data, sim, Sigma, Sigma_YY_inv) {
  ## Initialization
  ntaxa <- length(phylo$tip.label)
  conditional_law_X <- list(expectations = rep(NA,ntaxa-1), 
                            variances = rep(NA,ntaxa-1), 
                            covariances = rep(NA,ntaxa-2))
  ## Mean
  m_Y <- extract.simulate(sim, where="tips", what="expectations")
  m_Z <- extract.simulate(sim, where="nodes", what="expectations")
  conditional_law_X$optimal.values <- c(extract.simulate(sim,
                                                         where="tips",
                                                         what="optimal.values"),
                                        extract.simulate(sim,
                                                         where="nodes",
                                                         what="optimal.values")) # NULL if BM
  ## Variance Covariance
  Sigma_YZ <- extract.variance_covariance(Sigma, what="YZ")
  Sigma_ZZ <- extract.variance_covariance(Sigma, what="ZZ")
  temp <- Sigma_YZ %*% Sigma_YY_inv
  conditional_law_X$expectations <- c(Y_data, m_Z + temp%*%(Y_data-m_Y))
  conditional_variance_covariance <- Sigma_ZZ - temp%*%t(Sigma_YZ)
  conditional_law_X$variances <- c(rep(0, ntaxa), diag(conditional_variance_covariance))
  conditional_law_X$covariances <- c(rep(0, ntaxa), extract.covariance_parents(phylo, conditional_variance_covariance))
  return(conditional_law_X)
}


##
#' @title Extract sub-matrices of variance.
#'
#' @description
#' \code{extract.variance_covariance} return the adequate sub-matrix.
#' 
#' @param struct structural matrix of size (ntaxa+nNode)*p, result 
#' of function \code{compute_variance_covariance}
#' @param what: sub-matrix to be extracted:
#'                "YY" : sub-matrix of tips (p*ntaxa first lines and columns)
#'                "YZ" : sub matrix tips x nodes (p*nNodes last rows and p*ntaxa first columns)
#'                "ZZ" : sub matrix of nodes (p*nNodes last rows and columns)
#' 
#' @return sub-matrix of variance covariance.
#' 
##
extract.variance_covariance <- function(struct, what=c("YY","YZ","ZZ")){
  ntaxa <- attr(struct, "ntaxa")
  p <- attr(struct, "p")
  if (what=="YY") {
    return(struct[1:(p * ntaxa), 1:(p * ntaxa)])
  } else if (what=="YZ") {
    return(struct[(p * ntaxa + 1):(dim(struct)[1]), 1:(p * ntaxa)])
  } else if (what=="ZZ") {
    return(struct[(p * ntaxa + 1):(dim(struct)[1]),(p * ntaxa + 1):(dim(struct)[2])])
  }
}

##
# extract.covariance_parents (phylo, struct)
# PARAMETERS:
#            @phylo (tree)
#            @struct (matrix) structural matrix of size ntaxa+nNode, result of function compute_times_ca, compute_dist_phy or compute_variance_covariance
# RETURNS:
#            (vector) : for every node i, entry i-ntaxa of the vector is (i,pa(i)). For the root (ntaxa+1), entry 1 is NA
# DEPENDENCIES:
#            none
# PURPOSE:
#            Extract covariances needed
# NOTES:
#            none
# REVISIONS:
#            22/05/14 - Initial release
##
extract.covariance_parents<- function(phylo, struct){
  ntaxa <- length(phylo$tip.label)
  cov <- rep(NA,dim(phylo$edge)[1] - ntaxa + 1)
  for (i in (ntaxa+2):(dim(phylo$edge)[1]+1)) {
    pa <- getAncestor(phylo,i)
    cov[i-ntaxa] <- struct[i-ntaxa,pa-ntaxa]
  }
  return(cov)
}

##
#' @title Complete variance covariance matrix for BM
#'
#' @description
#' \code{compute_variance_covariance.BM} computes the (n+m)*p squared variance covariance
#' matrix of vec(X).
#'
#'  @param times_shared times of shared ancestry of all nodes and tips, result of function
#'  \code{compute_times_ca}
#'  @param params_old (list) : old parameters to be used in the E step
#' 
#' @return matrix of variance covariance for the BM
#' 
##
compute_variance_covariance.BM <- function(times_shared, params_old, ...) {
  p <- nrow(params_old$shifts$values)
  if (p == 1){ # Dimension 1 (next would also work, but slightly faster)
    J <- matrix(1, nrow = dim(times_shared)[1], ncol = dim(times_shared)[2])
    varr <- params_old$variance * times_shared
    if (params_old$root.state$random) {
      varr <- varr + params_old$root.state$var.root * J
    }
  } else {
    varr <- kronecker(times_shared, params_old$variance)
    if (params_old$root.state$random) {
      varr <- varr + kronecker(diag(1, dim(times_shared)), 
                               params_old$root.state$var.root)
    }
  }
  attr(varr, "p") <- p
  attr(varr, "ntaxa") <- attr(params_old, "ntaxa")
  return(varr)
}

compute_variance_covariance.OU <- function(times_shared, distances_phylo, params_old, ...) {
  alpha <- params_old$selection.strength
  sigma2 <- params_old$variance
  Var <- sigma2/(2*alpha) * (1 - exp(- 2 * alpha * times_shared)) * exp(- alpha * distances_phylo)
  if (!params_old$root.state$random) {
    return(Var)
  } else if (params_old$root.state$stationary.root) {
    return(sigma2/(2*alpha) * exp(- alpha * distances_phylo))
  } else {
    times_nodes <- list(diag(times_shared))
    sum_times <- do.call('rbind',rep(times_nodes,length(diag(times_shared)))) + do.call('cbind',rep(times_nodes,length(diag(times_shared))))
    gamma2 <- params_old$root.state$var.root
    return( gamma2 * exp(- alpha * sum_times) + Var)
  }
}

##
#' @title Compute moments of params_old
#'
#' @description
#' \code{compute_mean_variance.simple} computes the quantities needed to compute
#' mean and variance matrix with parameters params_old.
#'
#' @details
#' This function is used by functions \code{compute_E.simple} and 
#' \code{compute_log_likelihood.simple}.
#'
#' @param phylo Input tree.
#' @param times_shared (matrix) : times of shared ancestry, result of function 
#' \code{compute_times_ca}
#' @param distances_phylo (matrix) : phylogenetics distance, result of function 
#' \code{compute_dist_phy}
#' @param process a two letter string indicating the process to consider
#' @param params_old a list of parameters to be used in the computations
#' 
#' @return sim (list) : result of funtion \code{simulate} with the appropriate
#'  parameters
#' @return Sigma matrix of variance covariance, result of function 
#' \code{compute_variance_covariance}
#' @return Sigma_YY_inv inverse of vairance matrix of the data
#' 29/09/14 - Initial release
##
compute_mean_variance.simple <- function (phylo,
                                          times_shared,
                                          distances_phylo,
                                          process=c("BM","OU"),
                                          params_old, ...) {
  ## Choose process 
  process  <- match.arg(process)
  compute_variance_covariance  <- switch(process, 
                                         BM = compute_variance_covariance.BM,
                                         OU = compute_variance_covariance.OU)
  ## Mean
  sim <- simulate(phylo = phylo, 
                  process = process,
                  p = attr(params_old, "p"),
                  root.state = params_old$root.state, 
                  shifts = params_old$shifts, 
                  variance = params_old$variance, 
                  optimal.value = params_old$optimal.value, 
                  selection.strength = params_old$selection.strength)
  ## Variance Covariance
  Sigma <- compute_variance_covariance(times_shared = times_shared, 
                                       distances_phylo = distances_phylo,
                                       params_old = params_old)
  Sigma_YY <- extract.variance_covariance(Sigma, what="YY")
  Sigma_YY_inv <- solve(Sigma_YY)
  return(list(sim = sim, Sigma = Sigma, Sigma_YY_inv = Sigma_YY_inv))
}

#####################################################################
## Functions to compute the log likelihood
#####################################################################

##
#' @title Squared Mahalanobis Distance
#'
#' @description
#' \code{compute_mahalanobis_distance.simple} computes the squared mahalanobis distance 
#' between the data and mean at tips  of the data in the simple case where the inverse
#' of the variance matrix is given.
#'
#' @details
#' This function takes parameters sim and Sigma_YY_inv from  
#' \code{compute_mean_variance.simple}. It uses function \code{extract.simulate}
#' to extract the needed quantities from these objects.
#'
#' @param phylo Input tree.
#' @param Y_data : vector indicating the data at the tips.
#' @param sim (list) : result of function \code{simulate}.
#' @param Sigma_YY_inv : invert of the variance-covariance matrix of the data.
#' 
#' @return squared Mahalanobis distance between data and mean at the tips.
##
compute_mahalanobis_distance.simple <- function(phylo, Y_data, sim, Sigma_YY_inv){
  ntaxa <- length(phylo$tip.label)
  m_Y <- extract.simulate(sim, where="tips", what="expectations")
  m_Y <- as.vector(m_Y)
  Y_data <- as.vector(Y_data)
  MD <- t(Y_data - m_Y)%*%Sigma_YY_inv%*%(Y_data - m_Y)
  return(MD)
}


##
#' @title Log Likelihood
#'
#' @description
#' \code{compute_log_likelihood.simple} computes the log-likelihood of the data 
#' in the simple case where the inverse of the variance matrix is given.
#'
#' @details
#' This function takes parameters sim, Sigma and Sigma_YY_inv from  
#' \code{compute_mean_variance.simple}. It uses functions 
#' \code{extract.variance_covariance}, \code{extract.covariance_parents}, and
#'  \code{extract.simulate} to extract the needed quantities from these objects.
#'
#' @param phylo Input tree.
#' @param Y_data : vector indicating the data at the tips
#' @param sim (list) : result of function \code{simulate}
#' @param Sigma : variance-covariance matrix, result of function \code{compute_variance_covariance}
#' @param Sigma_YY_inv : invert of the variance-covariance matrix of the data
#' 
#' @return log likelihood of the data
#' 
#' 29/09/14 - Initial release
##
compute_log_likelihood.simple <- function(phylo, Y_data, sim, Sigma, Sigma_YY_inv){
  ntaxa <- length(phylo$tip.label)
  Sigma_YY <- extract.variance_covariance(Sigma, what="YY")
  logdetSigma_YY <- determinant(Sigma_YY, logarithm = TRUE)$modulus
  m_Y <- extract.simulate(sim, where="tips", what="expectations")
  LL <- ntaxa * log(2*pi) + logdetSigma_YY
  m_Y_vec <- as.vector(m_Y)
  Y_data_vec <- as.vector(Y_data)
  LL <- LL + t(Y_data_vec - m_Y_vec) %*% Sigma_YY_inv %*% (Y_data_vec - m_Y_vec)
  return(-LL/2)
}

compute_entropy.simple <- function(Sigma, Sigma_YY_inv){
  Sigma_YZ <- extract.variance_covariance(Sigma, what="YZ")
  Sigma_ZZ <- extract.variance_covariance(Sigma, what="ZZ")
  conditional_variance_covariance <- Sigma_ZZ - Sigma_YZ%*%Sigma_YY_inv%*%t(Sigma_YZ)
  N <- dim(Sigma_ZZ)[1]
  logdet_conditional_variance_covariance <- determinant(conditional_variance_covariance, logarithm = TRUE)$modulus
  return(N/2*log(2*pi*exp(1)) + 1/2*logdet_conditional_variance_covariance)
}

compute_log_likelihood_with_entropy.simple <- function(CLL, H){
  return(CLL + H)
}