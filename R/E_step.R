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
#'                   "expectation" : matrix of size p x (ntaxa+nNodes), with ntaxa
#' fisrt columns set to Y_data (tips), and from ntaxa+1 to conditional expectation
#' of the nodes conditionned to the tips E[Z_j|Y]
#'                   "variances" : array of size p x p x (ntaxa+nNodes) with ntaxa first 
#' matrices of zeros (tips) and conditional variance of the nodes conditionned to the tips 
#' Var[Z_j|Y]
#'                  "covariances" : array of size p x p x (ntaxa+nNodes) with ntaxa first 
#' matrices of zeros (tips) and conditional covariance of the nodes and their parents
#' conditionned to the tips Cov[Z_j,Z_pa(j)|Y], with NA for the root.
#'                   "optimal.values" : matrix of size p x ntaxa+nNodes of optimal
#' values beta(t_j)
#' 
##
compute_E.simple <- function (phylo, Y_data_vec, sim, Sigma, Sigma_YY_chol_inv,
                              missing, masque_data) {
  ## Initialization
  ntaxa <- length(phylo$tip.label)
  nNodes <- phylo$Nnode
  p <- dim(sim)[1]
  nMiss <- sum(missing)
  # index_missing <- (nNodes * p + 1):(nNodes * p + nMiss)
  index_missing <- c(rep(TRUE, nMiss), rep(FALSE, nNodes * p))
  conditional_law_X <- list(expectations = matrix(NA, p, ntaxa + nNodes), 
                            variances = array(NA, c(p, p, ntaxa + nNodes)), 
                            covariances = matrix(NA, c(p, p, ntaxa + nNodes)))
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
  Sigma_YZ <- extract.variance_covariance(Sigma, what="YZ", masque_data)
  Sigma_ZZ <- extract.variance_covariance(Sigma, what="ZZ", masque_data)
#  temp <- Sigma_YZ %*% Sigma_YY_inv
  temp <- Sigma_YZ %*% Sigma_YY_chol_inv
  # Y_data_vec <- as.vector(Y_data)
  m_Y_vec <- as.vector(m_Y)[!missing]
  m_Z_vec <- c(as.vector(m_Y)[missing], as.vector(m_Z))
  # Conditionnal expectation of unkonwn values
  exp_Z_vec <- m_Z_vec + temp %*% crossprod(Sigma_YY_chol_inv,
                                            (Y_data_vec - m_Y_vec))
  expcond <- rep(NA, (nNodes + ntaxa) * p)
  # Data
  expcond[masque_data] <- Y_data_vec
  # Missing Data and Nodes
  expcond[!masque_data] <- as.vector(exp_Z_vec)
  conditional_law_X$expectations <- matrix(expcond, nrow = p)
#  conditional_variance_covariance <- Sigma_ZZ - temp %*% t(Sigma_YZ)
  ## Variances
  conditional_variance_covariance <- Sigma_ZZ - tcrossprod(temp)
  attr(conditional_variance_covariance, "p_dim") <- p
  conditional_variance_covariance_nodes <- conditional_variance_covariance[!index_missing, !index_missing]
  attr(conditional_variance_covariance_nodes, "p_dim") <- p
  # Data tips
  var_tips <- array(0, c(p, p, ntaxa))
  cov_tips <- array(0, c(p, p, ntaxa))
  if (nMiss > 0){
    conditional_variance_covariance_tips <- conditional_variance_covariance[index_missing, index_missing, drop = FALSE]
    conditional_variance_covariance_tips_nodes <- conditional_variance_covariance[!index_missing, index_missing, drop = FALSE]
    
    missing_tips <- (which(missing) - 1) %/% p + 1
    missing_tips_uniques <- unique(missing_tips)
    par_missing_tips <- getAncestors(phylo, missing_tips_uniques)
    par_missing_tips <- par_missing_tips - ntaxa
    par_missing_tips <- sapply(par_missing_tips,
                               function(z) (p * (z - 1) + 1):(p * z))
    par_missing_tips <- matrix(par_missing_tips, nrow = p)
    missing_chars <- (which(missing) - 1) %% p + 1
    grpes_missing <- sapply(1:ntaxa, function(z) missing_tips == z)
    for (i in 1:length(missing_tips_uniques)){
      tip <- missing_tips_uniques[i]
      tipgrp <- grpes_missing[, tip]
      var_tips[missing_chars[tipgrp], missing_chars[tipgrp], tip] <- as.matrix(conditional_variance_covariance_tips[tipgrp, tipgrp])
      cov_tips[, missing_chars[tipgrp], tip] <- as.matrix(conditional_variance_covariance_tips_nodes[par_missing_tips[, i], tipgrp])
    }
  }
  # Nodes - varariances
  var_nodes <- extract.variance_nodes(phylo,
                                      conditional_variance_covariance_nodes)
  conditional_law_X$variances <- array(c(var_tips,
                                         var_nodes), c(p, p, ntaxa + nNodes))
  # Nodes - covariances
  cov_nodes <- extract.covariance_parents(phylo,
                                          conditional_variance_covariance_nodes)
  conditional_law_X$covariances <- array(c(cov_tips,
                                           cov_nodes), c(p, p, ntaxa + nNodes))
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
#' @param missing; missing values of Y_data
#' 
#' @return sub-matrix of variance covariance.
#' 
##
extract.variance_covariance <- function(struct, what=c("YY","YZ","ZZ"),
                                        masque_data = c(rep(TRUE, attr(struct, "ntaxa") * attr(struct, "p_dim")), rep(FALSE, (dim(struct)[1] - attr(struct, "ntaxa")) * attr(struct, "p_dim")))){
  ntaxa <- attr(struct, "ntaxa")
  p <- attr(struct, "p_dim")
  if (what=="YY") {
    return(struct[masque_data, masque_data])
  } else if (what=="YZ") {
    return(struct[!masque_data, masque_data])
  } else if (what=="ZZ") {
    return(struct[!masque_data, !masque_data])
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
extract.covariance_parents <- function(phylo, struct){
  ntaxa <- length(phylo$tip.label)
  p <- attr(struct, "p_dim")
  m <- dim(phylo$edge)[1] - ntaxa + 1
  # cov1 <- array(NA, c(p, p, m))
  daughters <- phylo$edge[, 2]
  parents <- phylo$edge[, 1]
  masque <- daughters > ntaxa + 1
  daughters <- daughters[masque]
  parents <- parents[masque]
  range_node <- function(node){
    tmp <- sapply(node, function(z) ((z - ntaxa - 1) * p + 1):((z - ntaxa) * p))
    return(as.vector(tmp))
  }
  arr <- struct[range_node(daughters), range_node(parents)] # + struct[range_node(parents), range_node(daughters)]
  masque <- rep(c(rep(rep(c(T, F),
                          c(p, p * (m - 2))),
                      p - 1),
                  rep(c(T, F),
                      c(p, p * (m - 1)))),
                m)[1:(p*m)^2]
  arr <- array(as.matrix(arr)[masque], c(p, p, m - 1))
  arr[, , daughters - ntaxa - 1] <- arr
  arr <- array(c(rep(NA, p*p), arr), c(p, p, m))
#   for (i in (ntaxa + 2):(dim(phylo$edge)[1] + 1)) {
#     pa <- getAncestor(phylo,i)
#     range_i <- ((i - ntaxa - 1) * p + 1):((i - ntaxa) * p)
#     range_pa <- ((pa - ntaxa - 1) * p + 1):((pa - ntaxa) * p)
#     cov1[1:p, 1:p, i - ntaxa] <- as.matrix(struct[range_i, range_pa] + struct[range_pa, range_i])
#   }
  return(arr)
}

extract.variance_nodes <- function(phylo, struct){
  ntaxa <- length(phylo$tip.label)
  p <- attr(struct, "p_dim")
  m <- dim(phylo$edge)[1] - ntaxa + 1
#   fun <- function(i){
#     range <- ((i - ntaxa - 1) * p + 1):((i - ntaxa) * p)
#     return(struct[range, range])
#   }
#   aa <- sapply((ntaxa + 1):(dim(phylo$edge)[1] + 1), fun, simplify = "array")
#   return(aa)
  ##
#   varr1 <- array(NA, c(p, p, m))
#   for (i in (ntaxa + 1):(dim(phylo$edge)[1] + 1)) {
#     range_i <- ((i - ntaxa - 1) * p + 1):((i - ntaxa) * p)
#     varr1[1:p, 1:p, i - ntaxa] <- as.matrix(struct[range_i, range_i])
#   }
#   return(varr)
  ##
  masque <- rep(c(rep(rep(c(T, F),
                          c(p, p * (m-1))),
                      p - 1),
                  rep(c(T, F),
                      c(p, p * m))),
                m)[1:(p*m)^2]
  varr <- array(as.matrix(struct)[masque], c(p, p, m))
  return(varr)
}

##
#' @title Get variance matrix of a node
#'
#' @description
#' \code{get_variance_node} returns the conditional variance of a node, or the conditional
#' corariance of a node and its parent.
#' 
#' @param vars matrix of size p x p*(ntaxa+nNodes) result of funtion \code{compute_E.simple},
#' entry "variances" or "covariances".
#' @param node for which to extract the matrix.
#' 
#' @return sub-matrix of variance for the node.
#' 
##

get_variance_node <- function(node, vars){
  return(vars[, , node])
#   p <- nrow(vars)
#   nN <- ncol(vars) / p
#   range <- ((node - 1) * p + 1):(node * p)
#   return(vars[1:p, range, drop = F])
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
    varr <- as.vector(params_old$variance) * times_shared
    if (params_old$root.state$random) {
      J <- matrix(1, nrow = dim(times_shared)[1], ncol = dim(times_shared)[2])
      varr <- varr + as.vector(params_old$root.state$var.root) * J
    }
  } else {
    varr <- kronecker(times_shared, params_old$variance)
    if (params_old$root.state$random) {
      varr <- varr + kronecker(as(diag(1, dim(times_shared)), "symmetricMatrix"), 
                               params_old$root.state$var.root)
    }
  }
  attr(varr, "p_dim") <- p
  attr(varr, "ntaxa") <- attr(params_old, "ntaxa")
  return(varr)
}

compute_variance_covariance.scOU <- function(times_shared, distances_phylo, params_old, ...) {
  p <- nrow(params_old$shifts$values)
  if (is.null(p)) p <- 1
  alpha <- as.vector(params_old$selection.strength)
  sigma2 <- params_old$variance
  Var <-  kronecker((1 - exp(- 2 * alpha * times_shared)) * exp(- alpha * distances_phylo),
                                                               sigma2 / (2 * alpha))
  if (!params_old$root.state$random) {
    varr <- Var
  } else if (params_old$root.state$stationary.root) {
    varr <- kronecker(exp(- alpha * distances_phylo),
                      sigma2 / (2 * alpha))
  } else {
    times_nodes <- list(diag(times_shared))
    sum_times <- do.call('rbind',
                         rep(times_nodes, length(diag(times_shared))))
    sum_times <- sum_times + do.call('cbind',
                                     rep(times_nodes, length(diag(times_shared))))
    gamma2 <- params_old$root.state$var.root
    varr <- kronecker(exp(- alpha * sum_times), gamma2) + Var
  }
  attr(varr, "p_dim") <- p
  attr(varr, "ntaxa") <- attr(params_old, "ntaxa")
  return(varr)
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
#' #@return Sigma_YY_inv inverse of vairance matrix of the data
#' @return Sigma_YY_chol_inv invert of cholesky matrix of Sigma_YY:
#'  (Sigma_YY)^(-1) = tcrossprod(Sigma_YY_chol_inv)
#' 29/09/14 - Initial release
##
compute_mean_variance.simple <- function (phylo,
                                          times_shared,
                                          distances_phylo,
                                          process=c("BM", "OU", "rBM", "scOU"),
                                          params_old,
                                          masque_data, ...) {
  ## Choose process 
  process  <- match.arg(process)
  compute_variance_covariance  <- switch(process, 
                                         BM = compute_variance_covariance.BM,
                                         OU = compute_variance_covariance.scOU,
                                         scOU = compute_variance_covariance.scOU)
  ## Mean
  sim <- simulate(phylo = phylo, 
                  process = process,
                  p = attr(params_old, "p_dim"),
                  root.state = params_old$root.state, 
                  shifts = params_old$shifts, 
                  variance = params_old$variance, 
                  optimal.value = params_old$optimal.value, 
                  selection.strength = params_old$selection.strength)
  ## Variance Covariance
  Sigma <- compute_variance_covariance(times_shared = times_shared, 
                                       distances_phylo = distances_phylo,
                                       params_old = params_old)
  Sigma_YY <- extract.variance_covariance(Sigma, what="YY", masque_data = masque_data)
  Sigma_YY_chol <- chol(Sigma_YY)
  Sigma_YY_chol_inv <- backsolve(Sigma_YY_chol, diag(ncol(Sigma_YY_chol)))
  #Sigma_YY_inv <- chol2inv(Sigma_YY_chol)
  return(list(sim = sim, Sigma = Sigma, Sigma_YY_chol_inv = Sigma_YY_chol_inv))
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
compute_mahalanobis_distance.simple <- function(phylo, Y_data_vec, sim,
                                                Sigma_YY_chol_inv, missing){
  ntaxa <- length(phylo$tip.label)
  m_Y <- extract.simulate(sim, where="tips", what="expectations")
  m_Y <- as.vector(m_Y)[!missing]
#  MD <- t(Y_data - m_Y)%*%Sigma_YY_inv%*%(Y_data - m_Y)
  MD <- tcrossprod(t(Y_data_vec - m_Y) %*% Sigma_YY_chol_inv)
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
compute_log_likelihood.simple <- function(phylo, Y_data_vec, sim,
                                          Sigma, Sigma_YY_chol_inv,
                                          missing, masque_data){
  ntaxa <- length(phylo$tip.label)
  Sigma_YY <- extract.variance_covariance(Sigma, what="YY", masque_data)
  logdetSigma_YY <- determinant(Sigma_YY, logarithm = TRUE)$modulus
  m_Y <- extract.simulate(sim, where="tips", what="expectations")
  LL <- ntaxa * log(2*pi) + logdetSigma_YY
  m_Y_vec <- as.vector(m_Y)[!missing]
#  LL <- LL + t(Y_data_vec - m_Y_vec) %*% Sigma_YY_inv %*% (Y_data_vec - m_Y_vec)
  LL <- LL + tcrossprod(t(Y_data_vec - m_Y_vec) %*% Sigma_YY_chol_inv)
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