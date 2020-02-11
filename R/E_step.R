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
#'  \code{extract_simulate_internal} to extract the needed quantities from these objects.
#'
#' @param phylo Input tree.
#' @param Y_data : vector indicating the data at the tips
#' @param sim (list) : result of function \code{simulate}
#' @param Sigma : variance-covariance matrix, result of function \code{compute_variance_covariance}
#' @param Sigma_YY_inv : invert of the variance-covariance matrix of the data
#' 
#' @return conditional_law_X (list) : list of conditional statistics :
#'                   "expectation" : matrix of size p x (ntaxa+Nnode), with ntaxa
#' first columns set to Y_data (tips), and from ntaxa+1 to conditional expectation
#' of the nodes conditioned to the tips E[Z_j|Y]
#'                   "variances" : array of size p x p x (ntaxa+Nnode) with ntaxa first 
#' matrices of zeros (tips) and conditional variance of the nodes conditioned to the tips 
#' Var[Z_j|Y]
#'                  "covariances" : array of size p x p x (ntaxa+Nnode) with ntaxa first 
#' matrices of zeros (tips) and conditional covariance of the nodes and their parents
#' conditioned to the tips Cov[Z_j,Z_pa(j)|Y], with NA for the root.
#'                   "optimal.values" : matrix of size p x ntaxa+Nnode of optimal
#' values beta(t_j)
#' 
#' @keywords internal
##

compute_E.simple <- function(phylo,
                             times_shared,
                             distances_phylo,
                             process,
                             params_old,
                             masque_data = c(rep(TRUE, attr(params_old, "p_dim") * length(phylo$tip.label)),
                                             rep(FALSE, attr(params_old, "p_dim") * phylo$Nnode)),
                             F_moments,
                             Y_data_vec_known,
                             miss = rep(FALSE, attr(params_old, "p_dim") * length(phylo$tip.label)),
                             Y_data,
                             U_tree, ...){
  moments <- compute_mean_variance.simple(phylo = phylo,
                                          times_shared = times_shared,
                                          distances_phylo = distances_phylo,
                                          process = process,
                                          params_old = params_old,
                                          masque_data = masque_data,
                                          F_moments = F_moments,
                                          U_tree = U_tree)
  log_likelihood_old <- compute_log_likelihood.simple(phylo = phylo,
                                                      Y_data_vec = Y_data_vec_known,
                                                      sim = moments$sim,
                                                      Sigma = moments$Sigma,
                                                      Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                                      miss = miss, 
                                                      masque_data = masque_data,
                                                      C_YY = F_moments$C_YY,
                                                      Y_data = Y_data,
                                                      C_YY_chol_inv = F_moments$C_YY_chol_inv,
                                                      R = params_old$variance)
  conditional_law_X <- compute_cond_law.simple(phylo = phylo,
                                               Y_data_vec = Y_data_vec_known,
                                               sim = moments$sim,
                                               Sigma = moments$Sigma,
                                               Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                               miss = miss,
                                               masque_data = masque_data,
                                               F_means = F_moments$F_means,
                                               F_vars = F_moments$F_vars,
                                               R = params_old$variance,
                                               Y_data = Y_data)
  return(list(log_likelihood_old = log_likelihood_old,
              conditional_law_X = conditional_law_X))
}

compute_cond_law.simple <- function (phylo, Y_data_vec, sim,
                                     Sigma, Sigma_YY_chol_inv,
                                     miss = rep(FALSE, dim(sim)[1] * length(phylo$tip.label)),
                                     masque_data = c(rep(TRUE, dim(sim)[1] * length(phylo$tip.label)),
                                                     rep(FALSE, dim(sim)[1] * phylo$Nnode)),
                                     ...) {
  ## Initialization
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  p <- dim(sim)[1]
  nMiss <- sum(miss)
  # index_missing <- (Nnode * p + 1):(Nnode * p + nMiss)
  index_missing <- c(rep(TRUE, nMiss), rep(FALSE, Nnode * p))
  conditional_law_X <- list(expectations = matrix(NA, p, ntaxa + Nnode), 
                            variances = array(NA, c(p, p, ntaxa + Nnode)), 
                            covariances = matrix(NA, c(p, p, ntaxa + Nnode)))
  ## Mean
  m_Y <- extract_simulate_internal(sim, where="tips", what="expectations")
  m_Z <- extract_simulate_internal(sim, where="nodes", what="expectations")
  conditional_law_X$optimal.values <- cbind(extract_simulate_internal(sim,
                                                             where="tips",
                                                             what="optimal.values"),
                                            extract_simulate_internal(sim,
                                                             where="nodes",
                                                             what="optimal.values")) # NULL if BM
  ## Variance Covariance
  Sigma_YZ <- extract.variance_covariance(Sigma, what="YZ", masque_data)
  Sigma_ZZ <- extract.variance_covariance(Sigma, what="ZZ", masque_data)
#  temp <- Sigma_YZ %*% Sigma_YY_inv
  temp <- Sigma_YZ %*% Sigma_YY_chol_inv
  # Y_data_vec <- as.vector(Y_data)
  m_Y_vec <- as.vector(m_Y)[!miss]
  m_Z_vec <- c(as.vector(m_Y)[miss], as.vector(m_Z))
  # Conditional expectation of unkonwn values
  exp_Z_vec <- m_Z_vec + temp %*% crossprod(Sigma_YY_chol_inv,
                                            (Y_data_vec - m_Y_vec))
  expcond <- rep(NA, (Nnode + ntaxa) * p)
  # Data
  expcond[masque_data] <- Y_data_vec
  # Missing Data and Nodes
  expcond[!masque_data] <- as.vector(exp_Z_vec)
  conditional_law_X$expectations <- matrix(expcond, nrow = p)
#  conditional_variance_covariance <- Sigma_ZZ - temp %*% t(Sigma_YZ)
  ## Variances
  conditional_variance_covariance <- Sigma_ZZ - Matrix::tcrossprod(temp)
  attr(conditional_variance_covariance, "p_dim") <- p
  conditional_variance_covariance_nodes <- conditional_variance_covariance[!index_missing, !index_missing]
  attr(conditional_variance_covariance_nodes, "p_dim") <- p
  # Data tips
  var_tips <- array(0, c(p, p, ntaxa))
  cov_tips <- array(0, c(p, p, ntaxa))
  if (nMiss > 0){
    conditional_variance_covariance_tips <- conditional_variance_covariance[index_missing, index_missing, drop = FALSE]
    conditional_variance_covariance_tips_nodes <- conditional_variance_covariance[!index_missing, index_missing, drop = FALSE]
    
    missing_tips <- (which(miss) - 1) %/% p + 1
    missing_tips_uniques <- unique(missing_tips)
    par_missing_tips <- getAncestors(phylo, missing_tips_uniques)
    par_missing_tips <- par_missing_tips - ntaxa
    par_missing_tips <- sapply(par_missing_tips,
                               function(z) (p * (z - 1) + 1):(p * z))
    par_missing_tips <- matrix(par_missing_tips, nrow = p)
    missing_chars <- (which(miss) - 1) %% p + 1
    grpes_missing <- matrix(sapply(1:ntaxa, function(z) missing_tips == z),
                            nrow = nMiss)
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
                                         var_nodes), c(p, p, ntaxa + Nnode))
  # Nodes - covariances
  cov_nodes <- extract.covariance_parents(phylo,
                                          conditional_variance_covariance_nodes)
  conditional_law_X$covariances <- array(c(cov_tips,
                                           cov_nodes), c(p, p, ntaxa + Nnode))
  return(conditional_law_X)
}

###############################################################################
## No Missing
###############################################################################
compute_E.simple.nomissing.BM <- function(phylo,
                                          times_shared,
                                          distances_phylo,
                                          process,
                                          params_old,
                                          masque_data = c(rep(TRUE, attr(params_old, "p_dim") * length(phylo$tip.label)),
                                                          rep(FALSE, attr(params_old, "p_dim") * phylo$Nnode)),
                                          F_moments,
                                          Y_data_vec_known,
                                          miss = rep(FALSE, attr(params_old, "p_dim") * length(phylo$tip.label)),
                                          Y_data,
                                          U_tree, ...){
  moments <- compute_mean_variance.simple.nomissing.BM(phylo = phylo,
                                                       times_shared = times_shared,
                                                       distances_phylo = distances_phylo,
                                                       process = process,
                                                       params_old = params_old,
                                                       masque_data = masque_data,
                                                       F_moments = F_moments,
                                                       U_tree = U_tree)
  log_likelihood_old <- compute_log_likelihood.simple.nomissing.BM(phylo = phylo,
                                                                   Y_data_vec = Y_data_vec_known,
                                                                   sim = moments$sim,
                                                                   Sigma = moments$Sigma,
                                                                   Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                                                   miss = miss, 
                                                                   masque_data = masque_data,
                                                                   C_YY = F_moments$C_YY,
                                                                   Y_data = Y_data,
                                                                   C_YY_chol_inv = F_moments$C_YY_chol_inv,
                                                                   R = params_old$variance)
  conditional_law_X <- compute_cond_law.simple.nomissing.BM(phylo = phylo,
                                                            Y_data_vec = Y_data_vec_known,
                                                            sim = moments$sim,
                                                            Sigma = moments$Sigma,
                                                            Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                                            miss = miss,
                                                            masque_data = masque_data,
                                                            F_means = F_moments$F_means,
                                                            F_vars = F_moments$F_vars,
                                                            R = params_old$variance,
                                                            Y_data = Y_data)
  return(list(log_likelihood_old = log_likelihood_old,
              conditional_law_X = conditional_law_X))
}


compute_cond_law.simple.nomissing.BM <- function (phylo, Y_data, sim,
                                                  F_means, F_vars, R, ...) {
  ## Initialization
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  p <- dim(sim)[1]
  conditional_law_X <- list(expectations = matrix(NA, p, ntaxa + Nnode), 
                            variances = array(NA, c(p, p, ntaxa + Nnode)), 
                            covariances = matrix(NA, c(p, p, ntaxa + Nnode)))
  ## Mean
  m_Y <- extract_simulate_internal(sim, where="tips", what="expectations")
  m_Z <- extract_simulate_internal(sim, where="nodes", what="expectations")
  # Conditional expectation of unkonwn values
  conditional_law_X$expectations <- m_Z + (Y_data- m_Y) %*% F_means
  conditional_law_X$expectations <- cbind(Y_data, conditional_law_X$expectations)
  conditional_law_X$expectations <- matrix(conditional_law_X$expectations, dim(conditional_law_X$expectations))
  ## Variances
  # Data tips
  var_tips <- array(0, c(p, p, ntaxa))
  cov_tips <- array(0, c(p, p, ntaxa))
  # Nodes - varariances
  R <- array(R, dim(R))
  var_nodes <- R %o% array(diag(F_vars))
  conditional_law_X$variances <- array(c(var_tips,
                                         var_nodes), c(p, p, ntaxa + Nnode))
  # Nodes - covariances
  daughters <- phylo$edge[phylo$edge[, 2] > ntaxa, 2] # Only the nodes
  parents <- phylo$edge[phylo$edge[, 2] > ntaxa, 1]
  cov_nodes <- diag(F_vars[daughters - ntaxa, parents - ntaxa])
  cov_nodes[daughters - ntaxa] <- cov_nodes
  cov_nodes[1] <- NA # Root is NA
  cov_nodes <- R %o% array(cov_nodes)
  conditional_law_X$covariances <- array(c(cov_tips,
                                           cov_nodes), c(p, p, ntaxa + Nnode))
  return(conditional_law_X)
}

##
#' @title Compute fixed moments for E step.
#'
#' @description
#' \code{compute_fixed_moments} compute the fixed matrices used in E step when
#' there is no missing data.
#' 
#' @param times_shared matrix of the tree, result of function
#'  \code{compute_times_ca}
#' 
#' @return F_means the matrix to use for actualization of means.
#' @return F_vars the matrix to use for the actualization of variances.
#' 
#' @keywords internal
#' 
##
compute_fixed_moments <- function(times_shared, ntaxa){
  masque_data <- c(rep(TRUE, ntaxa), rep(FALSE, dim(times_shared)[1] - ntaxa))
  C_YZ <- extract.variance_covariance(times_shared, what="YZ", masque_data)
  C_YY <- extract.variance_covariance(times_shared, what="YY", masque_data)
  C_ZZ <- extract.variance_covariance(times_shared, what="ZZ", masque_data)
  C_YY_chol <- chol(C_YY)
  C_YY_chol_inv <- backsolve(C_YY_chol, diag(ncol(C_YY_chol)))
  temp <- C_YZ %*% C_YY_chol_inv
  F_means <- tcrossprod(C_YY_chol_inv, temp)  # solve(C_YY) %*% t(C_YZ)
  F_vars <- C_ZZ - tcrossprod(temp)
  if (!isSymmetric(F_vars, tol = 1000 * .Machine$double.eps)){
    stop("Something went wrong, matrix F_vars sould be symmetric.")
  }
  F_vars <- forceSymmetric(F_vars)
  return(list(C_YY = C_YY, C_YY_chol_inv = C_YY_chol_inv,
              F_means = F_means, F_vars = F_vars))
}

##
#' @title Extract sub-matrices of variance.
#'
#' @description
#' \code{extract.variance_covariance} return the adequate sub-matrix.
#' 
#' @param struct structural matrix of size (ntaxa+Nnode)*p, result 
#' of function \code{compute_variance_covariance}
#' @param what: sub-matrix to be extracted:
#'                "YY" : sub-matrix of tips (p*ntaxa first lines and columns)
#'                "YZ" : sub matrix tips x nodes (p*Nnode last rows and p*ntaxa first columns)
#'                "ZZ" : sub matrix of nodes (p*Nnode last rows and columns)
#' @param miss; missing values of Y_data
#' 
#' @return sub-matrix of variance covariance.
#' 
#' @keywords internal
#' 
##
extract.variance_covariance <- function(struct, what=c("YY","YZ","ZZ"),
                                        masque_data = c(rep(TRUE, attr(struct, "ntaxa") * attr(struct, "p_dim")),
                                                        rep(FALSE, (dim(struct)[1] - attr(struct, "ntaxa")) * attr(struct, "p_dim")))){
  # ntaxa <- attr(struct, "ntaxa")
  # p <- attr(struct, "p_dim")
  if (what=="YY") {
    res <- struct[masque_data, masque_data]
    # res <- as(res, "dpoMatrix")
  } else if (what=="YZ") {
    res <- struct[!masque_data, masque_data]
  } else if (what=="ZZ") {
    res <- struct[!masque_data, !masque_data]
    # res <- as(res, "dpoMatrix")
  }
  return(res)
}

##
# extract.covariance_parents (phylo, struct)
# PARAMETERS:
#            @phylo (tree)
#            @struct (matrix) structural matrix of size ntaxa+Nnode, result of function compute_times_ca, compute_dist_phy or compute_variance_covariance
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
#' covariance of a node and its parent.
#' 
#' @param vars matrix of size p x p*(ntaxa+Nnode) result of function \code{compute_E.simple},
#' entry "variances" or "covariances".
#' @param node for which to extract the matrix.
#' 
#' @return sub-matrix of variance for the node.
#' 
#' @keywords internal
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
#' @param times_shared times of shared ancestry of all nodes and tips, result of function
#'  \code{compute_times_ca}
#' @param params_old (list) : old parameters to be used in the E step
#' 
#' @return matrix of variance covariance for the BM
#' 
#' @keywords internal
#' 
##
compute_variance_covariance.BM <- function(times_shared, params_old, ...) {
  p <- nrow(params_old$shifts$values)
  J <- matrix(1, nrow = dim(times_shared)[1], ncol = dim(times_shared)[2])
  if (p == 1){ # Dimension 1 (next would also work, but slightly faster)
    varr <- as.vector(params_old$variance) * times_shared
    if (params_old$root.state$random) {
      varr <- varr + as.vector(params_old$root.state$var.root) * J
    }
  } else {
    varr <- kronecker(times_shared, params_old$variance)
    if (params_old$root.state$random) {
      varr <- varr + kronecker(as(J, "symmetricMatrix"), params_old$root.state$var.root)
    }
  }
  attr(varr, "p_dim") <- p
  attr(varr, "ntaxa") <- attr(params_old, "ntaxa")
  return(varr)
}

##
#' @title Tips Variances for the BM
#'
#' @description
#' \code{compute_variance_block_diagonal.BM} computes the n p*p variance
#' matrices of each tip vector.
#'
#' @param times_shared times of shared ancestry of all nodes and tips, result 
#'  of function \code{compute_times_ca}
#' @param params_old (list) : old parameters to be used in the E step
#' 
#' @return p * p * ntaxa array with ntaxa variance matrices
#' 
#' @keywords internal
#' 
##
compute_variance_block_diagonal.BM <- function(times_shared,
                                               params_old,
                                               ntaxa, ...) {
  p <- dim(params_old$variance)[1]
  if (is.null(p)){
    stop("Variance should be a matrix.")
  }
  sigma2 <- as.matrix(params_old$variance)
  # random root if needed
  if (params_old$root.state$random){
    var_root <- as.matrix(params_old$root.state$var.root)
  } else {
    var_root <- matrix(0, p, p)
  }
  return(sigma2 %o% diag(times_shared)[1:ntaxa] + var_root %o% rep(1, ntaxa))
}

##
#' @title Matrix of tree-induced correlations for the BM
#'
#' @description
#' \code{compute_tree_correlations_matrix.BM} returns times_shared its provided argument.
#'
#' @param times_shared times of shared ancestry of all nodes and tips, result of function
#'  \code{compute_times_ca}
#' 
#' @return times_shared 
#' 
#' @keywords internal
##
compute_tree_correlations_matrix.BM <- function(times_shared, params_old, ...) {
  res <- times_shared
  attr(res, "p_dim") <- nrow(params_old$shifts$values)
  attr(res, "ntaxa") <- attr(params_old, "ntaxa")
  return(times_shared)
}

##
#' @title Complete variance covariance matrix for scOU
#'
#' @description
#' \code{compute_variance_covariance.scOU} computes the (n+m)*p squared variance covariance
#' matrix of vec(X).
#'
#' @param times_shared times of shared ancestry of all nodes and tips, result of function
#'  \code{compute_times_ca}
#' @param distances_phylo (matrix) : phylogenetic distance, result of function 
#' \code{compute_dist_phy}
#' @param params_old (list) : old parameters to be used in the E step
#' 
#' @return matrix of variance covariance for the scOU
#' 
#' @keywords internal
#' 
##
compute_variance_covariance.scOU <- function(times_shared, distances_phylo, params_old, ...) {
  p <- nrow(params_old$shifts$values)
  if (is.null(p)) p <- 1
  alpha <- as.vector(params_old$selection.strength)
  sigma2 <- params_old$variance
  S <- compute_tree_correlations_matrix.scOU(times_shared, distances_phylo, params_old, ...)
  varr <- kronecker(S, sigma2 / (2 * alpha))
  if (params_old$root.state$random && !params_old$root.state$stationary.root) {
    times_nodes <- list(Matrix::diag(times_shared))
    sum_times <- do.call('rbind',
                         rep(times_nodes, length(Matrix::diag(times_shared))))
    sum_times <- sum_times + do.call('cbind',
                                     rep(times_nodes, length(Matrix::diag(times_shared))))
    gamma2 <- params_old$root.state$var.root
    varr <- kronecker(exp(- alpha * sum_times), gamma2) + varr
  }
  # varr <- as(varr, "dpoMatrix")
  attr(varr, "p_dim") <- p
  attr(varr, "ntaxa") <- attr(params_old, "ntaxa")
  return(varr)
}

##
#' @title Matrix of tree-induced correlations for the scOU
#'
#' @description
#' \code{compute_tree_correlations_matrix.scOU} computes the (n+m)x(m+n) matrix of correlations
#' induced by the tree. It takes two cases in consideration: root fixed, or root in stationary
#' state.
#'
#' @param times_shared times of shared ancestry of all nodes and tips, result of function
#'  \code{compute_times_ca}
#' @param distances_phylo (matrix) : phylogenetic distance, result of function 
#' \code{compute_dist_phy}
#' @param params_old (list) : old parameters to be used in the E step
#' 
#' @return matrix of variance covariance for the scOU
#' 
#' @keywords internal
#' 
##
compute_tree_correlations_matrix.scOU <- function(times_shared, distances_phylo, params_old, ...) {
  p <- nrow(params_old$shifts$values)
  if (is.null(p)) p <- 1
  alpha <- as.vector(params_old$selection.strength)
  sigma2 <- params_old$variance
  if (params_old$root.state$stationary.root) {
    res <- exp(- alpha * distances_phylo)
  } else {
    res <-  (1 - exp(- 2 * alpha * times_shared)) * exp(- alpha * distances_phylo)
  }
  attr(res, "p_dim") <- p
  attr(res, "ntaxa") <- attr(params_old, "ntaxa")
  return(res)
}

##
#' @title Complete variance covariance matrix for OU
#'
#' @description
#' \code{compute_variance_covariance.OU} computes the (n+m)*p squared variance
#' covariance matrix of vec(X).
#'
#' @param times_shared times of shared ancestry of all nodes and tips, result 
#'  of function \code{compute_times_ca}
#' @param params_old (list) : old parameters to be used in the E step
#' 
#' @return matrix of variance covariance for the OU
#' 
#' @keywords internal
#' 
##
compute_variance_covariance.OU <- function(times_shared,
                                           distances_phylo,
                                           params_old, ...) {
  p <- dim(params_old$variance)[1]
  ntaxa <- dim(times_shared)[1]
  if (is.null(p)){
    return(compute_variance_covariance.scOU(times_shared,
                                            distances_phylo,
                                            params_old, ...))
  }
  alpha_mat <- params_old$selection.strength
  sigma2 <- params_old$variance
  varr <- matrix(NA, p*ntaxa, p*ntaxa)
  # random root if needed
  if (params_old$root.state$random){
    var_root <- params_old$root.state$var.root
  } else {
    var_root <- matrix(0, p, p)
  }
  eig_alpha <- eigen(alpha_mat, symmetric = TRUE)
  P <- eig_alpha$vectors
  sigma_trans <- t(P) %*% sigma2 %*% P
  var_root_trans <- t(P) %*% var_root %*% P
  vv <- matrix(NA, p, p)
  for (q in 1:p){
    for (r in 1:p){
      vv[q, r] <- 1/(eig_alpha$values[q] + eig_alpha$values[r])
    }
  }
  ee <- exp(eig_alpha$values)
  for (i in 1:ntaxa){
    for (j in 1:ntaxa){
      ti <- times_shared[i, i]
      tj <- times_shared[j, j]
      tij <- times_shared[i, j]
      cpij <- tcrossprod(ee^(-ti), ee^(-tj))
      temp <- vv * (tcrossprod(ee^tij) - 1)
      temp <- cpij * (temp * sigma_trans + var_root_trans)
      varr[((i-1)*(p)+1):(i*p), ((j-1)*(p)+1):(j*p)] <- as.matrix(P %*% temp %*% t(P))
#       temp <- vv * exp(-eig_alpha$values * ti) %*% t(exp(-eig_alpha$values * tj)) * (exp(eig_alpha$values * tij) %*% t(exp(eig_alpha$values * tij)) - 1)
#       varr[((i-1)*(p)+1):(i*p),((j-1)*(p)+1):(j*p)] <- as.matrix(P %*% (temp * sigma_trans + exp(-eig_alpha$values * ti) %*% t(exp(-eig_alpha$values * tj)) * var_root_trans) %*% t(P))
    }
  }
  attr(varr, "p_dim") <- p
  attr(varr, "ntaxa") <- attr(params_old, "ntaxa")
  return(varr)
}

##
#' @title Tips Variances for the OU
#'
#' @description
#' \code{compute_variance_block_diagonal.OU} computes the n p*p variance
#' matrices of each tip vector.
#'
#' @param times_shared times of shared ancestry of all nodes and tips, result 
#'  of function \code{compute_times_ca}
#' @param params_old (list) : old parameters to be used in the E step
#' 
#' @return p * p * ntaxa array with ntaxa variance matrices
#' 
#' @keywords internal
#' 
##
compute_variance_block_diagonal.OU <- function(times_shared,
                                               params_old,
                                               ntaxa, ...) {
  p <- dim(params_old$variance)[1]
  if (is.null(p)){
    stop("Variance should be a matrix.")
  }
  alpha_mat <- params_old$selection.strength
  sigma2 <- as.matrix(params_old$variance)
  # random root if needed
  if (params_old$root.state$random){
    var_root <- as.matrix(params_old$root.state$var.root)
  } else {
    var_root <- matrix(0, p, p)
  }
  eig_alpha <- eigen(alpha_mat)
  P <- eig_alpha$vectors
  P_inv <- solve(P)
  sigma_trans <- P_inv %*% sigma2 %*% t(P_inv)
  var_root_trans <- P_inv %*% var_root %*% t(P_inv)
  safe_exp <- function(a, b, t) {
    alpha <- a + b
    if (alpha == 0) return(t)
    return((exp(alpha * t) - 1) / alpha)
  }
  safe_exp_vec <- function(a, b, t) {
    mapply(function(x, y) safe_exp(x, y, t), a, b)
  }
  vv <- matrix(NA, p, p)
  for (q in 1:p){
    for (r in 1:p){
      vv[q, r] <- 1 / (eig_alpha$values[q] + eig_alpha$values[r])
    }
  }
  ee <- exp(eig_alpha$values)
  var1 <- sapply(diag(times_shared)[1:ntaxa],
                 function(ti) return(P %*% (tcrossprod(ee^(-ti), ee^(-ti)) * (outer(eig_alpha$values, eig_alpha$values, function(a, b) return(safe_exp_vec(a, b, ti))) * sigma_trans + var_root_trans)) %*% t(P)))
  return(array(var1, c(p, p, ntaxa)))
}

##
#' @title Complete variance covariance matrix for OU
#'
#' @description
#' \code{compute_variance_covariance.OU} computes the (n+m)*p squared variance
#' covariance matrix of vec(X).
#'
#' @param times_shared times of shared ancestry of all nodes and tips, result 
#'  of function \code{compute_times_ca}
#' @param params_old (list) : old parameters to be used in the E step
#' 
#' @return matrix of variance covariance for the OU
#' 
#' @keywords internal
#' 
##
compute_variance_covariance.OU.nonsym <- function(times_shared,
                                                  distances_phylo,
                                                  params_old, ...) {
  p <- dim(params_old$variance)[1]
  ntaxa <- dim(times_shared)[1]
  if (is.null(p)){
    return(compute_variance_covariance.scOU(times_shared,
                                            distances_phylo,
                                            params_old, ...))
  }
  alpha_mat <- params_old$selection.strength
  sigma2 <- params_old$variance
  varr <- matrix(NA, p*ntaxa, p*ntaxa)
  # random root if needed
  if (params_old$root.state$random){
    var_root <- params_old$root.state$var.root
  } else {
    var_root <- matrix(0, p, p)
  }
  eig_alpha <- eigen(alpha_mat)
  P <- eig_alpha$vectors
  P_inv <- solve(P)
  sigma_trans <- P_inv %*% sigma2 %*% t(P_inv)
  var_root_trans <- P_inv %*% var_root %*% t(P_inv)
  vv <- matrix(NA, p, p)
  for (q in 1:p){
    for (r in 1:p){
      vv[q, r] <- 1/(eig_alpha$values[q] + eig_alpha$values[r])
    }
  }
  ee <- exp(eig_alpha$values)
  for (i in 1:ntaxa){
    print(i)
    for (j in 1:ntaxa){
      ti <- times_shared[i, i]
      tj <- times_shared[j, j]
      tij <- times_shared[i, j]
      cpij <- tcrossprod(ee^(-ti), ee^(-tj))
      temp <- vv * (tcrossprod(ee^tij) - 1)
      temp <- cpij * (temp * sigma_trans + var_root_trans)
      varr[((i-1)*(p)+1):(i*p), ((j-1)*(p)+1):(j*p)] <- as.matrix(P %*% temp %*% t(P))
      #       temp <- vv * exp(-eig_alpha$values * ti) %*% t(exp(-eig_alpha$values * tj)) * (exp(eig_alpha$values * tij) %*% t(exp(eig_alpha$values * tij)) - 1)
      #       varr[((i-1)*(p)+1):(i*p),((j-1)*(p)+1):(j*p)] <- as.matrix(P %*% (temp * sigma_trans + exp(-eig_alpha$values * ti) %*% t(exp(-eig_alpha$values * tj)) * var_root_trans) %*% t(P))
    }
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
#' @param distances_phylo (matrix) : phylogenetic distance, result of function 
#' \code{compute_dist_phy}
#' @param process a two letter string indicating the process to consider
#' @param params_old a list of parameters to be used in the computations
#' 
#' @return sim (list) : result of function \code{simulate} with the appropriate
#'  parameters
#' @return Sigma matrix of variance covariance, result of function 
#' \code{compute_variance_covariance}
# #@return Sigma_YY_inv inverse of variance matrix of the data
#' @return Sigma_YY_chol_inv invert of Cholesky matrix of Sigma_YY:
#'  (Sigma_YY)^(-1) = tcrossprod(Sigma_YY_chol_inv)
#'  
#' @keywords internal
# 29/09/14 - Initial release
##
compute_mean_variance.simple <- function(phylo,
                                         times_shared,
                                         distances_phylo,
                                         process=c("BM", "OU", "rBM", "scOU"),
                                         params_old,
                                         masque_data = c(rep(TRUE, attr(params_old, "p_dim") * length(phylo$tip.label)),
                                                         rep(FALSE, attr(params_old, "p_dim") * phylo$Nnode)),
                                         sim = NULL, 
                                         U_tree = NULL, ...) {
  ## Choose process 
  process  <- match.arg(process)
  if (attr(params_old, "p_dim") == 1 && process == "OU") process <- "scOU" 
  compute_variance_covariance  <- switch(process, 
                                         BM = compute_variance_covariance.BM,
                                         OU = compute_variance_covariance.OU,
                                         scOU = compute_variance_covariance.scOU)
  ## Mean
  if (is.null(sim)){
    sim <- simulate_internal(phylo = phylo, 
                             process = process,
                             p = attr(params_old, "p_dim"),
                             root.state = params_old$root.state, 
                             shifts = params_old$shifts, 
                             variance = params_old$variance, 
                             optimal.value = params_old$optimal.value, 
                             selection.strength = params_old$selection.strength,
                             simulate_random = FALSE,
                             U_tree = U_tree)
  }
  ## Variance Covariance
  Sigma <- compute_variance_covariance(times_shared = times_shared, 
                                       distances_phylo = distances_phylo,
                                       params_old = params_old)
  Sigma_YY <- extract.variance_covariance(Sigma, what="YY", masque_data = masque_data)
  Sigma_YY_chol <- chol(Sigma_YY)
  Sigma_YY_chol_inv <- backsolve(Sigma_YY_chol, Matrix::diag(ncol(Sigma_YY_chol)))
  #Sigma_YY_inv <- chol2inv(Sigma_YY_chol)
  return(list(sim = sim, Sigma = Sigma, Sigma_YY_chol_inv = Sigma_YY_chol_inv))
}

compute_mean_variance.simple.nomissing.BM <- function (phylo,
                                                       process = "BM",
                                                       params_old,
                                                       F_moments,
                                                       U_tree = NULL, ...) {
  ## Mean
  sim <- simulate_internal(phylo = phylo, 
                           process = process,
                           p = attr(params_old, "p_dim"),
                           root.state = params_old$root.state, 
                           shifts = params_old$shifts, 
                           variance = params_old$variance, 
                           optimal.value = params_old$optimal.value, 
                           selection.strength = params_old$selection.strength,
                           simulate_random = FALSE,
                           U_tree = U_tree)
  return(list(sim = sim, C_YY = F_moments$C_YY,
              C_YY_chol_inv = F_moments$C_YY_chol_inv,
              F_means = F_moments$F_means, F_vars = F_moments$F_vars))
}

#####################################################################
## Functions to compute the log likelihood
#####################################################################

##
#' @title Residuals
#'
#' @description
#' \code{compute_residuals.simple} computes the residuals after the fit of the
#' data in the simple case where the inverse of the variance matrix is given.
#'
#' @details
#' This function takes parameters sim and Sigma_YY_inv from  
#' \code{compute_mean_variance.simple}. It uses function \code{extract_simulate_internal}
#' to extract the needed quantities from these objects.
#'
#' @param phylo Input tree.
#' @param Y_data : vector indicating the data at the tips.
#' @param sim (list) : result of function \code{simulate}.
#' @param Sigma_YY_inv : invert of the variance-covariance matrix of the data.
#' 
#' @keywords internal
#' 
# @return vector of residuals
##
compute_residuals.simple <- function(phylo, Y_data_vec, sim,
                                     Sigma_YY_chol_inv, miss){
  ntaxa <- length(phylo$tip.label)
  m_Y <- extract_simulate_internal(sim, where="tips", what="expectations")
  m_Y <- as.vector(m_Y)[!miss]
  resis <- Sigma_YY_chol_inv %*% (Y_data_vec - m_Y)
  return(resis)
}

##
#' @title Squared Mahalanobis Distance
#'
#' @description
#' \code{compute_mahalanobis_distance.simple} computes the squared Mahalanobis distance 
#' between the data and mean at tips  of the data in the simple case where the inverse
#' of the variance matrix is given.
#'
#' @details
#' This function takes parameters sim and Sigma_YY_inv from  
#' \code{compute_mean_variance.simple}. It uses function \code{extract_simulate_internal}
#' to extract the needed quantities from these objects.
#'
#' @param phylo Input tree.
#' @param Y_data : vector indicating the data at the tips.
#' @param sim (list) : result of function \code{simulate}.
#' @param Sigma_YY_inv : invert of the variance-covariance matrix of the data.
#' 
#' @return squared Mahalanobis distance between data and mean at the tips.
#' 
#' @keywords internal
##
compute_mahalanobis_distance.simple <- function(phylo, Y_data_vec, sim,
                                                Sigma_YY_chol_inv,
                                                miss = rep(FALSE, dim(sim)[1] * length(phylo$tip.label)),
                                                ...){
  ntaxa <- length(phylo$tip.label)
  m_Y <- extract_simulate_internal(sim, where="tips", what="expectations")
  m_Y <- as.vector(m_Y)[!miss]
#  MD <- t(Y_data - m_Y)%*%Sigma_YY_inv%*%(Y_data - m_Y)
  MD <- tcrossprod(t(Y_data_vec - m_Y) %*% Sigma_YY_chol_inv)
  return(MD)
}

compute_mahalanobis_distance.simple.nomissing.BM <- function(phylo, Y_data, sim,
                                                             C_YY_chol_inv, R, ...){
  ntaxa <- length(phylo$tip.label)
  m_Y <- extract_simulate_internal(sim, where="tips", what="expectations")
#  MD <- t(Y_data - m_Y)%*%Sigma_YY_inv%*%(Y_data - m_Y)
  R_chol <- t(chol(R)) # R = R_chol %*% t(R_chol)
  R_chol_inv <- t(backsolve(t(R_chol), diag(ncol(R_chol))))
  MD <- R_chol_inv %*% (Y_data - m_Y) %*% C_YY_chol_inv
  return(sum(MD^2))
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
#'  \code{extract_simulate_internal} to extract the needed quantities from these objects.
#'
#' @param phylo Input tree.
#' @param Y_data : vector indicating the data at the tips
#' @param sim (list) : result of function \code{simulate}
#' @param Sigma : variance-covariance matrix, result of function \code{compute_variance_covariance}
#' @param Sigma_YY_inv : invert of the variance-covariance matrix of the data
#' 
#' @return log likelihood of the data
#' 
#' @keywords internal
#' 
# 29/09/14 - Initial release
##
compute_log_likelihood.simple <- function(phylo, Y_data_vec, sim,
                                          Sigma, Sigma_YY_chol_inv,
                                          miss = rep(FALSE, dim(sim)[1] * length(phylo$tip.label)),
                                          masque_data = c(rep(TRUE, dim(sim)[1] * length(phylo$tip.label)),
                                                          rep(FALSE, dim(sim)[1] * phylo$Nnode)), ...){
  # ntaxa <- length(phylo$tip.label)
  Sigma_YY <- extract.variance_covariance(Sigma, what="YY", masque_data)
  logdetSigma_YY <- Matrix::determinant(Sigma_YY, logarithm = TRUE)$modulus
  m_Y <- extract_simulate_internal(sim, where="tips", what="expectations")
  LL <- length(Y_data_vec) * log(2*pi) + logdetSigma_YY
  m_Y_vec <- as.vector(m_Y)[!miss]
#  LL <- LL + t(Y_data_vec - m_Y_vec) %*% Sigma_YY_inv %*% (Y_data_vec - m_Y_vec)
  LL <- LL + tcrossprod(t(Y_data_vec - m_Y_vec) %*% Sigma_YY_chol_inv)
  return(-LL/2)
}

compute_log_likelihood.simple.nomissing.BM <- function(phylo, Y_data, sim,
                                                       C_YY, C_YY_chol_inv, R, ...){
  # ntaxa <- length(phylo$tip.label)
  logdetC_YY <- determinant(C_YY, logarithm = TRUE)$modulus
  logdetR <- determinant(R, logarithm = TRUE)$modulus
  m_Y <- extract_simulate_internal(sim, where="tips", what="expectations")
  LL <- prod(dim(Y_data)) * log(2*pi) + dim(R)[1] * logdetC_YY + dim(C_YY)[1] * logdetR
  LL <- LL + compute_mahalanobis_distance.simple.nomissing.BM(phylo, Y_data, sim,
                                                              C_YY_chol_inv, R)
  return(-LL/2)
}

# compute_log_det.simple <- function(phylo, Y_data_vec, sim,
#                                    Sigma, Sigma_YY_chol_inv,
#                                    miss, masque_data){
#   ntaxa <- length(phylo$tip.label)
#   Sigma_YY <- extract.variance_covariance(Sigma, what="YY", masque_data)
#   logdetSigma_YY <- determinant(Sigma_YY, logarithm = TRUE)$modulus
#   return( - (ntaxa * log(2*pi) + logdetSigma_YY) / 2)
# }

# compute_log_maha.simple <- function(phylo, Y_data_vec, sim,
#                                     Sigma, Sigma_YY_chol_inv,
#                                     miss, masque_data){
#   ntaxa <- length(phylo$tip.label)
#   m_Y <- extract_simulate_internal(sim, where="tips", what="expectations")
#   m_Y_vec <- as.vector(m_Y)[!miss]
#   LL <- tcrossprod(t(Y_data_vec - m_Y_vec) %*% Sigma_YY_chol_inv)
#   return(-LL/2)
# }

# compute_entropy.simple <- function(Sigma, Sigma_YY_inv){
#   Sigma_YZ <- extract.variance_covariance(Sigma, what="YZ")
#   Sigma_ZZ <- extract.variance_covariance(Sigma, what="ZZ")
#   conditional_variance_covariance <- Sigma_ZZ - Sigma_YZ%*%Sigma_YY_inv%*%t(Sigma_YZ)
#   N <- dim(Sigma_ZZ)[1]
#   logdet_conditional_variance_covariance <- determinant(conditional_variance_covariance, logarithm = TRUE)$modulus
#   return(N/2*log(2*pi*exp(1)) + 1/2*logdet_conditional_variance_covariance)
# }

# compute_log_likelihood_with_entropy.simple <- function(CLL, H){
#   return(CLL + H)
# }

# ## Likelihood and Mahalanobis computation for independent traits
# lik_maha_ind <- function(phylo,
#                          times_shared,
#                          distances_phylo,
#                          process,
#                          independent,
#                          params_old,
#                          Y_data,
#                          masque_data = c(rep(TRUE, dim(Y_data)[1] * length(phylo$tip.label)),
#                                          rep(FALSE, dim(Y_data)[1] * phylo$Nnode)),
#                          F_moments = NULL,
#                          Y_data_vec_known = as.vector(Y_data),
#                          miss = rep(FALSE, dim(Y_data)[1] * length(phylo$tip.label))){
#   if (independent){
#     # if independent, params_old is a list of p params
#     params_old <- split_params_independent(params_old)
#     masque_data_matr <- matrix(masque_data,
#                                ncol = length(phylo$tip.label) + phylo$Nnode)
#     miss_matr <- matrix(miss,
#                         ncol = length(phylo$tip.label))
#     res <- vector(mode = "list", length = length(params_old))
#     for (i in 1:length(res)){
#       res[[i]] <- lik_maha_ind(phylo = phylo,
#                                times_shared = times_shared,
#                                distances_phylo = distances_phylo,
#                                process = process,
#                                params_old = params_old[[i]],
#                                masque_data = masque_data_matr[i, ],
#                                F_moments = F_moments,
#                                independent = FALSE,
#                                Y_data_vec_known = Y_data[i, !miss_matr[i, ]],
#                                miss = miss_matr[i, ],
#                                Y_data = Y_data[i, , drop = F])
#     }
#     return(res)
#   }
#   moments <- compute_mean_variance.simple(phylo = phylo,
#                                           times_shared = times_shared,
#                                           distances_phylo = distances_phylo,
#                                           process = process,
#                                           params_old = params_old,
#                                           masque_data = masque_data,
#                                           F_moments = F_moments)
#   log_likelihood <- compute_log_likelihood.simple(phylo = phylo,
#                                                   Y_data_vec = Y_data_vec_known,
#                                                   sim = moments$sim,
#                                                   Sigma = moments$Sigma,
#                                                   Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
#                                                   miss = miss, 
#                                                   masque_data = masque_data,
#                                                   C_YY = F_moments$C_YY,
#                                                   Y_data = Y_data,
#                                                   C_YY_chol_inv = F_moments$C_YY_chol_inv,
#                                                   R = params_old$variance)
#   ## Compute Mahalanobis norm between data and mean at tips
#   maha_data_mean <- compute_mahalanobis_distance.simple(phylo = phylo,
#                                                         Y_data_vec = Y_data_vec_known,
#                                                         sim = moments$sim,
#                                                         Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
#                                                         miss = miss,
#                                                         Y_data = Y_data,
#                                                         C_YY_chol_inv = F_moments$C_YY_chol_inv,
#                                                         R = params_old$variance)
#   return(list(log_likelihood = log_likelihood,
#               maha_data_mean = maha_data_mean))
# }

###############################################################################
## Upward_downward
###############################################################################

#' @useDynLib PhylogeneticEM, .registration=TRUE
#' @importFrom Rcpp evalCpp
# @import RcppArmadillo

compute_E.upward_downward <- function(phylo,
                                      Y_data,
                                      process,
                                      params_old,
                                      U_tree, ...){
  Delta <- shifts.list_to_matrix(phylo, params_old$shifts)
  params_old$root.state$var.root <- as.matrix(params_old$root.state$var.root)
  if (process == "BM"){
    params_old$variance <- as.matrix(params_old$variance)
    return(upward_downward_BM(Y_data, phylo$edge, Delta,
                              params_old$variance, phylo$edge.length,
                              params_old$root.state))
  } else if (process == "OU"){
    Beta1 <- tcrossprod(Delta, U_tree) + params_old$optimal.value
    Beta <- Beta1[, phylo$edge[, 2], drop = F] # re-order by edge
    if (params_old$root.state$stationary.root){
      Stationary_Var <- as.matrix(params_old$root.state$var.root)
    } else {
      Stationary_Var <- as.matrix(compute_stationary_variance(params_old$variance, params_old$selection.strength))
    }
    params_old$selection.strength <- as.matrix(params_old$selection.strength)
    res <- upward_downward_OU(Y_data, phylo$edge,
                              Beta, Stationary_Var,
                              phylo$edge.length, params_old$selection.strength,
                              params_old$root.state)
    res$conditional_law_X$optimal.values <- Beta1
    return(res)
    
  } else if (process == "scOU"){
    Beta1 <- tcrossprod(Delta, U_tree) + params_old$optimal.value
    Beta <- Beta1[, phylo$edge[, 2], drop = F] # re-order by edge
    if (params_old$root.state$stationary.root){
      Stationary_Var <- as.matrix(params_old$root.state$var.root)
    } else {
      if (params_old$selection.strength[1] < 0){
        Stationary_Var <- -as.matrix(compute_stationary_variance(params_old$variance, -params_old$selection.strength))
      } else {
        Stationary_Var <- as.matrix(compute_stationary_variance(params_old$variance, params_old$selection.strength))
      }
    }
    Alpha <- params_old$selection.strength * diag(rep(1, ncol(Stationary_Var)))
    res <- upward_downward_OU(Y_data, phylo$edge,
                              Beta, Stationary_Var,
                              phylo$edge.length, Alpha,
                              params_old$root.state)
    res$conditional_law_X$optimal.values <- Beta1
    return(res)
  }
}

##
#' @title Log Likelihood of a fitted object
#'
#' @description
#' \code{log_likelihood} computes the log likelihood of some parameters.
#'
#' @param x an object of class \code{\link{params_process}} or \code{\link{PhyloEM}}.
#' @param Y_data matrix of data at the tips, size p x ntaxa. Each line is a
#' trait, and each column is a tip. The column names are checked against the
#' tip names of the tree.
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
# @param U_tree (optional) full incidence matrix of the tree, result of function
#' \code{\link{incidence.matrix.full}}. Can be specified to avoid extra computations.
#' @param ... for a \code{PhyloEM} object, further arguments to be passed on to
#' \code{\link{params_process.PhyloEM}} (to choose which parameters to extract from
#' the results, see documentation of this function).
#'     
#' @return
#' The log likelihood of the data with the provided parameters on the tree.
#'  
#' @seealso \code{\link{params_process}}, \code{\link{PhyloEM}}
#' 
#' @export
#' 
##
log_likelihood <- function(x, ...) UseMethod("log_likelihood")

##
#' @describeIn log_likelihood \code{\link{params_process}} object
#' @export
##
log_likelihood.params_process <- function(x,
                                          Y_data,
                                          phylo, ...){
  phy <- reorder(phylo, order = "postorder")
  # Trace edges
  x$shifts$edges <- correspondanceEdges(edges = x$shifts$edges,
                                        from = phylo, to = phy)
  Delta <- shifts.list_to_matrix(phy, x$shifts)
  x$root.state$var.root <- as.matrix(x$root.state$var.root)
  ## BM Process
  if (x$process == "BM"){
    x$variance <- as.matrix(x$variance)
    return(log_likelihood_BM(Y_data, phy$edge, Delta,
                              x$variance, phy$edge.length,
                              x$root.state))
    ## OU process
  } else if (x$process == "OU"){
    U_tree <- incidence.matrix.full(phy)
    Beta1 <- tcrossprod(Delta, U_tree) + x$optimal.value
    Beta <- Beta1[, phy$edge[, 2], drop = F] # re-order by edge
    if (x$root.state$stationary.root){
      Stationary_Var <- as.matrix(x$root.state$var.root)
    } else {
      Stationary_Var <- as.matrix(compute_stationary_variance(x$variance, x$selection.strength))
    }
    x$selection.strength <- as.matrix(x$selection.strength)
    res <- log_likelihood_OU(Y_data, phy$edge,
                              Beta, Stationary_Var,
                              phy$edge.length, x$selection.strength,
                              x$root.state)
    return(res)
    ## scOU process
  } else if (x$process == "scOU"){
    U_tree <- incidence.matrix.full(phy)
    Beta1 <- tcrossprod(Delta, U_tree) + x$optimal.value
    Beta <- Beta1[, phy$edge[, 2], drop = F] # re-order by edge
    if (x$root.state$stationary.root){
      Stationary_Var <- as.matrix(x$root.state$var.root)
    } else {
      if (x$selection.strength[1] < 0){
        Stationary_Var <- -as.matrix(compute_stationary_variance(x$variance,
                                                                 -unique(diag(x$selection.strength))))
      } else {
        Stationary_Var <- as.matrix(compute_stationary_variance(x$variance, 
                                                                unique(diag(x$selection.strength))))      
      }
    }
    Alpha <- x$selection.strength * diag(rep(1, ncol(Stationary_Var)))
    res <- log_likelihood_OU(Y_data, phy$edge,
                             Beta, Stationary_Var,
                             phy$edge.length, Alpha,
                             x$root.state)
    return(res)
  }
}

##
#' @describeIn log_likelihood \code{\link{PhyloEM}} object
#' @export
##
log_likelihood.PhyloEM <- function(x, ...){
  
  params <- params_process(x, ...)
  
  return(log_likelihood.params_process(params,
                                       x$Y_data,
                                       x$phylo))
}

##
#' @title Residuals of a fitted object
#'
#' @description
#' \code{residuals} computes the residuals of some parameters.
#'
#' @param object an object of class \code{\link{params_process}} or \code{\link{PhyloEM}}.
# @param Y_data matrix of data at the tips, size p x ntaxa. Each line is a
#' trait, and each column is a tip. The column names are checked against the
#' tip names of the tree.
# @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
# @param U_tree (optional) full incidence matrix of the tree, result of function
#' \code{\link{incidence.matrix.full}}. Can be specified to avoid extra computations.
#' @param ... for a \code{PhyloEM} object, further arguments to be passed on to
#' \code{\link{params_process.PhyloEM}} (to choose which parameters to extract from
#' the results, see documentation of this function).
#'     
#' @return
#' The log likelihood of the data with the provided parameters on the tree.
#'  
#' @seealso \code{\link{params_process}}, \code{\link{PhyloEM}}
#' 
#' @export
##
residuals.PhyloEM <- function(object, ...){
  
  m_Y <- imputed_traits(object, trait = 1:object$p,
                        where = c("tips"),
                        what = c("expectations"), ...)
  
  return(object$Y_data - m_Y)
}

###############################################################################
## Wrapper (independent case)
###############################################################################
##
#' @title Wrapper for E step in EM
#'
#' @description
#' \code{wrapper_E_step} is used in the EM algorithm. It calls itself
#' recursively in case of independent parameters.
#' 
#' @keywords internal
##
wrapper_E_step <- function(phylo,
                           times_shared,
                           distances_phylo,
                           process,
                           params_old,
                           masque_data,
                           F_moments,
                           independent,
                           Y_data_vec_known,
                           miss,
                           Y_data,
                           U_tree,
                           compute_E){
  if (independent){
    # if independent, params_old is a list of p params
    masque_data_matr <- matrix(masque_data,
                               ncol = length(phylo$tip.label) + phylo$Nnode)
    miss_matr <- matrix(miss,
                        ncol = length(phylo$tip.label))
    res <- vector(mode = "list", length = length(params_old))
    for (i in 1:length(res)){
      res[[i]] <- wrapper_E_step(phylo = phylo,
                                 times_shared = times_shared,
                                 distances_phylo = distances_phylo,
                                 process = process,
                                 params_old = params_old[[i]],
                                 masque_data = masque_data_matr[i, ],
                                 F_moments = F_moments,
                                 independent = FALSE,
                                 Y_data_vec_known = Y_data[i, !miss_matr[i, ]],
                                 miss = miss_matr[i, ],
                                 Y_data = Y_data[i, , drop = F],
                                 U_tree = U_tree,
                                 compute_E = compute_E)
    }
    return(res)
  }
  return(compute_E(phylo = phylo,
                   times_shared = times_shared,
                   distances_phylo = distances_phylo,
                   process = process,
                   params_old = params_old,
                   masque_data = masque_data,
                   F_moments = F_moments,
                   independent = FALSE,
                   Y_data_vec_known = Y_data_vec_known,
                   miss = miss,
                   Y_data = Y_data,
                   U_tree = U_tree))
}
