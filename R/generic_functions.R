# {General functions}
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

########################################################
# Here are some general functions used in all the files.
########################################################
##
# getAncestor ( phy, x )
# PARAMETERS:
# @phy (tree) Input tree
# @x (int) Number of a node in the tree
# RETURNS:
# (int) Number of parental node of node x in tree phy
# DEPENDENCIES:
# none
# PURPOSE:
# Get the ancestor of node x
# NOTES:
# none
# REVISIONS:
# 16/05/14 - Initial release
##
getAncestor <- function(phy, x){
  if (x == Ntip(phy) + 1) return(NA)
  i <- which(phy$edge[, 2] == x)
  return(phy$edge[i, 1])
}
getAncestors <- function(phy, x){
  i <- match(x, phy$edge[, 2])
  return(phy$edge[i, 1])
}
replaceInList <- function (x, FUN, ...) {
  if (is.list(x)) {
    for (i in seq_along(x)) {
      x[i] <- list(replaceInList(x[[i]], FUN, ...))
    }
    x
  }
  else FUN(x, ...)
}

##
#' @title Correspondence between edges numbers
#'
#' @description
#' \code{correspondenceEdges} takes edges numbers on an input tree, and gives back
#' their corresponding numbers on the output tree.
#'
#' @param edges vector of index of edges in the tree "from"
#' @param from initial input tree (format "\code{phylo}")
#' @param to aimed output tree (format "\code{phylo}")
#'
#' @return vector of index of edges in the tree "to"
#'
#' @export
#26/05/14
##
correspondenceEdges <- function(edges, from, to){
  mm <- match(from$edge[, 2], to$edge[, 2])
  newEdges <- mm[edges]
  return(newEdges)
}

correspondanceEdges <- correspondenceEdges

##
# compute_times_ca (phy)
# PARAMETERS:
# @phy (tree) input tree
# RETURNS:
# (matrix) : entry (i,j) of the matrix is t_ij, the time of shared ancestry between nodes i and j
# DEPENDENCIES:
# (node.depth.edgelength, mrca)
# PURPOSE:
# Compute t_ij
# NOTES:
# none
# REVISIONS:
# 22/05/14 - Initial release
##
##
#' @title Common Ancestors Times
#'
#' @description
#' \code{compute_times_ca} computes the times t_ij between the root and the common
#' ancestor of two tips i, j.
#' 
#' @details
#' This function relies on \code{ape} functions
#' \code{\link[ape]{node.depth.edgelength}} and \code{\link[ape]{mrca}}.
#'
#' @param phy a phylogenetic tree of class \code{\link[ape]{phylo}}.
#'
#' @return a matrix of times of shared evolution, ordered as the tips of the
#' tree. The matrix is of type \code{\link[Matrix]{symmetricMatrix-class}}.
#'
#' @export
##
compute_times_ca <- function(phy) {
  times <- ape::node.depth.edgelength(phy)
  prac <- ape::mrca(phy,full=TRUE)
  times_ca <- matrix(times[prac],dim(prac))
  # attr(times_ca, "ntaxa") <- length(phy$tip.label)
#   times_ca <- phy$root.edge + times_ca # Add the root length
  return(as(times_ca, "symmetricMatrix"))
}

##
# compute_dist_phy (phy)
# PARAMETERS:
# @phy (tree) input tree
# RETURNS:
# (matrix) : entry (i,j) of the matrix is d_ij, the phylogenetic distance between nodes i and j
# DEPENDENCIES:
# (dist.nodes)
# PURPOSE:
# Compute d_ij
# NOTES:
# none
# REVISIONS:
# 22/05/14 - Initial release
##
##
#' @title Phylogenetic Distances
#'
#' @description
#' \code{compute_dist_phy} computes the phylogenetic distances d_ij between all the
#' tips i, j.
#' 
#' @details
#' This function relies on \code{ape} function
#' \code{\link[ape]{dist.nodes}}.
#'
#' @param phy a phylogenetic tree of class \code{\link[ape]{phylo}}.
#'
#' @return a matrix of phylogenetic distances, ordered as the tips of the
#' tree. The matrix is of type \code{\link[Matrix]{symmetricMatrix-class}}.
#'
#' @export
##
compute_dist_phy <- function(phy) {
  dist_phy <- ape::dist.nodes(phy)
  attr(dist_phy, "ntaxa") <- length(phy$tip.label)
  return(as(dist_phy, "symmetricMatrix"))
}

scale.tree <- function(phylo){
  if (!is.ultrametric(phylo)) stop("The tree is not ultrametric")
  ntaxa <- length(phylo$tip.label)
  height <- min(ape::node.depth.edgelength(phylo)[1:ntaxa]) - .Machine$double.eps^0.5# take the min so that any error is above 1
  phylo$edge.length <- phylo$edge.length/height
  return(phylo)
}

###############################################################################
## Functions to wander on the tree
###############################################################################
##
#' @title Generic recursion down the tree.
#'
#' @description
#' \code{recursionDown} uses the function \code{updateDown} to compute
#' daughters rows of matrix param.
#' @details
#' This function is to be used in other more complex function that need to
#' update a quantity from the root to the tips of a tree. Note that the
#' input tree must be in cladewise order.
#'
#' @param phy Input tree, in cladewise order.
#' @param params Matrix of parameters to update by the recursion
#' @param updateDown Function to be used for the update
#' @param ... Arguments to be used by the function updateDown
#'
#' @return Matrix of parameters updated.
#' 
#' @keywords internal
#'
##
recursionDown <- function(phy, params, updateDown, ...) {
  if (attr(phy,"order") != "cladewise") stop("The tree must be in cladewise order")
  ## Choose function to subset
  if (hasArg(subset_node)){
    subset_node <- list(...)[["subset_node"]]
  } else {
    subset_node <- subset_node.default
  }
  if (hasArg(allocate_subset_node)){
    allocate_subset_node <- list(...)[["allocate_subset_node"]]
  } else {
    allocate_subset_node <- allocate_subset_node.default
  }
  ## Tree recursion from root to tips
  for (e in 1:nrow(phy$edge)) {
    edge <- phy$edge[e, ]
    length <- phy$edge.length[e]
    parent <- edge[1]
    daughter <- edge[2]
    params <- allocate_subset_node(daughter, params,
                                   updateDown(edgeNbr = e,
                                              ancestral = subset_node(parent, params),
                                              length = length, ...))
  }
  return(params)
}

allocate_subset_node.default <- function(node, matrix, value){
  matrix[node, ] <- value
  return(matrix)
}

subset_node.default <- function(node, matrix){
  return(matrix[node, ])
}

##
# recursionUp ( phy, params, updateUp, ... )
# PARAMETERS:
# @phy (tree) Input tree, in postorder order
# @params (matrix) Matrix of parameters to update by the recursion
# @updateUp (function) Function to be used for the update
# RETURNS:
# (matrix) Matrix of parameters updated
# DEPENDENCIES:
# none
# PURPOSE:
# Do the recursion from the tips to the root. params is updated row after row.
# NOTES:
# The input tree must be in postorder order
# REVISIONS:
# 21/05/14 - Initial release
##
recursionUp <- function(phy, params, updateUp, ...){
  if (attr(phy,"order") != "postorder") stop("The tree must be in postorder order")
  ## Tree recursion
  e <- 1
  while (e <= nrow(phy$edge)) {
    edge <- phy$edge[e, ]
    parent <- edge[1]
    ii <- which(phy$edge[,1]==parent)
    daughters <- phy$edge[ii,2]
    params[parent,] <- updateUp(edgesNbr=ii,
                                daughters=daughters,
                                daughtersParams = params[daughters,,drop=F],
                                parent = parent, ...)
    e <- ii[length(ii)]+1
  }
  return(params)
}

recursionUp_list <- function(phy, params, updateUp, ...){
  if (attr(phy,"order") != "postorder") stop("The tree must be in postorder order")
  ## Tree recursion
  e <- 1
  while (e <= nrow(phy$edge)) {
    edge <- phy$edge[e, ]
    parent <- edge[1]
    ii <- which(phy$edge[,1]==parent)
    daughters <- phy$edge[ii,2]
    params[[parent]] <- updateUp(edgesNbr=ii,
                                daughters=daughters,
                                daughtersParams = params[daughters],
                                parent = parent, ...)
    e <- ii[length(ii)]+1
  }
  return(params)
}

###############################################################################
## Functions to generate trees with fixed topologies
###############################################################################
##
# rtree.sym (n)
# PARAMETERS:
# @n (int) tree with 2^n tips
# RETURNS:
# (tree) A symetric tree with 2^n tips
# DEPENDENCIES:
# read.tree
# PURPOSE:
# Generate a symetric tree
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
##
rtree.sym <- function(n){
  tree <- "A"
  for (k in 1:n) {
    tree <- paste("(", tree, ",", tree, ")", sep="")
  }
  return(read.tree(text=paste(tree, ";", sep="")))
}

##
# rtree.comb (n)
# PARAMETERS:
# @n (int) tree with n tips
# RETURNS:
# (tree) A comb-like tree with n tips
# DEPENDENCIES:
# read.tree
# PURPOSE:
# Generate a comb-like tree
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
##
rtree.comb <- function(n){
  if (n == 1) return(read.tree(text="(A);"))
  tree <- "A"
  for (k in 2:n) {
    tree <- paste("(A,", tree, ")", sep="")
  }
  return(read.tree(text=paste(tree, ";", sep="")))
}

###############################################################################
## Functions to test the parameters of the processes
###############################################################################
##
#' @title Check selection strength
#'
#' @description
#' \code{check.selection.strength} checks, if the process is an OU, if the
#' selection strength is not too low, in which case the process is replaced
#' with a BM.
#'
#' @details
#' This function return a process in a form of a character. If the entry
#' process is "BM", then nothing is done. If it is "OU", then the verification
#' of the selection strength is done.
#'
#' @param process : "BM" or "OU"
#' @param selection.strength the selection strength parameter (if OU)
#' @param eps the tolerance for the selection strength
#'
#' @return character : "BM" or "OU"
#'
#' @keywords internal
#16/06/14 - Initial release
##
check.selection.strength <- function(process, selection.strength = NA,
                                     eps = 10^(-6), ...){
  if (process == "BM") {
    return("BM")
  } else if (sum(abs(selection.strength)) < eps) {
    warning(paste("The selection strength is too low (L1-norm<", eps, "), process is considered to be a simple Brownian Motion", sep=""))
    return("BM")
  } else if (any(Re(eigen(selection.strength)$values) < eps)) {
    warning("All the eigen values of the selection strengh do not have a strictly positive real part. That might cause some issue. Proceed with caution.")
  }
  return(process)
}

##
#' @title Test state of root.
#'
#' @description
#' \code{test.root.state} test whether the parameters of root.state given
#' by the user are coherent. If not, it returns a new corrected list to
#' define root.state.
#'
#' @details
#' To test coherence, the following priorities are applied:
#'  random > stationary.root > values.root = exp.root = var.root
#'
#' @param root.state A list giving the root state
#' @param process "BM", "OU" or "scOU"
#' @param ... parameters of the process (if OU)
#'
#' @return Coherent list root.state.
#'
#' @keywords internal
# 28/05/14 - Initial release
##
test.root.state <- function(root.state, process=c("BM", "OU", "scOU", "OUBM"), ...) {
  process <- match.arg(process)
  process <- check.selection.strength(process, ...)
  if (process == "BM") {
    return(test.root.state.BM(root.state))
  } else if (process %in% c("OU", "scOU", "OUBM")) {
    return(test.root.state.OU(root.state, ...))
  }
}

test.root.state.BM <- function(root.state, ...) {
  if (!is.null(root.state$stationary.root) && root.state$stationary.root){
    warning("The BM does not have a stationary state. root.state$stationary.root is set to NULL")
    root.state$stationary.root <- NULL
  }
  if (root.state$random && !anyNA(root.state$value.root)) {
    warning("As root state is supposed random, its value is not defined and set to NA")
    root.state$value.root <- NA
    root.state$var.root <- as(root.state$var.root, "dpoMatrix")
  }
  if (!root.state$random && (!anyNA(root.state$exp.root) || !anyNA(root.state$exp.root))) {
    warning("As root state is supposed fixed, its expectation and variance are not defined and set to NA")
    root.state$exp.root <- NA
    root.state$var.root <- NA
  }
  return(root.state)
}

test.root.state.OU <- function(root.state, process, variance, selection.strength, optimal.value, ...) {
  if (root.state$random && !anyNA(root.state$value.root)) {
    warning("As root state is supposed random, its value is not defined and set to NA")
    root.state$value.root <- NA
    root.state$var.root <- as(root.state$var.root, "dpoMatrix")
  }
  if (!root.state$random && (!anyNA(root.state$exp.root) || !anyNA(root.state$exp.root))) {
    warning("As root state is supposed fixed, its expectation and variance are not defined and set to NA")
    root.state$exp.root <- NA
    root.state$var.root <- NA
  }
  if (is.null(root.state$stationary.root)) {
    warning("root.state$stationary.root was not defined, and is now set to its default value")
    if (root.state$random){
      root.state$stationary.root <- TRUE 
    } else {
      root.state$stationary.root <- FALSE
    }
  }
  if (!root.state$random && root.state$stationary.root) {
    warning("As root state is supposed fixed, the root cannot be at its stationary state. root.state$stationary.root is set to FALSE")
    root.state$stationary.root <- FALSE
  }
#   if (root.state$stationary.root &&
#         (!isTRUE(all.equal(root.state$exp.root, optimal.value)) ||
#            !isTRUE(all.equal(root.state$var.root, variance/(2 * selection.strength))))) {
#     warning("As root is supposed to be at stationary distribution, mu=beta and gamma2=sigma2/(2*alpha)")
#     root.state$exp.root <- optimal.value
#     root.state$var.root <- variance/(2 * selection.strength)
#   }
  root.state <- coherence_stationary_case(root.state, optimal.value,
                                           variance, selection.strength)
  return(root.state)
}

coherence_stationary_case <- function(root.state, optimal.value,
                                       variance, selection.strength){
  if (!root.state$stationary.root){
    return(root.state) ## Do nothing
  } else {
    if (!isTRUE(all.equal(root.state$exp.root, optimal.value))){
      root.state$exp.root <- optimal.value
      warning("As root is supposed to be in stationary case, root expectation was set to be equal to optimal value.")
    }
    
    root_var_expected <- compute_stationary_variance(variance, selection.strength)
    if(!isTRUE(all.equal(root.state$var.root, root_var_expected))){
      root.state$var.root <- as(root_var_expected, "dpoMatrix")
      warning("As the root is supposed to be in stationary state, root variance Gamma was set to: vec(Gamma) = (A kro_plus A)^{-1}vec(R).")
    }
    return(root.state)
  }
}

##
#' @title Compute the stationary variance matrix 
#'
#' @description
#' \code{compute_stationary_variance} computes the stationary variance matrix of
#' an OU process.
#'
#' @param variance the variance (rate matrix) of the process.
#' @param selection.strength the selection strength (alpha) matrix of the 
#' process.
#'
#' @return A positive definite Matrix of class \code{\link[Matrix]{dpoMatrix-class}}.
#' 
#' @export
##
compute_stationary_variance <- function(variance, selection.strength){
  if (is.null(selection.strength)) return(NA)
  if (length(as.vector(selection.strength)) == 1){
    vv <- as.matrix(variance) / (2 * selection.strength)
    vv <- Matrix(vv)
    vv <- as(vv, "dpoMatrix")
    return(vv)
  } else if (isDiagonal(selection.strength)) {
    dd <- diag(selection.strength)
    vv <- variance / outer(dd, dd, "+")
    vv <- as(vv, "dpoMatrix")
  } else {
    variance_vec <- as.vector(variance)
    kro_sum_A <- kronecker_sum(selection.strength, selection.strength)
    kro_sum_A_inv <- solve(kro_sum_A)
    root_var_vec <- kro_sum_A_inv %*% variance_vec
    gamma <- matrix(root_var_vec, dim(variance))
    if (!isSymmetric(gamma, tol = (.Machine$double.eps)^(0.7))) stop("Error in computation of stationary variance: matrix computed was not symmetric.")
    gamma <- symmpart(gamma)
    gamma <- nearPD(gamma)
    return(gamma$mat)
    # return(as(gamma, "symmetricMatrix"))
    # return(as(gamma, "dsyMatrix"))
    # return(as(gamma, "dpoMatrix"))
  }
}

compute_variance_from_stationary <- function(var.root, selection.strength){
  if (dim(var.root)[1] == 1){
    return(var.root * (2 * selection.strength))
  } else {
    var.root_vec <- as.vector(var.root)
    kro_sum_A <- kronecker_sum(selection.strength, selection.strength)
    variance_vec <- kro_sum_A%*%var.root_vec
    return(matrix(variance_vec, dim(var.root)))
  }
}

kronecker_sum <- function(M, N){
  if (!is.matrix(M) || !is.matrix(N))
    stop("Entries of Kronecker sum must be matrices")
  if ((length(dim(M)) != 2) || (length(dim(N)) != 2)) 
    stop("Entries of Kronecker sum must be matrices")
  if ((dim(M)[1] != dim(M)[2]) || (dim(N)[1] != dim(N)[2]))
    stop("Entries of Kronecker sum must be squared matrice.")
  m <- dim(M)[1]; Im <- diag(1, m, m)
  n <- dim(N)[1]; In <- diag(1, n, n)
  return(kronecker(M, Im) + kronecker(In, N))
}

##
# @title Log Likelihood of a model
#
# @description
# \code{likelihood.OU} computes the likelihhod of the data given a model. This
# is a non-efficient debugging function.
#
# @details
# This function uses functions \code{compute_mean_variance.simple}, \code{compute_times_ca}, \code{compute_dist_phy}, \code{compute_log_likelihood.simple}
#
# @param Y the vector of the data at the tips
# @param phylo a phylogenetic tree
# @param params list of parameters with the correct structure
#
# @return boolean
# 
# @keywords internal
#02/10/14 - Initial release
##
# log_likelihood.OU <- function(Y, phylo, params, ...) {
#   moments <- compute_mean_variance.simple(phylo = phylo,
#                                           times_shared = compute_times_ca(phylo),
#                                           distances_phylo = compute_dist_phy(phylo),
#                                           process = "OU",
#                                           params_old = params, ...)
#   LL <- compute_log_likelihood.simple(phylo = phylo,
#                                       Y_data = Y,
#                                       sim = moments$sim,
#                                       Sigma = moments$Sigma,
#                                       Sigma_YY_inv = moments$Sigma_YY_inv)
#   return(LL)
# }

##
#' @title Check dimensions of the parameters
#'
#' @description
#' \code{check_dimensions} checks dimensions of the parameters. 
#' If wrong, throw an error.
#'
#' @param p dimension of the trait simulated
#' @param root.state (list) state of the root, with:
#'     random : random state (TRUE) or deterministic state (FALSE)
#'     value.root : if deterministic, value of the character at the root
#'     exp.root : if random, expectation of the character at the root
#'     var.root : if random, variance of the character at the root
#' @param shifts (list) position and values of the shifts :
#'     edges : vector of the K id of edges where the shifts are
#'     values : matrix p x K of values of the shifts on the edges (one column = one shift)
#'     relativeTimes : vector of dimension K of relative time of the shift from the
#'     parent node of edges
#' @param variance variance-covariance matrix size p x p 
#' @param selection.strength matrix of selection strength size p x p (OU)
#' @param optimal.value vector of p optimal values at the root (OU)
#'     
#' @return Nothing
#' 
#' @keywords internal
#' 
# 25/08/15 - Multivariate
##

check_dimensions <- function(p,
                             root.state, shifts, variance,
                             selection.strength = NULL, optimal.value = NULL){
  root.state <- check_dimensions.root.state(p, root.state)
  #if (!is.null(unlist(shifts)))
  shifts <- check_dimensions.shifts(p, shifts)
  variance <- check_dimensions.matrix(p, p, variance, "variance")
  variance <- as(variance, "dpoMatrix")
  if (!is.null(selection.strength))
    if (is.vector(selection.strength) && length(selection.strength) == p){
      selection.strength <- diag(selection.strength, ncol = length(selection.strength))
    }
    if (is.vector(selection.strength) && length(selection.strength) == 1){
      selection.strength <- diag(rep(selection.strength, p))
    }
    selection.strength <- check_dimensions.matrix(p, p, selection.strength, "selection strength")
  if (!is.null(optimal.value))
    optimal.value <- check_dimensions.vector(p, optimal.value, "optimal value")
  
  return(params = list(root.state = root.state,
                       shifts = shifts,
                       variance = variance,
                       selection.strength = selection.strength,
                       optimal.value = optimal.value))
}

check_dimensions.matrix <- function(p, q, matrix, name = "matrix"){
  if (is.null(matrix)) matrix <- matrix(0, p, q)
  if (p == 1){
    if (is.vector(matrix) && length(matrix) != q) 
      stop(paste0(matrix, " should be a scalar in dimension q = ", q, "."))
    dim(matrix) <- c(1, q)
  }
  if (is.vector(matrix) || !all(dim(matrix) == c(p, q))) 
    stop(paste0("Dimensions of ", matrix, " matrix do not match"))
  return(matrix)
}

check_dimensions.vector <- function(p, v, name = "vector"){
  if (!is.vector(v)) stop(paste0(name, " should be a vector."))
  if (length(v) != p) 
    stop(paste0("Dimensions of ", name, " do not match"))
  return(v)
}

check_dimensions.root.state <- function(p, root.state){
  if (root.state$random){
    root.state$exp.root <- check_dimensions.vector(p, root.state$exp.root, "Root Expectation")
    root.state$var.root <- check_dimensions.matrix(p, p, root.state$var.root, "root variance")
    root.state$var.root <- as(root.state$var.root, "dpoMatrix")
  } else {
    root.state$value.root <- check_dimensions.vector(p, root.state$value.root, "Root Value")
  }
  return(root.state)
}

check_dimensions.shifts <- function(p, shifts){
  K <- length(shifts$edges)
  shifts$values <- check_dimensions.matrix(p, K, shifts$values, "shifts values")
  if (sum(shifts$relativeTimes) == 0) # If all zero, re-formate
    shifts$relativeTimes <- rep(0, K)
  shifts$relativeTimes <- check_dimensions.vector(K, shifts$relativeTimes, "shifts relative Times")
  return(shifts)
}

##
#' @title Find a reasonable grid for alpha
#' 
#' @description Grid so that 
#' 2*ln(2)*quantile(d_ij)/factor_up_alpha < t_1/2 < factor_down_alpha * ln(2) * h_tree,
#' with t_1/2 the phylogenetic half life: t_1/2 = log(2)/alpha.
#' Ensures that for alpha_min, it is almost a BM, and for alpha_max,
#' almost all the tips are decorrelated.
#' 
#' @details
#' If \code{quantile_low_distance=0}, then \code{quantile(d_ij)=min(d_ij)}, and, for any
#' two tips i,j, the correlation between i and j is bounded by exp(-factor_up_alpha/2).
#' Those values of alpha will be used for the re-scaling of the tree, which has an 
#' exponential term in exp(2*alpha*h). The function makes sure that this number is
#' below the maximal float allowed (equals to \code{.Machine$double.xmax}).
#'
#' @param phy phylogenetic tree of class "\code{phylo}"
#' @param alpha fixed vector of alpha values if already known. Default to NULL.
#' @param nbr_alpha the number of elements in the grid
#' @param factor_up_alpha factor for up scalability
#' @param factor_down_alpha factor for down scalability
#' @param quantile_low_distance quantile for min distance
#' @param log_transform whether to take a log scale for the spacing of alpha
#' values. Default to TRUE.
#' @param allow_negative whether to allow negative values for alpha (Early Burst).
#' See documentation of \code{\link{PhyloEM}} for more details. Default to FALSE.
#' @param ... not used.
#'     
#' @return A grid of alpha values
#' 
#' @seealso \code{\link{transform_branch_length}}, \code{\link{.Machine}}
#' 
#' @export
#' 
##
find_grid_alpha <- function(phy, alpha = NULL,
                            nbr_alpha = 10,
                            factor_up_alpha = 2,
                            factor_down_alpha = 3,
                            quantile_low_distance = 0.0001,
                            log_transform = TRUE, 
                            allow_negative = FALSE, ...){
  if (!is.null(alpha)) return(alpha)
  dtips <- cophenetic(phy)
  d_min <- quantile(dtips[dtips > 0], quantile_low_distance)
  h_tree <- node.depth.edgelength(phy)[1]
  alpha_min <- 1 / (factor_down_alpha * h_tree)
  alpha_max <- factor_up_alpha / (2 * d_min)
  alpha_max_machine <- log(.Machine$double.xmax^0.975)/(2*h_tree)
  if (alpha_max > alpha_max_machine){
    warning("The chosen alpha_max was above the machine precision. Taking alpha_max as the largest possible on this machine.")
    alpha_max <- alpha_max_machine
  }
  if (allow_negative) nbr_alpha <- nbr_alpha %/% 2
  if (log_transform){
    alpha_grid <- exp(seq(log(alpha_min), log(alpha_max), length.out = nbr_alpha))
  } else {
    alpha_grid <- seq(alpha_min, alpha_max, length.out = nbr_alpha)
  }
  if (allow_negative){
    alpha_min_neg_machine <- log(.Machine$double.eps^0.9)/(2*h_tree)
    if (log_transform){
      alpha_grid <- c(-exp(seq(log(-alpha_min_neg_machine), log(alpha_min), length.out = nbr_alpha)), 0, alpha_grid)
    } else {
      alpha_grid <- c(seq(alpha_min_neg_machine, -alpha_min, length.out = nbr_alpha), 0, alpha_grid)
    }
  }
  return(alpha_grid)
}

##
#' @title Check range of alpha
#' 
#' @description Check that the chosen values of alpha are not too large
#' or too small, in order to avoid numerical instabilities.
#'
#' @param alpha a vector of alpha values.
#' @param h_tree the total height of the tree.
#' 
#' @keywords internal
#' 
##
check_range_alpha <- function(alpha, h_tree){
  alpha_max_machine <- log(.Machine$double.xmax^0.98)/(2*h_tree)
  alpha_min_machine <- log(.Machine$double.eps^0.98)/(2*h_tree)
  if (any(alpha > alpha_max_machine)) {
    stop(paste0("The value for the selection strengh you took is too big, and will lead to numerical instabilities. Please consider using a value below ", alpha_max_machine))
  }
  if (any(alpha < alpha_min_machine)) {
    stop(paste0("The value for the selection strengh you took is too small, and will lead to numerical instabilities. Please consider using a value above ", alpha_min_machine))
  }
}

##
#' @title Transform branch length for a re-scaled BM
#' 
#' @description Re-scale the branch length of the tree so that a BM running
#' on the new tree produces the same observations at the tips than an OU with
#' parameter alpha.
#'
#' @param phylo A phylogenetic tree of class \code{\link[ape]{phylo}}, with branch
#' lengths.
#' @param alp Value of the selection strength.
#'     
#' @return phylo The same phylogenetic tree, with transformed branch lengths.
#' 
#' @export
# 25/08/15 - Multivariate
##
transform_branch_length <- function(phylo, alp){
  if (alp == 0){
    return(phylo)
  } else {
    nodes_depth <- node.depth.edgelength(phylo)
    h_tree <- nodes_depth[1]
    fun <- function(z){
      return(1 / (2 * alp) * exp(- 2 * alp * h_tree) * (exp(2 * alp * nodes_depth[z[2]]) - exp(2 * alp * nodes_depth[z[1]])))
    }
    ## Root edge if exists
    phylo$edge.length <- apply(phylo$edge, 1, fun)
    if (!is.null(phylo$root.edge)){
      phylo$root.edge <- 1 / (2 * alp) * exp(- 2 * alp * h_tree) * phylo$root.edge
    }
    return(phylo)
  }
}

##
#' @title Scale variance and selection strength from a linear transform
#' 
#' @description Used for process equivalencies on re-scaled trees.
#'
#' @param params Parameters list
#' @param f Factor of the linear transform. If t' = f * t, the function takes
#' parameters from phylo' back to phylo.
#'     
#' @return re-scaled parameters
#' 
#' @keywords internal
#' 
##
scale_params <- function(params, f){
  if (!is.null(params$variance)) params$variance <- f * params$variance
  if (!is.null(params$selection.strength)) params$selection.strength <- f * params$selection.strength
  return(params)
}

##
#' @title Split independent parameters into a list of parameters
#' 
#' @description \code{split_params_independent} split a params object for a 
#' process with p independent traits into p params objects.
#' The reverse operation is done by \code{merge_params_independent}
#'
#' @param params parameters
#'     
#' @return A list of p parameters
#' 
#' @keywords internal
#' 
##
split_params_independent <- function(params){
  p <- dim(params$variance)[1]
  params_split <- vector(mode = "list", length = p)
  for (l in 1:p){
    params_split[[l]] <- params
    params_split[[l]]$variance <- params$variance[l, l]
    params_split[[l]]$selection.strength <- params$selection.strength[l, l]
    if (!is.null(params$shifts$edges)){
      params_split[[l]]$shifts$values <- params$shifts$values[l, ]
    }
    if (!anyNA(params$root.state$value.root)){
      params_split[[l]]$root.state$value.root <- params$root.state$value.root[l]
    }
    if (!anyNA(params$root.state$exp.root)){
      params_split[[l]]$root.state$exp.root <- params$root.state$exp.root[l]
    }
    if (!anyNA(params$root.state$var.root)){
      params_split[[l]]$root.state$var.root <- params$root.state$var.root[l, l]
    }
    params_split[[l]]$optimal.value <- params$optimal.value[l]
    params_split[[l]] <- check_dimensions(1,
                                          params_split[[l]]$root.state,
                                          params_split[[l]]$shifts,
                                          params_split[[l]]$variance,
                                          params_split[[l]]$selection.strength,
                                          params_split[[l]]$optimal.value)
    if (!is.null(attr(params, "p_dim"))) attr(params_split[[l]], "p_dim") <- 1
  }
  return(params_split)
}

##
#' @title Merge a list of independent parameters into into one parameter
#' 
#' @description \code{merge_params_independent} merges a list of p params
#' objects into one param object of dimension p
#' The reverse operation is done by \code{split_params_independent}
#'
#' @param params_split a list of parameters
#'     
#' @return A parameter object
#' 
#' @keywords internal
#' 
##
merge_params_independent <- function(params_split){
  p <- length(params_split)
  params <- params_split[[1]]
  if (p > 1){
    params$variance <- diag(sapply(params_split, function(z) return(as.vector(z$variance))))
    if (!is.null(params$selection.strength)){
      params$selection.strength <- diag(sapply(params_split, function(z) return(z$selection.strength)))
    }
    if (length(params$shifts$edges) > 1){
      params$shifts$values <- t(sapply(params_split, function(z) return(z$shifts$values)))
    } else if (length(params$shifts$edges) == 1) {
      params$shifts$values <- sapply(params_split, function(z) return(z$shifts$values))
      dim(params$shifts$values) <- c(p,1)
    } else {
      params$shifts$values <- matrix(0, p, 0)
    }
    if (!anyNA(params$root.state$value.root)){
      params$root.state$value.root <- sapply(params_split, function(z) return(z$root.state$value.root))
    }
    if (!anyNA(params$root.state$exp.root)){
      params$root.state$exp.root <- sapply(params_split, function(z) return(z$root.state$exp.root))
    }
    if (!anyNA(as.vector(params$root.state$var.root))){
      params$root.state$var.root <- diag(sapply(params_split, function(z) return(as.vector(z$root.state$var.root))))
    }
    if (!is.null(params$optimal.value)){
      params$optimal.value <- sapply(params_split, function(z) return(z$optimal.value))
    }
  }
  params <- check_dimensions(p,
                             params$root.state,
                             params$shifts,
                             params$variance,
                             params$selection.strength,
                             params$optimal.value)
  if (!is.null(attr(params_split[[1]], "p_dim"))) attr(params, "p_dim") <- p
  return(params)
}