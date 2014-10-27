# {General functions}
# Copyright (C) {2014}  {SR, MM, PB}
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
#            @phy (tree) Input tree
#            @x (int) Number of a node in the tree
# RETURNS:
#            (int) Number of parental node of node x in tree phy
# DEPENDENCIES:
#            none
# PURPOSE:
#            Get the ancestor of node x
# NOTES:
#            none
# REVISIONS:
#            16/05/14 - Initial release
##
getAncestor <- function(phy, x){
  if (x == Ntip(phy) + 1) return(NA)
  i <- which(phy$edge[, 2] == x)
  return(phy$edge[i, 1])
}

replaceInList <- function (x, FUN, ...)  {
  if (is.list(x)) {
    for (i in seq_along(x)) {
      x[i] <- list(replaceInList(x[[i]], FUN, ...))
    }
    x
  }
  else FUN(x, ...)
}

##
# correspondanceEdges (edges,from,to)
# PARAMETERS:
#            @edges (vector) Vector of index of edges in the tree "from"
#            @from (tree) Initial tree
#            @to (tree) Destination tree, should be the same as tree "from", but with a different parametrization
# RETURNS:
#            (vector) Vector of index of edges in the tree "to"
# DEPENDENCIES:
#            none
# PURPOSE:
#            If the parametrization of the tree is changed (when put to "cladewise"), make the same changes on the edges where the shifts occur.
# NOTES:
#            none
# REVISIONS:
#            20/05/14 - Initial release
#            26/05/14 - Simplification using match
##
correspondanceEdges <- function(edges,from,to){
  mm <- match(from$edge[,2],to$edge[,2])
  newEdges <- mm[edges]
  return(newEdges)
}

##
# compute_times_ca (phy)
# PARAMETERS:
#            @phy (tree) imput tree
# RETURNS:
#            (matrix) : entry (i,j) of the matrix is t_ij, the time of shared ancestry between nodes i and j
# DEPENDENCIES:
#            (node.depth.edgelength, mrca)
# PURPOSE:
#            Compute t_ij
# NOTES:
#            none
# REVISIONS:
#            22/05/14 - Initial release
##
compute_times_ca <- function(phy) {
  times <- node.depth.edgelength(phy)
  prac <- mrca(phy,full=TRUE)
  times_ca <- matrix(times[prac],dim(prac))
  attr(times_ca, "ntaxa")  <- length(phy$tip.label)
  return(times_ca)
}

##
# compute_dist_phy (phy)
# PARAMETERS:
#            @phy (tree) imput tree
# RETURNS:
#            (matrix) : entry (i,j) of the matrix is d_ij, the phylogenetic distance between nodes i and j
# DEPENDENCIES:
#            (dist.nodes)
# PURPOSE:
#            Compute d_ij
# NOTES:
#            none
# REVISIONS:
#            22/05/14 - Initial release
##
compute_dist_phy <- function(phy) {
  dist_phy <- dist.nodes(phy)
  attr(dist_phy, "ntaxa")  <- length(phy$tip.label)
  return(dist_phy)
}

###############################################################################
## Functions to wander on the tree
###############################################################################
##
#' @title Generic recursion down the tree.
#'
#' @description
#' \code{recursionDown} uses the function \code{\link{updateDown}} to compute
#' daughters rows of matrix param.

#' @details
#' This functin is to be used in other more complex function that need to 
#' update a quantity from the root to the tips of a tree. Note that the 
#' input tree must be in claddewise order.
#'
#' @param phy Input tree, in cladewise order.
#' @param params Matrix of parameters to update by the recursion
#' @param updateDown Function to be used for the update
#' @param ... Arguments to be used by the function updateDown
#' 
#' @return Matrix of parameters updated.
#'
#'21/05/14 - Initial release
##
recursionDown <- function(phy, params, updateDown, ...) {
  if (attr(phy,"order") != "cladewise") stop("The tree must be in cladewise order")
  ## Tree recursion from root to tips
  for (e in 1:nrow(phy$edge)) {
    edge  <- phy$edge[e, ]
    length <- phy$edge.length[e]
    parent  <- edge[1]
    daughter  <- edge[2]
    params[daughter,] <- updateDown(edgeNbr = e,
                                    ancestral = params[parent,],
                                    length = length, ...)
  }
  return(params)
}

##
# recursionUp ( phy, params, updateUp, ... )
# PARAMETERS:
#            @phy (tree) Input tree, in postorder order
#            @params (matrix) Matrix of parameters to update by the recursion
#            @updateUp (function) Function to be used for the update
# RETURNS:
#            (matrix)  Matrix of parameters updated
# DEPENDENCIES:
#            none
# PURPOSE:
#            Do the recursion from the tips to the root. params is updated row after row.
# NOTES:
#            The input tree must be in postorder order
# REVISIONS:
#            21/05/14 - Initial release
##
recursionUp <- function(phy, params, updateUp, ...){
  if (attr(phy,"order") != "postorder") stop("The tree must be in postorder order")
  ## Tree recursion
  e <- 1
  while (e <= nrow(phy$edge)) {
    edge  <- phy$edge[e, ]
    parent  <- edge[1]
    ii <- which(phy$edge[,1]==parent)
    daughters  <- phy$edge[ii,2]
    params[parent,] <- updateUp(edgesNbr=ii,
                                daughters=daughters,
                                daughtersParams = params[daughters,,drop=F],
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
#            @n (int) tree with 2^n tips
# RETURNS:
#            (tree)  A symetric tree with 2^n tips
# DEPENDENCIES:
#            read.tree
# PURPOSE:
#            Generate a symetric tree
# NOTES:
#            none
# REVISIONS:
#            26/05/14 - Initial release
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
#            @n (int) tree with n tips
# RETURNS:
#            (tree)  A comb-like tree with n tips
# DEPENDENCIES:
#            read.tree
# PURPOSE:
#            Generate a comb-like tree
# NOTES:
#            none
# REVISIONS:
#            26/05/14 - Initial release
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
#' @param eps the tolerence for the selection strength
#' 
#' @return character : "BM" or "OU"
#'
#'16/06/14 - Initial release
##
check.selection.strength <- function(process, selection.strength=NA, eps=10^(-6), ...){
  if (process == "BM") {
    return("BM")
  } else if (abs(selection.strength) < eps) {
    warning(paste("The selection strength is too low (alpha<", eps, "), process is considered to be a simple Brownian Motion", sep=""))
    return("BM")
  } else {
    return("OU")
  }
}

##
#' @title Test state of root.
#'
#' @description
#' \code{test.root.state} test wether the parameters of root.state given
#' by the user are coherent. If not, it returns a new corrected list to
#' define root.state.
#'
#' @details
#' To test coherence, the following priorities are applied : random > stationary.root > values.root = exp.root = var.root
#'
#' @param root.state A list giving the root state
#' @param process "BM" or "OU"
#' @param ... parameters of the process (if OU)
#' 
#' @return Coherent list root.state.
#'
#' 28/05/14 - Initial release
##
test.root.state <- function(root.state, process=c("BM","OU"), ...) {
  process <- match.arg(process)
  process <- check.selection.strength(process, ...)
  if (process == "BM") {
    return(test.root.state.BM(root.state))
  } else if (process == "OU") {
    return(test.root.state.OU(root.state, ...))
  }
}

test.root.state.BM <- function(root.state, ...) {
  if (!is.null(root.state$stationary.root) && root.state$stationary.root){
    warning("The BM does not have a stationary state. root.state$stationary.root is set to NULL")
    root.state$stationary.root <- NULL
  }
  if (root.state$random && !is.na(root.state$value.root)) {
    warning("As root state is supposed random, its value is not defined and set to NA")
    root.state$value.root <- NA
  }
  if (!root.state$random && (!is.na(root.state$exp.root) || !is.na(root.state$exp.root))) {
    warning("As root state is supposed fixed, its expectation and variance are not defined and set to NA")
    root.state$exp.root <- NA
    root.state$var.root <- NA
  }
  return(root.state)
}

test.root.state.OU <- function(root.state, process, variance, selection.strength, optimal.value, ...) {
  if (root.state$random && !is.na(root.state$value.root)) {
    warning("As root state is supposed random, its value is not defined and set to NA")
    root.state$value.root <- NA
  }
  if (!root.state$random && (!is.na(root.state$exp.root) || !is.na(root.state$exp.root))) {
    warning("As root state is supposed fixed, its expectation and variance are not defined and set to NA")
    root.state$exp.root <- NA
    root.state$var.root <- NA
  }
  if (is.null(root.state$stationary.root)) {
    warning("root.state$stationary.root was not defined, and is now set to its default value (TRUE)")
    root.state$stationary.root <- TRUE
  }
  if (!root.state$random && root.state$stationary.root) {
    warning("As root state is supposed fixed, the root cannot be at its stationary state. root.state$stationary.root is set to FALSE")
    root.state$stationary.root <- FALSE
  }
  if (root.state$stationary.root &&
        (!isTRUE(all.equal(root.state$exp.root, optimal.value)) ||
           !isTRUE(all.equal(root.state$var.root, variance/(2 * selection.strength))))) {
    warning("As root is supposed to be at stationary distribution, mu=beta and gamma2=sigma2/(2*alpha)")
    root.state$exp.root <- optimal.value
    root.state$var.root <- variance/(2 * selection.strength)
  }
  return(root.state)
}

##
#' @title Log Likelihood of a model
#'
#' @description
#' \code{likelihood.OU} computes the likelihhod of the data given a model. This
#' is a non-efficient debugging function.
#'
#' @details
#' This function uses functions \code{compute_mean_variance.simple}, \code{compute_times_ca}, \code{compute_dist_phy}, \code{compute_log_likelihood.simple}
#'
#' @param Y the vector of the data at the tips
#' @param phylo a phylogenetic tree
#' @param params list of parameters with the correct structure
#' 
#' @return boolean
#'
#'02/10/14 - Initial release
##

log_likelihood.OU <- function(Y, phylo, params, ...) {
  moments <- compute_mean_variance.simple(phylo = phylo,
                                          times_shared = compute_times_ca(phylo),
                                          distances_phylo = compute_dist_phy(phylo),
                                          process = "OU",
                                          params_old = params, ...)
  LL <- compute_log_likelihood.simple(phylo = phylo,
                                      Y_data = Y, 
                                      sim = moments$sim, 
                                      Sigma = moments$Sigma, 
                                      Sigma_YY_inv = moments$Sigma_YY_inv)
  return(LL)
}
