# {Parcimony Number}
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
## Here is a function to compute the number of partimonious allocation of shifts
## on the tree, given a clustering of the tips.
## Dependencies : generic_functions.R
###############################################################################
##
# parcimonyNumber (phylo,clusters=rep(1,length(phy$tip.label)))
# PARAMETERS:
# @phylo (tree) imput tree
# @clusters (vector) : vector indicating the clusters of each tip
# RETURNS:
# (matrix) matrix with ncluster columns and Nnodes+ntaxa rows. Each row i contains the numbers of parcimonious reconstruction of the subtree under the node i, if starting with state j (column)
# DEPENDENCIES:
# init.parcimonyNumber, update.parcimonyNumber (, extrac.parcimonyNumber)
# PURPOSE:
# Find the number of parcimonious repartition of the shift that can lead to a given clustering
# NOTES:
# none
# REVISIONS:
# 19/05/14 - Initial release
# 21/05/14 - Gestion of case only one cluster. Extraction of the recursion.
##
parcimonyNumber <- function(phylo,clusters=rep(1,length(phylo$tip.label))){
  phy <- reorder(phylo,"postorder")
  ntaxa <- length(phy$tip.label)
  ## Initialization
  nbrReconstructions <- init.parcimonyNumber(phy,clusters)
  ## Tree recursion
  nbrReconstructions <- recursionUp(phy, nbrReconstructions, update.parcimonyNumber)
  attr(nbrReconstructions, "ntaxa") <- ntaxa
  return(nbrReconstructions)
}

##
# init.parcimonyNumber (phy,clusters)
# PARAMETERS:
# @(phy,clusters) see note above
# RETURNS:
# (matrix) matrix with ncluster columns and Nnodes+ntaxa rows. For each tip i, row i is the vector Ind(i=j). All other rows are set to NAs
# DEPENDENCIES:
# none
# PURPOSE:
# Initialise the matrix for function parcimonyNumber
# NOTES:
# none
# REVISIONS:
# 19/05/14 - Initial release
##
init.parcimonyNumber <- function(phy,clusters){
  ntaxa <- length(phy$tip.label)
  clus <- unique(clusters)
  nclus <- length(clus)
  nbrReconstructions <- matrix(NA,nrow=1 + nrow(phy$edge),ncol=nclus)
  for (i in 1:ntaxa){
    nbrReconstructions[i,] <- as.numeric(clusters[i]==clus)
  }
  return(nbrReconstructions)
}

##
# update.parcimonyNumber (daughtersNbr)
# PARAMETERS:
# @daughtersNbr (matrix) matrix with ncluster columns. Each row contains the result of parcimonyNumber for the children of a given node
# RETURNS:
# (vector) vector containing the number of parcimonious reconstructions from the current node if starting with cluster j
# DEPENDENCIES:
# none
# PURPOSE:
# Update
# NOTES:
# none
# REVISIONS:
# 19/05/14 - Initial release
##
update.parcimonyNumber <- function(daughtersParams, ...){
  nclus <- dim(daughtersParams)[2]
  nbrAdm <- rep(0,nclus)
  isNul <- daughtersParams==0
  costs <- colSums(1-isNul)
  Krond <- which(costs==max(costs))
  sums <- rowSums(daughtersParams)*isNul
  for (k in Krond){
    nbrAdm[k] <- prod(daughtersParams[,k]+sums[,k])
  }
  return(nbrAdm)
}

##
# extract.parcimonyNumber (nbrReconstructions,node=attr(nbrReconstructions, "ntaxa")+1)
# PARAMETERS:
# @nbrReconstructions (matrix) return of the function parcimonyNumber
# @node (int) : node from which to start (default value : root)
# RETURNS:
# (int) Number of parcimonious reconstructions for the subtree starting at node node.
# DEPENDENCIES:
# parcimonyNumber
# PURPOSE:
# Extract the values wanted from the raw result of function parcimonyNumber
# NOTES:
# none
# REVISIONS:
# 19/05/14 - Initial release
##
extract.parcimonyNumber <- function(nbrReconstructions,node=attr(nbrReconstructions, "ntaxa")+1){
  if (is.vector(nbrReconstructions)){
    return(1)
  } else{
    return(sum(nbrReconstructions[node,]))
  }
}

###############################################################################
## Find the clustering of the tips, given the shifts
###############################################################################

##
#' @title Tips descendants of nodes.
#'
#' @description
#' \code{enumerate_tips_under_edges} gives, for each edge of the tree, the labels
#' of the tips that have this edge as an ancestor.
#'
#' @details
#' This function uses function \code{prop.part} from package \code{ape}.
#'
#' @param tree phylogenetic tree
#' 
#' @return list of size n+m-1, entry i is the vector of tips bellow edge i.
#'
##
enumerate_tips_under_edges <- function (tree) {
  ntaxa <- length(tree$tip.label)
  temp <- prop.part(tree)
  subtree.list <- vector("list", nrow(tree$edge))
  for (i in 1:nrow(tree$edge)) {
    node <- tree$edge[i, 2]
    if (node > ntaxa) {
      subtree.list[[i]] <- temp[[node -ntaxa]]
    } else {
      subtree.list[[i]] <- node
    }
  }
  return(subtree.list)
}

##
#' @title Clusters of the tips corresponding to a list of shifts.
#'
#' @description
#' \code{clusters_from_shifts} take a vector of shifts edges, and gives the
#' clustering of the tips induced by them.
#'
#' @details
#' By default, this function uses \code{enumerate_tips_under_edges} to compute 
#' the list of tips under each edge.
#'
#' @param tree phylogenetic tree
#' @param edges a vector of edges of the tree, where the shifts are
#' @param part.list a list giving the descendant tips of each edge
#' 
#' @return list of size n+m-1, entry i is the vector of tips bellow edge i.
#'
##
clusters_from_shifts <- function (tree, edges, part.list = enumerate_tips_under_edges(tree)) {
  ntaxa <- length(tree$tip.label)
  part <- rep(0, ntaxa)
  edges <- sort(edges)
  for (i in 1:length(edges)) {
    part[part.list[[edges[i]]]] <- i
  }
  return(part)
}