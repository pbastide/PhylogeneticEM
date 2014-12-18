# {Partitions Number}
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
## Here is a function to compute the number of partitions of the tips with a
## given number of groups.
## Dependencies : generic_functions.R
###############################################################################

#########################
## Main functions
#########################
##
# partitionsNumber (phylo, npart)
# PARAMETERS:
# @phylo (tree) imput tree
# @npart (int) : number of partitions of the tips allowed
# RETURNS:
# (matrix) matrix with Nnodes+ntaxa rows and 2*npart columns. Each column contains two vectors : for k=1:npart it contains the number of partitions with k groups compatible with the tree and the shift process; and for k=(npart+1):2*npart, it contains the number of "marqued" partitions with (k-npart) groups compatible with the tree and the shift process.
# DEPENDENCIES:
# init.partitionsNumber, update.partitionsNumber.bin (, extract.partitionsNumber)
# PURPOSE:
# Find the number of partitions with npart groups of the tips that are compatible with the structure of the tree and the shift process
# NOTES:
# Only works for rooted, binary trees.
# REVISIONS:
# 26/05/14 - Initial release
##
partitionsNumber <- function(phylo, npart){
  if (!is.rooted(phylo)) stop("The tree must be rooted !")
  if (is.binary.tree(phylo)){
    update.partitionsNumber <- update.partitionsNumber.bin
  } else {
    update.partitionsNumber <- update.partitionsNumber.gen
  }
  phy <- reorder(phylo,"postorder")
  ntaxa <- length(phy$tip.label)
  ## Initialization
  nbrCompatiblePartitions <- init.partitionsNumber(phy,npart)
  ## Tree recursion
  nbrCompatiblePartitions <- recursionUp(phy, nbrCompatiblePartitions, update.partitionsNumber)
  attr(nbrCompatiblePartitions, "ntaxa") <- ntaxa
  attr(nbrCompatiblePartitions, "npart") <- npart
  return(nbrCompatiblePartitions)
}

##
# init.partitionsNumber (phy,clusters)
# PARAMETERS:
# @(phy,npart) see note above
# RETURNS:
# (matrix) matrix with Nnodes+ntaxa rows and 2*npart columns. All rows from 1 to ntaxa are set to 0, exept for columns 1 and npart+1, set to one. All rows from ntaxa to the end are set to NAs
# DEPENDENCIES:
# none
# PURPOSE:
# Initialise the matrix for function partitionsNumber
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
##
init.partitionsNumber <- function(phy,npart){
  ntaxa <- length(phy$tip.label)
  nbrCompatiblePartitions <- matrix(NA, nrow = 1 + nrow(phy$edge), ncol = 2*npart)
  nbrCompatiblePartitions[1:ntaxa, ] <- 0
  nbrCompatiblePartitions[1:ntaxa, 1] <- rep(1,ntaxa)
  nbrCompatiblePartitions[1:ntaxa, npart + 1] <- rep(1,ntaxa)
  return(nbrCompatiblePartitions)
}

##
# update.partitionsNumber.bin (daughtersParams, ...)
# PARAMETERS:
# @daughtersParams (matrix) matrix with 2*npart columns. Each row contains the result of partitionsNumber for the children of a given node
# RETURNS:
# (vector) vector of size 2*npart. For k=1:npart it contains the number of partitions with k groups compatible with the sub tree starting at the current node ; and for k=(npart+1):2*npart, it contains the number of "marqued" partitions with (k-npart) groups compatible with the sub tree starting at the current node.
# DEPENDENCIES:
# none
# PURPOSE:
# Update
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
##
update.partitionsNumber.bin <- function(daughtersParams, ...){
  npart <- dim(daughtersParams)[2]%/%2
  if (npart == 1) return( c(1,1) )
  nbrComp <- rep(NA, npart)
  nbrMarkedComp <- rep(NA, npart)
  nbrComp[1] <- 1; nbrMarkedComp[1] <- 1 # Always 1 way for one color
  left_nbrComp <- daughtersParams[1, 1:npart]
  right_nbrComp <- daughtersParams[2, 1:npart]
  left_nbrMarkedComp <- daughtersParams[1, (npart+1):(2*npart)]
  right_nbrMarkedComp <- daughtersParams[2, (npart+1):(2*npart)]
  for (k in 2:npart){
    nbrComp[k] <- sum( left_nbrComp[1:(k-1)] * rev(right_nbrComp[1:(k-1)]) ) + sum( left_nbrMarkedComp[1:k] * rev(right_nbrMarkedComp[1:k]) )
    nbrMarkedComp[k] <- sum( left_nbrComp[1:(k-1)] * rev(right_nbrMarkedComp[1:(k-1)]) ) + sum( right_nbrComp[1:(k-1)] * rev(left_nbrMarkedComp[1:(k-1)]) ) + sum( left_nbrMarkedComp[1:k] * rev(right_nbrMarkedComp[1:k]) )
  }
  nbrAdm <- c(nbrComp, nbrMarkedComp)
  return(nbrAdm)
}

##
#' @title Update formula in the general case
#'
#' @description
#' \code{update.partitionsNumber.gen} apply the actualisation formula to get
#' Nk and Ak of a node given its daughters.
#'
#' @details
#' Uses functions \code{sum.partitions} and \code{sum.simplex} to compute
#' the needed sums.
#'
#' @param daughtersParams : matrix with 2*npart columns. Each row contains
#' the result of partitionsNumber for the children of a given node
#'
#' @return vector of size 2*npart. For k=1:npart it contains the number of
#' partitions with k groups compatible with the sub tree starting at the
#' current node ; and for k=(npart+1):2*npart, it contains the number of
#' "marqued" partitions with (k-npart) groups compatible with the sub tree
#' starting at the current node.
#'
#'05/06/14 - Initial release
##
update.partitionsNumber.gen <- function(daughtersParams, ...){
  npart <- dim(daughtersParams)[2]%/%2
  if (npart == 1) return( c(1,1) )
  nbrComp <- rep(NA, npart)
  nbrMarkedComp <- rep(NA, npart)
  nbrComp[1] <- 1; nbrMarkedComp[1] <- 1 # Always 1 way for one color
  N_daughters <- daughtersParams[, 1:npart]
  A_daughters <- daughtersParams[, (npart+1):(2*npart)]
  p <- dim(N_daughters)[1]
  for (K in 2:npart){
    N <- N_daughters[,1:K]
    A <- A_daughters[,1:K]
    nbrComp[K] <- sum.simplex(N, K, p) + sum.partitions(A, N, K, p, 2)
    nbrMarkedComp[K] <- sum.partitions(A, N, K, p, 1)
  }
  nbrAdm <- c(nbrComp, nbrMarkedComp)
  return(nbrAdm)
}

##
# extract.partitionsNumber (nbrCompatiblePartitions, node=attr(nbrCompatiblePartitions, "ntaxa")+1, npart=attr(nbrCompatiblePartitions, "npart"), marqued=FALSE)
# PARAMETERS:
# @nbrCompatiblePartitions (matrix) return of the function partitionsNumber
# @node (int) : node from which to start (default value : root)
# @npart (int) : number of groups to be maximized (default value : npart from the previous partitionsNumber call)
# @marqued (bool) if true, return the number of marqued partitions
# RETURNS:
# (int) Number of partitions with npart groups compatible with the tree and the shift process, for sub-tree starting at node node.
# DEPENDENCIES:
# partitionsNumber
# PURPOSE:
# Extract the values wanted from the raw result of function partitionsNumber
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
# 27/05/14 - Add "marqued"
##
extract.partitionsNumber <- function(nbrCompatiblePartitions, node=attr(nbrCompatiblePartitions, "ntaxa")+1, npart=attr(nbrCompatiblePartitions, "npart"), marqued=FALSE){
  if (marqued) return(nbrCompatiblePartitions[node,attr(nbrCompatiblePartitions, "npart")+npart])
  return(nbrCompatiblePartitions[node,npart])
}
##################################
## Some small technical functions
##################################
##
#' @title Product of elements of a matrix
#'
#' @description
#' \code{prod.index} return the procuct of choosen elements of a matrix.
#'
#' @details
#' This function is to be used in \code{sum.simplex} to be applied to all the
#' elements given by \code{xsimplex}.
#' Performs the product : X[1,Id[1]]*X[2,Id[2]]*...*X[p,Id[p]] if all the
#' elements of Id are positive. Otherwise, just return 0.
#'
#' @param X a martix with p rows and K column. Each row contains the number of
#' partition in 1<=k<=K groups for one of the p children of a given node
#' @param Id a vector of length p, result of the function \code{xsimplex}.
#'
#' @return double : the result of the product.
#'
#'05/06/14 - Initial release
##
prod.index <- function(X,Id){
  if (0 %in% Id) return(0) # Indice hors limites
  return(prod(diag(X[,Id])))
}

##
#' @title Sum on a simplex
#'
#' @description
#' \code{sum.simplex} returns the sum on k_1+...+k_p=K, k_i>0 of the products
#' of NN[i,k_i].
#'
#' @details
#' This function uses \code{xsimplex} to perform the product of NN[i,k_i] for
#' all combination of k_i such that k_1+...+k_p=K, k_i>0, using function
#' \code{prod.index}. Then sum all the products.
#'
#' @param NN a martix with p rows and K column. Each row contains the number of
#' partition in 1<=k<=K groups for one of the p children of a given node
#' @param K an integer. Tne number of groups wanted.
#' @param p an integer. The number of daughters of a node.
#'
#' @return double : the result of the sum.
#'
#'05/06/14 - Initial release
##
sum.simplex <- function (NN, K, p) {
  return(sum(xsimplex(p, K, fun=prod.index, simplify=TRUE, X=NN)))
}

##
#' @title Sum on a simplex
#'
#' @description
#' \code{sum.prod.comb} returns the sum on k_1+...+k_p=K+|I|-1, k_i>0 of the
#' products of prod(A[i,k_i], i in I)*prod(N[i,k_i], i not in I).
#'
#' @details
#' This function uses \code{sum.simplex} to perform the wanted sum on a ad-hoc
#' matrix, combination of rows of A and N.
#'
#' @param I a vector of integers representing a subset of [1,p]
#' @param A a martix with p rows and K column. Each row contains the number of
#' marqued partition in 1<=k<=K groups for one of the p children of a
#' given node.
#' @param N a martix with p rows and K column. Each row contains the number of
#' partition in 1<=k<=K groups for one of the p children of a given node
#' @param K an integer. Tne number of groups wanted.
#' @param p an integer. The number of daughters of a node.
#'
#' @return double : the result of the sum.
#'
#'05/06/14 - Initial release
##
sum.prod.comb <- function(I, A, N, K, p){
  KK <- K+length(I)-1
  AN <- matrix(NA, nrow=p, ncol=K)
  AN[I,] <- A[I,]
  AN[-I,] <- N[-I,]
  return(sum.simplex(AN, KK, p))
}

##
#' @title Sum on subsets of a given cardinal.
#'
#' @description
#' \code{sum.partitions.cardFixed} returns the sum on I subset of [1,p], |I|
#' fixed, of the sums computed by \code{sum.prod.comb}.
#'
#' @details
#' This function uses \code{combn} to enumarate all the subsets I of [1,p] of
#' a given cardinal, and then performs the wanted sum on these subsets.
#'
#' @param A a martix with p rows and K column. Each row contains the number of
#' marqued partition in 1<=k<=K groups for one of the p children of a
#' given node.
#' @param N a martix with p rows and K column. Each row contains the number of
#' partition in 1<=k<=K groups for one of the p children of a given node
#' @param K an integer. Tne number of groups wanted.
#' @param p an integer. The number of daughters of a node.
#' @param cardI an integer. The cardinal of the subset wanted.
#'
#' @return double : the result of the sum.
#'
#'05/06/14 - Initial release
##
sum.partitions.cardFixed <- function(A, N, K, p, cardI){
  return(sum(combn(p, cardI, fun=sum.prod.comb, simplify=TRUE, A=A, N=N, K=K, p=p)))
}

##
#' @title Sum on all subsets.
#'
#' @description
#' \code{sum.partitions} returns the sum on all I subset of [1,p], with |I|>m.
#'
#' @details
#' This function applies \code{sum.partitions.cardFixed} to all integer between
#' m and p, and sum the results.
#'
#' @param A a martix with p rows and K column. Each row contains the number of
#' marqued partition in 1<=k<=K groups for one of the p children of a
#' given node.
#' @param N a martix with p rows and K column. Each row contains the number of
#' partition in 1<=k<=K groups for one of the p children of a given node
#' @param K an integer. Tne number of groups wanted.
#' @param p an integer. The number of daughters of a node.
#' @param m an integer. The minimum cadinal of a subset allowed.
#'
#' @return double : the result of the sum.
#'
#'05/06/14 - Initial release
##
sum.partitions <- function(A, N, K, p, m) {
  return(sum(sapply(m:p, function(x) sum.partitions.cardFixed(A,N,K,p,x))))
}

###############################################################################
## Complexity Break Point
###############################################################################

complexity_break_point <- function(n){
  return(floor((10 * n - 13 - sqrt(20 * n^2 - 20 * n + 9)) / 10) + 1)
}