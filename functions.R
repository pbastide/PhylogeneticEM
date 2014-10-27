# {Test}
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

#
{ ##
# Purpose:           Define functions used in main.R
# Called From:       main.R
# Author:            MM, PB
# Notes:             --
# Revisions:         Created: 16/05/14
#                    Last change: 28/05/14 by PB
} ##
#

{ ##
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
} ##
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

## Simulate ###################################################################
###############################################################################

{ ##
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
} ##
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

{ ##
# simulate (phylo, root.state = list(random=FALSE,value.root,exp.root,var.root), process = c("BM", "OU"), shifts = list(edges=NULL,values=NULL,relativeTimes=NULL), ...)
# PARAMETERS:
#            @phylo (tree) Input tree
#            @root.state (list) state of the root :
#                 random : random state (TRUE) or deterministic state (FALSE)
#                 value.root : if deterministic, value of the character at the root
#                 exp.root : if random, expectation of the character at the root
#                 var.root : if random, variance of the character at the root
#            @process (string) Random process to simulate. Possible values :
#                 "BM" : Brownian Motion
#                 "OU" : Ornstein-Uhlenbeck
#            @shifts (list) position and values of the shifts :
#                 edges : vector of id of edges where the shfts are
#                 values : vector of same dimension of values of the shifts on the edges
#                 relativeTimes : vector of same dimension of relative time of the shift form the parent node of edges
#            @... other parameters : parameters of the process :
#                 for BM : variance
#                 for OU : variance, optimal.value, selection.strength
# RETURNS:
#            (matrix) Each row i contains a vector of values computed for node i :
#                       BM : simulated state, expected value
#                       OU : simulated state, expected value, optimal value
# DEPENDENCIES:
#            init, update, correspondanceEdges (, extract.simulate), check.selection.strength
# PURPOSE:
#            Simulate a random process on the tree phylo. Return the simulated and expected  values of all tips and nodes.
# NOTES:
#            none
# REVISIONS:
#            16/05/14 - Initial release
#            20/05/14 - Gestion of edges (function correspondanceEdges)
#            21/05/14 - extraction of the recursion
#            16/06/14 - check.selection.strength
} ##
simulate <- function(phylo, process = c("BM", "OU"), root.state = list(random=FALSE, stationary.root=FALSE, value.root, exp.root, var.root), shifts = list(edges=NULL,values=NULL,relativeTimes=NULL), eps=10^(-6), selection.strength=NULL, variance=NULL, optimal.value=NULL) {
  phy <- reorder(phylo, order = "cladewise")
  ## Trace edges
  shifts_ordered <- shifts
  shifts_ordered$edges <- correspondanceEdges(edges=shifts$edges,from=phylo,to=phy)
  ## Set branch stochastic process
  process <- match.arg(process)
  process <- check.selection.strength(process=process, selection.strength=selection.strength, eps=eps) # if OU, check if selection.strength is not too low.
  init  <- switch(process, 
                    BM = init.simulate.BM, 
                    OU = init.simulate.OU)
  updateDown  <- switch(process, 
                        BM = update.simulate.BM, 
                        OU = update.simulate.OU)
  ## Check root
  root.state <- test.root.state(root.state=root.state, 
                                process=process, 
                                eps=eps, 
                                variance=variance,
                                selection.strength=selection.strength,
                                optimal.value=optimal.value)
  ## Initialisation and setting root state
  ntaxa <- length(phy$tip.label)
  paramSimu <- init(phy=phy,
                    root.state=root.state, 
                    optimal.value=optimal.value)
  ## Tree recursion
  paramSimu <- recursionDown(phy=phy, 
                             params=paramSimu, 
                             updateDown=updateDown, 
                             shifts=shifts_ordered, 
                             variance=variance,
                             eps=eps, 
                             selection.strength=selection.strength)
  attr(paramSimu, "ntaxa")  <- ntaxa
  return(paramSimu)
}

{ ##
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
} ##
correspondanceEdges <- function(edges,from,to){
  mm <- match(from$edge[,2],to$edge[,2])
  newEdges <- mm[edges]
  return(newEdges)
}

{ ##
# init.simulate.StateAndExp (phy, root.state = list(random=FALSE,value.root,exp.root,var.root))
# PARAMETERS:
#            @phy (tree) Input tree
#            @root.state (list) state of the root : (see note above)
# RETURNS:
#            (matrix) Matrix with as many rows as the cumulated number of tips and nodes in the tree phy, and 2 columns. Each row 1<=i<=ntaxa contains a vector of values computed for tip i. All other rows are set to NAs
# DEPENDENCIES:
#            none
# PURPOSE:
#            Initialize the matrix used in simulate, with simulated and expected values of the root node.
# NOTES:
#            none
# REVISIONS:
#            16/05/14 - Initial release
} ##
init.simulate.StateAndExp <- function(phy,root.state){
  ntaxa <- length(phy$tip.label)
  paramSimu <- matrix(NA,nrow=1 + nrow(phy$edge),ncol=2)
  if (!root.state$random) { # The root is not random
    paramSimu[ntaxa + 1,]  <- c(root.state$value.root,
                                root.state$value.root)
  } else { # The value of the root is random N(exp.root,var.root)
    paramSimu[ntaxa + 1,]  <- c(root.state$exp.root + sqrt(root.state$var.root)*rnorm(1),
                                root.state$exp.root)
  }
  return(paramSimu)
}

{ ##
# init.simulate.BM (phy, root.state = list(random=FALSE,value.root,exp.root,var.root))
# PARAMETERS:
#            @phy (tree) Input tree
#            @root.state (list) state of the root (see note above)
#            @... other parameters : parameters of the process (see note above)
# RETURNS:
#            (matrix) Matrix with as many rows as the cumulated number of tips and nodes in the tree phy, and 2 columns. Row number ntaxa+1 (root) contains a vector of values computed for the root node. All other rows are set to NAs
# DEPENDENCIES:
#            initStateAndExp
# PURPOSE:
#            Initialize the matrix used in simulate, with simulated and expected values of the root node.
# NOTES:
#            none
# REVISIONS:
#            16/05/14 - Initial release
} ##
init.simulate.BM <- function(phy,root.state,...){
  return(init.simulate.StateAndExp(phy,root.state))
}

{ ##
# init.simulate.OU (phy, root.state = list(random=FALSE,value.root,exp.root,var.root))
# PARAMETERS:
#            @phy (tree) Input tree
#            @root.state (list) state of the root (see note above)
#            @... other parameters : parameters of the process (see note above)
# RETURNS:
#            (matrix) Matrix with as many rows as the cumulated number of tips and nodes in the tree phy, and 2 columns. Row number ntaxa+1 (root) contains a vector of values computed for the root node. All other rows are set to NAs
# DEPENDENCIES:
#            initStateAndExp
# PURPOSE:
#            Initialize the matrix used in simulate, with simulated and expected values of the root node, and optimal value of the root node 
# NOTES:
#            none
# REVISIONS:
#            16/05/14 - Initial release
} ##
init.simulate.OU <- function(phy, root.state, optimal.value, ...){
  paramSimu <- init.simulate.StateAndExp(phy,root.state)
  ntaxa <- length(phy$tip.label)
  beta <- rep(NA, length = 1 + nrow(phy$edge)) # selection strength
  beta[ntaxa + 1] <- optimal.value
  paramSimu <- cbind(paramSimu,beta)
}

{ ##
# update.simulate.BM (edgeNbr, ancestral, length, shifts, variance, ...)
# PARAMETERS:
#            @edgeNbr (int) Number of the edge considered
#            @ancestral (vector) Computed vector for the parental node
#            @length (double) Length of the current edge
#            @shifts (list) position and values of the shifts (see note above)
#            @variance (double) variance of the BM
# RETURNS:
#            (vector) simulated state, expected value for the daughter node
# DEPENDENCIES:
#            none
# PURPOSE:
#            Update the process in the case of the BM
# NOTES:
#            none
# REVISIONS:
#            16/05/14 - Initial release
} ##
update.simulate.BM <- function(edgeNbr, ancestral, length, shifts, variance, ...){
  shiftsIndex <- which(shifts$edges==edgeNbr) # If no shifts = NULL, and sum = 0
  return(c( ancestral[1] + sum(shifts$values[shiftsIndex]) + sqrt(length*variance) * rnorm(1),
            ancestral[2] + sum(shifts$values[shiftsIndex]) ) )
}

{ ##
# update.simulate.OU (edgeNbr, ancestral, length, shifts, variance, selection.strength, eps=10^(-6), ...)
# PARAMETERS:
#            @(edgeNbr, ancestral, length, shifts, variance, selection.strength) : see note above
#            @eps (double) : if the selection strenght is smaller than eps, simulate according to a BM instead of an OU
# RETURNS:
#            (vector) simulated state, expected value, optimal value for the daughter node
# DEPENDENCIES:
#            updateBM
# PURPOSE:
#            Update the process in the case of the OU
# NOTES:
#            none
# REVISIONS:
#            16/05/14 - Initial release
} ##
update.simulate.OU <- function(edgeNbr, ancestral, length, shifts, variance, selection.strength, ...){
  shiftsIndex <- which(shifts$edges==edgeNbr) # If no shifts = NULL, and sum = 0
  beta <- ancestral[3] + sum(shifts$values[shiftsIndex])
  ee <- exp(-selection.strength*length)
  ss <- sum(shifts$values[shiftsIndex]*( 1-exp( -selection.strength * length * (1-shifts$relativeTimes[shiftsIndex]) ) ))
  SimExp <- c( ancestral[3]*(1-ee) + ancestral[1]*ee + ss + sqrt(variance*(1-ee^2)/(2*selection.strength))*rnorm(1),
               ancestral[3]*(1-ee) + ancestral[2]*ee + ss )
  return(c(SimExp,beta))
}

{ ##
# extract.simulate (paramSimu, where=c("tips","nodes"), what=c("states","expectations"))
# PARAMETERS:
#            @paramSimu (matrix) return of the function simulate
#            @where (string) : where to extract the values : at the "tips" or the internal "nodes" ?
#            @what (string) : which value to extract : the simulated "states" or the "expectations" ?
# RETURNS:
#            (vector) values choosen for nodes/tips choosen
# DEPENDENCIES:
#            simulate
# PURPOSE:
#            Extract the values wanted from the raw result of function simulate.
# NOTES:
#            none
# REVISIONS:
#            16/05/14 - Initial release
#            28/05/14 - Add optimal.values
#            02/06/14 - Case where optimal.value is asked for a BM (retrun NULL)
} ##
extract.simulate <- function(paramSimu, where=c("tips","nodes"), what=c("states","expectations","optimal.values")){
  where <- match.arg(where)
  what <- match.arg(what)
  ntaxa <- attr(paramSimu,"ntaxa")
  if (where=="tips") {
    rows <- 1:ntaxa
  } else if (where=="nodes") {
    rows <- (ntaxa+1):nrow(paramSimu)
  }
  if (what=="states") {
    col <- 1
  } else if (what=="expectations") {
    col <- 2
  } else if (what=="optimal.values") {
    col <- 3
  }
  if (col > dim(paramSimu)[2] ){
    return(NULL) # Case of optimal.value asked for a BM
  } else {
    return(paramSimu[rows,col])
  }
}

{ ##
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
} ##
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

{ ##
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
} ##
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


### ParcimonyNumber############################################################
###############################################################################

{ ##
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
} ##
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

{ ##
# parcimonyNumber (phylo,clusters=rep(1,length(phy$tip.label)))
# PARAMETERS:
#            @phylo (tree) imput tree
#            @clusters (vector) : vector indicating the clusters of each tip
# RETURNS:
#            (matrix) matrix with ncluster columns and Nnodes+ntaxa rows. Each row i contains the numbers of parcimonious reconstruction of the subtree under the node i, if starting with state j (column)
# DEPENDENCIES:
#            init.parcimonyNumber, update.parcimonyNumber (, extrac.parcimonyNumber)
# PURPOSE:
#            Find the number of parcimonious repartition of the shift that can lead to a given clustering
# NOTES:
#            none
# REVISIONS:
#            19/05/14 - Initial release
#            21/05/14 - Gestion of case only one cluster. Extraction of the recursion.
} ##
parcimonyNumber <- function(phylo,clusters=rep(1,length(phy$tip.label))){
  phy <- reorder(phylo,"postorder")
  ntaxa <- length(phy$tip.label)
  ## Initialization
  nbrReconstructions <- init.parcimonyNumber(phy,clusters)
  ## Tree recursion
  nbrReconstructions <- recursionUp(phy, nbrReconstructions, update.parcimonyNumber)
  attr(nbrReconstructions, "ntaxa")  <- ntaxa
  return(nbrReconstructions)
}

{ ##
# init.parcimonyNumber (phy,clusters)
# PARAMETERS:
#            @(phy,clusters) see note above
# RETURNS:
#            (matrix) matrix with ncluster columns and Nnodes+ntaxa rows. For each tip i, row i is the vector Ind(i=j). All other rows are set to NAs
# DEPENDENCIES:
#            none
# PURPOSE:
#            Initialise the matrix for function parcimonyNumber
# NOTES:
#            none
# REVISIONS:
#            19/05/14 - Initial release
} ##
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

{ ##
# update.parcimonyNumber (daughtersNbr)
# PARAMETERS:
#            @daughtersNbr (matrix) matrix with ncluster columns. Each row contains the result of parcimonyNumber for the children of a given node 
# RETURNS:
#            (vector) vector containing the number of parcimonious reconstructions from the current node if starting with cluster j
# DEPENDENCIES:
#            none
# PURPOSE:
#            Update
# NOTES:
#            none
# REVISIONS:
#            19/05/14 - Initial release
} ##
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

{ ##
# extract.parcimonyNumber (nbrReconstructions,node=attr(nbrReconstructions, "ntaxa")+1)
# PARAMETERS:
#            @nbrReconstructions (matrix) return of the function parcimonyNumber
#            @node (int) : node from which to start (default value : root)
# RETURNS:
#            (int) Number of parcimonious reconstructions for the subtree starting at node node.
# DEPENDENCIES:
#            parcimonyNumber
# PURPOSE:
#            Extract the values wanted from the raw result of function parcimonyNumber
# NOTES:
#            none
# REVISIONS:
#            19/05/14 - Initial release
} ##
extract.parcimonyNumber <- function(nbrReconstructions,node=attr(nbrReconstructions, "ntaxa")+1){
  if (is.vector(nbrReconstructions)){
    return(1)
  } else{
    return(sum(nbrReconstructions[node,]))
  }
}

## partitionsNumber ###########################################################
###############################################################################

{ ##
# partitionsNumber (phylo, npart)
# PARAMETERS:
#            @phylo (tree) imput tree
#            @npart (int) : number of partitions of the tips allowed
# RETURNS:
#            (matrix) matrix with Nnodes+ntaxa rows and 2*npart columns. Each column contains two vectors : for k=1:npart it contains the number of partitions with k groups compatible with the tree and the shift process; and for k=(npart+1):2*npart, it contains the number of "marqued" partitions with (k-npart) groups compatible with the tree and the shift process.
# DEPENDENCIES:
#            init.partitionsNumber, update.partitionsNumber.bin (, extract.partitionsNumber)
# PURPOSE:
#            Find the number of partitions with npart groups of the tips that are compatible with the structure of the tree and the shift process
# NOTES:
#            Only works for rooted, binary trees.
# REVISIONS:
#            26/05/14 - Initial release
} ##
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
  attr(nbrCompatiblePartitions, "ntaxa")  <- ntaxa
  attr(nbrCompatiblePartitions, "npart")  <- npart
  return(nbrCompatiblePartitions)
}

{ ##
# init.partitionsNumber (phy,clusters)
# PARAMETERS:
#            @(phy,npart) see note above
# RETURNS:
#            (matrix) matrix with Nnodes+ntaxa rows and 2*npart columns. All rows from 1 to ntaxa are set to 0, exept for columns 1 and npart+1, set to one. All rows from ntaxa to the end are set to NAs
# DEPENDENCIES:
#            none
# PURPOSE:
#            Initialise the matrix for function partitionsNumber
# NOTES:
#            none
# REVISIONS:
#            26/05/14 - Initial release
} ##
init.partitionsNumber <- function(phy,npart){
  ntaxa <- length(phy$tip.label)
  nbrCompatiblePartitions <- matrix(NA, nrow = 1 + nrow(phy$edge), ncol = 2*npart)
  nbrCompatiblePartitions[1:ntaxa, ] <- 0
  nbrCompatiblePartitions[1:ntaxa, 1] <- rep(1,ntaxa)
  nbrCompatiblePartitions[1:ntaxa, npart + 1] <- rep(1,ntaxa)
  return(nbrCompatiblePartitions)
}

{ ##
# update.partitionsNumber.bin (daughtersParams, ...)
# PARAMETERS:
#            @daughtersParams (matrix) matrix with 2*npart columns. Each row contains the result of partitionsNumber for the children of a given node 
# RETURNS:
#            (vector) vector of size 2*npart. For k=1:npart it contains the number of partitions with k groups compatible with the sub tree starting at the current node ; and for k=(npart+1):2*npart, it contains the number of "marqued" partitions with (k-npart) groups compatible with the sub tree starting at the current node.
# DEPENDENCIES:
#            none
# PURPOSE:
#            Update
# NOTES:
#            none
# REVISIONS:
#            26/05/14 - Initial release
} ##
update.partitionsNumber.bin <- function(daughtersParams, ...){
  npart <- dim(daughtersParams)[2]%/%2
  if (npart == 1)  return( c(1,1) )
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

{ ##
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
} ##
update.partitionsNumber.gen <- function(daughtersParams, ...){
  npart <- dim(daughtersParams)[2]%/%2
  if (npart == 1)  return( c(1,1) )
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

{ ##
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
} ##
prod.index <- function(X,Id){
  if (0 %in% Id) return(0) # Indice hors limites
  return(prod(diag(X[,Id])))
}

{ ##
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
} ##
sum.simplex <- function (NN, K, p) {
  return(sum(xsimplex(p, K, fun=prod.index, simplify=TRUE, X=NN)))
}

{ ##
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
} ##
sum.prod.comb <- function(I, A, N, K, p){
  KK <- K+length(I)-1
  AN <- matrix(NA, nrow=p, ncol=K)
  AN[I,] <- A[I,]
  AN[-I,] <- N[-I,]
  return(sum.simplex(AN, KK, p))
}

{ ##
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
} ##
sum.partitions.cardFixed <- function(A, N, K, p, cardI){
  return(sum(combn(p, cardI, fun=sum.prod.comb, simplify=TRUE, A=A, N=N, K=K, p=p)))
}

{ ##
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
} ##
sum.partitions <- function(A, N, K, p, m) {
  return(sum(sapply(m:p, function(x) sum.partitions.cardFixed(A,N,K,p,x))))
}

{ ##
# extract.partitionsNumber (nbrCompatiblePartitions, node=attr(nbrCompatiblePartitions, "ntaxa")+1, npart=attr(nbrCompatiblePartitions, "npart"), marqued=FALSE)
# PARAMETERS:
#            @nbrCompatiblePartitions (matrix) return of the function partitionsNumber
#            @node (int) : node from which to start (default value : root)
#            @npart (int) : number of groups to be maximized (default value : npart from the previous partitionsNumber call)
#            @marqued (bool) if true, return the number of marqued partitions
# RETURNS:
#            (int)  Number of partitions with npart groups compatible with the tree and the shift process, for sub-tree starting at node node.
# DEPENDENCIES:
#            partitionsNumber
# PURPOSE:
#            Extract the values wanted from the raw result of function partitionsNumber
# NOTES:
#            none
# REVISIONS:
#            26/05/14 - Initial release
#            27/05/14 - Add "marqued"
} ##
extract.partitionsNumber <- function(nbrCompatiblePartitions, node=attr(nbrCompatiblePartitions, "ntaxa")+1, npart=attr(nbrCompatiblePartitions, "npart"), marqued=FALSE){
  if (marqued) return(nbrCompatiblePartitions[node,attr(nbrCompatiblePartitions, "npart")+npart])
  return(nbrCompatiblePartitions[node,npart])
}

{ ##
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
} ##
rtree.sym <- function(n){
  tree <- "A"
  for (k in 1:n) {
    tree <- paste("(", tree, ",", tree, ")", sep="")
  }
  return(read.tree(text=paste(tree, ";", sep="")))
}

{ ##
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
} ##
rtree.comb <- function(n){
  if (n == 1) return(read.tree(text="(A);"))
  tree <- "A"
  for (k in 2:n) {
    tree <- paste("(A,", tree, ")", sep="")
  }
  return(read.tree(text=paste(tree, ";", sep="")))
}

## estimateEM #################################################################
###############################################################################

{ ##
# estimateEM (phylo, Y_data, process=c("BM","OU"), tol=10^(-5),  method.variance=c("simple"), method.init=c("default"), nbr_of_shifts=0, ...)
# PARAMETERS:
#            @phylo (tree) imput tree
#            @Y_data (vector) : vector indicating the data at the tips
#            @process (string) Random process to simulate. Possible values :
#                 "BM" : Brownian Motion
#                 "OU" : Ornstein-Uhlenbeck
#            @tol (list) : tolerance for shutoff. TO BE DEFINED
#            @method.variance (string) : method to compute the variance :
#                 "simple" : invert the matrix Sigma_YY with function "solve" of R
#            @method.init : method to initialize the parameters
#                 "default" : put the parameters to a pre-defined value
#            @nbr_of_shifts : number of shifts wanted for the inference
#            @specialCase : boolean. If true, the inference of the parameters is done in the special case for the OU : alpha known, shifts on nodes, root in stationnay state
# RETURNS:
#            (list) list of parameters for the model fitted
# DEPENDENCIES:
#            init.EM, compute_E, compute_M, compute_times_ca, compute_dist_phy, shutoff.EM, update.parcimonyNumber (, ...)
# PURPOSE:
#            Run the EM algorithm
# NOTES:
#            Under construction
# REVISIONS:
#            21/05/14 - Initial release (non fonctionnal)
#            22/05/14 - Minimal "working" release
#            02/06/14 - OU in the special case
#            10/06/14 - Test of divergence
} ##
estimateEM <- function(phylo, 
                       Y_data, 
                       process = c("BM","OU"), 
                       tol=list(variance = 10^(-5), 
                                value.root = 10^(-5), 
                                exp.root = 10^(-5), 
                                var.root = 10^(-5),
                                selection.strength = 10^(-5)),  
                       Nbr_It_Max = 500, 
                       method.variance = c("simple"), 
                       method.init = c("default", "lasso"),
                       method.init.alpha = c("default", "estimation"),
                       nbr_of_shifts = 0,
                       random.root = TRUE,
                       stationnary.root = TRUE,
                       shifts_at_nodes = TRUE,
                       alpha_known = FALSE,
                       eps = 10^(-3),
                       known.selection.strength = 1,
                       init.selection.strength = 1,
                       max_selection.strength = 100,
                       use_sigma_for_lasso = TRUE,
                       max_triplet_number = 10000,
                       min_params=list(variance = 10^(-3), 
                                       value.root = -10^(3), 
                                       exp.root = -10^(3), 
                                       var.root = 10^(-3),
                                       selection.strength = 10^(-3)),
                       max_params=list(variance = 10^(3), 
                                       value.root = 10^(3), 
                                       exp.root = 10^(3), 
                                       var.root = 10^(3),
                                       selection.strength = 10^(3)),
                       var.init.root = 1,
                       methods.segmentation = c("max_costs_0", "lasso", "same_shifts", "same_shifts_same_values", "best_single_move"), ...){
  ## Check consistancy
  if (alpha_known && missing(known.selection.strength)) stop("The selection strength alpha is supposed to be known, but is not specified. Please add an argument known.selection.strength to the call of the function.")
  ## Choose process
  process <- match.arg(process)
  process <- check.selection.strength(process, known.selection.strength, eps)
  # specialCase <- stationnary.root && shifts_at_nodes && alpha_known
  compute_M  <- switch(process, 
                       BM = compute_M.BM,
                       OU = compute_M.OU(stationnary.root, shifts_at_nodes, alpha_known))
  shutoff.EM  <- switch(process, 
                        BM = shutoff.EM.BM,
                        OU = shutoff.EM.OU(stationnary.root, shifts_at_nodes, alpha_known))
  is.finite.params  <- switch(process, 
                              BM = is.finite.params.BM,
                              OU = is.finite.params.OU(stationnary.root, shifts_at_nodes, alpha_known))
  is.in.ranges.params  <- switch(process, 
                              BM = is.in.ranges.params.BM,
                              OU = is.in.ranges.params.OU(stationnary.root, shifts_at_nodes, alpha_known))
  compute_MaxCompleteLogLik <- switch(process, 
                                      BM = compute_MaxCompleteLogLik.BM,
                                      OU = compute_MaxCompleteLogLik.OU(stationnary.root, shifts_at_nodes))
  conditional_expectation_log_likelihood <- switch(process, 
                                                   BM = conditional_expectation_log_likelihood.BM,
                                                   OU = conditional_expectation_log_likelihood.OU(stationnary.root, shifts_at_nodes))
  ## init alpha
  method.init.alpha  <- match.arg(method.init.alpha)
  if (!stationnary.root && (method.init.alpha == "estimation")){
    method.init.alpha <- "default"
    warning("The estimation initialization of alpha does only work when the root is stationnary. The initialization is set to the default one.")
  }
  init.alpha<- switch(process, 
                      BM = init.alpha.BM,
                      OU = init.alpha.OU)
  init.alpha.gamma<- switch(process, 
                            BM = init.alpha.gamma.BM,
                            OU = init.alpha.gamma.OU)
  ## Choose method
  method.variance  <- match.arg(method.variance)
  compute_E  <- switch(method.variance, 
                       simple = compute_E.simple)
  compute_mean_variance  <- switch(method.variance, 
                                   simple = compute_mean_variance.simple)
  compute_log_likelihood  <- switch(method.variance, 
                                    simple = compute_log_likelihood.simple)
  method.init  <- match.arg(method.init)
  init.EM  <- switch(method.init, 
                     default = init.EM.default(process),
                     lasso = init.EM.lasso)
  method.init.alpha  <- match.arg(method.init.alpha)
  methods.segmentation <- match.arg(methods.segmentation, several.ok = TRUE)
  ## Initialization
  ntaxa <- length(phylo$tip.label)
  times_shared <- compute_times_ca(phylo)
  distances_phylo <- compute_dist_phy(phylo)
  init.a.g <- init.alpha.gamma(method.init.alpha)(phylo = phylo,
                            Y_data = Y_data,
                            nbr_of_shifts = nbr_of_shifts,
                            distances_phylo = distances_phylo,
                            init.selection.strength = init.selection.strength,
                            max_triplet_number = max_triplet_number,
                            known.selection.strength = known.selection.strength,
                            alpha_known = alpha_known,
                            init.var.root = var.init.root)
  init.var.root <- init.a.g$gamma_0
  if (!alpha_known) {
    init.selection.strength <- init.a.g$alpha_0
  } else {
    init.selection.strength <- known.selection.strength
  }
  params_init <- init.EM(phylo = phylo,
                         Y_data = Y_data,
                         process = process, 
                         times_shared = times_shared, 
                         distances_phylo = distances_phylo, 
                         nbr_of_shifts = nbr_of_shifts, 
                         selection.strength.init = init.selection.strength, 
                         random.init = random.root,
                         stationnary.root.init = stationnary.root,
                         use_sigma = use_sigma_for_lasso,
                         method.init.alpha = method.init.alpha,
                         var.root.init = init.var.root,
                         ...)
  params <- params_init
  params$root.state <- test.root.state(root.state=params$root.state, 
                                  process=process, 
                                  optimal.value=params$optimal.value,
                                  variance=params$variance, 
                                  selection.strength=params$selection.strength)
  attr(params, "ntaxa")  <- ntaxa
  params_old <- NULL
  ## Iteration
  Nbr_It <- 0
  params_history <- vector("list")#, Nbr_It_Max)
#   CLL_history <- NULL
  number_new_shifts <- NULL
  while ( Nbr_It == 0 || # Initialisation
            ( !shutoff.EM(params_old,params,tol) && # Shutoff
              is.in.ranges.params(params, min=min_params, max=max_params) && # Divergence ?
              Nbr_It < Nbr_It_Max ) ) { # Nbr of iteration
    ## Actualization
    Nbr_It <- Nbr_It + 1
    params_old <- params
    ## Log likelihood
    moments <- compute_mean_variance(phylo = phylo,
                                     times_shared = times_shared,
                                     distances_phylo = distances_phylo,
                                     process = process,
                                     params_old = params_old)
    log_likelihood <- compute_log_likelihood(phylo = phylo,
                                             Y_data = Y_data,
                                             sim = moments$sim,
                                             Sigma = moments$Sigma,
                                             Sigma_YY_inv = moments$Sigma_YY_inv)
    attr(params_old, "log_likelihood") <- log_likelihood
    ## E step
    conditional_law_X <- compute_E(phylo = phylo,
                                   Y_data = Y_data,
                                   sim = moments$sim,
                                   Sigma = moments$Sigma,
                                   Sigma_YY_inv = moments$Sigma_YY_inv)
    ## Log likelihood as the sum of conditional + entropy
#     H <- compute_entropy.simple(moments$Sigma, moments$Sigma_YY_inv)
#     CLL <- conditional_expectation_log_likelihood(phylo = phylo,
#                                           conditional_law_X = conditional_law_X, 
#                                           sigma2 = params_old$variance,
#                                           mu = params_old$root.state$exp.root,
#                                           shifts = params_old$shifts,
#                                           alpha = params_old$selection.strength)
#     log_likelihood_bis <- compute_log_likelihood_with_entropy.simple(CLL, H)
#     attr(params_old, "log_likelihood_bis") <- log_likelihood_bis
    ## Store params for history
    params_history[[paste(Nbr_It - 1, sep="")]] <- params_old
    ## M step
    params <- compute_M(phylo = phylo, 
                        Y_data = Y_data, 
                        conditional_law_X = conditional_law_X, 
                        nbr_of_shifts = nbr_of_shifts, 
                        random.root = random.root,
                        known.selection.strength = known.selection.strength,
                        alpha_old = params_old$selection.strength,
                        max_selection.strength = max_selection.strength,
                        eps = eps,
                        methods.segmentation = methods.segmentation,
                        beta_0_old = params_old$optimal.value,
                        shifts_old = params_old$shifts)
    attr(params, "ntaxa")  <- ntaxa
    ## Number of shifts that changed position ?
    number_new_shifts <- c(number_new_shifts,
                           sum(!(params$shifts$edges %in% params_old$shifts$edges)))
    ## Check that the M step rised the conditional expectation of the completed likelihood
#     CLL_old <- conditional_expectation_log_likelihood(phylo = phylo,
#                                     conditional_law_X = conditional_law_X, 
#                                     sigma2 = params_old$variance,
#                                     mu = params_old$root.state$exp.root,
#                                     shifts = params_old$shifts,
#                                     alpha = params_old$selection.strength)
#     names(CLL_old) <- "CCL_old"
#     CLL_new <- conditional_expectation_log_likelihood(phylo = phylo,
#                                     conditional_law_X = conditional_law_X, 
#                                     sigma2 = params$variance,
#                                     mu = params$root.state$exp.root,
#                                     shifts = params$shifts,
#                                     alpha = params$selection.strength)
#     names(CLL_new) <- "CCL_new"
#     if (CLL_old > CLL_new) { warning("The conditional expectation of the completed log likelihood decreased after step M !")}
#     CLL_history <- cbind(CLL_history, c(CLL_old, CLL_new))
#     attr(params, "MaxCompleteLogLik") <- CLL_new
  }
  result <- list(params = params, 
                 ReconstructedNodesStates = conditional_law_X$expectations[(ntaxa+1):length(conditional_law_X$expectations)], 
                 params_old = params_old, 
                 params_init = params_init,
                 params_history = params_history,
                 number_new_shifts = number_new_shifts)
#                  CLL_history = CLL_history

  attr(result, "Nbr_It") <- Nbr_It
  attr(result, "Divergence") <- !is.in.ranges.params(result$params, min=min_params, max=max_params) # TRUE if has diverged
  if (Nbr_It == Nbr_It_Max) warning(paste("The maximum number of iterations (Nbr_It_Max = ",Nbr_It_Max,") was reached.",sep=""))
  return(result)
}

{ ##
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
} ##
init.EM.default <- function(process){
  if (process=="BM"){
    return(init.EM.default.BM)
  } else if (process=="OU"){
    return(init.EM.default.OU)
  }
}

init.EM.default.BM <- function(variance.init=1, random.init=TRUE, value.root.init=0, exp.root.init=1, var.root.init=1, edges.init=NULL, values.init=NULL, relativeTimes.init=NULL, ...) {
  if (random.init) {
    value.root.init <- NA
  } else {
    exp.root.init <- NA
    var.root.init <- NA
  }
  params_init=list(variance=variance.init,
                   root.state=list(random=random.init,
                                   value.root=value.root.init,
                                   exp.root=exp.root.init,
                                   var.root=var.root.init),
                   shifts=list(edges=edges.init,
                               values=values.init,
                               relativeTimes=relativeTimes.init))
  return(params_init)
}

init.EM.default.OU <- function(variance.init=1, random.init=TRUE, stationary.root.init=TRUE, value.root.init=1, exp.root.init=1, var.root.init=1, edges.init=NULL, values.init=NULL, relativeTimes.init=NULL, selection.strength.init=1, optimal.value.init=0, ...) {
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
  params_init=list(variance=variance.init,
                   root.state=list(random=random.init,
                                   stationary.root=stationary.root.init,
                                   value.root=value.root.init,
                                   exp.root=exp.root.init,
                                   var.root=var.root.init),
                   shifts=list(edges=edges.init,
                               values=values.init,
                               relativeTimes=relativeTimes.init),
                   selection.strength=selection.strength.init,
                   optimal.value=optimal.value.init)
  return(params_init)
}

{ ##
#' @title Do a lasso regression with the number of non-zero variables fixed.
#'
#' @description
#' \code{lasso_regression_K_fixed} does the following regression :
#' ||Yp-Xp.delta|| + lambda |delta|_1 using the function \code{glmnet} of 
#' package \code{glmnet}, where delta is a vector representing the shifts 
#' occuring on the branches. It does a gauss lasso regression using function 
#' \code{lm} on top of it. This function is used in functions 
#' \code{init.EM.lasso}, \code{segmentation.OU.specialCase.lasso}, ...
#'
#' @details
#' lambda is choosen so that delta has the right number of non zero components.
#' If not possible, either temporaly raise the number of shifts and then select
#' only the shifts with the highest modulus, or if not possible, throw an error.
#'
#' @param Yp (transformed) data
#' @param Xp (transformed) matrix of regression
#' @param K number of non-zero components allowed
#' 
#' @return E0.gauss the intercept (value at the root)
#' @return shifts.gauss the list of shifts found on the branches
#'
#'06/10/14 - Initial release
} ##
lasso_regression_K_fixed <- function (Yp, Xp, K, intercept.penalty = FALSE ) {
  ## Penalty on the first coordinate = intercept : force first cooerdinate to be null
  excl <- NULL
  if (intercept.penalty) excl <- c(1)
  ## fit
  fit <- glmnet(x = 0 + Xp, y = Yp, alpha = 1, nlambda = 50000, dfmax = K + 1, exclude = excl)
  ## Find the lambda that gives the right number of ruptures
  K_2 <- K
  ## Check that lambda goes far enought
  if (K_2 > max(fit$df)) {
    fit <- glmnet(x = 0 + Xp, y = Yp, alpha = 1, dfmax = K + 1, exclude = excl)
  }
  ## If the right lambda does not exists, raise the number of shifts
  while (sum(fit$df == K_2) == 0 && K_2 < 500) {
    warning("During lasso regression, could not find the right lambda for the number of shifts K. Temporarly raised it to do the lasso regression, and furnishing the K largest coefficients.")
    K_2 <- K_2 + 1
  }
  ## If could not find the right lambda, do a default initialization
  if (sum(fit$df == K_2) == 0) {
    stop("Lasso Initialisation fail : could not find a satisfying number of shifts.")
  } else {
    delta <- coef(fit, s = fit$lambda[min(which(fit$df == K_2))])
    E0 <- delta[1]; # Intercept
    delta <- delta[-1];
    ## Gauss lasso
    projection <- which(delta != 0)
    Xproj <- 0 + Xp[, projection]
    fit.gauss <- lm(Yp ~ Xproj)
    delta.gauss <- rep(0, dim(Xp)[2])
    E0.gauss <- coef(fit.gauss)[1]; names(E0.gauss) <- NULL
    delta.gauss[projection] <- coef(fit.gauss)[-1]
    ## If we had to raise the number of shifts, go back to the initial number, taking the K largest shifts
    edges <- order(-abs(delta.gauss))[1:K]
    delta.gauss.final <- rep(0, length(delta.gauss))
    delta.gauss.final[edges] <- delta.gauss[edges]
    shifts.gauss <- shifts.vector_to_list(delta.gauss.final);
    return(list(E0.gauss = E0.gauss, shifts.gauss = shifts.gauss))
  }
}

{ ##
#' @title Initialisation of the shifts using Lasso.
#'
#' @description
#' \code{init.EM.lasso} does the following regression :
#' ||Y_data-T.delta||_(Sigma_YY^(-1)) + lambda |delta|_1
#' using the function \code{glmnet} of package \code{glmnet}, throught function
#' \code{lasso_regression_K_fixed}. T is the incidence matrix of the tree, and 
#' delta the vectorial representation of the shifts (see functions 
#' \code{incidence.matrix} and \code{shifts.list_to_vector} for further details).
#'
#' @details
#' A cholesky decomposition of function Sigma_YY^(-1) is used.
#' lambda is choosen so that delta has the right number of non zero components.
#'
#' @param Y_data data at the tips.
#' @param times_shared (matrix) : times of shared ancestry, result of function 
#' \code{compute_times_ca}.
#' @param distances_phylo (matrix) : phylogenetics distance, result of function 
#' \code{compute_dist_phy}
#' @param nbr_of_shifts number of shifts used in the EM algorithm
#' 
#' @return params_init the list of initial parameters to be used, in the right
#'  format.
#'
#'18/06/14 - Initial release
#'06/10/14 - Externalization of function lasso
} ##
init.EM.lasso <- function(phylo, Y_data, process, times_shared, distances_phylo, nbr_of_shifts, use_sigma=TRUE, variance.init=1, random.init=TRUE, value.root.init=0, exp.root.init=1, var.root.init=1, edges.init=NULL, values.init=NULL, relativeTimes.init=NULL, selection.strength.init=1, optimal.value.init=0, ...) {
  init.EM.default <- init.EM.default(process)
  ## Choose the norm :
  if (use_sigma) {
    # Choose process
    compute_variance_covariance  <- switch(process, 
                                           BM = compute_variance_covariance.BM,
                                           OU = compute_variance_covariance.OU)
    # Initialize Sigma with default parameters
    params.default <- init.EM.default(selection.strength.init=selection.strength.init, random.init=random.init, stationnary.root.init=stationnary.root.init, var.root.init = var.root.init)
    Sigma <- compute_variance_covariance(times_shared=times_shared, distances_phylo=distances_phylo, params_old=params.default)
    Sigma_YY <- extract.variance_covariance(Sigma, what="YY")
    # Cholesky
    Sig_chol <- chol(Sigma_YY)
    Sig_chol_inv <- t(solve(Sig_chol)) # Sigma_YY_inv = t(Sig_chol_inv)%*%Sig_chol_inv
    # Transform Y_data and T
    Tr <- incidence.matrix(phylo)
    Tp <- Sig_chol_inv%*%Tr
    Yp <- Sig_chol_inv%*%Y_data
  } else {
    # Return untransformed Y_data and T
    Tp <- incidence.matrix(phylo)
    Yp <- Y_data
  }
  ## Fit
  fit <- try(lasso_regression_K_fixed(Yp = Yp, Xp = Tp, K = nbr_of_shifts))
  if (inherits(fit, "try-error")) {
    warning("Lasso Initialisation fail : could not find a satisfying number of shifts. Proceeding to a default initialization.")
    return(init.EM.default(selection.strength.init=selection.strength.init, random.init=random.init, stationnary.root.init=stationnary.root.init, ...))
  } else { 
    E0.gauss <- fit$E0.gauss
    shifts.gauss <- fit$shifts.gauss
    params_init <- init.EM.default(value.root.init=E0.gauss[1], 
                                   exp.root.init=E0.gauss[1], 
                                   optimal.value.init=E0.gauss[1],
                                   edges.init=shifts.gauss$edges, 
                                   values.init=shifts.gauss$values, 
                                   relativeTimes.init=shifts.gauss$relativeTimes, 
                                   selection.strength.init=selection.strength.init, 
                                   random.init = random.init, 
                                   var.root.init = var.root.init, ...)
    return(params_init)
  }
#   fit.sig <- glmnet(x=0+Tp, y=Yp, alpha=1, nlambda=50000, dfmax=nbr_of_shifts+1)
#   ## Find the lambda that gives the right number of ruptures
#   nbr_of_shifts_2 <- nbr_of_shifts
#   ## Check that lambda goes far enought
#   if (nbr_of_shifts_2 > max(fit.sig$df)) {
#     fit.sig <- glmnet(x=0+Tp, y=Yp, alpha=1, dfmax=nbr_of_shifts+1)
#   }
#   ## If the right lambda does not exists, raise the number of shifts
#   while (sum(fit.sig$df==nbr_of_shifts_2) == 0 && nbr_of_shifts_2 < 500) {
#     nbr_of_shifts_2 <- nbr_of_shifts_2 + 1
#   }
#   ## If coudl not find the right lambda, do a default initialization
#   if (sum(fit.sig$df==nbr_of_shifts_2) == 0) {
#     warning("Lasso Initialisation fail : could not find a satisfying number of shifts. Proceeding to a default initialization.")
#     return(init.EM.default(selection.strength.init=selection.strength.init, random.init=random.init, stationnary.root.init=stationnary.root.init, ...))
#   } else {
#     delta.sig <- coef(fit.sig, s=fit.sig$lambda[min(which(fit.sig$df==nbr_of_shifts_2))])
#     E0.sig <- delta.sig[1]; # Intercept
#     delta.sig <- delta.sig[-1];
#     ## Gauss lasso
#     projection.sig <- which(delta.sig != 0)
#     Tproj <- 0 + Tp[, projection.sig]
#     fit.sig.gauss <- lm(Yp ~ Tproj)
#     delta.sig.gauss <- rep(0, dim(Tp)[2])
#     E0.sig.gauss <- coef(fit.sig.gauss)[1]; names(E0.sig.gauss) <- NULL
#     delta.sig.gauss[projection.sig] <- coef(fit.sig.gauss)[-1]
#     ## If we had to raise the number of shifts, go back to the initial number, taking the K largest shifts
#     edges <- order(-abs(delta.sig.gauss))[1:nbr_of_shifts]
#     delta.sig.gauss.final <- rep(0, length(delta.sig.gauss))
#     delta.sig.gauss.final[edges] <- delta.sig.gauss[edges]
#     shifts.init.sig.gauss <- shifts.vector_to_list(delta.sig.gauss.final);
#     ## Initialize parameters
#     params_init <- init.EM.default(value.root.init=E0.sig.gauss[1], 
#                                    exp.root.init=E0.sig.gauss[1], 
#                                    optimal.value.init=E0.sig.gauss[1],
#                                    edges.init=shifts.init.sig.gauss$edges, 
#                                    values.init=shifts.init.sig.gauss$values, 
#                                    relativeTimes.init=shifts.init.sig.gauss$relativeTimes, 
#                                    selection.strength.init=selection.strength.init, 
#                                    random.init = random.init, 
#                                    var.root.init = var.root.init, ...)
#     return(params_init)
#   }
}

init.alpha.BM <- function(method.init.alpha){
  return(function(...) return(0))
}

init.alpha.OU <- function(method.init.alpha){
  return(switch(method.init.alpha, 
                default = init.alpha.default,
                estimation = init.alpha.estimation))
}

init.alpha.gamma.BM <- function(method.init.alpha){
  return(function(var.init.root, ...) return(list(alpha_0 = 0,
                                                  gamma_0 = var.init.root)))
}

init.alpha.gamma.OU <- function(method.init.alpha){
  return(switch(method.init.alpha, 
                default = init.alpha.gamma.default,
                estimation = init.alpha.gamma.estimation))
}

init.alpha.default <- function(init.selection.strength, known.selection.strength, alpha_known, ...){
  if (alpha_known) {
    return(known.selection.strength)
  } else {
    return(init.selection.strength)
  }
}

{ ##
#' @title Initialisation the selection strength alpha using robust estimation
#'
#' @description
#' \code{init.alpha.estimation} fits (Y_i-Y_j)^2 ~ gamma^2(1-exp(-alpha*d_ij))
#' for all couples of tips (i,j) that have the same mean, i.e than are not
#' separated by a shift. Shifts are initialized thanks to a lasso
#' (function \code{init.EM.lasso}).
#'
#' @details
#' Function \code{nlrob} is used for the robust fit.
#'
#' @param phylo phylogenetic tree.
#' @param Y_data data at the tips.
#' @param nbr_of_shifts : number of shifts wanted
#' @param distances_phylo (matrix) : phylogenetics distance, result of function 
#' \code{compute_dist_phy}
#' @param nbr_of_shifts number of shifts used in the EM algorithm
#' 
#' @return params_init the list of initial parameters to be used, in the right
#'  format.
#'
#'10/07/14 - Initial release
} ##

init.alpha.estimation <- function(phylo, Y_data, nbr_of_shifts, distances_phylo, max_triplet_number, ...){
  ## Initialize a vector with the group of each tip
  tips_groups <- rep(0, length(phylo$tip.label))
  names(tips_groups) <- phylo$tip.label
  ## Initialize shifts by a lasso without sigma
  if (nbr_of_shifts > 0) {
    lasso <- init.EM.lasso(phylo=phylo, Y_data=Y_data, process="OU", nbr_of_shifts=nbr_of_shifts, use_sigma=FALSE)
    ## Roeorder phylo and trace edges
    phy <- reorder(phylo, order = "cladewise")
    edges_shifts <- correspondanceEdges(edges=lasso$shifts$edges,from=phylo,to=phy)
    ## Set groups of tips (one group = all the tips under a given shift)
    Tr <- incidence.matrix(phy)
    for (ed in order(edges_shifts)) { # Do it in order so that edges with higher numbers erase groups (edges closer from the tips)
      ed_sh <- edges_shifts[ed]
      tips_groups[phy$tip.label[Tr[,ed_sh]]] <- ed
    }
  } else {
    edges_shifts <- NULL
  }
  ## For each group, take all the triplets of tips to estimate the covariance sigma_ij
  cor_hat <- NULL # estimations from trilpets of pairs corelations
  square_diff <- NULL # (Y_i-Y_j)^2
  dists <- NULL # corresponding phylogenetic distances between pairs
  hat_gam <- rep(0, length(edges_shifts)+1)
  for (grp in 0:length(edges_shifts)) {
    tips <- which(tips_groups==grp)
    hat_gam[grp+1] <- var(Y_data[tips])
    #     if (length(tips) > 2){
    # #       ## Genrate all combinations
    # #       Combinations <- combn(x=tips, m=3)
    # #       Ncomb <- ncol(Combinations)
    # #       ## If too many of them, sample from them
    # #       if (Ncomb > max_triplet_number) {
    # #         Combinations <- Combinations[,sample(Ncomb, max_triplet_number)]
    # #       }
    #       Ncomb <- choose(length(tips), 3)
    #       ## If too many of them, sample from them (with replacement)
    #       if (Ncomb > max_triplet_number) {
    #         Combinations <- replicate(max_triplet_number, sample(tips, 3))
    #       } else {
    #         Combinations <- combn(x=tips, m=3)
    #       }
    #       temp_cor <- apply(X=Combinations, MARGIN=2, FUN=estimate_covariance_from_triplet, Y_data=Y_data, distances_phylo=distances_phylo)
    #       temp_dist <- apply(X=Combinations, MARGIN=2, FUN=compute_dist_triplet, distances_phylo=distances_phylo)
    #       temp_cor <- as.vector(temp_cor); temp_dist <- as.vector(temp_dist)
    #       cor_hat <- c(cor_hat, temp_cor)
    #       dists <- c(dists, temp_dist)
    #     }
    #   }
    if (length(tips) > 1){
      Z <- outer(Y_data[tips], Y_data[tips], function(x,y){x-y} )
      square_diff <- c(square_diff, (Z[upper.tri(Z)])^2)
      Z <- distances_phylo[tips,tips]
      dists <- c(dists, Z[upper.tri(Z)])
    }
  }
#     ## Empirical mean of estimators for corelations estimated several times
#     un_dists <- unique(dists)
#     un_cor_hat <- rep(NA, length(un_dists))
#   }
#   for (dd in 1:length(un_dists)){
#     un_cor_hat[dd] <- mean(cor_hat[dists==un_dists[dd]])
#   }
#   pos_cor_hat <- cor_hat[cor_hat>=0]
#   dists_pos <- dists[cor_hat>=0]
  ## Fit sigma_ij against gamma*exp(-alpha*d_ij)
#  fit <- nls(pos_cor_hat ~ gam*exp(-alpha*dists_pos), start=list(gam=mean(hat_gam), alpha=1))
#  fit <- nls(square_diff ~ gam*(1-exp(-alpha*dists)), start=list(gam=mean(hat_gam), alpha=1))
  df <- data.frame(square_diff=square_diff, dists=dists)
  fit.rob <- try(nlrob(square_diff ~ gam*(1-exp(-alpha*dists)), data=df, start=list(gam=mean(hat_gam, na.rm=TRUE), alpha=1)))
  if (inherits(fit.rob, "try-error")) {
      warning("Robust estimation of alpha failed")
      return(init.alpha.default(...))
  } else { 
      return(unname(coef(fit.rob)["alpha"]))
  }
}

compute_dist_triplet <- function(distances_phylo,v) {
  phylo_dist <- c(distances_phylo[v[1],v[2]], distances_phylo[v[2],v[3]], distances_phylo[v[1],v[3]])
  return(phylo_dist)
}

estimate_covariance_from_triplet <- function(Y_data, distances_phylo, v){
  Y_i <- Y_data[v[1]]; Y_j <- Y_data[v[2]]; Y_k <- Y_data[v[3]];
  hat_gamma <- var(Y_data[v])
  cor <- hat_gamma - (1/9) * ((Y_i + Y_j + Y_k)^2 + (Y_i - Y_j)^2 + (Y_j - Y_k)^2 + (Y_i - Y_k)^2)
  hat_sigmas <- c(Y_i*Y_j, Y_j*Y_k, Y_i*Y_k) + cor
  return(hat_sigmas)
}

init.alpha.gamma.default <- function(init.selection.strength, known.selection.strength, alpha_known, init.var.root, ...){
  return(list(alpha_0 = init.alpha.default(init.selection.strength, known.selection.strength, alpha_known),
              gamma_0 = init.var.root))
}

init.alpha.gamma.estimation <- function(phylo, Y_data, nbr_of_shifts, distances_phylo, max_triplet_number, ...){
  ## Initialize a vector with the group of each tip
  tips_groups <- rep(0, length(phylo$tip.label))
  names(tips_groups) <- phylo$tip.label
  ## Initialize shifts by a lasso without sigma
  if (nbr_of_shifts > 0) {
    lasso <- init.EM.lasso(phylo=phylo, Y_data=Y_data, process="OU", nbr_of_shifts=nbr_of_shifts, use_sigma=FALSE)
    ## Roeorder phylo and trace edges
    phy <- reorder(phylo, order = "cladewise")
    edges_shifts <- correspondanceEdges(edges=lasso$shifts$edges,from=phylo,to=phy)
    ## Set groups of tips (one group = all the tips under a given shift)
    Tr <- incidence.matrix(phy)
    for (ed in order(edges_shifts)) { # Do it in order so that edges with higher numbers erase groups (edges closer from the tips)
      ed_sh <- edges_shifts[ed]
      tips_groups[phy$tip.label[Tr[,ed_sh]]] <- ed
    }
  } else {
    edges_shifts <- NULL
  }
  ## For each group, take all the triplets of tips to estimate the covariance sigma_ij
  cor_hat <- NULL # estimations from trilpets of pairs corelations
  square_diff <- NULL # (Y_i-Y_j)^2
  dists <- NULL # corresponding phylogenetic distances between pairs
  hat_gam <- rep(0, length(edges_shifts)+1)
  for (grp in 0:length(edges_shifts)) {
    tips <- which(tips_groups==grp)
    hat_gam[grp+1] <- var(Y_data[tips])
    if (length(tips) > 1){
      Z <- outer(Y_data[tips], Y_data[tips], function(x,y){x-y} )
      square_diff <- c(square_diff, (Z[upper.tri(Z)])^2)
      Z <- distances_phylo[tips,tips]
      dists <- c(dists, Z[upper.tri(Z)])
    }
  }
  gamma_0 <- mean(hat_gam, na.rm=TRUE)
  df <- data.frame(square_diff=square_diff, dists=dists)
  fit.rob <- try(nlrob(square_diff ~ gam*(1-exp(-alpha*dists)), data=df, start=list(gam = gamma_0, alpha=1)))
  if (inherits(fit.rob, "try-error")) {
    warning("Robust estimation of alpha failed")
    return(list(alpha_0 = init.alpha.gamma.default(...)$alpha_0,
                gamma_0 = gamma_0))
  } else { 
    return(list(alpha_0 = unname(coef(fit.rob)["alpha"]), 
                gamma_0 = gamma_0))
  }
}

{ ##
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
} ##
compute_mean_variance.simple <- function (phylo, times_shared, distances_phylo, process=c("BM","OU"), params_old, ...) {
  ## Choose process 
  process  <- match.arg(process)
  compute_variance_covariance  <- switch(process, 
                                         BM = compute_variance_covariance.BM,
                                         OU = compute_variance_covariance.OU)
  ## Mean
  sim <- simulate(phylo = phylo, 
                  process = process, 
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

{ ##
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
} ##
compute_E.simple <- function (phylo, Y_data, sim, Sigma, Sigma_YY_inv) {
  ## Initialization
  ntaxa <- length(phylo$tip.label)
  conditional_law_X <- list(expectations = rep(NA,ntaxa-1), 
                            variances = rep(NA,ntaxa-1), 
                            covariances = rep(NA,ntaxa-2))
  ## Mean
  m_Y <- extract.simulate(sim, where="tips", what="expectations")
  m_Z <- extract.simulate(sim, where="nodes", what="expectations")
  conditional_law_X$optimal.values <- c(extract.simulate(sim, where="tips", what="optimal.values"), extract.simulate(sim, where="nodes", what="optimal.values")) # NULL if BM
  ## Variance Covariance
  Sigma_YZ <- extract.variance_covariance(Sigma, what="YZ")
  Sigma_ZZ <- extract.variance_covariance(Sigma, what="ZZ")
  temp <- Sigma_YZ%*%Sigma_YY_inv
  conditional_law_X$expectations <- c(Y_data, m_Z + temp%*%(Y_data-m_Y))
  conditional_variance_covariance <- Sigma_ZZ - temp%*%t(Sigma_YZ)
  conditional_law_X$variances <- c(rep(0, ntaxa), diag(conditional_variance_covariance))
  conditional_law_X$covariances <- c(rep(0, ntaxa), extract.covariance_parents(phylo, conditional_variance_covariance))
  return(conditional_law_X)
}

{ ##
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
} ##
compute_log_likelihood.simple <- function(phylo, Y_data, sim, Sigma, Sigma_YY_inv){
  ntaxa <- length(phylo$tip.label)
  Sigma_YY <- extract.variance_covariance(Sigma, what="YY")
  logdetSigma_YY <- determinant(Sigma_YY, logarithm = TRUE)$modulus
  m_Y <- extract.simulate(sim, where="tips", what="expectations")
  LL <- ntaxa * log(2*pi) + logdetSigma_YY
  LL <- LL + t(Y_data - m_Y)%*%Sigma_YY_inv%*%(Y_data-m_Y)
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

{ ##
# # compute_E.simple (phylo, times_shared, distances_phylo, Y_data, process=c("BM","OU"), params_old, ...)
# # PARAMETERS:
# #            @phylo (tree) imput tree
# #            @Y_data (vector) : vector indicating the data at the tips
# #            @process (string) Random process to simulate.
# #            @params_old (list) : old parameters to be used in the E step
# #            @times_shared (matrix) : times of shared ancestry, result of function compute_times_ca (see note below)
# #            @distances_phylo (matrix) : phylogenetics distance, result of function compute_dist_phy (see note below)
# # RETURNS:
# #            @conditional_law_X (list) : list of conditionnal statistics :
# #                   "expectation" : vector of length ntaxa+nNodes, with ntaxa fisrt values set to Y_data (tips), and from ntaxa+1 of conditional expectation of the nodes conditionned to the tips E[Z_j|Y]
# #                   "variances" : vector of length ntaxa+nNodes with ntaxa 0 (tips) and conditional variance of the nodes conditionned to the tips Var[Z_j|Y]
# #                   "covariances" : vector of length ntaxa+nNodes with ntaxa 0 (tips) and conditional covariance of the nodes and their parents conditionned to the tips Cov[Z_j,Z_pa(j)|Y], with NA for the root.
# #                   "optimal.values" : vector of length ntaxa+nNodes of optimal values beta(t_j)
# # DEPENDENCIES:
# #            compute_times_ca, compute_dist_phy, compute_variance_covariance, extract.variance_covariance, extract.covariance_parents, simulate, extract.simulate
# # PURPOSE:
# #            E step in the simple case
# # NOTES:
# #            none
# # REVISIONS:
# #            21/05/14 - Initial release (non fonctionnal)
# #            22/05/14 - Minimal "working" release
# #            02/06/14 - Add case OU
# } ##
# compute_E.simple <- function (phylo, times_shared, distances_phylo, Y_data, process=c("BM","OU"), params_old, ...) {
#   ## Choose process 
#   process  <- match.arg(process)
#   compute_variance_covariance  <- switch(process, 
#                                          BM = compute_variance_covariance.BM,
#                                          OU = compute_variance_covariance.OU)
#   ## Initialization
#   ntaxa <- length(phylo$tip.label)
#   conditional_law_X <- list(expectations=rep(NA,ntaxa-1), variances=rep(NA,ntaxa-1), covariances=rep(NA,ntaxa-2))
#   ## Mean
#   sim <- simulate(phylo=phylo, 
#                   process=process, 
#                   root.state = params_old$root.state, 
#                   shifts = params_old$shifts, 
#                   variance=params_old$variance, 
#                   optimal.value=params_old$optimal.value, 
#                   selection.strength=params_old$selection.strength)
#   m_Y <- extract.simulate(sim, where="tips", what="expectations")
#   m_Z <- extract.simulate(sim, where="nodes", what="expectations")
#   conditional_law_X$optimal.values <- c(extract.simulate(sim, where="tips", what="optimal.values"), extract.simulate(sim, where="nodes", what="optimal.values")) # NULL if BM
#   ## Variance Covariance
#   Sigma <- compute_variance_covariance(times_shared = times_shared, 
#                                        distances_phylo = distances_phylo,
#                                        params_old = params_old)
#   Sigma_YY <- extract.variance_covariance(Sigma, what="YY")
#   Sigma_YZ <- extract.variance_covariance(Sigma, what="YZ")
#   Sigma_ZZ <- extract.variance_covariance(Sigma, what="ZZ")
#   ## Conditionnal mean and variance covariance
#   Sigma_YY_inv <- solve(Sigma_YY)
#   temp <- Sigma_YZ%*%Sigma_YY_inv
#   conditional_law_X$expectations <- c(Y_data, m_Z + temp%*%(Y_data-m_Y))
#   conditional_variance_covariance <- Sigma_ZZ - temp%*%t(Sigma_YZ)
#   conditional_law_X$variances <- c(rep(0, ntaxa), diag(conditional_variance_covariance))
#   conditional_law_X$covariances <- c(rep(0, ntaxa), extract.covariance_parents(phylo, conditional_variance_covariance))
#   return(conditional_law_X)
}

{ ##
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
} ##
compute_times_ca <- function(phy) {
  times <- node.depth.edgelength(phy)
  prac <- mrca(phy,full=TRUE)
  times_ca <- matrix(times[prac],dim(prac))
  attr(times_ca, "ntaxa")  <- length(phy$tip.label)
  return(times_ca)
}

{ ##
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
} ##
compute_dist_phy <- function(phy) {
  dist_phy <- dist.nodes(phy)
  attr(dist_phy, "ntaxa")  <- length(phy$tip.label)
  return(dist_phy)
}

{ ##
# extract.variance_covariance (struct, what=c("YY","YZ","ZZ"))
# PARAMETERS:
#            @struct (matrix) structural matrix of size ntaxa+nNode, result of function compute_times_ca, compute_dist_phy or compute_variance_covariance
#            @what (string) what to extract :
#                 "YY" : sub-matrix of tips (ntaxa first lines and columns)
#                 "YZ" : sub matrix tips x nodes (nNodes last rows and ntaxa first columns)
#                 "ZZ" : sub matrix of nodes (nNodes last rows and columns)
# RETURNS:
#            (matrix) : sub-matrix of the entry matrix corresponding to the wanted values
# DEPENDENCIES:
#            none
# PURPOSE:
#            Extract the right sub matrix
# NOTES:
#            none
# REVISIONS:
#            22/05/14 - Initial release
} ##
extract.variance_covariance <- function(struct, what=c("YY","YZ","ZZ")){
  ntaxa <- attr(struct, "ntaxa")
  if (what=="YY") {
    return(struct[1:ntaxa,1:ntaxa])
  } else if (what=="YZ") {
    return(struct[(ntaxa+1):(dim(struct)[1]),1:ntaxa])
  } else if (what=="ZZ") {
    return(struct[(ntaxa+1):(dim(struct)[1]),(ntaxa+1):(dim(struct)[2])])
  }
}

{ ##
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
} ##
extract.covariance_parents<- function(phylo, struct){
  ntaxa <- length(phylo$tip.label)
  cov <- rep(NA,dim(phylo$edge)[1] - ntaxa + 1)
  for (i in (ntaxa+2):(dim(phylo$edge)[1]+1)) {
    pa <- getAncestor(phylo,i)
    cov[i-ntaxa] <- struct[i-ntaxa,pa-ntaxa]
  }
  return(cov)
}

{ ##
# compute_variance_covariance.BM (times_shared, params_old, ...) 
# PARAMETERS:
#            @times_shared (matrix) : times of shared ancestry, result of function compute_times_ca (see note above)
#            @params_old (list) : old parameters to be used in the E step
# RETURNS:
#            (matrix) : matrix of variance covariance for the BM
# DEPENDENCIES:
#            compute_times_ca
# PURPOSE:
#            Compute variance covariance matrix in the case of the BM
# NOTES:
#            none
# REVISIONS:
#            22/05/14 - Initial release
} ##
compute_variance_covariance.BM <- function(times_shared, params_old, ...) {
  J <- matrix(1, nrow=dim(times_shared)[1], ncol=dim(times_shared)[2])
  if (params_old$root.state$random) {
    return(params_old$root.state$var.root * J + params_old$variance * times_shared)
  } else {
    return(params_old$variance * times_shared)
  }
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

{ ##
# compute_M (phylo, process, Y_data, conditional_law_X, nbr_of_shifts)
# PARAMETERS:
#            @phylo (tree) imput tree
#            @process (string) Random process to simulate.
#            @Y_data (vector) : vector indicating the data at the tips
#            @conditional_law_X (list) result of compute_E (see note above)
#            @nbr_of_shifts : number of shifts wanted for the inference
# RETURNS:
#            (list) list of parameters for the model fitted
# DEPENDENCIES:
#            compute_E, init.EM.default, compute_diff_exp, compute_var_diff
# PURPOSE:
#            M step
# NOTES:
#            none
# REVISIONS:
#            22/05/14 - Initial release
} ##

# compute_M <- function(process) {
#   if (process=="BM"){
#     return(compute_M.BM)
#   }
# }

compute_M.BM <- function(phylo, Y_data, conditional_law_X, nbr_of_shifts, random.root, ...) {
  ## Initialization
  ntaxa <- length(phylo$tip.label)
  params <- init.EM.default.BM(random.init = random.root, ...)
  ## Actualization of the root
  if (params$root.state$random) {
    params$root.state$value.root <- NA
    params$root.state$exp.root <- conditional_law_X$expectations[ntaxa+1]
    params$root.state$var.root <- conditional_law_X$variances[ntaxa+1]
  } else {
    params$root.state$value.root <- conditional_law_X$expectations[ntaxa+1]
    params$root.state$exp.root <- NA
    params$root.state$var.root <- NA
  }
  ## Segmentation
  diff_exp <- compute_diff_exp.BM(phylo=phylo, 
                                  conditional_law_X=conditional_law_X)
  costs0 <- 1/(phylo$edge.length) * diff_exp^2
  seg <- segmentation.BM(nbr_of_shifts=nbr_of_shifts, 
                         costs0=costs0,
                         diff_exp=diff_exp)
  params$shifts <- seg$shifts
  edges_max <- seg$edges_max
  ## Variance
  var_diff <- compute_var_diff.BM(phylo=phylo, 
                                  conditional_law_X=conditional_law_X)
  params$variance <- compute_var_M.BM(phylo=phylo,
                                      var_diff=var_diff,
                                      costs0=costs0,
                                      edges_max=edges_max)
  return(params)
}

compute_M.OU <- function(stationnary.root, shifts_at_nodes, alpha_known){
  if (stationnary.root && shifts_at_nodes && alpha_known) {
    return(compute_M.OU.specialCase)
  } else if (stationnary.root && shifts_at_nodes) {
    return(compute_M.OU.stationnary.root_AND_shifts_at_nodes)
  } else {
    stop("The EM algorithm for the OU is only defined (for the moment) for a stationnary root and shifts at nodes !")
  }
}

compute_M.OU.specialCase <- function(phylo, Y_data, conditional_law_X, nbr_of_shifts, known.selection.strength, methods.segmentation, beta_0_old, shifts_old, ...){
  ## Initialization
  ntaxa <- length(phylo$tip.label)
  params <- init.EM.default.OU(selection.strength.init=known.selection.strength,
                               stationnary.root.init=TRUE,
                               random.init=TRUE,
                               ...)
  ## Choose method(s) for segmentation
  segmentation.OU.specialCase <- function(method.segmentation){
    segmentation <- switch(method.segmentation, 
                           max_costs_0 = segmentation.OU.specialCase.max_costs_0,
                           lasso = segmentation.OU.specialCase.lasso,
                           same_shifts = segmentation.OU.specialCase.same_shifts,
                           same_shifts_same_values = segmentation.OU.specialCase.same_shifts_same_values,
                           best_single_move = segmentation.OU.specialCase.best_single_move)
    return(segmentation(phylo = phylo, 
                        nbr_of_shifts = nbr_of_shifts, 
                        conditional_law_X = conditional_law_X, 
                        selection.strength = known.selection.strength,
                        beta_0_old = beta_0_old,
                        shifts_old = shifts_old))
  }
  ## Segmentation
  segs <- sapply(methods.segmentation, segmentation.OU.specialCase, simplify = FALSE)
 #edges_max <- seg$edges_max
  ## Variance
  var_diff <- compute_var_diff.OU(phylo=phylo, 
                                  conditional_law_X=conditional_law_X, 
                                  selection.strength=known.selection.strength)
  # Function to compute gamma2 for all parameters obtained by all segmentations
  compute_var_M <- function(method.segmentation){
    return(compute_var_M.OU.specialCase(phylo=phylo, 
                                        var_diff=var_diff, 
                                        costs=segs[[method.segmentation]]$costs, 
                                        selection.strength=known.selection.strength,
                                        conditional_root_variance=unname(conditional_law_X$variances[ntaxa+1])))
  }
  var.roots <- sapply(methods.segmentation, compute_var_M)
  ## Compute objective function for each set of parameters, and choose the best one
  cond_exp_log_lik <- function(method.segmentation){
    return(unname(conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes(phylo = phylo,
                conditional_law_X = conditional_law_X, 
                sigma2 = 2 * known.selection.strength * var.roots[method.segmentation],
                mu = segs[[method.segmentation]]$beta_0,
                shifts = segs[[method.segmentation]]$shifts,
                alpha = known.selection.strength)))
  }
  obj_funcs <- sapply(methods.segmentation, cond_exp_log_lik)
  best.method.seg <- which.max(obj_funcs)
  ## Actualize paremters with the ones found by the best segmentation method
  params$root.state$var.root <- unname(var.roots[best.method.seg])
  params$variance <- 2 * known.selection.strength * params$root.state$var.root
  params$shifts <- segs[[best.method.seg]]$shifts
  params$root.state$exp.root <- segs[[best.method.seg]]$beta_0
  params$optimal.value <- params$root.state$exp.root
  attr(params, "segmentation_algorithm_used") <- names(best.method.seg)
  return(params)
}

compute_M.OU.stationnary.root_AND_shifts_at_nodes <- function(phylo, Y_data, conditional_law_X, nbr_of_shifts, alpha_old, max_selection.strength, eps, methods.segmentation, beta_0_old = beta_0_old, shifts_old = shifts_old, ...){
  ## Estimate all parameters with alpha of the previous step
  params <- compute_M.OU.specialCase(phylo = phylo, 
                                     Y_data = Y_data, 
                                     conditional_law_X = conditional_law_X,
                                     nbr_of_shifts = nbr_of_shifts,
                                     known.selection.strength = alpha_old,
                                     methods.segmentation = methods.segmentation,
                                     beta_0_old = beta_0_old,
                                     shifts_old = shifts_old)
  ## Estimate new alpha
  params$selection.strength <- estimate.alpha(phylo = phylo,
                                              conditional_law_X = conditional_law_X, 
                                              sigma2 = params$variance,
                                              mu = params$root.state$exp.root,
                                              shifts = params$shifts,
                                              alpha_old = alpha_old,
                                              max_selection.strength = max_selection.strength)
  ## Change value of the root (stationnary) accordingly
#  params$root.state <- suppressWarnings(test.root.state(root.state=params$root.state, process="OU", variance=params$variance, selection.strength=params$selection.strength, optimal.value=params$optimal.value))
# ATTENTION : Changement dans la manière de mettre à jour
   params$variance <- 2 * params$selection.strength * params$root.state$var.root
  return(params)
}

{ ##
#' @title Function to estimate alpha
#'
#' @description
#' \code{optimize} to maximize function
#' \code{conditional_expectation_log_likelihood.OU} in alpha. The interval is
#' set to [alpha_old/2, 2*alpha_old], supposing that the previous guess of 
#' alpha_old is not far from reality.
#'
#' @details
#' This function uses functions \code{compute_var_diff.OU} 
#' and \code{compute_diff_exp.OU} in the process. Carefull : only works if the
#' root is stationnary, and shifts at nodes.
#'
#' @param phylo Input tree.
#' @param conditional_law_X result of function \code{compute_E.OU}
#' @param sigma2 variance of params
#' @param mu mean of the root state
#' @param shifts list of shifts on the tree
#' @param alpha_old previous estimation of the selection strength
#' @param max_selection.strength the maximal value of alpha authorized by
#'  the user
#' 
#' @return double : estimation of alpha
#'
#'09/07/14 - Initial release
#'02/10/14 - Take newly estimated shifts into consideration
} ##
estimate.alpha <- function(phylo, conditional_law_X, sigma2, mu, shifts, alpha_old, max_selection.strength){
  #opt <- optimize(R_function, phylo=phylo, conditional_law_X=conditional_law_X, sigma2=sigma2, mu=mu, interval=c(alpha_old/2, alpha_old*2))
  betas <- compute_betas(phylo = phylo, optimal.value = mu, shifts = shifts)
  opt2 <- optimize(conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes.from_betas, 
                   phylo = phylo, conditional_law_X = conditional_law_X, 
                   sigma2 = sigma2, 
                   mu = mu, 
                   betas = betas,
                   interval = c(max(0.01, alpha_old/2), 
                                min(max_selection.strength, alpha_old*2)), 
                   maximum = TRUE)
  # Check that function is concave on alpha found
  #hess <- fdHess(opt2$maximum, fun=conditional_expectation_log_likelihood.OU, phylo=phylo, conditional_law_X=conditional_law_X, sigma2=sigma2, mu=mu)
  #if (hess$Hessian > 0) warning("The function to maximize is not concave in the maximum fund.")
  return(opt2$maximum)
}

{ ##
#' @title Function to be minimized in alpha
#'
#' @description
#' \code{R_function} is the oposite of the conditional expectation of the 
#' completed log-likelihood.
#'
#' @details
#' This function uses functions \code{compute_var_diff.OU} 
#' and \code{compute_diff_exp.OU} in the process. Carefull : only works if the
#' root is stationnary, and shifts at nodes.
#'
#' @param phylo Input tree.
#' @param conditional_law_X result of function \code{compute_E.OU}
#' @param sigma2 variance of params
#' @param alpha selection strength to be estimated
#' 
#' @return double : value of the function
#'
#'09/07/14 - Initial release
} ##
R_function <- function(phylo, conditional_law_X, sigma2, mu, alpha){
  ntaxa <- length(phylo$tip.label)
  nNodes <- phylo$Nnode
  ## Terms in log
  ee <- exp(- alpha * phylo$edge.length )
  K_1 <- conditional_law_X$variances[ntaxa+1] + (conditional_law_X$expectations[ntaxa+1] - mu)^2
  R_func <- K_1 * alpha / sigma2 - (nNodes + ntaxa - 1)*log(alpha)/2 + sum(log(1-ee^2))/2
  ## Terms with the variance
  var_diff <- compute_var_diff.OU(phylo=phylo, 
                                  conditional_law_X=conditional_law_X, 
                                  selection.strength=alpha)
  R_func <- R_func + alpha*sum((1 - ee^2)^(-1) * var_diff)/sigma2
  ## Terms with the expectation
  diff_exp <- compute_diff_exp.OU(phylo=phylo, 
                                  conditional_law_X=conditional_law_X, 
                                  selection.strength=alpha)
  daughters <- phylo$edge[,2]
  betas <- conditional_law_X$optimal.values[daughters]
  R_func <- R_func + alpha*sum((1 - ee^2)^(-1) * (diff_exp - betas * (1-ee))^2)/sigma2
  return(R_func)
}


conditional_expectation_log_likelihood.OU <- function(stationnary.root, shifts_at_nodes, alpha_known){
  if (stationnary.root && shifts_at_nodes) {
    return(conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes)
  } else {
    stop("The EM algorithm for the OU is only defined (for the moment) for a stationnary root and shifts at nodes !")
  }
}

{ ##
#' @title Expectation conditional to the tips of the completed log-likelihood.
#' This is an old version.
#'
#' @description
#' \code{conditional_expectation_log_likelihood.OU} computes the expectation
#'  conditional to the tips of the completed log-likelihood, given the 
#'  first and second order moments of the nodes, and the parameters of the 
#'  OU process.
#'
#' @details
#' This function uses functions \code{compute_var_diff.OU} 
#' and \code{compute_diff_exp.OU} in the process. Carefull : only works if the
#' root is stationnary, and shifts at nodes.
#'
#' @param phylo Input tree.
#' @param conditional_law_X result of function \code{compute_E.OU}, containing
#' first and second moments.
#' @param sigma2 variance of the OU.
#' @param mu mean of the root.
#' @param alpha selection strength of the OU.
#' 
#' @return double : value of the function
#'
#'09/07/14 - Initial release
} ##
conditional_expectation_log_likelihood.OU.OLD <- function(phylo, conditional_law_X, sigma2, mu, alpha){
  ntaxa <- length(phylo$tip.label)
  nNodes <- phylo$Nnode
  ## Constante
  cst <- -1/2 * ((nNodes + ntaxa) * log(2*pi))
  LogLik <- cst
  ## Terms in log
  ee <- exp(- alpha * phylo$edge.length )
  LogLik <- LogLik - (nNodes + ntaxa - 1)*log(sigma2/(2*alpha))/2 - sum(log(1-ee^2))/2
  ## Terms with the variance
  K_1 <- conditional_law_X$variances[ntaxa+1] + (conditional_law_X$expectations[ntaxa+1] - mu)^2
  var_diff <- compute_var_diff.OU(phylo=phylo, 
                                  conditional_law_X=conditional_law_X, 
                                  selection.strength=alpha)
  LogLik <- LogLik - (alpha / sigma2) * (K_1 + sum((1 - ee^2)^(-1) * var_diff))
  ## Terms with the expectation
  diff_exp <- compute_diff_exp.OU(phylo=phylo, 
                                  conditional_law_X=conditional_law_X, 
                                  selection.strength=alpha)
  daughters <- phylo$edge[,2]
  betas <- conditional_law_X$optimal.values[daughters]
  LogLik <- LogLik - (alpha / sigma2) * sum((1 - ee^2)^(-1) * (diff_exp - betas * (1-ee))^2)
  return(LogLik)
}

{ ##
#' @title Expectation conditional to the tips of the completed log-likelihood.
#'
#' @description
#' \code{conditional_expectation_log_likelihood.OU} computes the expectation
#'  conditional to the tips of the completed log-likelihood, given the 
#'  first and second order moments of the nodes, and the parameters of the 
#'  OU process.
#'
#' @details
#' This function uses functions \code{compute_var_diff.OU} 
#' and \code{compute_diff_exp.OU} in the process. Carefull : only works if the
#' root is stationnary, and shifts at nodes.
#'
#' @param phylo Input tree.
#' @param conditional_law_X result of function \code{compute_E.OU}, containing
#' first and second moments.
#' @param sigma2 variance of the OU.
#' @param mu mean of the root.
#' @param alpha selection strength of the OU.
#' 
#' @return double : value of the function
#'
#'09/07/14 - Initial release
#'02/10/14 - take new shifts in consideration
} ##
conditional_expectation_log_likelihood.OU.stationary_root_shifts_at_nodes <- function(phylo, conditional_law_X, sigma2, mu, shifts, alpha){
  ntaxa <- length(phylo$tip.label)
  nNodes <- phylo$Nnode
  ## Constante
  cst <- -1/2 * ((nNodes + ntaxa) * log(2*pi))
  LogLik <- cst
  ## Terms in log
  ee <- exp(- alpha * phylo$edge.length )
  LogLik <- LogLik - (nNodes + ntaxa - 1)*log(sigma2/(2*alpha))/2 - sum(log(1-ee^2))/2
  ## Terms with the variance
  K_1 <- conditional_law_X$variances[ntaxa+1] + (conditional_law_X$expectations[ntaxa+1] - mu)^2
  var_diff <- compute_var_diff.OU(phylo=phylo, 
                                  conditional_law_X=conditional_law_X, 
                                  selection.strength=alpha)
  LogLik <- LogLik - (alpha / sigma2) * (K_1 + sum((1 - ee^2)^(-1) * var_diff))
  ## Terms with the expectation
  diff_exp <- compute_diff_exp.OU(phylo=phylo, 
                                  conditional_law_X=conditional_law_X, 
                                  selection.strength=alpha)
  parents <- phylo$edge[,1]
  daughters <- phylo$edge[,2]
  betas <- conditional_law_X$optimal.values[parents]
  delta <- shifts.list_to_vector(phylo, shifts) # vector of shifts (branches)
  LogLik <- LogLik - (alpha / sigma2) * sum((1 - ee^2)^(-1) * (diff_exp - (betas + delta) * (1-ee))^2)
  return(LogLik)
}

conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes <- function(phylo, conditional_law_X, sigma2, mu, shifts, alpha){
  betas <- compute_betas(phylo = phylo, optimal.value = mu, shifts = shifts)
  return(conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes.from_betas(phylo, conditional_law_X, sigma2, mu, betas, alpha))
}

conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes.from_betas <- function(phylo, conditional_law_X, sigma2, mu, betas, alpha){
  ntaxa <- length(phylo$tip.label)
  nNodes <- phylo$Nnode
  ## Constante
  cst <- -1/2 * ((nNodes + ntaxa) * log(2*pi))
  LogLik <- cst
  ## Terms in log
  ee <- exp(- alpha * phylo$edge.length )
  LogLik <- LogLik - (nNodes + ntaxa)*log(sigma2/(2*alpha))/2 - sum(log(1-ee^2))/2
  ## Terms with the variance
  K_1 <- conditional_law_X$variances[ntaxa+1] + (conditional_law_X$expectations[ntaxa+1] - mu)^2
  var_diff <- compute_var_diff.OU(phylo=phylo, 
                                  conditional_law_X=conditional_law_X, 
                                  selection.strength=alpha)
  LogLik <- LogLik - (alpha / sigma2) * (K_1 + sum((1 - ee^2)^(-1) * var_diff))
  ## Terms with the expectation
  diff_exp <- compute_diff_exp.OU(phylo = phylo, 
                                  conditional_law_X = conditional_law_X, 
                                  selection.strength = alpha)
  parents <- phylo$edge[,1]
  daughters <- phylo$edge[,2]
  betas <- betas[daughters]
  LogLik <- LogLik - (alpha / sigma2) * sum((1 - ee^2)^(-1) * (diff_exp - betas * (1-ee))^2)
  return(LogLik)
}


{ ##
#' @title Segmentation in the BM case
#'
#' @description
#' \code{segmentation.BM} performs the segmentation algorithm described.
#'
#' @details
#'  This function takes the largest values of costs0, and make them null,
#'  thanks to delta and tau.
#'
#' @param nbr_of_shifts Number of shifts on the phylogeny allowed
#' @param costs0 Cost of each edge
#' @param diff_exp Difference of expectations
#' 
#' @return List containing : edges_max : array of nbr_of_shifts edges where
#' costs0 is maximal.
#'                           shifts:list containing the computed tau and delta
#'
#'02/06/14 - Initial release
} ##
segmentation.BM <- function(nbr_of_shifts, costs0, diff_exp){
  if (nbr_of_shifts>0) {
    edges_max <- order(-costs0)[1:nbr_of_shifts]
    edges <- edges_max
    values <- diff_exp[edges_max]
    relativeTimes <- rep(0,length(edges_max))
    shifts=list(edges=edges, values=values, relativeTimes=relativeTimes)
    return(list(edges_max=edges_max, shifts=shifts))
  } else {
    edges_max <- length(costs0)+1
    return(list(edges_max=edges_max, shifts=NULL))
  }
}

{ ##
#' @title Segmentation in the OU special case, algo 1
#'
#' @description
#' \code{segmentation.OU.specialCase.max_costs_0} performs the first segmentation 
#' algorithm described.
#'
#' @details
#'  This function takes the largest values of costs0, and make them null,
#'  thanks to delta and tau.
#'
#' @param phylo a phylogenetic tree
#' @param nbr_of_shifts Number of shifts on the phylogeny allowed
#' @param conditional_law_X moments of the conditional law of X given Y, result
#' of function \code{compute_M.OU.specialCase}
#' @param selection.strength the selection strength
#' 
#' @return List containing : beta_0 : the optimal value at the root
#'                           shifts : list containing the computed tau and delta
#'                           costs : vector of costs
#'
#'10/06/14 - Initial release
#'06/10/14 - Change name to include other algorithms
} ##
segmentation.OU.specialCase.max_costs_0 <- function(phylo, nbr_of_shifts, conditional_law_X, selection.strength, ...){
  ntaxa <- length(phylo$tip.label)
  ## Computation of mu=beta0
  beta_0 <- conditional_law_X$expectations[ntaxa+1]
  ## Computation of costs
  diff_exp <- compute_diff_exp.OU(phylo = phylo, 
                                  conditional_law_X = conditional_law_X, 
                                  selection.strength = selection.strength)
  parents <- phylo$edge[,1]
  betas <- conditional_law_X$optimal.values
  ee <- exp(- selection.strength * phylo$edge.length )
  costs0 <- (1 - ee^2)^(-1) * (diff_exp - betas[parents] * (1-ee))^2
  ## Segmentation per se
  if (nbr_of_shifts > 0) {
    # Max of costs
    edges_max <- order(-costs0)[1:nbr_of_shifts]
    edges <- edges_max
    # Put them to 0
    ee <- exp(- selection.strength * phylo$edge.length )
    parents <- phylo$edge[edges_max,1]
    values <- (1 - ee[edges_max])^(-1) * diff_exp[edges_max] - betas[parents]
    relativeTimes <- rep(0,length(edges_max))
    # Define shifts
    shifts <- list(edges = edges, values = values, relativeTimes = relativeTimes)
    # Compute new costs
    costs <- costs0
    costs[edges] <- 0
    return(list(beta_0 = beta_0, edges_max = edges_max, shifts = shifts, costs = costs))
  } else {
    edges_max <- length(costs0) + 1
    return(list(beta_0 = beta_0, edges_max = edges_max, shifts = NULL, costs = costs0))
  }
}

{ ##
#' @title Segmentation in the OU special case, using lasso regression
#'
#' @description
#' \code{segmentation.OU.specialCase.lasso} performs the segmentation using a 
#' lasso regression to select for the edges where the shifts are added.
#'
#' @details
#'  This function re-write the sum of costs to be minimized as a least squares 
#'  regression problem, and uses a lasso regression to solve it. It uses
#'  functions \code{incidence.matrix.full} to express the problem as a 
#'  linear model.
#'  
#' @param phylo a phylogenetic tree
#' @param nbr_of_shifts Number of shifts on the phylogeny allowed
#' @param conditional_law_X moments of the conditional law of X given Y, result
#' of function \code{compute_M.OU.specialCase}
#' @param selection.strength the selection strength
#' 
#' @return List containing : beta_0 : the optimal value at the root
#'                           shifts : list containing the computed tau and delta
#'                           costs : vector of costs
#'                           
#'06/10/14 - Initial release
} ##
segmentation.OU.specialCase.lasso <- function(phylo, nbr_of_shifts, conditional_law_X, selection.strength, ...){
  ntaxa <- length(phylo$tip.label)
  ## Computation of answer matrix D
  diff_exp <- compute_diff_exp.OU(phylo = phylo, 
                                  conditional_law_X = conditional_law_X, 
                                  selection.strength = selection.strength)
  parents <- phylo$edge[,1]
  ee <- exp(- selection.strength * phylo$edge.length )
  D <- (1 - ee^2)^(-1/2) * diff_exp
  D <- c(conditional_law_X$expectations[ntaxa+1], D)
  ## Regression matrix : modified incidence matrix
  U <- incidence.matrix.full(phylo)
  U <- cbind(rep(1, dim(U)[1]), U)
  A <- diag(c(1,sqrt((1-ee)/(1+ee))))
  Xp <- A%*%U
  ## Segmentation per se
  if (nbr_of_shifts > 0) {
    # Lasso regression
    fit <- lasso_regression_K_fixed(Yp = D, Xp = Xp, K = nbr_of_shifts, intercept.penalty = TRUE)
    # Define shifts
    shifts <- fit$shifts.gauss
    shifts$edges <- shifts$edges - 1
    # Define mu = beta_0
    beta_0 <- fit$E0.gauss
    # Compute new costs
    Delta <- shifts.list_to_vector(phylo, shifts)
    Delta <- c(beta_0, Delta)
    costs <- (D - Xp%*%Delta)^2
    return(list(beta_0 = beta_0, shifts = shifts, costs = costs))
  } else {
   #edges_max <- length(costs0) + 1
    beta_0 <- conditional_law_X$expectations[ntaxa+1]
    Delta <- c(beta_0, rep(0, length(phylo$edge)))
    costs <- (D - Xp%*%Delta)^2
    return(list(beta_0 = beta_0, shifts = NULL, costs = costs))
  }
}

{ ##
#' @title Segmentation in the OU special case, conserving the same shifts.
#'
#' @description
#' \code{segmentation.OU.specialCase.same_shifts_same_values} keeps the same
#' parameters and compute the quantities needed. It is here to ensure that we do
#' at least better than the previous step.
#'
#' @details
#'  This function takes the old shifts parameters, and compute costs with the new
#'  moments.
#'
#' @param phylo a phylogenetic tree
#' @param conditional_law_X moments of the conditional law of X given Y, result
#' of function \code{compute_M.OU.specialCase}
#' @param selection.strength the selection strength
#' @param beta_0_old the previous (and concerved) value of beta_0
#' @param shifts_old the previous (and concerved) list of shifts
#' 
#' @return List containing : beta_0 : the optimal value at the root
#'                           shifts : list containing the computed tau and delta
#'                           costs : vector of costs
#'
#'06/10/14 - Initial release
} ##
segmentation.OU.specialCase.same_shifts_same_values <- function(phylo, conditional_law_X, selection.strength, beta_0_old, shifts_old, ...){
  ntaxa <- length(phylo$tip.label)
  ## Computation of costs
  diff_exp <- compute_diff_exp.OU(phylo = phylo, 
                                  conditional_law_X = conditional_law_X, 
                                  selection.strength = selection.strength)
  betas <- compute_betas(phylo = phylo, optimal.value = beta_0_old, shifts = shifts_old)
  daughters <- phylo$edge[,2]
  ee <- exp(- selection.strength * phylo$edge.length )
  costs <- (1 - ee^2)^(-1) * (diff_exp - betas[daughters] * (1-ee))^2
  costs <- c((conditional_law_X$expectations[ntaxa+1] - beta_0_old)^2, costs)
  return(list(beta_0 = beta_0_old, shifts = shifts_old, costs = costs))
}

{ ##
#' @title Segmentation in the OU special case, conserving the same shifts
#'  position.
#'
#' @description
#' \code{segmentation.OU.specialCase.same_shifts} keeps the same shifts position,
#' and optimize the sum of costs using function 
#' \code{optimize_costs_given_shift_position.OU.specialCase}. 
#'
#' @details
#'  This is the best move if keeping the previous shifts positions.
#'
#' @param phylo a phylogenetic tree
#' @param conditional_law_X moments of the conditional law of X given Y, result
#' of function \code{compute_M.OU.specialCase}
#' @param selection.strength the selection strength
#' @param shifts_old the previous list of shifts (only position is used)
#' 
#' @return List containing : beta_0 : the optimal value at the root
#'                           shifts : list containing the computed tau and delta
#'                           costs : vector of costs
#'
#'06/10/14 - Initial release
} ##
segmentation.OU.specialCase.same_shifts <- function(phylo, conditional_law_X, selection.strength, shifts_old, ...){
  return(optimize_costs_given_shift_position.OU.specialCase(phylo = phylo,
                                        conditional_law_X = conditional_law_X,
                                        selection.strength = selection.strength,
                                        shifts_edges = shifts_old$edges))
}

segmentation.OU.specialCase.best_single_move <- function(phylo, conditional_law_X, selection.strength, shifts_old, ...){
  # Construct vector of all allowed combinations of shifts (variation from a base scenario)
  shifts_edges <- shifts_old$edges
  K <- length(shifts_edges)
  nEdges <- dim(phylo$edge)[1]
  allowed_moves <- which(!(1:nEdges %in% shifts_edges))
  scenarii <- t(matrix(shifts_edges, (nEdges - K) * K + 1, nrow = K))
  for (i in 1:K) {
    scenarii[1 + ((i-1)*(nEdges - K)+1):(i*(nEdges - K)), i] <- allowed_moves
  }
  # Function to be applyed to each row
  fun <- function(sh_ed){
    seg <- optimize_costs_given_shift_position.OU.specialCase(phylo = phylo,
                                        conditional_law_X = conditional_law_X,
                                        selection.strength = selection.strength,
                                        shifts_edges = sh_ed)
    totalCost <- sum(seg$costs)
    return(list(seg = seg, 
                totalCost = totalCost))
  }
  # Apply to each row and take the minimal total cost
  allSegs <- apply(scenarii, 1, fun)
  dd <- do.call(rbind, allSegs)
  min_conf <- which.min(dd[,"totalCost"])
  return(dd[[min_conf]])
}

{ ##
#' @title Minimization of the sum of costs, given the shift position.
#'
#' @description
#' \code{optimize_costs_given_shift_position.OU.specialCase} minimize the sum of
#' costs when the shift position is fixed.
#'
#' @details
#' This function find the regimes of each node using function 
#' \code{allocate_regimes_from_shifts} and optimize the sum of costs, computed 
#' using function \code{compute_diff_exp.OU}, in the values of the optimal values
#' betas (using a close formula). It then goes back to a shift expression of the
#' problem using function \code{compute_shifts_from_betas}.
#'
#' @param phylo a phylogenetic tree
#' @param conditional_law_X moments of the conditional law of X given Y, result
#' of function \code{compute_M.OU.specialCase}
#' @param selection.strength the selection strength
#' @param shifts_edges the vector of the position of the shifts on the tree
#' 
#' @return List containing : beta_0 : the optimal value at the root
#'                           shifts : list containing the computed tau and delta
#'                           costs : vector of costs
#'
#'15/10/14 - Initial release
} ##
optimize_costs_given_shift_position.OU.specialCase <- function(phylo, conditional_law_X, selection.strength, shifts_edges, ...){
  ntaxa <- length(phylo$tip.label)
  ## Computation values of betas
  diff_exp <- compute_diff_exp.OU(phylo = phylo, 
                                  conditional_law_X = conditional_law_X, 
                                  selection.strength = selection.strength)
  # Regimes of the branches from alod positions of shifts
  regimes <- allocate_regimes_from_shifts(phylo, shifts_edges)
  # Quantities needed
  daughters <- phylo$edge[,2]
  ee <- exp(- selection.strength * phylo$edge.length )
  numerateur <- diff_exp / (1 + ee)
  denominateur <- (1 - ee) / (1 + ee)
  # Root branch
  nEdges <- dim(phylo$edge)[1]
  numerateur[nEdges + 1]  <- conditional_law_X$expectations[ntaxa+1]
  denominateur[nEdges + 1]  <- 1
  # Function to compute betas on the regimes
  fun <- function(reg){
    edg <- which(phylo$edge[,2] %in% which(regimes == reg))
    if (reg == 0) edg <- c(edg, nEdges+1)
    return((sum(numerateur[edg])) / (sum(denominateur[edg])))
  }
  regimes_values <- 0:(length(unique(regimes))-1)
  beta_values <- sapply(regimes_values, fun)
  ## Computation of corresponding shifts values
  betas <- beta_values[regimes + 1]
  shifts <- compute_shifts_from_betas(phylo, betas)
  ## Computation of costs
  costs <- (1 - ee^2)^(-1) * (diff_exp - betas[daughters] * (1-ee))^2
  costs <- c((conditional_law_X$expectations[ntaxa+1] - beta_values[1])^2, costs)
  return(list(beta_0 = beta_values[1], shifts = shifts, costs = costs))
}

{ ##
# compute_diff_exp (phylo, Y_data, conditional_law_X)
# PARAMETERS:
#            @phylo (tree) imput tree
#            @Y_data (vector) : vector indicating the data at the tips
#            @conditional_law_X (list) result of compute_E (see note above)
# RETURNS:
#            (vector) entry i is Y_i-E[Z_pa(i)|Y] if i is a tip, and E[Z_i|Y]-E[Z_pa(i)|Y] if i is a node
# DEPENDENCIES:
#            none
# PURPOSE:
#            compute differences of conditionnal expectations
# NOTES:
#            none
# REVISIONS:
#            22/05/14 - Initial release
} ##
compute_diff_exp.BM <- function(phylo, conditional_law_X) {
  diff_exp <- rep(NA, nrow(phylo$edge))
  daughters <- phylo$edge[,2]
  parents <- phylo$edge[,1]
  diff_exp <- conditional_law_X$expectations[daughters] - conditional_law_X$expectations[parents]
  return(diff_exp)
}

compute_diff_exp.OU <- function(phylo, conditional_law_X, selection.strength) {
  diff_exp <- rep(NA, nrow(phylo$edge))
  daughters <- phylo$edge[,2]
  parents <- phylo$edge[,1]
  diff_exp <- conditional_law_X$expectations[daughters] - exp(-selection.strength * phylo$edge.length) * conditional_law_X$expectations[parents]
  return(diff_exp)
}

{ ##
# compute_var_diff (phylo, conditional_law_X)
# PARAMETERS:
#            @phylo (tree) imput tree
#            @conditional_law_X (list) result of compute_E (see note above)
# RETURNS:
#            (vector) entry i is Var[Z_pa(i)|Y] if i is a tip, and Var[Z_i-Z_pa(i)|Y] if i is a node
# DEPENDENCIES:
#            none
# PURPOSE:
#            compute conditionnal variances of differences
# NOTES:
#            none
# REVISIONS:
#            22/05/14 - Initial release
} ##
compute_var_diff.BM <- function(phylo, conditional_law_X) {
  var_diff <- rep(NA, nrow(phylo$edge))
  daughters <- phylo$edge[,2]
  parents <- phylo$edge[,1]
  var_diff <- conditional_law_X$variances[daughters] + conditional_law_X$variances[parents] - 2 * conditional_law_X$covariances[daughters]
  return(var_diff)
}

compute_var_diff.OU <- function(phylo, conditional_law_X, selection.strength) {
  var_diff <- rep(NA, nrow(phylo$edge))
  daughters <- phylo$edge[,2]
  parents <- phylo$edge[,1]
  ee <- exp(- selection.strength * phylo$edge.length)
  var_diff <- conditional_law_X$variances[daughters] + ee^2 * conditional_law_X$variances[parents] - 2 * ee * conditional_law_X$covariances[daughters]
  return(var_diff)
}

{ ##
#' @title Computation of the variance.
#'
#' @description
#' \code{compute_var_M.BM} finds the variance that is the maximum of likelihood
#'
#' @details
#' Given the variances, the costs0 and the edges where the shifts occurs, the
#' computation of the maximum of likelihood in the variance is simple.
#'
#' @param phylo Tree
#' @param var_diff variances of diferences
#' @param costs0 Cost of each edge
#' @param edges_max Edges where the shifts occur
#' 
#' @return a double : the computed variance
#'
#'02/06/14 - Initial release
} ##
compute_var_M.BM <- function(phylo, var_diff, costs0, edges_max){
  ntaxa <- length(phylo$tip.label)
  nNodes <- phylo$Nnode
  return(1/(ntaxa + nNodes-1) * ( sum(1/(phylo$edge.length) * var_diff) + sum( costs0[-edges_max]) ))
}

compute_var_M.OU.specialCase <- function(phylo, var_diff, costs, selection.strength, conditional_root_variance){
  ntaxa <- length(phylo$tip.label)
  nNodes <- phylo$Nnode
  ee <- exp(- selection.strength * phylo$edge.length )
  #return(1/(ntaxa + nNodes) * ( conditional_root_variance + sum((1 - ee^2)^(-1) * var_diff ) + sum( costs0[-edges_max]) ))
  return(1/(ntaxa + nNodes) * ( conditional_root_variance + sum((1 - ee^2)^(-1) * var_diff ) + sum(costs) ))
}

{ ##
# shutoff.EM.BM (params_old,params,tol)
# PARAMETERS:
#            @params_old (list) : old parameters
#            @params (list) : new parameters
#            @tol (list) : tolerance
# RETURNS: 
#            (bool) should we shutoff ?
# DEPENDENCIES:
#            none
# PURPOSE:
#            shutoff?
# NOTES:
#            TO BE DEFINED
# REVISIONS:
#            22/05/14 - Initial release
} ##
shutoff.EM.BM <- function(params_old, params, tol) {
  if (params_old$root.state$random) {
    return(shutoff.EM.BM.randroot(params_old,params,tol))
  } else {
    return(shutoff.EM.BM.fixedroot(params_old,params,tol))
  }
}

shutoff.EM.BM.randroot <- function(params_old, params, tol){
  if (abs(params_old$variance-params$variance)<tol$variance &&
      abs(params_old$root.state$exp.root-params$root.state$exp.root)<tol$exp.root &&
      abs(params_old$root.state$var.root-params$root.state$var.root)<tol$var.root) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

shutoff.EM.BM.fixedroot <- function(params_old, params, tol){
  if (abs(params_old$variance-params$variance)<tol$variance &&
      abs(params_old$root.state$value.root-params$root.state$value.root)<tol$value.root) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

shutoff.EM.OU <- function(stationnary.root, shifts_at_nodes, alpha_known){
  if (stationnary.root && shifts_at_nodes && alpha_known) {
    return(shutoff.EM.OU.specialCase)
  } else if (stationnary.root && shifts_at_nodes) {
    return(shutoff.EM.OU.stationnary.root_AND_shifts_at_nodes)
  } else {
    stop("The EM algorithm for the OU is only defined (for the moment) for a stationnary root and shifts at nodes !")
  }
}

shutoff.EM.OU.specialCase <- function(params_old, params, tol){
  if (abs(params_old$variance-params$variance)<tol$variance &&
        abs(params_old$root.state$exp.root-params$root.state$exp.root)<tol$exp.root &&
        abs(params_old$root.state$var.root-params$root.state$var.root)<tol$var.root) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

shutoff.EM.OU.stationnary.root_AND_shifts_at_nodes <- function(params_old, params, tol){
  if (abs(params_old$variance-params$variance)<tol$variance &&
        abs(params_old$root.state$exp.root-params$root.state$exp.root)<tol$exp.root &&
        abs(params_old$root.state$var.root-params$root.state$var.root)<tol$var.root &&
        abs(params_old$selection.strength-params$selection.strength)<tol$selection.strength) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

{ ##
#' @title Check wether parameters are finite.
#'
#' @description
#' \code{is.finite.params} checks whether calculated parameters in the EM are
#' finite or not.
#'
#' @details
#' This function is used to test the convergence of the algorithm.
#'
#' @param params list of parameters with the correct structure
#' 
#' @return boolean
#'
#'10/06/14 - Initial release
} ##

is.finite.params.BM <- function(params) {
  if (params$root.state$random) {
    return(is.finite.params.BM.randroot(params))
  } else {
    return(is.finite.params.BM.fixedroot(params))
  }
}

is.finite.params.BM.randroot <- function(params) {
  if (is.finite(params$variance) &&
        is.finite(params$root.state$exp.root) &&
        is.finite(params$root.state$var.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.finite.params.BM.fixedroot <- function(params) {
  if ( is.finite(params$variance) &&
       is.finite(params$root.state$value.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.finite.params.OU <- function(stationnary.root, shifts_at_nodes, alpha_known){
  if (stationnary.root && shifts_at_nodes && alpha_known) {
    return(is.finite.params.OU.specialCase)
  } else if (stationnary.root && shifts_at_nodes) {
    return(is.finite.params.OU.stationnary.root_AND_shifts_at_nodes)
  } else {
    stop("The EM algorithm for the OU is only defined (for the moment) for a stationnary root and shifts at nodes !")
  }
}

is.finite.params.OU.specialCase <- function(params) {
  if (is.finite(params$variance) &&
        is.finite(params$root.state$exp.root) &&
        is.finite(params$root.state$var.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.finite.params.OU.stationnary.root_AND_shifts_at_nodes <- function(params) {
  if (is.finite(params$variance) &&
        is.finite(params$root.state$exp.root) &&
        is.finite(params$root.state$var.root) &&
        is.finite(params$selection.strength)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

{ ##
#' @title Check wether parameters are in ranges.
#'
#' @description
#' \code{is.in.ranges.params} checks whether calculated parameters in the EM are
#' in the defined ranges.
#'
#' @details
#' This function is used to test the convergence of the algorithm.
#'
#' @param params list of parameters with the correct structure
#' @param min list of minimum values for the parameters
#' @param max list of maximum values for the parameters
#' 
#' @return boolean
#'
#'16/07/14 - Initial release
} ##

is.in.ranges <- function(p, min, max){
  if (p < min || p > max){
    return(FALSE)
  } else {
    return(TRUE)
  }
  
}

is.in.ranges.params.BM <- function(params, min, max) {
  if (params$root.state$random) {
    return(is.in.ranges.params.BM.randroot(params, min, max))
  } else {
    return(is.in.ranges.params.BM.fixedroot(params, min, max))
  }
}

is.in.ranges.params.BM.randroot <- function(params, min, max) {
  if (is.in.ranges(params$variance, min$variance, max$variance) &&
        is.in.ranges(params$root.state$exp.root, min$exp.root, max$exp.root) &&
        is.in.ranges(params$root.state$var.root, min$var.root, max$var.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.in.ranges.params.BM.fixedroot <- function(params, min, max) {
  if ( is.in.ranges(params$variance, min$variance, max$variance) &&
         is.in.ranges(params$root.state$value.root, min$value.root, max$value.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.in.ranges.params.OU <- function(stationnary.root, shifts_at_nodes, alpha_known){
  if (stationnary.root && shifts_at_nodes && alpha_known) {
    return(is.in.ranges.params.OU.specialCase)
  } else if (stationnary.root && shifts_at_nodes) {
    return(is.in.ranges.params.OU.stationnary.root_AND_shifts_at_nodes)
  } else {
    stop("The EM algorithm for the OU is only defined (for the moment) for a stationnary root and shifts at nodes !")
  }
}

is.in.ranges.params.OU.specialCase <- function(params, min, max) {
  if (is.in.ranges(params$variance, min$variance, max$variance) &&
        is.in.ranges(params$root.state$exp.root, min$exp.root, max$exp.root) &&
        is.in.ranges(params$root.state$var.root, min$var.root, max$var.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.in.ranges.params.OU.stationnary.root_AND_shifts_at_nodes <- function(params, min, max) {
  if (is.in.ranges(params$variance, min$variance, max$variance) &&
        is.in.ranges(params$root.state$exp.root, min$exp.root, max$exp.root) &&
        is.in.ranges(params$root.state$var.root, min$var.root, max$var.root) &&
        is.in.ranges(params$selection.strength, min$selection.strength, max$selection.strength)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

{ ##
# compute_MaxCompleteLogLik.BM (phylo, params)
# PARAMETERS:
#            @params (list) : parameters
# RETURNS: 
#            (double) : Value of the maximum of E[log(p_theta(X))|Y] given parameters theta
# DEPENDENCIES:
#            none
# PURPOSE:
#            Compute the value of the maximum
# NOTES:
#            
# REVISIONS:
#            24/05/14 - Initial release
#            02/06/14 - Change constant
} ##
compute_MaxCompleteLogLik.BM <- function(phylo, params) {
  ntaxa <- attr(params, "ntaxa")
  nNodes <- phylo$Nnode
  cst <- -1/2 * ((nNodes + ntaxa) * log(2*pi)) - sum(log(phylo$edge.length))
  if (params$root.state$random) {
    return( cst - 1/2 * (log(params$root.state$var.root) + (nNodes) * log(params$variance) + 1) )
  } else {
    return(cst - 1/2 * ((nNodes) * log(params$variance) + 1))
  }
}

compute_MaxCompleteLogLik.OU <- function(stationnary.root, shifts_at_nodes){
  if (stationnary.root && shifts_at_nodes) {
    return(compute_MaxCompleteLogLik.OU.stationnary.root_AND_shifts_at_nodes)
  } else {
    stop("The EM algorithm for the OU is only defined (for the moment) for a stationnary root and shifts at nodes !")
  }
}

compute_MaxCompleteLogLik.OU.stationnary.root_AND_shifts_at_nodes <- function(phylo, params) {
  ntaxa <- attr(params, "ntaxa")
  nNodes <- phylo$Nnode
  ee <- exp(- params$selection.strength * phylo$edge.length)
  cst <- -1/2 * ((nNodes + ntaxa) * log(2*pi)) - sum(log(1-ee^2))
  if (params$root.state$random) {
    return( cst - 1/2 * (log(params$root.state$var.root) + (nNodes) * log(params$variance) + 1) )
  } else {
    return(cst - 1/2 * ((nNodes) * log(params$root.state$var.root) + 1))
  }
}

compute_LogLikCondToShifts <- function(Y_data, m_Y, detSigmaYY, Sigma_YY_inv){
  ntaxa <- length(m_Y)
  return(n*log(2*pi) + log(detSigmaYY) + t(Y_data - m_Y)%*%(Y_data - m_Y))
}

compute_LogLikCondToShifts <- function(phylo, Y_data, process, shifts.edges, random.root, stationary.root=NULL, alpha){
  ntaxa <- length(Y_data)
  shifts <- compute_shiftsFromAlpha(alpha, shifts.edges)
  root.state = test.root.state(list(random=random.root, 
                                    stationary.root=stationary.root, 
                                    value.root=theta_0, 
                                    exp.root=theta_0, 
                                    var.root= ))
  sim <- simulate(phylo=phylo, 
                  process=process, 
                  root.state = root.state, 
                  shifts = shifts, 
                  variance=0, 
                  optimal.value=theta_1, 
                  selection.strength=alpha)
  m_Y <- extract.simulate(sim, where="tips", what="expectations")
  return(n*log(2*pi) + log(detSigmaYY) + t(Y_data - m_Y)%*%(Y_data - m_Y))
}

compute_LogLikdata <- function(phylo, Y_data, process, params, Sigma_YY){
  ntaxa <- length(Y_data)
  sim <- simulate(phylo = phylo, 
                  process = process, 
                  root.state = params$root.state, 
                  shifts = params$shifts, 
                  variance = params$variance, 
                  optimal.value = params$optimal.value, 
                  selection.strength = params$selection.strength)
  m_Y <- extract.simulate(sim, where="tips", what="expectations")
  return(n*log(2*pi) + log(detSigmaYY) + t(Y_data - m_Y)%*%(Y_data - m_Y))
}

{ ##
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
} ##

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

## make.name ##################################################################
###############################################################################

{ ##
# make.name (Name, TreeType, Y.state, Z.state, phylo, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, ...)
# PARAMETERS:
#            @process (string) : which process ?
#            @paramsSimu (list) : parameters used for th simulation
#            @paramsEstimate (list) : parameters found by an (optional) estimation from the data
#            @estimate (bool) : wether the data is issued from an estimation or a direct simulation
#            @...
# RETURNS:
#            (string) standardized name for the data
# DEPENDENCIES:
#            catch.ProcessParams, catch.TolParams
# PURPOSE:
#            Generate a standardized name for the data
# NOTES:
#            none
# REVISIONS:
#            26/05/14 - Initial release
} ##
make.name <- function(process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, params_algo_EM=NULL, ...) {
  ## Choose Process
  catch.ProcessParams  <- switch(process, 
                                 BM = catch.ProcessParams.BM, 
                                 OU = catch.ProcessParams.OU)
  catch.TolParams  <- switch(process, 
                             BM = catch.TolParams.BM, 
                             OU = catch.TolParams.OU)
  ## Define File name
  RootState <- paste("_root-rand=",paramsSimu$root.state$random,"_val-root=",paramsSimu$root.state$value.root,"_exp-root=",paramsSimu$root.state$exp.root,"_var-root=",paramsSimu$root.state$var.root,sep="")
  if (is.null(paramsEstimate$shifts$edges)){
    ShiftsState <- "_no-shift"
  } else {
    ShiftsState <- paste("_shifts-edges=", paste(paramsSimu$shifts$edges,collapse="-"), "_shifts-val=", paste(paramsSimu$shifts$values, collapse="-"), "_shifts-T=", paste(paramsSimu$shifts$relativeTimes,collapse="-"), sep="")
  }
  ProcessParams <- catch.ProcessParams(paramsSimu)
  ## Parametrers of the EM if relevent
  if (estimate) {
    TolParams <- "" # catch.TolParams(params_algo_EM)
    EstimParams <- paste(TolParams, "_process-used=", params_algo_EM$process, "_met-variance=", params_algo_EM$method.variance, "_met-init=", params_algo_EM$method.init, "_nbr-shifts=", params_algo_EM$nbr_of_shifts, sep="")
  } else {
    EstimParams <- ""
  }
  return(paste(ProcessParams, RootState, ShiftsState, EstimParams, sep=""))
}

{ ##
# plot.process (Name, TreeType, Y.state, Z.state, phylo, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, ...)
# PARAMETERS:
#            @Name (string) begening of the name
#            @TreeType (string) : type pf tree
#            @Y.state (vector) : value of the traits at the tips
#            @Z.state (vector) : value of the traits at the internal nodes
#            @process (string) : which process ?
#            @paramsSimu (list) : parameters used for th simulation
#            @paramsEstimate (list) : parameters found by an (optional) estimation from the data
#            @estimate (bool) : wether the data is issued from an estimation or a direct simulation
#            @directory (string) : where to store the plot
#            @...
# RETURNS:
#            (pdf) a pdf plot of the process
# DEPENDENCIES:
#            catch.LegendProcess, make.name
# PURPOSE:
#            Generate a pdf file with a plot of the process
# NOTES:
#            none
# REVISIONS:
#            26/05/14 - Initial release
#            06/06.14 - Add directory
} ##
plot.process <- function(Name, TreeType, Y.state, Z.state, phylo, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, position.legend="bottomleft", directory, params_algo_EM=NULL) {
  ## Define File Name
  FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, params_algo_EM)
  ## Define legend
  LegendProcess <- catch.LegendProcess(process, paramsEstimate, estimate, params_algo_EM)
  LegendRoot <- c(#paste("Random Root = ",round(paramsEstimate$root.state$random,2),sep=""),
                  #paste("Root Value (if not random) = ",round(paramsEstimate$root.state$value.root,2),sep=""),
                  paste("Root expectation = ",round(paramsEstimate$root.state$exp.root,2),sep=""),
                  paste("Root variance = ",round(paramsEstimate$root.state$var.root,2), sep=""))
  ## Plot
  FileName <- paste(Name, TreeType, FileName, sep="")
  pdf(paste(directory, FileName, ".pdf",sep=""), height=10,width=20)
  plot(phylo, show.tip.label = FALSE)
  tiplabels(pch = 19, cex = abs(Y.state)/mean(abs(Y.state)), col = ifelse(Y.state >= 0, "orangered", "lightblue"))
  nodelabels(pch = 19, cex = abs(Z.state)/mean(abs(Z.state)), col = ifelse(Z.state >= 0, "orangered", "lightblue"))
  if ( !is.null(paramsEstimate$shifts$edges) ) {
    edgelabels(text=round(paramsEstimate$shifts$values,2), edge=paramsEstimate$shifts$edges, bg="chocolate4", cex=2.5)
  }
  legend(paste(position.legend),legend=c(LegendProcess,LegendRoot), col="black", cex = 2)
  dev.off()
}

save.process <- function(Name, TreeType, XX, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, directory, ...) {
  ## Define File Name
  FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, ...)
  FileName <- paste(Name, TreeType, FileName, sep="")
  save(XX, paramsSimu, paramsEstimate, file=paste(directory, FileName, ".RData",sep=""))
}

write.table.results <- function(Name, TreeType, res, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, directory, ...) {
  ## Define File Name
  FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, ...)
  FileName <- paste(Name, TreeType, FileName, sep="")
  write.table(res, paste(directory, FileName, ".csv",sep=""))
}

load.process <- function(Name, TreeType, XX, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, directory, ...) {
  ## Define File Name
  FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, ...)
  FileName <- paste(Name, TreeType, FileName, sep="")
  load(file=paste(directory, FileName, ".RData",sep=""))
}

{ ##
# catch.ProcessParams (paramsSimu)
# PARAMETERS:
#            @paramsSimu (list) : parameters used for th simulation
# RETURNS:
#            (string) string containing all the informations on the parameters used for the simulation
# DEPENDENCIES:
#            none
# PURPOSE:
#            generate adequate string
# NOTES:
#            none
# REVISIONS:
#            26/05/14 - Initial release
} ##
catch.ProcessParams.BM <- function(paramsSimu){
  return(paste("_process=BM","_var=",paramsSimu$variance,sep=""))
}

catch.ProcessParams.OU <- function(paramsSimu){
  return(paste("_process=OU","_var=",paramsSimu$variance,"_opt-val=",paramsSimu$optimal.value,"_sel-strength=",paramsSimu$selection.strength,sep=""))
}

{ ##
# catch.LegendProcess (paramsEstimate)
# PARAMETERS:
#            @paramsEstimate (list) : parameters estimated from the data
# RETURNS:
#            (vector) strings containing all the informations on the parameters used for the simulation
# DEPENDENCIES:
#            none
# PURPOSE:
#            generate adequate vector of string for a legend
# NOTES:
#            none
# REVISIONS:
#            26/05/14 - Initial release
} ##
catch.LegendProcess <- function(process, paramsEstimate, estimate, params_algo_EM=NULL) {
  if (estimate) {
    proc <- params_algo_EM$process
  } else {
    proc <- process
  }
  if (proc=="BM") {
    return(catch.LegendProcess.BM(paramsEstimate))
  } else if (proc=="OU") {
    return(catch.LegendProcess.OU(paramsEstimate))
  }
}
catch.LegendProcess.BM <- function(paramsEstimate){
  return( c("Process : BM",
            paste("Process Variance = ", round(paramsEstimate$variance,2), sep="")) )
}

catch.LegendProcess.OU <- function(paramsEstimate){
  return( c("Process : OU",
            paste("Process Variance = ", round(paramsEstimate$variance,2),sep=""),
            paste("Beta_0 = ", round(paramsEstimate$optimal.value,2), sep=""),
            paste("Selection Strength = ", round(paramsEstimate$selection.strength,2), sep="")) )
}

{ ##
# catch.TolParams (params_algo_EM)
# PARAMETERS:
#            @params_algo_EM (list) : parameters of the EM algorithm used
# RETURNS:
#            (string) string containing all the informations on the parameters used for the tolerence of the parameters in the estimation
# DEPENDENCIES:
#            none
# PURPOSE:
#            generate adequate string
# NOTES:
#            none
# REVISIONS:
#            26/05/14 - Initial release
} ##
catch.TolParams.BM <- function(params_algo_EM){
  return(paste("_tol-variance=", params_algo_EM$tol$variance,
               "_tol-exp-root=", params_algo_EM$tol$exp.root,
               "_tol-var-root=", params_algo_EM$tol$var.root,sep=""))
}

catch.TolParams.OU <- function(params_algo_EM){
  return(paste("_tol-variance=", params_algo_EM$tol$variance,
               "_tol-exp-root=", params_algo_EM$tol$exp.root,
               "_tol-var-root=", params_algo_EM$tol$var.root,
               "_tol-optim-value=", params_algo_EM$tol$optim.value,
               "_tol-selection-strength=", params_algo_EM$tol$selection.strength,sep=""))
}

## incidence.matrix ###########################################################
###############################################################################

{ ##
#' @title Incidence matrix of a tree.
#'
#' @description
#' \code{incidence.matrix} computes the incidence matrix T of a tree : for a 
#' lineage i and a branch b, T[i,b]=1 if b is in the lineage i, and 0 
#' otherwise.
#'
#' @details
#' This function uses the general up tree recursion function \code{recursionUp}
#' with the initialization function \code{init.incidence.matrix} and update
#' function \code{update.incidence.matrix}.
#'
#' @param phylo Input tree.
#' 
#' @return Matrix of incidence.
#'
#'17/06/14 - Initial release
} ##
incidence.matrix <- function(phylo){
  ## Reorder and keep track
  phy <- reorder(phylo, order = "postorder")
  cor <- correspondanceEdges(edges=1:nrow(phy$edge), from=phylo, to=phy)
  ## Init and recurence
  T <- init.incidence.matrix(phy)
  T <- recursionUp(phy, T, update.incidence.matrix)
  ## Take for each node its parenting branch
  daughters <- phy$edge[,2]
  T <- T[daughters,]
  ## Return to original tree
  T <- T[cor, ]
  return(t(T))
}

{ ##
#' @title Initialization for incidence matrix
#'
#' @description
#' \code{init.incidence.matrix} initialize the matrix updated in
#' \code{update.incidence.matrix} for the computation of the incidence matrix
#' in \code{incidence.matrix}.
#'
#' @details
#' The initialized matrix has ntaxa column and nNodes rows. Each node
#' represent its parental branch. A row corresponding to a tip i is initialized
#' to a vector of zeros, with only entry i equal to one. (Branch ending at 
#' tip i is only in the i^th lineage)
#'
#' @param phy Input tree.
#' 
#' @return Matrix with nNodes rows and ntaxa column.
#'
#'17/06/14 - Initial release
} ##
init.incidence.matrix <- function(phy){
  ntaxa <- length(phy$tip.label)
  T <- matrix(NA, nrow = 1 + nrow(phy$edge), ncol = ntaxa)
  for (i in 1:ntaxa){
    T[i,] <- 1:ntaxa == i
  }
  return(T)
}

{ ##
#' @title Update function for incidence matrix
#'
#' @description
#' \code{update.incidence.matrix} updates the matrix initialized in 
#' \code{init.incidence.matrix} for the computation of the incidence matrix
#' in \code{incidence.matrix}.
#'
#' @details
#' A node belongs to all the lineages of its daughters.
#'
#' @param daughtersParams : rows of updated matrix corresponding to the
#' daughters of the current node.
#' 
#' @return Vector of length ntaxa, indicating to which lineages the branch
#' above the current node belongs to.
#'
#'17/06/14 - Initial release
} ##
update.incidence.matrix <- function(daughtersParams, ...){
  inc <- colSums(daughtersParams)
  return(inc>0)
}

{ ##
#' @title Incidence matrix of a tree.
#'
#' @description
#' \code{incidence.matrix.full} computes the incidence matrix T of a tree : for a 
#' node i and a branch b, T[i,b]=1 if b is in the lineage i, and 0 
#' otherwise.
#'
#' @details
#' This function uses the general up tree recursion function \code{recursionUp}
#' with the initialization function \code{init.incidence.matrix.full} and update
#' function \code{update.incidence.matrix.full}.
#'
#' @param phylo Input tree.
#' 
#' @return Matrix of incidence, size ntaxa + nNodes
#'
#'06/10/14 - Initial release
} ##
incidence.matrix.full <- function(phylo){
  ## Reorder and keep track
  phy <- reorder(phylo, order = "postorder")
  cor <- correspondanceEdges(edges=1:nrow(phy$edge), from=phylo, to=phy)
  ## Init and recurence
  U <- init.incidence.matrix.full(phy)
  U <- recursionUp(phy, U, updateUp = update.incidence.matrix.full)
  ## Take for each node its parenting branch
  daughters <- phy$edge[,2]
  U <- U[daughters,]
  ## Return to original tree
  U <- U[cor, ]
  return(t(U))
}

{ ##
#' @title Initialization for incidence matrix (full tree)
#'
#' @description
#' \code{init.incidence.matrix.full} initialize the matrix updated in
#' \code{update.incidence.matrix.full} for the computation of the incidence
#'  matrix of the full tree in \code{incidence.matrix.full}.
#'
#' @details
#' The initialized matrix is squared of size ntaxa + nNodes. Each node
#' represent its parental branch. A row corresponding to a tip i is initialized
#' to a vector of zeros, with only entry i equal to one. (Branch ending at 
#' tip i is only in the i^th lineage)
#'
#' @param phy Input tree.
#' 
#' @return Matrix of size ntaxa + nNodes.
#'
#'06/10/14 - Initial release
} ##
init.incidence.matrix.full <- function(phy){
  ntaxa <- length(phy$tip.label)
  nNodes <- phy$Nnode
  U <- matrix(NA, nrow = nNodes, ncol = ntaxa + nNodes)
  T_1 <- matrix(0, nrow = ntaxa, ncol = ntaxa + nNodes)
  U <- rbind(T_1, U)
  diag(U) <- rep(1, ntaxa + nNodes)
  return(U)
}

{ ##
#' @title Update function for incidence matrix
#'
#' @description
#' \code{update.incidence.matrix.full} updates the matrix initialized in 
#' \code{init.incidence.matrix.full} for the computation of the incidence matrix
#' in \code{incidence.matrix.full}.
#'
#' @details
#' A node belongs to all the lineages of its daughters.
#'
#' @param daughtersParams : rows of updated matrix corresponding to the
#' daughters of the current node.
#' 
#' @return Vector of length ntaxa + nNodes, indicating to which lineages the 
#' branch above the current node belongs to.
#'
#'06/10/14 - Initial release
} ##
update.incidence.matrix.full <- function(daughtersParams, parent, ...){
  inc <- colSums(daughtersParams)
  inc[parent] <- 1
  return(inc > 0)
}

{ ##
#' @title Compute the vector of shifts.
#'
#' @description
#' \code{shifts.list_to_vector} takes the list description of the shifts 
#' to give the vectorial representation of the shifts : the b th element of 
#' the vector has the value of the shift occuring on that branch b.
#'
#' @details
#'
#' @param phy Input tree.
#' @param shifts : list description of the shifts : shifts$edges, shifts$values
#' 
#' @return Vector of length nbranch.
#' 
#' @seealso \code{shifts.vector_to_list}
#' 
#'17/06/14 - Initial release
} ##
shifts.list_to_vector <- function(phy, shifts){
  delta <- rep(0, nrow(phy$edge))
  delta[shifts$edges] <- shifts$values
  return(delta)
}

{ ##
#' @title Compute the list of shifts.
#'
#' @description
#' \code{shifts.vector_to_list} takes the vectorial description of the shifts 
#' to create the list description of the shifts.
#'
#' @details
#'
#' @param delta : vector description of the shift.
#' 
#' @return Vector of length nbranch.
#' 
#' @seealso \code{shifts.list_to_vector}
#' 
#'17/06/14 - Initial release
} ##
shifts.vector_to_list <- function(delta){
  edsh <- which(delta != 0)
  shifts <- list(edges = edsh, values = delta[edsh], relativeTimes = rep(0, length(edsh)))
  return(shifts)
}

{ ##
#' @title Initialisation for the computation of the optimal values
#'
#' @description
#' \code{init.compute_betas} initialize the vector of optimal values at nodes and
#' tips, with the value at the root.
#'
#' @details
#' This function is used in function \code{compute_betas} and is designed to 
#' furnish function \code{update.compute_betas} with the right structure of data.
#'
#' @param phy Input tree.
#' @param optimal.value the optimal value at the root of the tree
#' 
#' @return Matrix of size (nNodes + ntaxa)x1 of NAs, with the optiaml value
#'  at the root.
#'
#'06/10/14 - Initial release
} ##
init.compute_betas <- function(phy, optimal.value, ...){
  ntaxa <- length(phy$tip.label)
  beta <- matrix(nrow = 1 + nrow(phy$edge), ncol = 1) # selection strength
  beta[ntaxa + 1,] <- optimal.value
  return(beta)
}

{ ##
#' @title Update function ofr optimal value computation
#'
#' @description
#' \code{update.compute_betas} computes the optimal value at a daughter node, 
#' knowing the optimal value at the parent node and the vector of shifts.
#'
#' @details
#' This function is used in function \code{compute_betas} and is designed to 
#' furnish function \code{recursionDown} with the right structure of data.
#'
#' @param edgeNbr : Number of the edge considered
#' @param ancestral : Computed vector for the parental node
#' @param shifts position and values of the shifts 
#' 
#' @return Updated matrix of size (nNodes + ntaxa)x1.
#'
#'06/10/14 - Initial release
} ##
update.compute_betas <- function(edgeNbr, ancestral, shifts, ...){
  shiftsIndex <- which(shifts$edges == edgeNbr) #If no shifts = NULL, and sum = 0
  beta <- ancestral + sum(shifts$values[shiftsIndex])
  return(beta)
}


{ ##
#' @title Computation of the optimal values at nodes and tips.
#'
#' @description
#' \code{compute_betas} computes the optimal values at the nodes and tips of the
#' tree, given the value at the root and the list of shifts occuring in the tree.
#'
#' @details
#' This function uses function \code{recursionDown} for recursion on the tree, 
#' with \code{init.compute_betas} for the initialization of the vector, and 
#' \code{update.compute_betas} for its actualization.
#'
#' @param phylo : imput tree
#' @param optimal.value the optimal value at the root of the tree
#' @param shifts position and values of the shifts 
#' 
#' @return Vector of size (ntaxa + nNodes) of the ptimal values at the tips
#' of the tree.
#'
#'06/10/14 - Initial release
} ##
compute_betas <- function(phylo, optimal.value, shifts){
  phy <- reorder(phylo, order = "cladewise")
  ## Trace edges
  shifts_ordered <- shifts
  shifts_ordered$edges <- correspondanceEdges(edges=shifts$edges,from=phylo,to=phy)
  betas <- init.compute_betas(phy, optimal.value)
  betas <- recursionDown(phy, betas, update.compute_betas, shifts_ordered)
  return(betas)
}

{ ##
#' @title Initialisation for the allocation of shifts.
#'
#' @description
#' \code{init.allocate_regimes_from_shifts} initialize the vector of regimes 
#' at nodes and tips, with the regime (0) at the root.
#'
#' @details
#' This function is used in function \code{allocate_regimes_from_shifts} and is 
#' designed to furnish function \code{update.allocate_regimes_from_shifts} with 
#' the right structure of data.
#'
#' @param phy Input tree.
#' 
#' @return Matrix of size (nNodes + ntaxa)x1 of NAs, with the 0 at the root.
#'
#'10/10/14 - Initial release
} ##
init.allocate_regimes_from_shifts <- function(phy, ...){
  ntaxa <- length(phy$tip.label)
  regimes <- matrix(nrow = 1 + nrow(phy$edge), ncol = 1) 
  regimes[ntaxa + 1,] <- 0
  return(regimes)
}

{ ##
#' @title Update function for regime allocation.
#'
#' @description
#' \code{update.allocate_regimes_from_shifts} computes the regime of a daughter 
#' node, knowing the regime at the parent node and the vector of shifts positions
#'
#' @details
#' This function is used in function \code{allocate_regimes_from_shifts} and is 
#' designed to furnish function \code{recursionDown} with the right structure of
#'  data.
#'
#' @param edgeNbr : Number of the edge considered
#' @param ancestral : regime of the parent node
#' @param shifts_edges positions on edges
#' 
#' @return regime of the daughter node.
#'
#'10/10/14 - Initial release
} ##
update.allocate_regimes_from_shifts <- function(edgeNbr, ancestral, shifts_edges, ...){
  shiftsIndex <- which(shifts_edges == edgeNbr) # If no shifts = integer(0)
  if (length(shiftsIndex) == 0){
    return(ancestral) # No shift : keep same regime
  } else {
    return(shiftsIndex) # Shift : new regime numbered by the position in the list
  }
}


{ ##
#' @title Allocation of regimes to nodes.
#'
#' @description
#' \code{allocate_regimes_from_shifts} allocate a number (from 0 to the number 
#' of shifts) to each node, corresponding to its regime : all nodes below shift 
#' i are numbered by i.
#' 
#' @details
#' This function uses function \code{recursionDown} for recursion on the tree, 
#' with \code{init.allocate_regimes_from_shifts} for the initialization of the
#' vector, and \code{update.allocate_regimes_from_shifts} for its actualization.
#'
#' @param phylo : imput tree
#' @param shifts_edges edges were the shifts are
#' 
#' @return Vector of size (ntaxa + nNodes) of the regimes of each node and tip.
#'
#'10/10/14 - Initial release
} ##
allocate_regimes_from_shifts <- function(phylo, shifts_edges){
  phy <- reorder(phylo, order = "cladewise")
  ## Trace edges
  shifts_ordered <- correspondanceEdges(edges=shifts_edges, from=phylo, to=phy)
  regimes <- init.allocate_regimes_from_shifts(phy, shifts_ordered)
  regimes <- recursionDown(phy, regimes, update.allocate_regimes_from_shifts, shifts_ordered)
  return(regimes)
}


{ ##
#' @title Allocation of shifts to edges
#'
#' @description
#' \code{allocate_shifts_from_regimes} returns the position of the shifts induced
#' by the allocation of the regimes. Only works in an "infinite site" model.
#' 
#' @details
#' This function uses function fun on each row of matrix of edges.
#'
#' @param phylo : imput tree
#' @param regimes : vector of size (ntaxa + nNodes) of the regimes of each node
#' and tip.
#' 
#' @return vector of edges numbers where the shifts are
#'
#'10/10/14 - Initial release
} ##
allocate_shifts_from_regimes <- function(phylo, regimes){
  fun <- function(edge){
    if (regimes[edge[1]] == regimes[edge[2]]) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
  shifts_edges <- apply(phylo$edge, 1, fun)
  return(which(shifts_edges))
}

{ ##
#' @title Computation of shifts from the vector of optimal values
#'
#' @description
#' \code{compute_shifts_from_betas} computes the list of shifts corresponding
#' to the vector of optimal values on nodes.
#' 
#' @details
#' This function uses function fun on each row of matrix of edges.
#'
#' @param phylo : imput tree
#' @param betas : vector of size (ntaxa + nNodes) of the optimal values at each
#' node and tip.
#' 
#' @return list of shifts
#'
#'13/10/14 - Initial release
} ##
compute_shifts_from_betas <- function(phylo, betas){
  fun <- function(edge){
    if (betas[edge[1]] == betas[edge[2]]) {
      return(c(0, 0))
    } else {
      return(c(1, betas[edge[2]] - betas[edge[1]]))
    }
  }
  shifts_ed <- apply(phylo$edge, 1, fun)
  edges <- which(shifts_ed[1,] == 1)
  shifts <- list(edges = edges,
                 values = shifts_ed[2, edges],
                 relativeTimes = rep(0, length(edges)))
  return(shifts)
}

## plot.history ###############################################################
###############################################################################

list_to_table.history <- function(params_history) {
  ll <- unlist(sapply(params_history, function(x) attr(x, "log_likelihood")[1]))
  ll_bis <- unlist(sapply(params_history, function(x) attr(x, "log_likelihood_bis")[1]))
  method <- unlist(sapply(params_history, function(x) attr(x, "segmentation_algorithm_used")))
  nbr_of_shifts <- length(params_history[['1']]$shifts$edges)
  params_history[['0']] <- replaceInList(params_history[['0']], function(x) if(is.null(x))rep(0,nbr_of_shifts) else x)
  params_history <- lapply(params_history, unlist)
  history <- do.call(cbind, params_history)
  history <- rbind(history, log_likelihood = c(ll, NA), log_likelihood_bis = c(ll_bis, NA), segmentation_algorithm = c(method, NA, NA))
  return(history)
}

write.table.history <- function(history, params_algo_EM, PATH, ...) {
  ## Define File Name
  name <- paste(PATH, "history_parameters", "_init=", params_algo_EM$method.init, "_initalpha=", params_algo_EM$method.init.alpha, "_nbrofshifts=", params_algo_EM$nbr_of_shifts, ".csv", sep="")
  write.csv2(history, name, ...)
}

plot.history.OU.stationnary <- function(params_history, paramsSimu, PATH, params_algo_EM, name){
  history <- list_to_table.history(params_history)
  params_simu  <-  unlist(paramsSimu)[names(history[,1])]
  pdf(paste(PATH, "history_parameters", "_init=", params_algo_EM$method.init, "_initalpha=", params_algo_EM$method.init.alpha, "_nbrofshifts=", params_algo_EM$nbr_of_shifts, ".pdf", sep=""), width = 12, height = 8)
  ## Create grid
  pushViewport(viewport(layout = grid.layout(2+1, 3, heights = unit(c(1,5,5), "null"))))
  ## Title of the page
  grid.text(paste("Initialization : ", params_algo_EM$method.init, "\n", "Alpha Initialization : ", params_algo_EM$method.init.alpha, sep=""), vp = vplayout(1,1:3))
  ## Continuous parameters
  row <- c(2,2,3,3); col <- c(1,2,1,2)
  params_to_plot <- c("variance", "selection.strength", "root.state.var.root", "optimal.value")
  params_to_plot_legend <- c(expression(sigma^2), expression(alpha), expression(gamma^2), expression(beta[0]))
  for (s in 1:4) {
    # Choose right score
    df <- history[params_to_plot[s],]
    # Plot
    p <- qplot(seq_along(df)-1, df, xlab="iterations", ylab=params_to_plot_legend[s])
    p <- p + geom_hline(yintercept=params_simu[params_to_plot[s]])
    print(p, vp=vplayout(row[s],col[s]))
  }
  ## Plot with differents colours for each shift
  # Values
  df_val <- as.data.frame(t(history[grepl('shifts.values', names(history[,1])),, drop=F]))
  df_val_long <- melt(df_val, variable.name = "shift", value.name="shift.value")
  df_val_long$iterations <- rep(seq_along(df_val[,1])-1, dim(df_val)[2])
  df_val_long$shift <- as.factor(rep(seq_along(df_val[1,]), each=dim(df_val)[1]))
  df_val_long$true <- rep(params_simu[grepl('shifts.values', names(history[,1]))], each=dim(df_val)[1])
  p <- ggplot(data=df_val_long, aes(x=iterations, y=shift.value, colour = shift)) + geom_point()
  p <- p + geom_hline(aes(yintercept=true, colour=shift), data=df_val_long)
  p <- p + ylab(expression(delta))
  print(p, vp=vplayout(2,3))
  # Edges
  df_ed <- as.data.frame(t(history[grepl('shifts.edges', names(history[,1])),, drop=F]))
  df_ed_long <- melt(df_ed, variable.name = "shift", value.name="shift.edge")
  df_ed_long$iterations <- rep(seq_along(df_ed[,1])-1, dim(df_ed)[2])
  df_ed_long$shift <- as.factor(rep(seq_along(df_ed[1,]), each=dim(df_ed)[1]))
  df_ed_long$true <- rep(params_simu[grepl('shifts.edges', names(history[,1]))], each=dim(df_ed)[1])
  p <- ggplot(data=df_ed_long, aes(x=iterations, y=shift.edge, colour = shift)) + geom_point()
  p <- p + geom_hline(aes(yintercept=true, colour=shift), data=df_ed_long)
  p <- p + ylab(expression(tau))
  print(p, vp=vplayout(3,3))
  dev.off()
}

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
dtoc <- function(x) gsub(".", ",", x, fixed=TRUE)
