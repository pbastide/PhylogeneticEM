# {Simulate}
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
## Here are functions used to simulate a process on the tree
## Dependencies : generic_functions.R
###############################################################################
##
# simulate (phylo, root.state = list(random=FALSE,value.root,exp.root,var.root), process = c("BM", "OU"), shifts = list(edges=NULL,values=NULL,relativeTimes=NULL), ...)
# PARAMETERS:
# @phylo (tree) Input tree
# @root.state (list) state of the root :
# random : random state (TRUE) or deterministic state (FALSE)
# value.root : if deterministic, value of the character at the root
# exp.root : if random, expectation of the character at the root
# var.root : if random, variance of the character at the root
# @process (string) Random process to simulate. Possible values :
# "BM" : Brownian Motion
# "OU" : Ornstein-Uhlenbeck
# @shifts (list) position and values of the shifts :
# edges : vector of id of edges where the shfts are
# values : vector of same dimension of values of the shifts on the edges
# relativeTimes : vector of same dimension of relative time of the shift form the parent node of edges
# @... other parameters : parameters of the process :
# for BM : variance
# for OU : variance, optimal.value, selection.strength
# RETURNS:
# (matrix) Each row i contains a vector of values computed for node i :
# BM : simulated state, expected value
# OU : simulated state, expected value, optimal value
# DEPENDENCIES:
# init, update, correspondanceEdges (, extract.simulate), check.selection.strength
# PURPOSE:
# Simulate a random process on the tree phylo. Return the simulated and expected values of all tips and nodes.
# NOTES:
# none
# REVISIONS:
# 16/05/14 - Initial release
# 20/05/14 - Gestion of edges (function correspondanceEdges)
# 21/05/14 - extraction of the recursion
# 16/06/14 - check.selection.strength
##
simulate <- function(phylo, process = c("BM", "OU"), root.state = list(random=FALSE, stationary.root=FALSE, value.root, exp.root, var.root), shifts = list(edges=NULL,values=NULL,relativeTimes=NULL), eps=10^(-6), selection.strength=NULL, variance=NULL, optimal.value=NULL) {
  phy <- reorder(phylo, order = "cladewise")
  ## Trace edges
  shifts_ordered <- shifts
  shifts_ordered$edges <- correspondanceEdges(edges=shifts$edges,from=phylo,to=phy)
  ## Set branch stochastic process
  process <- match.arg(process)
  process <- check.selection.strength(process=process, selection.strength=selection.strength, eps=eps) # if OU, check if selection.strength is not too low.
  init <- switch(process,
                 BM = init.simulate.BM,
                 OU = init.simulate.OU)
  updateDown <- switch(process,
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
  attr(paramSimu, "ntaxa") <- ntaxa
  return(paramSimu)
}
##
# init.simulate.StateAndExp (phy, root.state = list(random=FALSE,value.root,exp.root,var.root))
# PARAMETERS:
# @phy (tree) Input tree
# @root.state (list) state of the root : (see note above)
# RETURNS:
# (matrix) Matrix with as many rows as the cumulated number of tips and nodes in the tree phy, and 2 columns. Each row 1<=i<=ntaxa contains a vector of values computed for tip i. All other rows are set to NAs
# DEPENDENCIES:
# none
# PURPOSE:
# Initialize the matrix used in simulate, with simulated and expected values of the root node.
# NOTES:
# none
# REVISIONS:
# 16/05/14 - Initial release
##
init.simulate.StateAndExp <- function(phy,root.state){
  ntaxa <- length(phy$tip.label)
  paramSimu <- matrix(NA,nrow=1 + nrow(phy$edge),ncol=2)
  if (!root.state$random) { # The root is not random
    paramSimu[ntaxa + 1,] <- c(root.state$value.root,
                               root.state$value.root)
  } else { # The value of the root is random N(exp.root,var.root)
    paramSimu[ntaxa + 1,] <- c(root.state$exp.root + sqrt(root.state$var.root)*rnorm(1),
                               root.state$exp.root)
  }
  return(paramSimu)
}
##
# init.simulate.BM (phy, root.state = list(random=FALSE,value.root,exp.root,var.root))
# PARAMETERS:
# @phy (tree) Input tree
# @root.state (list) state of the root (see note above)
# @... other parameters : parameters of the process (see note above)
# RETURNS:
# (matrix) Matrix with as many rows as the cumulated number of tips and nodes in the tree phy, and 2 columns. Row number ntaxa+1 (root) contains a vector of values computed for the root node. All other rows are set to NAs
# DEPENDENCIES:
# initStateAndExp
# PURPOSE:
# Initialize the matrix used in simulate, with simulated and expected values of the root node.
# NOTES:
# none
# REVISIONS:
# 16/05/14 - Initial release
##
init.simulate.BM <- function(phy,root.state,...){
  return(init.simulate.StateAndExp(phy,root.state))
}
##
# init.simulate.OU (phy, root.state = list(random=FALSE,value.root,exp.root,var.root))
# PARAMETERS:
# @phy (tree) Input tree
# @root.state (list) state of the root (see note above)
# @... other parameters : parameters of the process (see note above)
# RETURNS:
# (matrix) Matrix with as many rows as the cumulated number of tips and nodes in the tree phy, and 2 columns. Row number ntaxa+1 (root) contains a vector of values computed for the root node. All other rows are set to NAs
# DEPENDENCIES:
# initStateAndExp
# PURPOSE:
# Initialize the matrix used in simulate, with simulated and expected values of the root node, and optimal value of the root node
# NOTES:
# none
# REVISIONS:
# 16/05/14 - Initial release
##
init.simulate.OU <- function(phy, root.state, optimal.value, ...){
  paramSimu <- init.simulate.StateAndExp(phy,root.state)
  ntaxa <- length(phy$tip.label)
  beta <- rep(NA, length = 1 + nrow(phy$edge)) # selection strength
  beta[ntaxa + 1] <- optimal.value
  paramSimu <- cbind(paramSimu,beta)
}
##
# update.simulate.BM (edgeNbr, ancestral, length, shifts, variance, ...)
# PARAMETERS:
# @edgeNbr (int) Number of the edge considered
# @ancestral (vector) Computed vector for the parental node
# @length (double) Length of the current edge
# @shifts (list) position and values of the shifts (see note above)
# @variance (double) variance of the BM
# RETURNS:
# (vector) simulated state, expected value for the daughter node
# DEPENDENCIES:
# none
# PURPOSE:
# Update the process in the case of the BM
# NOTES:
# none
# REVISIONS:
# 16/05/14 - Initial release
##
update.simulate.BM <- function(edgeNbr, ancestral, length, shifts, variance, ...){
  shiftsIndex <- which(shifts$edges==edgeNbr) # If no shifts = NULL, and sum = 0
  return(c( ancestral[1] + sum(shifts$values[shiftsIndex]) + sqrt(length*variance) * rnorm(1),
            ancestral[2] + sum(shifts$values[shiftsIndex]) ) )
}
##
# update.simulate.OU (edgeNbr, ancestral, length, shifts, variance, selection.strength, eps=10^(-6), ...)
# PARAMETERS:
# @(edgeNbr, ancestral, length, shifts, variance, selection.strength) : see note above
# @eps (double) : if the selection strenght is smaller than eps, simulate according to a BM instead of an OU
# RETURNS:
# (vector) simulated state, expected value, optimal value for the daughter node
# DEPENDENCIES:
# updateBM
# PURPOSE:
# Update the process in the case of the OU
# NOTES:
# none
# REVISIONS:
# 16/05/14 - Initial release
##
update.simulate.OU <- function(edgeNbr, ancestral, length, shifts, variance, selection.strength, ...){
  shiftsIndex <- which(shifts$edges==edgeNbr) # If no shifts = NULL, and sum = 0
  beta <- ancestral[3] + sum(shifts$values[shiftsIndex])
  ee <- exp(-selection.strength*length)
  ss <- sum(shifts$values[shiftsIndex]*( 1-exp( -selection.strength * length * (1-shifts$relativeTimes[shiftsIndex]) ) ))
  SimExp <- c( ancestral[3]*(1-ee) + ancestral[1]*ee + ss + sqrt(variance*(1-ee^2)/(2*selection.strength))*rnorm(1),
               ancestral[3]*(1-ee) + ancestral[2]*ee + ss )
  return(c(SimExp,beta))
}
##
# extract.simulate (paramSimu, where=c("tips","nodes"), what=c("states","expectations"))
# PARAMETERS:
# @paramSimu (matrix) return of the function simulate
# @where (string) : where to extract the values : at the "tips" or the internal "nodes" ?
# @what (string) : which value to extract : the simulated "states" or the "expectations" ?
# RETURNS:
# (vector) values choosen for nodes/tips choosen
# DEPENDENCIES:
# simulate
# PURPOSE:
# Extract the values wanted from the raw result of function simulate.
# NOTES:
# none
# REVISIONS:
# 16/05/14 - Initial release
# 28/05/14 - Add optimal.values
# 02/06/14 - Case where optimal.value is asked for a BM (retrun NULL)
##
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