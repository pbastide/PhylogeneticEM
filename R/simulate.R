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
#' @title Simulate a Stochastic Process on a tree
#'
#' @description
#' \code{simulate} simulate a stochastic process on a tree.
#'
#' @param phylo: Input tree.
#' @param process: choose BM or OU process
#' @param p: dimention of the trait simulated
#' @param root.state (list): state of the root, with:
#'     random : random state (TRUE) or deterministic state (FALSE)
#'     value.root : if deterministic, value of the character at the root
#'     exp.root : if random, expectation of the character at the root
#'     var.root : if random, variance of the character at the root (pxp matrix)
#' @param shifts (list) position and values of the shifts :
#'     edges : vector of the K id of edges where the shifts are
#'     values : matrix p x K of values of the shifts on the edges (one column = one shift)
#'     relativeTimes : vector of dimension K of relative time of the shift from the
#'     parentnode of edges
#' @param eps: tolerance for the value of the norm 1 of the selection strength matrix for OU
#' @param variance: variance-covariance matrix size p x p 
#' @param selection.strenght: matrix of selection strength size p x p (OU)
#' @param optimal.value: vector of p optimal values at the root (OU)
#'     
#' @return paramSimu: array p x nNodes x 2 (BM) or p x nNodes x 3 (OU). For each trait t, 
#' 1 <= t <= p, paramSimu[t, , ] has tree columns, containing respectively the simulated state,
#'  expected value and optimal value for all the nodes.
#' 
#' 16/05/14 - Initial release
#' 20/05/14 - Gestion of edges (function correspondanceEdges)
#' 21/05/14 - extraction of the recursion
#' 16/06/14 - check.selection.strength
#' 24/08/15 - Multivariate
##

simulate <- function(phylo,
                     process = c("BM", "OU", "scOU", "StudentOU"),
                     p = 1,
                     # independent = FALSE,
                     root.state = list(random = FALSE, 
                                       stationary.root = FALSE, 
                                       value.root, 
                                       exp.root,
                                       var.root),
                     shifts = list(edges = NULL,
                                   values = NULL,
                                   relativeTimes = NULL),
                     eps=10^(-6),
                     selection.strength = NULL,
                     variance = NULL,
                     optimal.value = NULL,
                     checks = TRUE,
                     simulate_random = TRUE,
                     U_tree = NULL,
                     df = 1) {
  # library(MASS)
  ntaxa <- length(phylo$tip.label)
  ## Set branch stochastic process
  process <- match.arg(process)
  if (process == "scOU"){
    # Use a OU (more efficient things can be done)
    # process <- "OU"
    # Check selection strength
#      # scalar
#     if (is.null(dim(selection.strength))){
#       selection.strength <- selection.strength *  diag(rep(1, p))
#     }
     # Matrix provided
    if (!is.null(dim(selection.strength))){
      zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
        if (length(x) == 1) return(TRUE)
        x <- range(x) / mean(x)
        isTRUE(all.equal(x[1], x[2], tolerance = tol))
      }
      if (!all(selection.strength[!diag(nrow(selection.strength))] == 0)){
        stop("Process is said to be scalar OU, but selection strengh matrix is not diagonal.")
      }
      if (!zero_range(diag(selection.strength))){
        stop("Process is said to be scalar OU, but selection strengh matrix is diagonal, but not scalar.")
      }
      selection.strength <- selection.strength[1, 1]
    }
    ## Check Dimensions
    if (checks){
      parameters <- check_dimensions(p = p, root.state = root.state,
                                     shifts = shifts, variance = variance,
                                     selection.strength = NULL, 
                                     optimal.value = optimal.value)
      root.state <- parameters$root.state
      shifts <- parameters$shifts
      variance <- as(parameters$variance, "dpoMatrix")
      # selection.strength <- parameters$selection.strength
      optimal.value <- parameters$optimal.value
    }
  } else {
    ## Check Dimensions
    if (checks){
      parameters <- check_dimensions(p = p, root.state = root.state,
                                     shifts = shifts, variance = variance,
                                     selection.strength = selection.strength, 
                                     optimal.value = optimal.value)
      root.state <- parameters$root.state
      shifts <- parameters$shifts
      variance <- as(parameters$variance, "dpoMatrix")
      selection.strength <- parameters$selection.strength
      optimal.value <- parameters$optimal.value
    }
  }
  ## Optimal values
  if (is.null(optimal.value) && process %in% c("OU", "scOU")){
    stop("Optimal values for the OU simulation must be specified.")
  }
  ## Special case BM
  if ((process == "BM") && !simulate_random){
    paramSimu <- compute_expectations.BM(phylo,
                                         root.state = root.state,
                                         shifts = shifts,
                                         U_tree = U_tree)
    attr(paramSimu, "ntaxa") <- ntaxa
    return(paramSimu)
  }
  ## Reorder tree
  phy <- reorder(phylo, order = "cladewise")
  # Trace edges
  shifts_ordered <- shifts
  shifts_ordered$edges <- correspondanceEdges(edges = shifts$edges,
                                              from = phylo, to = phy)
  ## Set branch stochastic process
  process <- match.arg(process)
  process <- check.selection.strength(process = process,
                                      selection.strength = selection.strength,
                                      eps = eps) # if OU, check if selection.strength is not too low.
  init <- switch(process,
                 BM = init.simulate.BM,
                 OU = init.simulate.OU,
                 scOU = init.simulate.OU,
                 StudentOU = init.simulate.OU)
  updateDown <- switch(process,
                       BM = update.simulate.BM,
                       OU = update.simulate.OU,
                       scOU = update.simulate.scOU,
                       StudentOU = update.simulate.StudentOU)
  ## Check root
  if (checks) {
    root.state <- test.root.state(root.state = root.state,
                                  process = process,
                                  eps = eps,
                                  variance = variance,
                                  selection.strength = selection.strength,
                                  optimal.value = optimal.value)
  }
  ## Initialisation and setting root state
  paramSimu <- init(phy = phy,
                    p = p,
                    root.state = root.state,
                    optimal.value = optimal.value,
                    simulate_random = simulate_random)
  if (process %in% c("scOU", "OU")){
    if (root.state$stationary.root){
      stationary_variance <- root.state$var.root
    } else {
      stationary_variance <- compute_stationary_variance(variance, selection.strength)
    }
  } else {
    stationary_variance <- NA
  }
  ## Tree recursion
  paramSimu <- recursionDown(phy = phy,
                             params = paramSimu,
                             updateDown = updateDown,
                             subset_node = subset_node.simulate,
                             allocate_subset_node = allocate_subset_node.simulate,
                             shifts = shifts_ordered,
                             variance = variance,
                             eps = eps,
                             selection.strength = selection.strength,
                             stationary_variance = stationary_variance,
                             simulate_random = simulate_random,
                             df = df)
  attr(paramSimu, "ntaxa") <- ntaxa
  return(paramSimu)
}

##
#' @title Iteration allocation
#'
#' @description
#' \code{allocate_subset_node.simulate} slice the data correctly and allocate the right part.
#' To be used in function UpdateDown.
#'
#' @param node: node on which to slice
#' @param array: structure to be sliced
#' @param value: value to be attributed to the slice
#'     
#' @return array: array p x nNodes x 2 (BM), with slice corresponding to node filled with value
#'  
##

allocate_subset_node.simulate <- function(node, array, value){
  array[, node, ] <- value
  return(array)
}

subset_node.simulate <- function(node, array){
  return(array[, node, , drop = F])
}


##
#' @title Initialize state and expectation matrices
#'
#' @param phy: Input tree.
#' @param p: dimension of the trait simulated
#' @param root.state (list): state of the root, with:
#'     random : random state (TRUE) or deterministic state (FALSE)
#'     value.root : if deterministic, value of the character at the root
#'     exp.root : if random, expectation of the character at the root
#'     var.root : if random, variance of the character at the root (pxp matrix)
#'     
#' @return paramSimu: array p x nNodes x 2 (BM), filled with NAs.
#' Slice paramSimu[, ntaxa + 1, ] (array p x 2) is initialized with simulated states and root
#' expectarions for all the traits.
#'  
##

init.simulate.StateAndExp <- function(phy, p, root.state, simulate_random){
  ntaxa <- length(phy$tip.label)
  paramSimu <- array(NA, dim = c(p, 1 + nrow(phy$edge), 2))
  if (!root.state$random) { # The root is not random
    paramSimu <- allocate_subset_node.simulate(ntaxa + 1, paramSimu,
                                               cbind(root.state$value.root,
                                                     root.state$value.root))
  } else { # The value of the root is random N(exp.root, var.root)
    if (simulate_random){
      sim_rand <- MASS::mvrnorm(1, mu = root.state$exp.root,
                                Sigma = root.state$var.root)
    } else {
      sim_rand <- root.state$exp.root
    }
    paramSimu <- allocate_subset_node.simulate(ntaxa + 1, paramSimu,
                                               cbind(sim_rand,
                                                     root.state$exp.root))
  }
  return(paramSimu)
}


##
#' @title Initialize BM
#'
#' @param phy: Input tree.
#' @param root.state (list): state of the root, with:
#'     random : random state (TRUE) or deterministic state (FALSE)
#'     value.root : if deterministic, value of the character at the root
#'     exp.root : if random, expectation of the character at the root
#'     var.root : if random, variance of the character at the root (pxp matrix)
#'     
#' @return paramSimu: array p x nNodes x 2 (BM), filled with NAs.
#' Slice paramSimu[, ntaxa + 1, ] (array p x 2) is initialized with simulated states and root
#' expectarions for all the traits.
#'  
##

init.simulate.BM <- function(phy, p, root.state, simulate_random, ...){
  return(init.simulate.StateAndExp(phy, p, root.state, simulate_random))
}

##
#' @title Initialize state and expectation matrices
#'
#' @param phy: Input tree.
#' @param p: dimension of the trait simulated
#' @param root.state (list): state of the root, with:
#'     random : random state (TRUE) or deterministic state (FALSE)
#'     value.root : if deterministic, value of the character at the root
#'     exp.root : if random, expectation of the character at the root
#'     var.root : if random, variance of the character at the root (pxp matrix)
#'     
#' @return paramSimu: array p x nNodes x 3, filled with NAs.
#' Slice paramSimu[, ntaxa + 1, ] (array p x 3) is initialized with simulated states, root
#' expectations, and optimal values for all the traits.
#'  
##

init.simulate.OU <- function(phy, p, root.state, optimal.value,
                             simulate_random, ...){
  paramSimu <- array(NA, dim = c(p, 1 + nrow(phy$edge), 3))
  paramSimu[, , c(1,2)] <- init.simulate.StateAndExp(phy, p, root.state,
                                                     simulate_random)
  ntaxa <- length(phy$tip.label)
  beta <- array(NA, dim = c(p, 1 + nrow(phy$edge))) # selection strength
  beta[ , ntaxa + 1] <- optimal.value
  paramSimu[, , 3] <- beta
  return(paramSimu)
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
update.simulate.BM <- function(edgeNbr, ancestral, length, shifts, variance,
                               simulate_random, ...){
  shiftsIndex <- which(shifts$edges == edgeNbr) # If no shifts = NULL, and sum = 0
  shiftsValues <- rowSums(shifts$values[, shiftsIndex, drop = F])
  if (simulate_random){
    sim_value = MASS::mvrnorm(1, mu = shiftsValues, Sigma = length*variance)
  } else {
    sim_value = shiftsValues
  }
  return(cbind(ancestral[, , 1, drop = F] + sim_value,
               ancestral[, , 2, drop = F] + shiftsValues))
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
update.simulate.OU <- function(edgeNbr, ancestral,
                               length, shifts, selection.strength,
                               stationary_variance,
                               simulate_random, ...){
  shiftsIndex <- which(shifts$edges == edgeNbr) # If no shifts = NULL, and sum = 0
  if (length(shiftsIndex) == 0){
    r <- 0
  } else {
    r <- shifts$relativeTimes[shiftsIndex]
  }
  beta <- ancestral[, , 3] + rowSums(shifts$values[, shiftsIndex, drop = F])
  if (r == 0){
    ee <- expm(-selection.strength * length)
    I <- diag(1, dim(selection.strength))
    plus_exp <- (I - ee) %*% beta
  } else {
    ee_d <- expm(-selection.strength * length * (1-r))
    ee_p <- expm(-selection.strength * length * r)
    ee <- expm(-selection.strength * length)
    I <- diag(1, dim(selection.strength))
    plus_exp <- (I - ee_d) %*% beta + (I - ee_p) %*% ancestral[ , , 3]
  }
  Exp <- ee %*% ancestral[ , , 2] + plus_exp
  Exp <- as.matrix(Exp)
  if (simulate_random){
    Sim <- MASS::mvrnorm(1,
                         mu = ee %*% ancestral[ , , 1] + plus_exp,
                         Sigma = stationary_variance - ee %*% stationary_variance %*% t(ee))
  } else {
    Sim <- Exp
  }
  child <- ancestral
  p <- dim(ancestral)[1]
  child[, , 1] <- array(Sim, dim = c(p, 1, 1))
  child[, , 2] <- array(Exp, dim = c(p, 1, 1))
  child[, , 3] <- array(beta, dim = c(p, 1, 1))
  return(child)
}

update.simulate.scOU <- function(edgeNbr, ancestral,
                               length, shifts, selection.strength,
                               stationary_variance,
                               simulate_random, ...){
  shiftsIndex <- which(shifts$edges == edgeNbr) # If no shifts = NULL, and sum = 0
  if (length(shiftsIndex) == 0){
    r <- 0
  } else {
    r <- shifts$relativeTimes[shiftsIndex]
  }
  beta <- ancestral[, , 3] + rowSums(shifts$values[, shiftsIndex, drop = F])
  if (r == 0){
    ee <- exp(-selection.strength * length)
    plus_exp <- (1 - ee) * beta
  } else {
    ee_d <- exp(-selection.strength * length * (1-r))
    ee_p <- exp(-selection.strength * length * r)
    ee <- exp(-selection.strength * length)
    plus_exp <- (1 - ee_d) * beta + (1 - ee_p) * ancestral[ , , 3]
  }
  Exp <- ee * ancestral[ , , 2] + plus_exp
  Exp <- as.matrix(Exp)
  if (simulate_random){
    Sim <- MASS::mvrnorm(1,
                         mu = ee * ancestral[ , , 1] + plus_exp,
                         Sigma = (1 - ee^2) * stationary_variance)
  } else {
    Sim <- Exp
  }
  child <- ancestral
  p <- dim(ancestral)[1]
  child[, , 1] <- array(Sim, dim = c(p, 1, 1))
  child[, , 2] <- array(Exp, dim = c(p, 1, 1))
  child[, , 3] <- array(beta, dim = c(p, 1, 1))
  return(child)
}

update.simulate.StudentOU <- function(edgeNbr, ancestral, length, shifts, variance, selection.strength, df, ...){
  shiftsIndex <- which(shifts$edges==edgeNbr) # If no shifts = NULL, and sum = 0
  beta <- ancestral[3] + sum(shifts$values[shiftsIndex])
  ee <- exp(-selection.strength*length)
  ss <- sum(shifts$values[shiftsIndex]*( 1-exp( -selection.strength * length * (1-shifts$relativeTimes[shiftsIndex]) ) ))
  SimExp <- c( ancestral[3]*(1-ee) + ancestral[1]*ee + ss + sqrt(variance*(1-ee^2)/(2*selection.strength))*rt(1, df),
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
extract.simulate <- function(paramSimu,
                             where=c("tips", "nodes"),
                             what=c("states", "expectations", "optimal.values")){
  where <- match.arg(where)
  what <- match.arg(what)
  ntaxa <- attr(paramSimu,"ntaxa")
  if (where=="tips") {
    rows <- 1:ntaxa
  } else if (where=="nodes") {
    rows <- (ntaxa+1):dim(paramSimu)[2]
  }
  if (what=="states") {
    col <- 1
  } else if (what=="expectations") {
    col <- 2
  } else if (what=="optimal.values") {
    col <- 3
  }
  if (col > dim(paramSimu)[3] ){
    return(NULL) # Case of optimal.value asked for a BM
  } else {
    return(matrix(paramSimu[, rows, col], nrow = dim(paramSimu)[1]))
  }
}

##
#' @title Compute the expected states of a BM
#'
#' @description
#' \code{compute_expectations.BM} use the matrix formulation to compute the 
#' expected values at all the nodes.
#'
#' @param phylo: Input tree.
#' @param root.state (list): state of the root, with:
#'     random : random state (TRUE) or deterministic state (FALSE)
#'     value.root : if deterministic, value of the character at the root
#'     exp.root : if random, expectation of the character at the root
#'     var.root : if random, variance of the character at the root (pxp matrix)
#' @param shifts (list) position and values of the shifts :
#'     edges : vector of the K id of edges where the shifts are
#'     values : matrix p x K of values of the shifts on the edges (one column = one shift)
#'     relativeTimes : vector of dimension K of relative time of the shift from the
#'     parentnode of edges
#'     
#' @return paramSimu: array p x nNodes x 2 (BM). For each trait t, 1 <= t <= p,
#'  paramSimu[t, , ] has two columns, both containing the expected values for
#'  all the nodes.
##
compute_expectations.BM <- function(phylo, root.state, shifts, U_tree = NULL){
  if (is.null(U_tree)) U_tree <- incidence.matrix.full(phylo)
  Delta <- shifts.list_to_matrix(phylo, shifts)
  if (root.state$random){
    root_value = root.state$exp.root
  } else {
    root_value = root.state$value.root
  }
  paramSimu <- tcrossprod(Delta, U_tree) + root_value
  # Put in the right format (for extraction)
  paramSimu <- array(rep(paramSimu, 2), c(dim(paramSimu), 2))
  return(paramSimu)
}