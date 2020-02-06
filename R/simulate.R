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
#' @param x an object of class \code{\link{params_process}} or \code{\link{PhyloEM}}.
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param checks whether to check the entry parameters for consistency. Default 
#' to TRUE.
#' @param simulate_random set to FALSE if only the expected values are needed
#' (and not the random sample). Default to TRUE.
#' @param U_tree optional, full incidence matrix of the tree, result of function
#' \code{\link{incidence.matrix.full}}. Can be specified to avoid extra computations.
#' @param times_shared optional, times of shared ancestry of all nodes and tips,
#' result of function \code{\link{compute_times_ca}}. Can be specified to avoid extra
#' computations.
# @param nsim Unused.
# @param seed Unused.
#' @param ... for a \code{PhyloEM} object, further arguments to be passed on to
#' \code{\link{params_process.PhyloEM}} (to choose which parameters to extract from
#' the results, see documentation of this function).
#'     
#' @return An S3 object of class \code{simul_process}. This contains:
#' \describe{
#'  \item{sim_traits}{an array with dimensions p x Nnode x 2 (BM)
#'  or p x Nnode x 3 (OU). For each trait t, 1 <= t <= p, sim_traits[t, , ] has
#'  tree columns, containing respectively the simulated state,
#'  expected value and optimal value for all the nodes.}
#'  \item{phylo}{the phylogenetic tree used for the simulations (class \code{phylo}).}
#'  \item{params}{the parameters used for the simulations
#'  (class \code{params_proces}).}
#'  }
#'  
#' @seealso \code{\link{params_process}}, \code{\link{PhyloEM}}
#' 
#' @export
#' 
##
simul_process <- function(x, ...) UseMethod("simul_process")

##
#' @describeIn simul_process \code{\link{params_process}} object
#' @export
##
simul_process.params_process <- function(x, 
                                         phylo, simulate_random = TRUE,
                                         checks = TRUE,
                                         U_tree = NULL,
                                         times_shared = NULL, ...){
  
  if (x$process == "BM"){ ## Just to be safe
    x$selection.strength <- NULL
    x$optimal.value <- NULL
  }
  
  sim <- simulate_internal(phylo,
                           process = x$process,
                           p = ncol(x$variance),
                           root.state = x$root.state,
                           shifts = x$shifts,
                           eps = 10^(-6),
                           selection.strength = x$selection.strength,
                           variance = x$variance,
                           optimal.value = x$optimal.value,
                           checks = checks,
                           simulate_random = simulate_random,
                           U_tree = U_tree,
                           times_shared = times_shared,
                           df = x$df)
  
  res <- list(sim_traits = sim,
              phylo = phylo,
              params = x)
  class(res) <- "simul_process"
  return(res)
}

##
#' @describeIn simul_process \code{\link{PhyloEM}} object
#' @export
##
simul_process.PhyloEM <- function(x, 
                                  simulate_random = TRUE,
                                  checks = TRUE,
                                  U_tree = NULL,
                                  times_shared = NULL, ...){
  
  params <- params_process(x, ...)
  
  return(simul_process.params_process(params, 
                                      x$phylo,
                                      simulate_random = simulate_random,
                                      checks = checks,
                                      U_tree = U_tree,
                                      times_shared = times_shared))
}

##
#' ##
#' @title Extraction of simulated traits
#'
#' @description
#' \code{extract.simul_process} takes an object of class "\code{simul_process}",
#' result of function \code{\link{simul_process}}, and extracts the traits values, 
#' expectations or optimal values at the tips or the internal nodes.
#'
#' @param x an object of class "\code{simul_process}", result of function
#' \code{\link{simul_process}}.
#' @param where one of "tips" (the default) or "nodes". Where to extract the results.
#' @param what one of "states" (the default), "expectation", or "optimal.values".
#' @param ... unused
#' 
#' @return A matrix giving the selected quantities at the selected nodes or tips. If
#' the tips or nods are labeled, then the colnames of the matrix are set accordingly.
#' 
#' @seealso \code{\link{simul_process}}
#' 
#' @export
##
extract.simul_process <- function(x,
                                  where = c("tips", "nodes"),
                                  what = c("states", "expectations", "optimal.values"),
                                  ...) {
  res <- extract_simulate_internal(x$sim_traits, where = where, what = what)
  if (where == "tips"){
    colnames(res) <- x$phylo$tip.label
  } else if (where == "nodes"){
    colnames(res) <- x$phylo$node.label
  }
  return(res)
}

##
#' @title Plot for class \code{simul_process}
#'
#' @description
#' This function takes an object of class \code{params_process}, and plots them along
#' with some data at the tips of the tree.
#'
#' @param x an object of class \code{params_process}.
#' @param phylo a phylogenetic tree.
#' @param data a matrix of data at the tips of the tree. Must have p rows and
#' ntaxa columns. If these are simulated, use the \code{\link{extract.simul_process}}
#' function.
#' @param ancestral_states if \code{plot_ancestral_states=TRUE}, the ancestral states 
#' must be specified. If these are simulated, use the
#' \code{\link{extract.simul_process}} function.
#' @inheritParams plot.PhyloEM
#' 
#' @return
#' NULL
#' 
#' @seealso \code{\link{simul_process}}, \code{\link{plot.PhyloEM}},
#' \code{\link{params_BM}}, \code{\link{params_OU}}
#' 
#' @export
#'

plot.params_process <- function(x,
                                phylo,
                                data = NULL,
                                traits,
                                automatic_colors = TRUE,
                                color_characters = "black",
                                color_edges = "black",
                                plot_ancestral_states = FALSE,
                                ancestral_states = NULL,
                                imposed_scale,
                                ancestral_cex = 2,
                                ancestral_pch = 19,
                                value_in_box = FALSE,
                                ancestral_as_shift = FALSE,
                                shifts_cex = 0.6,
                                shifts_bg = "chocolate4",
                                root_bg = "chocolate4",
                                shifts_adj = 0,
                                root_adj = 1,
                                color_shifts_regimes = FALSE,
                                regime_boxes = FALSE,
                                alpha_border = 70,
                                show.tip.label = FALSE,
                                label_cex = 0.5,
                                label_offset = 0,
                                axis_cex = 0.7,
                                edge.width = 1,
                                margin_plot = NULL,
                                gray_scale = FALSE,
                                ...){
  if (missing(traits)) traits <- 1:ncol(x$variance)
  ## Checking consistency
  if (plot_ancestral_states && length(traits) > 1) stop("Ancestral state plotting is only allowed for one single trait. Please select the trait you would like to plot with argument 'traits' (see documentation).")
  if (plot_ancestral_states && is.null(ancestral_states)){
    stop("The ancestral traits must be specified. Use function 'extract' if you simulated them.")
  }
  if (value_in_box && length(traits) > 1){
    stop("Showing the shifts values on the tree is only allowed for one single trait. Please select the trait you would like to plot with argument 'traits' (see documentation).")
  }
  if (!is.null(data)){
    if (nrow(data) < length(traits)){
      stop("The data matrix must have as many rows as the number of traits to be plotted.")
    }
    if (ncol(data) != length(phylo$tip.label)){
      stop("The data matrix must have as many columns as the number of tips in the tree.")
    } 
  }
  
  # ## Save curent par
  # .pardefault <- par(no.readonly = T)
  # on.exit(par(.pardefault), add = TRUE)
  
  ## parameters
  params <- x
  
  # If on trait, select relevant quantities
  if (length(traits) == 1){
    if (length(as.vector(params$selection.strength)) == 1) params$selection.strength <- diag(rep(params$selection.strength, ncol(x$variance)))
    params <- split_params_independent(params)
    params <- params[[traits]]
  }
  
  Y_state <- data
  if (missing(imposed_scale)) imposed_scale <- Y_state
  
  ## Plotting
  plot.data.process.actual(Y.state = Y_state[traits, , drop = FALSE],
                           phylo = phylo,
                           params = params,
                           process = x$process,
                           miss = is.na(Y_state[traits, , drop = FALSE]),
                           imposed_scale = imposed_scale,
                           root_adj = root_adj,
                           shifts_adj = shifts_adj,
                           shifts_bg = shifts_bg,
                           root_bg = root_bg,
                           quant.root = 0.25,
                           color_characters = color_characters,
                           color_edges = color_edges,
                           edge.width = edge.width,
                           automatic_colors = automatic_colors,
                           regime_boxes = regime_boxes,
                           alpha_border = alpha_border,
                           value_in_box = value_in_box,
                           shifts_cex = shifts_cex,
                           axis_cex = axis_cex,
                           margin_plot = margin_plot,
                           color_shifts_regimes = color_shifts_regimes,
                           # shifts_regimes = shifts_regimes,
                           plot_ancestral_states = plot_ancestral_states,
                           ancestral_states = ancestral_states,
                           # imposed_scale.nodes = imposed_scale.nodes,
                           ancestral_cex = ancestral_cex,
                           ancestral_pch = ancestral_pch,
                           label_cex = label_cex,
                           show.tip.label = show.tip.label,
                           # underscore = underscore,
                           # label.offset = label.offset,
                           ancestral_as_shift = ancestral_as_shift,
                           gray_scale = gray_scale,
                           ...)
}

##
#' @title Simulate a Stochastic Process on a tree
#'
#' @description
#' \code{simulate_internal} simulate a stochastic process on a tree.
#'
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param process The model used for the simulation. One of "BM" (for a full BM
#' model, univariate or multivariate); "OU" (for a full OU model, univariate or
#' multivariate); or "scOU" (for a "scalar OU" model).
#' @param p Dimension of the simulated trait
#' @param root.state List describing the state of the root, with:
#' \describe{
#'  \item{random}{random state (TRUE) or deterministic state (FALSE)}
#'  \item{value.root}{if deterministic, value of the character at the root}
#'  \item{exp.root}{if random, expectation of the character at the root}
#'  \item{var.root}{if random, variance of the character at the root (pxp matrix)}
#'  }
#' @param shifts List with position and values of the shifts :
#' \describe{
#'  \item{edges}{vector of the K id of edges where the shifts are}
#'  \item{values}{matrix p x K of values of the shifts on the edges
#'   (one column = one shift)}
#'  \item{relativeTimes}{vector of dimension K of relative time of the shift from the parent node of edges}
#'  }
#' @param eps Tolerance for the value of the norm 1 of the selection strength matrix for OU
#' @param variance Variance-covariance matrix size p x p 
#' @param selection.strength Matrix of selection strength size p x p (OU)
#' @param optimal.value Vector of p optimal values at the root (OU)
#' @param checks whether to check the entry parameters for consistency. Default 
#' to TRUE.
#' @param simulate_random set to FALSE if only the expected values are needed
#' (and not the random sample). Default to TRUE.
#' @param df if the process is "StudentOU", the number of degree of freedom of
#' the chosen student law. default to 1.
#' @param U_tree optional, full incidence matrix of the tree, result of function
#' \code{incidence.matrix.full}.
#' @param times_shared optional, times of shared ancestry of all nodes and tips,
#' result of function \code{\link{compute_times_ca}}. Can be specified to avoid extra
#' computations.
#'     
#' @return paramSimu An array with dimensions p x Nnode x 2 (BM)
#'  or p x Nnode x 3 (OU). For each trait t, 1 <= t <= p, paramSimu[t, , ] has
#'  tree columns, containing respectively the simulated state,
#'  expected value and optimal value for all the nodes.
#'  
#' @keywords internal
#' 
# 16/05/14 - Initial release
# 20/05/14 - Gestion of edges (function correspondanceEdges)
# 21/05/14 - extraction of the recursion
# 16/06/14 - check.selection.strength
# 24/08/15 - Multivariate
##
simulate_internal <- function(phylo,
                              process = c("BM", "OU", "scOU", "OUBM", "StudentOU"),
                              p = 1,
                              # independent = FALSE,
                              root.state = list(random = FALSE, 
                                                stationary.root = FALSE, 
                                                value.root = NA, 
                                                exp.root = NA,
                                                var.root = NA),
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
                              times_shared = NULL,
                              df = 1) {
  # library(MASS)
  ntaxa <- length(phylo$tip.label)
  ## Set branch stochastic process
  process <- match.arg(process)
  if (process == "OUBM") {
    if (!isDiagonal(selection.strength)){
      stop("The OUBM simulation is only implemented for a diagonal matrix.")
    }
    if (!is.null(shifts$edges)) {
      stop("The OUBM simulation is only implemented for a process with no shift.")
    }
  }
  if (process == "scOU"){
     # Matrix provided
    if (!is.null(dim(selection.strength))){
      zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
        if (length(x) == 1) return(TRUE)
        x <- range(x) / mean(x)
        isTRUE(all.equal(x[1], x[2], tolerance = tol))
      }
      if (!isDiagonal(selection.strength)){
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
  if (is.null(optimal.value) && process %in% c("OU", "scOU", "OUBM")){
    stop("Optimal values for the OU simulation must be specified.")
  }
  ## Special case BM
  if (!simulate_random){
    if (process == "BM") {
      paramSimu <- compute_expectations.BM(phylo,
                                           root.state = root.state,
                                           shifts = shifts,
                                           U_tree = U_tree)
      attr(paramSimu, "ntaxa") <- ntaxa
      return(paramSimu) 
    } else if (process == "scOU"){
      paramSimu <- compute_expectations.scOU(phylo,
                                             root.state = root.state,
                                             shifts = shifts,
                                             alpha = selection.strength,
                                             U_tree = U_tree,
                                             times_shared = times_shared)
      attr(paramSimu, "ntaxa") <- ntaxa
      return(paramSimu) 
    }
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
                 StudentOU = init.simulate.OU,
                 OUBM = init.simulate.OU)
  updateDown <- switch(process,
                       BM = update.simulate.BM,
                       OU = update.simulate.OU,
                       scOU = update.simulate.scOU,
                       StudentOU = update.simulate.StudentOU,
                       OUBM = update.simulate.OUBM)
  ## Check root
  if (checks) {
    root.state <- test.root.state(root.state = root.state,
                                  process = process,
                                  eps = eps,
                                  variance = variance,
                                  selection.strength = selection.strength,
                                  optimal.value = optimal.value)
  }
  ## Initialization and setting root state
  paramSimu <- init(phy = phy,
                    p = p,
                    root.state = root.state,
                    optimal.value = optimal.value,
                    simulate_random = simulate_random)
  if (process %in% c("scOU", "OU", "OUBM")){
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
# extract.simulate (paramSimu, where=c("tips","nodes"), what=c("states","expectations"))
# PARAMETERS:
# @paramSimu (matrix) return of the function simulate
# @where (string) : where to extract the values : at the "tips" or the internal "nodes" ?
# @what (string) : which value to extract : the simulated "states" or the "expectations" ?
# RETURNS:
# (vector) values chosen for nodes/tips chosen
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
extract_simulate_internal <- function(paramSimu,
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
#' @return array: array p x Nnode x 2 (BM), with slice corresponding to node filled with value
#'  
#' @keywords internal
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
#' @description Function used in \code{\link{simulate}} for BM/OU initializations.
#'
#' @param phy: Input tree.
#' @param p: dimension of the trait simulated
#' @param root.state (list): state of the root, with:
#'     random : random state (TRUE) or deterministic state (FALSE)
#'     value.root : if deterministic, value of the character at the root
#'     exp.root : if random, expectation of the character at the root
#'     var.root : if random, variance of the character at the root (pxp matrix)
#'     
#' @return paramSimu: array p x Nnode x 2 (BM), filled with NAs.
#' Slice paramSimu[, ntaxa + 1, ] (array p x 2) is initialized with simulated states and root
#' expectations for all the traits.
#' 
#' @keywords internal
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
#' @description Function used in \code{\link{simulate}} for BM initialization.
#'
#' @param phy Input tree.
#' @param root.state (list) State of the root, with:
#'     random : random state (TRUE) or deterministic state (FALSE)
#'     value.root : if deterministic, value of the character at the root
#'     exp.root : if random, expectation of the character at the root
#'     var.root : if random, variance of the character at the root (pxp matrix)
#'     
#' @return paramSimu Array p x Nnode x 2 (BM), filled with NAs.
#' Slice paramSimu[, ntaxa + 1, ] (array p x 2) is initialized with simulated
#' states and root expectations for all the traits.
#' 
#' @keywords internal
#'  
##

init.simulate.BM <- function(phy, p, root.state, simulate_random, ...){
  return(init.simulate.StateAndExp(phy, p, root.state, simulate_random))
}

##
#' @title Initialize state and expectation matrices
#' 
#' @description Function used in \code{\link{simulate}} for OU initialization.
#'
#' @param phy: Input tree.
#' @param p: dimension of the trait simulated
#' @param root.state (list): state of the root, with:
#'     random : random state (TRUE) or deterministic state (FALSE)
#'     value.root : if deterministic, value of the character at the root
#'     exp.root : if random, expectation of the character at the root
#'     var.root : if random, variance of the character at the root (pxp matrix)
#'     
#' @return paramSimu: array p x Nnode x 3, filled with NAs.
#' Slice paramSimu[, ntaxa + 1, ] (array p x 3) is initialized with simulated states, root
#' expectations, and optimal values for all the traits.
#'  
#' @keywords internal
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
# @eps (double) : if the selection strength is smaller than eps, simulate according to a BM instead of an OU
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

update.simulate.OUBM <- function(edgeNbr, ancestral,
                                 length, shifts, selection.strength,
                                 stationary_variance,
                                 variance,
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
    Sigma_sim <- as.matrix(stationary_variance - ee %*% stationary_variance %*% t(ee))
    Sigma_sim[is.na(Sigma_sim)] <- length * variance[is.na(Sigma_sim)]
    Sim <- MASS::mvrnorm(1,
                         mu = ee %*% ancestral[ , , 1] + plus_exp,
                         Sigma = Sigma_sim)
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
#'     parent node of edges
#'     
#' @return paramSimu: array p x Nnode x 2 (BM). For each trait t, 1 <= t <= p,
#'  paramSimu[t, , ] has two columns, both containing the expected values for
#'  all the nodes.
#'  
#' @keywords internal
#'  
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

##
#' @title Compute the expected states of a scOU
#'
#' @description
#' \code{compute_expectations.scOU} use the matrix formulation to compute the 
#' expected values at all the nodes. Assumes a stationary root.
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
#'     parent node of edges
#'     
#' @return paramSimu: array p x Nnode x 2 (BM). For each trait t, 1 <= t <= p,
#'  paramSimu[t, , ] has two columns, both containing the expected values for
#'  all the nodes.
#'  
#' @keywords internal
#'  
##
compute_expectations.scOU <- function(phylo, root.state, shifts, alpha, 
                                      U_tree = NULL, times_shared = NULL){
  if (is.null(U_tree)) U_tree <- incidence.matrix.full(phylo)
  if (is.null(times_shared)) times_shared <- compute_times_ca(phylo)
  A <- diag(exp(-alpha * diag(times_shared)))
  B <- diag(exp(alpha * diag(times_shared[phylo$edge[, 1], phylo$edge[, 1]])))
  Delta <- shifts.list_to_matrix(phylo, shifts)
  if (root.state$random){
    root_value = root.state$exp.root
  } else {
    root_value = root.state$value.root
  }
  paramSimu <- tcrossprod(Delta, U_tree - A%*%U_tree%*%B) + root_value
  # Put in the right format (for extraction)
  paramSimu <- array(rep(paramSimu, 2), c(dim(paramSimu), 3))
  ## Computes optimal values
  opt_val <- tcrossprod(Delta, U_tree) + root_value
  paramSimu[, , 3] <- opt_val
  return(paramSimu)
}