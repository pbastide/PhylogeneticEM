# {Estimate EM}
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
## Here is the complete implementation of the EM.
## Dependencies : generic_functions.R
## : simulate.R
## : shifts_manipulations.R
## : init_EM.R
## : E_step.R
## : M_step.R
## : shutoff.R
###############################################################################

# @import stats
# @import utils
# @import graphics
#' @import methods
#' @import ape
#' @import Matrix
NULL

#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @importFrom stats coef cophenetic cov lm lm.fit mad median na.omit nls.control
#' optimize quantile rbinom reorder residuals rnorm rt var
#' @importFrom graphics axis close.screen legend par plot rect screen segments
#' split.screen strwidth text
#' @importFrom grDevices col2rgb gray gray.colors palette rainbow rgb
#' @importFrom utils capture.output combn setTxtProgressBar txtProgressBar
NULL
###############################################################################
## estimateEM
###############################################################################

##
# estimateEM (phylo, Y_data, process=c("BM","OU"), tol=10^(-5),  method.variance=c("simple"), method.init=c("default"), nbr_of_shifts=0, ...)
# PARAMETERS:
#            @phylo (tree) input tree
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
#            init.EM, compute_E, compute_M, compute_times_ca, compute_dist_phy, shutoff.EM, update.parsimonyNumber (, ...)
# PURPOSE:
#            Run the EM algorithm
# NOTES:
#            Under construction
# REVISIONS:
#            21/05/14 - Initial release (non fonctionnal)
#            22/05/14 - Minimal "working" release
#            02/06/14 - OU in the special case
#            10/06/14 - Test of divergence
##
##
#' @title Perform One EM
#'
#' @description
#' \code{EstimateEM} performs one EM for one given number of shifts. It is called
#' from function \code{\link{PhyloEM}}. Its use is mostly internal, and most user
#' should not need it.
#'
#' @details
#' See documentation of \code{\link{PhyloEM}} for further details.
#' All the parameters monitoring the EM (like \code{tol_EM}, \code{Nbr_It_Max}, etc.)
#' can be called from \code{PhyloEM}.
#' 
#' @inheritParams PhyloEM
#' @param Y_data_imp (optional) imputed data if previously computed, same format as
#' \code{Y_data}. Mostly here for internal calls.
#' @param tol_EM the tolerance for the convergence of the parameters. A named list, with
#' items:
#' \describe{
#' \item{variance}{default to 10^(-2)}
#' \item{value.root}{default to 10^(-2)}
#' \item{exp.root}{default to 10^(-2)}
#' \item{var.root}{default to 10^(-2)}
#' \item{selection.strength}{default to 10^(-2)}
#' \item{normalized_half_life}{default to 10^(-2)}
#' \item{log_likelihood}{default to 10^(-2)}
#' }
#' @param Nbr_It_Max the maximal number of iterations of the EM allowed. Default to
#' 500 iterations.
#' @param nbr_of_shifts the number of shifts allowed.
#' @param alpha_known is the selection strength assumed to be known ?
#' Default to FALSE.
#' @param eps tolerance on the selection strength value before switching to a BM.
#' Default to 10^(-3).
#' @param known.selection.strength if \code{alpha_known=TRUE}, the value of the
#' known selection strength.
#' @param init.selection.strength (optional) a starting point for the selection
#' strength value.
#' @param max_selection.strength the maximal value allowed of the selection strength.
#' Default to 100.
#' @param use_sigma_for_lasso whether to use the first estimation of the variance
#' matrix in the lasso regression. Default to TRUE.
#' @param max_triplet_number for the initialization of the selection strength value
#' (when estimated), the maximal number of triplets of tips to be considered.
#' @param min_params a named list containing the minimum allowed values for the
#' parameters. If the estimation is smaller, then the EM stops, and is considered to
#' be divergent. Default values:
#' \describe{
#' \item{variance}{default to 0}
#' \item{value.root}{default to -10^(5)}
#' \item{exp.root}{default to -10^(5)}
#' \item{var.root}{default to 0}
#' \item{selection.strength}{default to 0}
#' }
#' @param max_params a named list containing the maximum allowed values for the
#' parameters. If the estimation is larger, then the EM stops, and is considered to
#' be divergent. Default values:
#' \describe{
#' \item{variance}{default to 10^(5)}
#' \item{value.root}{default to 10^(5)}
#' \item{exp.root}{default to 10^(5)}
#' \item{var.root}{default to 10^(5)}
#' \item{selection.strength}{default to 10^(5)}
#' }
#' @param var.init.root optional initialization value for the variance of the root.
#' @param variance.init optional initialization value for the variance.
#' @param times_shared (optional) times of shared ancestry of all nodes and tips,
#' result of function \code{\link{compute_times_ca}}
#' @param distances_phylo (optional) phylogenetic distances, result of function 
#' \code{\link{compute_dist_phy}}.
#' @param subtree.list (optional) tips descendants of all the edges, result of
#' function \code{\link{enumerate_tips_under_edges}}.
#' @param T_tree (optional) matrix of incidence of the tree, result of function 
#' \code{\link{incidence.matrix}}.
#' @param U_tree (optional) full matrix of incidence of the tree, result of function
#' \code{\link{incidence.matrix.full}}.
#' @param h_tree (optional) total height of the tree.
#' @param F_moments (optional, internal)
#' @param tol_half_life should the tolerance criterion be applied to the
#' phylogenetic half life (TRUE, default) or to the raw selection strength ?
#' @param warning_several_solutions whether to issue a warning if several equivalent
#' solutions are found (default to TRUE).
#' @param convergence_mode one of "relative" (the default) or "absolute". Should the
#' tolerance be applied to the raw parameters, or to the renormalized ones ?
#' @param check_convergence_likelihood should the likelihood be taken into
#' consideration for convergence assessment ? (default to TRUE).
#' @param method.OUsun Method to be used in univariate OU. One of "rescale" 
#' (rescale the tree to fit a BM) or "raw" (directly use an OU, only available for
#' univariate processes).
#' @param sBM_variance Is the root of the BM supposed to be random and
#' "stationary"? Used for BM equivalent computations. Default to FALSE.
#' @param allow_negative whether to allow negative values for alpha (Early Burst).
#' See documentation of \code{\link{PhyloEM}} for more details. Default to FALSE.
#' @param trait_correlation_threshold the trait correlation threshold to stop the analysis. Default to 0.9.
#' 
#' @return
#' An object of class \code{EstimateEM}.
#' 
#' @seealso \code{\link{PhyloEM}}
#' 
#' @export
#'
##
estimateEM <- function(phylo, 
                       Y_data, 
                       Y_data_imp = Y_data,
                       process = c("BM", "OU", "scOU", "rBM"), 
                       independent = FALSE,
                       tol_EM = list(variance = 10^(-2), 
                                     value.root = 10^(-2), 
                                     exp.root = 10^(-2), 
                                     var.root = 10^(-2),
                                     selection.strength = 10^(-2),
                                     normalized_half_life = 10^(-2),
                                     log_likelihood = 10^(-2)),  
                       Nbr_It_Max = 500, 
                       method.variance = c("simple", "upward_downward"), 
                       method.init = c("default", "lasso"),
                       method.init.alpha = c("default", "estimation"),
                       method.init.alpha.estimation = c("regression", 
                                                        "regression.MM", 
                                                        "median"),
                       nbr_of_shifts = 0,
                       random.root = TRUE,
                       stationary.root = TRUE,
                       alpha_known = FALSE,
                       eps = 10^(-3),
                       known.selection.strength = 1,
                       init.selection.strength = 1,
                       max_selection.strength = 100,
                       use_sigma_for_lasso = TRUE,
                       max_triplet_number = 10000,
                       min_params=list(variance = 0, 
                                       value.root = -10^(5), 
                                       exp.root = -10^(5), 
                                       var.root = 0,
                                       selection.strength = 0),
                       max_params=list(variance = 10^(5), 
                                       value.root = 10^(5), 
                                       exp.root = 10^(5), 
                                       var.root = 10^(5),
                                       selection.strength = 10^(5)),
                       var.init.root = diag(1, nrow(Y_data)),
                       variance.init = diag(1, nrow(Y_data), nrow(Y_data)),
                       methods.segmentation = c(#"max_costs_0", 
                                                "lasso", 
                                                "same_shifts", 
                                                #"same_shifts_same_values",
                                                "best_single_move"),
                                                #"lasso_one_move"),
                       check.tips.names = FALSE,
                       times_shared = NULL, # These can be specified to save time
                       distances_phylo = NULL, 
                       subtree.list = NULL,
                       T_tree = NULL, 
                       U_tree = NULL,
                       h_tree = NULL,
                       F_moments = NULL,
                       tol_half_life = TRUE,
                       warning_several_solutions = TRUE,
                       convergence_mode = c("relative", "absolute"),
                       check_convergence_likelihood = TRUE,
                       sBM_variance = FALSE,
                       method.OUsun = c("rescale", "raw"),
                       # impute_init_Rphylopars = FALSE,
                       K_lag_init = 0,
                       allow_negative = FALSE,
                       trait_correlation_threshold = 0.9,
                       ...){
  
  ntaxa <- length(phylo$tip.label)
  ## Check that the vector of data is in the correct order and dimensions ####
  Y_data <- check_data(phylo, Y_data, check.tips.names)
  trait_correlation_threshold <- check_correlations(Y_data, trait_correlation_threshold)
  ## Find dimension
  p <- nrow(Y_data)
  
  ## Root edge of the tree ###################################################
  if (is.null(phylo$root.edge)) phylo$root.edge <- 0
  
  ########## Check consistancy ################################################
  if (alpha_known && missing(known.selection.strength)) stop("The selection strength alpha is supposed to be known, but is not specified. Please add an argument known.selection.strength to the call of the function.")
  #  known.selection.strength <- check_dimensions.matrix(p, p, known.selection.strength, "known.selection.strength")
  if (independent && 
      alpha_known &&
      !missing(known.selection.strength) && 
      (length(known.selection.strength) != p)){
    warning("The vector of known selection strength provided has not the correct dimension (should be of length p). It will be recycled.")
    known.selection.strength <- rep(known.selection.strength, p)[1:p]
  }
  if (independent && 
      # !missing(init.selection.strength) && 
      (length(init.selection.strength) != p)){
    if (!missing(init.selection.strength)) {
      warning("The vector of init selection strength provided has not the correct dimension (should be of length p). It will be recycled.")
    }
    init.selection.strength <- rep(init.selection.strength, p)[1:p]
  }
  if (sBM_variance){
    if (!random.root){
      warning("Process sBM assumes a random root. Switching to a random root.")
      random.root <- TRUE
    }
    if (phylo$root.edge == 0){
      stop("Something went wrong: I need to know the length of the root edge for the sBM.")
    }
  }
  
  ########## Choose process ###################################################
  process <- match.arg(process)
  original_process <- process
  temp <- choose_process_EM(process = process, p = p,
                            random.root = random.root,
                            stationary.root = stationary.root,
                            alpha_known = alpha_known,
                            known.selection.strength = known.selection.strength,
                            eps = eps,
                            sBM_variance = sBM_variance,
                            method.OUsun = method.OUsun,
                            independent = independent,
                            allow_negative = allow_negative)
  process <- temp$process
  transform_scOU <- temp$transform_scOU # Transform back to get an OU ?
  rescale_tree <- temp$rescale_tree # Rescale the tree ?
  sBM_variance <- temp$sBM_variance
  
  ########## Missing Data #####################################################
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  p <- nrow(Y_data)
  miss <- as.vector(is.na(Y_data))
  Flag_Missing <- any(miss)
  Y_data_vec <- as.vector(Y_data)
  Y_data_vec_known <- as.vector(Y_data[!miss])
  # Vectorized Data Mask
  masque_data <- rep(FALSE, (ntaxa + Nnode) * p)
  masque_data[1:(p*ntaxa)] <- !miss
  
  ########## Choose functions #################################################
  # specialCase <- stationary.root && shifts_at_nodes && alpha_known
  shifts_at_nodes <- TRUE
  compute_M  <- switch(process, 
                       BM = compute_M.BM,
                       OU = compute_M.OU(stationary.root, shifts_at_nodes, alpha_known))
  shutoff.EM  <- switch(process, 
                        BM = shutoff.EM.BM,
                        OU = shutoff.EM.OU(stationary.root, shifts_at_nodes, 
                                           alpha_known, tol_half_life))
  has_converged  <- switch(convergence_mode[1], 
                           relative = has_converged_relative,
                           absolute = has_converged_absolute)
  # is.finite.params  <- switch(process, 
  #                             BM = is.finite.params.BM,
  #                             OU = is.finite.params.OU(stationary.root, shifts_at_nodes, alpha_known))
  is.in.ranges.params  <- switch(process, 
                                 BM = is.in.ranges.params.BM,
                                 OU = is.in.ranges.params.OU(stationary.root, shifts_at_nodes, alpha_known))
  #   compute_MaxCompleteLogLik <- switch(process, 
  #                                       BM = compute_MaxCompleteLogLik.BM,
  #                                       OU = compute_MaxCompleteLogLik.OU(stationary.root, shifts_at_nodes))
  #   conditional_expectation_log_likelihood <- switch(process, 
  #                                                    BM = conditional_expectation_log_likelihood.BM,
  #                                                    OU = conditional_expectation_log_likelihood.OU(stationary.root, shifts_at_nodes))
  
  ########## init alpha #######################################################
  method.init.alpha  <- match.arg(method.init.alpha)
  if (!stationary.root && (method.init.alpha == "estimation")){
    method.init.alpha <- "default"
    warning("The estimation initialization of alpha does only work when the root is stationary. The initialization is set to the default one.")
  }
  # init.alpha <- switch(process, 
  #                      BM = init.alpha.BM,
  #                      OU = init.alpha.OU)
  init.alpha.gamma<- switch(process, 
                            BM = init.alpha.gamma.BM,
                            OU = init.alpha.gamma.OU)
  
  ########## Moments computation method #######################################
  method.variance  <- match.arg(method.variance)
  if (process == "BM" && !Flag_Missing && method.variance == "simple"){
    method.variance <- "simple.nomissing.BM"
  }
  compute_E  <- switch(method.variance, 
                       simple = compute_E.simple,
                       simple.nomissing.BM = compute_E.simple.nomissing.BM,
                       upward_downward = compute_E.upward_downward)
  # compute_mean_variance  <- switch(method.variance, 
  #                                  simple = compute_mean_variance.simple,
  #                                  simple.nomissing.BM = compute_mean_variance.simple.nomissing.BM)
  # compute_log_likelihood  <- switch(method.variance, 
  #                                   simple = compute_log_likelihood.simple,
  #                                   simple.nomissing.BM = compute_log_likelihood.simple.nomissing.BM)
  # compute_mahalanobis_distance  <- switch(method.variance, 
  #                                         simple = compute_mahalanobis_distance.simple,
  #                                         simple.nomissing.BM = compute_mahalanobis_distance.simple.nomissing.BM)
  
  ########## Initialization Method ############################################
  method.init  <- match.arg(method.init)
  # Lasso initialization for OU only works for stationary root
  #   if (!stationary.root && (method.init == "lasso")){
  #     method.init <- "default"
  #     warning("The lasso initialization of alpha does only work when the root is stationary. The initialization is set to the default one.")
  #   }
  init.EM  <- switch(method.init, 
                     default = init.EM.default(process),
                     lasso = init.EM.lasso)
  method.init.alpha  <- match.arg(method.init.alpha)
  methods.segmentation <- match.arg(methods.segmentation, several.ok = TRUE)
  # if (independent){
  #   tmp <- methods.segmentation %in% c("lasso_one_move")
  #   if (any(tmp)){
  #     warning("Lasso segmentation methods are not implemented for multivariate independent OU. Removing these methods from the list.")
  #     methods.segmentation <- methods.segmentation[!tmp]
  #     if (length(methods.segmentation) == 0) {
  #       warning("The list of segmentations methods was empty. Adding the best single move method.")
  #       methods.segmentation <- "best_single_move"
  #     }
  #   }
  # }
  method.init.alpha.estimation  <- match.arg(method.init.alpha.estimation,
                                             several.ok = TRUE)
  
  ########## Fixed Quantities #################################################
  ntaxa <- length(phylo$tip.label)
  ## Transform the branch lengths if needed
  phy_original <- phylo
  if (rescale_tree){
    if (sBM_variance) phylo$root.edge <- 1
    phylo <- transform_branch_length(phylo, known.selection.strength)
  }
  if (is.null(times_shared)) times_shared <- compute_times_ca(phylo)
  if (is.null(distances_phylo)) distances_phylo <- compute_dist_phy(phylo)
  if (is.null(subtree.list)) subtree.list <- enumerate_tips_under_edges(phylo)
  if (is.null(T_tree)) T_tree <- incidence.matrix(phylo)
  if (is.null(U_tree)) U_tree <- incidence.matrix.full(phylo)
  if (is.null(h_tree)) h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
  if ((is.null(F_moments)) 
      && (!Flag_Missing) 
      && (process == "BM")
      && (method.variance != "upward_downward")){
    # Add root edge to the branch lengths (root assumed fixed by default)
    root_edge_length <- 0
    if (!is.null(phylo$root.edge)) root_edge_length <- phylo$root.edge
    F_moments <- compute_fixed_moments(times_shared + root_edge_length, ntaxa)
  } 
  
  ########## Re-scale tree to 100 #############################################
  
  # factor_rescale <- min(phy_original$edge.length) / min(phylo$edge.length) # min new = min old
  factor_rescale <- 1 / h_tree # total height to 1
  # factor_rescale <- 1
  
  h_tree <- factor_rescale * h_tree
  times_shared <- factor_rescale * times_shared
  distances_phylo <- factor_rescale * distances_phylo
  phylo$edge.length <- factor_rescale * phylo$edge.length
  phylo$root.edge <- factor_rescale * phylo$root.edge
  
  if (!Flag_Missing && process == "BM"){
    F_moments$C_YY = F_moments$C_YY * factor_rescale
    F_moments$C_YY_chol_inv = F_moments$C_YY_chol_inv / sqrt(factor_rescale)
    F_moments$F_vars = F_moments$F_vars * factor_rescale
  }
  
  known.selection.strength <-  known.selection.strength / factor_rescale
  init.selection.strength <- init.selection.strength / factor_rescale
  variance.init <- variance.init / factor_rescale

  ########## Initialization of alpha and Variance #############################
  init.a.g <- init.alpha.gamma(method.init.alpha)(phylo = phylo,
                                                  Y_data = Y_data,
                                                  nbr_of_shifts = nbr_of_shifts,
                                                  times_shared = times_shared,
                                                  distances_phylo = distances_phylo,
                                                  init.selection.strength = init.selection.strength,
                                                  max_triplet_number = max_triplet_number,
                                                  known.selection.strength = known.selection.strength,
                                                  alpha_known = alpha_known,
                                                  init.var.root = var.init.root,
                                                  method.init.alpha.estimation = method.init.alpha.estimation,
                                                  tol_EM = tol_EM,
                                                  h_tree = h_tree,
                                                  random.root = random.root,
                                                  T_tree = T_tree,
                                                  subtree.list = subtree.list,
                                                  miss = miss,
                                                  masque_data = masque_data,
                                                  independent = independent,
                                                  ...)
  if (is.null(init.a.g$gamma_0)){
    init.var.root <- NULL
  } else {
    init.var.root <- colMeans(init.a.g$gamma_0, na.rm = TRUE)
    init.var.root <- diag(init.var.root, ncol = length(init.var.root))
  }
  if (process == "OU"){
    if(!alpha_known) {
      if ((sum(is.finite(init.a.g$alpha_0)) != 0)){
        ## Only if not all NAs ot infinite
        tmp <- colMeans(init.a.g$alpha_0, na.rm = TRUE)
        if (anyNA(tmp)){
          warning("Estimating selection strength failed for some traits. Replacing corresponding values of alpha by default ones.")
          tmp[is.na(tmp)] <- init.selection.strength[is.na(tmp)]
        }
        init.selection.strength <- tmp
      }
    } else {
      init.selection.strength <- known.selection.strength
    }
  }
  if (process == "BM"){
    init.selection.strength <- rep(0, p)
  }
  
  ## Init of Rate matrix for BM
  if (process == "BM" && (!random.root || (random.root && sBM_variance))){
    # if (impute_init_Rphylopars && any(is.na(Y_data_imp))){
    #   # Impute the data only once, if needed.
    #   message("Imputing data for lasso initialization.")
    #   Y_data_imp <- impute.data.Rphylopars(phylo, Y_data, process, random.root)
    # }
    variance.init <- init.variance.BM.estimation(phylo = phylo, 
                                                 Y_data = Y_data, 
                                                 Y_data_imp = Y_data_imp,
                                                 Y_data_vec_known = Y_data_vec_known,
                                                 nbr_of_shifts = nbr_of_shifts, 
                                                 times_shared = times_shared,
                                                 distances_phylo = distances_phylo, 
                                                 h_tree = h_tree,
                                                 random.root = random.root,
                                                 T_tree = T_tree,
                                                 subtree.list = subtree.list,
                                                 miss = miss,
                                                 selection.strength.init = init.selection.strength,
                                                 # impute_init_Rphylopars = impute_init_Rphylopars,
                                                 masque_data = masque_data,
                                                 ...)
  }
  
  ########## Initialization of all parameters #################################
  params_init <- init.EM(phylo = phylo,
                         Y_data = Y_data,
                         Y_data_imp = Y_data_imp,
                         Y_data_vec_known = Y_data_vec_known,
                         process = process, 
                         times_shared = times_shared, 
                         distances_phylo = distances_phylo, 
                         nbr_of_shifts = nbr_of_shifts, 
                         K_lag_init = K_lag_init,
                         selection.strength.init = init.selection.strength, 
                         random.init = random.root,
                         stationary.root.init = stationary.root,
                         use_sigma = use_sigma_for_lasso,
                         method.init.alpha = method.init.alpha,
                         var.root.init = init.var.root,
                         T_tree = T_tree,
                         subtree.list = subtree.list,
                         miss = miss,
                         variance.init = variance.init,
                         sBM_variance = sBM_variance,
                         # impute_init_Rphylopars = impute_init_Rphylopars,
                         masque_data = masque_data,
                         independent = independent,
                         ...)
  params <- params_init
  params$root.state <- test.root.state(root.state = params$root.state, 
                                       process = process, 
                                       optimal.value = params$optimal.value,
                                       variance = params$variance, 
                                       selection.strength = params$selection.strength)
  attr(params, "ntaxa")  <- ntaxa
  attr(params, "p_dim")  <- p
  params_old <- NULL
  
  ########## Iterations #######################################################
  Nbr_It <- 0
  params_history <- vector("list")#, Nbr_It_Max)
  #   CLL_history <- NULL
  number_new_shifts <- NULL
  CV_log_lik <- FALSE
  while ( Nbr_It == 0 || # initialization
          K_lag_init > 0 || 
          ( !(CV_log_lik && # CV of log-Likelihood ?
              shutoff.EM(params_old, params, tol_EM, has_converged, h_tree)) && # Shutoff
            is.in.ranges.params(params, min = min_params, max = max_params) && #Divergence?
            Nbr_It < Nbr_It_Max ) ) { # Nbr of iteration
    ## Actualization
    Nbr_It <- Nbr_It + 1
    params_old <- params
    if (K_lag_init > 0) K_lag_init <- K_lag_init - K_lag_init # reduce lag progressivelly
    ########## Log Likelihood and E step ####################################
    # Check convergence of loglik ?
    if (check_convergence_likelihood && Nbr_It > 1){
      log_likelihood_old_old <- log_likelihood_old
    }    
    # Independent ?
    if (independent){
      params_old <- split_params_independent(params_old)
    }
    ## result
    temp <- wrapper_E_step(phylo = phylo,
                           times_shared = times_shared,
                           distances_phylo = distances_phylo,
                           process = process,
                           params_old = params_old,
                           masque_data = masque_data,
                           F_moments = F_moments,
                           independent = independent,
                           Y_data_vec_known = Y_data_vec_known,
                           miss = miss,
                           Y_data = Y_data,
                           U_tree = U_tree,
                           # compute_mean_variance = compute_mean_variance,
                           # compute_log_likelihood = compute_log_likelihood,
                           # compute_mahalanobis_distance = compute_mahalanobis_distance,
                           compute_E = compute_E)
    ## Format result if independent
    if (independent){
      log_likelihood_old <- sum(sapply(temp, function(z) return(z$log_likelihood_old)))
      # Store params for history
      params_history[[paste(Nbr_It - 1, sep="")]] <- merge_params_independent(params_old)
      attr(params_history[[paste(Nbr_It - 1, sep="")]], "log_likelihood") <- log_likelihood_old
      # attr(params_history[[paste(Nbr_It - 1, sep="")]], "mahalanobis_distance_data_mean") <- sum(sapply(temp, function(z) return(z$maha_data_mean)))
      # Conditional law
      conditional_law_X <- lapply(temp, function(z) return(z$conditional_law_X))
    } else {
      log_likelihood_old <- as.vector(temp$log_likelihood_old)
      attr(params_old, "log_likelihood") <- log_likelihood_old
      # attr(params_old, "mahalanobis_distance_data_mean") <- as.vector(temp$maha_data_mean)
      # Store params for history
      params_history[[paste(Nbr_It - 1, sep="")]] <- params_old
      # Conditional law
      conditional_law_X <- temp$conditional_law_X
    }
    rm(temp)
    if (check_convergence_likelihood && Nbr_It > 1){
      CV_log_lik <- has_converged_absolute(log_likelihood_old_old,
                                           log_likelihood_old,
                                           tol_EM$log_likelihood)
    }
    
    #   ## Log likelihood as the sum of conditional + entropy
    #         H <- compute_entropy.simple(moments$Sigma, moments$Sigma_YY_inv)
    #         CLL <- conditional_expectation_log_likelihood(phylo = phylo,
    #                                               conditional_law_X = conditional_law_X, 
    #                                               sigma2 = params_old$variance,
    #                                               mu = params_old$root.state$exp.root,
    #                                               shifts = params_old$shifts,
    #                                               alpha = params_old$selection.strength)
    #         log_likelihood_bis <- compute_log_likelihood_with_entropy.simple(CLL, H)
    #         attr(params_old, "log_likelihood_bis") <- log_likelihood_bis
    
    
    ########## M step #########################################################
    if (independent){
      # correct format for selection strength case independent
      alpha_old = sapply(params_old, function(z) z$selection.strength)
    } else {
      alpha_old = params_old$selection.strength
    }
    params <- compute_M(phylo = phylo, 
                        Y_data = Y_data, 
                        conditional_law_X = conditional_law_X, 
                        nbr_of_shifts = nbr_of_shifts + K_lag_init, 
                        random.root = random.root,
                        known.selection.strength = known.selection.strength,
                        alpha_old = alpha_old,
                        max_selection.strength = max_selection.strength,
                        eps = eps,
                        methods.segmentation = methods.segmentation,
                        beta_0_old = params_history[[paste(Nbr_It - 1, sep="")]]$optimal.value,
                        shifts_old = params_history[[paste(Nbr_It - 1, sep="")]]$shifts,
                        variance_old = params_history[[paste(Nbr_It - 1, sep="")]]$variance,
                        mu_old = params_history[[paste(Nbr_It - 1, sep="")]]$root.state$value.root,
                        subtree.list = subtree.list,
                        sBM_variance = sBM_variance,
                        params_old = params_old)
    # If independent, go back to merged parameters.
    if (independent){
      if (p > 1) params <- merge_params_independent(params)
      params_old <- params_history[[paste(Nbr_It - 1, sep="")]]
    }
    attr(params, "ntaxa")  <- ntaxa
    attr(params, "p_dim")  <- p
    ## Number of shifts that changed position ?
    if (independent){
      number_new_shifts <- c(number_new_shifts,
                             sum(!(params$shifts$edges %in% params_old[[1]]$shifts$edges)))
    } else {
      number_new_shifts <- c(number_new_shifts,
                             sum(!(params$shifts$edges %in% params_old$shifts$edges))) 
    }
  }
  
  ########## Scale back parameters to original tree ###########################
  
  params <- scale_params(params, factor_rescale)
  
  h_tree <- h_tree / factor_rescale
  times_shared <- times_shared / factor_rescale
  distances_phylo <- distances_phylo / factor_rescale
  phylo$edge.length <- phylo$edge.length / factor_rescale
  phylo$root.edge <- phylo$root.edge / factor_rescale
  
  if (!Flag_Missing && process == "BM"){
    F_moments$C_YY = F_moments$C_YY / factor_rescale
    F_moments$C_YY_chol_inv = F_moments$C_YY_chol_inv * sqrt(factor_rescale)
    F_moments$F_vars = F_moments$F_vars / factor_rescale
  }
  
  known.selection.strength <-  known.selection.strength * factor_rescale
  init.selection.strength <- init.selection.strength * factor_rescale
  variance.init <- variance.init * factor_rescale
  
  ########## Go back to OU parameters if needed ###############################
  if (transform_scOU){
    ## Go back to original tree and process
    phylo <- phy_original
    process <- suppressWarnings(check.selection.strength(original_process,
                                                         known.selection.strength,
                                                         eps))
    times_shared <- compute_times_ca(phy_original)
    distances_phylo <- compute_dist_phy(phy_original)
    ## Compute equivalent parameters
    params_scOU <- go_back_to_original_process(phy_original = phy_original,
                                               known.selection.strength = known.selection.strength,
                                               sBM_variance = sBM_variance,
                                               params = params)
  } else {
    params_scOU <- params # If a BM, params_scOU = params
  }
  
  ########## Compute scores and ancestral states for final parameters ##########
  if ((original_process %in% c("OU", "scOU"))
      && (method.variance == "simple.nomissing.BM")){
    ## Go back to simple method if switched to a different one.
    compute_E <- compute_E.simple
    # compute_mean_variance  <- compute_mean_variance.simple
    # compute_log_likelihood  <- compute_log_likelihood.simple
    # compute_mahalanobis_distance  <- compute_mahalanobis_distance.simple
  }

  if (independent){
    params_scOU <- split_params_independent(params_scOU)
  }

  temp <- wrapper_E_step(phylo = phylo,
                         times_shared = times_shared,
                         distances_phylo = distances_phylo,
                         process = process,
                         params_old = params_scOU,
                         masque_data = masque_data,
                         F_moments = F_moments,
                         independent = independent,
                         Y_data_vec_known = Y_data_vec_known,
                         miss = miss,
                         Y_data = Y_data,
                         U_tree = U_tree,
                         # compute_mean_variance = compute_mean_variance,
                         # compute_log_likelihood = compute_log_likelihood,
                         # compute_mahalanobis_distance = compute_mahalanobis_distance,
                         compute_E = compute_E)
  ## Format results
  if (independent){ # Independent
    ## Likelihood and Mahalanobis of last parameters
    params_scOU <- merge_params_independent(params_scOU)
    attr(params_scOU, "log_likelihood") <- sum(sapply(temp, function(z) return(z$log_likelihood_old)))
    # attr(params_scOU, "mahalanobis_distance_data_mean") <- sum(sapply(temp, function(z) return(z$maha_data_mean)))
    params_history[[paste(Nbr_It, sep="")]] <- params_scOU
    ## "Ancestral States Reconstruction"
    condlaw <- lapply(temp, function(z) return(z$conditional_law_X))
    condlaw <- do.call(rbind, condlaw)
    conditional_law_X <- vector("list")
    conditional_law_X$expectations <- do.call(rbind, condlaw[, "expectations"])
    conditional_law_X$optimal.values <- do.call(rbind, condlaw[, "optimal.values"])
    conditional_law_X$variances <- do.call(rbind, condlaw[, "variances"])
    conditional_law_X$covariances <- do.call(rbind, condlaw[, "covariances"])
    if (p > 1){
      conditional_law_X$variances <- plyr::aaply(conditional_law_X$variances, 2,
                                                 diag)
      conditional_law_X$variances <- aperm(conditional_law_X$variances, c(2, 3, 1))
      conditional_law_X$covariances <- plyr::aaply(conditional_law_X$covariances, 2, diag)
      conditional_law_X$covariances <- aperm(conditional_law_X$covariances, c(2, 3, 1))
    } else {
      conditional_law_X$variances <- array(conditional_law_X$variances,
                                           c(1, 1, ncol(conditional_law_X$variances)))
      conditional_law_X$covariances <- array(conditional_law_X$covariances,
                                             c(1, 1, ncol(conditional_law_X$covariances)))
    }
    rm(condlaw)
    ## Mean at tips with estimated parameters
    # m_Y_estim <- extract_simulate_internal(tmpsim, where="tips", what="expectations")
    # m_Y_estim <- lapply(temp, function(z) extract_simulate_internal(z$moments$sim, where="tips", what="expectations"))
    # m_Y_estim <- do.call(rbind, m_Y_estim)
  } else { ## NOT independent
    ## Likelihood and Mahalanobis of last parameters
    attr(params_scOU, "log_likelihood") <- as.vector(temp$log_likelihood_old)
    # attr(params_scOU, "mahalanobis_distance_data_mean") <- as.vector(temp$maha_data_mean)
    params_history[[paste(Nbr_It, sep="")]] <- params_scOU
    ## "Ancestral States Reconstruction"
    conditional_law_X <- temp$conditional_law_X
    ## Mean at tips with estimated parameters
    # m_Y_estim <- temp$conditional_law_X$expectation[, 1:ntaxa]
    # m_Y_estim <- extract_simulate_internal(tmpsim, where="tips", what="expectations")
  }
  rm(temp)
  tmpsim <- simulate_internal(phylo = phylo,
                              process = process,
                              p = attr(params_scOU, "p_dim"),
                              root.state = params_scOU$root.state,
                              shifts = params_scOU$shifts,
                              variance = params_scOU$variance,
                              optimal.value = params_scOU$optimal.value,
                              selection.strength = params_scOU$selection.strength,
                              simulate_random = FALSE,
                              U_tree = U_tree)
  ## Mean at tips with estimated parameters
  m_Y_estim <- extract_simulate_internal(tmpsim, where="tips", what="expectations")
  rm(tmpsim)
  
  ########## Number of equivalent solutions ###################################
  clusters <- clusters_from_shifts(phylo, params$shifts$edges,
                                   part.list = subtree.list)
  Neq <- extract.parsimonyNumber(parsimonyNumber(phylo, clusters))
  if (Neq > 1 && warning_several_solutions) message("There are some equivalent solutions to the solution found.")
  attr(params_scOU, "Neq") <- Neq
  
  ### lsq ####
  if (!independent){
    lsq <- sum(diag(params_scOU$variance))
  } else {
    pp <- split_params_independent(params_scOU)
    lsq <- sapply(pp, function(z) sum(diag(z$variance)))
    lsq <- sum(lsq)
  }
  
  ########## Result  ##########################################################
  conditional_law_X$expectations <- matrix(conditional_law_X$expectations, nrow = p)
  result <- list(params = params_scOU, # Return untransformed parameters as default
                 params_raw = params,
                 ReconstructedNodesStates = conditional_law_X$expectations[ , (ntaxa+1):ncol(conditional_law_X$expectations), drop = FALSE],
                 ReconstructedTipsStates = conditional_law_X$expectations[ , 1:ntaxa, drop = FALSE],
                 ReconstructedNodesVariances = conditional_law_X$variances[ , , (ntaxa+1):ncol(conditional_law_X$expectations), drop = FALSE],
                 ReconstructedTipsVariances = conditional_law_X$variances[ , , 1:ntaxa, drop = FALSE],
                 m_Y_estim = m_Y_estim,
                 params_old = params_old, 
                 params_init = params_init,
                 alpha_0 = init.a.g$alpha_0,
                 gamma_0 = init.a.g$gamma_0,
                 params_history = params_history,
                 number_new_shifts = number_new_shifts,
                 number_equivalent_solutions = Neq,
                 least_squares_raw = sum((Y_data - m_Y_estim)^2, na.rm = TRUE),
                 least_squares = lsq)
  #                  CLL_history = CLL_history
  if (transform_scOU) result$params_scOU <-  params_scOU
  ## Handle convergence
  attr(result, "Nbr_It") <- Nbr_It
  attr(result, "Divergence") <- !is.in.ranges.params(result$params_raw,
                                                     min = min_params,
                                                     max = max_params) # TRUE if has diverged (use raw parameters for diagnostic)
  if (Nbr_It == Nbr_It_Max) warning(paste("The maximum number of iterations (Nbr_It_Max = ",Nbr_It_Max,") was reached.",sep=""))
  return(result)
}

###############################################################################
## PhyloEM
###############################################################################

##
#' @title Model Estimation with Detection of Shifts
#'
#' @description
#' \code{PhyloEM} is the main function of the package. It uses maximum likelihood
#' methods to fit a BM or an OU process for several traits evolving along a
#' phylogenetic tree, with automatic shift detection on the branches of the tree.
#' This function can handle missing data.
#'
#' @details
#' Several models can be used:
#' \itemize{
#' \item BM with fixed root, univariate or multivariate.
#' \item OU with fixed or stationary root, univariate or multivariate.
#' }
#' For the OU in the multivariate setting, two assumptions can be made:
#' \itemize{
#' \item Independent traits. This amounts to diagonal rate and selection matrices.
#' \item "Scalar OU" (scOU): the rate matrix can be full, but the selection 
#' strength matrix is assumed to be scalar, i.e. all the traits are supposed to
#' go to their optimum values with the same speed.
#' }
#' 
#' Note that the "scalar OU" model can also be seen as a re-scaling of the tree.
#' The selection strength parameter alpha can then be interpreted as a measure
#' of the "phylogenetic signal":
#' \itemize{
#' \item If alpha is close to 0, then the process is similar to a BM on the original tree,
#' and the signal is strong.
#' \item If alpha is large, then the re-scaled tree is similar to a star-tree,
#' and the signal is weak.
#' }
#' When there are no shifts, and the root is taken to be constant, this
#' model is actually equivalent to an AC model (Uyeda et al. 2015).
#' With this interpretation in mind, one might want to explore 
#' negative values of alpha, in order to fit a DC (or Early Burst) model.
#' With no shift and a fixed root, the same proof shows that the scOU
#' with alpha negative is equivalent to the DC model. There are two
#' strong caveats in doing that.
#' \itemize{
#' \item The interpretation of the OU as modeling the dynamic of a trait
#' undergoing stabilizing selection is lost. In this case, the scOU can only
#' be seen as a re-scaling of the tree, similar to Pagel's delta.
#' \item The values of the "optimal values", and of the shifts on them, cannot
#' be interpreted as such (the process is actually going away from this values,
#' instead of being attracted). When looking at these values, one should only use
#' the un-normalized values happening of the underlying BM. You can extract those
#' using the \code{\link{params_process}} function with \code{rBM = TRUE}.
#' }
#'
#' @param phylo A phylogenetic tree of class \code{phylo} 
#' (from package \code{\link{ape}}).
#' @param Y_data Matrix of data at the tips, size p x ntaxa. Each line is a
#' trait, and each column is a tip. The column names are checked against the
#' tip names of the tree.
#' @param process The model used for the fit. One of "BM" (for a full BM model, 
#' univariate or multivariate); "OU" (for an OU with independent traits, 
#' univariate or multivariate); or "scOU" (for a "scalar OU" model, see details).
#' @param check_postorder Re-order the tree in post-order. If the Upward-Downward
#' algorithm is used, the tree need to be in post-order. Default to TRUE if the
#' upward-downward is used, otherwise automatically set to FALSE.
#' @param independent Are the trait assumed to be independent from one another?
#' Default to FALSE. OU in a multivariate setting only works if TRUE.
#' @param K_max The maximum number of shifts to be considered. Default to 
#' \eqn{max(|\sqrt ntaxa|, 10)}.
#' @param use_previous Should the initialization for K+1 shifts use the 
#' estimation for $K$ shifts already obtained? Default to FALSE.
#' @param order Should the estimations be done for K increasing (TRUE) or K
#' decreasing (FALSE)? If use_previous=FALSE, this has no influence, except if one
#' initialization fails. Default to TRUE.
#' @param method.selection Method selection to be used. Several ones can be
#' used at the same time. One of "LINselect" for the Baraud Giraud Huet LINselect 
#' method; "DDSE" for the Slope Heuristic or "Djump" for the Jump Heuristic, last
#' two based the Birgé Massart method.
#' @param C.BM1 Multiplying constant to be used for the BigeMassart1 method.
#' Need to be positive. Default to 0.1.
#' @param C.BM2 Multiplying constant to be used for the BigeMassart2 method.
#' Default to 2.5.
#' @param C.LINselect Multiplying constant to be used for the LINselect method.
#' Need to be greater than 1. Default to 1.1.
#' @param method.variance Algorithm to be used for the moments computations at the
#' E step. One of "simple" for the naive method; of "upward_downward" for the 
#' Upward Downward method (usually faster). Default to "upward_downward".
#' @param method.init The initialization method. One of "lasso" for the LASSO
#' base initialization method; or "default" for user-specified initialization
#' values. Default to "lasso".
#' @param method.init.alpha For OU model, initialization method for the selection
#' strength alpha. One of "estimation" for a cherry-based initialization, using
#' \code{\link[robustbase]{nlrob}}; or "default" for user-specified 
#' initialization values. Default to "estimation".
#' @param method.init.alpha.estimation If method.init.alpha="estimation",
#' choice of the estimation(s) methods to be used. Choices among "regression",
#' (method="M" is passed to \code{\link[robustbase]{nlrob}}); "regression.MM"
#' (method="MM" is passed to \code{\link[robustbase]{nlrob}}) or "median"
#' (\code{\link[robustbase]{nlrob}} is not used, a simple median is taken).
#' Default to all of them.
#' @param methods.segmentation For OU, method(s) used at the M step to find new
#' candidate shifts positions. Choices among "lasso" for a LASSO-based algorithm;
#' and "best_single_move" for a one-move at a time based heuristic. Default to 
#' both of them. Using only "lasso" might speed up the function a lot.
#' @param alpha_grid whether to use a grid for alpha values. Default to TRUE. This
#' is the only available method for scOU. This method is not available for OU with
#' multivariate traits. OU with univariate traits can take both TRUE or FALSE. If
#' TRUE, a grid based on the branch length of the tree is automatically computed,
#' using function \code{\link{find_grid_alpha}}.
#' @param nbr_alpha If \code{alpha_grid=TRUE}, the number of alpha values on the
#' grid. Default to 10.
#' @param random.root whether the root is assumed to be random (TRUE) of fixed
#' (FALSE). Default to TRUE
#' @param stationary.root whether the root is assumed to be in the stationary 
#' state. Default to TRUE.
#' @param alpha If the estimation is done with a fixed alpha (either known, or
#' on a grid), the possible value for alpha. Default to NULL.
#' @param check.tips.names whether to check the tips names of the tree against
#' the column names of the data. Default to TRUE.
#' @param progress.bar whether to display a progress bar of the computations.
#' Default to TRUE.
#' @param estimates The result of a previous run of this same function. This
#' function can be re-run for other model election method. Default to NULL.
#' @param save_step If alpha_grid=TRUE, whether to save the intermediate results
#' for each value of alpha (in a temporary file). Useful for long computations.
#' Default to FALSE.
# @param sBM_variance DEPRECATED. Used for BM equivalent computations. 
# Default to FALSE.
#' @param rescale_OU For the Univariate OU, should the tree be re-scaled to use a BM ? 
#' This can speed up the computations a lot. However, it can make it harder for the EM to 
#' explore the space of parameters, and hence lead to a sub-optimal solution.
#' Default to TRUE.
#' @param parallel_alpha If alpha_grid=TRUE, whether to run the 
#' estimations with different values of alpha on separate cores. Default to 
#' FALSE. If TRUE, the log is written as a temporary file.
#' @param Ncores If parallel_alpha=TRUE, number of cores to be used.
# @param exportFunctions DEPRECATED. TO BE REMOVED.
# @param impute_init_Rphylopars whether to use 
# \code{\link[Rphylopars]{Rphylopars-package}} for initialization. 
# Default to FALSE.
#' @param K_lag_init Number of extra shifts to be considered at the initialization
#' step. Increases the accuracy, but can make computations quite slow of taken
#' too high. Default to 5.
#' @param light_result if TRUE (the default), the object returned is made light,
#' without easily computable quantities. If FALSE, the object can be very heavy, but
#' its subsequent manipulations can be faster (especially for plotting).
#' @param tol_tree tolerance to consider a branch length significantly greater than zero, or
#' two lineages lengths to be different, when checking for ultrametry. 
#' (Default to .Machine$double.eps^0.5). See \code{\link{is.ultrametric}} and \code{\link{di2multi}}.
#' @param allow_negative whether to allow negative values for alpha (Early Burst).
#' See details. Default to FALSE.
#' @param option_is.ultrametric option for \code{\link{is.ultrametric}} check. Default to 1.
#' @param trait_correlation_threshold the trait correlation threshold to stop the analysis. Default to 0.9.
#' @param ... Further arguments to be passed to \code{\link{estimateEM}}, including
#' tolerance parameters for stopping criteria, maximal number of iterations, etc.
#' 
#' 
#' @return
#' An object of class \code{PhyloEM}. Relevant quantities can be extracted from it 
#' using helper functions \code{\link{params_process.PhyloEM}},
#' \code{\link{imputed_traits.PhyloEM}}
#' 
#' @seealso \code{\link{plot.PhyloEM}}, \code{\link{params_process.PhyloEM}},
#' \code{\link{imputed_traits.PhyloEM}}
#' 
#' @examples
#' \dontrun{
#' ## Load Data
#' data(monkeys)
#' ## Run method
#' # Note: use more alpha values for better results.
#' res <- PhyloEM(Y_data = monkeys$dat,        ## data
#'                phylo = monkeys$phy,         ## phylogeny
#'                process = "scOU",            ## scalar OU
#'                random.root = TRUE,          ## root is stationary
#'                stationary.root = TRUE,
#'                K_max = 10,                  ## maximal number of shifts
#'                nbr_alpha = 4,               ## number of alpha values
#'                parallel_alpha = TRUE,       ## parallelize on alpha values
#'                Ncores = 2)
#' ## Plot selected solution (LINselect)
#' plot(res) # three shifts
#' ## Plot selected solution (DDSE)
#' plot(res, method.selection = "DDSE") # no shift
#' ## Extract and solution with 5 shifts
#' params_5 <- params_process(res, K = 5)
#' plot(res, params = params_5)
#' ## Show all equivalent solutions
#' eq_sol <- equivalent_shifts(monkeys$phy, params_5)
#' plot(eq_sol)
#' }
#' 
#' @export
#'
##
# @return summary a data frame with K_max lines, and columns:
#    - alpha_estim the estimated selection strength
#    - gamma_estim the estimated root variance
#    - beta_0_estim the estimated value of root optimum
#    - EM_steps number of iterations needed before convergence
#    - DV_estim has the EM diverged ?
#    - CV_estim has the EM converged ?
#    - log_likelihood log likelihood of the data using the estimated parameters
#    - mahalanobis_distance_data_mean the Mahalanobis distance between the data
# and the estimated means at the tips
#    - least_squares the Mahalanobis distance, renormalized by gamma^2: 
# mahalanobis_distance_data_mean * gamma_estim.
#    - mean_number_new_shifts the mean number of shifts that changed over the 
# iterations of the EM
#    - number_equivalent_solutions the number of equivalent solutions to 
# the solution found.
#    - K_try the number of shifts allowed.
#    - complexity the complexity for K_try
#    - time the CPU time needed.
# @return params a list of inferred parameters for each EM.

PhyloEM <- function(phylo, Y_data, process = c("BM", "OU", "scOU", "rBM"),
                    check_postorder = TRUE,
                    independent = FALSE,
                    K_max = max(floor(sqrt(length(phylo$tip.label))), 10),
                    use_previous = FALSE,
                    order = TRUE,
                    method.selection = c("LINselect", "DDSE", "Djump"),
                    C.BM1 = 0.1, C.BM2 = 2.5, C.LINselect = 1.1,
                    method.variance = c("upward_downward", "simple"),
                    method.init = "lasso",
                    method.init.alpha = "estimation",
                    method.init.alpha.estimation = c("regression", "regression.MM", "median"), 
                    methods.segmentation = c("lasso", "best_single_move"),
                    alpha_grid = TRUE,
                    nbr_alpha = 10,
                    random.root = TRUE,
                    stationary.root = random.root,
                    alpha = NULL,
                    check.tips.names = TRUE,
                    progress.bar = TRUE,
                    estimates = NULL,
                    save_step = FALSE,
                    # sBM_variance = FALSE,
                    rescale_OU = TRUE,
                    parallel_alpha = FALSE,
                    Ncores = 3,
                    # exportFunctions = ls(),
                    # impute_init_Rphylopars = FALSE,
                    K_lag_init = 5,
                    light_result = TRUE,
                    tol_tree = .Machine$double.eps^0.5,
                    allow_negative = FALSE,
                    option_is.ultrametric = 1,
                    trait_correlation_threshold = 0.9,
                    ...){
  ## Required packages
  # library(doParallel)
  # library(foreach)
  # library(ape)
  # library(glmnet) # For Lasso initialization
  # library(robustbase) # For robust fitting of alpha
  ## Check the tree  ##########################################################
  process <- match.arg(process)
  if ((process != "BM") && !is.ultrametric(phylo, tol = tol_tree, option = option_is.ultrametric)) stop("The tree must be ultrametric.")
  if (any(abs(phylo$edge.length) < tol_tree)){
    stop("The tree has zero-length branches.
         Please use `ape::di2multi` function to transform the zero-length branches into ploytomies.")
  }
  if (any(phylo$edge.length < 0)){
    stop("The tree has negative branch lengths. This is not allowed.")
  }
  phylo_given <- phylo
  method.variance  <- match.arg(method.variance)
  if (method.variance == "simple") check_postorder <- FALSE
  if (check_postorder){
    phylo_original_order <- phylo
    phylo <- reorder(phylo, "postorder") 
  }
  ## Check that the vector of data is in the correct order and dimensions #####
  Y_data <- check_data(phylo, Y_data, check.tips.names)
  trait_correlation_threshold <- check_correlations(Y_data, trait_correlation_threshold)
  p <- nrow(Y_data)
  ntaxa <- length(phylo$tip.label)
  ## Independent traits #######################################################
  method.OUsun <- "rescale"
  if (p == 1) {
    if (!rescale_OU) method.OUsun <- "raw"
    independent <- FALSE
  }
  if (independent && missing(alpha_grid) && missing(nbr_alpha)) alpha_grid <- FALSE
  ## Adaptations to the BM ####################################################
  if (process == "BM"){
    if (independent){
      if (p > 1) warning("The independent option is not available for the BM. The traits are supposed to be correlated.")
      independent <- FALSE
    }
    alpha_grid <- TRUE
    alpha <- 0
  }
  ## Model Selection ##########################################################
  method.selection  <- match.arg(method.selection,
                                 choices = c("LINselect", "DDSE", "Djump",
                                             "BirgeMassart1", "BirgeMassart2",
                                             "BGH", "BGHuni", "BGHlsq", "BGHml",
                                             "BGHlsqraw", "BGHmlraw",
                                             "pBIC", "pBIC_l1ou"),
                                 several.ok = TRUE)
  method.selection <- expand_method_selection(method.selection)
  if (p > 1){
    method.selection <- method.selection[method.selection != "BGH"]
    method.selection <- method.selection[method.selection != "BGHuni"]
    # warning("BGH is not implemented for multivariate data.")
  }
  if (p == 1){
    if ("BGH" %in% method.selection){
      method.selection[method.selection == "BGH"] <- "BGHuni"
    }
    method.selection <- method.selection[method.selection != "BGHlsq"]
    method.selection <- method.selection[method.selection != "BGHml"]
    method.selection <- method.selection[method.selection != "BGHlsqraw"]
    method.selection <- method.selection[method.selection != "BGHmlraw"]
  }
  if ("BirgeMassart1" %in% method.selection || "BirgeMassart2" %in% method.selection){
    if (K_max < 10){
      warning("Slope and Jump heuristics need at least 10 observations. Consider choosing K_max >= 10, or an other model selection.")
      method.selection <- method.selection[method.selection != "BirgeMassart1"]
      method.selection <- method.selection[method.selection != "BirgeMassart2"]
    }
    # library(capushe) 
  }
  if (!alpha_grid){
    method.selection <- method.selection[method.selection != "BGHlsq"]
    method.selection <- method.selection[method.selection != "BGHlsqraw"]
  }
  if (length(method.selection) == 0) stop("No selection method were selected or suited to the problem (see relevant warnings). Please fix before carying on.")
  
  ## Inference per se
  if (!is.null(estimates)){ # If the user already has the estimates
    X <- estimates ## Get estimations from presiously computed results
    check_postorder <- FALSE
  } else if (alpha_grid) { # For a grid estimation of alpha
    X <- PhyloEM_grid_alpha(phylo = phylo,
                            Y_data = Y_data,
                            process = process, 
                            independent = independent, 
                            K_max = K_max,
                            use_previous = use_previous, 
                            order = order, 
                            method.selection = method.selection, 
                            C.BM1 = C.BM1, 
                            C.BM2 = C.BM2, 
                            C.LINselect = C.LINselect, 
                            method.variance = method.variance, 
                            method.init = method.init, 
                            method.init.alpha = "default", 
                            method.init.alpha.estimation = method.init.alpha.estimation, 
                            methods.segmentation = methods.segmentation, 
                            alpha_known = TRUE, 
                            random.root = random.root, 
                            stationary.root = stationary.root, 
                            alpha = alpha, 
                            nbr_alpha = nbr_alpha,
                            check.tips.names = check.tips.names, 
                            progress.bar = progress.bar, 
                            estimates = estimates, 
                            save_step = save_step, 
                            # sBM_variance = sBM_variance, 
                            method.OUsun = method.OUsun, 
                            parallel_alpha = parallel_alpha, 
                            Ncores = Ncores, 
                            # exportFunctions = exportFunctions, 
                            # impute_init_Rphylopars = impute_init_Rphylopars, 
                            K_lag_init = K_lag_init,
                            light_result = light_result,
                            allow_negative = allow_negative,
                            trait_correlation_threshold = trait_correlation_threshold,
                            ...)
  } else { # For an in-loop estimation of alpha (independent = TRUE)
    if ((p > 1) && !independent){
      stop("Estimation of alpha outside of a grid is only implemented for independent traits. Please consider either use a grid for alpha values (parameter alpha_grid = TRUE), or independent traits (parameter independent = TRUE). See documentation.")
    }
    X <- PhyloEM_alpha_estim(phylo = phylo,
                             Y_data = Y_data,
                             process = process, 
                             independent = TRUE, 
                             K_max = K_max,
                             use_previous = use_previous, 
                             order = order, 
                             method.selection = method.selection, 
                             C.BM1 = C.BM1, 
                             C.BM2 = C.BM2, 
                             C.LINselect = C.LINselect, 
                             method.variance = method.variance, 
                             method.init = method.init, 
                             method.init.alpha = method.init.alpha, 
                             method.init.alpha.estimation = method.init.alpha.estimation, 
                             methods.segmentation = methods.segmentation, 
                             alpha_known = FALSE, 
                             random.root = random.root, 
                             stationary.root = stationary.root, 
                             alpha = NULL, 
                             check.tips.names = check.tips.names, 
                             progress.bar = progress.bar, 
                             estimates = estimates, 
                             save_step = save_step, 
                             method.OUsun = "raw", 
                             # impute_init_Rphylopars = impute_init_Rphylopars, 
                             K_lag_init = K_lag_init,
                             light_result = light_result,
                             trait_correlation_threshold = trait_correlation_threshold,
                             ...)
  }
  
  ## Return to original order if needed
  if (check_postorder){
    X <- return_to_original_order(X, phylo_original_order, phylo)
  }
  
  ## Save some paramaters
  X$phylo <- phylo_given
  X$p <- p
  X$process <- process
  X$times_shared <- compute_times_ca(X$phylo)
  X$distances_phylo <- compute_dist_phy(X$phylo)
  X$subtree.list <- enumerate_tips_under_edges(X$phylo)
  X$T_tree <- incidence.matrix(X$phylo)
  X$U_tree <- incidence.matrix.full(X$phylo)
  X$h_tree <- max(diag(as.matrix(X$times_shared))[1:ntaxa])
  X$Y_data <- Y_data
  X$light_result <- light_result
  
  ## Model Selection
  model_selection <- function(one.method.selection){
    mod_sel  <- switch(one.method.selection, 
                       BirgeMassart1 = model_selection_BM1,
                       BirgeMassart2 = model_selection_BM2,
                       BGHuni = model_selection_BGH,
                       BGHlsq = model_selection_BGH_leastsquares,
                       BGHml = model_selection_BGH_ml,
                       BGHmlraw = model_selection_BGH_mlraw,
                       BGHlsqraw = model_selection_BGH_leastsquares_raw,
                       pBIC = model_selection_pBIC,
                       pBIC_l1ou = model_selection_pBIC_l1ou)
    selection <- try(mod_sel(X, ntaxa = ncol(Y_data),
                             C.BM1 = C.BM1, C.BM2 = C.BM2, C.LINselect = C.LINselect,
                             tree = phylo_given, independent = independent,
                             T_tree = X$T_tree, times_shared = X$times_shared, 
                             distances_phylo = X$distances_phylo,
                             process = X$process, Y_data = X$Y_data))
    if (inherits(selection, "try-error")){
      warning(paste0("Model Selection ",  one.method.selection, " failled"))
    } else if (one.method.selection == "BGHlsq") {
      X$alpha_min <- selection
    } else if (one.method.selection == "BGHlsqraw") {
      X$alpha_min_raw <- selection
    } else {
      X$alpha_max <- selection
    }
    return(X)
  }
  ## Selection(s)
  for (meth.sel in method.selection){
    X <- model_selection(meth.sel)
  }
  
  ## Class and return
  class(X) <- "PhyloEM"
  return(X)
}

###############################################################################
## Handling Functions
###############################################################################

##
#' @title Parameter estimates
#'
#' @description
#' \code{params} takes an object of class \code{\link{PhyloEM}}, and returns the 
#' inferred parameters of the process.
#'
#' @param x an object of class \code{\link{PhyloEM}}
#' @param method.selection (optional) the method selection to be used.
#' One of "LINselect", "DDSE", "Djump". Default to "LINselect".
#' @param K (optional) an integer giving the number of shifts for which to retrieve
#' the parameters. Default to NULL (automatically selected number of shifts, see
#' \code{method.selection} argument).
#' @param alpha (optional) a value of alpha for which to retrieve the parameters. Can
#' be an (un-ambiguous) estimation of the true value. If
#' specified, then \code{K} must be precised too. Default to NULL (automatically
#' selected value, see \code{method.selection} argument).
#' @param rBM (optional) if TRUE, and if the process is "scOU", returns the raw
#' parameters of the BM on the re-scaled tree. Default to FALSE, except if 
#' the selection strength is negative (see doc of \code{\link{PhyloEM}} for
#' an explanation of this particular case).
#' @param init (optional) if TRUE, gives the parameters from the initialization of
#' the EM. Default to FALSE. This has no effect if \code{K} is not specified.
#' @param ... unused.
#' 
#' @return
#' An object of class \code{\link{params_process}}.
#' 
#' @seealso \code{\link{PhyloEM}}, \code{\link{imputed_traits.PhyloEM}}
#' 
# params_process(res)
# params_process(res, K = 3)
# params_process(res, K = 3, alpha = 0.33)
# params_process(res, K = 3, alpha = 3.12, rBM = TRUE)
# params_process(res, K = 3, alpha = 3.12, init = TRUE)
#' 
#' @export
#'
##
params_process.PhyloEM <- function(x, method.selection = NULL,
                                   K = NULL, alpha = NULL, rBM = FALSE,
                                   init = FALSE, ...){
  ## Select a given K
  if (is.null(K) && !is.null(alpha)) stop("If you specify alpha, you must also provide K.")
  if (!is.null(K)){
    if (!(K %in% x$K_try)){
      stop(paste0("The value of K: ", K, " was not found in the fitted object."))
    }
  ## Select a given alpha
    if (!is.null(alpha)){
      if (alpha == 0){
        tmp <- which("alpha_0" == names(x))
      } else {
        tmp <- grep(alpha, names(x))
      }
      if (length(tmp) == 0){
        stop(paste0("The value of alpha: ", alpha, " was not found in the fitted object."))
      } else if (length(tmp) > 1){
        stop(paste0("The value of alpha: ", alpha, " is ambiguous (several possible)."))
      } else {
        alpha_name <- names(x)[tmp]
      }
    } else {
      alpha_name <- "alpha_max"
    }
    if (init){
      res <- x[[alpha_name]]$params_init_estim[[paste0(K)]]
    } else {
      res <- x[[alpha_name]]$params_estim[[paste0(K)]] 
    }
  } else {
  ## Take the selected parameters (default)
    m_sel <- get_method_selection(x, method.selection = method.selection)
    res <- extract_params(x, m_sel[1], m_sel[2])
  }
  ## Case scOU with negative value
  if ((length(as.vector(res$selection.strength)) == 1)
        && (res$selection.strength < 0)
        && !rBM) {
    warning("The 'selection strength' is negative. One should only look at the un-normalized values of the shifts. To do so, please call this function using 'rBM = TRUE'.")
  }
  ## Return to rBM parameters if needed
  if (rBM){
    if (is.null(res$selection.strength) || sum(res$selection.strength) < .Machine$double.eps^(1/2)){
      ## Already a BM
      res <- res
      message("The process already was a BM.")
    } else if (length(as.vector(res$selection.strength)) > 1){
      ## Not an scOU
      stop("There is only an equivalent rBM process for an univariate OU or a multivariate scOU. Seems that you are using a full multivariate OU.")
    } else {
      res <- compute_raw_parameters(x$phylo, res) 
    }
  }
  ## Check dimensions
  tmp <- check_dimensions(x$p,
                          res$root.state,
                          res$shifts,
                          res$variance,
                          res$selection.strength,
                          res$optimal.value)
  res$root.state <- tmp$root.state
  res$shifts <- tmp$shifts
  res$variance <- tmp$variance
  res$selection.strength <- tmp$selection.strength
  res$optimal.value <- tmp$optimal.value
  ## Process
  res$process <- x$process
  if (rBM || (!is.null(alpha) && alpha == 0)) res$process <- "BM"
  ## Check root state
  res$root.state <- test.root.state(res$root.state,
                                    res$process,
                                    variance = res$variance,
                                    selection.strength = res$selection.strength,
                                    optimal.value = res$optimal.value)
  if (!is.null(rownames(x$Y_data))) res <- name_params(res, rownames(x$Y_data))
  res$variance <- as(res$variance, "dpoMatrix")
  class(res) <- "params_process"
  if (attr(res, "Neq") > 1){
    warning("There are several equivalent solutions for this shift position.")
  }
  return(res)
}

name_params <- function(res, names) {
  ## root state
  if (isnonnullna(res$root.state$value.root)) names(res$root.state$value.root) <- names
  if (isnonnullna(res$root.state$exp.root)) names(res$root.state$exp.root) <- names
  res$root.state$var.root <- name_matrix(res$root.state$var.root, names)
  ## shifts
  if (isnonnullna(res$shifts$values)) rownames(res$shifts$values) <- names
  ## variance
  res$variance <- name_matrix(res$variance, names)
  ## selection strength
  res$selection.strength <- name_matrix(res$selection.strength, names)
  ## optimal values
  if (isnonnullna(res$optimal.value)) names(res$optimal.value) <- names
  return(res)
}

isnonnullna <- function(x) {
  return(!is.null(x) && !any(is.na(x)))
}

name_matrix <- function(M, names) {
  if (isnonnullna(M)) {
    if (!is.null(attr(class(M), "package")) && attr(class(M), "package") == "Matrix") {
      dimnames(M) <- rep.int(list(names), 2L)
    } else {
      colnames(M) <- rownames(M) <- names
    }
  }
  return(M)
}

extract_params <- function(x, method, alpha_str){
  if (!is.null(x[[alpha_str]][[method]])){
    res <- x[[alpha_str]][[method]]$params_select
  } else {
    stop(paste0(method, "method was not used in the fit. Please use function 'model_selection' to add this criterion."))
  }
}

get_method_selection <- function(x, method.selection = NULL) {
  if (is.null(method.selection)){
    ## Take the selected parameters (default)
    if (x$p == 1){
      method.selection <- "BGHuni"
    } else {
      for (method in c("BGHml", "BGHlsq", "DDSE_BM1", "Djump_BM1",
                       "BGHmlraw", "pBIC")){
        if (method == "BGHlsq"){
          alpha_str <- "alpha_min"
        } else if (method == "BGHlsqraw"){
          alpha_str <- "alpha_min_raw"
        } else {
          alpha_str <- "alpha_max"
        }
        if (!is.null(x[[alpha_str]][[method]])){
          method.selection <- method
          break
        }
      }
      if (is.null(method.selection)) stop("No model selection procedure was found !")
    }
  } else {
    method.selection <- match.arg(method.selection,
                                  choices = c("LINselect", "DDSE", "Djump", "pBIC",
                                              "BGHlsq", "BGHml",
                                              "BGHlsqraw", "BGHmlraw",
                                              "BGH", "BGHuni",
                                              "log_likelihood")) 
    if (method.selection == "BGH") method.selection <- "BGHuni"
    if (method.selection == "LINselect") method.selection <- "BGHml"
  }
  # Result
  if (method.selection %in% c("DDSE", "DDSE_BM1")){
    res <- c("DDSE_BM1", "alpha_max", "DDSE", "max")
  }
  if (method.selection %in% c("Djump", "Djump_BM1")){
    res <- c("Djump_BM1", "alpha_max", "Djump", "max")
  }
  if (method.selection == "pBIC"){
    res <- c("pBIC", "alpha_max", "pBIC", "min")
  } 
  if (method.selection == "BGHuni"){
    res <- c("BGHuni", "alpha_max", "BGHuni", "max")
  } 
  if (method.selection == "BGHml"){
    res <- c("BGHml", "alpha_max", "BGHml", "min")
  } 
  if (method.selection == "BGHlsq"){
    res <- c("BGHlsq", "alpha_min", "BGHlsq", "min")
  } 
  if (method.selection == "BGHmlraw"){
    res <- c("BGHmlraw", "alpha_max", "BGHmlraw", "min")
  } 
  if (method.selection == "BGHlsqraw"){
    res <- c("BGHlsqraw", "alpha_min_raw", "BGHlsqraw", "min")
  }
  if (method.selection == "log_likelihood"){
    res <- c("log_likelihood", "alpha_max", "log likelihood", "no_max_min")
  }
  
  return(res)
}

##
#' @title Merge fits from independent runs of PhyloEM.
#'
#' @description
#' \code{merge_rotations} takes several fits from \code{\link{PhyloEM}}, and
#' merge them according to the best score (maximum likelihood or least squares).
#' For each number of shifts, 
#' The datasets needs to be equal up to a rotation. This is tested thanks to a QR
#' decomposition, see function \code{\link{find_rotation}}.
#'
#' @param ... objects of class \code{\link{PhyloEM}} fitted on datasets that are equal up to a rotation.
#' @param method.selection (optional) selection method to be applied to the merged fit. 
#' See \code{\link{params_process.PhyloEM}}.
#' @param tol (optional) relative numerical tolerance. See \code{\link{find_rotation}}.
#' 
#' @examples
#' \dontrun{
#' ## Load Data
#' data(monkeys)
#' ## Run method
#' # Note: use more alpha values for better results.
#' res <- PhyloEM(Y_data = monkeys$dat,        ## data
#'                phylo = monkeys$phy,         ## phylogeny
#'                process = "scOU",            ## scalar OU
#'                random.root = TRUE,          ## root is stationary
#'                stationary.root = TRUE,
#'                K_max = 10,                  ## maximal number of shifts
#'                nbr_alpha = 4,               ## number of alpha values
#'                parallel_alpha = TRUE,       ## parallelize on alpha values
#'                Ncores = 2)
#' ## Rotate dataset
#' rot <- matrix(c(cos(pi/4), -sin(pi/4), sin(pi/4), cos(pi/4)), nrow= 2, ncol = 2)
#' Yrot <- t(rot) %*% monkeys$dat
#' rownames(Yrot) <- rownames(monkeys$dat)
#' ## Fit rotated dataset
#' # Note: use more alpha values for better results.
#' res_rot <- PhyloEM(Y_data = Yrot,               ## rotated data
#'                    phylo = monkeys$phy,         
#'                    process = "scOU",            
#'                    random.root = TRUE,          
#'                    stationary.root = TRUE,
#'                    K_max = 10,                  
#'                    nbr_alpha = 4,               
#'                    parallel_alpha = TRUE,       
#'                    Ncores = 2)
#' ## Merge the two
#' res_merge <- merge_rotations(res, res_rot)
#' ## Plot the selected result
#' plot(res_merge)
#' ## Plot the model selection criterion
#' plot_criterion(res_merge)
#' }
#' 
#' @return
#' An object of class \code{\link{PhyloEM}}, result of the merge.
#' 
#' @export
#'
##
merge_rotations <- function(...,  method.selection = NULL, tol = NULL) {
  ress <- list(...)
  ## Basic tests
  if (length(ress) <= 1) stop("There should be at least 2 results to merge.")
  rots <- lapply(ress, function(x) find_rotation(ress[[1]], x, tol = tol))
  ## Criterion
  m_sels <- sapply(ress, get_method_selection)
  m_sel_1 <- m_sels[, 1]
  if (any(apply(m_sels, 2, function(x) x[1] != m_sel_1[1]))) {
    stop("The same criterion should be used for all the PhyloEM objects.")
  }
  ## Find maximum likelihood
  ll_max <- apply(sapply(ress, function(x) x$alpha_max$results_summary[["log_likelihood"]]), 1, which.max)
  lsq_min <- apply(sapply(ress, function(x) x$alpha_min$results_summary[["least_squares"]]), 1, which.min)
  lsq_min_raw <- apply(sapply(ress, function(x) x$alpha_min_raw$results_summary[["least_squares_raw"]]), 1, which.min)
  ## Make new result
  resMerge <- ress[[1]]
  resMerge[grep("alpha_[[:digit:]]", names(resMerge))] <- NULL
  for (K_t in 0:(length(ll_max)-1)) {
    resMerge <- merge_rotate(resMerge, ress, rots, "alpha_max", ll_max, K_t)
    resMerge <- merge_rotate(resMerge, ress, rots, "alpha_min", lsq_min, K_t)
    resMerge <- merge_rotate(resMerge, ress, rots, "alpha_min_raw", lsq_min_raw, K_t)
  }
  resMerge <- model_selection(resMerge, method.selection = m_sel_1)
  return(resMerge)
}

##
#' @title Test for rotation invariant datasets
#'
#' @description
#' \code{find_rotation} takes two fits from from \code{\link{PhyloEM}},
#' and test if their datasets are equal up to a rotation.
#'
#' @param res1 an object of class \code{\link{PhyloEM}}.
#' @param res2 an object of class \code{\link{PhyloEM}}.
#' @param tol relative numerical tolerance. Default to \code{.Machine$double.eps^(0.5)}.
#' 
#' 
#' @return
#' If appropriate, the rotation matrix rot such that dat1 = rot %*% dat2.
#' 
#' @export
#'
##
find_rotation <- function(res1, res2, tol = NULL) {
  dat1 <- t(res1$Y_data)
  dat2 <- t(res2$Y_data)
  # Basic checks
  if (any(dim(dat1) != dim(dat2))) stop("The datasets should have the same dimension.")
  # NAs
  if (anyNA(dat1)) {
    if (any(!(rowSums(is.na(dat1)) %in% c(0, ncol(dat1))))) stop("Rotations can only be applied to datasets that have entire species missing (i.e. entire columns in `Y_data`).")
    if (any(!(rowSums(is.na(dat2)) %in% c(0, ncol(dat2))))) stop("Rotations can only be applied to datasets that have entire species missing (i.e. entire columns in `Y_data`).")
    if (any(!(is.na(dat1) == is.na(dat2)))) stop("The two datasets used in the analyses do not have the same missing data.")
    dat1 <- na.omit(dat1)
    dat2 <- na.omit(dat2)
  }
  # Fit
  fit_12 <- lm.fit(dat1, dat2)
  # Linearly dependent ?
  if (is.null(tol)) tol <- .Machine$double.eps^(0.5)
  tolMax <- tol * sum(abs((dat1)))
  testQR <- (sum(abs(fit_12$residuals)) <= tolMax)
  if (!testQR) stop("The datasets are not linearly mapped.")
  # Rotation ?
  rot <- fit_12$coefficients
  testRot <- isTRUE(all.equal(as.vector(t(rot) %*% rot), c(1, 0, 0, 1)))
  if (!testRot) stop("The datasets are not linked by a rotation.")
  return(unname(rot))
}

merge_rotate <- function(resMerge, ress, rots, alpha_str, ll_mm, K_t) {
  res_mm <- ress[[ll_mm[K_t + 1]]][[alpha_str]]
  rot_mm <- rots[[ll_mm[K_t + 1]]]
  resMerge[[alpha_str]]$results_summary[K_t + 1, ] <- res_mm$results_summary[K_t + 1, ]
  resMerge[[alpha_str]]$params_estim[[paste(K_t)]] <- rotate_params(res_mm$params_estim[[paste(K_t)]], rot_mm)
  resMerge[[alpha_str]]$params_init_estim[[paste(K_t)]] <- rotate_params(res_mm$params_init_estim[[paste(K_t)]], rot_mm)
  resMerge[[alpha_str]]$edge.quality[[paste(K_t)]] <- res_mm$edge.quality[[paste(K_t)]]
  if (!resMerge$light_result){
    resMerge[[alpha_str]]$Yhat[[paste(K_t)]] <- rot_mm %*% res_mm$Yhat[[paste(K_t)]]
    resMerge[[alpha_str]]$Zhat[[paste(K_t)]] <- rot_mm %*% res_mm$Zhat[[paste(K_t)]]
    resMerge[[alpha_str]]$Yvar[[paste(K_t)]] <- res_mm$Yvar[[paste(K_t)]]
    resMerge[[alpha_str]]$Zvar[[paste(K_t)]] <- res_mm$Zvar[[paste(K_t)]]
    resMerge[[alpha_str]]$m_Y_estim[[paste(K_t)]] <- rot_mm %*% res_mm$m_Y_estim[[paste(K_t)]] 
  }
  return(resMerge)
}

rotate_params <- function(params, rot) {
  rot_params <- params
  # Vectors
  rot_params$shifts$values <- rot %*% params$shifts$values
  vec <- params$optimal.value
  if (!(is.null(vec) || (length(vec) == 1 && is.na(vec)))) {
    rot_params$optimal.value <- as.vector(rot %*% vec)
  }
  vec <- params$root.state$value.root
  if (!(is.null(vec) || (length(vec) == 1 && is.na(vec)))) {
    rot_params$root.state$value.root <- as.vector(rot %*% vec)
  }
  vec <- params$root.state$exp.root
  if (!(is.null(vec) || (length(vec) == 1 && is.na(vec)))) {
    rot_params$root.state$exp.root <- as.vector(rot %*% vec)
  }
  vec <- params$lambda
  if (!(is.null(vec) || (length(vec) == 1 && is.na(vec)))) {
    rot_params$lambda <- as.vector(rot %*% vec)
  }
  # Variances
  rot_params$variance <- forceSymmetric(rot %*% params$variance %*% t(rot))
  if (!(length(params$root.state$var.root) == 1 && is.na(params$root.state$var.root))) {
    rot_params$root.state$var.root <- forceSymmetric(rot %*% params$root.state$var.root %*% t(rot))
  }
  return(rot_params)
}

##
#' @title Ancestral State Reconstruction
#'
#' @description
#' \code{imputed_traits.PhyloEM} takes an object of class \code{\link{PhyloEM}},
#' and returns the imputed traits values, either at the internal nodes (ancestral
#' state reconstruction) or at the tips (data imputation)
#'
#' @param x an object of class \code{\link{PhyloEM}}.
#' @param trait an integer giving the trait to extract. Default to 1.
#' @param save_all if TRUE, arguments \code{where} and \code{what} are ignored, and
#' all the moments are kept for further extraction with the same function, specifying
#' the argument \code{reconstructed_states}. Default to FALSE.
#' @param where either "nodes" for ancestral state reconstruction, or "tips" for
#' data imputation.
#' @param what the quantity to retrieve. Either the imputed traits (default), their
#' conditional variances, or the simple expectations under the selected process.
#' @param params (optional) some user-specified parameters.
#' Must be of class \code{\link{params_process}}. If left blank, they are extracted
#' using the \code{method.selection} argument (see below).
#' @param method.selection (optional) the method selection to be used.
#' One of "LINselect", "DDSE", "Djump". Default to "LINselect".
#' @param reconstructed_states if the reconstructed states have already been
#' computed (by a previous call of the function, with \code{save_all=TRUE}),
#' they can be passed on here (avoids multiple computations of the E step).
#' @param ... further arguments to be passed on to
#' \code{\link{params_process.PhyloEM}}
#' 
#' 
#' @return
#' A matrix or array with the computed quantities.
#' 
#' @seealso \code{\link{params_process.PhyloEM}}, \code{\link{PhyloEM}}
#' 
# imputed_traits(res_new)
# imputed_traits(res_new, K = 10, alpha = 3)
#' 
#' @export
#'
##
imputed_traits <- function(x, ...) UseMethod("imputed_traits")

##
#' @describeIn imputed_traits \code{\link{PhyloEM}} object
#' @export
##
imputed_traits.PhyloEM <- function(x, trait = 1,
                                   save_all = FALSE,
                                   where = c("nodes", "tips"),
                                   what = c("imputed", "variances", "expectations"),
                                   params = NULL,
                                   method.selection = NULL,
                                   reconstructed_states = NULL,
                                   ...){
  ## Computes all the moments if needed
  if (is.null(reconstructed_states)){
    if (save_all) what <- c("imputed", "variances", "expectations")
    reconstructed_states <- compute_ancestral_traits(x, params, method.selection, what, ...)
  }
  
  ## Stop here if save_all=TRUE
  if (save_all) return(reconstructed_states)
  
  ## Else, extract the right moments
  where <- match.arg(where)
  what <- match.arg(what)
  if (where == "nodes"){
    if (what == "imputed"){
      res <- reconstructed_states$Zhat[trait, , drop = F]
    } else if (what == "variances"){
      res <- reconstructed_states$Zvar[trait, trait, , drop = F]
    } else if (what == "expectations"){
      res <- reconstructed_states$m_Z_estim[trait, , drop = F]
    }
  } else if (where == "tips"){
    if (what == "imputed"){
      res <- reconstructed_states$Yhat[trait, , drop = F]
    } else if (what == "variances"){
      res <- reconstructed_states$Yvar[trait, trait, , drop = F]
    } else if (what == "expectations"){
      res <- reconstructed_states$m_Y_estim[trait, , drop = F]
    }
  }
  
  return(res)
}

compute_ancestral_traits <- function(x,
                                     params,
                                     method.selection,
                                     what = c("imputed", "variances", "expectations"),
                                     ...){
  
  ## parameters
  if (is.null(params)){
    params <- params_process(x, method.selection, ...)
  } else {
    if (!inherits(params, "params_process")) {
      stop("The user specified parameters must be of class 'params_process'.")
    }
  }
  
  ## Heavy results
  if (!x$light_result && x$process == "scOU"){
    K <- length(params$shifts$edges)
    alpha <- unique(diag(params$selection.strength))
    tmp <- x[[paste0("alpha_", alpha)]]
    res <- list(
      m_Y_estim = tmp$m_Y_estim[[paste0(K)]],
      m_Z_estim = NULL,
      Zhat = tmp$Zhat[[paste0(K)]],
      Yhat = tmp$Yhat[[paste0(K)]],
      Zvar = tmp$Zvar[[paste0(K)]],
      Yvar = tmp$Yvar[[paste0(K)]]
    )
    return(res)
  }

  ## Needed quatities
  ntaxa <- length(x$phylo$tip.label)
  miss <- as.vector(is.na(x$Y_data))
  Y_data_vec_known <- as.vector(x$Y_data[!miss])
  masque_data <- rep(FALSE, (ntaxa + x$phylo$Nnode) * x$p)
  masque_data[1:(x$p * ntaxa)] <- !miss
  
  ## what to do
  what <- match.arg(what, several.ok = TRUE)
  res <- vector(mode = "list")
  
  ## Compute the expectations
  if ("expectations" %in% what){
    tmpsim <- simulate_internal(phylo = x$phylo, 
                                process = x$process,
                                p = x$p,
                                root.state = params$root.state, 
                                shifts = params$shifts, 
                                variance = params$variance, 
                                optimal.value = params$optimal.value, 
                                selection.strength = params$selection.strength,
                                simulate_random = FALSE,
                                U_tree = x$U_tree)
    res$m_Y_estim <- extract_simulate_internal(tmpsim,
                                               where = "tips",
                                               what = "expectations")
    res$m_Z_estim <- extract_simulate_internal(tmpsim,
                                               where = "nodes",
                                               what = "expectations")
  }
  
  ## Post order
  if ("imputed" %in% what || "variances" %in% what){
    phy <- reorder(x$phylo, "postorder")
    params$shifts$edges <- correspondanceEdges(edges = params$shifts$edges,
                                               from = x$phylo, to = phy)
    U_tree <- x$U_tree[, correspondanceEdges(edges = 1:nrow(phy$edge),
                                             from = phy, to = x$phylo)]
    
    ## E step
    temp <- wrapper_E_step(phylo = phy,
                           times_shared = NULL,
                           distances_phylo = NULL,
                           process = x$process,
                           params_old = params,
                           masque_data = masque_data,
                           F_moments = NULL,
                           independent = FALSE,
                           Y_data_vec_known = Y_data_vec_known,
                           miss = miss,
                           Y_data = x$Y_data,
                           U_tree = U_tree,
                           compute_E = compute_E.upward_downward)
    
    ## Checking for consistency
    if (!isTRUE(all.equal(as.vector(temp$log_likelihood_old),
                          attr(params, "log_likelihood"),
                          tol = .Machine$double.eps ^ 0.2))){
      warning(paste0("For K = ", length(params$shifts$edges), ", the log_likelihood of the transformed parameters on the er-scaled tree is different from the log_likelihood of the parameters on the original tree, with a tolerance of ", .Machine$double.eps ^ 0.2, "."))
      # stop("Something went wrong: log likelihood of supposedly equivalent parameters are not equal.")
    }
    res$Zhat <- temp$conditional_law_X$expectations[ , (ntaxa+1):ncol(temp$conditional_law_X$expectations), drop = FALSE]
    res$Yhat <- temp$conditional_law_X$expectations[ , 1:ntaxa, drop = FALSE]
    res$Zvar <- temp$conditional_law_X$variances[ , , (ntaxa+1):ncol(temp$conditional_law_X$expectations), drop = FALSE]
    res$Yvar <- temp$conditional_law_X$variances[ , , 1:ntaxa, drop = FALSE]
  }
  
  ## Result
  return(res)
}

##
#' @export
#' @method print PhyloEM
##
print.PhyloEM <- function(x, ...){
  cat("Result of the PhyloEM algorithm.\n")
  cat("Selected parameters by the default method:")
  print(params_process.PhyloEM(x))
  cat("\n\nSee help to see all plotting and handling functions.")
}

##
#' @title Make the result lighter
#'
#' @description
#' \code{enlight.PhyloEM} takes an object of class \code{\link{PhyloEM}},
#' and returns the same object, without saving the quantities that can be easily
#' re-computed using function \code{\link{imputed_traits.PhyloEM}}.
#' 
#' @details 
#' The resulting object can be much lighter, saving a lot of memory space, but each
#' call to the function \code{\link{imputed_traits.PhyloEM}} will be longer. As 
#' function \code{\link{plot.PhyloEM}} relies on this function, this makes the 
#' plotting also longer.
#' This has the same effect as setting the option "\code{light_result=TRUE}" in the
#' call of \code{\link{PhyloEM}}.
#'
#' @param x an object of class \code{\link{PhyloEM}}.
#' 
#' @return
#' Same as entry, lighter.
#' 
#' @seealso \code{\link{PhyloEM}}, \code{\link{imputed_traits.PhyloEM}},
#' \code{\link{plot.PhyloEM}}
#' 
#' @export
#'
##
enlight <- function(x) UseMethod("enlight")

##
#' @describeIn enlight \code{\link{PhyloEM}} object
#' @export
##
enlight.PhyloEM <- function(x){
  if (x$light_result){
    message("The object was already in light format.")
    return(x)
  }
  ## delete imputed quantities
  lres <- x
  put_imput_null <- function(y){
    y$Zhat <- NULL
    y$Yhat <- NULL
    y$Zvar <- NULL
    y$Yvar <- NULL
    y$m_Y_estim <- NULL
    return(y)
  }
  inds <- grep("alpha_", names(lres))
  lres[inds] <- lapply(lres[inds], put_imput_null)
  ## Put light_result to FALSE
  lres$light_result <- TRUE
  return(lres)
}

###############################################################################
## PhyloEM_grid_alpha
###############################################################################

PhyloEM_grid_alpha <- function(phylo, Y_data, process = c("BM", "OU", "scOU", "rBM"),
                               independent = FALSE,
                               K_max, use_previous = TRUE,
                               order = TRUE,
                               method.selection = c("BirgeMassart1", "BirgeMassart2",
                                                    "BGHuni", "pBIC", "pBIC_l1ou",
                                                    "BGHlsq", "BGHml",
                                                    "BGHlsqraw", "BGHmlraw"),
                               C.BM1 = 0.1, C.BM2 = 2.5, C.LINselect = 1.1,
                               method.variance = "simple",
                               method.init = "default",
                               method.init.alpha = "default",
                               method.init.alpha.estimation = c("regression",
                                                                "regression.MM",
                                                                "median"), 
                               methods.segmentation = c("lasso", "best_single_move"),
                               alpha_known = TRUE,
                               random.root = FALSE,
                               stationary.root = FALSE,
                               alpha = NULL,
                               nbr_alpha = nbr_alpha,
                               check.tips.names = FALSE,
                               progress.bar = TRUE,
                               estimates = NULL,
                               save_step = FALSE,
                               sBM_variance = FALSE,
                               method.OUsun = "rescale",
                               parallel_alpha = FALSE,
                               Ncores = 3,
                               # exportFunctions = ls(),
                               # impute_init_Rphylopars = FALSE,
                               K_lag_init = 0,
                               light_result = TRUE,
                               allow_negative = FALSE,
                               trait_correlation_threshold = 0.9,
                               ...){
  # reqpckg <- c("ape", "glmnet", "robustbase")
  reqpckg <- c("PhylogeneticEM")
  ntaxa <- length(phylo$tip.label)
  ## Save Original process and tree
  process <- match.arg(process)
  process_original <- process
  original_phy <- phylo
  
  ## Fixed Quantities
  times_shared_original <- compute_times_ca(phylo)
  distances_phylo_original <- compute_dist_phy(phylo)
  subtree.list_original <- enumerate_tips_under_edges(phylo)
  h_tree_original <- max(diag(as.matrix(times_shared_original))[1:ntaxa])
  T_tree = incidence.matrix(phylo)
  U_tree = incidence.matrix.full(phylo)
  
  ## Missing informations
  miss <- as.vector(is.na(Y_data))
  Y_data_vec <- as.vector(Y_data)
  Y_data_vec_known <- as.vector(Y_data[!miss])
  # Vectorized Data Mask
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  p <- nrow(Y_data)
  masque_data <- rep(FALSE, (ntaxa + Nnode) * p)
  masque_data[1:(p*ntaxa)] <- !miss
  
  ## Compute alpha
  if (process == "BM") {
    alpha <- 0
  } else {
    alpha <- find_grid_alpha(phylo, alpha, nbr_alpha = nbr_alpha, allow_negative = allow_negative, ...)
    if (stationary.root) alpha <- alpha[alpha != 0]
  }
  ## Check alpha for numerical instabilities
  check_range_alpha(alpha, h_tree_original)
  ## Loop on alpha
  estimate_alpha_several_K <- function(alp, 
                                       original_phy, Y_data,
                                       process_original,
                                       process,
                                       independent,
                                       K_max, 
                                       use_previous,
                                       order,
                                       method.variance,
                                       method.init,
                                       method.init.alpha,
                                       method.init.alpha.estimation, 
                                       methods.segmentation,
                                       alpha_known,
                                       random.root,
                                       stationary.root,
                                       sBM_variance,
                                       method.OUsun,
                                       # impute_init_Rphylopars,
                                       p,
                                       ntaxa,
                                       progress.bar,
                                       times_shared_original,
                                       distances_phylo_original,
                                       subtree.list_original,
                                       h_tree_original,
                                       T_tree,
                                       U_tree,
                                       K_lag_init,
                                       light_result,
                                       allow_negative,
                                       trait_correlation_threshold,
                                       ...){
    if(progress.bar){
      message(paste0("Alpha ", alp))
    }
    temp <- choose_process_EM(process = process_original,
                              p = p,
                              random.root = random.root,
                              stationary.root = stationary.root,
                              alpha_known = alpha_known,
                              known.selection.strength = alp,
                              sBM_variance = sBM_variance,
                              method.OUsun = method.OUsun,
                              independent = independent,
                              allow_negative = allow_negative)
    
    rescale_tree <- temp$rescale_tree # Rescale the tree ?
    transform_scOU <- temp$transform_scOU # Re-transform parameters back ?
    sBM_variance <- temp$sBM_variance
    ## Transform branch lengths if needed
    phylo <- original_phy
    if (sBM_variance){ # process sBM : need a root branch.
      phylo$root.edge <- 1
    }
    if (rescale_tree) {
      phylo <- transform_branch_length(phylo, alp)
    }
    ## Fixed quantities
    if (rescale_tree){
      times_shared <- compute_times_ca(phylo)
      distances_phylo <- compute_dist_phy(phylo)
      subtree.list <- enumerate_tips_under_edges(phylo)
      h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
    } else {
      times_shared <- times_shared_original
      distances_phylo <- distances_phylo_original
      subtree.list <- subtree.list_original
      h_tree <- h_tree_original
    }
    ## Fixed Quantities if no missing data
    Flag_Missing <- any(is.na(Y_data)) # TRUE if some missing values
    if ((!Flag_Missing) && (method.variance != "upward_downward")){
      # Add root edge to the branch lengths (root assumed fixed by default)
      root_edge_length <- 0
      if (!is.null(phylo$root.edge)) root_edge_length <- phylo$root.edge
      F_moments <- compute_fixed_moments(times_shared + root_edge_length, ntaxa)
    } else {
      F_moments = NULL
    }
    ## Impute data if needed
    # if (!order && !impute_init_Rphylopars && any(is.na(Y_data))){
    #   warning("There are some missing values, and the inference is not done by increasing values of shifts, so they cannot be inferred. Using Rphylopars for the initialization (impute_init_Rphylopars = TRUE)")
    #   impute_init_Rphylopars <- TRUE
    # }
    Y_data_imp <- Y_data
    # if (any(is.na(Y_data_imp))
    #     && impute_init_Rphylopars
    #     && temp$process == "BM"){
    #   ## Re-scale tree to unit height
    #   factor_rescale <- 1 / h_tree # total height to 1
    #   phylo_temp <- phylo
    #   phylo_temp$edge.length <- factor_rescale * phylo$edge.length
    #   phylo_temp$root.edge <- factor_rescale * phylo$root.edge
    #   Y_data_imp <- try(impute.data.Rphylopars(phylo_temp,
    #                                            Y_data,
    #                                            temp$process,
    #                                            random.root))
    #   if (inherits(Y_data_imp, "try-error")) { # If fails, replace with mean of the trait
    #     Y_data_imp <- Y_data
    #     for (j in 1:(dim(Y_data_imp)[1])){
    #       Y_data_imp[j, is.na(Y_data_imp[j, ])] <- mean(Y_data_imp[j, ], na.rm = TRUE)
    #     }
    #   }
    #   rm(phylo_temp)
    # }
    ## Estimations
    X <- Phylo_EM_sequencial(phylo = phylo,
                             Y_data = Y_data,
                             Y_data_imp = Y_data_imp,
                             process = temp$process,
                             independent = independent,
                             K_max = K_max,
                             # curent = X,
                             use_previous = use_previous,
                             order = order,
                             method.variance = method.variance,
                             method.init = method.init,
                             method.init.alpha = method.init.alpha,
                             method.init.alpha.estimation = method.init.alpha.estimation, 
                             methods.segmentation = methods.segmentation,
                             alpha_known = alpha_known,
                             random.root = random.root,
                             stationary.root = stationary.root,
                             alp = alp,
                             check.tips.names = check.tips.names,
                             progress.bar = progress.bar,
                             times_shared = times_shared,
                             distances_phylo = distances_phylo,
                             subtree.list = subtree.list,
                             T_tree = T_tree,
                             U_tree = U_tree,
                             h_tree = h_tree,
                             F_moments = F_moments,
                             save_step = save_step,
                             sBM_variance = sBM_variance,
                             method.OUsun = method.OUsun,
                             # impute_init_Rphylopars = impute_init_Rphylopars,
                             K_lag_init = K_lag_init,
                             light_result = light_result,
                             allow_negative = allow_negative,
                             trait_correlation_threshold = trait_correlation_threshold,
                             ...)
    ## Trnasform back parameters to OU if needed
    if (transform_scOU){
      ## Compute equivalent parameters
      fun1 <- function(params){
        params_scOU <- go_back_to_original_process(phy_original = original_phy,
                                                   known.selection.strength = alp,
                                                   sBM_variance = sBM_variance,
                                                   params = params)
      }
      X$params_estim <- lapply(X$params_estim, fun1)
      rm(fun1)
      ## Normalize least squares
      X$results_summary$least_squares <- X$results_summary$least_squares / (2 * abs(alp))
      ## Ancestral state reconstruction
      if (!light_result){
        compute_E  <- switch(method.variance,
                             simple = compute_E.simple,
                             upward_downward = compute_E.upward_downward)
        fun2 <- function(params_scOU){
          temp <- wrapper_E_step(phylo = original_phy,
                                 times_shared = times_shared_original,
                                 distances_phylo = distances_phylo_original,
                                 process = process_original,
                                 params_old = params_scOU,
                                 masque_data = masque_data,
                                 F_moments = NULL,
                                 independent = FALSE,
                                 Y_data_vec_known = Y_data_vec_known,
                                 miss = miss,
                                 Y_data = Y_data,
                                 U_tree = U_tree,
                                 # compute_mean_variance = compute_mean_variance.simple,
                                 # compute_log_likelihood = compute_log_likelihood.simple,
                                 # compute_mahalanobis_distance = compute_mahalanobis_distance.simple,
                                 compute_E = compute_E)
          if (!isTRUE(all.equal(as.vector(temp$log_likelihood_old),
                                attr(params_scOU, "log_likelihood"),
                                tol = .Machine$double.eps ^ 0.2))){
            warning(paste0("For K = ", length(params_scOU$shifts$edges), ", the log_likelihood of the transformed parameters on the er-scaled tree is different from the log_likelihood of the parameters on the original tree, with a tolerance of ", .Machine$double.eps ^ 0.2, "."))
            # stop("Something went wrong: log likelihood of supposedly equivalent parameters are not equal.")
          }
          tmpsim <- simulate_internal(phylo = original_phy,
                                      process = process_original,
                                      p = attr(params_scOU, "p_dim"),
                                      root.state = params_scOU$root.state,
                                      shifts = params_scOU$shifts,
                                      variance = params_scOU$variance,
                                      optimal.value = params_scOU$optimal.value,
                                      selection.strength = params_scOU$selection.strength,
                                      simulate_random = FALSE,
                                      U_tree = U_tree)
          return(list(Zhat = temp$conditional_law_X$expectations[ , (ntaxa+1):ncol(temp$conditional_law_X$expectations)],
                      Yhat = temp$conditional_law_X$expectations[ , 1:ntaxa],
                      Zvar = temp$conditional_law_X$variances[ , , (ntaxa+1):ncol(temp$conditional_law_X$expectations)],
                      Yvar = temp$conditional_law_X$variances[ , , 1:ntaxa],
                      m_Y_estim = extract_simulate_internal(tmpsim,
                                                            where="tips",
                                                            what="expectations")))
        }
        temp_list <- lapply(X$params_estim, fun2)
        rm(fun2)
        ## "Ancestral States Reconstruction"
        X$Zhat <- lapply(temp_list, function(z) z$Zhat)
        X$Yhat <- lapply(temp_list, function(z) z$Yhat)
        X$Zvar <- lapply(temp_list, function(z) z$Zvar)
        X$Yvar <- lapply(temp_list, function(z) z$Yvar)
        X$m_Y_estim <- lapply(temp_list, function(z) z$m_Y_estim)
        rm(temp_list) 
      }
    }
    return(X)
  }
  a_greek <- NULL
  if (parallel_alpha){
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is needed for parallel computation (option 'parallel_alpha = TRUE'). Please install this package, or set the option to 'FALSE'.",
           call. = FALSE)
    }
    cl <- parallel::makeCluster(Ncores, outfile = "")
                                # outfile = tempfile(pattern = "log_file_dopar_"))
    doParallel::registerDoParallel(cl)
    X <- foreach::foreach(a_greek = alpha, .packages = reqpckg) %dopar%
    {
      estimate_alpha_several_K(alp = a_greek,
                               original_phy = original_phy, Y_data = Y_data,
                               process_original = process_original,
                               process = process,
                               independent = independent,
                               K_max = K_max, 
                               use_previous = use_previous,
                               order = order,
                               method.variance = method.variance,
                               method.init = method.init,
                               method.init.alpha = method.init.alpha,
                               method.init.alpha.estimation = method.init.alpha.estimation, 
                               methods.segmentation = methods.segmentation,
                               alpha_known = alpha_known,
                               random.root = random.root,
                               stationary.root = stationary.root,
                               sBM_variance = sBM_variance,
                               method.OUsun = method.OUsun,
                               # impute_init_Rphylopars = impute_init_Rphylopars,
                               p = p,
                               ntaxa = ntaxa,
                               progress.bar = progress.bar,
                               times_shared_original = times_shared_original,
                               distances_phylo_original = distances_phylo_original,
                               subtree.list_original = subtree.list_original,
                               h_tree_original = h_tree_original,
                               T_tree = T_tree,
                               U_tree = U_tree,
                               K_lag_init = K_lag_init,
                               light_result = light_result,
                               allow_negative = allow_negative,
                               trait_correlation_threshold = trait_correlation_threshold,
                               ...)
    }
    parallel::stopCluster(cl)
  } else {
    X <- foreach::foreach(a_greek = alpha, .packages = reqpckg) %do%
    {
      estimate_alpha_several_K(alp = a_greek,
                               original_phy = original_phy, Y_data = Y_data,
                               process_original = process_original,
                               process = process,
                               independent = independent,
                               K_max = K_max, 
                               use_previous = use_previous,
                               order = order,
                               method.variance = method.variance,
                               method.init = method.init,
                               method.init.alpha = method.init.alpha,
                               method.init.alpha.estimation = method.init.alpha.estimation, 
                               methods.segmentation = methods.segmentation,
                               alpha_known = alpha_known,
                               random.root = random.root,
                               stationary.root = stationary.root,
                               sBM_variance = sBM_variance,
                               method.OUsun = method.OUsun,
                               # impute_init_Rphylopars = impute_init_Rphylopars,
                               p = p,
                               ntaxa = ntaxa,
                               progress.bar = progress.bar,
                               times_shared_original = times_shared_original,
                               distances_phylo_original = distances_phylo_original,
                               subtree.list_original = subtree.list_original,
                               h_tree_original = h_tree_original,
                               T_tree = T_tree,
                               U_tree = U_tree,
                               K_lag_init = K_lag_init,
                               light_result = light_result,
                               allow_negative = allow_negative,
                               trait_correlation_threshold = trait_correlation_threshold,
                               ...)
    }
  }
  
  ## Format Output
  names(X) <- paste0("alpha_", alpha)
  X$Y_data <- Y_data
  X$K_try <- 0:K_max
  X$ntaxa <- ntaxa

  ## Select max solution for each K
  X <- merge_max_grid_alpha(X, alpha, light_result)
  # if ("BGHlsq" %in% method.selection){
  X <- merge_min_grid_alpha(X, light_result) 
  # }
  # if ("BGHlsqraw" %in% method.selection){
  X <- merge_min_grid_alpha(X, light_result, raw = TRUE) 
  # }
  return(X)
}

###############################################################################
## PhyloEM_alpha_estim
###############################################################################

PhyloEM_alpha_estim <- function(phylo, Y_data, process = c("BM", "OU", "scOU", "rBM"),
                                independent = TRUE,
                                K_max, use_previous = TRUE,
                                order = TRUE,
                                method.selection = c("BirgeMassart1", "BirgeMassart2", "BGHuni", "pBIC", "pBIC_l1ou", "BGHlsq", "BGHml", "BGHlsqraw", "BGHmlraw"),
                                C.BM1 = 0.1, C.BM2 = 2.5, C.LINselect = 1.1,
                                method.variance = c("simple", "upward_downward"),
                                method.init = "default",
                                method.init.alpha = "default",
                                method.init.alpha.estimation = c("regression", "regression.MM", "median"), 
                                methods.segmentation = c("lasso", "best_single_move"),
                                alpha_known = FALSE,
                                random.root = FALSE,
                                stationary.root = FALSE,
                                alpha = NULL,
                                check.tips.names = FALSE,
                                progress.bar = TRUE,
                                estimates = NULL,
                                save_step = FALSE,
                                method.OUsun = "raw",
                                # impute_init_Rphylopars = FALSE,
                                K_lag_init = 0,
                                light_result = TRUE,
                                trait_correlation_threshold = 0.9,
                                ...){
  ## Fixed quantities
  ntaxa <- length(phylo$tip.label)
  times_shared <- compute_times_ca(phylo)
  distances_phylo <- compute_dist_phy(phylo)
  subtree.list <- enumerate_tips_under_edges(phylo)
  T_tree <- incidence.matrix(phylo)
  U_tree <- incidence.matrix.full(phylo)
  h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
  
  Y_data_imp <- Y_data
  #     if (any(is.na(Y_data_imp))
  #         && impute_init_Rphylopars
  #         && temp$process == "BM"){
  #       ## Re-scale tree to unit height
  #       factor_rescale <- 1 / h_tree # total height to 1
  #       phylo_temp <- phylo
  #       phylo_temp$edge.length <- factor_rescale * phylo$edge.length
  #       phylo_temp$root.edge <- factor_rescale * phylo$root.edge
  #       Y_data_imp <- try(impute.data.Rphylopars(phylo_temp,
  #                                                Y_data,
  #                                                temp$process,
  #                                                random.init))
  #       if (inherits(Y_data_imp, "try-error")) { # If fails, replace with mean of the trait
  #         Y_data_imp <- Y_data
  #         for (j in 1:(dim(Y_data_imp)[1])){
  #           Y_data_imp[j, is.na(Y_data_imp[j, ])] <- mean(Y_data_imp[j, ], na.rm = TRUE)
  #         }
  #       }
  #       rm(phylo_temp)
  #     }
  ## Estimations
  X <- Phylo_EM_sequencial(phylo = phylo,
                           Y_data = Y_data,
                           Y_data_imp = Y_data_imp,
                           process = process,
                           independent = independent,
                           K_max = K_max,
                           # curent = X,
                           use_previous = use_previous,
                           order = order,
                           method.variance = method.variance,
                           method.init = method.init,
                           method.init.alpha = method.init.alpha,
                           method.init.alpha.estimation = method.init.alpha.estimation, 
                           methods.segmentation = methods.segmentation,
                           alpha_known = alpha_known,
                           random.root = random.root,
                           stationary.root = stationary.root,
                           alp = NULL,
                           check.tips.names = check.tips.names,
                           progress.bar = progress.bar,
                           times_shared = times_shared,
                           distances_phylo = distances_phylo,
                           subtree.list = subtree.list,
                           T_tree = T_tree,
                           U_tree = U_tree,
                           h_tree = h_tree,
                           save_step = save_step,
                           method.OUsun = method.OUsun,
                           # impute_init_Rphylopars = impute_init_Rphylopars,
                           K_lag_init = K_lag_init,
                           light_result = light_result,
                           trait_correlation_threshold = trait_correlation_threshold,
                           ...)
  
  ## Format Output
  X2 <- vector("list", 4)
  names(X2) <- c("alpha_max", "Y_data", "K_try", "ntaxa")
  X2$alpha_max <- X
  X2$Y_data <- Y_data
  X2$K_try <- 0:K_max
  X2$ntaxa <- ntaxa
  
  return(X2)
}

###############################################################################
## Phylo_EM_sequencial
###############################################################################

Phylo_EM_sequencial <- function(phylo, Y_data,
                                Y_data_imp,
                                process,
                                independent,
                                K_max,
                                #                               curent = list(Y_data = Y_data,
                                #                                             K_try = 0:K_max,
                                #                                             ntaxa = length(phylo$tip.label)),
                                use_previous = TRUE,
                                order = TRUE,
                                method.variance = c("simple", "upward_downward"),
                                method.init = "default",
                                method.init.alpha = "default",
                                method.init.alpha.estimation = c("regression",
                                                                 "regression.MM",
                                                                 "median"), 
                                methods.segmentation = c("lasso", "best_single_move"),
                                alpha_known = FALSE,
                                random.root = TRUE,
                                stationary.root = TRUE,
                                alp = NULL,
                                check.tips.names = FALSE,
                                progress.bar = TRUE,
                                times_shared = NULL,
                                distances_phylo = NULL,
                                subtree.list = NULL,
                                T_tree = NULL,
                                U_tree = NULL,
                                h_tree = NULL,
                                F_moments = NULL,
                                save_step = FALSE,
                                sBM_variance = FALSE,
                                method.OUsun = "rescale", 
                                # impute_init_Rphylopars = FALSE,
                                K_lag_init = 0,
                                light_result = TRUE,
                                allow_negative = FALSE,
                                trait_correlation_threshold = trait_correlation_threshold,
                                ...){
  p <- nrow(Y_data)
  ntaxa <- length(phylo$tip.label)
  ## First estim
  if (order){
    K_first <- 0
    K_last <- K_max
    next_it <- function(K_t) { return(K_t + 1) }
    prev_it <- function(K_t) { return(K_t - 1) }
  } else {
    K_first <- K_max
    K_last <- 0
    next_it <- function(K_t) { return(K_t - 1) }
    prev_it <- function(K_t) { return(K_t + 1) }
  }
  ## Set up
  XX <- vector('list', K_max + 1)
  names(XX) <- 0:K_max
  ## Progress Bar
  if(progress.bar){
    pb <- txtProgressBar(min = 0, max = K_max + 1, style = 3)
  }
  ## First
  XX[[paste0(K_first)]] <- estimateEM_wrapper_scratch(phylo = phylo,
                                                      Y_data = Y_data,
                                                      Y_data_imp = Y_data_imp,
                                                      process = process,
                                                      independent = independent,
                                                      K_t = K_first,
                                                      method.variance = method.variance,
                                                      random.root = random.root,
                                                      stationary.root = stationary.root,
                                                      alpha_known = alpha_known,
                                                      alpha = alp,
                                                      method.init = method.init,
                                                      method.init.alpha = method.init.alpha,
                                                      methods.segmentation = methods.segmentation,
                                                      times_shared = times_shared, 
                                                      distances_phylo = distances_phylo,
                                                      subtree.list = subtree.list,
                                                      T_tree = T_tree,
                                                      U_tree = U_tree,
                                                      h_tree = h_tree,
                                                      F_moments = F_moments,
                                                      warning_several_solutions = FALSE,
                                                      sBM_variance = sBM_variance,
                                                      method.OUsun = method.OUsun,
                                                      # impute_init_Rphylopars = impute_init_Rphylopars,
                                                      K_lag_init = K_lag_init,
                                                      allow_negative = allow_negative,
                                                      trait_correlation_threshold = trait_correlation_threshold,
                                                      ...)
  if (K_first == 0 && any(is.na(Y_data))){
    Y_data_imp <- XX[["0"]]$Yhat
  }
  pp <- check_dimensions(p,
                         XX[[paste0(K_first)]]$params$root.state,
                         XX[[paste0(K_first)]]$params$shifts,
                         XX[[paste0(K_first)]]$params$variance,
                         XX[[paste0(K_first)]]$params$selection.strength,
                         XX[[paste0(K_first)]]$params$optimal.value)
  XX[[paste0(K_first)]]$params$root.state <- pp$root.state
  XX[[paste0(K_first)]]$params$shifts <- pp$shifts
  XX[[paste0(K_first)]]$params$variance <- pp$variance
  XX[[paste0(K_first)]]$params$selection.strength <- pp$selection.strength
  XX[[paste0(K_first)]]$params$optimal.value <- pp$optimal.value
  # update progress bar
  if(progress.bar) setTxtProgressBar(pb, 1); counter <- 2;
  ## Iterations
  if (K_first != K_last){
    for (K_t in (next_it(K_first)):(K_last)){
      XX[[paste0(K_t)]] <- estimateEM_wrapper(use_previous)(phylo = phylo,
                                                            Y_data = Y_data,
                                                            Y_data_imp = Y_data_imp,
                                                            process = process,
                                                            independent = independent,
                                                            K_t = K_t,
                                                            prev = XX[[paste0(prev_it(K_t))]],
                                                            method.variance = method.variance,
                                                            random.root = random.root,
                                                            stationary.root = stationary.root,
                                                            alpha_known = alpha_known,
                                                            alpha = alp,
                                                            method.init = method.init,
                                                            method.init.alpha = method.init.alpha,
                                                            methods.segmentation = methods.segmentation,
                                                            times_shared = times_shared, 
                                                            distances_phylo = distances_phylo,
                                                            subtree.list = subtree.list,
                                                            T_tree = T_tree, 
                                                            U_tree = U_tree,
                                                            h_tree = h_tree,
                                                            F_moments = F_moments,
                                                            warning_several_solutions = FALSE,
                                                            sBM_variance = sBM_variance,
                                                            method.OUsun = method.OUsun,
                                                            # impute_init_Rphylopars = impute_init_Rphylopars,
                                                            K_lag_init = K_lag_init,
                                                            allow_negative = allow_negative,
                                                            trait_correlation_threshold = trait_correlation_threshold,
                                                            ...)
      pp <- check_dimensions(p,
                             XX[[paste0(K_t)]]$params$root.state,
                             XX[[paste0(K_t)]]$params$shifts,
                             XX[[paste0(K_t)]]$params$variance,
                             XX[[paste0(K_t)]]$params$selection.strength,
                             XX[[paste0(K_t)]]$params$optimal.value)
      XX[[paste0(K_t)]]$params$root.state <- pp$root.state
      XX[[paste0(K_t)]]$params$shifts <- pp$shifts
      XX[[paste0(K_t)]]$params$variance <- pp$variance
      XX[[paste0(K_t)]]$params$selection.strength <- pp$selection.strength
      XX[[paste0(K_t)]]$params$optimal.value <- pp$optimal.value
      if(progress.bar) setTxtProgressBar(pb, counter); counter <- counter + 1
    }
  }
  ## Format results and return
  res <- format_output_several_K_single(XX, light_result)
  if (save_step) save(res,
                      file = tempfile(pattern = paste0("Alpha=", alp, "_"),
                                      fileext = c(".RData")))
  return(res)
}

###############################################################################
## merge_max_grid_alpha
###############################################################################

merge_max_grid_alpha <- function(X, alpha, light_result = TRUE){
  summary_all <- X[[paste0("alpha_", alpha[1])]]$results_summary
  for (alp in alpha[-1]){
    summary_all <- rbind(summary_all,
                         X[[paste0("alpha_", alp)]]$results_summary)
  }
  summary_all$alpha_name <- rep(alpha, each = length(X$K_try))
  X$alpha_max$results_summary <- matrix(NA, nrow = length(X$K_try),
                                        ncol = ncol(summary_all))
  colnames(X$alpha_max$results_summary) <- colnames(summary_all)
  X$alpha_max$edge.quality <- vector(length = length(X$K_try), mode = "list")
  for (K_t in X$K_try){
    max_sum <- summary_all[summary_all$K_try == K_t, ]
    max_sum <- max_sum[max_sum$log_likelihood == max(max_sum$log_likelihood), ]
# subset(subset(summary_all, K_try == K_t), log_likelihood == max(log_likelihood))
    res_max <- X[[paste0("alpha_", max_sum$alpha_name)]]
    params <- res_max$params_estim[[paste(K_t)]]
    X$alpha_max$results_summary[K_t + 1, ] <- as.vector(unname(as.matrix(max_sum)))
    X$alpha_max$params_estim[[paste(K_t)]] <- params
    # X$alpha_max$params_raw[[paste(K_t)]] <- res_max$params_raw[[paste(K_t)]]
    X$alpha_max$params_init_estim[[paste(K_t)]] <- res_max$params_init_estim[[paste(K_t)]]
    X$alpha_max$edge.quality[[paste(K_t)]] <- res_max$edge.quality[[paste(K_t)]]
    if (!light_result){
      X$alpha_max$Yhat[[paste(K_t)]] <- res_max$Yhat[[paste(K_t)]]
      X$alpha_max$Zhat[[paste(K_t)]] <- res_max$Zhat[[paste(K_t)]]
      X$alpha_max$Yvar[[paste(K_t)]] <- res_max$Yvar[[paste(K_t)]]
      X$alpha_max$Zvar[[paste(K_t)]] <- res_max$Zvar[[paste(K_t)]]
      X$alpha_max$m_Y_estim[[paste(K_t)]] <- res_max$m_Y_estim[[paste(K_t)]]
    }
  }
  X$alpha_max$results_summary <- as.data.frame(X$alpha_max$results_summary)
  return(X)
}

# add_lsq <- function(X){
#   nums <- grep("alpha_-?[[:digit:]]", names(X))
#   for (i in nums){
#     X[[i]]$results_summary$least_squares <- sapply(X[[i]]$m_Y_estim,
#                                                    function(z) sum((X$Y_data - z)^2))
#     # X[[i]]$results_summary$least_squares <- sapply(X[[i]]$params_estim, function(z) sum(diag(z$variance)))
#   }
#   return(X)
# }

merge_min_grid_alpha <- function(X, light_result = TRUE, raw = FALSE){
  nums <- grep("alpha_-?[[:digit:]]", names(X))
  summary_all <- X[[nums[1]]]$results_summary
  for (i in nums[-1]){
    summary_all <- rbind(summary_all,
                         X[[i]]$results_summary)
  }
  summary_all$alpha_name <- rep(as.numeric(sub("alpha_", "", names(X)[nums])),
                                each = length(X$K_try))
  if (raw){
    lsq_str <- "least_squares_raw"
    alpha_min_str <- "alpha_min_raw"
  } else {
    lsq_str <- "least_squares"
    alpha_min_str <- "alpha_min"
  }
  X[[alpha_min_str]]$results_summary <- matrix(NA, nrow = length(X$K_try),
                                        ncol = ncol(summary_all))
  colnames(X[[alpha_min_str]]$results_summary) <- colnames(summary_all)
  X[[alpha_min_str]]$edge.quality <- vector(length = length(X$K_try), mode = "list")
  for (K_t in X$K_try){
    min_sum <- summary_all[summary_all$K_try == K_t, ]
    min_sum <- min_sum[which.min(min_sum[[lsq_str]]), ]
    # subset(subset(summary_all, K_try == K_t), log_likelihood == min(log_likelihood))
    res_min <- X[[match(paste0("alpha_", min_sum$alpha_name), names(X))]]
    params <- res_min$params_estim[[paste(K_t)]]
    X[[alpha_min_str]]$results_summary[K_t + 1, ] <- as.vector(unname(as.matrix(min_sum)))
    X[[alpha_min_str]]$params_estim[[paste(K_t)]] <- params
    # X[[alpha_min_str]]$params_raw[[paste(K_t)]] <- res_min$params_raw[[paste(K_t)]]
    X[[alpha_min_str]]$params_init_estim[[paste(K_t)]] <- res_min$params_init_estim[[paste(K_t)]]
    X[[alpha_min_str]]$edge.quality[[paste(K_t)]] <- res_min$edge.quality[[paste(K_t)]]
    if (!light_result){
      X[[alpha_min_str]]$Yhat[[paste(K_t)]] <- res_min$Yhat[[paste(K_t)]]
      X[[alpha_min_str]]$Zhat[[paste(K_t)]] <- res_min$Zhat[[paste(K_t)]]
      X[[alpha_min_str]]$Yvar[[paste(K_t)]] <- res_min$Yvar[[paste(K_t)]]
      X[[alpha_min_str]]$Zvar[[paste(K_t)]] <- res_min$Zvar[[paste(K_t)]]
      X[[alpha_min_str]]$m_Y_estim[[paste(K_t)]] <- res_min$m_Y_estim[[paste(K_t)]] 
    }
  }
  X[[alpha_min_str]]$results_summary <- as.data.frame(X[[alpha_min_str]]$results_summary)
  return(X)
}

###############################################################################
## estimateEM_wrapper
###############################################################################

estimateEM_wrapper <- function(use_previous){
  if (use_previous) return(estimateEM_wrapper_previous)
  return(estimateEM_wrapper_scratch)
}

estimateEM_wrapper_previous <- function(phylo, Y_data,
                                        Y_data_imp,
                                        process,
                                        independent = independent,
                                        K_t, prev,
                                        method.variance,
                                        random.root, stationary.root,
                                        alpha_known, alpha,
                                        methods.segmentation,
                                        method.init,
                                        method.init.alpha,
                                        sBM_variance,
                                        method.OUsun,
                                        # impute_init_Rphylopars,
                                        K_lag_init,
                                        allow_negative,
                                        trait_correlation_threshold,
                                        ...){
  tt <- system.time(results_estim_EM <- estimateEM(phylo = phylo, 
                                                   Y_data = Y_data, 
                                                   Y_data_imp = Y_data_imp,
                                                   process = process, 
                                                   independent = independent,
                                                   method.variance = method.variance, 
                                                   method.init = method.init,
                                                   method.init.alpha = method.init.alpha,
                                                   nbr_of_shifts = K_t,
                                                   random.root = random.root,
                                                   stationary.root = stationary.root,
                                                   alpha_known = alpha_known,
                                                   known.selection.strength = alpha,
                                                   init.selection.strength = prev$alpha_estim,
                                                   var.init.root = prev$params_raw$root.state$var.root,
                                                   exp.root.init = prev$params_raw$root.state$exp.root,
                                                   variance.init = prev$params_raw$variance,
                                                   value.root.init = prev$params_raw$root.state$value.root,
                                                   edges.init = prev$params_raw$shifts$edges,
                                                   values.init = prev$params_raw$shifts$values,
                                                   #relativeTimes.init = prev$params_raw$shifts$relativeTimes,
                                                   methods.segmentation = methods.segmentation,
                                                   sBM_variance = sBM_variance,
                                                   method.OUsun = method.OUsun,
                                                   # impute_init_Rphylopars = impute_init_Rphylopars,
                                                   K_lag_init = K_lag_init,
                                                   allow_negative = allow_negative,
                                                   trait_correlation_threshold = trait_correlation_threshold,
                                                   ...))
  return(format_output(results_estim_EM, phylo, tt))
}

estimateEM_wrapper_scratch <- function(phylo, Y_data,
                                       Y_data_imp,
                                       process,
                                       independent,
                                       K_t,
                                       method.variance,
                                       random.root, stationary.root,
                                       alpha_known, alpha,
                                       methods.segmentation,
                                       method.init,
                                       method.init.alpha,
                                       sBM_variance,
                                       method.OUsun,
                                       # impute_init_Rphylopars,
                                       K_lag_init,
                                       allow_negative,
                                       trait_correlation_threshold,
                                       ...){
  tt <- system.time(results_estim_EM <- estimateEM(phylo = phylo, 
                                                   Y_data = Y_data, 
                                                   Y_data_imp = Y_data_imp,
                                                   process = process, 
                                                   independent = independent,
                                                   method.variance = method.variance, 
                                                   method.init = method.init,
                                                   method.init.alpha = method.init.alpha,
                                                   nbr_of_shifts = K_t,
                                                   random.root = random.root,
                                                   stationary.root = stationary.root,
                                                   alpha_known = alpha_known,
                                                   known.selection.strength = alpha,
                                                   methods.segmentation = methods.segmentation,
                                                   sBM_variance = sBM_variance,
                                                   method.OUsun = method.OUsun,
                                                   # impute_init_Rphylopars = impute_init_Rphylopars,
                                                   K_lag_init = K_lag_init,
                                                   allow_negative = allow_negative,
                                                   trait_correlation_threshold = trait_correlation_threshold,
                                                   ...))
  return(format_output(results_estim_EM, phylo, tt))
}

###############################################################################
## Helper Functions
###############################################################################

##
#' @title Test the format of data entry.
#'
#' @description
#' \code{check_data} tests if the data matrix has the right format, and if it is
#' correctly ordered to match the tips names.
#'
#' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
#' @param Y_data matrix of data at the tips (pxntaxa)
#' @param check.tips.names (bool) whether to check the tips names or not
# @param trait_correlation_threshold threshold for trait correlation. Default to 0.9.
#' 
#' @return Y_data a re-ordered matrix of data (if necessary)
#' 
#' @keywords internal
#'
##

check_data <- function(phylo, Y_data, check.tips.names){
  if (is.vector(Y_data)){
    p <- 1
    Y_data <- matrix(Y_data, 1, length(Y_data))
  } else {
    p <- nrow(Y_data)
  }
  if (ncol(Y_data) != length(phylo$tip.label)){
    stop("The data matrix should have as many columns as the number of taxa (p x ntaxa).")
  }
  if (check.tips.names){
    if((is.null(phylo$tip.label) || is.null(colnames(Y_data)))){
      warning("The columns of data matrix and/or the tips of the phylogeny are not named. Could not check for consistency : please make sure that you gave them in the right order.")
    } else {
      if (!all(phylo$tip.label == colnames(Y_data))){
        correspondances <- match(phylo$tip.label, colnames(Y_data))
        if (length(unique(correspondances)) != length(phylo$tip.label)){
          stop("The names of the column data matrix do not match the tip labels.")
        }
        warning("The vector of data was not sorted in the correct order, when compared with the tips label. I am re-ordering the vector of data.")
        Y_data <- Y_data[ , correspondances, drop = FALSE]
      }
    }
  }
  return(as.matrix(Y_data))
}

## Check that correlations between traits are not too high.
check_correlations <- function(Y_data, cor_th = 0.9) {
  cor_data <- stats::cor(t(Y_data))
  cor_data[lower.tri(cor_data, diag = T)] <- 0
  ind_cor <- which(abs(cor_data) >= cor_th, arr.ind = TRUE)
  if (length(ind_cor) > 0) {
    for (x in 1:nrow(ind_cor)) {
      traits_cor <- rownames(Y_data)[ind_cor[x, ]]
      message(paste0("Traits ", traits_cor[1], " and ", traits_cor[2],
                     " have a very high correlation of ",
                     cor_data[ind_cor[x, 1], ind_cor[x, 2]], ".\n",
                     "This correlation measure does not take shifs into account, and could be an artifact of grouped data, in which case this message can be ignored.\n",
                     "If it is not induced by shifts, it might however induce some numerical errors.\n",
                     "If problems arrise, please consider reducing the dimension of the dataset, e.g. using a pre-processing PCA."))
    }
  }
  return(1.0)
}

choose_process_EM <- function(process, p, random.root, stationary.root,
                              alpha_known,
                              known.selection.strength = 1, eps = 10^(-3),
                              sBM_variance = FALSE,
                              method.OUsun = "rescale",
                              independent = FALSE,
                              allow_negative = FALSE){
  ## Reduce Process
  transform_scOU <- FALSE # Should we re-transform back the parameters to get an OU ?
  rescale_tree <- FALSE # Should we re-scale the tree ?
  if (process == "OU"){
    if (p > 1 && !independent) {
      stop("The EM algorithm for shift detection for the general OU is not implemented in the multivariate case. Please consider choosing one of the two following assumptions: (1) the selection strength matrix A is scalar (A = a*I), and the rate matrix R is full (choose process = scOU) or (2) both matices A and R are diagonal, i.e. all the trait are independent (choose independent = TRUE).")
    }
    if (p == 1) {
      process <- "scOU"
    }
  }
  if (process == "scOU"){
    if ((!is.null(known.selection.strength)) && known.selection.strength < 0 && !allow_negative){
      stop("The 'selection strength' you gave is negative. This might not be what you want to do. See manual for the interpretation of such a process. If you really want a negative alpha, please set 'allow_negative=TRUE'.")
    }
    if (random.root){
      if ((!is.null(known.selection.strength)) && known.selection.strength < 0){
        stop("The scalar OU with negative selection strength cannot have a random root. Please try again with 'random_root=FALSE'.")
      }
      if ((method.OUsun == "rescale")){
        if (!alpha_known){
          stop("The re-scaled scalar OU is only implemented for known selection strength. Please consider using a grid.")
        }
        sBM_variance <- TRUE
        if (!stationary.root){
          warning("The scalar OU process with a general random root is not implemented for the multivariate case. The root is taken to be in the stationary state.")
        }
        transform_scOU <- TRUE
        rescale_tree <- TRUE
        process <- "BM"
      } else {
        if (p > 1){
          stop("Only the re-scale method is available for a multivariate scOU.")
        }
        process <- "OU"
      }
    } else {
      if (!alpha_known){
        stop("The re-scaled scalar OU is only implemented for known selection strength. Please consider using a grid.")
      }
      transform_scOU <- TRUE
      rescale_tree <- TRUE
      process <- "BM"
    }
  }
  if ((process == "OU") && alpha_known){
    process <- check.selection.strength(process, known.selection.strength, eps)
  }
  if ((process == "BM") && random.root && !sBM_variance){
    warning("The root parameters cannot be estimated for a BM with a random root. Switching to a fixed root.")
    random.root <- FALSE
  }
  if (process == "rBM"){
    rescale_tree <- TRUE
    process <- "BM"
  }
  return(list(process = process,
              transform_scOU = transform_scOU,
              rescale_tree = rescale_tree,
              sBM_variance = sBM_variance))
}

return_to_original_order <- function(X, phy_o, phy_r){
  reorder_shifts <- function(params){
    if (length(params$shifts$edges) > 0){
      params$shifts$edges <- correspondanceEdges(edges = params$shifts$edges,
                                                 from = phy_r, to = phy_o)
    }
    return(params)
  }
  reorder_edges <- function(edge_qual){
    if (!is.null(names(edge_qual))){
      names(edge_qual) <- correspondanceEdges(edges = as.numeric(names(edge_qual)),
                                              from = phy_r, to = phy_o)
    }
    return(edge_qual)
  }
  reorder_all_shifts <- function(AA){
    AA$params_estim <- lapply(AA$params_estim, reorder_shifts)
    # AA$params_raw <- lapply(AA$params_raw, reorder_shifts)
    AA$params_init_estim <- lapply(AA$params_init_estim, reorder_shifts)
    AA$edge.quality <- lapply(AA$edge.quality, reorder_edges)
    return(AA)
  }
  X[grep("alpha_", names(X))] <- lapply(X[grep("alpha_", names(X))], reorder_all_shifts)
  return(X)
}

##
#' @title Scale the parameters back to the original process
#'
#' @description
#' \code{go_back_to_original_process} takes the inferred parameters with a BM
#' on a rescaled tree, and gives back the equivalent parameters of the OU on 
#' the original process.
#'
#' @param phy_original the original phylogenetic tree
#' @param known.selection.strength the known selection strength of the original
#' OU.
#' @param sBM_variance boolean. Is the root random ?
#' @param params the inferred parameters of the BM on the re-scaled tree.
#' 
#' 
#' @return params_scOU the equivalent parameters of the OU on the original tree.
#' 
#' @keywords internal
#'
##
go_back_to_original_process <- function (phy_original,
                                         known.selection.strength,
                                         sBM_variance, params) {
  params_scOU <- params
  ## lambda parameter
  if (sBM_variance){
    params_scOU$lambda <- params$root.state$exp.root
  } else {
    params_scOU$lambda <- params$root.state$value.root
  }
  ## Default: beta_0 = mu = lambda
  params_scOU$optimal.value <- params_scOU$lambda
  ## shifts values
  params_scOU$shifts <- transform_shifts_values(params$shifts,
                                                from = 0,
                                                to = known.selection.strength,
                                                phylo = phy_original)
  params_scOU$selection.strength <- known.selection.strength
  params_scOU$root.state$stationary.root <- FALSE
  if (sBM_variance){
    params_scOU$root.state$var.root <- params_scOU$variance / (2 * params_scOU$selection.strength)
    params_scOU$root.state$stationary.root <- TRUE
  }
  return(params_scOU)
}

compute_raw_parameters <- function (phy_original,
                                    params) {
  params_rBM <- params
  ## lambda parameter
  params_rBM$lambda <- NULL
  params_rBM$optimal.value <- NULL
  ## shifts values
  params_rBM$shifts <- transform_shifts_values(params$shifts,
                                               from = params$selection.strength,
                                               to = 0,
                                               phylo = phy_original)
  params_rBM$selection.strength <- NULL
  params_rBM$root.state$stationary.root <- FALSE
  return(params_rBM)
}

##
#' @title Run the EM for several values of K
#'
#' @description
#' \code{estimateEM_several_K.OUsr} uses function \code{estimateEM} on the data, 
#' for all values of K between 0 and K_max.
#'
#' @details
#' The EM is first launched for K=0, with alpha and gamma estimated. The
#' estimated values of alpha, gamma and beta_0 found by this first EM are then
#' used as initialization parameters for all the other runs of the EM for other
#' K.
#' The EMs are parallelized thanks to packages \code{foreach} and 
#' \code{doParallel}.
#' WARNING : this code only work of OU with stationary root, on an ultrametric
#' tree.
#' 
#'
#' @param results_estim_EM output of function \code{estimateEM}
#' @param time to run the function
#' 
#' @return summary a data frame with K_max lines, and columns:
#'    - alpha_estim the estimated selection strength
#'    - gamma_estim the estimated root variance
#'    - beta_0_estim the estimated value of root optimum
#'    - EM_steps number of iterations needed before convergence
#'    - DV_estim has the EM diverged ?
#'    - CV_estim has the EM converged ?
#'    - log_likelihood log likelihood of the data using the estimated parameters
#'    - mahalanobis_distance_data_mean the Mahalanobis distance between the data
#' and the estimated means at the tips
#'    - least_squares the Mahalanobis distance, renormalized by gamma^2: 
#' mahalanobis_distance_data_mean * gamma_estim.
#'    - mean_number_new_shifts the mean number of shifts that changed over the 
#' iterations of the EM
#'    - number_equivalent_solutions the number of equivalent solutions to 
#' the solution found.
#'    - K_try the number of shifts allowed.
#'    - complexity the complexity for K_try
#'    - time the CPU time needed.
#' @return params a list of inferred parameters
#' @return params_init a list of initial parameters
#' @return alpha_0 initial values of alpha
#' @return gamma_0 initial values of gamma
#' @return Zhat reconstructed node states
#' @return m_Y_estim reconstructed tip states
#' @return edge.quality for each edge, relative number of iterations in which they
#'  were present.
#'  
#' @keywords internal
#'
##
format_output <- function(results_estim_EM, phylo, time = NA){
  params <- results_estim_EM$params
  params_raw <- results_estim_EM$params_raw
  params_init <- results_estim_EM$params_history[['0']]
  X <- NULL
  X$params <- params
  X$params_raw <- params_raw
  X$params_init <- params_init
  X$alpha_0 <- results_estim_EM$alpha_0
  if (!is.null(X$alpha_0)) names(X$alpha_0) <- paste0("alpha_0_", names(X$alpha_0))
  X$gamma_0 <- results_estim_EM$gamma_0
  # if (!is.null(X$gamma_0)) names(X$gamma_0) <- paste0("gamma_0_", names(X$gamma_0))
  X$Zhat <- results_estim_EM$ReconstructedNodesStates
  X$Zvar <- results_estim_EM$ReconstructedNodesVariances
  X$Yhat <- results_estim_EM$ReconstructedTipsStates
  X$Yvar <- results_estim_EM$ReconstructedTipsVariances
  X$m_Y_estim <- results_estim_EM$m_Y_estim
  #  X$raw_results <- results_estim_EM
  if (is.null(params$selection.strength)){# Handle BM case
    params$selection.strength <- NA
    params_init$selection.strength <- NA
  }
  K_t <- length(params$shifts$edges)
  X$summary <- data.frame(
    ## Estimated Parameters
    # "alpha_estim" = params$selection.strength,
    # "gamma_estim" = params$root.state$var.root,
    # "beta_0_estim" = params$root.state$exp.root,
    "log_likelihood" = attr(params, "log_likelihood")[1],
    "least_squares" = results_estim_EM$least_squares,
    "least_squares_raw" = results_estim_EM$least_squares_raw,
    # "mahalanobis_distance_data_mean" = attr(params, "mahalanobis_distance_data_mean"),
    #"least_squares" = attr(params, "mahalanobis_distance_data_mean") * params$root.state$var.root,
    ## Convergence Monitoring Quantities
    "EM_steps" = attr(results_estim_EM, "Nbr_It"),
    "DV_estim" = attr(results_estim_EM, "Divergence"),
    "CV_estim" = (attr(results_estim_EM, "Nbr_It") != 1000) && !attr(results_estim_EM, "Divergence"),
    "mean_number_new_shifts" = mean(results_estim_EM$number_new_shifts),
    ## Other useful informations
    "number_equivalent_solutions" = results_estim_EM$number_equivalent_solutions,
    "K_try" = K_t,
    "complexity" = extract.partitionsNumber(partitionsNumber(phylo, K_t + 1)),
    "time" = time["elapsed"],
    ## Initial Estimated Parameters
    # "alpha_estim_init" = params_init$selection.strength,
    # "gamma_estim_init" = params_init$root.state$var.root,
    # "beta_0_estim_init" = params_init$root.state$exp.root,
    "log_likelihood_init" = attr(params_init, "log_likelihood")[1]
    # "mahalanobis_distance_data_mean_init" = attr(params_init, "mahalanobis_distance_data_mean")
    #"least_squares_init" = attr(params_init, "mahalanobis_distance_data_mean") * params_init$root.state$var.root
  )
  #   X$summary <- as.data.frame(c(X$summary, X$alpha_0))
  #   X$summary <- as.data.frame(c(X$summary, X$gamma_0))
  ## shifts to be kept for init.
  kpsh <- order(-colSums(params_init$shifts$values))[0:length(params$shifts$edges)]
  results_estim_EM$params_history[['0']]$shifts$edges <- params_init$shifts$edges[kpsh]
  results_estim_EM$params_history[['0']]$shifts$values <- params_init$shifts$values[kpsh]
  results_estim_EM$params_history[['0']]$shifts$relativeTimes <- params_init$shifts$relativeTimes[kpsh]
  ## Compute edge quality
  extract.edges <- function(x) {
    z <- unlist(lapply(x, function(y) y$shifts$edges))
    if (!is.null(z)) z <- matrix(z, nrow = K_t)      
    return(z)
  }
  selected.edges <- extract.edges(results_estim_EM$params_history)
  compute.quality <- function(i) {
    res <- 0 + (selected.edges == i)
    return(mean(colSums(res)))
  }
  if (K_t != 0){
    edge.quality <- unlist(lapply(params$shifts$edges, compute.quality))
    names(edge.quality) <- params$shifts$edges
  } else {
    edge.quality <- NA
  }
  X$edge.quality <- edge.quality
  return(X)
}

# format_output_several_K <- function(res_sev_K, out, alpha = "estimated"){
#   alpha <- paste0("alpha_", alpha)
#   out[[alpha]] <- format_output_several_K_single(res_sev_K)
#   return(out)
# }

format_output_several_K_single <- function(res_sev_K, light_result = TRUE){
  dd <- do.call(rbind, res_sev_K)
  df <- do.call(rbind, dd[ , "summary"])
  df <- as.data.frame(df)
  ## Results
  res <- vector("list")
  res$results_summary <- df
  res$params_estim <- dd[, "params"]
  # res$params_raw <- dd[, "params_raw"]
  res$params_init_estim <- dd[, "params_init"]
  res$alpha_0 <- dd[,"alpha_0" == colnames(dd)]
  if (!light_result){
    res$Zhat <- dd[, "Zhat"]
    res$Yhat <- dd[, "Yhat"]
    res$Zvar <- dd[, "Zvar"]
    res$Yvar <- dd[, "Yvar"]
    res$m_Y_estim <- dd[, "m_Y_estim"]
    if (length(res$params_estim) == 1){
      nn <- paste0(length(res$params_estim[[1]]$shifts$edges))
      names(res$Zhat) <- nn
      names(res$Yhat) <- nn
      names(res$Zvar) <- nn
      names(res$Yvar) <- nn
      names(res$m_Y_estim) <- nn
    }
  }
  res$edge.quality <- dd[, "edge.quality"]
  if (length(res$params_estim) == 1){
    nn <- paste0(length(res$params_estim[[1]]$shifts$edges))
    names(res$params_estim) <- nn
    names(res$params_init_estim) <- nn
    if (length(res$alpha_0) > 0) names(res$alpha_0) <- nn
    names(res$edge.quality) <- nn
  }
  return(res)
}

# ##@title Maximal number of shifts allowed
# 
# @description
# \code{compute_K_max} computes the quantity
# min(floor(kappa * ntaxa / (2 + log(2) + log(ntaxa))), ntaxa - 7))
# that is the maximal dimension allowed to get theoretical garenties during
# the selection model, when using the procedure defined by Baraud et al (2009)
# 
# @details
# See Baraud et al (2009)
# 
# @param ntaxa the number of tips
# @param kappa a real strictly bellow 1.
# 
# @return K_max the maximal number of shifts allowed.
# 
# @keywords internal
# 
# #
# compute_K_max <- function(ntaxa, kappa = 0.9){
#   if (kappa >= 1) stop("For K_max computation, one must have kappa < 1")
#   return(min(floor(kappa * ntaxa / (2 + log(2) + log(ntaxa))), ntaxa - 7))
# }

expand_method_selection <- function(method.selection){
  if ("Djump" %in% method.selection){
    method.selection <- add_method_selection("BirgeMassart1", method.selection)
    method.selection <- method.selection[method.selection != "Djump"]
  }
  if ("DDSE" %in% method.selection){
    method.selection <- add_method_selection("BirgeMassart1", method.selection)
    method.selection <- method.selection[method.selection != "DDSE"]
  }
  if ("LINselect" %in% method.selection){
    for (meth in c("BGHuni", "BGHlsq", "BGHml", "BGHlsqraw", "BGHmlraw")){
      method.selection <- add_method_selection(meth, method.selection)
    }
    method.selection <- method.selection[method.selection != "LINselect"]
  }
  return(method.selection)
}

add_method_selection <- function(meth, method.selection){
  if (!(meth %in% method.selection)){
    method.selection <- c(method.selection, meth)
  }
  return(method.selection)
}