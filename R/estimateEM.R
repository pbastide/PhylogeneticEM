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

##
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
estimateEM <- function(phylo, 
                       Y_data, 
                       process = c("BM","OU"), 
                       tol=list(variance = 10^(-5), 
                                value.root = 10^(-5), 
                                exp.root = 10^(-5), 
                                var.root = 10^(-5),
                                selection.strength = 10^(-5),
                                normalized_half_life = 10^(-5)),  
                       Nbr_It_Max = 500, 
                       method.variance = c("simple"), 
                       method.init = c("default", "lasso"),
                       method.init.alpha = c("default", "estimation"),
                       method.init.alpha.estimation = c("regression", 
                                                        "regression.MM", 
                                                        "median"),
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
                       methods.segmentation = c("max_costs_0", 
                                                "lasso", 
                                                "same_shifts", 
                                                "same_shifts_same_values",
                                                "best_single_move", 
                                                "lasso_one_move"),
                       check.tips.names = FALSE,
                       times_shared = NULL, # These can be specified to save time
                       distances_phylo = NULL, 
                       subtree.list = NULL,
                       T_tree = NULL, 
                       h_tree = NULL,
                       tol_half_life = TRUE, ...){
  
  ## Check consistancy #########################################
  if (alpha_known && missing(known.selection.strength)) stop("The selection strength alpha is supposed to be known, but is not specified. Please add an argument known.selection.strength to the call of the function.")
  
  ## Choose process #########################################
  process <- match.arg(process)
  if ((process == "OU") && alpha_known){
    process <- check.selection.strength(process, known.selection.strength, eps)
    random.root = FALSE
  }
  if ((process == "BM") && random.root){
    warning("The Process is BM with random root : model selection won't work.")
  }
  # specialCase <- stationnary.root && shifts_at_nodes && alpha_known
  compute_M  <- switch(process, 
                       BM = compute_M.BM,
                       OU = compute_M.OU(stationnary.root, shifts_at_nodes, alpha_known))
  shutoff.EM  <- switch(process, 
                        BM = shutoff.EM.BM,
                        OU = shutoff.EM.OU(stationnary.root, shifts_at_nodes, alpha_known, tol_half_life))
  is.finite.params  <- switch(process, 
                              BM = is.finite.params.BM,
                              OU = is.finite.params.OU(stationnary.root, shifts_at_nodes, alpha_known))
  is.in.ranges.params  <- switch(process, 
                                 BM = is.in.ranges.params.BM,
                                 OU = is.in.ranges.params.OU(stationnary.root, shifts_at_nodes, alpha_known))
#   compute_MaxCompleteLogLik <- switch(process, 
#                                       BM = compute_MaxCompleteLogLik.BM,
#                                       OU = compute_MaxCompleteLogLik.OU(stationnary.root, shifts_at_nodes))
#   conditional_expectation_log_likelihood <- switch(process, 
#                                                    BM = conditional_expectation_log_likelihood.BM,
#                                                    OU = conditional_expectation_log_likelihood.OU(stationnary.root, shifts_at_nodes))

  ## init alpha #########################################
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
  compute_mahalanobis_distance  <- switch(method.variance, 
                                          simple = compute_mahalanobis_distance.simple)

  ## Iniialization Method #########################################
  method.init  <- match.arg(method.init)
  # Lasso initialization for OU only works for stationnary root
  if (!stationnary.root && (method.init == "lasso")){
    method.init <- "default"
    warning("The lasso initialization of alpha does only work when the root is stationnary. The initialization is set to the default one.")
  }
  init.EM  <- switch(method.init, 
                     default = init.EM.default(process),
                     lasso = init.EM.lasso)
  method.init.alpha  <- match.arg(method.init.alpha)
  methods.segmentation <- match.arg(methods.segmentation, several.ok = TRUE)
  method.init.alpha.estimation  <- match.arg(method.init.alpha.estimation, several.ok = TRUE)

  ## Fixed Quantities #########################################
  ntaxa <- length(phylo$tip.label)
  if (is.null(times_shared)) times_shared <- compute_times_ca(phylo)
  if (is.null(distances_phylo)) distances_phylo <- compute_dist_phy(phylo)
  if (is.null(subtree.list)) subtree.list <- enumerate_tips_under_edges(phylo)
  if (is.null(T_tree)) T_tree <- incidence.matrix(phylo)
  if (is.null(h_tree)) h_tree <- max(diag(times_shared)[1:ntaxa])

  ## Check that the vector of data is in the correct order and dimensions ################
  Y_data <- check_data(phylo, Y_data, check.tips.names)
  ## Find dimension
  p <- nrow(Y_data)
  ## Initialization
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
                                                  tol = tol,
                                                  h_tree = h_tree)
  init.var.root <- mean(init.a.g$gamma_0[is.finite(init.a.g$gamma_0)])
  if (process == "OU"){
    if(!alpha_known) {
      if ((sum(is.finite(init.a.g$alpha_0)) != 0)){
        ## Only if not all NAs ot infinite
        init.selection.strength <- mean(init.a.g$alpha_0[is.finite(init.a.g$alpha_0)])
      }
    } else {
      init.selection.strength <- known.selection.strength
    }
  }
#   # Always start with some shifts, in case of default initialisation (if number of shifts different from 0)
#   if (!exists("edges.init") || is.null(edges.init)){
#     if (nbr_of_shifts != 0){
#       init_edges <- sample_shifts_edges(phylo, nbr_of_shifts, part.list = subtree.list)
#     }
#     else {
#       init_edges <- NULL
#     }
#   } else {
#     init_edges <- edges.init
#   }
# Initialization per se
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
                         T_tree = T_tree,
                         subtree.list = subtree.list,
                         ...)
  params <- params_init
  params$root.state <- test.root.state(root.state = params$root.state, 
                                       process = process, 
                                       optimal.value = params$optimal.value,
                                       variance = params$variance, 
                                       selection.strength = params$selection.strength)
  attr(params, "ntaxa")  <- ntaxa
  attr(params, "p")  <- p
  params_old <- NULL
  ## Iteration
  Nbr_It <- 0
  params_history <- vector("list")#, Nbr_It_Max)
  #   CLL_history <- NULL
  number_new_shifts <- NULL
  while ( Nbr_It == 0 || # Initialisation
            ( !shutoff.EM(params_old, params, tol, h_tree) && # Shutoff
                is.in.ranges.params(params, min = min_params, max = max_params) && #Divergence?
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
    ## Compute Mahalanobis norm between data and mean at tips
    maha_data_mean <- compute_mahalanobis_distance(phylo = phylo,
                                                   Y_data = Y_data,
                                                   sim = moments$sim,
                                                   Sigma_YY_inv = moments$Sigma_YY_inv)
    attr(params_old, "mahalanobis_distance_data_mean") <- maha_data_mean
    ## E step
    conditional_law_X <- compute_E(phylo = phylo,
                                   Y_data = Y_data,
                                   sim = moments$sim,
                                   Sigma = moments$Sigma,
                                   Sigma_YY_inv = moments$Sigma_YY_inv)
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
                        shifts_old = params_old$shifts,
                        variance_old = params_old$variance,
                        subtree.list = subtree.list)
    attr(params, "ntaxa")  <- ntaxa
    attr(params, "p")  <- p
    ## Number of shifts that changed position ?
    number_new_shifts <- c(number_new_shifts,
                           sum(!(params$shifts$edges %in% params_old$shifts$edges)))
    #   ## Check that the M step rised the conditional expectation of the completed likelihood
    #         CLL_old <- conditional_expectation_log_likelihood(phylo = phylo,
    #                                         conditional_law_X = conditional_law_X, 
    #                                         sigma2 = params_old$variance,
    #                                         mu = params_old$root.state$exp.root,
    #                                         shifts = params_old$shifts,
    #                                         alpha = params_old$selection.strength)
    #         names(CLL_old) <- "CCL_old"
    #         CLL_new <- conditional_expectation_log_likelihood(phylo = phylo,
    #                                         conditional_law_X = conditional_law_X, 
    #                                         sigma2 = params$variance,
    #                                         mu = params$root.state$exp.root,
    #                                         shifts = params$shifts,
    #                                         alpha = params$selection.strength)
    #         names(CLL_new) <- "CCL_new"
    #         if (CLL_old > CLL_new) { warning("The conditional expectation of the completed log likelihood decreased after step M !")}
    #         CLL_history <- cbind(CLL_history, c(CLL_old, CLL_new))
    #         attr(params, "MaxCompleteLogLik") <- CLL_new
  }
  ## Compute log-likelihood for final parameters
  moments <- compute_mean_variance(phylo = phylo,
                                   times_shared = times_shared,
                                   distances_phylo = distances_phylo,
                                   process = process,
                                   params_old = params)
  log_likelihood <- compute_log_likelihood(phylo = phylo,
                                           Y_data = Y_data,
                                           sim = moments$sim,
                                           Sigma = moments$Sigma,
                                           Sigma_YY_inv = moments$Sigma_YY_inv)
  attr(params, "log_likelihood") <- log_likelihood
  params_history[[paste(Nbr_It, sep="")]] <- params
  ## Compute Mahalanobis norm between data and mean at tips
  maha_data_mean <- compute_mahalanobis_distance(phylo = phylo,
                                                 Y_data = Y_data,
                                                 sim = moments$sim,
                                                 Sigma_YY_inv = moments$Sigma_YY_inv)
  attr(params, "mahalanobis_distance_data_mean") <- maha_data_mean
  ## Mean at tips with estimated parameters
  m_Y_estim <- extract.simulate(moments$sim, where="tips", what="expectations")
  ## Number of equivalent solutions
  clusters <- clusters_from_shifts_ism(phylo, params$shifts$edges, part.list = subtree.list)
  Neq <- extract.parsimonyNumber(parsimonyNumber(phylo, clusters))
  if (Neq > 1) message("There are some equivalent solutions to the solution found.")
  ## Result
  result <- list(params = params, 
                 ReconstructedNodesStates = conditional_law_X$expectations[(ntaxa+1):length(conditional_law_X$expectations)],
                 ReconstructedTipsStates = m_Y_estim,
                 params_old = params_old, 
                 params_init = params_init,
                 alpha_0 = init.a.g$alpha_0,
                 gamma_0 = init.a.g$gamma_0,
                 params_history = params_history,
                 number_new_shifts = number_new_shifts,
                 number_equivalent_solutions = Neq)
  #                  CLL_history = CLL_history
  
  ## Handle convergence
  attr(result, "Nbr_It") <- Nbr_It
  attr(result, "Divergence") <- !is.in.ranges.params(result$params, min=min_params, max=max_params) # TRUE if has diverged
  if (Nbr_It == Nbr_It_Max) warning(paste("The maximum number of iterations (Nbr_It_Max = ",Nbr_It_Max,") was reached.",sep=""))
  return(result)
}

##
#' @title Maximal number of shifts allowed
#'
#' @description
#' \code{compute_K_max} computes the quantity 
#' min(floor(kappa * ntaxa / (2 + log(2) + log(ntaxa))), ntaxa - 7))
#' that is the maximal dimention allowed to get theoretical garenties during
#' the selection model, when using the procedure defined by Baraud et al (2009)
#'
#' @details
#' See Baraud et al (2009)
#'
#' @param ntaxa the number of tips
#' @param kappa a real strictly bellow 1.
#' 
#' @return K_max the maximal number of shifts allowed.
#'
##
compute_K_max <- function(ntaxa, kappa = 0.9){
  if (kappa >= 1) stop("For K_max computation, one must have kappa < 1")
  return(min(floor(kappa * ntaxa / (2 + log(2) + log(ntaxa))), ntaxa - 7))
}

##
#' @title Run the EM for several values of K
#'
#' @description
#' \code{estimateEM_several_K.OUsr} uses function \code{estimateEM} on the data, 
#' for all values of K between 0 and K_max.
#'
#' @details
#' The EM is fisrt launched for K=0, with alpha and gamma estimated. The
#' estimated values of alpha, gamma and beta_0 found by this fisrt EM are then
#' used as initialisation parameters for all the other runs of the EM for other
#' K.
#' The EMs are parralelized thanks to packages \code{foreach} and 
#' \code{doParallel}.
#' WARNING : this code only work of OU with stationnary root, on an ultrametric
#' tree.
#' 
#'
#' @param phylo a phylogenetic tree
#' @param Y_data vector of data at the tips
#' @param K_max the maximal number of shifts allowed. By default, computed with 
#' function \code{compute_K_max}.
#' @param ... other arguments to pass to \coed{estimateEM}.
#' 
#' @return summary a data frame with K_max lines, and columns:
#'    - alpha_estim the estimated selection strength
#'    - gamma_estim the estimated root variance
#'    - beta_0_estim the estimated value of root optumum
#'    - EM_steps number of iterations needed before convergence
#'    - DV_estim has the EM diverged ?
#'    - CV_estim has the EM converged ?
#'    - log_likelihood log likelihood of the data using the estimated parameters
#'    - mahalanobis_distance_data_mean the mahalanobis distance between the data
#' and the estimated means at the tips
#'    - least_squares the mahalanobis distance, renormalized by gamma^2: 
#' mahalanobis_distance_data_mean * gamma_estim.
#'    - mean_number_new_shifts the mean number of shifts that changed over the 
#' iterations of the EM
#'    - number_equivalent_solutions the number of equivalent solutions to 
#' the solution found.
#'    - K_try the number of shifts allowed.
#'    - complexity the complexity for K_try
#'    - time the CPU time needed.
#' @return params a list of infered parameters for each EM.
#'
##

estimateEM_several_K.OUsr <- function(phylo, 
                                 Y_data, 
                                 K_max = compute_K_max(length(tree$tip.label), 0.9), 
                                 ...){
  ## Fixed quantities
  times_shared <- compute_times_ca(phylo)
  distances_phylo <- compute_dist_phy(phylo)
  subtree.list <- enumerate_tips_under_edges(phylo)
  T_tree <- incidence.matrix(phylo)
  ## With no shift
  X_0 <- estimation_wrapper.OUsr(0, 
                                 phylo = phylo, 
                                 Y_data = Y_data, 
                                 times_shared = times_shared, 
                                 distances_phylo = distances_phylo,
                                 subtree.list = subtree.list,
                                 T_tree = T_tree, ...)
  beta_0_noshift  <-  X_0$summary$beta_0_estim
  gamma_noshift <- X_0$summary$gamma_estim
  alpha_noshift <- X_0$summary$alpha_estim
  ## Parallel computation using results as initialization
  cl <- makeCluster(Ncores)
  registerDoParallel(cl)
  reqpckg <- c("ape", "quadrupen", "robustbase")
  estimations <- foreach(i = 1:K_max, 
                         .packages = reqpckg, .export=ls(envir=globalenv())) %dopar% {
    estimation_wrapper.OUsr(i, phylo = phylo, Y_data = Y_data,
                            times_shared = times_shared, 
                            distances_phylo = distances_phylo,
                            subtree.list = subtree.list,
                            T_tree = T_tree,
                            method.init.alpha = "default",
                            exp.root.init = beta_0_noshift,
                            var.init.root = gamma_noshift,
                            init.selection.strength = alpha_noshift, ...)
  }
  stopCluster(cl)
  ## Join 0 and rest together
  names(estimations) <- 1:K_max
  X_0_bis <- list('0' = X_0)
  estimations <- c(X_0_bis, estimations)
  dd <- do.call(rbind, estimations)
  df <- do.call(rbind, dd[ , "summary"])
  df <- as.data.frame(df)
  return(list(summary = df, 
              parameters = dd[, "params"]))
}

##
#' @title A wrapper for estimateEM for OU with stationnary root.
#'
#' @description
#' \code{estimation_wrapper.OUsr} call estimateEM with a set of parameters well
#' suited for the OU with stationnary root.
#' It is used in \code{estimateEM_several_K.OUsr}.
#'
#' @param K_t the number of shifts allowed
#' @param phylo a phylogenetic tree
#' @param Y_data vector of data at the tips
#' @param alpha_known a boolean
#' @param alpha the value of alpha if known
#' @param method.init.alpha the initialization method for alpha
#' @param ... other arguments to pass to \coed{estimateEM}.
#' 
#' @return summary a data frame with columns:
#'    - alpha_estim the estimated selection strength
#'    - gamma_estim the estimated root variance
#'    - beta_0_estim the estimated value of root optumum
#'    - EM_steps number of iterations needed before convergence
#'    - DV_estim has the EM diverged ?
#'    - CV_estim has the EM converged ?
#'    - log_likelihood log likelihood of the data using the estimated parameters
#'    - mahalanobis_distance_data_mean the mahalanobis distance between the data
#' and the estimated means at the tips
#'    - least_squares the mahalanobis distance, renormalized by gamma^2: 
#' mahalanobis_distance_data_mean * gamma_estim.
#'    - mean_number_new_shifts the mean number of shifts that changed over the 
#' iterations of the EM
#'    - number_equivalent_solutions the number of equivalent solutions to 
#' the solution found.
#'    - K_try the number of shifts allowed.
#'    - complexity the complexity for K_try
#'    - time the CPU time needed.
#' @return params a list of infered parameters
#' @return params_init a list of initial parameters
#' @return Zhat the reconstructed node states
#' @return edge.quality the quality of each selected edge
#' @return raw_results complete result of \code{estimateEM}
#'
##

estimation_wrapper.OUsr <- function(K_t, phylo, Y_data,
                                    alpha_known = FALSE, alpha = 0,
                                    method.init = "lasso",
                                    method.init.alpha = "estimation",
                                    method.init.alpha.estimation = c("regression",
                                                                     "regression.MM",
                                                                     "median"),
                                    Nbr_It_Max = 1000,
                                    tol_h_l = 10^(-2),
                                    ...) {
  time <- system.time(
    results_estim_EM <- estimateEM(phylo = phylo, 
                                   Y_data = Y_data, 
                                   tol = list(variance=10^(-4), 
                                              value.root=10^(-4), 
                                              exp.root=10^(-4), 
                                              var.root=10^(-4), 
                                              selection.strength=10^(-3),
                                              normalized_half_life = tol_h_l), 
                                   process = "OU", 
                                   method.variance = "simple", 
                                   method.init = method.init,
                                   method.init.alpha = method.init.alpha,
                                   method.init.alpha.estimation = method.init.alpha.estimation,
                                   Nbr_It_Max = Nbr_It_Max, 
                                   nbr_of_shifts = K_t, 
                                   alpha_known = alpha_known, ##
                                   known.selection.strength = alpha,
                                   min_params=list(variance = 10^(-4), 
                                                   value.root = -10^(4), 
                                                   exp.root = -10^(4), 
                                                   var.root = 10^(-4),
                                                   selection.strength = 10^(-4)),
                                   max_params=list(variance = 10^(4), 
                                                   value.root = 10^(4), 
                                                   exp.root = 10^(4), 
                                                   var.root = 10^(4),
                                                   selection.strength = 10^(4)),
                                   methods.segmentation = c("lasso", "best_single_move"), ...)
  )
  params <- results_estim_EM$params
  params_init <- results_estim_EM$params_history['0']$'0'
  X <- NULL
  X$params = params
  X$params_init <- params_init
  X$alpha_0 <- results_estim_EM$alpha_0
  names(X$alpha_0) <- paste0("alpha_0_", names(X$alpha_0))
  X$gamma_0 <- results_estim_EM$gamma_0
  names(X$gamma_0) <- paste0("gamma_0_", names(X$gamma_0))
  X$Zhat <- results_estim_EM$ReconstructedNodesStates
  X$m_Y_estim <- results_estim_EM$ReconstructedTipsStates
#  X$raw_results <- results_estim_EM
  if (is.null(params$selection.strength)){# Handle BM case
    params$selection.strength <- NA
    params_init$selection.strength <- NA
  }
  X$summary <- data.frame(
    ## Estimated Parameters
    "alpha_estim" = params$selection.strength,
    "gamma_estim" = params$root.state$var.root,
    "beta_0_estim" = params$root.state$exp.root,
    "log_likelihood" = attr(params, "log_likelihood")[1],
    "mahalanobis_distance_data_mean" = attr(params, "mahalanobis_distance_data_mean"),
    "least_squares" = attr(params, "mahalanobis_distance_data_mean") * params$root.state$var.root,
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
     "alpha_estim_init" = params_init$selection.strength,
     "gamma_estim_init" = params_init$root.state$var.root,
     "beta_0_estim_init" = params_init$root.state$exp.root,
     "log_likelihood_init" = attr(params_init, "log_likelihood")[1],
     "mahalanobis_distance_data_mean_init" = attr(params_init, "mahalanobis_distance_data_mean"),
     "least_squares_init" = attr(params_init, "mahalanobis_distance_data_mean") * params_init$root.state$var.root
  )
  X$summary <- as.data.frame(c(X$summary, X$alpha_0))
  X$summary <- as.data.frame(c(X$summary, X$gamma_0))
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


##
#' @title Test the format of data entry.
#'
#' @description
#' \code{check_data} tests if the data matrix has the right format, and if it is correctly
#' ordered to match the tips names.
#'
#' @param phylo a phylogenetic tree
#' @param Y_data matrix of data at the tips (pxntaxa)
#' @param check.tips.names (bool) wether to check the tips names or not
#' 
#' @return Y_data a re-ordered matrix of data (if necessary)
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
        Y_data <- Y_data[ , correspondances]
      }
    }
  }
  return(as.matrix(Y_data))
}