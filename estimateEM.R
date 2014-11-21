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
#   compute_MaxCompleteLogLik <- switch(process, 
#                                       BM = compute_MaxCompleteLogLik.BM,
#                                       OU = compute_MaxCompleteLogLik.OU(stationnary.root, shifts_at_nodes))
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
  ## Iniialization Method
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
  ## Fixed Quantities
  ntaxa <- length(phylo$tip.label)
  times_shared <- compute_times_ca(phylo)
  distances_phylo <- compute_dist_phy(phylo)
#  t_tree <-  min(node.depth.edgelength(phylo)[1:ntaxa])
  subtree.list <- enumerate_tips_under_edges(phylo)
  T_tree <- incidence.matrix(phylo)
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
                         T_tree = T_tree,
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
                        subtree.list = subtree.list)
    attr(params, "ntaxa")  <- ntaxa
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
  ## Number of equivalent solutions
  clusters <- clusters_from_shifts(phylo, params$shifts$edges)
  Neq <- extract.parsimonyNumber(parsimonyNumber(phylo, clusters))
  if (Neq > 1) message("There are some equivalent solutions to the solution found.")
  ## Result
  result <- list(params = params, 
                 ReconstructedNodesStates = conditional_law_X$expectations[(ntaxa+1):length(conditional_law_X$expectations)], 
                 params_old = params_old, 
                 params_init = params_init,
                 params_history = params_history,
                 number_new_shifts = number_new_shifts,
                 number_equivalent_solutions = Neq)
  #                  CLL_history = CLL_history
  
  attr(result, "Nbr_It") <- Nbr_It
  attr(result, "Divergence") <- !is.in.ranges.params(result$params, min=min_params, max=max_params) # TRUE if has diverged
  if (Nbr_It == Nbr_It_Max) warning(paste("The maximum number of iterations (Nbr_It_Max = ",Nbr_It_Max,") was reached.",sep=""))
  return(result)
}