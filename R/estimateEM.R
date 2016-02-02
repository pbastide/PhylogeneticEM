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
                       Y_data_imp = Y_data,
                       process = c("BM", "OU", "scOU", "rBM"), 
                       tol = list(variance = 10^(-3), 
                                  value.root = 10^(-3), 
                                  exp.root = 10^(-3), 
                                  var.root = 10^(-3),
                                  selection.strength = 10^(-3),
                                  normalized_half_life = 10^(-3),
                                  log_likelihood = 10^(-3)),  
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
                       min_params=list(variance = 10^(-5), 
                                       value.root = -10^(5), 
                                       exp.root = -10^(5), 
                                       var.root = 10^(-5),
                                       selection.strength = 10^(-5)),
                       max_params=list(variance = 10^(5), 
                                       value.root = 10^(5), 
                                       exp.root = 10^(5), 
                                       var.root = 10^(5),
                                       selection.strength = 10^(5)),
                       var.init.root = diag(1, nrow(Y_data)),
                       variance.init = diag(1, nrow(Y_data), nrow(Y_data)),
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
                       tol_half_life = TRUE,
                       warning_several_solutions = TRUE,
                       convergence_mode = c("relative", "absolute"),
                       check_convergence_likelihood = TRUE,
                       sBM_variance = FALSE,
                       method.OUsun = c("rescale", "raw"),
                       impute_init_Rphylopars = TRUE,
                       ...){
  
  ntaxa <- length(phylo$tip.label)
  ## Check that the vector of data is in the correct order and dimensions ################
  Y_data <- check_data(phylo, Y_data, check.tips.names)
  ## Find dimension
  p <- nrow(Y_data)
  
  ## Root edge of the tree ###############################################################
  if (is.null(phylo$root.edge)) phylo$root.edge <- 0
  
  ########## Check consistancy ################################################
  if (alpha_known && missing(known.selection.strength)) stop("The selection strength alpha is supposed to be known, but is not specified. Please add an argument known.selection.strength to the call of the function.")
#  known.selection.strength <- check_dimensions.matrix(p, p, known.selection.strength, "known.selection.strength")
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
                            stationnary.root = stationnary.root,
                            alpha_known = alpha_known,
                            known.selection.strength = known.selection.strength,
                            eps = eps,
                            sBM_variance = sBM_variance,
                            method.OUsun = method.OUsun)
  process <- temp$process
  transform_scOU <- temp$transform_scOU # Transform back to get an OU ?
  rescale_tree <- temp$rescale_tree # Rescale the tree ?
  sBM_variance <- temp$sBM_variance
  if (sBM_variance) phylo$root.edge <- 1
  
  ########## Choose functions #################################################
  # specialCase <- stationnary.root && shifts_at_nodes && alpha_known
  compute_M  <- switch(process, 
                       BM = compute_M.BM,
                       OU = compute_M.OU(stationnary.root, shifts_at_nodes, alpha_known))
  shutoff.EM  <- switch(process, 
                        BM = shutoff.EM.BM,
                        OU = shutoff.EM.OU(stationnary.root, shifts_at_nodes, alpha_known, tol_half_life))
  has_converged  <- switch(convergence_mode[1], 
                           relative = has_converged_relative,
                           absolute = has_converged_absolute)
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

  ########## init alpha #######################################################
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
  
  ########## Moments computation method #######################################
  method.variance  <- match.arg(method.variance)
  compute_E  <- switch(method.variance, 
                       simple = compute_E.simple)
  compute_mean_variance  <- switch(method.variance, 
                                   simple = compute_mean_variance.simple)
  compute_log_likelihood  <- switch(method.variance, 
                                    simple = compute_log_likelihood.simple)
  compute_mahalanobis_distance  <- switch(method.variance, 
                                          simple = compute_mahalanobis_distance.simple)

  ########## Initialization Method ############################################
  method.init  <- match.arg(method.init)
  # Lasso initialization for OU only works for stationnary root
#   if (!stationnary.root && (method.init == "lasso")){
#     method.init <- "default"
#     warning("The lasso initialization of alpha does only work when the root is stationnary. The initialization is set to the default one.")
#   }
  init.EM  <- switch(method.init, 
                     default = init.EM.default(process),
                     lasso = init.EM.lasso)
  method.init.alpha  <- match.arg(method.init.alpha)
  methods.segmentation <- match.arg(methods.segmentation, several.ok = TRUE)
  method.init.alpha.estimation  <- match.arg(method.init.alpha.estimation, several.ok = TRUE)
  
  ########## Fixed Quantities #################################################
  ntaxa <- length(phylo$tip.label)
  ## Transform the branch lengths if needed
  phy_original <- phylo
  if (rescale_tree){
    phylo <- transform_branch_length(phylo, known.selection.strength)
  }
  if (is.null(times_shared)) times_shared <- compute_times_ca(phylo)
  if (is.null(distances_phylo)) distances_phylo <- compute_dist_phy(phylo)
  if (is.null(subtree.list)) subtree.list <- enumerate_tips_under_edges(phylo)
  if (is.null(T_tree)) T_tree <- incidence.matrix(phylo)
  if (is.null(h_tree)) h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
  
  ########## Re-scale tree to 100 #############################################
  
  # factor_rescale <- min(phy_original$edge.length) / min(phylo$edge.length) # min new = min old
  factor_rescale <- 1 / h_tree # total height to 1
  # factor_rescale <- 1
  
  h_tree <- factor_rescale * h_tree
  times_shared <- factor_rescale * times_shared
  distances_phylo <- factor_rescale * distances_phylo
  phylo$edge.length <- factor_rescale * phylo$edge.length
  phylo$root.edge <- factor_rescale * phylo$root.edge
  
  known.selection.strength <-  known.selection.strength / factor_rescale
  init.selection.strength <- init.selection.strength / factor_rescale
  variance.init <- variance.init / factor_rescale
  
  ########## Missing Data #####################################################
  ntaxa <- length(phylo$tip.label)
  nNodes <- phylo$Nnode
  p <- nrow(Y_data)
  missing <- as.vector(is.na(Y_data))
  Y_data_vec <- as.vector(Y_data)
  Y_data_vec_known <- as.vector(Y_data[!missing])
  # Vectorized Data Mask
  masque_data <- rep(FALSE, (ntaxa + nNodes) * p)
  masque_data[1:(p*ntaxa)] <- !missing
  
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
                                                  tol = tol,
                                                  h_tree = h_tree,
                                                  random.root = random.root,
                                                  T_tree = T_tree,
                                                  h_tree = h_tree,
                                                  subtree.list = subtree.list,
                                                  missing = missing,
                                                  ...)
  init.var.root <- init.a.g$gamma_0 # mean(init.a.g$gamma_0[is.finite(init.a.g$gamma_0)])
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
  if (process == "BM"){
    init.selection.strength <- 0
  }
  
  ## Init of Rate matrix for BM
  if (process == "BM" && (!random.root || (random.root && sBM_variance))){
    variance.init <- init.variance.BM.estimation(phylo = phylo, 
                                                 Y_data = Y_data, 
                                                 Y_data_imp = Y_data_imp,
                                                 nbr_of_shifts = nbr_of_shifts, 
                                                 times_shared = times_shared,
                                                 distances_phylo = distances_phylo, 
                                                 h_tree = h_tree,
                                                 random.root = random.root,
                                                 T_tree = T_tree,
                                                 subtree.list = subtree.list,
                                                 missing = missing,
                                                 selection.strength.init = init.selection.strength,
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
                         selection.strength.init = init.selection.strength, 
                         random.init = random.root,
                         stationnary.root.init = stationnary.root,
                         use_sigma = use_sigma_for_lasso,
                         method.init.alpha = method.init.alpha,
                         var.root.init = init.var.root,
                         T_tree = T_tree,
                         subtree.list = subtree.list,
                         missing = missing,
                         variance.init = variance.init,
                         sBM_variance = sBM_variance,
                         impute_init_Rphylopars = impute_init_Rphylopars,
                         masque_data = masque_data,
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
  while ( Nbr_It == 0 || # Initialisation
            ( !(CV_log_lik && # CV of log-Likelihood ?
              shutoff.EM(params_old, params, tol, has_converged, h_tree)) && # Shutoff
                is.in.ranges.params(params, min = min_params, max = max_params) && #Divergence?
                Nbr_It < Nbr_It_Max ) ) { # Nbr of iteration
    ## Actualization
    Nbr_It <- Nbr_It + 1
    params_old <- params
    ########## Log Likelihood #################################################
    # Check convergence of loglik ?
    if (check_convergence_likelihood && Nbr_It > 1){
      log_likelihood_old_old <- log_likelihood_old
    }
    moments <- compute_mean_variance(phylo = phylo,
                                     times_shared = times_shared,
                                     distances_phylo = distances_phylo,
                                     process = process,
                                     params_old = params_old,
                                     masque_data = masque_data)
    log_likelihood_old <- compute_log_likelihood(phylo = phylo,
                                                 Y_data_vec = Y_data_vec_known,
                                                 sim = moments$sim,
                                                 Sigma = moments$Sigma,
                                                 Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                                 missing = missing, 
                                                 masque_data = masque_data)
    attr(params_old, "log_likelihood") <- log_likelihood_old
    if (check_convergence_likelihood && Nbr_It > 1){
      CV_log_lik <- has_converged_absolute(log_likelihood_old_old, log_likelihood_old, tol$log_likelihood)
    }
    ## Compute Mahalanobis norm between data and mean at tips
    maha_data_mean <- compute_mahalanobis_distance(phylo = phylo,
                                                   Y_data_vec = Y_data_vec_known,
                                                   sim = moments$sim,
                                                   Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                                   missing = missing)
    attr(params_old, "mahalanobis_distance_data_mean") <- maha_data_mean
    
    ########## E step #########################################################
    conditional_law_X <- compute_E(phylo = phylo,
                                   Y_data_vec = Y_data_vec_known,
                                   sim = moments$sim,
                                   Sigma = moments$Sigma,
                                   Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                   missing = missing,
                                   masque_data = masque_data)
    rm(moments)
    if (process == "OU"){
      if (p > 1) stop("Multivariate OU not yet implemented.")
      conditional_law_X$expectations <- as.vector(conditional_law_X$expectations)
      conditional_law_X$variances <- as.vector(conditional_law_X$variances)
      conditional_law_X$covariances <- as.vector(conditional_law_X$covariances)
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
    ## Store params for history
    params_history[[paste(Nbr_It - 1, sep="")]] <- params_old
    
    ########## M step #########################################################
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
                        mu_old = params_old$root.state$value.root,
                        subtree.list = subtree.list,
                        sBM_variance = sBM_variance)
    attr(params, "ntaxa")  <- ntaxa
    attr(params, "p_dim")  <- p
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
  
  ########## Scale back parameters to original tree ###########################
  
  params <- scale_params(params, factor_rescale)
  
  h_tree <- h_tree / factor_rescale
  times_shared <- times_shared / factor_rescale
  distances_phylo <- distances_phylo / factor_rescale
  phylo$edge.length <- phylo$edge.length / factor_rescale
  phylo$root.edge <- phylo$root.edge / factor_rescale
  
  known.selection.strength <-  known.selection.strength * factor_rescale
  init.selection.strength <- init.selection.strength * factor_rescale
  variance.init <- variance.init * factor_rescale
  
  ########## Go back to OU parameters if needed ###############################
  params_scOU <- params # If a BM, params_scOU = params
  if (transform_scOU){
    ## Go back to original tree and process
    phylo <- phy_original
    times_shared <- compute_times_ca(phy_original)
    distances_phylo <- compute_dist_phy(phy_original)
    process <- suppressWarnings(check.selection.strength(original_process,
                                                         known.selection.strength,
                                                         eps))
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
                                                  phylo = phylo)
    params_scOU$selection.strength <- known.selection.strength
    params_scOU$root.state$stationary.root <- FALSE
    if (sBM_variance){
      params_scOU$root.state$var.root <- params_scOU$variance / (2 * params_scOU$selection.strength)
      params_scOU$root.state$stationary.root <- TRUE
    }
  }
  
  ########## Compute scores and ancestral states for final parameters ##########
  ## Compute log-likelihood for final parameters
  moments <- compute_mean_variance(phylo = phylo,
                                   times_shared = times_shared,
                                   distances_phylo = distances_phylo,
                                   process = process,
                                   params_old = params_scOU,
                                   masque_data = masque_data)
  log_likelihood <- compute_log_likelihood(phylo = phylo,
                                           Y_data_vec = Y_data_vec_known,
                                           sim = moments$sim,
                                           Sigma = moments$Sigma,
                                           Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                           missing = missing,
                                           masque_data = masque_data)
  attr(params_scOU, "log_likelihood") <- log_likelihood
  params_history[[paste(Nbr_It, sep="")]] <- params_scOU
  ## Compute Mahalanobis norm between data and mean at tips
  maha_data_mean <- compute_mahalanobis_distance(phylo = phylo,
                                                 Y_data_vec = Y_data_vec_known,
                                                 sim = moments$sim,
                                                 Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                                 missing = missing)
  attr(params_scOU, "mahalanobis_distance_data_mean") <- maha_data_mean
  ## "Ancestral state reconstruction"
  conditional_law_X <- compute_E(phylo = phylo,
                                 Y_data_vec = Y_data_vec_known,
                                 sim = moments$sim,
                                 Sigma = moments$Sigma,
                                 Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                 missing = missing,
                                 masque_data = masque_data)
  ## Mean at tips with estimated parameters
  m_Y_estim <- extract.simulate(moments$sim, where="tips", what="expectations")
  ## Number of equivalent solutions
  clusters <- clusters_from_shifts_ism(phylo, params$shifts$edges, part.list = subtree.list)
  Neq <- extract.parsimonyNumber(parsimonyNumber(phylo, clusters))
  if (Neq > 1 && warning_several_solutions) message("There are some equivalent solutions to the solution found.")
  attr(params_scOU, "Neq") <- Neq
  
  ########## Result  ##########################################################
  conditional_law_X$expectations <- matrix(conditional_law_X$expectations, nrow = p)
  result <- list(params = params_scOU, # Return untransformed parameters as default
                 params_raw = params,
                 ReconstructedNodesStates = conditional_law_X$expectations[ , (ntaxa+1):ncol(conditional_law_X$expectations)],
                 ReconstructedTipsStates = conditional_law_X$expectations[ , 1:ntaxa],
                 ReconstructedNodesVariances = conditional_law_X$variances[ , , (ntaxa+1):ncol(conditional_law_X$expectations)],
                 ReconstructedTipsVariances = conditional_law_X$variances[ , , 1:ntaxa],
                 m_Y_estim = m_Y_estim,
                 params_old = params_old, 
                 params_init = params_init,
                 alpha_0 = init.a.g$alpha_0,
                 gamma_0 = init.a.g$gamma_0,
                 params_history = params_history,
                 number_new_shifts = number_new_shifts,
                 number_equivalent_solutions = Neq)
  #                  CLL_history = CLL_history
  if (transform_scOU) result$params_scOU <-  params_scOU
  ## Handle convergence
  attr(result, "Nbr_It") <- Nbr_It
  attr(result, "Divergence") <- !is.in.ranges.params(result$params, min=min_params, max=max_params) # TRUE if has diverged
  if (Nbr_It == Nbr_It_Max) warning(paste("The maximum number of iterations (Nbr_It_Max = ",Nbr_It_Max,") was reached.",sep=""))
  return(result)
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
#' @param results_estim_EM output of function \code{estimateEM}
#' @param time to run the function
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
#' @return params a list of infered parameters
#' @return params_init a list of initial paramters
#' @return alpha_0 initial values of alpha
#' @return gmma_0 initial values of gamma
#' @return Zhat reconstructed node states
#' @return m_Y_estim reconstructed tip states
#' @return edge.quality for each edge, relative number of iterations in which they
#'  were present.
#'
##
format_output <- function(results_estim_EM, phylo, time = NA){
  params <- results_estim_EM$params
  params_raw <- results_estim_EM$params_raw
  params_init <- results_estim_EM$params_history['0']$'0'
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
    "alpha_estim" = params$selection.strength,
    # "gamma_estim" = params$root.state$var.root,
    "beta_0_estim" = params$root.state$exp.root,
    "log_likelihood" = attr(params, "log_likelihood")[1],
    "mahalanobis_distance_data_mean" = attr(params, "mahalanobis_distance_data_mean"),
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
    "beta_0_estim_init" = params_init$root.state$exp.root,
    "log_likelihood_init" = attr(params_init, "log_likelihood")[1],
    "mahalanobis_distance_data_mean_init" = attr(params_init, "mahalanobis_distance_data_mean")
    #"least_squares_init" = attr(params_init, "mahalanobis_distance_data_mean") * params_init$root.state$var.root
  )
#   X$summary <- as.data.frame(c(X$summary, X$alpha_0))
#   X$summary <- as.data.frame(c(X$summary, X$gamma_0))
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

format_output_several_K <- function(res_sev_K, out, alpha = "estimated"){
  alpha <- paste0("alpha_", alpha)
  out[[alpha]] <- format_output_several_K_single(res_sev_K)
  return(out)
}

format_output_several_K_single <- function(res_sev_K){
  dd <- do.call(rbind, res_sev_K)
  df <- do.call(rbind, dd[ , "summary"])
  df <- as.data.frame(df)
  ## Results
  res <- vector("list")
  res$results_summary <- df
  res$params_estim <- dd[, "params"]
  res$params_raw <- dd[, "params_raw"]
  res$params_init_estim <- dd[, "params_init"]
  res$alpha_0 <- dd[,"alpha_0" == colnames(dd)]
  res$Zhat <- dd[, "Zhat"]
  res$Yhat <- dd[, "Yhat"]
  res$Zvar <- dd[, "Zvar"]
  res$Yvar <- dd[, "Yvar"]
  res$m_Y_estim <- dd[, "m_Y_estim"]
  res$edge.quality <- dd[, "edge.quality"]
  return(res)
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

PhyloEM <- function(phylo, Y_data, process = c("BM", "OU", "scOU", "rBM"),
                    K_max, use_previous = TRUE,
                    order = TRUE,
                    method.selection = c("BirgeMassart1", "BirgeMassart2", "BGH"),
                    C.BM1 = 0.1, C.BM2 = 2.5, C.BGH = 1.1,
                    method.variance = "simple",
                    method.init = "default",
                    method.init.alpha = "default",
                    method.init.alpha.estimation = c("regression", "regression.MM", "median"), 
                    methods.segmentation = c("lasso", "best_single_move"),
                    alpha_known = TRUE,
                    random.root = FALSE,
                    stationnary.root = FALSE,
                    alpha = NULL,
                    check.tips.names = FALSE,
                    progress.bar = TRUE,
                    estimates = NULL,
                    save_step = FALSE,
                    sBM_variance = FALSE,
                    method.OUsun = "rescale",
                    parallel_alpha = FALSE,
                    Ncores = 3,
                    exportFunctions = ls(),
                    impute_init_Rphylopars = TRUE,
                    ...){
  ## Required packages
  library(doParallel)
  library(foreach)
  library(ape)
  library(glmnet) # For Lasso initialization
  library(robustbase) # For robust fitting of alpha
  reqpckg <- c("ape", "glmnet", "robustbase")
  ## Check the tree
  if (!is.ultrametric(phylo)) stop("The tree must be ultrametric.")
  ## Check that the vector of data is in the correct order and dimensions ################
  Y_data <- check_data(phylo, Y_data, check.tips.names)
  p <- nrow(Y_data)
  ntaxa <- length(phylo$tip.label)
  ## Model Selection
  method.selection  <- match.arg(method.selection, several.ok = TRUE)
  if (p > 1 && "BGH" %in% method.selection){
    method.selection <- method.selection[-which(method.selection == "BGH")]
    if (length(method.selection) == 0) stop("BGH is not implemented for multivariate data.")
    warning("BGH is not implemented for multivariate data.")
  }
  if (method.selection == "BirgeMassart1" || method.selection == "BirgeMassart2") library(capushe)
  
  ## Save Original process and tree
  process <- match.arg(process)
  process_original <- process
  original_phy <- phylo
  
  ## Compute alpha
  if (process == "BM") {
    alpha <- 0
  } else {
    alpha <- find_grid_alpha(phylo, alpha, ...)
  }
  
  if (!is.null(estimates)){
    ## Set Up
    X <- estimates ## Get estimations from presiously computed results
  } else {
    ## Loop on alpha
    estimate_alpha_several_K <- function(alp, ...){
      ## Process
#       process <- match.arg(process)
#       process_original <- process
#       original_phy <- phylo
      temp <- choose_process_EM(process = process_original,
                                p = p,
                                random.root = random.root,
                                stationnary.root = stationnary.root,
                                alpha_known = alpha_known,
                                known.selection.strength = alp,
                                sBM_variance = sBM_variance,
                                method.OUsun = method.OUsun)
      
      rescale_tree <- temp$rescale_tree # Rescale the tree ?
      transform_scOU <- temp$transform_scOU # Re-transform parameters back ?
      sBM_variance <- temp$sBM_variance
      if (sBM_variance){ # process sBM : need a root branch.
        original_phy$root.edge <- 1
      }
      ## Transform branch lengths if needed
      phylo <- original_phy
      if (rescale_tree) {
        phylo <- transform_branch_length(phylo, alp)
      }
      ## Fixed quantities
      times_shared <- compute_times_ca(phylo)
      distances_phylo <- compute_dist_phy(phylo)
      subtree.list <- enumerate_tips_under_edges(phylo)
      T_tree <- incidence.matrix(phylo)
      h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
      ## Impute data if needed
      Y_data_imp <- Y_data
      if (impute_init_Rphylopars && temp$process == "BM"){
        ## Re-scale tree to unit height
        factor_rescale <- 1 / h_tree # total height to 1
        phylo_temp <- phylo
        phylo_temp$edge.length <- factor_rescale * phylo$edge.length
        phylo_temp$root.edge <- factor_rescale * phylo$root.edge
        Y_data_imp <- try(impute.data.Rphylopars(phylo_temp,
                                                 Y_data,
                                                 temp$process,
                                                 random.init))
        if (inherits(Y_data_imp, "try-error")) { # If fails, replace with mean of the trait
          Y_data_imp <- Y_data
          for (j in 1:(dim(Y_data_imp)[1])){
            Y_data_imp[j, is.na(Y_data_imp[j, ])] <- mean(Y_data_imp[j, ], na.rm = TRUE)
          }
        }
        rm(phylo_temp)
      }
      ## Estimations
      X <- Phylo_EM_sequencial(phylo = original_phy,
                               Y_data = Y_data,
                               Y_data_imp = Y_data_imp,
                               process = process,
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
                               stationnary.root = stationnary.root,
                               alp = alp,
                               check.tips.names = check.tips.names,
                               progress.bar = progress.bar,
                               times_shared = times_shared,
                               distances_phylo = distances_phylo,
                               subtree.list = subtree.list,
                               T_tree = T_tree,
                               h_tree = h_tree,
                               save_step = save_step,
                               sBM_variance = sBM_variance,
                               method.OUsun = method.OUsun,
                               impute_init_Rphylopars = impute_init_Rphylopars,
                               ...)
    }
    if (parallel_alpha){
      cl <- makeCluster(Ncores)
      registerDoParallel(cl)
      X <- foreach(alp = alpha, .packages = reqpckg, .export = exportFunctions) %dopar%
      {
        estimate_alpha_several_K(alp)
      }
      stopCluster(cl)
    } else {
      X <- foreach(alp = alpha, .packages = reqpckg) %do%
      {
        estimate_alpha_several_K(alp)
      }
    }
    names(X) <- paste0("alpha_", alpha)
    X$Y_data <- Y_data
    X$K_try <- 0:K_max
    X$ntaxa <- ntaxa
  }
  ## Select max solution for each K
  X <- merge_max_grid_alpha(X, alpha)
  ## Model Selection
  model_selection <- function(one.method.selection){
    mod_sel  <- switch(one.method.selection, 
                       BirgeMassart1 = model_selection_BM1,
                       BirgeMassart2 = model_selection_BM2,
                       BGH = model_selection_BGH)
    selection <- try(mod_sel(X$alpha_max, ntaxa = ncol(Y_data),
                             C.BM1 = C.BM1, C.BM2 = C.BM2, C.BGH = C.BGH))
    if (inherits(selection, "try-error")){
      warning(paste0("Model Selection ",  one.method.selection, " failled"))
    } else {
      X$alpha_max <- selection
    }
    return(X)
  }
  ## Selection(s)
  for (meth.sel in method.selection){
    X <- model_selection(meth.sel)
  }
  return(X)
}

Phylo_EM_sequencial <- function(phylo, Y_data,
                                Y_data_imp,
                                process, K_max,
                                #                                 curent = list(Y_data = Y_data,
                                #                                               K_try = 0:K_max,
                                #                                               ntaxa = length(phylo$tip.label)),
                                use_previous = TRUE,
                                order = TRUE,
                                method.variance = "simple",
                                method.init = "default",
                                method.init.alpha = "default",
                                method.init.alpha.estimation = c("regression",
                                                                 "regression.MM", "median"), 
                                methods.segmentation = c("lasso", "best_single_move"),
                                alpha_known = FALSE,
                                random.root = TRUE,
                                stationnary.root = TRUE,
                                alp = alp,
                                check.tips.names = FALSE,
                                progress.bar = TRUE,
                                times_shared = NULL,
                                distances_phylo = NULL,
                                subtree.list = NULL,
                                T_tree = NULL,
                                h_tree = NULL,
                                save_step = TRUE,
                                sBM_variance = FALSE,
                                method.OUsun = "rescale", 
                                impute_init_Rphylopars = TRUE,
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
    message(paste0("Alpha ", alp))
    pb <- txtProgressBar(min = 0, max = K_max + 1, style = 3)
  }
  ## First
  XX[[paste0(K_first)]] <- estimateEM_wrapper_scratch(phylo = phylo,
                                                      Y_data = Y_data,
                                                      Y_data_imp = Y_data_imp,
                                                      process = process,
                                                      K_t = K_first,
                                                      method.variance = method.variance,
                                                      random.root = random.root,
                                                      stationnary.root = stationnary.root,
                                                      alpha_known = alpha_known,
                                                      alpha = alp,
                                                      method.init = method.init,
                                                      method.init.alpha = method.init.alpha,
                                                      methods.segmentation = methods.segmentation,
                                                      times_shared = times_shared, 
                                                      distances_phylo = distances_phylo,
                                                      subtree.list = subtree.list,
                                                      T_tree = T_tree,
                                                      h_tree = h_tree,
                                                      warning_several_solutions = FALSE,
                                                      sBM_variance = sBM_variance,
                                                      method.OUsun = method.OUsun,
                                                      impute_init_Rphylopars = impute_init_Rphylopars,
                                                      ...)
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
  for (K_t in (next_it(K_first)):(K_last)){
    XX[[paste0(K_t)]] <- estimateEM_wrapper(use_previous)(phylo = phylo,
                                                          Y_data = Y_data,
                                                          Y_data_imp = Y_data_imp,
                                                          process = process,
                                                          K_t = K_t,
                                                          prev = XX[[paste0(prev_it(K_t))]],
                                                          method.variance = method.variance,
                                                          random.root = random.root,
                                                          stationnary.root = stationnary.root,
                                                          alpha_known = alpha_known,
                                                          alpha = alp,
                                                          method.init = method.init,
                                                          method.init.alpha = method.init.alpha,
                                                          methods.segmentation = methods.segmentation,
                                                          times_shared = times_shared, 
                                                          distances_phylo = distances_phylo,
                                                          subtree.list = subtree.list,
                                                          T_tree = T_tree, 
                                                          h_tree = h_tree,
                                                          warning_several_solutions = FALSE,
                                                          sBM_variance = sBM_variance,
                                                          method.OUsun = method.OUsun,
                                                          impute_init_Rphylopars = impute_init_Rphylopars,
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
  ## Format results and return
  res <- format_output_several_K_single(XX)
  if (save_step) save(res, file = paste0("Tmp_", alp, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
  return(res)
}

merge_max_grid_alpha <- function(X, alpha){
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
    max_sum <- subset(subset(summary_all, K_try == K_t), log_likelihood == max(log_likelihood))
    res_max <- X[[paste0("alpha_", max_sum$alpha_name)]]
    params <- res_max$params_estim[[paste(K_t)]]
    X$alpha_max$results_summary[K_t + 1, ] <- as.vector(unname(as.matrix(max_sum)))
    X$alpha_max$params_estim[[paste(K_t)]] <- params
    X$alpha_max$params_raw[[paste(K_t)]] <- res_max$params_raw[[paste(K_t)]]
    X$alpha_max$params_init_estim[[paste(K_t)]] <- res_max$params_init_estim[[paste(K_t)]]
    X$alpha_max$Yhat[[paste(K_t)]] <- res_max$Yhat[[paste(K_t)]]
    X$alpha_max$Zhat[[paste(K_t)]] <- res_max$Zhat[[paste(K_t)]]
    X$alpha_max$Yvar[[paste(K_t)]] <- res_max$Yvar[[paste(K_t)]]
    X$alpha_max$Zvar[[paste(K_t)]] <- res_max$Zvar[[paste(K_t)]]
    X$alpha_max$edge.quality[[paste(K_t)]] <- res_max$edge.quality[[paste(K_t)]]
    X$alpha_max$m_Y_estim[[paste(K_t)]] <- res_max$m_Y_estim[[paste(K_t)]]
  }
  X$alpha_max$results_summary <- as.data.frame(X$alpha_max$results_summary)
  return(X)
}

estimateEM_wrapper <- function(use_previous){
  if (use_previous) return(estimateEM_wrapper_previous)
  return(estimateEM_wrapper_scratch)
}

estimateEM_wrapper_previous <- function(phylo, Y_data,
                                        Y_data_imp,
                                        process, K_t, prev,
                                        method.variance,
                                        random.root, stationnary.root,
                                        alpha_known, alpha,
                                        methods.segmentation,
                                        method.init,
                                        method.init.alpha,
                                        sBM_variance,
                                        method.OUsun,
                                        impute_init_Rphylopars,
                                        ...){
  tt <- system.time(results_estim_EM <- estimateEM(phylo = phylo, 
                                                   Y_data = Y_data, 
                                                   Y_data_imp = Y_data_imp,
                                                   process = process, 
                                                   method.variance = method.variance, 
                                                   method.init = method.init,
                                                   method.init.alpha = method.init.alpha,
                                                   nbr_of_shifts = K_t,
                                                   random.root = random.root,
                                                   stationnary.root = stationnary.root,
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
                                                   impute_init_Rphylopars = impute_init_Rphylopars,
                                                   ...))
  return(format_output(results_estim_EM, phylo, tt))
}

estimateEM_wrapper_scratch <- function(phylo, Y_data,
                                       Y_data_imp,
                                       process, K_t,
                                       method.variance,
                                       random.root, stationnary.root,
                                       alpha_known, alpha,
                                       methods.segmentation,
                                       method.init,
                                       method.init.alpha,
                                       sBM_variance,
                                       method.OUsun,
                                       impute_init_Rphylopars,
                                       ...){
  tt <- system.time(results_estim_EM <- estimateEM(phylo = phylo, 
                                                   Y_data = Y_data, 
                                                   Y_data_imp = Y_data_imp,
                                                   process = process, 
                                                   method.variance = method.variance, 
                                                   method.init = method.init,
                                                   method.init.alpha = method.init.alpha,
                                                   nbr_of_shifts = K_t,
                                                   random.root = random.root,
                                                   stationnary.root = stationnary.root,
                                                   alpha_known = alpha_known,
                                                   known.selection.strength = alpha,
                                                   methods.segmentation = methods.segmentation,
                                                   sBM_variance = sBM_variance,
                                                   method.OUsun = method.OUsun,
                                                   impute_init_Rphylopars = impute_init_Rphylopars,
                                                   ...))
  return(format_output(results_estim_EM, phylo, tt))
}

PhyloEM_core <- function(phylo, Y_data, process = c("BM", "OU", "scOU", "rBM"),
                         K_max, use_previous = TRUE,
                         order = TRUE,
                         method.variance = "simple",
                         method.init = "default",
                         method.init.alpha = "default",
                         method.init.alpha.estimation = c("regression", "regression.MM", "median"), 
                         methods.segmentation = c("lasso", "best_single_move"),
                         alpha_known = TRUE,
                         random.root = FALSE,
                         stationnary.root = FALSE,
                         alpha = NULL,
                         check.tips.names = FALSE,
                         progress.bar = TRUE,
                         save_step = FALSE,
                         sBM_variance = FALSE,
                         method.OUsun = "rescale",
                    ...){
  ## Check the tree
  if (!is.ultrametric(phylo)) stop("The tree must be ultrametric.")
  ## Check that the vector of data is in the correct order and dimensions ################
  Y_data <- check_data(phylo, Y_data, check.tips.names)
  p <- nrow(Y_data)
  ntaxa <- length(phylo$tip.label)
  ## Process
  process <- match.arg(process)
  process_original <- process
  original_phy <- phylo
  temp <- choose_process_EM(process = process, p = p,
                            random.root = random.root,
                            stationnary.root = stationnary.root,
                            alpha_known = alpha_known,
                            known.selection.strength = alpha,
                            sBM_variance = sBM_variance,
                            method.OUsun = method.OUsun)
  
  rescale_tree <- temp$rescale_tree # Rescale the tree ?
  transform_scOU <- temp$transform_scOU # Re-transform parameters back ?
  sBM_variance <- temp$sBM_variance
  if (sBM_variance){ # process sBM : need a root branch.
    original_phy$root.edge <- 1
  }
  
  ## Compute alpha
  if (process == "BM") {
    alpha <- 0
  } else {
    alpha <- find_grid_alpha(phylo, alpha, ...)
  }
  ## Estimates
  X <- list(Y_data = Y_data,
            K_try = 0:K_max,
            ntaxa = ntaxa)
  
  ## Loop on alpha
  for (alp in alpha){
    ## Transform branch lengths if needed
    phylo <- original_phy
    if (rescale_tree) {
      phylo <- transform_branch_length(phylo, alp)
    }
    ## Fixed quantities
    times_shared <- compute_times_ca(phylo)
    distances_phylo <- compute_dist_phy(phylo)
    subtree.list <- enumerate_tips_under_edges(phylo)
    T_tree <- incidence.matrix(phylo)
    h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
    ## Impute data if needed
    Y_data_imp <- Y_data
    if (temp$process == "BM"){
      Y_data_imp <- impute.data.Rphylopars(phylo, Y_data, temp$process, random.init)
    }
    ## Estimations
    X <- Phylo_EM_sequencial(phylo = original_phy,
                             Y_data = Y_data,
                             Y_data_imp = Y_data_imp,
                             process = process,
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
                             stationnary.root = stationnary.root,
                             alp = alp,
                             check.tips.names = check.tips.names,
                             progress.bar = progress.bar,
                             times_shared = times_shared,
                             distances_phylo = distances_phylo,
                             subtree.list = subtree.list,
                             T_tree = T_tree,
                             h_tree = h_tree,
                             save_step = save_step,
                             sBM_variance = sBM_variance,
                             method.OUsun = method.OUsun,
                             ...)
  }
  ## Select max solution for each K
  X <- merge_max_grid_alpha(X, alpha)
  return(X)
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
  X <- format_output(results_estim_EM)
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

choose_process_EM <- function(process, p, random.root, stationnary.root, alpha_known,
                              known.selection.strength = 1, eps = 10^(-3),
                              sBM_variance = FALSE,
                              method.OUsun = "rescale"){
  ## Reduce Process
  transform_scOU <- FALSE # Should we re-transform back the parameters to get an OU ?
  rescale_tree <- FALSE # Should we re-scale the tree ?
  if (process == "OU"){
    if (p > 1) {
      warning("The general OU is not implemented in the multivariate case. Switching to the scalar OU (scOU).")
    }
    process <- "scOU"
  }
  if (process == "scOU"){
    if (random.root){
      if ((method.OUsun == "rescale")){
        if (!alpha_known){
          stop("The re-scaled scalar OU is only implemented for known selection strength. Please consider using a grid.")
        }
        sBM_variance <- TRUE
        if (!stationnary.root){
          warning("The scalar OU process with a general random root is not implemented for the multivariate case. The root is taken to be in the stationnary state.")
        }
        transform_scOU <- TRUE
        rescale_tree <- TRUE
        process <- "BM"
      } else {
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