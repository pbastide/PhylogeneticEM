# {Initialisations of the EM}
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
## Here are the functions to initialize the EM.
## Dependencies : generic_functions.R
## : simulate.R
## : shifts_manipulations.R
###############################################################################

######################################################
## Defaults initialisations
######################################################
##
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
##
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


#########################################################
## Lasso initializations
#########################################################

##
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
##
lasso_regression_K_fixed.glmnet <- function (Yp, Xp, K, intercept.penalty = FALSE ) {
  ## Penalty on the first coordinate = intercept : force first cooerdinate to be null
  excl <- NULL
  if (intercept.penalty) excl <- c(1)
  ## fit
  fit <- glmnet(x = 0 + Xp, y = Yp, alpha = 1, exclude = excl)
  ## Find the lambda that gives the right number of ruptures
  # Check that lambda goes far enought
  if (K > max(fit$df)) {
    fit <- glmnet(x = 0 + Xp, y = Yp, alpha = 1, lambda.min.ratio = 0, exclude = excl)
  }
  if (K > max(fit$df)) {
    stop("Lasso regression failed. There are too many variables.")
  }
  ## If the right lambda does not exists, find it.
  count <- 0
  while (sum(fit$df == K) == 0 && count < 500) {
    count <- count + 1
    K_inf <- K-1
    while ((sum(K_inf == fit$df) == 0) && (K_inf >= 0)) {
      K_inf <- K_inf - 1
    }
    lambda_inf <- fit$lambda[tail(which(K_inf == fit$df), n=1)]
    K_sup <- K + 1
    while ((sum(K_sup == fit$df) == 0) && (K_sup <= max(fit$df))) {
      K_sup <- K_sup + 1
    }
    lambda_sup <- fit$lambda[head(which(K_sup == fit$df), n=1)]
    lambda <- seq(from = lambda_inf, to = lambda_sup, length.out = 100)
    fit <- glmnet(x = 0 + Xp, y = Yp, alpha = 1, lambda = lambda, exclude = excl)
  }
  ## If the right lambda does not exists, raise the number of shifts
  K_2 <- K
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
    # If lm fails to find some coeeficients, put them to 0
    delta.gauss[is.na(delta.gauss)] <- 0.1
    ## If we had to raise the number of shifts, go back to the initial number, taking the K largest shifts
    edges <- order(-abs(delta.gauss))[1:K]
    delta.gauss.final <- rep(0, length(delta.gauss))
    delta.gauss.final[edges] <- delta.gauss[edges]
    shifts.gauss <- shifts.vector_to_list(delta.gauss.final);
    return(list(E0.gauss = E0.gauss, shifts.gauss = shifts.gauss))
  }
}

lasso_regression_K_fixed <- function (Yp, Xp, K, root = NULL) {
  ## Root is the intercept, should be excluded from varaiable selection
  # In that case, project Yp on the orthogonal of the root
  if (!is.null(root)){
    L <- Xp[,root]
    norme_L <- drop(crossprod(L))
    Xp_noroot <- Xp[ , -root, drop = FALSE]
    Xp_orth <- Xp_noroot - (tcrossprod(L)%*%Xp_noroot)/norme_L
    Yp_orth <- Yp - crossprod(Yp, L)/(norme_L)*L
    intercept <- FALSE
  } else {
    Xp_orth <- Xp
    Yp_orth <- Yp
    intercept <- TRUE
  }
  ## fit
  fit <- elastic.net(x = 0 + Xp_orth, y = Yp_orth, lambda2 = 0, nlambda1 = 500, intercept = intercept)
  df <- rowSums(fit@active.set)
  ## Find the lambda that gives the right number of ruptures
  # Check that lambda goes far enought
  if (K > max(df)) {
    fit <- elastic.net(x = 0 + Xp_orth, y = Yp_orth, lambda2 = 0, min.ratio = 10^(-10), intercept = intercept)
    df <- rowSums(fit@active.set)
  }
  if (K > max(df)) {
    stop("Lasso regression failed. There are too many variables.")
  }
  ## If the right lambda does not exists, find it.
  count <- 0
  while (sum(df == K) == 0 && count < 500) {
    count <- count + 1
    K_inf <- K-1
    while ((sum(K_inf == df) == 0) && (K_inf >= 0)) {
      K_inf <- K_inf - 1
    }
    lambda_inf <- fit@lambda1[tail(which(K_inf == df), n = 1)]
    K_sup <- K + 1
    while ((sum(K_sup == df) == 0) && (K_sup <= max(df))) {
      K_sup <- K_sup + 1
    }
    lambda_sup <- fit@lambda1[head(which(K_sup == df), n = 1)]
    lambda <- seq(from = lambda_inf, to = lambda_sup, length.out = 100)
    fit <- elastic.net(x = 0 + Xp_orth, y = Yp_orth, lambda1 = lambda, lambda2 = 0, intercept = intercept)
    df <- rowSums(fit@active.set)
  }
  ## If the right lambda does not exists, raise the number of shifts
  K_2 <- K
  while (sum(df == K_2) == 0 && K_2 < 500) {
    warning("During lasso regression, could not find the right lambda for the number of shifts K. Temporarly raised it to do the lasso regression, and furnishing the K largest coefficients.")
    K_2 <- K_2 + 1
  }
  ## If could not find the right lambda, do a default initialization
  if (sum(df == K_2) == 0) {
    stop("Lasso Initialisation failed : could not find a satisfying number of shifts.")
  }
  ## Select the row with the right number of coefficients
  index <- min(which(df == K_2))
  delta <- fit@coefficients[index,]
  # If we put aside the root, replace it in the coefficients
  if (!is.null(root)){
    #deltabis <- unname(c(delta[1:(root - 1)], 0, tail(delta, n = max(0, length(delta) - root + 1))))
    delta <- append(delta, 0, after = root - 1)
  }
  # Check that the matrix is of full rank
  projection <- which(delta != 0)
  Xproj <- Xp[ , projection, drop = FALSE]
  if (dim(Xproj)[2] != qr(Xproj)$rank) {
    stop("The selected variables do not produce a full rank regression matrix !")
  }
  ## If we had to raise the number of shifts, go back to the initial number, taking the K largest shifts
  edges <- order(-abs(delta))[1:K]
  delta.bis <- rep(0, length(delta))
  delta.bis[edges] <- delta[edges]
  ## Gauss lasso
  return(compute_gauss_lasso(Yp, Xp, delta.bis, root))
}

##
#' @title Do a lm on top of a lasso regression.
#'
#' @description
#' \code{compute_gauss_lasso} takes the variables selected by a lasso procedure, and
#' uses them to do a simple linear least square regression. Function used is
#' \code{lm} for non-transformed data (root = NULL), and \code{lm.fit} for
#' transformed data (root = an integer).
#'
#' @details
#' Depending on the value of root, the behaviour is different. If root is null, then
#' we fit a linear regression with an intercept. If root is equal to an integer,
#' then the "intercept" column of the matrix Xp (that has possibly been trough a 
#' multiplication with a cholesky decomposition of the variance) is included, rather
#' than the intercept.
#'
#' @param Yp (transformed) data
#' @param Xp (transformed) matrix of regression
#' @param delta regression coefficients obtained with a lasso regression
#' @param root the position of the root (intercept) in delta
#' 
#' @return E0.gauss the intercept (value at the root)
#' @return shifts.gauss the list of shifts found on the branches
#' @return residuals the residuals of the regression
#'
##
compute_gauss_lasso <- function (Yp, Xp, delta, root) {
  projection <- which(delta != 0)
  if (is.null(root)) { # If no one is excluded, "real" intercept
    Xproj <- 0 + Xp[, projection, drop = FALSE]
    fit.gauss <- lm(Yp ~ Xproj)
    delta.gauss <- rep(0, dim(Xp)[2])
    E0.gauss <- coef(fit.gauss)[1]; names(E0.gauss) <- NULL
    delta.gauss[projection] <- coef(fit.gauss)[-1]
  } else { # take intercept (root) into consideration
    Xproj <- 0 + Xp[, c(root, projection), drop = FALSE]
    fit.gauss <- lm.fit(x = Xproj, y = Yp)
    delta.gauss <- rep(0, dim(Xp)[2])
    E0.gauss <- coef(fit.gauss)[1]; names(E0.gauss) <- NULL
    delta.gauss[projection] <- coef(fit.gauss)[-1]
  }
  # If lm fails to find some coeficients
  if (anyNA(delta.gauss)) {
    warning("There were some NA in the lm fit for the Gauss Lasso. These were replaced with the values obtained from the Lasso. This is not the optimal solution.")
    delta.gauss[is.na(delta.gauss)] <- delta[is.na(delta.gauss)]
  }
  # Result
  shifts.gauss <- shifts.vector_to_list(delta.gauss);
  return(list(E0.gauss = E0.gauss, 
              shifts.gauss = shifts.gauss,
              residuals = residuals(fit.gauss)))
}

##
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
##
init.EM.lasso <- function(phylo, Y_data, process, times_shared = compute_times_ca(phylo), distances_phylo, nbr_of_shifts, use_sigma=TRUE, variance.init=1, random.init=TRUE, value.root.init=0, exp.root.init=1, var.root.init=1, edges.init=NULL, values.init=NULL, relativeTimes.init=NULL, selection.strength.init=1, optimal.value.init=0, T_tree = incidence.matrix(phylo), ...) {
  ntaxa <- length(phylo$tip.label)
  init.EM.default <- init.EM.default(process)
  ## Actualization of incidence matrix
  Tr <- T_tree
  ac_tree <- incidence_matrix_actualization_factors(tree = phylo, 
                                                    selection.strength = selection.strength.init,
                                                    times_shared = times_shared)
  Tr <- T_tree * ac_tree
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
    Tr <- cbind(Tr, rep(1, dim(Tr)[1])) # Here we use hypothesis : stationnary root.
    Tp <- Sig_chol_inv%*%Tr
    Yp <- Sig_chol_inv%*%Y_data
    fit <- try(lasso_regression_K_fixed(Yp = Yp, Xp = Tp, K = nbr_of_shifts, root = dim(Tr)[2]))
  } else {
    # Return untransformed Y_data and T
    Tp <- Tr
    Yp <- Y_data
    fit <- try(lasso_regression_K_fixed(Yp = Yp, Xp = Tp, K = nbr_of_shifts))
  }
  ## Fit
  if (inherits(fit, "try-error")) {
    warning("Lasso Initialisation fail : could not find a satisfying number of shifts. Proceeding to a default initialization.")
    return(init.EM.default(selection.strength.init = selection.strength.init, 
                           random.init = random.init, 
                           stationnary.root.init = stationnary.root.init, ...))
  } else { 
    E0.gauss <- fit$E0.gauss
    shifts.gauss <- fit$shifts.gauss
#     ## If OU, apply the correct factor to shifts
#     if (process == "OU" && !is.null(times_shared) && !is.null(shifts.gauss$values)){
#       parents <- phylo$edge[shifts.gauss$edges,1]
#       factors <- compute_actualization_factors(selection.strength = selection.strength.init, 
#                                                t_tree = t_tree, 
#                                                times_shared = times_shared, 
#                                                parents = parents)
#       shifts.gauss$values <- shifts.gauss$values/factors
#     }
    params_init <- init.EM.default(value.root.init = E0.gauss[1], 
                                   exp.root.init = E0.gauss[1], 
                                   optimal.value.init = E0.gauss[1],
                                   edges.init = shifts.gauss$edges, 
                                   values.init = shifts.gauss$values, 
                                   relativeTimes.init = shifts.gauss$relativeTimes, 
                                   selection.strength.init =selection.strength.init, 
                                   random.init = random.init, 
                                   var.root.init = var.root.init, ...)
    return(params_init)
  }
}

compute_actualization_factors <- function(selection.strength, 
                                                t_tree, 
                                                times_shared, 
                                                parents){
  factors <- exp(- selection.strength * (t_tree - diag(times_shared[parents, parents, drop = FALSE])))
  return(1-factors)
}

###############################################
## Robust initialisation of alpha (and gamma)
###############################################

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

##
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
##

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

init.alpha.gamma.estimation <- function(phylo, Y_data, nbr_of_shifts, distances_phylo, max_triplet_number, alpha_known, ...){
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
  if (alpha_known) {
    return(list(alpha_0 = init.alpha.gamma.default(alpha_known, ...)$alpha_0,
                gamma_0 = gamma_0))
  } else {
    df <- data.frame(square_diff=square_diff, dists=dists)
    fit.rob <- try(nlrob(square_diff ~ gam*(1-exp(-alpha*dists)), data=df, start=list(gam = gamma_0, alpha=1)))
    if (inherits(fit.rob, "try-error")) {
      warning("Robust estimation of alpha failed")
      return(list(alpha_0 = init.alpha.gamma.default(alpha_known, ...)$alpha_0,
                  gamma_0 = gamma_0))
    } else { 
      return(list(alpha_0 = unname(coef(fit.rob)["alpha"]), 
                  gamma_0 = gamma_0))
    }
  }
}

