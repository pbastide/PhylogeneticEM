# {M step}
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
## Here are the functions to compute the M step.
## Dependencies : generic_functions.R
## : simulate.R
## : shifts_manipulations.R
###############################################################################

##
# compute_M (phylo, process, Y_data, conditional_law_X, nbr_of_shifts)
# PARAMETERS:
#            @phylo (tree) input tree
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
##

# compute_M <- function(process) {
#   if (process=="BM"){
#     return(compute_M.BM)
#   }
# }

compute_M.BM <- function(phylo,
                         Y_data,
                         conditional_law_X,
                         nbr_of_shifts,
                         random.root,
                         variance_old,
                         mu_old,
                         sBM_variance, ...) {
  ## Initialization
  ntaxa <- length(phylo$tip.label)
  p <- nrow(Y_data)
  params <- init.EM.default.BM(random.init = random.root, p = nrow(Y_data), ...)
  ## Non-random root : diff_exp for children of root is just expectation 
  if (!random.root){
    conditional_law_X$expectations[, ntaxa + 1] <- rep(0, p)
  }
  ## Segmentation
  diff_exp <- compute_diff_exp.BM(phylo = phylo, 
                                  conditional_law_X = conditional_law_X)
  ## Deal with numerical imprecision
  diff_exp[abs(diff_exp) < .Machine$double.eps ^ 0.5] <- 0
  ## Transform diff_exp and lengths if root fixed
  trans <- compute_transformed_diff_exp_0.BM(phylo, random.root,
                                             diff_exp, mu_old)
  diff_exp_trans <- trans$diff_exp
  lengths_ed <- trans$lengths_ed
  costs0 <- compute_costs_0.simple(phylo = phylo,
                                   diff_exp = diff_exp_trans,
                                   variance_old = variance_old,
                                   lengths_ed = lengths_ed)
  seg <- segmentation.BM(nbr_of_shifts = nbr_of_shifts, 
                         costs0 = costs0,
                         diff_exp = diff_exp_trans)
  params$shifts <- check_dimensions.shifts(p, seg$shifts)
  edges_max <- seg$edges_max
  ## Variance and root parameters
  var_diff <- compute_var_diff.BM(phylo = phylo, 
                                  conditional_law_X = conditional_law_X)
  if (sBM_variance){
    ## Variance then root
    params$variance <- compute_var_M.sBM(phylo = phylo,
                                         var_diff = var_diff,
                                         diff_exp = diff_exp,
                                         edges_max = edges_max,
                                         conditional_law_X = conditional_law_X)
    params <- update.root.sBM(phylo, params, conditional_law_X, seg, ntaxa)
  } else {
    ## Root then variance
    params <- update.root.BM(phylo, params, conditional_law_X, seg, ntaxa)
    params$variance <- compute_var_M.BM(phylo = phylo,
                                        var_diff = var_diff,
                                        diff_exp = diff_exp,
                                        edges_max = edges_max,
                                        random.root = random.root,
                                        mu = params$root.state$value.root)
  }
  ## Variance
  return(params)
}

update.root.BM <- function (phylo, params, conditional_law_X, seg, ntaxa) {
  if (params$root.state$random) {
    params$root.state$value.root <- NA
    params$root.state$exp.root <- conditional_law_X$expectations[ , ntaxa + 1]
    params$root.state$var.root <- get_variance_node(ntaxa + 1,
                                                    conditional_law_X$variances)
  } else {
    params$root.state$value.root <- compute_root_value.BM(phylo,
                                                          conditional_law_X$expectations,
                                                          seg$shifts)
    params$root.state$exp.root <- NA
    params$root.state$var.root <- NA
  }
  return(params)
}

update.root.sBM <- function (phylo, params, conditional_law_X, seg, ntaxa) {
  if (!params$root.state$random) {
    stop("Something went wrong: in M step, root is said fixed, but trying to fit a sBM.")
    }
    params$root.state$value.root <- NA
    params$root.state$exp.root <- conditional_law_X$expectations[ , ntaxa + 1]
    params$root.state$var.root <- params$variance * phylo$root.edge
  return(params)
}

compute_costs_0.simple <- function(phylo, diff_exp, variance_old, lengths_ed){
  p <- nrow(diff_exp)
  if (p == 1){
    costs0 <- diff_exp^2
  } else {
    ## Compute inverse of variance
    R_chol <- chol(variance_old)
    R_chol_inv <- backsolve(R_chol, diag(ncol(R_chol)))
#    R_inv <- solve(variance_old)
    ## Compute Mahalanobis norms
    costs0 <- rep(NA, ncol(diff_exp))
    for (i in 1:ncol(diff_exp)){
      costs0[i] <- tcrossprod(t(diff_exp[, i]) %*% R_chol_inv)
#      costs0[i] <- t(diff_exp[, i]) %*% R_inv %*% diff_exp[, i]
    }
  }
  return(1/lengths_ed * costs0)
}

compute_transformed_diff_exp_0.BM <- function(phylo, random.root, diff_exp, mu_old){
  lengths_ed <- phylo$edge.length
  if(!random.root){
    ntaxa = length(phylo$tip.label)
    parents <- phylo$edge[,1]
    root_edges <- which(parents == ntaxa + 1)
    if (length(root_edges) == 2){ # Root has only two descendants
#      diff_exp[, root_edges[1]] <- conditional_law_X$expectations[ , daughters[root_edges[1]], drop = F] - conditional_law_X$expectations[ , daughters[root_edges[2]], drop = F]
      diff_exp[, root_edges[1]] <- diff_exp[, root_edges[1]] - diff_exp[, root_edges[2]]
      diff_exp[ , root_edges[2]] <- 0
      lengths_ed[root_edges[1]] <- sum(lengths_ed[root_edges])
    } else {
      diff_exp[, root_edges] <- diff_exp[, root_edges] - mu_old
    }
  }
  return(list(diff_exp = diff_exp,
              lengths_ed = lengths_ed))
}

compute_root_value.BM <- function(phylo,
                                  expectations,
                                  shifts){
  ntaxa = length(phylo$tip.label)
  p <- nrow(expectations)
  Nedge <- nrow(phylo$edge)
  daughters <- phylo$edge[ , 2]
  parents <- phylo$edge[ , 1]
  root_edges <- which(parents == ntaxa + 1)
  deltas <- matrix(0, p, Nedge)
  deltas <- shifts.list_to_matrix(phylo, shifts, p)
  expe_root <- expectations[ , daughters[root_edges], drop = F]
  shifts_root <- deltas[, root_edges, drop = F]
  if (all(shifts_root == 0) 
      && isTRUE(all.equal(expe_root[, 1], expe_root[, 2], tolerance = 1e-20))){
    # If no shift and identical, avoid numerical errors
    return(expe_root[, 1])
  } else if (isTRUE(all.equal(expe_root[, 1] - shifts_root [, 1],
                              expe_root[, 2]))){
    # If one shift in binary case fixed root
    return(expe_root[, 2])
  } else {
    mu <- rowSums(sweep((expe_root - shifts_root), 2, 1/(phylo$edge.length[root_edges]), '*'))
    mu <- mu / sum(1/phylo$edge.length[root_edges])
    return(mu)
  }
}

compute_M.OU <- function(stationary.root, shifts_at_nodes, alpha_known){
  if (stationary.root && shifts_at_nodes && alpha_known) {
    return(compute_M.OU.specialCase)
  } else if (stationary.root && shifts_at_nodes) {
    return(compute_M.OU.stationary.root_AND_shifts_at_nodes)
  } else {
    stop("The EM algorithm for the univariate raw OU is only defined (for the moment) for a stationary root and shifts at nodes !")
  }
}

compute_M.OU.specialCase <- function(phylo, Y_data, conditional_law_X,
                                     nbr_of_shifts, known.selection.strength,
                                     methods.segmentation, beta_0_old,
                                     shifts_old, subtree.list, params_old, ...){
  ## Initialization
  p <- nrow(Y_data)
  ntaxa <- length(phylo$tip.label)
  params <- init.EM.default.OU(selection.strength.init=known.selection.strength,
                               p = p,
                               stationary.root.init=TRUE,
                               random.init=TRUE,
                               ...)
  if (p > 1){
    params <- split_params_independent(params)
  }
  ## Comutation of regression matrix idoine
  D <- vector("list", p)
  Xp <- vector("list", p)
  for (l in 1:p){
    regMat <- compute_regression_matrices(phylo = phylo,
                                          conditional_law_X = conditional_law_X[[l]],
                                          selection.strength = known.selection.strength[l])
    D[[l]] <- regMat$D
    D[[l]] <- matrix(D[[l]], nrow = 1)
    Xp[[l]] <- regMat$Xp
  }
  ## Choose method(s) for segmentation
  segmentation.OU.specialCase <- function(method.segmentation){
    segmentation <- switch(method.segmentation, 
                           #max_costs_0 = segmentation.OU.specialCase.max_costs_0,
                           lasso = segmentation.OU.specialCase.lasso,
                           #lasso_one_move = segmentation.OU.specialCase.lasso_one_move,
                           same_shifts = segmentation.OU.specialCase.same_shifts,
                           #same_shifts_same_values = segmentation.OU.specialCase.same_shifts_same_values,
                           best_single_move = segmentation.OU.specialCase.best_single_move)
    return(segmentation(phylo = phylo, 
                        nbr_of_shifts = nbr_of_shifts, 
                        conditional_law_X = conditional_law_X, 
                        selection.strength = known.selection.strength,
                        beta_0_old = beta_0_old,
                        shifts_old = shifts_old,
                        D = D,
                        Xp = Xp,
                        params_old = params_old))
  }
  ## Segmentation
  segs <- sapply(methods.segmentation, segmentation.OU.specialCase,
                 simplify = FALSE)
  ## If everyone failed, keep the same shifts
  if (all(sapply(segs, function(z) (length(z$costs) > 0 && is.infinite(z$costs))
                                  || any(is.infinite(z[[1]]$costs))))){
    methods.segmentation <- "same_shifts"
    segs <- sapply(methods.segmentation, segmentation.OU.specialCase, simplify = FALSE)
  }
  #edges_max <- seg$edges_max
  ## Variance
  if (p == 1){
    var_diff <- compute_var_diff.OU(phylo = phylo, 
                                    conditional_law_X = conditional_law_X[[1]], 
                                    selection.strength = known.selection.strength)
    # Function to compute gamma2 for all parameters obtained by all segmentations
    compute_var_M <- function(method.segmentation){
      return(compute_var_M.OU.specialCase(phylo=phylo, 
                                          var_diff=var_diff, 
                                          costs=segs[[method.segmentation]][[1]]$costs, 
                                          selection.strength=known.selection.strength,
                                          conditional_root_variance=unname(conditional_law_X[[1]]$variances[ntaxa+1])))
    }
    var.roots <- sapply(methods.segmentation, compute_var_M)
    # cond log lik
    cond_exp_log_lik <- function(method.segmentation){
      return(conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes(phylo = phylo,
                                                                                                   conditional_law_X = conditional_law_X[[1]], 
                                                                                                   sigma2 = 2 * known.selection.strength * var.roots[method.segmentation],
                                                                                                   mu = segs[[method.segmentation]][[1]]$beta_0,
                                                                                                   shifts = segs[[method.segmentation]][[1]]$shifts,
                                                                                                   alpha = known.selection.strength))
    }
    obj_funcs <- sapply(methods.segmentation, cond_exp_log_lik)
    obj_funcs <- matrix(obj_funcs, ncol = length(methods.segmentation))
  } else {
    var_diff <- vector("list", p)
    var.roots <- vector("list", p)
    obj_funcs <- matrix(NA, nrow = p, ncol = length(methods.segmentation))
    for (l in 1:p){
      var_diff[[l]] <- compute_var_diff.OU(phylo = phylo, 
                                           conditional_law_X = conditional_law_X[[l]], 
                                           selection.strength = known.selection.strength[l])
      # Function to compute gamma2 for all parameters obtained by all segmentations
      compute_var_M <- function(method.segmentation){
        return(compute_var_M.OU.specialCase(phylo = phylo, 
                                            var_diff = var_diff[[l]], 
                                            costs = segs[[method.segmentation]][[l]]$costs, 
                                            selection.strength = known.selection.strength[l],
                                            conditional_root_variance = unname(conditional_law_X[[l]]$variances[ntaxa+1])))
      }
      var.roots[[l]] <- sapply(methods.segmentation, compute_var_M)
      # cond log lik
      cond_exp_log_lik <- function(method.segmentation){
        return(conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes(phylo = phylo,
                                                                                                     conditional_law_X = conditional_law_X[[l]], 
                                                                                                     sigma2 = 2 * known.selection.strength[l] * var.roots[[l]][method.segmentation],
                                                                                                     mu = segs[[method.segmentation]][[l]]$beta_0,
                                                                                                     shifts = segs[[method.segmentation]][[l]]$shifts,
                                                                                                     alpha = known.selection.strength[l]))
      }
      obj_funcs[l, ] <- sapply(methods.segmentation, cond_exp_log_lik)
    }
  }
  ## Compute objective function for each set of parameters, and choose the best one
  best.method.seg <- which.max(colSums(obj_funcs))
  sh_ed <- segs[[best.method.seg]]$shifts$edges
  if (p > 1) sh_ed <- segs[[best.method.seg]][[1]]$shifts$edges
  ## Take the best method, that provides a parsimonious solution
  while(!check_parsimony(phylo,
                         sh_ed,
                         subtree.list)
          && any(is.finite(obj_funcs))){
    obj_funcs[, best.method.seg] <- rep(-Inf, p)
    best.method.seg <- which.max(colSums(obj_funcs))
    sh_ed <- segs[[best.method.seg]]$shifts$edges
    if (p > 1) sh_ed <- segs[[best.method.seg]][[1]]$shifts$edges
  }
  ## If no solution is parsimonious, keep the same one.
  if (prod(is.infinite(obj_funcs)) == 1){
    warning("Could not find any parsimonious solution at the M Step. Keeping the same shifts.")
    if (is.null(shifts_old$edges)){
      warning("Had to use the same_shift method, but with no old shift. Taking random ones. This is probably not what you intented to do.")
      shifts_olds$edges <- sample_shifts_edges(phylo, nbr_of_shifts,
                                               part.list = subtree.list)
    }
    segs <- sapply("same_shifts", segmentation.OU.specialCase, simplify = FALSE)
    if (p == 1){
      var.roots <- sapply("same_shifts", compute_var_M)
    } else {
      for (l in 1:p){
        compute_var_M <- function(method.segmentation){
          return(compute_var_M.OU.specialCase(phylo = phylo, 
                                              var_diff = var_diff[[l]], 
                                              costs = segs[[method.segmentation]][[l]]$costs, 
                                              selection.strength = known.selection.strength[l],
                                              conditional_root_variance = unname(conditional_law_X[[l]]$variances[ntaxa+1])))
        }
        var.roots[[l]] <- sapply(methods.segmentation, compute_var_M)
      }
    }
    best.method.seg <- 1
  }
  ## Actualize paremters with the ones found by the best segmentation method
  if (p == 1){
    params$root.state$var.root <- unname(var.roots[best.method.seg])
    params$variance <- 2 * known.selection.strength * params$root.state$var.root
    params$shifts <- segs[[best.method.seg]][[1]]$shifts
    params$root.state$exp.root <- segs[[best.method.seg]][[1]]$beta_0
    params$optimal.value <- params$root.state$exp.root
    attr(params, "segmentation_algorithm_used") <- names(best.method.seg)
    ## Dimensions
    pp <- check_dimensions(1,
                           params$root.state,
                           params$shifts,
                           params$variance,
                           params$selection.strength,
                           params$optimal.value)
    params$root.state <- pp$root.state
    params$shifts <- pp$shifts
    params$variance <- pp$variance
    params$selection.strength <- pp$selection.strength
    params$optimal.value <- pp$optimal.value
    return(params)
  } else {
    for (l in 1:p){
      params[[l]]$root.state$var.root <- unname(var.roots[[l]][best.method.seg])
      params[[l]]$variance <- 2 * known.selection.strength[l] * params[[l]]$root.state$var.root
      params[[l]]$shifts <- segs[[best.method.seg]][[l]]$shifts
      params[[l]]$root.state$exp.root <- segs[[best.method.seg]][[l]]$beta_0
      params[[l]]$optimal.value <- params[[l]]$root.state$exp.root
      attr(params[[l]], "segmentation_algorithm_used") <- methods.segmentation[best.method.seg]
      ## Dimensions
      pp <- check_dimensions(1,
                             params[[l]]$root.state,
                             params[[l]]$shifts,
                             params[[l]]$variance,
                             params[[l]]$selection.strength,
                             params[[l]]$optimal.value)
      params[[l]]$root.state <- pp$root.state
      params[[l]]$shifts <- pp$shifts
      params[[l]]$variance <- pp$variance
      params[[l]]$selection.strength <- pp$selection.strength
      params[[l]]$optimal.value <- pp$optimal.value
    }
      return(params)
  }
}

compute_M.OU.stationary.root_AND_shifts_at_nodes <- function(phylo,
                                                    Y_data,
                                                    conditional_law_X,
                                                    nbr_of_shifts,
                                                    alpha_old,
                                                    max_selection.strength,
                                                    eps,
                                                    methods.segmentation,
                                                    beta_0_old,
                                                    shifts_old,
                                                    subtree.list,
                                                    params_old, ...){
  ## Estimate all parameters with alpha of the previous step
  params <- compute_M.OU.specialCase(phylo = phylo, 
                                     Y_data = Y_data, 
                                     conditional_law_X = conditional_law_X,
                                     nbr_of_shifts = nbr_of_shifts,
                                     known.selection.strength = alpha_old,
                                     methods.segmentation = methods.segmentation,
                                     beta_0_old = beta_0_old,
                                     shifts_old = shifts_old, 
                                     subtree.list = subtree.list,
                                     params_old = params_old)
  ## Estimate new alpha
  p <- nrow(Y_data)
  if (p > 1){
    for (l in 1:p){
      params[[l]]$selection.strength <- estimate.alpha(phylo = phylo,
                                                       conditional_law_X = conditional_law_X[[l]], 
                                                       sigma2 = params[[l]]$variance,
                                                       mu = params[[l]]$root.state$exp.root,
                                                       shifts = params[[l]]$shifts,
                                                       alpha_old = alpha_old[l],
                                                       max_selection.strength = max_selection.strength)
      ## Change value of the root (stationary) accordingly
      params[[l]]$variance <- 2 * params[[l]]$selection.strength * params[[l]]$root.state$var.root
    }
  } else {
    params$selection.strength <- estimate.alpha(phylo = phylo,
                                                conditional_law_X = conditional_law_X[[1]], 
                                                sigma2 = params$variance,
                                                mu = params$root.state$exp.root,
                                                shifts = params$shifts,
                                                alpha_old = alpha_old,
                                                max_selection.strength = max_selection.strength)
    params$selection.strength <- matrix(params$selection.strength, 1, 1)
    ## Change value of the root (stationary) accordingly
    params$variance <- 2 * params$selection.strength * params$root.state$var.root
  }
  return(params)
}


################################################################
## Some technical functions
################################################################

##
#' @title Compute differences of expectations between node and parent.
#'
#' @description
#' \code{compute_diff_exp} compute the differences of conditional expectations between all the
#' nodes and their parents.
#'
#' @param phylo a phylogenetic tree
#' @param conditional_law_X result of function \code{compute_E}
#' 
#' @return matrix p x Nedge containing, for each edge e finishing at node i,
#' the quantity E[Z_i|Y]-E[Z_pa(i)|Y].
#' 
#' @keywords internal
#'
##
compute_diff_exp.BM <- function(phylo, conditional_law_X) {
  diff_exp <- matrix(NA, nrow(conditional_law_X$expectations), nrow(phylo$edge))
  daughters <- phylo$edge[ , 2]
  parents <- phylo$edge[ , 1]
  diff_exp <- conditional_law_X$expectations[ , daughters, drop = F] - conditional_law_X$expectations[ , parents, drop = F]
  return(diff_exp)
}

compute_diff_exp.OU <- function(phylo, conditional_law_X, selection.strength) {
  diff_exp <- rep(NA, nrow(phylo$edge))
  daughters <- phylo$edge[,2]
  parents <- phylo$edge[,1]
  diff_exp <- conditional_law_X$expectations[daughters] - exp(-selection.strength * phylo$edge.length) * conditional_law_X$expectations[parents]
  # diff_exp <- conditional_law_X$expectations[, daughters, drop = F] - exp(-selection.strength * phylo$edge.length) * conditional_law_X$expectations[, parents, drop = F]
  return(diff_exp)
}

##
#' @title Compute variances of differences between nodes and parents.
#'
#' @description
#' \code{compute_var_diff} computes variances of differences between all the
#' nodes and their parents.
#'
#' @param phylo a phylogenetic tree
#' @param conditional_law_X result of function \code{compute_E}
#' 
#' @return array p x p x Nedge containing, for each edge e finishing at node i,
#' the quantity Var[Z_i-Z_pa(i)|Y].
#' 
#' @keywords internal
#'
##
compute_var_diff.BM <- function(phylo, conditional_law_X) {
  p <- nrow(conditional_law_X$expectations)
  Nedge <- nrow(phylo$edge)
  var_diff <- array(NA, c(p, p, Nedge))
  daughters <- phylo$edge[,2]
  parents <- phylo$edge[,1]
  for (e in 1:Nedge){
    # range_e <- ((e - 1) * p + 1):(e * p)
    dauvar <- get_variance_node(daughters[e], conditional_law_X$variances)
    parvar <- get_variance_node(parents[e], conditional_law_X$variances)
    daucov <- get_variance_node(daughters[e], conditional_law_X$covariances)
    var_diff[1:p, 1:p, e] <- dauvar + parvar - daucov - t(daucov)
    # tmp <- dauvar + parvar - daucov - t(daucov)
#     if (!isSymmetric(tmp, tol = 10000 * .Machine$double.eps)){
#       stop("var_diff[", e, "] should be symmetric.")
#     }
    # var_diff[1:p, 1:p, e] <- as.matrix(forceSymmetric(tmp))
  }
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

##
#' @title Compute weighted sum of var_diff
#'
#' @description
#' \code{compute_sum_var_diff} computes sum_{e edge} ell_j * Var[X_j - X_pa(j) | Y]
#'
#' @param phylo a phylogenetic tree
#' @param var_diff result of function \code{compute_var_diff.BM}
#' 
#' @return matrix p x p
#' 
#' @keywords internal
#'
##
compute_sum_var_diff <- function(phylo, var_diff){
  p <- dim(var_diff)[1]
  if (p == 1){
    return(Matrix(sum(var_diff * 1/phylo$edge.length)))
  } else {
    Nedge <- dim(var_diff)[3]
    vv <- sweep(var_diff, MARGIN = 3, STATS = 1/phylo$edge.length,
                FUN = '*', check.margin = FALSE)
    # vv <- var_diff %*% diag(1/rep(phylo$edge.length, each = p)) # mult each column by length
    # arr <- array(vv, dim = c(p, p, Nedge))
    res <- apply(vv, 1, rowSums)
    if (!isSymmetric(res, tol = .Machine$double.eps^0.5)){
      stop("Sum of variances should be symmetric. It is not.")
    }
    return(forceSymmetric(res))
  }
}

##
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
#' @param var_diff variances of differences result of function \code{compute_var_diff.BM}
#' @param diff_exp differences of expectations result of function \code{compute_diff_exp.BM}
#' @param edges_max Edges where the shifts occur result of function \code{segmentation.BM}
#' 
#' @return a p x p matrix : the computed variance
#' 
#' @keywords internal
##
compute_var_M.BM <- function(phylo, var_diff, diff_exp, edges_max, random.root, mu){
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  p <- nrow(var_diff)
  varr <- compute_sum_var_diff(phylo, var_diff)
  if (!random.root){ ## If root not random, substract root value
    root_edges <- which(phylo$edge[,1] == ntaxa + 1)
    diff_exp[, root_edges] <- diff_exp[, root_edges] - mu
  }
  expp <- as(tcrossprod(sweep(diff_exp[, -edges_max, drop = F], 2,
                              sqrt(1/phylo$edge.length[-edges_max]), '*')),
             "symmetricMatrix")
  return(1/(ntaxa + Nnode - 1) * (varr + expp))
}

compute_var_M.sBM <- function(phylo, var_diff, diff_exp, edges_max, conditional_law_X){
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  p <- nrow(var_diff)
  varr <- compute_sum_var_diff(phylo, var_diff)
  varr <- varr + 1 / phylo$root.edge * get_variance_node(ntaxa + 1,
                                   conditional_law_X$variances)
  expp <- as(tcrossprod(sweep(diff_exp[, -edges_max, drop = F], 2,
                              sqrt(1/phylo$edge.length[-edges_max]), '*')),
             "symmetricMatrix")
  # expp <- expp + 1 / phylo$root.edge * tcrossprod(conditional_law_X$expectations[ , ntaxa + 1])
  return(1/(ntaxa + Nnode) * (varr + expp))
}

compute_var_M.OU.specialCase <- function(phylo, var_diff, costs, selection.strength, conditional_root_variance){
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  ee <- exp(- selection.strength * phylo$edge.length )
  #return(1/(ntaxa + Nnode) * ( conditional_root_variance + sum((1 - ee^2)^(-1) * var_diff ) + sum( costs0[-edges_max]) ))
  return(1/(ntaxa + Nnode) * (conditional_root_variance + sum((1 - ee^2)^(-1) * var_diff) + sum(costs)))
}

################################################################
## Estimation of alpha
################################################################

##
#' @title Function to estimate alpha
#'
#' @description
#' \code{optimize} to maximize the
#' conditional expectation log likelihood in alpha. The interval is
#' set to [alpha_old/2, 2*alpha_old], supposing that the previous guess of 
#' alpha_old is not far from reality.
#'
#' @details
#' This function uses functions \code{compute_var_diff.OU} 
#' and \code{compute_diff_exp.OU} in the process. Careful : only works if the
#' root is stationary, and shifts at nodes.
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
#' @keywords internal
#'
#09/07/14 - Initial release
#02/10/14 - Take newly estimated shifts into consideration
##
estimate.alpha <- function(phylo,
                           conditional_law_X,
                           sigma2,
                           mu,
                           shifts,
                           alpha_old,
                           max_selection.strength){
  #opt <- optimize(R_function, phylo=phylo, conditional_law_X=conditional_law_X, sigma2=sigma2, mu=mu, interval=c(alpha_old/2, alpha_old*2))
  betas <- compute_betas_from_shifts(phylo = phylo, optimal.value = mu, shifts = shifts)
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

#########################################################
## Objective function computation
#########################################################

# conditional_expectation_log_likelihood.OU <- function(stationary.root, shifts_at_nodes, alpha_known){
#   if (stationary.root && shifts_at_nodes) {
#     return(conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes)
#   } else {
#     stop("The EM algorithm for the OU is only defined (for the moment) for a stationary root and shifts at nodes !")
#   }
# }

# #' ##@title Expectation conditional to the tips of the completed log-likelihood.
# #' This is an old version.
# #' 
# #' @description
# #' \code{conditional_expectation_log_likelihood.OU} computes the expectation
# #'  conditional to the tips of the completed log-likelihood, given the
# #'  first and second order moments of the nodes, and the parameters of the
# #'  OU process.
# #' 
# #' @details
# #' This function uses functions \code{compute_var_diff.OU}
# #' and \code{compute_diff_exp.OU} in the process. Careful : only works if the
# #' root is stationary, and shifts at nodes.
# #' 
# #' @param phylo Input tree.
# #' @param conditional_law_X result of function \code{compute_E.OU}, containing
# #' first and second moments.
# #' @param sigma2 variance of the OU.
# #' @param mu mean of the root.
# #' @param alpha selection strength of the OU.
# #' 
# #' @return double : value of the function
# #' 
# #' @keywords internal
# #' 
# #' #09/07/14 - Initial releas
# ##
# conditional_expectation_log_likelihood.OU.OLD <- function(phylo, conditional_law_X, sigma2, mu, alpha){
#   ntaxa <- length(phylo$tip.label)
#   Nnode <- phylo$Nnode
#   ## Constante
#   cst <- -1/2 * ((Nnode + ntaxa) * log(2*pi))
#   LogLik <- cst
#   ## Terms in log
#   ee <- exp(- alpha * phylo$edge.length )
#   LogLik <- LogLik - (Nnode + ntaxa - 1)*log(sigma2/(2*alpha))/2 - sum(log(1-ee^2))/2
#   ## Terms with the variance
#   K_1 <- conditional_law_X$variances[ntaxa+1] + (conditional_law_X$expectations[ntaxa+1] - mu)^2
#   var_diff <- compute_var_diff.OU(phylo=phylo,
#                                   conditional_law_X=conditional_law_X,
#                                   selection.strength=alpha)
#   LogLik <- LogLik - (alpha / sigma2) * (K_1 + sum((1 - ee^2)^(-1) * var_diff))
#   ## Terms with the expectation
#   diff_exp <- compute_diff_exp.OU(phylo=phylo,
#                                   conditional_law_X=conditional_law_X,
#                                   selection.strength=alpha)
#   daughters <- phylo$edge[,2]
#   betas <- conditional_law_X$optimal.values[daughters]
#   LogLik <- LogLik - (alpha / sigma2) * sum((1 - ee^2)^(-1) * (diff_exp - betas * (1-ee))^2)
#   return(LogLik)
# }

###
## @title Expectation conditional to the tips of the completed log-likelihood.
##
## @description
## \code{conditional_expectation_log_likelihood.OU} computes the expectation
##  conditional to the tips of the completed log-likelihood, given the 
##  first and second order moments of the nodes, and the parameters of the 
##  OU process.
##
## @details
## This function uses functions \code{compute_var_diff.OU} 
## and \code{compute_diff_exp.OU} in the process. Careful : only works if the
## root is stationary, and shifts at nodes.
##
## @param phylo Input tree.
## @param conditional_law_X result of function \code{compute_E.OU}, containing
## first and second moments.
## @param sigma2 variance of the OU.
## @param mu mean of the root.
## @param alpha selection strength of the OU.
## 
## @return double : value of the function
## 
## @keywords internal
##
##09/07/14 - Initial release
##02/10/14 - take new shifts in consideration
###
# 
# # conditional_expectation_log_likelihood.OU.stationary_root_shifts_at_nodes <- function(phylo, conditional_law_X, sigma2, mu, shifts, alpha){
# #   ntaxa <- length(phylo$tip.label)
# #   Nnode <- phylo$Nnode
# #   ## Constante
# #   cst <- -1/2 * ((Nnode + ntaxa) * log(2*pi))
# #   LogLik <- cst
# #   ## Terms in log
# #   ee <- exp(- alpha * phylo$edge.length )
# #   LogLik <- LogLik - (Nnode + ntaxa - 1)*log(sigma2/(2*alpha))/2 - sum(log(1-ee^2))/2
# #   ## Terms with the variance
# #   K_1 <- conditional_law_X$variances[ntaxa+1] + (conditional_law_X$expectations[ntaxa+1] - mu)^2
# #   var_diff <- compute_var_diff.OU(phylo=phylo, 
# #                                   conditional_law_X=conditional_law_X, 
# #                                   selection.strength=alpha)
# #   LogLik <- LogLik - (alpha / sigma2) * (K_1 + sum((1 - ee^2)^(-1) * var_diff))
# #   ## Terms with the expectation
# #   diff_exp <- compute_diff_exp.OU(phylo=phylo, 
# #                                   conditional_law_X=conditional_law_X, 
# #                                   selection.strength=alpha)
# #   parents <- phylo$edge[,1]
# #   daughters <- phylo$edge[,2]
# #   betas <- conditional_law_X$optimal.values[parents]
# #   delta <- shifts.list_to_vector(phylo, shifts) # vector of shifts (branches)
# #   LogLik <- LogLik - (alpha / sigma2) * sum((1 - ee^2)^(-1) * (diff_exp - (betas + delta) * (1-ee))^2)
# #   return(LogLik)
# # }

conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes <- function(phylo, conditional_law_X, sigma2, mu, shifts, alpha){
  betas <- compute_betas_from_shifts(phylo = phylo, optimal.value = mu, shifts = shifts)
  return(conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes.from_betas(phylo, conditional_law_X, sigma2, mu, betas, alpha))
}

conditional_expectation_log_likelihood_real_shifts.OU.stationary_root_shifts_at_nodes.from_betas <- function(phylo, conditional_law_X, sigma2, mu, betas, alpha){
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  ## Constante
  cst <- -1/2 * ((Nnode + ntaxa) * log(2*pi))
  LogLik <- cst
  ## Terms in log
  ee <- exp(- alpha * phylo$edge.length )
  LogLik <- LogLik - (Nnode + ntaxa)*log(sigma2/(2*alpha))/2 - sum(log(1-ee^2))/2
  ## Terms with the variance
  K_1 <- as.vector(conditional_law_X$variances)[ntaxa+1] + (as.vector(conditional_law_X$expectations)[ntaxa+1] - mu)^2
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
  return(unname(as.vector(LogLik)))
}

#######################################################################
## Segmentation algorithms
#######################################################################

##
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
#' @keywords internal
#02/06/14 - Initial release
##
segmentation.BM <- function(nbr_of_shifts, costs0, diff_exp){
  if (nbr_of_shifts > 0) {
    edges_max <- order(-costs0)[1:nbr_of_shifts]
    edges <- edges_max
    values <- diff_exp[ , edges_max, drop = F]
    relativeTimes <- rep(0,length(edges_max))
    shifts <- list(edges = edges,
                   values = values,
                   relativeTimes = relativeTimes)
    return(list(edges_max = edges_max,
                shifts = shifts))
  } else {
    edges_max <- length(costs0) + 1
    return(list(edges_max = edges_max,
                shifts = NULL))
  }
}

# ##
# #' @title Segmentation in the OU special case, algo 1
# #'
# #' @description
# #' \code{segmentation.OU.specialCase.max_costs_0} performs the first segmentation 
# #' algorithm described.
# #'
# #' @details
# #'  This function takes the largest values of costs0, and make them null,
# #'  thanks to delta and tau.
# #'
# #' @param phylo a phylogenetic tree
# #' @param nbr_of_shifts Number of shifts on the phylogeny allowed
# #' @param conditional_law_X moments of the conditional law of X given Y, result
# #' of function \code{compute_M.OU.specialCase}
# #' @param selection.strength the selection strength
# #' 
# #' @return List containing : beta_0 : the optimal value at the root
# #'                           shifts : list containing the computed tau and delta
# #'                           costs : vector of costs
# #'
# #' @keywords internal
# #10/06/14 - Initial release
# #06/10/14 - Change name to include other algorithms
# ##
# segmentation.OU.specialCase.max_costs_0 <- function(phylo, nbr_of_shifts,
#                                                     conditional_law_X,
#                                                     selection.strength, ...){
#   ntaxa <- length(phylo$tip.label)
#   ## Computation of mu=beta0
#   beta_0 <- conditional_law_X$expectations[ntaxa + 1]
#   ## Computation of costs
#   diff_exp <- compute_diff_exp.OU(phylo = phylo, 
#                                   conditional_law_X = conditional_law_X, 
#                                   selection.strength = selection.strength)
#   parents <- phylo$edge[, 1]
#   betas <- conditional_law_X$optimal.values
#   ee <- exp(- selection.strength * phylo$edge.length )
#   costs0 <- (1 - ee^2)^(-1) * (diff_exp - betas[parents] * (1-ee))^2
#   ## Segmentation per se
#   if (nbr_of_shifts > 0) {
#     # Max of costs
#     edges_max <- order(-costs0)[1:nbr_of_shifts]
#     edges <- edges_max
#     # Put them to 0
#     ee <- exp(- selection.strength * phylo$edge.length )
#     parents <- phylo$edge[edges_max,1]
#     values <- (1 - ee[edges_max])^(-1) * diff_exp[edges_max] - betas[parents]
#     relativeTimes <- rep(0,length(edges_max))
#     # Define shifts
#     shifts <- list(edges = edges, values = values, relativeTimes = relativeTimes)
#     # Compute new costs
#     costs <- costs0
#     costs[edges] <- 0
#     return(list(beta_0 = beta_0, edges_max = edges_max, shifts = shifts, costs = costs))
#   } else {
#     edges_max <- length(costs0) + 1
#     return(list(beta_0 = beta_0, edges_max = edges_max, shifts = NULL, costs = costs0))
#   }
# }

##
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
#' @keywords internal                           
#06/10/14 - Initial release
##
segmentation.OU.specialCase.lasso <- function(phylo, nbr_of_shifts, D, Xp, 
                                              penscale = rep(1, (nrow(phylo$edge) + 1)), ...){
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  ## Computation of answer matrix D : already done by now.
  if(is.list(D)){ # p independent traits
    p <- length(D)
    Dvec <- as.vector(do.call(cbind, D))
    Xkro <- bdiag(Xp)
    ## Re-order by branches
    traittoedge <- as.vector(sapply(1:((nrow(phylo$edge) + 1)),
                                    function(z) ((0:(p-1)) * (nrow(phylo$edge) + 1) + z)))
    Dvec <- Dvec[traittoedge]
    Xkro <- as.matrix(Xkro[traittoedge, traittoedge])
    ## Segmentation per se
    group <- rep(1:(ntaxa + Nnode), each = length(D))
    # Lasso regression
    fit <- try(lasso_regression_K_fixed.gglasso(Yvec = Dvec, Xkro = Xkro,
                                                K = nbr_of_shifts,
                                                root = ntaxa + Nnode,
                                                penscale = penscale,
                                                group = group,
                                                p_dim = p))
  } else { # Only one trait
    fit <- try(lasso_regression_K_fixed.glmnet_multivariate(Yp = D, Xp = Xp,
                                                            K = nbr_of_shifts,
                                                            root = ntaxa+Nnode,
                                                            penscale = penscale))
  }
  if (inherits(fit, "try-error")) {
    warning("At M step, Lasso regression failed.")
    if (is.list(D)){
      p = length(D)
      ret <- vector("list", p)
      for (l in 1:p){
        ret[[l]] <- list(beta_0 = 0, shifts = NULL, costs = Inf)
      }
      return(ret)
    } else {
      return(list(beta_0 = 0, shifts = NULL, costs = Inf))
    }
  } else {
    # Define shifts
    shifts <- shifts.matrix_to_list(t(fit$delta.gauss))
    # Define mu = beta_0
    beta_0 <- fit$E0.gauss
    # Compute new costs
    costs <- fit$residuals^2
    if (is.list(D)){ # Split independent traits
      p <- length(D)
      # re-order costs
      edgetotrait <- as.vector(sapply(1:p,
                                      function(z) ((0:(nrow(phylo$edge))) * p + z)))
      costs <- costs[edgetotrait]
      ret <- vector("list", p)
      for (l in 1:p){
        ret[[l]] <- list(beta_0 = beta_0[l],
                         shifts = list(edges = shifts$edges,
                                       values = shifts$values[l, ],
                                       relativeTimes = shifts$relativeTimes),
                         costs = costs[(l-1) * (nrow(phylo$edge) + 1) + 1:(nrow(phylo$edge) + 1)])
      }
      return(ret)
    } else { # One single trait
      return(list(beta_0 = beta_0, shifts = shifts, costs = costs))
    }
  }
}

compute_regression_matrices <- function(phylo, conditional_law_X, selection.strength, ...){
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  ## Computation of answer matrix D
  diff_exp <- compute_diff_exp.OU(phylo = phylo, 
                                  conditional_law_X = conditional_law_X, 
                                  selection.strength = selection.strength)
  daughters <- phylo$edge[,2]
  ee <- exp(- selection.strength * phylo$edge.length )
  D <- (1 - ee^2)^(-1/2) * diff_exp
  D <- D[match(1:(ntaxa + Nnode), phylo$edge[,2])]
  D[ntaxa + 1] <- conditional_law_X$expectations[ntaxa+1]
  ## Regression matrix : modified incidence matrix
  U <- incidence.matrix.full(phylo)
  U <- cbind(U, rep(1, dim(U)[1]))
  #A <- diag(c(sqrt((1-ee)/(1+ee)),1))
  A <- sqrt((1-ee)/(1+ee))
  A <- A[match(1:(ntaxa + Nnode), phylo$edge[,2])]
  A[ntaxa + 1] <- 1
  A <- diag(A)
  Xp <- A%*%U
  return(list(D = D, Xp = Xp))
}

# segmentation.OU.specialCase.lasso_one_move <- function(phylo, shifts_old, nbr_of_shifts, D, Xp, ...){
#   ## If no shifts, there is no such thing as a "single move"
#   if (is.null(shifts_old$edges)){
#     return(list(beta_0 = 0, shifts = shifts_old, costs = Inf))
#   }
#   ntaxa <- length(phylo$tip.label)
#   Nnode <- phylo$Nnode
#   ## do not penalise K-1 shifts
#   pens <- rep(1, ncol(Xp))
#   penscales <- matrix(1, nrow = nbr_of_shifts, ncol = ncol(Xp))
#   for (i in 1:nbr_of_shifts){
#     shs <- shifts_old$edges[-i]
#     penscales[i, shs] <- 0
#   }
#   ## Computation of answer matrix D : already done by now.
#   ## Segmentation per se
#   # Lasso regressions
#   fun <- function(penscale){
#     segmentation.OU.specialCase.lasso(phylo = phylo, 
#                                       nbr_of_shifts = nbr_of_shifts,
#                                       D = D,
#                                       Xp = Xp,
#                                       penscale = penscale)
#   }
#   segs <- apply(penscales, 1, fun)
#   costs <- sapply(segs, function(z) return(sum(z$cost)))
#   return(segs[[which.min(costs)]])
# }

# ##
# #' @title Segmentation in the OU special case, conserving the same shifts.
# #'
# #' @description
# #' \code{segmentation.OU.specialCase.same_shifts_same_values} keeps the same
# #' parameters and compute the quantities needed. It is here to ensure that we do
# #' at least better than the previous step.
# #'
# #' @details
# #'  This function takes the old shifts parameters, and compute costs with the new
# #'  moments.
# #'
# #' @param phylo a phylogenetic tree
# #' @param conditional_law_X moments of the conditional law of X given Y, result
# #' of function \code{compute_M.OU.specialCase}
# #' @param selection.strength the selection strength
# #' @param beta_0_old the previous (and concerved) value of beta_0
# #' @param shifts_old the previous (and concerved) list of shifts
# #' 
# #' @return List containing : beta_0 : the optimal value at the root
# #'                           shifts : list containing the computed tau and delta
# #'                           costs : vector of costs
# #'
# #' @keywords internal
# #' 
# #06/10/14 - Initial release
# ##
# segmentation.OU.specialCase.same_shifts_same_values <- function(phylo, conditional_law_X, selection.strength, beta_0_old, shifts_old, ...){
#   ntaxa <- length(phylo$tip.label)
#   browser()
#   ## Computation of costs
#   diff_exp <- compute_diff_exp.OU(phylo = phylo, 
#                                   conditional_law_X = conditional_law_X, 
#                                   selection.strength = selection.strength)
#   betas <- compute_betas_from_shifts(phylo = phylo, optimal.value = beta_0_old, shifts = shifts_old)
#   daughters <- phylo$edge[,2]
#   ee <- exp(- selection.strength * phylo$edge.length )
#   costs <- (1 - ee^2)^(-1) * (diff_exp - betas[daughters] * (1-ee))^2
#   costs <- c((conditional_law_X$expectations[ntaxa+1] - beta_0_old)^2, costs)
#   return(list(beta_0 = beta_0_old, shifts = shifts_old, costs = costs))
# }

##
#' @title Segmentation in the OU special case, conserving the same shifts
#'  position.
#'
#' @description
#' \code{segmentation.OU.specialCase.same_shifts} keeps the same shifts position,
#' and optimize the sum of costs using function 
#' \code{best_scenario}
# \code{optimize_costs_given_shift_position.OU.specialCase}. 
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
#' @keywords internal
#06/10/14 - Initial release
##
segmentation.OU.specialCase.same_shifts <- function(phylo, shifts_old, D, Xp, ...){
  edges <- shifts_old$edges
  if (is.null(edges)) edges <- 0
  dim(edges) <- c(1, length(edges))
  root <- length(phylo$tip.label) + phylo$Nnode
  return(best_scenario(as.vector(D), Xp, edges, root))
#   return(optimize_costs_given_shift_position.OU.specialCase(phylo = phylo,
#                                                             conditional_law_X = conditional_law_X,
#                                                             selection.strength = selection.strength,
#                                                             shifts_edges = shifts_old$edges))
}

# segmentation.OU.specialCase.best_single_move.old <- function(phylo, conditional_law_X, selection.strength, shifts_old, ...){
#   # Construct vector of all allowed combinations of shifts (variation from a base scenario)
#   shifts_edges <- shifts_old$edges
#   K <- length(shifts_edges)
#   Nedge <- dim(phylo$edge)[1]
#   allowed_moves <- which(!(1:Nedge %in% shifts_edges))
#   scenarii <- t(matrix(shifts_edges, (Nedge - K) * K + 1, nrow = K))
#   for (i in 1:K) {
#     scenarii[1 + ((i-1)*(Nedge - K)+1):(i*(Nedge - K)), i] <- allowed_moves
#   }
#   # Function to be applyed to each row
#   fun <- function(sh_ed){
#     seg <- optimize_costs_given_shift_position.OU.specialCase(phylo = phylo,
#                                                               conditional_law_X = conditional_law_X,
#                                                               selection.strength = selection.strength,
#                                                               shifts_edges = sh_ed)
#     totalCost <- sum(seg$costs)
#     return(list(seg = seg, 
#                 totalCost = totalCost))
#   }
#   # Apply to each row and take the minimal total cost
#   allSegs <- apply(scenarii, 1, fun)
#   dd <- do.call(rbind, allSegs)
#   min_conf <- which.min(dd[,"totalCost"])
#   return(dd[[min_conf]])
# }

##
# @title Minimization of the sum of costs, given the shift position.
#
# @description
# \code{optimize_costs_given_shift_position.OU.specialCase} minimize the sum of
# costs when the shift position is fixed.
#
# @details
# This function find the regimes of each node using function 
# \code{allocate_regimes_from_shifts} and optimize the sum of costs, computed 
# using function \code{compute_diff_exp.OU}, in the values of the optimal values
# betas (using a close formula). It then goes back to a shift expression of the
# problem using function \code{compute_shifts_from_betas}.
#
# @param phylo a phylogenetic tree
# @param conditional_law_X moments of the conditional law of X given Y, result
# of function \code{compute_M.OU.specialCase}
# @param selection.strength the selection strength
# @param shifts_edges the vector of the position of the shifts on the tree
# 
# @return List containing : beta_0 : the optimal value at the root
#                           shifts : list containing the computed tau and delta
#                           costs : vector of costs
#
# @keywords internal
# 
#15/10/14 - Initial release
##
# optimize_costs_given_shift_position.OU.specialCase <- function(phylo, conditional_law_X, selection.strength, shifts_edges, ...){
#   ntaxa <- length(phylo$tip.label)
#   ## Computation values of betas
#   diff_exp <- compute_diff_exp.OU(phylo = phylo, 
#                                   conditional_law_X = conditional_law_X, 
#                                   selection.strength = selection.strength)
#   # Regimes of the branches from alod positions of shifts
#   regimes <- allocate_regimes_from_shifts(phylo, shifts_edges)
#   # Quantities needed
#   daughters <- phylo$edge[,2]
#   ee <- exp(- selection.strength * phylo$edge.length )
#   numerateur <- diff_exp / (1 + ee)
#   denominateur <- (1 - ee) / (1 + ee)
#   # Root branch
#   Nedge <- dim(phylo$edge)[1]
#   numerateur[Nedge + 1]  <- conditional_law_X$expectations[ntaxa+1]
#   denominateur[Nedge + 1]  <- 1
#   # Function to compute betas on the regimes
#   fun <- function(reg){
#     edg <- which(phylo$edge[,2] %in% which(regimes == reg))
#     if (reg == 0) edg <- c(edg, Nedge+1)
#     return((sum(numerateur[edg])) / (sum(denominateur[edg])))
#   }
#   regimes_values <- 0:(length(unique(regimes))-1)
#   beta_values <- sapply(regimes_values, fun)
#   ## Computation of corresponding shifts values
#   betas <- beta_values[regimes + 1]
#   shifts <- compute_shifts_from_betas(phylo, betas)
#   ## Computation of costs
#   costs <- (1 - ee^2)^(-1) * (diff_exp - betas[daughters] * (1-ee))^2
#   costs <- c((conditional_law_X$expectations[ntaxa+1] - beta_values[1])^2, costs)
#   return(list(beta_0 = beta_values[1], shifts = shifts, costs = costs))
# }

segmentation.OU.specialCase.best_single_move <- function(phylo, shifts_old, D, Xp, params_old, ...){
  if (is.list(D)){ # case p > 1, independent
    p <- length(D)
    shifts_old <- params_old[[1]]$shifts
  }
  ## If no shifts, there is no such thing as a "single move"
  if (is.null(shifts_old$edges)){
    if (is.list(D)){
      p = length(D)
      ret <- vector("list", p)
      for (l in 1:p){
        ret[[l]] <- list(beta_0 = 0, shifts = NULL, costs = Inf)
      }
      return(ret)
    } else {
      return(list(beta_0 = 0, shifts = shifts_old, costs = Inf))
    }
  }
  ## Construct scenarii
  ntaxa <- length(phylo$tip.label)
  Nnode <- phylo$Nnode
  root <- ntaxa + Nnode
  # Construct vector of all allowed combinations of shifts (variation from a base scenario)
  shifts_edges <- shifts_old$edges
  K <- length(shifts_edges)
  Nedge <- dim(phylo$edge)[1]
  allowed_moves <- which(!(1:Nedge %in% shifts_edges))
  scenarii <- t(matrix(shifts_edges, (Nedge - K) * K + 1, nrow = K))
  for (i in 1:K) {
    scenarii[1 + ((i-1)*(Nedge - K)+1):(i*(Nedge - K)), i] <- allowed_moves
  }
  ## Choose the best scenario
  return(best_scenario(D, Xp, scenarii, root))
}

test_all_scenarii <- function(D, Xp, scenarii, root){
  # Function to be applyed to each row
  fun <- function(sh_ed){
    fit.lm <- lm.fit(x = Xp[, c(root, sh_ed), drop = FALSE], y = D)
    squared_res <- residuals(fit.lm)^2
    return(list(shifts_edges = sh_ed,
                coefs = unname(coef(fit.lm)),
                costs = squared_res,
                totalCost = sum(squared_res)))
  }
  # Apply to each row and take the minimal total cost
  allSegs <- apply(scenarii, 1, fun)
  return(allSegs)
}

best_scenario <- function (D, Xp, scenarii, root) {
  if (!is.list(D)){
    allSegs <- test_all_scenarii(as.vector(D), Xp, scenarii, root)
    dd <- do.call(rbind, allSegs)
    min_conf <- which.min(dd[,"totalCost"])
    # Go back to good format
    res <- dd[min_conf,]
    delta <- rep(0, dim(Xp)[2])
    delta[res$shifts_edges] <- res$coef[-1]
    shifts <- shifts.vector_to_list(delta)
    return(list(beta_0 = res$coef[1], shifts = shifts, costs = res$costs))
  } else {
    p <- length(D)
    # apply previous method to each dimension
    allSegslist <- vector("list", p)
    for (l in 1:p) {
      allSegslist[[l]] <- test_all_scenarii(as.vector(D[[l]]), Xp[[l]],
                                            scenarii, root)
    }
    # Find lowest total cost
    totalCostsAll <- vector("double", dim(scenarii)[1])
    for (i in 1:length(totalCostsAll)){
      tmp <- lapply(allSegslist, function(z) return(z[[i]]))
      totalCostsAll[i] <- sum(sapply(tmp, function(z) return(z$totalCost)))
    }
    min_conf <- which.min(totalCostsAll)
    # Go back to good format
    ret <- vector("list", p)
    for (l in 1:p){
      dd <- do.call(rbind, allSegslist[[l]])
      res <- dd[min_conf,]
      delta <- rep(0, dim(Xp[[l]])[2])
      delta[res$shifts_edges] <- res$coef[-1]
      shifts <- shifts.vector_to_list(delta)
      ret[[l]] <- list(beta_0 = res$coef[1], shifts = shifts, costs = res$costs)
    }
    return(ret)
  }
}

# best_scenario <- function (D, Xp, scenarii, root) {
#   # Function to be applyed to each row
#   fun <- function(sh_ed){
#     fit.lm <- lm.fit(x = Xp[, c(root, sh_ed), drop = FALSE], y = D)
#     squared_res <- residuals(fit.lm)^2
#     return(list(shifts_edges = sh_ed,
#                 coefs = unname(coef(fit.lm)),
#                 costs = squared_res,
#                 totalCost = sum(squared_res)))
#   }
#   # Apply to each row and take the minimal total cost
#   allSegs <- apply(scenarii, 1, fun)
#   dd <- do.call(rbind, allSegs)
#   min_conf <- which.min(dd[,"totalCost"])
#   # Go back to good format
#   res <- dd[min_conf,]
#   delta <- rep(0, dim(Xp)[2])
#   delta[res$shifts_edges] <- res$coef[-1]
#   shifts <- shifts.vector_to_list(delta)
#   return(list(beta_0 = res$coef[1], shifts = shifts, costs = res$costs))
# }
