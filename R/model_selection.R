# {model selection}
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
## Model Selection related functions.
###############################################################################

##
#' @title Penalty function type Birgé-Massart 1
#'
#' @description
#' \code{penalty_BirgeMassart_shape1} is the penalty shape defined by : 
#' pen_shape = (sqrt(K) + sqrt(2 * K * L_K))^2 with sum(exp(- K * L_K)) < infty :
#' L_K = B + 1/K * log(model_complexity).
#'
#' @details
#' See Birgé Massart (2001).
#' Must be applied to least-square criterion.
#' This penalty should be calibrated using the slope heuristic.
#'
#' @param K the number of shifts
#' @param p the dimension of the data
#' @param model_complexity the complexity of the set of models with dimension K
#' @param B a non-negative constant. Default is 0.1 
#' (as suggested in Cleynen & Lebarbier 2015)
#' 
#' @return value of the penalty
#' 
#' @seealso \code{\link{penalty_BaraudGiraudHuet_likelihood}},
#' \code{\link{penalty_BirgeMassart_shape2}}
#' 
#' @keywords internal
#'
##
penalty_BirgeMassart_shape1 <- function(K, p, model_complexity, B = 0.1){
  if (B <= 0) stop("Constant B in penalty shape 1 must be non-negative.")
  #return((sqrt(K) + sqrt(2 * B * K + 2 * log(model_complexity)))^2)
  return((K + 1) * p * (1 + sqrt(2) * sqrt(B + 1/((K + 1) * p) * log(model_complexity)))^2)
}

model_selection_BM1 <- function(res, C.BM1, ...){
  res <- res$alpha_max
  p <- nrow(res$params_estim$`0`$variance)
  pen_shape <- penalty_BirgeMassart_shape1(res$results_summary$K_try,
                                           p,
                                           res$results_summary$complexity,
                                           C.BM1)
  res <- model_selection_capushe(res, pen_shape, "BM1")
  return(res)
}

model_selection_capushe <- function(res, pen_shape, name){
  ## Format data for capushe
  data_capushe <- data.frame(names = res$results_summary$K_try, 
                             pen_shape = pen_shape,
                             complexity = res$results_summary$complexity,
                             contrast = -res$results_summary$log_likelihood)
  ## Capushe
  cap_res <- capushe::capushe(data_capushe)
  ## Assign results
  res[[paste0("capushe_output", name)]] <- cap_res
  res$results_summary[[paste0("pen_shape", name)]]  <- pen_shape
  # DDSE
  pen_DDSE <- 2*cap_res@DDSE@interval$interval["max"]*pen_shape
  crit_DDSE <- data_capushe$contrast + pen_DDSE
  res <- assign_results_model_selection(res, pen_DDSE, crit_DDSE, paste0("DDSE_", name))
  # Djump
  pen_Djump <- cap_res@Djump@ModelHat$Kopt*pen_shape
  crit_Djump <- data_capushe$contrast + pen_Djump
  res <- assign_results_model_selection(res, pen_Djump, crit_Djump, paste0("Djump_", name))
  return(res)
}

# assign_selected_model_capushe <- function(res, cap_res){
#   res$results_summary$K_select <- as.numeric(cap_res@DDSE@model)
#   if (cap_res@DDSE@model != cap_res@Djump@model){
#     res$results_summary$K_select_DDSE <- as.numeric(cap_res@DDSE@model)
#     res$results_summary$K_select_Djump <- as.numeric(cap_res@Djump@model) 
#   }
#   res$results_summary$pen_shape
#   return(res)
# }

##
#' @title Penalty function type Birgé-Massart 2
#'
#' @description
#' \code{penalty_BirgeMassart_shape2} is the penalty shape defined by : 
#' pen_shape = C*K_try + log(model_complexity).
#' It dominates the penalty defined by \code{penalty_BirgeMassart_shape1}.
#'
#' @details
#' See Birgé Massart (2001).
#' Must be applied to least-square criterion.
#' This penalty should be calibrated using the slope heuristic.
#'
#' @param K the number of shifts
#' @param p the dimension of the data
#' @param model_complexity the complexity of the set of models with dimension K.
#' @param C a non-negative constant. Default is 2.5 
#' (as suggested in Lebarbier 2005)
#' 
#' @return value of the penalty.
#' 
#' @seealso \code{\link{penalty_BirgeMassart_shape1}},
#' \code{\link{penalty_BaraudGiraudHuet_likelihood}}
#' 
#' @keywords internal
#' 
##

penalty_BirgeMassart_shape2 <- function(K, p, model_complexity, C = 2.5){
  if (C <= 0) stop("Constant C in penalty shape 2 must be non-negative.")
  #return(C*K + log(model_complexity))
  return(C * (K + 1) + log(model_complexity))
}

model_selection_BM2 <- function(res, C.BM2, ...){
  res <- res$alpha_max
  p <- nrow(res$params_estim$`0`$variance)
  pen_shape <- penalty_BirgeMassart_shape2(res$results_summary$K_try,
                                           p,
                                           res$results_summary$complexity,
                                           C.BM2)
  res <- model_selection_capushe(res, pen_shape, "BM2")
  return(res)
}

#############################################
## LINselect
#############################################
##
#' @title Penalty function type Baraud Giraud Huet.
#'
#' @description
#' \code{penalty_BaraudGiraudHuet_likelihood} is the penalty defined by : 
#' pen' = ntaxa * log(1 + pen/(ntaxa - K)) with 
#' pen = C * (ntaxa - K)/(ntaxa - K - 1) * EDkhi[K + 1; ntaxa - K - 1; exp(-Delta_K)/(K + 1)]
#' and Delta = log(model_complexity) + log(K + 1) 
#' such that sum(exp(-Delta_K)) < infty.
#'
#' @details
#' See Baraud Giraud Huet (2009, 2011).
#' Must be applied to log-likelihood criterion.
#' Function pen is computed using function \code{penalty} from package
#' \code{LINselect}.
#'
#' @param K the dimension of the model.
#' @param model_complexity the complexity of the set of models with dimension K.
#' @param ntaxa the number of tips.
#' @param C a constant, C > 1. Default is C = 1.1
#' (as suggested in Baraud Giraud Huet (2009))
#' 
#' @return value of the penalty.
#' 
#' @seealso \code{\link{penalty_BirgeMassart_shape1}},
#' \code{\link{penalty_BirgeMassart_shape2}}
#' 
#' @keywords internal
#'
##

penalty_BaraudGiraudHuet_likelihood <- function(K, model_complexity, ntaxa, 
                                                C = 1.1){
  Delta <- log(model_complexity) + log(K + 2)
  Delta <- c(0, Delta) # We start with dimension 1 (K=0)
  res <- LINselect::penalty(Delta, n = ntaxa, p = 2 * ntaxa - 2, K = C)
  res <- res[-1]
  return(ntaxa * log(1 + res/(ntaxa - K - 1)))
}

model_selection_BGH <- function(res, ntaxa, C.LINselect, ...){
  res <- res$alpha_max
  p <- nrow(res$params_estim$`0`$variance)
  ## Penalty
  pen <- 1/2 * penalty_BaraudGiraudHuet_likelihood(res$results_summary$K_try,
                                                   res$results_summary$complexity,
                                                   ntaxa,
                                                   C.LINselect)
  ## Criterion
  crit <- - res$results_summary$log_likelihood + pen
  ## Assign results
  res <- assign_results_model_selection(res, pen, crit, "BGHuni")
  return(res)
}

penalty_BaraudGiraudHuet_leastsquares <- function(K, model_complexity, ntaxa, 
                                                  C = 1.1){
  Delta <- log(model_complexity) + log(K + 2)
  Delta <- c(0, Delta) # We start with dimension 1 (K=0)
  res <- LINselect::penalty(Delta, n = ntaxa, p = 2 * ntaxa - 2, K = C)
  res <- res[-1]
  return(res / (ntaxa - K - 1))
}

model_selection_BGH_leastsquares <- function(res, ntaxa, C.LINselect, ...){
  # res <- add_lsq(res)
  # res <- merge_min_grid_alpha(res)
  res <- res$alpha_min
  p <- nrow(res$params_estim$`0`$variance)
  ## Penalty
  pen <- penalty_BaraudGiraudHuet_leastsquares(res$results_summary$K_try,
                                               res$results_summary$complexity,
                                               ntaxa,
                                               C.LINselect)
  ## least squares
  # lsq <- sapply(res$params_estim, function(z) sum(diag(z$variance)))
  lsq <- res$results_summary$least_squares
  crit <- ntaxa * lsq * (1 + pen)
  ## Assign results
  res <- assign_results_model_selection(res, pen, crit, "BGHlsq")
  return(res)
}

model_selection_BGH_ml <- function(res, ntaxa, C.LINselect, ...){
  # res <- add_lsq(res)
  # res <- merge_min_grid_alpha(res)
  res <- res$alpha_max
  p <- nrow(res$params_estim$`0`$variance)
  ## Penalty
  pen <- penalty_BaraudGiraudHuet_leastsquares(res$results_summary$K_try,
                                               res$results_summary$complexity,
                                               ntaxa,
                                               C.LINselect)
  ## least squares
  # lsq <- sapply(res$params_estim, function(z) sum(diag(z$variance)))
  lsq <- res$results_summary$least_squares
  crit <- ntaxa * lsq * (1 + pen)
  ## Assign results
  res <- assign_results_model_selection(res, pen, crit, "BGHml")
  return(res)
}

model_selection_BGH_mlraw <- function(res, ntaxa, C.LINselect, ...){
  # res <- add_lsq(res)
  # res <- merge_min_grid_alpha(res)
  res <- res$alpha_max
  p <- nrow(res$params_estim$`0`$variance)
  ## Penalty
  pen <- penalty_BaraudGiraudHuet_leastsquares(res$results_summary$K_try,
                                               res$results_summary$complexity,
                                               ntaxa,
                                               C.LINselect)
  ## least squares
  # lsq <- sapply(res$params_estim, function(z) sum(diag(z$variance)))
  lsq <- res$results_summary$least_squares_raw
  crit <- ntaxa * lsq * (1 + pen)
  ## Assign results
  res <- assign_results_model_selection(res, pen, crit, "BGHmlraw")
  return(res)
}

model_selection_BGH_leastsquares_raw <- function(res, ntaxa, C.LINselect, ...){
  # res <- add_lsq(res)
  # res <- merge_min_grid_alpha(res)
  res <- res$alpha_min_raw
  p <- nrow(res$params_estim$`0`$variance)
  ## Penalty
  pen <- penalty_BaraudGiraudHuet_leastsquares(res$results_summary$K_try,
                                               res$results_summary$complexity,
                                               ntaxa,
                                               C.LINselect)
  ## least squares
  # lsq <- sapply(res$params_estim, function(z) sum(diag(z$variance)))
  lsq <- res$results_summary$least_squares_raw
  crit <- ntaxa * lsq * (1 + pen)
  ## Assign results
  res <- assign_results_model_selection(res, pen, crit, "BGHlsqraw")
  return(res)
}

#############################################
## pBIC
#############################################
##
#' @title Penalty function type pBIC
#'
#' @description
#' \code{penalty_pBIC_scalarOU} is the pBIC.
#'
#' @param K the dimension of the model.
#' @param model_complexity the complexity of the set of models with dimension K.
#' @param ntaxa the number of tips.
#' 
#' @return value of the penalty.
#' 
#' @seealso \code{\link{penalty_BirgeMassart_shape1}},
#' \code{\link{penalty_BirgeMassart_shape2}}
#' 
#' @keywords internal
#'
##

penalty_pBIC <- function(all_params, model_complexity, independent, tree, times_shared,
                         distances_phylo, T_tree, p, K, ntaxa, process){ # nocov start
  if (independent){
    penalty_pBIC_unit <- penalty_pBIC_independent
  } else {
    penalty_pBIC_unit <- penalty_pBIC_scalarOU
  }
  ## Computations
  pen <- sapply(all_params, penalty_pBIC_unit,
                tree, times_shared, distances_phylo, T_tree, process)
  ## Model Complexity term
  pen <- pen + 2 * log(model_complexity)
  ## Constant term
  Cst <- p * (p+1)/2 * log(ntaxa - K - 1)
  # Cst <- Cst - p * (K + 1) * log(2 * pi)
  # Cst <- Cst - p*(p+1)/2 * log(2*pi) + 2*p * log(2)
  pen <- pen + Cst
  ## Alpha parameter
  if (independent){
    pen <- pen + p * log(ntaxa)
  } else {
    pen <- pen + log(ntaxa) 
  }
  return(pen)
} # nocov end

penalty_pBIC_scalarOU <- function(params, tree, times_shared,
                                  distances_phylo, T_tree, process){ # nocov start
  if (length(as.vector(params$selection.strength)) > 1){
    stop("pBIC only works for scalar OU.")
  }
  ntaxa <- length(tree$tip.label)
  K <- length(params$shifts$edges)
  p <- ncol(params$variance)
  ## Variance term
  # pen <- (p - K) * determinant(params$variance, logarithm = TRUE)$modulus
  # pen <- as.vector(pen)
  ## Model Design term
  # Correlation matrix
  compute_tree_correlations_matrix  <- switch(process, 
                                              BM = compute_tree_correlations_matrix.BM,
                                              OU = compute_tree_correlations_matrix.scOU,
                                              scOU = compute_tree_correlations_matrix.scOU)
  C <- compute_tree_correlations_matrix(times_shared, distances_phylo, params)
  C <- extract_variance_covariance(C, what="YY",
                                   masque_data = c(rep(TRUE, ntaxa),
                                                   rep(FALSE, dim(C)[1] - ntaxa)))
  C <- 1/(2*as.vector(params$selection.strength)) * C
  C_inv <- solve(C)
  # Tree Matrix
  ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                                                    selection.strength = as.vector(params$selection.strength),
                                                    times_shared = times_shared)
  Tr <- T_tree * ac_tree
  Tr <- Tr[, params$shifts$edges, drop = F]
  Tr <- cbind(Tr, rep(1, dim(Tr)[1]))
  # Det
  pen <- p * determinant(t(Tr) %*% C_inv %*% Tr, logarithm = TRUE)$modulus
  pen <- as.vector(pen)
  return(pen)
} # nocov end

penalty_pBIC_independent <- function(params, model_complexity, tree,
                                     times_shared, distances_phylo,
                                     T_tree, process){ # nocov start
  ntaxa <- length(tree$tip.label)
  K <- length(params$shifts$edges)
  p <- ncol(params$variance)
  ## Variance term
  pen <- (p - K) * determinant(params$variance, logarithm = TRUE)$modulus
  pen <- as.vector(pen)
  ## Model Design term
  # Correlation matrix
  compute_tree_correlations_matrix  <- switch(process, 
                                              BM = compute_tree_correlations_matrix.BM,
                                              OU = compute_tree_correlations_matrix.scOU,
                                              scOU = compute_tree_correlations_matrix.scOU)
  C <- compute_tree_correlations_matrix(times_shared, distances_phylo, params)
  C <- extract_variance_covariance(C, what="YY",
                                   masque_data = c(rep(TRUE, ntaxa),
                                                   rep(FALSE, dim(C)[1] - ntaxa)))
  C_inv <- solve(C)
  alphas <- diag(params$selection.strength)
  for (alp in alphas){
    # Tree Matrix
    ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                                                      selection.strength = alp,
                                                      times_shared = times_shared)
    Tr <- T_tree * ac_tree
    Tr <- Tr[, params$shifts$edges, drop = F]
    Tr <- cbind(Tr, rep(1, dim(Tr)[1]))
    # Det
    pen <- pen + determinant(t(Tr) %*% C_inv %*% Tr, logarithm = TRUE)
    pen <- as.vector(pen)
  }
  return(pen)
} # nocov end

model_selection_pBIC <- function(res, independent, tree, times_shared, 
                                 distances_phylo, T_tree, 
                                 ntaxa, process, ...){ # nocov start
  res <- res$alpha_max
  p <- nrow(res$params_estim$`0`$variance)
  ## Penalty
  pen <- 1/2 * penalty_pBIC(res$params_estim,
                            res$results_summary$complexity,
                            independent, tree,
                            times_shared, distances_phylo, T_tree,
                            p, res$results_summary$K_try, ntaxa, process)
  ## Criterion
  crit <- - res$results_summary$log_likelihood + pen
  ## Assign results
  res <- assign_results_model_selection(res, pen, crit, "pBIC")
  return(res)
} # nocov end

#############################################
## l1ou
#############################################
penalty_pBIC_l1ou <- function(all_params, model_complexity,
                              independent, tree, times_shared,
                              distances_phylo, T_tree, p, K, ntaxa,
                              process, Y_data){ # nocov start
  ## Computations
  pen <- sapply(all_params, penalty_pBIC_l1ou_unit,
                tree, times_shared, distances_phylo, T_tree, process, Y_data)
  ## Model Complexity term
  pen <- pen + 2 * log(model_complexity)
  return(pen)
} # nocov end

penalty_pBIC_l1ou_unit <- function(params, tree, times_shared, distances_phylo,
                                   T_tree, process, Y_data){ # nocov start
  ntaxa <- length(tree$tip.label)
  K <- length(params$shifts$edges)
  p <- ncol(params$variance)
  ## Variance term
  # browser()
  vars <- apply(Y_data, 1, var)
  vars <- vars - diag(params$variance)
  pen <- sum(log(abs(vars)))
  pen <- (K+1) * as.vector(pen)
  ## Model Design term
  # Correlation matrix
  compute_tree_correlations_matrix  <- switch(process, 
                                              BM = compute_tree_correlations_matrix.BM,
                                              OU = compute_tree_correlations_matrix.scOU,
                                              scOU = compute_tree_correlations_matrix.scOU)
  C <- compute_tree_correlations_matrix(times_shared, distances_phylo, params)
  C <- extract_variance_covariance(C, what="YY",
                                   masque_data = c(rep(TRUE, ntaxa),
                                                   rep(FALSE, dim(C)[1] - ntaxa)))
  C <- 1/(2*as.vector(params$selection.strength)) * C
  C_inv <- solve(C)
  # Tree Matrix
  ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                                                    selection.strength = as.vector(params$selection.strength),
                                                    times_shared = times_shared)
  Tr <- T_tree * ac_tree
  Tr <- Tr[, params$shifts$edges, drop = F]
  Tr <- cbind(Tr, rep(1, dim(Tr)[1]))
  # Det
  pen <- pen + p * determinant(t(Tr) %*% C_inv %*% Tr, logarithm = TRUE)$modulus
  pen <- as.vector(pen)
  # sigma and alpha
  pen <- pen + 2 * log(ntaxa)
  return(pen)
} # nocov end

model_selection_pBIC_l1ou <- function(res, independent, tree,
                                      times_shared, distances_phylo,
                                      T_tree, ntaxa, process, Y_data, ...){ # nocov start
  res <- res$alpha_max
  p <- nrow(res$params_estim$`0`$variance)
  ## Penalty
  pen <- 1/2 * penalty_pBIC_l1ou(res$params_estim,
                                 res$results_summary$complexity,
                                 independent, tree,
                                 times_shared, distances_phylo, T_tree,
                                 p, res$results_summary$K_try, ntaxa, process,
                                 Y_data)
  ## Criterion
  crit <- - res$results_summary$log_likelihood + pen
  ## Assign results
  res <- assign_results_model_selection(res, pen, crit, "pBIC_l1ou")
  return(res)
} # nocov end

#############################################
## Format results
#############################################

assign_results_model_selection <- function(res, pen, crit, name){
  ## Fill result summary
  res$results_summary[[paste0("pen_", name)]] <- pen
  res$results_summary[[paste0("crit_", name)]] <- crit
  K_select <- res$results_summary$K_try[which.min(crit)]
  res$results_summary[[paste0("K_select_", name)]] <- K_select
  res$K_select[[paste0("K_select_", name)]] <- K_select
  ## Extract parameters and reconstructions
  res[[paste0(name)]]$params_select <- res$params_estim[[paste(K_select)]]
  # res[[paste0(name)]]$params_raw <- res$params_raw[[paste(K_select)]]
  res[[paste0(name)]]$params_init_estim <- res$params_init_estim[[paste(K_select)]]
  res[[paste0(name)]]$results_summary <- res$results_summary[K_select + 1, ]
  # res[[paste0(name)]]$Yhat <- res$Yhat[[paste(K_select)]]
  # res[[paste0(name)]]$Zhat <- res$Zhat[[paste(K_select)]]
  # res[[paste0(name)]]$Yvar <- res$Yvar[[paste(K_select)]]
  # res[[paste0(name)]]$Zvar <- res$Zvar[[paste(K_select)]]
  # res[[paste0(name)]]$m_Y_estim <- res$m_Y_estim[[paste(K_select)]]
  res[[paste0(name)]]$edge.quality <- res$edge.quality[[paste(K_select)]]
  if (attr(res[[paste0(name)]]$params_select, "Neq") > 1) {
    message(paste0("There are some equivalent solutions to the set of shifts selected by the ",
                   name, " method."))
    }
  return(res)
}

#############################################
## Generic for model selection
#############################################
##
#' @title Model Selection of a fitted object
#'
#' @description
#' \code{model_selection} does the model selection on a fitted \code{\link{PhyloEM}} 
#' object, and returns the same fitted object.
#'
#' @param x a fitted \code{\link{PhyloEM}} object
#' @inheritParams PhyloEM
#'     
#' @return
#' The same object, but with a slot corresponding to the model selection used. See
#' function \code{\link{params_process.PhyloEM}} to retrieve the selected parameters.
#'  
#' @seealso \code{\link{PhyloEM}}, \code{\link{params_process.PhyloEM}},
#' \code{\link{imputed_traits.PhyloEM}}
#' 
#' @export
#' 
##
model_selection <- function(x, ...) UseMethod("model_selection")

##
#' @describeIn model_selection \code{\link{PhyloEM}} object
#' @export
##
model_selection.PhyloEM <- function(x,
                                    method.selection = c("LINselect", "DDSE", "Djump"),
                                    C.BM1 = 0.1, C.BM2 = 2.5, C.LINselect = 1.1,
                                    independent = FALSE, ...){
  mod_sel_unit <- function(one.method.selection){
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
    selection <- try(mod_sel(x, ntaxa = ncol(x$Y_data),
                             C.BM1 = C.BM1, C.BM2 = C.BM2, C.LINselect = C.LINselect,
                             tree = x$phylo, independent = independent,
                             T_tree = x$T_tree, times_shared = x$times_shared, 
                             distances_phylo = x$distances_phylo,
                             process = x$process, Y_data = x$Y_data))
    if (inherits(selection, "try-error")){
      warning(paste0("Model Selection ",  one.method.selection, " failled"))
    } else if (one.method.selection == "BGHlsq") {
      x$alpha_min <- selection
    } else if (one.method.selection == "BGHlsqraw") {
      x$alpha_min_raw <- selection
    } else {
      x$alpha_max <- selection
    }
    return(x)
  }
  method.selection  <- match.arg(method.selection,
                                 choices = c("LINselect", "DDSE", "Djump",
                                             "BirgeMassart1", "BirgeMassart2",
                                             "BGH", "BGHuni", "BGHlsq", "BGHml",
                                             "BGHlsqraw", "BGHmlraw",
                                             "pBIC", "pBIC_l1ou"),
                                 several.ok = TRUE)
  method.selection <- expand_method_selection(method.selection)
  if (x$p > 1){
    method.selection <- method.selection[method.selection != "BGH"]
    method.selection <- method.selection[method.selection != "BGHuni"]
    # warning("BGH is not implemented for multivariate data.")
  }
  if (x$p == 1){
    if ("BGH" %in% method.selection){
      method.selection[method.selection == "BGH"] <- "BGHuni"
    }
    method.selection <- method.selection[method.selection != "BGHlsq"]
    method.selection <- method.selection[method.selection != "BGHml"]
    method.selection <- method.selection[method.selection != "BGHlsqraw"]
    method.selection <- method.selection[method.selection != "BGHmlraw"]
  }
  if (length(method.selection) == 0) stop("No selection method were selected or suited to the problem (see relevant warnings).")
  for (meth.sel in method.selection){
    x <- mod_sel_unit(meth.sel)
  }
  return(x)
}