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
#' See Birge Massart (2001).
#' Must be applied to least-square crierion.
#' This penalty should be calibrated using the slope heuristic.
#'
#' @param K the number of shifts
#' @param p the dimention of the data
#' @param model_complexity the complexity of the set of models with dimention K
#' @param B a non-negative constant. Default is 0.1 
#' (as suggested in Cleymen Lebarbier 2015)
#' 
#' @return value of the penalty
#'
##
penalty_BirgeMassart_shape1 <- function(K, p, model_complexity, B = 0.1){
  if (B <= 0) stop("Constant B in penalty shape 1 must be non-negative.")
  #return((sqrt(K) + sqrt(2 * B * K + 2 * log(model_complexity)))^2)
  return((K + 1) * p * (1 + sqrt(2) * sqrt(B + 1/((K + 1) * p) * log(model_complexity)))^2)
}

model_selection_BM1 <- function(res, C.BM1, ...){
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
  cap_res <- capushe(data_capushe)
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

assign_selected_model_capushe <- function(res, cap_res){
  res$results_summary$K_select <- as.numeric(cap_res@DDSE@model)
  if (cap_res@DDSE@model != cap_res@Djump@model){
    res$results_summary$K_select_DDSE <- as.numeric(cap_res@DDSE@model)
    res$results_summary$K_select_Djump <- as.numeric(cap_res@Djump@model) 
  }
  res$results_summary$pen_shape
  return(res)
}

##
#' @title Penalty function type Birgé-Massart 2
#'
#' @description
#' \code{penalty_BirgeMassart_shape2} is the penalty shape defined by : 
#' pen_shape = C*K_try + log(model_complexity).
#' It dominates the penalty defined by \code{penalty_BirgeMassart_shape1}.
#'
#' @details
#' See Birge Massart (2001).
#' Must be applied to least-square crierion.
#' This penalty should be calibrated using the slope heuristic.
#'
#' @param K the number of shifts
#' @param p the dimension of the data
#' @param model_complexity the complexity of the set of models with dimention K.
#' @param C a non-negative constant. Default is 2.5 
#' (as suggested in Lebarbier 2005)
#' 
#' @return value of the penalty.
#'
##

penalty_BirgeMassart_shape2 <- function(K, p, model_complexity, C = 2.5){
  if (C <= 0) stop("Constant C in penalty shape 2 must be non-negative.")
  #return(C*K + log(model_complexity))
  return(C * (K + 1) + log(model_complexity))
}

model_selection_BM2 <- function(res, C.BM2, ...){
  p <- nrow(res$params_estim$`0`$variance)
  pen_shape <- penalty_BirgeMassart_shape1(res$results_summary$K_try,
                                           p,
                                           res$results_summary$complexity,
                                           C.BM2)
  res <- model_selection_capushe(res, pen_shape, "BM2")
  return(res)
}

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
#' Must be applied to log-likelihood crierion.
#' Function pen is computed using function \code{penalty} from package
#' \code{LINselect}.
#'
#' @param K the dimension of the model.
#' @param model_complexity the complexity of the set of models with dimention K.
#' @param ntaxa the number of tips.
#' @param C a constant, C > 1. Default is C = 1.1
#' (as suggested in Baraud Giraud Huet (2009))
#' 
#' @return value of the penalty.
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

model_selection_BGH <- function(res, ntaxa, C.BGH, ...){
  p <- nrow(res$params_estim$`0`$variance)
  ## Penalty
  pen <- 1/2 * penalty_BaraudGiraudHuet_likelihood(res$results_summary$K_try,
                                                   res$results_summary$complexity,
                                                   ntaxa,
                                                   C.BGH)
  ## Criterion
  crit <- - res$results_summary$log_likelihood + pen
  ## Assign results
  res <- assign_results_model_selection(res, pen, crit, "BGH")
  return(res)
}

assign_results_model_selection <- function(res, pen, crit, name){
  ## Fill result summary
  res$results_summary[[paste0("pen_", name)]] <- pen
  res$results_summary[[paste0("crit_", name)]] <- crit
  K_select <- res$results_summary$K_try[which.min(crit)]
  res$results_summary[[paste0("K_select_", name)]] <- K_select
  res$K_select[[paste0("K_select_", name)]] <- K_select
  ## Extract parameters and reconstructions
  res[[paste0(name)]]$params_select <- res$params_estim[[paste(K_select)]]
  res[[paste0(name)]]$Yhat <- res$Yhat[[paste(K_select)]]
  res[[paste0(name)]]$Zhat <- res$Zhat[[paste(K_select)]]
  res[[paste0(name)]]$Yvar <- res$Yvar[[paste(K_select)]]
  res[[paste0(name)]]$Zvar <- res$Zvar[[paste(K_select)]]
  if (attr(res[[paste0(name)]]$params_select, "Neq") > 1) {
    message(paste0("There are some equivalent solutions to the set of shifts selected by the ",
                   name, " method."))
    }
  return(res)
}