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
#' @param K the dimension of the model
#' @param model_complexity the complexity of the set of models with dimention K
#' @param B a non-negative constant. Default is 0.1 
#' (as suggested in Cleymen Lebarbier 2015)
#' 
#' @return value of the penalty
#'
##
penalty_BirgeMassart_shape1 <- function(K, model_complexity, B = 0.1){
  if (B <= 0) stop("Constant B in penalty shape 1 must be non-negative.")
  return((sqrt(K) + sqrt(2 * B * K + 2 * log(model_complexity)))^2)
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
#' @param K the dimension of the model.
#' @param model_complexity the complexity of the set of models with dimention K.
#' @param C a non-negative constant. Default is 2.5 
#' (as suggested in Lebarbier 2005)
#' 
#' @return value of the penalty.
#'
##

penalty_BirgeMassart_shape2 <- function(K, model_complexity, C = 2.5){
  if (C <= 0) stop("Constant C in penalty shape 2 must be non-negative.")
  return(C*K + log(model_complexity))
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
  Delta <- log(model_complexity) + log(K + 1)
  res <- LINselect::penalty(Delta, n = ntaxa, p = 2 * ntaxa - 2, K = C)
  return(ntaxa * log(1 + res/(ntaxa - K)))
}