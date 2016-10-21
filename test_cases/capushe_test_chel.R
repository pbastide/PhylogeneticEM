rm(list=ls())

library(ape)
library(glmnet) # For Lasso initialization
library(quadrupen)
#library(nlme) # For second derivative computation
#library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(ggplot2) # Plot
library(reshape2) # Plot
library(grid) # Plot
library(TreeSim)
library(plyr)
source("R/simulate.R")
source("R/estimateEM.R")
source("R/init_EM.R")
source("R/E_step.R")
source("R/M_step.R")
source("R/shutoff.R")
source("R/generic_functions.R")
source("R/shifts_manipulations.R")
source("R/plot_functions.R")
source("R/parsimonyNumber.R")

library(capushe)
# ## Easy case
# ntaxa <- 64 #256 #64
# K_true <- 5
# Ks <- 0:35
# data_type <- paste0("easy_ntaxa=", ntaxa, "K_true=", K_true)
# load(paste0("../Results/Likelihood_Plot/", data_type, "_estimation_K_in_", paste(Ks, collapse = "_"), ".RData"))

## Simulation
# ntaxa <- 64
# K_true <- 10
# Ks <- 1:50
# data_type <- paste0("random_ntaxa=", ntaxa, "K_true=", K_true)
# load(paste0("../Results/Likelihood_Plot/", data_type, "_estimation_K_in_", paste(Ks, collapse = "_"), ".RData"))

## Chelonia
# Ks <- 1:63
# data_type <- "chelonia"
# PATH <- "../Results/Chelonia/"
# load(paste0("../Results/Likelihood_Plot/", data_type, "_estimation_K_in_", paste(Ks, collapse = "_"), ".RData"))

## Chelonia bis
Ks <- 0:30
data_type <- "chelonia"
PATH <- "../Results/Chelonia/"
# alpha <- 0.15 # Pick a "reasonable" alpha
# data_type <- paste0(data_type, "_alpha=", alpha)
load(paste0("../Results/Likelihood_Plot/", data_type, "_estimation_K_in_", paste(Ks, collapse = "_"), ".RData"))

## Contrast function
# times_shared <- compute_times_ca(tree)
# distances_phylo <- compute_dist_phy(tree)
# T_tree <- incidence.matrix(tree)
# 
# least_squares.stationary_root <- function(tree, times_shared, distances_phylo, Y_data, alpha, shifts, beta_0){
#   ## Computation of m_Y
#   ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
#                                                     selection.strength = alpha,
#                                                     times_shared = times_shared)
#   T_tree_ac <- T_tree * ac_tree
#   delta <- shifts.list_to_vector(tree, shifts)
#   m_Y <- T_tree_ac %*% delta + beta_0
#   ## Matrix D
#   D <- exp(-alpha * distances_phylo[1:ntaxa,1:ntaxa])
#   Dinv <- solve(D)
#   return(t(Y_data - m_Y)%*%Dinv%*%(Y_data-m_Y))
# }
# 
# df$least_squares <- vector(length = length(Ks))
# for (K in Ks){
#   df$least_squares[K] <- least_squares.stationary_root(tree, times_shared, distances_phylo, data, df$alpha_estim[K], dd[,"shifts_estim"][[K]], df$beta_0_estim[K])
# }

#################################################################################
## Type Birgé Massart

## Penalties
pen1 <- function(K_try, model_complexity, c){
  return(c*K_try + log(model_complexity))
}

pen2 <- function(K_try, model_complexity, B){
  return((sqrt(K_try) + sqrt(2 * B * K_try + 2 * log(model_complexity)))^2)
}

B <- 0.1
plot(df$K_try, -df$least_squares)
plot(pen2(df$K_try, df$model_complexity, B), -df$least_squares)
plot(df$K_try, df$least_squares + pen2(df$K_try, df$model_complexity, B))
plot(df$K_try, pen2(df$K_try, df$model_complexity, B))
K_max <- 35
df_sub <- subset(df, K_try <= K_max)
pen_shape <- pen2(df_sub$K_try, df_sub$model_complexity, B)
data_capushe <- data.frame(names = df_sub$K_try, 
                           pen_shape = pen_shape,
                           complexity = df_sub$model_complexity,
                           contrast = df_sub$least_squares)
# Slope Heuristic
DDSE_results <- DDSE(data_capushe)
DDSE_results@model
DDSE_results@interval$interval
plot(DDSE_results)
# Dimension Jump
Djump_results <- Djump(data_capushe)
Djump_results
Djump_results@ModelHat$Kopt
plot(Djump_results)
# Comparison
plot(df_sub$K_try, df_sub$least_squares + 2*DDSE_results@interval$interval["max"] * pen_shape)
plot(df_sub$K_try, df_sub$least_squares + Djump_results@ModelHat$Kopt * pen_shape)
plot(df_sub$K_try, df_sub$least_squares + 0.2 *pen_shape)

## Plots
# True process
plot.process.actual(Y.state = extract_simulate_internal(XX, what="states", where="tips"),
                    Z.state = extract_simulate_internal(XX, what="states", where="nodes"),
                    phylo = tree, 
                    paramsEstimate = list(shifts = shifts, optimal.value = beta_0))
# DDSE
Kt <- as.integer(DDSE_results@model) + 1
plot.process.actual(Y.state = extract_simulate_internal(XX, what="states", where="tips"),
                    Z.state = estimations[[Kt]]$Zhat,
                    phylo = tree, 
                    paramsEstimate = list(shifts = estimations[[Kt]]$shifts, 
                                          optimal.value = estimations[[Kt]]$beta_0_estim))
# Djump
Kt <- as.integer(Djump_results@model) + 1
plot.process.actual(Y.state = extract_simulate_internal(XX, what="states", where="tips"),
                    Z.state = estimations[[Kt]]$Zhat,
                    phylo = tree, 
                    paramsEstimate = list(shifts = estimations[[Kt]]$shifts, 
                                          optimal.value = estimations[[Kt]]$beta_0_estim))
# True K
Kt <- 5 + 1
plot.process.actual(Y.state = extract_simulate_internal(XX, what="states", where="tips"),
                    Z.state = estimations[[Kt]]$Zhat,
                    phylo = tree, 
                    paramsEstimate = list(shifts = estimations[[Kt]]$shifts, 
                                          optimal.value = estimations[[Kt]]$beta_0_estim))

## Comparisons of selected model with stopping point.
Bs <- 0.1
Kt <- 11:35
DDSE_results <- matrix(NA, nrow = length(Bs), ncol = length(Kt))
Djump_results <- matrix(NA, nrow = length(Bs), ncol = length(Kt))
for (B in Bs){
  for (K_max in Kt){
    df_sub <- subset(df, K_try <= K_max)
    data_capushe <- data.frame(names = df_sub$K, 
                               pen_shape = pen2(df_sub$K_try, df_sub$model_complexity, B),
                               complexity = df_sub$model_complexity,
                               contrast = df_sub$least_squares)
    # Slope Heuristic
    DDSE_results[which(Bs == B), which(Kt == K_max)] <- DDSE(data_capushe)@model
    # Dimension Jump
    Djump_results[which(Bs == B), which(Kt == K_max)] <- Djump(data_capushe)@model
  }
}
colnames(DDSE_results) <- Kt
rownames(DDSE_results) <- Bs
colnames(Djump_results) <- Kt
rownames(Djump_results) <- Bs

#####################################################################
## Type Baraud et al

## penalty
pens1 <- function(K_try, model_complexity, C){
  Delta <- log(model_complexity) + log(K_try + 1)
  res <- LINselect::penalty(Delta, n = ntaxa, p = 2*ntaxa-2, K = C)
  return(ntaxa * log(1 + res/(ntaxa - K_try)))
}

pens2 <- function(K_try, model_complexity, C = 1){
  Delta <- log(model_complexity) + log(K_try + 1)
  res <- LINselect::penalty(Delta, n = ntaxa, p = 2*ntaxa-2, K = C)
  return(ntaxa * res/(ntaxa - K_try))
}

C <- 1.1
pen_shape1 <- 1/2 * pens1(df$K_try, df$model_complexity, C)
pen_shape2 <- 1/2 * pens2(df$K_try, df$model_complexity)

## Maximum number of shifts
kappa <- 0.9
Dmax <- min(floor(kappa * ntaxa / (2 + log(2) + log(ntaxa))), ntaxa - 7)

## Basic plots
plot(df$K_try, df$log_likelihood)

plot(pen_shape1, df$log_likelihood)
plot(df$K_try, -df$log_likelihood + pen_shape1)
plot(df$K_try, pen_shape1)
names(which.min(-df$log_likelihood + pen_shape1))

plot(pen_shape2, df$log_likelihood)
plot(df$K_try, -df$log_likelihood + pen_shape2)
plot(df$K_try, pen_shape2)

## Capushe for pen 2
K_max <- Dmax
df_sub <- subset(df, K_try <= K_max)
pen_shape_sub <- 1/2 * pens2(df_sub$K_try, df_sub$model_complexity)
data_capushe <- data.frame(names = df_sub$K_try, 
                           pen_shape = pen_shape_sub,
                           complexity = df_sub$model_complexity,
                           contrast = -df_sub$log_likelihood)
# Slope Heuristic
DDSE_results <- DDSE(data_capushe)
DDSE_results@model
DDSE_results@interval$interval
plot(DDSE_results)
# Dimension Jump
Djump_results <- Djump(data_capushe)
Djump_results
Djump_results@ModelHat$Kopt
plot(Djump_results)
# Comparison
plot(df_sub$K_try, -df_sub$log_likelihood + 2*DDSE_results@interval$interval["max"] * pen_shape_sub)
plot(df_sub$K_try, -df_sub$log_likelihood + Djump_results@ModelHat$Kopt * pen_shape_sub)
plot(df_sub$K_try, -df_sub$log_likelihood + pen_shape_sub)

