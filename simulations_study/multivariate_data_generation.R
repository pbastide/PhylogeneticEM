##############
## Parameters and functions
##############

require(doParallel)
require(foreach)
require(ape)
#require(glmnet) # For Lasso initialization
# require(quadrupen) # For Lasso initialization
# require(robustbase) # For robust fitting of alpha
library(TreeSim) # For simulation of the tree
library(Matrix)
#require(ggplot2) # Plot
#require(reshape2) # Plot
#require(grid) # Plot

## Source functions
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
source("R/partitionsNumber.R")

savedatafile = "../Results/Simulations_Multivariate/multivariate_simlist"

## Fixed parameters for simulations
process <- "OU"
p_base <- 4 # number of traits
beta_0 <- rep(0, p_base) # ancestral optimum
alpha_base <- 3 # selection strength
sigma_base <- 6 # sigma squared variance
gamma_base <- 1 # gamma squared stationnary variance
K_base <- 3 # number of shifts
ntaxa_base <- 160 # number of taxa
factor_shift_base <- 1 # multiplicative factor for shifts
r_base <- 0 # correlation coefficients
s_base <- 0 # non diagonal coefficient
NA_base <- 0 # % of NAs

## alpha grid
alpha_grid <- c(1, 2) # alpha varies with gamma squared fixed to 1

## Number of shifts for simulation
K_grid <- c(0, 3, 7, 11, 15) # number of shifts
factor_shift_grid <- c(0.5, 1, 1.5, 2) # multiplicative factor

## correlations
r_grid <- c(0.2, 0.4, 0.6) # non diagonal elements for A (r) or R (2r)

## Non scalar
s_grid <- c(0.5, 1, 1.5) # additive diagonal for two lowers of A

## NAs
NA_grid <- c(0.05, 0.1, 0.2, 0.5)

## Number of taxa
ntaxa_grid <- c(32, 64, 96, 128, 192, 215)

## replication depth (number of replicates per )
nrep <- 1:200

## The combination of simulation parameters
# Base
simparams_base <- expand.grid(alpha_base, gamma_base, K_base,
                              r_base, r_base, s_base,
                              factor_shift_base, ntaxa_base, NA_base, nrep,
                              "base")
# Alpha
simparams_alpha <- expand.grid(alpha_grid, gamma_base, K_base,
                               r_base, r_base, s_base,
                               factor_shift_base, ntaxa_base, NA_base, nrep,
                               "alpha_var")
# K AND Factors
simparams_K <- expand.grid(alpha_base, gamma_base, K_grid,
                           r_base, r_base, s_base,
                           factor_shift_grid, ntaxa_base, NA_base, nrep,
                           "K_var")
# r drift
simparams_r_drift <- expand.grid(alpha_base, gamma_base, K_base,
                                 r_grid, r_base, s_base,
                                 factor_shift_base, ntaxa_base, NA_base, nrep,
                                 "r_drift_var")
# r selection
simparams_r_selection <- expand.grid(alpha_base, gamma_base, K_base,
                                     r_base, r_grid, s_base,
                                     factor_shift_base, ntaxa_base, NA_base, nrep,
                                     "r_selection_var")
# s diagonal selection
simparams_s <- expand.grid(alpha_base, gamma_base, K_base,
                           r_base, r_base, s_grid,
                           factor_shift_base, ntaxa_base, NA_base, nrep,
                           "s_var")
# % of NAs
simparams_NA <- expand.grid(alpha_base, gamma_base, K_base,
                            r_base, r_base, s_base,
                            factor_shift_base, ntaxa_base, NA_grid, nrep,
                            "NA_var")
# ntaxa
simparams_ntaxa <- expand.grid(alpha_base, gamma_base, K_base,
                               r_base, r_base, s_base,
                               factor_shift_base, ntaxa_grid, NA_base, nrep,
                               "ntaxa_var")
simparams <- rbind(simparams_base,
                   simparams_alpha,
                   simparams_K,
                   simparams_r_drift,
                   simparams_r_selection,
                   simparams_s,
                   simparams_NA,
                   simparams_ntaxa
                  )
colnames(simparams) <- c("alpha", "gamma", "K", "rd", "rs", "s", "factor_shift",
                         "ntaxa", "NA_per", "nrep", "grp")
## Remove redondancies
# Base
mask <- (simparams$K == K_base) & (simparams$factor_shift == factor_shift_base) & (simparams$grp == "K_var")
simparams <- simparams[!mask, ]
# zero shift
mask <- (simparams$K == 0) & (simparams$factor_shift != 1) & (simparams$grp == "K_var")
simparams <- simparams[!mask, ] # 38*200 x 11

##############
## Generation of trees
##############
set.seed(17920904)
trees <- vector(mode = "list"); times_shared <- vector(mode = "list");
distances_phylo <- vector(mode = "list"); subtree.list <- vector(mode = "list");
T_tree <- vector(mode = "list"); K_try <- vector(mode = "list");
h_tree <- vector(mode = "list")
lambda <- 0.1
for (nta in c(128, 32, 64, 96, 160, 192, 215)){
  # Generate tree with nta taxa
  trees[[paste0(nta)]] <- sim.bd.taxa.age(n = nta, numbsim = 1, 
                                          lambda = lambda, mu = 0,
                                          age = 1, mrca = TRUE)[[1]]
  # Fixed tree quantities
  times_shared[[paste0(nta)]] <- compute_times_ca(trees[[paste0(nta)]])
  distances_phylo[[paste0(nta)]] <- compute_dist_phy(trees[[paste0(nta)]])
  subtree.list[[paste0(nta)]] <- enumerate_tips_under_edges(trees[[paste0(nta)]])
  T_tree[[paste0(nta)]] <- incidence.matrix(trees[[paste0(nta)]])
  h_tree[[paste0(nta)]] <- max(diag(times_shared[[paste0(nta)]])[1:nta])
  # Number of tries (depends on tree)
  K_try[[paste0(nta)]] <- 0:max(floor(sqrt(nta)), 10)
}

##############
## Generation of shifts
##############
plot(trees[["128"]], show.tip.label = FALSE); edgelabels(); tiplabels()

shifts_grid <- vector(mode = "list")

## 128 - 3
shifts_grid[["128_3"]] <- list(edges = c(9, 72, 209),
                               values=cbind(rep(2.4, p_base),
                                            rep(-2.1, p_base),
                                            rep(2.2, p_base)),
                               relativeTimes = 0)

# Means at the tips ?
Delta <- shifts.list_to_matrix(trees[["128"]], shifts_grid[["128_3"]])
W <- compute_actualization_matrix_ultrametric(trees[["128"]], alpha_base * diag(1, p_base, p_base))
vec_Y <- kronecker(T_tree[["128"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
X1.tips.exp.mat <- matrix(vec_Y, p_base, 128) + beta_0
unique(X1.tips.exp.mat[1, ])
plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
                         phylo = trees[["128"]], 
                         params = list(shifts = list(edges = shifts_grid[["128_3"]]$edges, values = shifts_grid[["128_3"]]$values[1, ])),
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
# Equivalent solutions ?
extract.parsimonyNumber(parsimonyNumber(trees[["128"]], clusters_from_shifts_ism(trees[["128"]], shifts_grid[["128_3"]]$edges)))

# ## 128 - 7
# shifts_grid[["128_7"]] <- list(edges = c(8, 72, 193,
#                                      117, 179, 33, 84),
#                            values = cbind(rep(2.1, p_base),
#                                           rep(-2.1, p_base),
#                                           rep(2.1, p_base),
#                                           rep(-2.1, p_base),
#                                           rep(-2.1, p_base),
#                                           rep(2.4, p_base),
#                                           rep(4.4, p_base)),
#                            relativeTimes = rep(0, p_base))
# 
# # Means at the tips ?
# Delta <- shifts.list_to_matrix(trees[["128"]], shifts_grid[["128_7"]])
# W <- compute_actualization_matrix_ultrametric(trees[["128"]], alpha_base * diag(1, p_base, p_base))
# vec_Y <- kronecker(T_tree[["128"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
# X1.tips.exp.mat <- matrix(vec_Y, p_base, ntaxa_base) + beta_0
# unique(X1.tips.exp.mat[1, ])
# plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
#                          phylo = trees[["128"]], 
#                          params = list(shifts = list(edges = shifts_grid[["128_7"]]$edges, values = shifts_grid[["128_7"]]$values[1, ])),
#                          adj.root = 0,
#                          automatic_colors = TRUE,
#                          margin_plot = NULL,
#                          cex = 2,
#                          bg_shifts = "azure2",
#                          bg_beta_0 = "azure2")
# # Equivalent solutions ?
# extract.parsimonyNumber(parsimonyNumber(trees[["128"]], clusters_from_shifts_ism(trees[["128"]], shifts_grid[["128_7"]]$edges)))
# 
# ## 128 - 11
# shifts_grid[["128_11"]] <- list(edges = c(8, 72, 193,
#                                       117, 179, 33, 84,
#                                       231, 210, 157, 127),
#                             values = cbind(rep(2.1, p_base),
#                                            rep(-2.1, p_base),
#                                            rep(2.1, p_base),
#                                            rep(-2.1, p_base),
#                                            rep(-2.1, p_base),
#                                            rep(2.4, p_base),
#                                            rep(4.4, p_base),
#                                            rep(-2.1, p_base),
#                                            rep(2.4, p_base),
#                                            rep(2.2, p_base),
#                                            rep(-2.5, p_base)),
#                             relativeTimes = rep(0, p_base))
# 
# # Means at the tips ?
# Delta <- shifts.list_to_matrix(trees[["128"]], shifts_grid[["128_11"]])
# W <- compute_actualization_matrix_ultrametric(trees[["128"]], alpha_base * diag(1, p_base, p_base))
# vec_Y <- kronecker(T_tree[["128"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
# X1.tips.exp.mat <- matrix(vec_Y, p_base, ntaxa_base) + beta_0
# unique(X1.tips.exp.mat[1, ])
# plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
#                          phylo = trees[["128"]], 
#                          params = list(shifts = list(edges = shifts_grid[["128_11"]]$edges, values = shifts_grid[["128_11"]]$values[1, ])),
#                          adj.root = 0,
#                          automatic_colors = TRUE,
#                          margin_plot = NULL,
#                          cex = 2,
#                          bg_shifts = "azure2",
#                          bg_beta_0 = "azure2")
# # Equivalent solutions ?
# extract.parsimonyNumber(parsimonyNumber(trees[["128"]], clusters_from_shifts_ism(trees[["128"]], shifts_grid[["128_11"]]$edges)))
# 
# ## 128 -  15
# shifts_grid[["128_15"]] <- list(edges = c(8, 72, 193,
#                                       117, 179, 33, 84,
#                                       231, 210, 157, 127,
#                                       20, 96, 199, 239),
#                             values = cbind(rep(2.1, p_base),
#                                            rep(-2.1, p_base),
#                                            rep(2.1, p_base),
#                                            rep(-2.1, p_base),
#                                            rep(-2.1, p_base),
#                                            rep(2.4, p_base),
#                                            rep(4.4, p_base),
#                                            rep(-2.1, p_base),
#                                            rep(2.4, p_base),
#                                            rep(2.2, p_base),
#                                            rep(-2.5, p_base),
#                                            rep(-5.2, p_base),
#                                            rep(2.6, p_base),
#                                            rep(-4.5, p_base),
#                                            rep(-2.5, p_base)),
#                             relativeTimes = rep(0, p_base))
# 
# # Means at the tips ?
# Delta <- shifts.list_to_matrix(trees[["128"]], shifts_grid[["128_15"]])
# W <- compute_actualization_matrix_ultrametric(trees[["128"]], alpha_base * diag(1, p_base, p_base))
# vec_Y <- kronecker(T_tree[["128"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
# X1.tips.exp.mat <- matrix(vec_Y, p_base, ntaxa_base) + beta_0
# unique(X1.tips.exp.mat[1, ])
# plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
#                          phylo = trees[["128"]], 
#                          params = list(shifts = list(edges = shifts_grid[["128_15"]]$edges, values = shifts_grid[["128_15"]]$values[1, ])),
#                          adj.root = 0,
#                          automatic_colors = TRUE,
#                          margin_plot = NULL,
#                          cex = 2,
#                          bg_shifts = "azure2",
#                          bg_beta_0 = "azure2")
# # Equivalent solutions ?
# extract.parsimonyNumber(parsimonyNumber(trees[["128"]], clusters_from_shifts_ism(trees[["128"]], shifts_grid[["128_15"]]$edges)))

## 32 -  3
plot(trees[["32"]], show.tip.label = FALSE); edgelabels(); tiplabels()
shifts_grid[["32_3"]] <- list(edges = c(6, 20, 48),
                                values = cbind(rep(2.4, p_base),
                                               rep(-2.4, p_base),
                                               rep(2.1, p_base)),
                                relativeTimes = rep(0, p_base))

# Means at the tips ?
Delta <- shifts.list_to_matrix(trees[["32"]], shifts_grid[["32_3"]])
W <- compute_actualization_matrix_ultrametric(trees[["32"]], alpha_base * diag(1, p_base, p_base))
vec_Y <- kronecker(T_tree[["32"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
X1.tips.exp.mat <- matrix(vec_Y, p_base, 32) + beta_0
unique(X1.tips.exp.mat[1, ])
plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
                         phylo = trees[["32"]], 
                         params = list(shifts = list(edges = shifts_grid[["32_3"]]$edges, values = shifts_grid[["32_3"]]$values[1, ])),
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
# Equivalent solutions ?
extract.parsimonyNumber(parsimonyNumber(trees[["32"]], clusters_from_shifts_ism(trees[["32"]], shifts_grid[["32_3"]]$edges)))

## 64 -  3
plot(trees[["64"]], show.tip.label = FALSE); edgelabels(); tiplabels()
shifts_grid[["64_3"]] <- list(edges = c(3, 44, 92),
                              values = cbind(rep(2.2, p_base),
                                             rep(-2.2, p_base),
                                             rep(2.2, p_base)),
                              relativeTimes = rep(0, p_base))

# Means at the tips ?
Delta <- shifts.list_to_matrix(trees[["64"]], shifts_grid[["64_3"]])
W <- compute_actualization_matrix_ultrametric(trees[["64"]], alpha_base * diag(1, p_base, p_base))
vec_Y <- kronecker(T_tree[["64"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
X1.tips.exp.mat <- matrix(vec_Y, p_base, 64) + beta_0
unique(X1.tips.exp.mat[1, ])
plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
                         phylo = trees[["64"]], 
                         params = list(shifts = list(edges = shifts_grid[["64_3"]]$edges, values = shifts_grid[["64_3"]]$values[1, ])),
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
# Equivalent solutions ?
extract.parsimonyNumber(parsimonyNumber(trees[["64"]], clusters_from_shifts_ism(trees[["64"]], shifts_grid[["64_3"]]$edges)))

## 96 -  3
plot(trees[["96"]], show.tip.label = FALSE); edgelabels(); tiplabels()
shifts_grid[["96_3"]] <- list(edges = c(48, 80, 116),
                              values = cbind(rep(2.2, p_base),
                                             rep(-2.2, p_base),
                                             rep(2.2, p_base)),
                              relativeTimes = rep(0, p_base))

# Means at the tips ?
Delta <- shifts.list_to_matrix(trees[["96"]], shifts_grid[["96_3"]])
W <- compute_actualization_matrix_ultrametric(trees[["96"]], alpha_base * diag(1, p_base, p_base))
vec_Y <- kronecker(T_tree[["96"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
X1.tips.exp.mat <- matrix(vec_Y, p_base, 96) + beta_0
unique(X1.tips.exp.mat[1, ])
plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
                         phylo = trees[["96"]], 
                         params = list(shifts = list(edges = shifts_grid[["96_3"]]$edges, values = shifts_grid[["96_3"]]$values[1, ])),
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
# Equivalent solutions ?
extract.parsimonyNumber(parsimonyNumber(trees[["96"]], clusters_from_shifts_ism(trees[["96"]], shifts_grid[["96_3"]]$edges)))

## 160 -  3
plot(trees[["160"]], show.tip.label = FALSE); edgelabels(); tiplabels()
shifts_grid[["160_3"]] <- list(edges = c(107, 62, 255),
                              values = cbind(rep(2.2, p_base),
                                             rep(-2.2, p_base),
                                             rep(2.2, p_base)),
                              relativeTimes = rep(0, p_base))

# Means at the tips ?
Delta <- shifts.list_to_matrix(trees[["160"]], shifts_grid[["160_3"]])
W <- compute_actualization_matrix_ultrametric(trees[["160"]], alpha_base * diag(1, p_base, p_base))
vec_Y <- kronecker(T_tree[["160"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
X1.tips.exp.mat <- matrix(vec_Y, p_base, 160) + beta_0
unique(X1.tips.exp.mat[1, ])
plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
                         phylo = trees[["160"]], 
                         params = list(shifts = list(edges = shifts_grid[["160_3"]]$edges, values = shifts_grid[["160_3"]]$values[1, ])),
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
# Equivalent solutions ?
extract.parsimonyNumber(parsimonyNumber(trees[["160"]], clusters_from_shifts_ism(trees[["160"]], shifts_grid[["160_3"]]$edges)))

## 160 - 7
shifts_grid[["160_7"]] <- list(edges = c(107, 62, 255,
                                         18, 204, 175, 276),
                               values = cbind(rep(2.2, p_base),
                                              rep(-2.2, p_base),
                                              rep(2.2, p_base),
                                              rep(2.2, p_base),
                                              rep(-2.2, p_base),
                                              rep(2.8, p_base),
                                              rep(2.4, p_base)),
                               relativeTimes = rep(0, p_base))

# Means at the tips ?
Delta <- shifts.list_to_matrix(trees[["160"]], shifts_grid[["160_7"]])
W <- compute_actualization_matrix_ultrametric(trees[["160"]], alpha_base * diag(1, p_base, p_base))
vec_Y <- kronecker(T_tree[["160"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
X1.tips.exp.mat <- matrix(vec_Y, p_base, 160) + beta_0
unique(X1.tips.exp.mat[1, ])
plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
                         phylo = trees[["160"]], 
                         params = list(shifts = list(edges = shifts_grid[["160_7"]]$edges, values = shifts_grid[["160_7"]]$values[1, ])),
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
# Equivalent solutions ?
extract.parsimonyNumber(parsimonyNumber(trees[["160"]], clusters_from_shifts_ism(trees[["160"]], shifts_grid[["160_7"]]$edges)))

## 160 - 11
shifts_grid[["160_11"]] <- list(edges = c(107, 62, 255,
                                         18, 204, 175, 276,
                                         145, 219, 314, 83),
                               values = cbind(rep(2.2, p_base),
                                              rep(-2.2, p_base),
                                              rep(2.2, p_base),
                                              rep(2.2, p_base),
                                              rep(-2.2, p_base),
                                              rep(2.8, p_base),
                                              rep(2.4, p_base),
                                              rep(2.8, p_base),
                                              rep(-3.1, p_base),
                                              rep(-3, p_base),
                                              rep(-2.4, p_base)),
                               relativeTimes = rep(0, p_base))

# Means at the tips ?
Delta <- shifts.list_to_matrix(trees[["160"]], shifts_grid[["160_11"]])
W <- compute_actualization_matrix_ultrametric(trees[["160"]], alpha_base * diag(1, p_base, p_base))
vec_Y <- kronecker(T_tree[["160"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
X1.tips.exp.mat <- matrix(vec_Y, p_base, 160) + beta_0
unique(X1.tips.exp.mat[1, ])
plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
                         phylo = trees[["160"]], 
                         params = list(shifts = list(edges = shifts_grid[["160_11"]]$edges, values = shifts_grid[["160_11"]]$values[1, ])),
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
# Equivalent solutions ?
extract.parsimonyNumber(parsimonyNumber(trees[["160"]], clusters_from_shifts_ism(trees[["160"]], shifts_grid[["160_11"]]$edges)))

## 160 -  15
shifts_grid[["160_15"]] <- list(edges = c(107, 62, 255,
                                          18, 204, 175, 276,
                                          145, 219, 314, 83,
                                          282, 265, 119, 47),
                                values = cbind(rep(2.2, p_base),
                                               rep(-2.2, p_base),
                                               rep(2.2, p_base),
                                               rep(2.2, p_base),
                                               rep(-2.2, p_base),
                                               rep(2.8, p_base),
                                               rep(2.4, p_base),
                                               rep(2.8, p_base),
                                               rep(-3.1, p_base),
                                               rep(-3, p_base),
                                               rep(-2.4, p_base),
                                               rep(3.2, p_base),
                                               rep(-4.5, p_base),
                                               rep(-4.8, p_base),
                                               rep(-3, p_base)),
                                relativeTimes = rep(0, p_base))

# Means at the tips ?
Delta <- shifts.list_to_matrix(trees[["160"]], shifts_grid[["160_15"]])
W <- compute_actualization_matrix_ultrametric(trees[["160"]], alpha_base * diag(1, p_base, p_base))
vec_Y <- kronecker(T_tree[["160"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
X1.tips.exp.mat <- matrix(vec_Y, p_base, 160) + beta_0
unique(X1.tips.exp.mat[1, ])
plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
                         phylo = trees[["160"]], 
                         params = list(shifts = list(edges = shifts_grid[["160_15"]]$edges, values = shifts_grid[["160_15"]]$values[1, ])),
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
# Equivalent solutions ?
extract.parsimonyNumber(parsimonyNumber(trees[["160"]], clusters_from_shifts_ism(trees[["160"]], shifts_grid[["160_15"]]$edges)))

## 192 -  3
plot(trees[["192"]], show.tip.label = FALSE); edgelabels(); tiplabels()
shifts_grid[["192_3"]] <- list(edges = c(57, 160, 307),
                               values = cbind(rep(2.2, p_base),
                                              rep(-2.2, p_base),
                                              rep(2.2, p_base)),
                               relativeTimes = rep(0, p_base))

# Means at the tips ?
Delta <- shifts.list_to_matrix(trees[["192"]], shifts_grid[["192_3"]])
W <- compute_actualization_matrix_ultrametric(trees[["192"]], alpha_base * diag(1, p_base, p_base))
vec_Y <- kronecker(T_tree[["192"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
X1.tips.exp.mat <- matrix(vec_Y, p_base, 192) + beta_0
unique(X1.tips.exp.mat[1, ])
plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
                         phylo = trees[["192"]], 
                         params = list(shifts = list(edges = shifts_grid[["192_3"]]$edges, values = shifts_grid[["192_3"]]$values[1, ])),
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
# Equivalent solutions ?
extract.parsimonyNumber(parsimonyNumber(trees[["192"]], clusters_from_shifts_ism(trees[["192"]], shifts_grid[["192_3"]]$edges)))

## 215 -  3
plot(trees[["215"]], show.tip.label = FALSE); edgelabels(); tiplabels()
shifts_grid[["215_3"]] <- list(edges = c(24, 234, 335),
                               values = cbind(rep(2.2, p_base),
                                              rep(-2.2, p_base),
                                              rep(2.2, p_base)),
                               relativeTimes = rep(0, p_base))

# Means at the tips ?
Delta <- shifts.list_to_matrix(trees[["215"]], shifts_grid[["215_3"]])
W <- compute_actualization_matrix_ultrametric(trees[["215"]], alpha_base * diag(1, p_base, p_base))
vec_Y <- kronecker(T_tree[["215"]], diag(1, p_base, p_base)) %*% W %*% as.vector(Delta)
X1.tips.exp.mat <- matrix(vec_Y, p_base, 215) + beta_0
unique(X1.tips.exp.mat[1, ])
plot.data.process.actual(Y.state = X1.tips.exp.mat[1, ],
                         phylo = trees[["215"]], 
                         params = list(shifts = list(edges = shifts_grid[["215_3"]]$edges, values = shifts_grid[["215_3"]]$values[1, ])),
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
# Equivalent solutions ?
extract.parsimonyNumber(parsimonyNumber(trees[["215"]], clusters_from_shifts_ism(trees[["215"]], shifts_grid[["215_3"]]$edges)))

## Clean up
rm(W, vec_Y, X1.tips.exp.mat, Delta)

##############
## Define date-stamp for file names
##############
datestamp_data <- format(Sys.time(), "%Y-%m-%d")

##############
## simulation function 
##############

moments_list <- vector("list") # keep the moments (avoid multiple computations)

# Return list of parameters + list of shifts + data at tips
datasetsim <- function(alpha, gamma, K, rd, rs, s, factor_shift,
                       ntaxa, NA_per, nrep, grp) {
  if (rs == 0 && s == 0){
   process_temp <- "scOU"
   alpha_mat <- alpha
  } else {
    process_temp <- "OU"
    alpha_mat <- diag(rep(alpha - rs, p_base)) + matrix(rs, p_base, p_base)
    alpha_mat <- alpha_mat + diag(c(rep(0, p_base / 2), rep(s, p_base / 2)))
  }
  tree <- trees[[paste0(ntaxa)]]
  if (K > 0){
    shifts <- shifts_grid[[paste0(ntaxa, "_", K)]]
    # Multiplicative factor 
    shifts$values <- factor_shift * shifts$values
    # Randomly multiply by +1 or -1 each trait
    plusminus <- sample(c(-1, 1), p_base, replace = TRUE)
    shifts$values <- shifts$values * plusminus
  } else {
    shifts <- NULL
  }
  var_mat <- diag(rep(2 * alpha * gamma - 2 * rd, p_base)) + matrix(2 * rd, p_base, p_base)
  var_mat <- as(var_mat, "symmetricMatrix")
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = beta_0,
                     var.root = compute_stationnary_variance(var_mat, alpha_mat))
  params <-  list(variance = var_mat,
                  root.state = root.state,
                  shifts = shifts,
                  selection.strength = alpha_mat, 
                  optimal.value = beta_0)
  XX <- simulate(phylo = tree,
                 process = process_temp,
                 p = p_base,
                 root.state = root.state, 
                 variance = var_mat,
                 shifts = shifts, 
                 selection.strength = alpha_mat, 
                 optimal.value = beta_0,
                 checks = TRUE)
  sim <- list(alpha = alpha,
              gamma = gamma,
              K = K,
              rd = rd,
              rs = rs,
              s = s,
              factor_shift = factor_shift,
              ntaxa = ntaxa,
              NA_per = NA_per,
              nrep = nrep,
              grp = grp,
              shifts = shifts,
              Y_true = extract.simulate(XX, what="states", where="tips"),
              Z_true = extract.simulate(XX, what = "states", where = "nodes"),
              m_Y_true = extract.simulate(XX, what="expectations", where="tips"))
  ## NAs
  sim$Y_data <- sim$Y_true
  if (NA_per > 0){
    nMiss <- floor(ntaxa * p_base * NA_per)
    miss <- sample(1:(p_base * ntaxa), nMiss, replace = FALSE)
    chars <- (miss - 1) %% p_base + 1
    tips <- (miss - 1) %/% p_base + 1
    for (i in 1:nMiss){
      sim$Y_data[chars[i], tips[i]] <- NA
    }
  }
  miss <- as.vector(is.na(sim$Y_data))
  Y_data_vec <- as.vector(sim$Y_data)
  Y_data_vec_known <- as.vector(sim$Y_data[!miss])
  # Vectorized Data Mask
  masque_data <- rep(FALSE, (sim$ntaxa + tree$Nnode) * p_base)
  masque_data[1:(p_base*sim$ntaxa)] <- !miss
  ## Compute true likelihood and difficulty of the problem
  attr(params, "p_dim") <- p_base
  sim$params_simu <- params
  name_config <- paste(alpha, gamma, K, rd, rs, s, factor_shift,
                       ntaxa, NA_per, sep = "_")
  if (is.null(moments_list[[name_config]])){
    moments_list[[name_config]] <<- compute_mean_variance.simple(phylo = tree,
                                    times_shared = times_shared[[paste0(ntaxa)]],
                                    distances_phylo = distances_phylo[[paste0(ntaxa)]],
                                    process = process_temp,
                                    params_old = params,
                                    masque_data = masque_data,
                                    sim = XX)
  }
  moments <- moments_list[[name_config]]
  moments$sim <- XX
  sim$log_likelihood.true <- compute_log_likelihood.simple(phylo = tree,
                                                           Y_data_vec = Y_data_vec_known,
                                                           sim = moments$sim,
                                                           Sigma = moments$Sigma,
                                                           Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                                           miss = miss,
                                                           masque_data = masque_data)
  ## Difficulty
  Sigma_YY_inv <- tcrossprod(moments$Sigma_YY_chol_inv)
  mu_0 <- (sum(Sigma_YY_inv))^(-1) * sum(Sigma_YY_inv %*% Y_data_vec_known)
  sim$difficulty <- as.vector(t(Y_data_vec_known - mu_0) %*% Sigma_YY_inv %*% (Y_data_vec_known - mu_0))
  # sim$difficulty <- as.vector(t(as.vector(sim$m_Y_true - mu_0)) %*% Sigma_YY_inv %*% as.vector(sim$m_Y_true - mu_0))
  sim$maha_data_mean <- compute_mahalanobis_distance.simple(phylo = tree,
                                                            Y_data_vec = Y_data_vec_known,
                                                            sim = moments$sim,
                                                            Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                                            miss = miss)
  return(sim)
}

##################################
## Simulations
##################################
## Set seed
set.seed(18051804)

## Sequencial simulations (for reproductability)
simlist <- foreach(i = 1:nrow(simparams)) %do% {
# simlist <- foreach(i = which(simparams$nrep == 1)) %do% {
# simlist <- foreach(i = c(5201, 5204)) %do% {
  sim <- datasetsim(alpha = simparams[i, "alpha"],
                    gamma = simparams[i, "gamma"],
                    K = simparams[i, "K"],
                    rd = simparams[i, "rd"],
                    rs = simparams[i, "rs"],
                    s = simparams[i, "s"],
                    factor_shift = simparams[i, "factor_shift"],
                    ntaxa = simparams[i, "ntaxa"],
                    NA_per = simparams[i, "NA_per"],
                    nrep= simparams[i, "nrep"],
                    grp = simparams[i, "grp"])
  sim$it <- i
  return(sim)
}

names(simlist) <- apply(simparams, 1, paste0, collapse = "_")

save.image(paste0(savedatafile, "_", datestamp_data, ".RData"))
