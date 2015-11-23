##########
## Chelonia EM
rm(list=ls())

PATH <- "../Results/Chelonia/"

# Dependencies of EM
library(ape)
library(quadrupen)
library(robustbase)
library(plyr)
# Tree Sim
library(TreeSim)
# Model selection
library(LINselect)
# Plot
library(reshape2)
library(ggplot2)

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
source("R/model_selection.R")

###################################################
## Import data
###################################################

## Import data set form package bayou
library(bayou)
data(chelonia)

tree <- chelonia$phy
data <- chelonia$dat
data <- data[match(tree$tip.label, names(data))]

K_max <- 20
data_type <- "chelonia"
ntaxa <- length(data)
K_true <- 16

## Random root ##############################################################################

## Re-scale tree to one
height_tree <- node.depth.edgelength(tree)[1]
tree$edge.length <- tree$edge.length / height_tree

alpha_grid <- find_grid_alpha(tree,
                              nbr_alpha = 20,
                              factor_up_alpha = 2,
                              factor_down_alpha = 4,
                              quantile_low_distance = 0.01,
                              log_transform = TRUE)

res <- PhyloEM(phylo = tree,
               Y_data = data,
               process = "OU", K_max = K_max,
               alpha_known = TRUE, alpha = alpha_grid[-1],
               random.root = TRUE, stationnary.root = TRUE,
               methods.segmentation = "lasso")

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_alpha_grid_stationnary_root.RData"))

params_select_stationnary <- res$alpha_max$BGH$params_select
total_time <- sum(sapply(res[-c(1, 2, 3, length(res))], function(z) z$results_summary$time))

OU_EM_stationary <- list(Nbr_shifts = length(params_select_stationnary$shifts$edges),
                         Nbr_regimes = length(params_select_stationnary$shifts$edges) + 1,
                         lnL = attr(params_select_stationnary, "log_likelihood"),
                         MlnL = NaN,
                         alpha = as.vector(params_select_stationnary$selection.strength),
                         half_life = log(2) / as.vector(params_select_stationnary$selection.strength),
                         sigma = as.vector(params_select_stationnary$variance),
                         gamma = as.vector(params_select_stationnary$variance / (2 * params_select_stationnary$selection.strength)),
                         time = total_time)

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_alpha_grid_stationnary_root.RData"))
save(params_select_stationnary, OU_EM_stationary,
     file = paste0(PATH, data_type, "_EM_stationnary_summary.RData"))

plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = res$alpha_max$BGH$params_select,
                         adj.root = 1.3,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 1.3,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2",
                         plot_ancestral_states = TRUE,
                         ancestral_states = res$alpha_max$BGH$Zhat,
                         no.margin = TRUE)

## Fixed Root ###########################################################################

## Re-scale tree to one
height_tree <- node.depth.edgelength(tree)[1]
tree$edge.length <- tree$edge.length / height_tree

alpha_grid <- find_grid_alpha(tree,
                              nbr_alpha = 50,
                              factor_up_alpha = 2,
                              factor_down_alpha = 4,
                              quantile_low_distance = 0.01,
                              log_transform = TRUE)

alpha_grid <- seq(0, 14, 0.2)

resb <- PhyloEM(phylo = tree,
               Y_data = data,
               use_previous = FALSE,
               method.init = "lasso",
               process = "OU", K_max = 10,
               alpha_known = TRUE, alpha = alpha_grid, random.root = FALSE,
               tol = list(variance = 10^(-2), 
                          value.root = 10^(-2)))

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_alpha_grid_fixed_root_rescaled_tree_lasso_init_R_init_regular_grid.RData"))

params_select_fixed <- resb$alpha_max$BGH$params_select
total_time <- sum(sapply(resb[-c(1, 2, 3, length(resb))], function(z) z$results_summary$time))

OU_EM_fixed <- list(Nbr_shifts = length(params_select_fixed$shifts$edges),
                    Nbr_regimes = length(params_select_fixed$shifts$edges) + 1,
                    lnL = attr(params_select_fixed, "log_likelihood"),
                    MlnL = NaN,
                    alpha = as.vector(params_select_fixed$selection.strength),
                    half_life = log(2) / as.vector(params_select_fixed$selection.strength),
                    sigma = as.vector(params_select_fixed$variance),
                    gamma = as.vector(params_select_fixed$variance / (2 * params_select_fixed$selection.strength)),
                    time = total_time)

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_alpha_grid_fixed_root_rescaled_tree_lasso_init_R_init_regular_grid.RData"))
save(params_select_fixed, OU_EM_fixed,
     file = paste0(PATH, data_type, "_EM_fixed_summary.RData"))

## Computation of log det and log maha for each solution
times_shared <- compute_times_ca(tree)
distances_phylo <- compute_dist_phy(tree)

compute_log_det_maha_BM <- function(param, alpha_temp){
  tree_trans <- transform_branch_length(tree, alpha_temp)
  times_shared_trans <- compute_times_ca(tree_trans)
  distances_phylo_trans <- compute_dist_phy(tree_trans)
  
  moments <- compute_mean_variance.simple(phylo = tree_trans,
                                          times_shared = times_shared_trans,
                                          distances_phylo = distances_phylo_trans,
                                          process = "BM",
                                          params_old = param,
                                          masque_data = c(rep(TRUE, ntaxa * 1),
                                                          rep(FALSE, tree$Nnode * 1)))
  
  log_det <- compute_log_det.simple(phylo = tree,
                                    Y_data_vec = as.vector(data),
                                    sim = moments$sim,
                                    Sigma = moments$Sigma,
                                    Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                    missing = rep(FALSE, ntaxa * 1),
                                    masque_data = c(rep(TRUE, ntaxa * 1),
                                                    rep(FALSE, tree$Nnode * 1)))
  
  log_maha <- compute_log_maha.simple(phylo = tree,
                                      Y_data_vec = as.vector(data),
                                      sim = moments$sim,
                                      Sigma = moments$Sigma,
                                      Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                      missing = rep(FALSE, ntaxa * 1),
                                      masque_data = c(rep(TRUE, ntaxa * 1),
                                                      rep(FALSE, tree$Nnode * 1)))
  
  return(c(log_det, log_maha, log_det + log_maha))
}

compute_log_det_maha_scOU <- function(param){
  moments <- compute_mean_variance.simple(phylo = tree,
                                          times_shared = times_shared,
                                          distances_phylo = distances_phylo,
                                          process = "scOU",
                                          params_old = param,
                                          masque_data = c(rep(TRUE, ntaxa * 1),
                                                          rep(FALSE, tree$Nnode * 1)))
  
  log_det <- compute_log_det.simple(phylo = tree,
                                    Y_data_vec = as.vector(data),
                                    sim = moments$sim,
                                    Sigma = moments$Sigma,
                                    Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                    missing = rep(FALSE, ntaxa * 1),
                                    masque_data = c(rep(TRUE, ntaxa * 1),
                                                    rep(FALSE, tree$Nnode * 1)))
  
  log_maha <- compute_log_maha.simple(phylo = tree,
                                      Y_data_vec = as.vector(data),
                                      sim = moments$sim,
                                      Sigma = moments$Sigma,
                                      Sigma_YY_chol_inv = moments$Sigma_YY_chol_inv,
                                      missing = rep(FALSE, ntaxa * 1),
                                      masque_data = c(rep(TRUE, ntaxa * 1),
                                                      rep(FALSE, tree$Nnode * 1)))
  
  return(c(log_det, log_maha, log_det + log_maha))
}

rresb <- resb[-c(1, 2, 3, 4, length(resb))]
log_lik_scOU_all <- sapply(rresb[c(1, seq(5, length(resb) - 1, 5))],
                           function(z) sapply(z$params_estim, compute_log_det_maha_scOU))
rownames(log_lik_scOU_all) <- rep(0:10, each = 3)

log_det_scOU <- log_lik_scOU_all[seq(1, ncol(log_lik_scOU_all), 3), ]
log_maha_scOU <- log_lik_scOU_all[seq(2, ncol(log_lik_scOU_all), 3), ]
log_lik_scOU <- log_lik_scOU_all[seq(3, ncol(log_lik_scOU_all), 3), ]

log_lik_BM_all <- sapply(rresb[c(1, seq(5, length(resb) - 1, 5))],
                           function(z){
                             alpha_temp <- z$params_estim$`0`$selection.strength
                             sapply(z$params_raw, compute_log_det_maha_BM, alpha_temp)
                           })
rownames(log_lik_BM_all) <- rep(0:10, each = 3)

log_det_BM <- log_lik_BM_all[seq(1, ncol(log_lik_BM_all), 3), ]
log_maha_BM <- log_lik_BM_all[seq(2, ncol(log_lik_BM_all), 3), ]
log_lik_BM <- log_lik_BM_all[seq(3, ncol(log_lik_BM_all), 3), ]

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_alpha_grid_fixed_root_rescaled_tree_lasso_init_R_init_regular_grid.RData"))

rbind(log_lik_scOU, log_lik_BM)

alpha_grid_bis <- seq(13, 14, 0.1)

resb_bis <- PhyloEM(phylo = tree,
                    Y_data = data,
                    use_previous = FALSE,
                    method.init = "lasso",
                    process = "OU", K_max = 10,
                    alpha_known = TRUE, alpha = alpha_grid_bis, random.root = FALSE,
                    tol = list(variance = 10^(-2), 
                               value.root = 10^(-2)))

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_alpha_grid_fixed_root_rescaled_tree_50_lasso_init.RData"))

resb$alpha_max$results_summary$alpha_name

plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = resb$alpha_max$BGH$params_select,
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = NULL,
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2",
                         plot_ancestral_states = TRUE,
                         ancestral_states = resb$alpha_max$BGH$Zhat)

rresb <- resb[-c(1, 2, 3, 15)]
sapply(rresb, function(z) z$results_summary$log_likelihood)[6, ]

rresb_bis <- resb_bis[-c(1, 2, 3, 15)]
sapply(rresb_bis, function(z) z$results_summary$log_likelihood)[, 11]

sapply(rresb, function(z) z$params_estim$`1`$shifts$edges)
sapply(rresb, function(z) z$params_estim$`1`$shifts$values)

## Alpha close to what found with random root
plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = resb$alpha_12.4$params_estim$`5`,
                         adj.root = 0,
                         automatic_colors = TRUE,
                         margin_plot = c(0, 0, 0, 0),
                         cex = 2,
                         bg_shifts = "azure2",
                         bg_beta_0 = "azure2")
                         # plot_ancestral_states = TRUE,
                         # ancestral_states = resb$alpha_12.4$Zhat$`5`)

resb$alpha_12.5303052395953$results_summary[c("log_likelihood_init", "log_likelihood")]

resb$alpha_12.5303052395953$params_estim$`5`$variance
resb$alpha_max$BGH$params_select$variance

resb$alpha_12.5303052395953$params_estim$`5`$root.state$value.root
resb$alpha_max$BGH$params_select$root.state$value.root


restest <- PhyloEM(phylo = tree,
                Y_data = data,
                process = "OU", K_max = 10,
                alpha_known = TRUE, alpha = alpha_grid[50], random.root = FALSE,
                tol = list(variance = 10^(-2), 
                           value.root = 10^(-2),
                           log_likelihood = 10^(-2)),
                Nbr_It_Max = 1000)
save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_alpha_fixed_root_rescaled_tree_high_alpha.RData"))

## Try mvMORPH with the shifts of high alpha
library(mvMORPH)
shifts_good <- resb$alpha_12.5303052395953$params_estim$`5`$shifts

tree_mapped <- shifts_to_simmap(tree, shifts_good$edges)
cols <- c("black", rainbow(5, start = 0, v = 0.5))
cols <- setNames(cols, 0:5)
plotSimmap(tree_mapped, colors = cols)

fit_mvOU <- mvOU(tree_mapped, data)
params_EM <- resb$alpha_12.5303052395953$params_estim$`5`

fit_mvOU$LogLik
attr(params_EM, "log_likelihood")

fit_mvOU$alpha
params_EM$selection.strength

fit_mvOU$sigma
params_EM$variance

fit_mvOU$theta
unique(compute_betas(tree, params_EM$optimal.value, params_EM$shifts))

## Other grid
alpha_grid_bis <- seq(11, 13, 0.1)

res_gridbis <- PhyloEM(phylo = tree,
                       Y_data = data,
                       process = "OU", K_max = 10,
                       alpha_known = TRUE, alpha = alpha_grid_bis, random.root = FALSE,
                       tol = list(variance = 10^(-2), 
                                  value.root = 10^(-2)))

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_alpha_grid_fixed_root_rescaled_tree_50_gridbis.RData"))

## BM ###########################################################################

res <- PhyloEM(phylo = tree,
               Y_data = data,
               process = "BM",
               K_max = K_max,
               random.root = FALSE)

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_BM.RData"))

plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = res$alpha_0.0568071425234231$params_estim$`8`,
                         adj.root = 1,
                         adj.nodes = 1,
                         #color_characters = colors_habitat,
                         #color_edges = colors_regimes,
                         regime_boxes = TRUE,
                         alpha.border = 70,
                         value_in_box = FALSE,
                         margin_plot = c(0.35, 0, 0, 0))