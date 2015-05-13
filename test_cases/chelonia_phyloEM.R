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
## Import or generate data
###################################################

## Import data set form package bayou
library(bayou)
data(chelonia)

tree <- chelonia$phy
data <- chelonia$dat
data <- data[match(tree$tip.label, names(data))]

# plot(tree, show.tip.label = FALSE)

## Fixed quantities
times_shared <- compute_times_ca(tree)
distances_phylo <- compute_dist_phy(tree)
subtree.list <- enumerate_tips_under_edges(tree)
T_tree <- incidence.matrix(tree)

K_max <- 20
data_type <- "chelonia"
ntaxa <- length(data)
K_true <- 16

# alpha <- 0.038 # Pick a "reasonable" alpha
# alpha <- 0.03204742 # K = 6 Baraud et al
# alpha <- 0.04470039 # K = 8 Birgé Massart
# data_type <- paste0(data_type, "_alpha=", alpha)

###############################################################################
## Alpha unknown
##############################################################################

estimations_several_K <- function(tree, data, K_max, order = TRUE){
  X <- list(Y_data = data,
            K_try = 0:K_max,
            ntaxa = ntaxa)
  ## Inference function, using previous as init
  fun <- function(K_t, sum_prev){
    return(estimation_wrapper.OUsr(K_t, 
                                   phylo = tree, 
                                   Y_data = data, 
                                   times_shared = times_shared, 
                                   distances_phylo = distances_phylo,
                                   subtree.list = subtree.list,
                                   T_tree = T_tree,
                                   alpha_known = FALSE,
                                   exp.root.init = sum_prev$beta_0_estim,
                                   var.init.root = sum_prev$gamma_estim,
                                   init.selection.strength = sum_prev$alpha_estim
                                   ))
  }
  ## Apply function for all K_try
  XX <- vector('list', K_max+1)
  names(XX) <- X$K_try
  # First estim
  if (order){
    K_first = 0
    K_last = K_max
    next_it <- function(K_t) { return(K_t + 1) }
    prev_it <- function(K_t) { return(K_t - 1) }
  } else {
    K_first = K_max
    K_last = 0
    next_it <- function(K_t) { return(K_t - 1) }
    prev_it <- function(K_t) { return(K_t + 1) }
  }
  XX[[paste0(K_first)]] <- estimation_wrapper.OUsr(K_first, 
                                 phylo = tree, 
                                 Y_data = data, 
                                 times_shared = times_shared, 
                                 distances_phylo = distances_phylo,
                                 subtree.list = subtree.list,
                                 T_tree = T_tree,
                                 alpha_known = FALSE)
  XX[[paste0(K_first)]]$summary <- XX[[paste0(K_first)]]$summary[,names(XX[[paste0(K_first)]]$summary)[-(21:length(names(XX[[paste0(K_first)]]$summary)))]]
  for (K_t in next_it(K_first):K_last){
    XX[[paste0(K_t)]] <- fun(K_t, XX[[paste0(prev_it(K_t))]]$summary)
    XX[[paste0(K_t)]]$summary <- XX[[paste0(K_t)]]$summary[,names(XX[[paste0(K_t)]]$summary)[-(21:length(names(XX[[paste0(K_t)]]$summary)))]]
  }
  ## Formate results
  dd <- do.call(rbind, XX)
  df <- do.call(rbind, dd[ , "summary"])
  df <- as.data.frame(df)
  df$alpha  <- X$alpha
  df$gamma  <- X$gamma
  df$K <- X$K
  df$n <- X$n
  df$ntaxa <- X$ntaxa
  df$grp <- X$grp
  df$log_likelihood_true <- X$log_likelihood.true[1]
  df$difficulty <- X$difficulty
  ## Results
  X$results_summary <- df
  X$params_estim <- dd[, "params"]
  X$params_init_estim <- dd[, "params_init"]
  X$alpha_0 <- dd[, "alpha_0"]
  X$Zhat <- dd[, "Zhat"]
  X$m_Y_estim <- dd[, "m_Y_estim"]
  X$edge.quality <- dd[, "edge.quality"]
  return(X)
}

simest_0_to_max <- estimations_several_K(tree, data, K_max, order = TRUE)
simest_max_to_0 <- estimations_several_K(tree, data, K_max, order = FALSE)

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_alpha_unknown.RData"))

###############################################################################
## Alpha known
##############################################################################
alpha_min <- min(c(simest_0_to_max$results_summary[, "alpha_estim"],
                   simest_max_to_0$results_summary[, "alpha_estim"]))
alpha_max <- max(c(simest_0_to_max$results_summary[, "alpha_estim"],
                   simest_max_to_0$results_summary[, "alpha_estim"]))
alpha_grid <- seq(from = alpha_min,
                  to = alpha_max,
                  length.out = 10)
alpha_grid <- c(alpha_grid, 1, 10)

estimations_several_K_alpha_known <- function(tree, data, K_max, known_alpha){
  X <- list(Y_data = data,
            K_try = 0:K_max,
            ntaxa = ntaxa)
  ## Inference function
  fun <- function(K_t, alpha){
    return(estimation_wrapper.OUsr(K_t, 
                                   phylo = tree, 
                                   Y_data = data, 
                                   times_shared = times_shared, 
                                   distances_phylo = distances_phylo,
                                   subtree.list = subtree.list,
                                   T_tree = T_tree,
                                   alpha_known = TRUE,
                                   alpha = alpha
    ))
  }
  ## Apply function for all K_try
  XX <- vector('list', K_max+1)
  names(XX) <- X$K_try
  for (K_t in X$K_try){
    XX[[paste0(K_t)]] <- fun(K_t, known_alpha)
    XX[[paste0(K_t)]]$summary <- XX[[paste0(K_t)]]$summary[,names(XX[[paste0(K_t)]]$summary)[-(21:length(names(XX[[paste0(K_t)]]$summary)))]]
  }
  ## Formate results
  dd <- do.call(rbind, XX)
  df <- do.call(rbind, dd[ , "summary"])
  df <- as.data.frame(df)
  df$alpha  <- X$alpha
  df$gamma  <- X$gamma
  df$K <- X$K
  df$n <- X$n
  df$ntaxa <- X$ntaxa
  df$grp <- X$grp
  df$log_likelihood_true <- X$log_likelihood.true[1]
  df$difficulty <- X$difficulty
  ## Results
  X$results_summary <- df
  X$params_estim <- dd[, "params"]
  X$params_init_estim <- dd[, "params_init"]
  X$alpha_0 <- dd[, "alpha_0"]
  X$Zhat <- dd[, "Zhat"]
  X$m_Y_estim <- dd[, "m_Y_estim"]
  X$edge.quality <- dd[, "edge.quality"]
  return(X)
}

simests_alpha_known <- vector("list", length(alpha_grid))
for (i in 1:length(alpha_grid)){
  known_alpha = alpha_grid[i]
  simests_alpha_known[[i]] <- estimations_several_K_alpha_known(tree, data, K_max, known_alpha)
}

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_several_tries.RData"))

###############################################################################
## Baraud Giraud Huet
###############################################################################
load(paste0(PATH, "migale/", data_type, "_estimation_K_max=", K_max, "_several_tries.RData"))

## Compute penalty and criteria for each set of parameters
penalty <- function(K_try, complexity, ntaxa){
  return(1/2 * penalty_BaraudGiraudHuet_likelihood(K_try, 
                                                   complexity, 
                                                   ntaxa, 
                                                   C = 1.1))
} 

criteria <- function(ll, pen){
  return(-ll + pen)
}

selected_K <- function(crit_ll, K_try){
  K_try[which.min(crit_ll)]
}

add_crit_and_pen <- function(z){
  z$pen_ll <- penalty(z$K_try,
                      z$complexity,
                      ntaxa)
  z$crit_ll <- criteria(z$log_likelihood,
                        z$pen_ll)
  z$K_select <- z$K_try[which.min(z$crit_ll)]
  z$ll_select <- z$log_likelihood[which.min(z$crit_ll)]
  return(z)
}

add_K_select_to_list <- function(simest){
  simest$results_summary <- add_crit_and_pen(simest$results_summary)
  simest$K_select <- unique(simest$results_summary$K_select)
  simest$params_select <- simest$params_estim[[paste(simest$K_select)]]
  return(simest)
}

fun <- function(z){
  z$results_summary$alpha <- signif(z$results_summary$alpha_estim, 2)
  return(z) 
}

simest_0_to_max <- add_K_select_to_list(simest_0_to_max)
simest_0_to_max$results_summary$alpha <- "K-1"
simest_max_to_0 <- add_K_select_to_list(simest_max_to_0)
simest_max_to_0$results_summary$alpha <- "K+1"
simests_alpha_known <- lapply(simests_alpha_known, add_K_select_to_list)
simests_alpha_known <- lapply(simests_alpha_known, fun)

simests_all <- simests_alpha_known
simests_all[["K-1"]] <- simest_0_to_max
simests_all[["K+1"]] <- simest_max_to_0

###################################################################
## Plots log likelihood
###################################################################
extract_data_frame <- function(simests){
  dd <- do.call(rbind, simests)
  results_summary <- as.data.frame(do.call(rbind, dd[,"results_summary"]))
  return(results_summary)
}

summary_alpha_known <- extract_data_frame(simests_all[c(1, 3, 5, 7, 10, 13, 14)])
p <- ggplot(summary_alpha_known, aes(x = K_try, y = log_likelihood, color = as.factor(alpha)))
p <- p + geom_line()
# p <- p + geom_point(data = crit_min, aes(x = K_try, y = crit_ll, size = 5))
p <- p + geom_point(aes(x = K_select, y = ll_select, color = as.factor(alpha), size = 5))
# p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = "K",
              y = "Log Likelihood")
p <- p + scale_size(name = "", labels = "Min")
p <- p + scale_color_discrete(name = expression(alpha))
#                               labels = unique(summary_alpha_known$alpha)
p <- p + scale_x_continuous(breaks = c(0, 5, 10, 15, 20))
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

summary_alpha_known <- extract_data_frame(simests_alpha_known)
p <- ggplot(summary_alpha_known, aes(x = as.factor(alpha_estim), y = K_select))
p <- p + geom_point()
# p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = expression(alpha),
              y = "K select")
p <- p + scale_x_discrete(labels = signif(unique(summary_alpha_known$alpha_estim), 2))
p <- p + scale_y_continuous(breaks = unique(summary_alpha_known$K_select))
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

summary_all <- extract_data_frame(simests_all)
p <- ggplot(summary_all, aes(x = K_try, y = time, color = as.factor(alpha)))
p <- p + geom_line()
# p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = "K",
              y = "Time (s)")
p <- p + scale_x_discrete(labels = signif(unique(summary_alpha_known$alpha_estim), 2))
p <- p + scale_y_continuous(breaks = unique(summary_alpha_known$K_select))
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# ll_plot <- melt(simest$results_summary[,c("K_try", "log_likelihood", "crit_ll", "pen_ll")],
#                 id.vars = "K_try",
#                 variable.name = "score",
#                 value.name = "value")
# ll_min_plot <- melt(crit_min[,c("K_try", "log_likelihood", "crit_ll", "pen_ll")],
#                     id.vars = "K_try", 
#                     variable.name = "score",
#                     value.name = "value")
# 
# ## Plot Likelihood
# p <- ggplot(ll_plot, aes(x = K_try, y = value, color = score))
# p <- p + geom_point()
# p <- p + geom_point(data = ll_min_plot, aes(x = K_try, y = value, color = score, size = 5))
# # p <- p + geom_vline(xintercept = K_true)
# p <- p + labs(x = "K",
#               y = "Log Likelihood")
# p <- p + theme_bw()
# p <- p + scale_size(name = "", labels = "Min")
# p <- p + scale_x_continuous(breaks = c(0, 1, 16))
# p <- p + theme(axis.text = element_text(size = 12),
#                strip.text = element_text(size = 12)
#                ##legend.position = c(0, 1),
#                ##legend.justification = c(0, 1)
# )
# p

##########################################################
## Plot Processes
##########################################################
nbrSol <- length(simests_all)
nbrLignes <- (nbrSol %/% 3) + 1
if (nbrSol %% 3 == 0) nbrLignes <- nbrLignes - 1
scr <- split.screen(c(nbrLignes, 3))
for (sol in 1:(nbrSol)) {
  ## Plot
  screen(scr[sol])
  plot.data.process.actual(Y.state = data,
                           phylo = tree, 
                           params = simests_all[[sol]]$params_select,
                           adj = 2,
                           automatic_colors = TRUE)
  legend("topleft", legend = unique(simests_all[[sol]]$results_summary$alpha), cex = 0.5)
}
close.screen(all.screens = TRUE)


###########################################################
# Summary results
###########################################################
K_select <- unique(simest$results_summary$K_select)
params_select <- simest$params_estim[[paste(K_select)]]
summary_select <- subset(simest$results_summary, K_try == K_select)

OU_EMselect <- list(Nbr_shifts = K_select,
                 Nbr_regimes = K_select + 1,
                 lnL = summary_select$log_likelihood,
                 MlnL = NA,
                 alpha = summary_select$alpha_estim,
                 half_life = log(2)/summary_select$alpha_estim,
                 sigma = 2*summary_select$alpha_estim*summary_select$gamma_estim,
                 gamma = summary_select$gamma_estim)

plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = params_select,
                         adj = 2,
                         automatic_colors = TRUE)

K_true <- 16
params_true <- simest$params_estim[[paste(K_true)]]
summary_true <- subset(simest$results_summary, K_try == K_true)

OU_EMtrue <- list(Nbr_shifts = K_true,
                    Nbr_regimes = K_true + 1,
                    lnL = summary_true$log_likelihood,
                    MlnL = NA,
                    alpha = summary_true$alpha_estim,
                    half_life = log(2)/summary_true$alpha_estim,
                    sigma = 2*summary_true$alpha_estim*summary_true$gamma_estim,
                    gamma = summary_true$gamma_estim)

plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = simest$params_estim[['2']],
                         adj = 2,
                         automatic_colors = TRUE)

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_reverse_order.RData"))

########################################################################
## Equivalent Solutions 
########################################################################
times_shared <- compute_times_ca(tree)
t_tree <-  min(node.depth.edgelength(tree)[1:ntaxa])
T_tree <- incidence.matrix(tree)
ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                                                  selection.strength = params_select$selection.strength,
                                                  times_shared = times_shared)
T_tree_ac <- T_tree * ac_tree

eq_shifts_edges_K_select <- equivalent_shifts_edges(tree, 
                                                    params_select$shifts$edges)
eq_shifts_values_K_select <- equivalent_shifts_values(tree,
                                                      shifts = params_select$shifts,
                                                      beta_0 = params_select$optimal.value,
                                                      eq_shifts_edges = eq_shifts_edges_K_select,
                                                      T_tree_ac = T_tree_ac)

plot_equivalent_shifts.actual(tree, eq_shifts_edges_K_select, eq_shifts_values_K_select, use.edge.length = FALSE, adj = 0)

###########################################################################
## Other Guesses
###########################################################################
simest2 <- estimation_wrapper.OUsr(1, 
                                   phylo = tree, 
                                   Y_data = data, 
                                   times_shared = times_shared, 
                                   distances_phylo = distances_phylo,
                                   subtree.list = subtree.list,
                                   T_tree = T_tree,
                                   alpha_known = FALSE)

plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = simest2$params,
                         adj = 2,
                         automatic_colors = TRUE)

simest3 <- estimation_wrapper.OUsr(2, 
                        phylo = tree, 
                        Y_data = data, 
                        times_shared = times_shared, 
                        distances_phylo = distances_phylo,
                        subtree.list = subtree.list,
                        T_tree = T_tree,
                        alpha_known = FALSE,
                        method.init.alpha = "default",
                        exp.root.init = simest$results_summary[21, "beta_0_estim"],
                        var.init.root = simest$results_summary[21, "gamma_estim"],
                        init.selection.strength = simest$results_summary[21, "alpha_estim"],
                        method.init = "default",
                        edges.init = c(77, 382),
                        relativeTimes.init = c(0,0))

plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = simest3$params,
                         adj = 2,
                         automatic_colors = TRUE)

simest4 <- estimation_wrapper.OUsr(2, 
                                   phylo = tree, 
                                   Y_data = data, 
                                   times_shared = times_shared, 
                                   distances_phylo = distances_phylo,
                                   subtree.list = subtree.list,
                                   T_tree = T_tree,
                                   alpha_known = FALSE,
                                   method.init.alpha = "default",
                                   optimal.value.init = simest$results_summary[3, "beta_0_estim"],
                                   var.init.root = simest$results_summary[3, "gamma_estim"],
                                   init.selection.strength = simest$results_summary[3, "alpha_estim"],
                                   method.init = "default",
                                   edges.init = c(47,382),
                                   relativeTimes.init = c(0,0),
                                   tol_h_l = 10^(-3))

simest5 <- estimation_wrapper.OUsr(2, 
                                   phylo = tree, 
                                   Y_data = data, 
                                   times_shared = times_shared, 
                                   distances_phylo = distances_phylo,
                                   subtree.list = subtree.list,
                                   T_tree = T_tree,
                                   alpha_known = TRUE,
                                   alpha = log(2)/17.6,
                                   method.init = "default",
                                   edges.init = c(47,382),
                                   relativeTimes.init = c(0,0))