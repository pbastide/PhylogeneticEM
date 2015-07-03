##########
## Test Cases EM

rm(list=ls())

PATH <- "../Results/Test_Cases/"

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
## Path quantities
###################################################
data_type <- "bimac" # chelonia bimac
K_max <- 5 # 20 5 

load(paste0(PATH, data_type, "_PhyloEM.RData"))

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

# Import data set form package geiger
library(geiger)
data(primates)
data_type <- "primates"
tree <- primates$phy
data <- primates$dat

# library(ouch)
# source("R/convert.R")
# # Anolis bimaculatus lizard size data
# data(bimac)
# data_type <- "bimac"
# bimac_ouch_tree <- with(bimac, ouchtree(node, ancestor, time, species))
# tree <- convert.pmc(bimac_ouch_tree)
# data <- log(bimac$size[23:45])
# names(data) <- bimac$species[23:45]

data <- data[match(tree$tip.label, names(data))]

plot(tree, show.tip.label = FALSE)

## Fixed quantities
times_shared <- compute_times_ca(tree)
distances_phylo <- compute_dist_phy(tree)
subtree.list <- enumerate_tips_under_edges(tree)
T_tree <- incidence.matrix(tree)

K_max <- 20
ntaxa <- length(data)

###############################################################################
## Alpha unknown
##############################################################################

estimations_several_K <- function(tree, data, K_max, order = TRUE, prev_init = TRUE){
  X <- list(Y_data = data,
            K_try = 0:K_max,
            ntaxa = ntaxa)
  if (prev_init){
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
                                     method.init.alpha = "default",
                                     exp.root.init = sum_prev$beta_0_estim,
                                     var.init.root = sum_prev$gamma_estim,
                                     init.selection.strength = sum_prev$alpha_estim
      ))
    }
  } else {
    fun <- function(K_t, ...){
      return(estimation_wrapper.OUsr(K_t, 
                                     phylo = tree, 
                                     Y_data = data, 
                                     times_shared = times_shared, 
                                     distances_phylo = distances_phylo,
                                     subtree.list = subtree.list,
                                     T_tree = T_tree,
                                     alpha_known = FALSE))
    }
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
    cat(paste(K_t))
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
simest_init_reg <- estimations_several_K(tree, data, K_max, prev_init = FALSE)

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
alpha_grid <- c(0, alpha_grid, 5)

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
  cat(paste(i))
  simests_alpha_known[[i]] <- estimations_several_K_alpha_known(tree, data, K_max, known_alpha)
}

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_several_tries.RData"))

###############################################################################
## Baraud Giraud Huet
###############################################################################
# load(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_alpha_unknown_init_reg.RData"))
# simest_init_reg <- simest_0_to_max
load(paste0(PATH, data_type, "_estimation_K_max=", K_max, "_several_tries.RData"))
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

alpha_type <- function(simest, bool){
  simest$results_summary$alpha_known <- bool
  return(simest)
}

fun <- function(z){
  z$results_summary$alpha <- signif(z$results_summary$alpha_estim, 2)
  return(z) 
}

simest_0_to_max <- add_K_select_to_list(simest_0_to_max)
simest_0_to_max$results_summary$alpha <- "K-1"
simest_0_to_max <- alpha_type(simest_0_to_max, FALSE)
simest_max_to_0 <- add_K_select_to_list(simest_max_to_0)
simest_max_to_0$results_summary$alpha <- "K+1"
simest_max_to_0 <- alpha_type(simest_max_to_0, FALSE)
simest_init_reg <- add_K_select_to_list(simest_init_reg)
simest_init_reg$results_summary$alpha <- "reg"
simest_init_reg <- alpha_type(simest_init_reg, FALSE)
simest_alpha_unknown <- list("K-1" = simest_0_to_max,
                             "K+1" = simest_max_to_0,
                             "reg" = simest_init_reg)
simests_alpha_known <- c(simests_alpha_known)#, simests_alpha_known_bis)
simests_alpha_known <- lapply(simests_alpha_known, add_K_select_to_list)
simests_alpha_known <- lapply(simests_alpha_known, fun)
simests_alpha_known <- lapply(simests_alpha_known, alpha_type, bool = TRUE)

simests_all <- vector("list")
simests_all <- simests_alpha_known
simests_all[["K-1"]] <- simest_0_to_max
simests_all[["K+1"]] <- simest_max_to_0
simests_all[["reg"]] <- simest_init_reg

# alpha_grid <- c(alpha_grid, -1, -2)

###################################################################
## Plots log likelihood
###################################################################
extract_data_frame <- function(simests){
  dd <- do.call(rbind, simests)
  results_summary <- as.data.frame(do.call(rbind, dd[,"results_summary"]))
  return(results_summary)
}

summary_all <- extract_data_frame(simests_all)
p <- ggplot(summary_all, aes(x = K_try, y = log_likelihood, color = as.factor(alpha), group = as.factor(alpha), linetype = alpha_known))
p <- p + geom_line()
# p <- p + geom_point(data = crit_min, aes(x = K_try, y = crit_ll, size = 5))
p <- p + geom_point(aes(x = K_select, y = ll_select, color = as.factor(alpha), size = 5))
# p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = "K",
              y = "Log Likelihood")
p <- p + scale_size(name = "", labels = "Min Criterion")
p <- p + scale_linetype(name = expression(alpha), 
                        labels = c("Estimated", "Fixed"))
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
p <- ggplot(summary_alpha_known, aes(x = alpha_estim, y = K_select))
p <- p + geom_point()
# p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = expression(alpha),
              y = "K select")
#p <- p + scale_x_discrete(labels = signif(unique(summary_alpha_known$alpha_estim), 2))
p <- p + scale_y_continuous(breaks = unique(summary_alpha_known$K_select))
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

summary_alpha_known <- extract_data_frame(simests_all)
p <- ggplot(summary_alpha_known, aes(x = K_try, y = alpha_estim,color = as.factor(alpha), group = as.factor(alpha)))
p <- p + geom_point()
# p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = "K",
              y = expression(alpha))
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

##########################################################
## Max Log Likelihood
##########################################################
select_grid <- function(simests, alpha_grid){
  summary_alpha_known <- extract_data_frame(simests)
  summary_max_ll <- ddply(summary_alpha_known,
                          .(K_try),
                          summarize,
                          log_likelihood = max(log_likelihood),
                          pen_ll = unique(pen_ll),
                          crit_ll = min(crit_ll),
                          alpha = "max_grid")
  K_select <- which.min(summary_max_ll$crit_ll) - 1
  summary_max_ll$K_select <- K_select
  summary_max_ll$ll_select <- summary_max_ll$log_likelihood[which.min(summary_max_ll$crit_ll)]
  index_select <- which(summary_alpha_known$crit_ll == min(summary_max_ll$crit_ll))
  index_select <- which(alpha_grid == summary_alpha_known[index_select, "alpha_estim"])
  params_select <- simests[[index_select]]$params_estim[[K_select + 1]]
  return(list(summary_max_ll = summary_max_ll,
              params_select = params_select,
              total_time = sum(summary_alpha_known$time),
              simest_select = simests[[index_select]]))
}

max_ll_select_grid <- select_grid(simests_all[2:length(alpha_grid)], alpha_grid[2:length(alpha_grid)])
params_select <- max_ll_select_grid$params_select

# # Chelonia only (smaller grid)
# max_ll_params_tres_grosse <- select_grid(simests_all[c(11, seq(21, 109, 20))], alpha_grid[c(11, seq(21, 109, 20))])
# params_select_tres_grosse <- max_ll_params_tres_grosse$params_select

plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = params_select,
                         automatic_colors = TRUE)

summary_grid_EM <- rbind(max_ll_select_grid$summary_max_ll,
                         extract_data_frame(simest_alpha_unknown)[,c("K_try",
                                                                "log_likelihood",
                                                                "pen_ll",
                                                                "crit_ll",
                                                                "alpha",
                                                                "K_select",
                                                                "ll_select")])
p <- ggplot(summary_grid_EM, aes(x = K_try, y = log_likelihood, color = as.factor(alpha), group = as.factor(alpha)))
p <- p + geom_line()
# p <- p + geom_point(data = crit_min, aes(x = K_try, y = crit_ll, size = 5))
p <- p + geom_point(aes(x = K_select, y = ll_select, color = as.factor(alpha), size = 5))
# p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = "K",
              y = "Log Likelihood")
p <- p + scale_size(name = "", labels = "Min")
p <- p + scale_color_discrete(name = expression(alpha))
#                               labels = unique(summary_alpha_known$alpha)
p <- p + scale_x_continuous(breaks = 0:5)
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

summary_all <- extract_data_frame(simests_all)
summary_alpha_known <- subset(summary_all, alpha_known == TRUE)
p <- ggplot(summary_alpha_known,
            aes(x = K_try, y = log_likelihood,
                color = alpha_estim, group = as.factor(alpha)))
p <- p + geom_line()
p <- p + geom_point(aes(x = K_select, y = ll_select, color = alpha_estim, size = 5))
## Plot Alpha Estim
p <- p + geom_line(data = subset(summary_all, alpha_known == FALSE),
                   mapping = aes(x = K_try, y = log_likelihood, group = as.factor(alpha)),
                   color = "lightgreen",
                   size = 1.5)
## Plot Maximum
p <- p + stat_summary(mapping = aes(x = K_try, y = log_likelihood, group = NULL),
                      data = summary_alpha_known,
                      fun.y = max,
                      geom="line", color = "red", size = 2, alpha = 0.7)
p <- p + geom_point(x = K_select_all,
                    y = summary_max_ll[K_select_all + 1, "log_likelihood"],
                    color = "red", size = 4)
p <- p + labs(x = "K",
              y = "Log Likelihood")
p <- p + scale_size(name = "", labels = "Min")
p <- p + scale_color_continuous(name = expression(alpha))
p <- p + scale_x_continuous(breaks = c(0, 5, 10, 15, 20))
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

##########################################################
## Plot Processes
##########################################################
simests_plot <- simests_all[c(11:20, 101)]
nbrSol <- length(simests_plot)
nbrLignes <- (nbrSol %/% 3) + 1
if (nbrSol %% 3 == 0) nbrLignes <- nbrLignes - 1
scr <- split.screen(c(nbrLignes, 3))
for (sol in 1:(nbrSol)) {
  ## Plot
  screen(scr[sol])
  plot.data.process.actual(Y.state = data,
                           phylo = tree, 
                           params = simests_plot[[sol]]$params_select,
                           adj = 2,
                           automatic_colors = TRUE)
  legend("topleft", legend = unique(simests_plot[[sol]]$results_summary$alpha), cex = 0.5)
}
close.screen(all.screens = TRUE)


###########################################################
# Summary results
###########################################################
# K_select <- unique(simest$results_summary$K_select)
# params_select <- simest$params_estim[[paste(K_select)]]
# summary_select <- subset(simest$results_summary, K_try == K_select)

# summary_alpha_known <- extract_data_frame(simests_all)
# K_select_all <- which.min(summary_max_ll$crit_ll) - 1
# index_select <- which(summary_alpha_known$crit_ll == min(summary_max_ll$crit_ll))
# index_select <- which(alpha_grid == summary_alpha_known[index_select, "alpha_estim"])
# params_select <- simests_all[[index_select]]$params_estim[[K_select_all]]
# summary_select <- simests_all[[index_select]]$results_summary[K_select_all +1, ]

## Bimac only (init fait son job)
params_select <- simest_init_reg$params_select

OU_EMselect <- list(Nbr_shifts = length(params_select$shifts$edges),
                    Nbr_regimes = length(params_select$shifts$edges) + 1,
                    lnL = attr(params_select, "log_likelihood")[1],
                    MlnL = NaN,
                    alpha = params_select$selection.strength,
                    half_life = log(2)/params_select$selection.strength,
                    sigma = params_select$variance,
                    gamma = params_select$root.state$var.root)#,
                    # time = sum(simest_init_reg$results_summary$time))

## Computation times
# Chelonia
OU_EMselect$time <- max_ll_params_tres_grosse$total_time
# Bimac
OU_EMselect$time <- sum(simest_init_reg$results_summary$time)

plot.data.process.actual(Y.state = data,
                         phylo = tree, 
                         params = params_select,
                         adj.root = 1,
                         adj.nodes = 0,
                         automatic_colors = TRUE)

# K_true <- 16
# # params_true <- simest$params_estim[[paste(K_true)]]
# # summary_true <- subset(simest$results_summary, K_try == K_true)
# index_true <- which(summary_alpha_known$crit_ll == summary_max_ll$crit_ll[K_true + 1])
# index_true <- which(alpha_grid == summary_alpha_known[index_true, "alpha_estim"])
# params_true <- simests_all[[index_true]]$params_estim[[K_true]]
# summary_true <- simests_all[[index_true]]$results_summary[K_true +1, ]
# 
# OU_EMtrue <- list(Nbr_shifts = K_true,
#                   Nbr_regimes = K_true + 1,
#                   lnL = summary_true$log_likelihood,
#                   MlnL = NA,
#                   alpha = summary_true$alpha_estim,
#                   half_life = log(2)/summary_true$alpha_estim,
#                   sigma = 2*summary_true$alpha_estim*summary_true$gamma_estim,
#                   gamma = summary_true$gamma_estim)
# 
# plot.data.process.actual(Y.state = data,
#                          phylo = tree, 
#                          params = params_true,
#                          adj = 2,
#                          automatic_colors = TRUE)

save.image(paste0(PATH, data_type, "_PhyloEM.RData"))
save(tree, data, alpha_grid, simest_init_reg, times_shared, K_max,
     #max_ll_params_tres_grosse, params_select_tres_grosse, # chelonia only
     params_select, OU_EMselect,
     file = paste0(PATH, data_type, "_PhyloEM_summary.RData"))

########################################################################
## Equivalent Solutions 
########################################################################
times_shared <- compute_times_ca(tree)
t_tree <-  min(node.depth.edgelength(tree)[1:ntaxa])
T_tree <- incidence.matrix(tree)
ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                                                  selection.strength = max_ll_select_grid$params_select$selection.strength,
                                                  times_shared = times_shared)
T_tree_ac <- T_tree * ac_tree

eq_shifts_edges_K_select <- equivalent_shifts_edges(tree, 
                                                    max_ll_select_grid$params_select$shifts$edges)
eq_shifts_values_K_select <- equivalent_shifts_values(tree,
                                                      shifts = max_ll_select_grid$params_select$shifts,
                                                      beta_0 = max_ll_select_grid$params_select$optimal.value,
                                                      eq_shifts_edges = eq_shifts_edges_K_select,
                                                      T_tree_ac = T_tree_ac)

plot_equivalent_shifts.actual(tree, eq_shifts_edges_K_select, eq_shifts_values_K_select, use.edge.length = FALSE, adj = 0)
