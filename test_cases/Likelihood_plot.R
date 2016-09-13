##########
## Test Case : Varying K in the inference.

rm(list=ls())

PATH <- "../Results/Likelihood_Plot/"

reqpckg <- c("ape", "quadrupen", "robustbase")

# Dependencies of EM
library(ape)
library(quadrupen)
library(robustbase)
library(plyr)
# Plot
library(ggplot2)
library(reshape2)
library(grid)
# Parallel Computation
library(doParallel)
library(foreach)
# Tree Sim
library(TreeSim)
# Model selection
library(LINselect)
library(capushe)

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
source("R/model_selection")

Ncores <- 3

###################################################
## Import or generate data
###################################################

## Import data set form package geiger
library(geiger)
data(chelonia)

tree <- chelonia$phy
data <- chelonia$dat
data <- data[match(tree$tip.label, names(data))]

plot(tree, show.tip.label = FALSE)

Ks <- 0:30
data_type <- "chelonia"
ntaxa <- length(data)
K_true <- 16

#alpha <- 0.15 # Pick a "reasonable" alpha
#data_type <- paste0(data_type, "_alpha=", alpha)

## Random data
# set.seed(20141211)
# ntaxa <- 64
# lambda <- 0.1
# tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = lambda, mu = 0, age = 1, mrca = TRUE)[[1]]
# plot(tree, show.tip.label = FALSE)
# 
# beta_0 <- 0
# gamma <- 0.1
# alpha <- 3
# 
# K_true <- 10
# shifts <- sample_shifts(tree, 18, K_true)
# 
# Ks <- 1:50
# data_type <- paste0("random_ntaxa=", ntaxa, "K_true=", K_true)
# 
# XX <- simulate(phylo = tree,
#                process = "OU",
#                root.state = list(random = TRUE,
#                                  stationary.root = TRUE,
#                                  value.root = NA,
#                                  exp.root = beta_0,
#                                  var.root = gamma), 
#                variance = 2*alpha*gamma,
#                shifts = shifts, 
#                selection.strength = alpha, 
#                optimal.value = beta_0)
# 
# data = extract.simulate(XX, what="states", where="tips")

## Easy data
# set.seed(20141211)
# ntaxa <- 64
# lambda <- 0.1
# tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = lambda, mu = 0, age = 1, mrca = TRUE)[[1]]
# plot(tree, show.tip.label = FALSE); edgelabels();
# 
# beta_0 <- 0
# gamma <- 0.1
# alpha <- 3
# 
# K_true <- 5
# shifts <- list(edges = c(61, 96, 87, 11, 30),
#                values = c(2, -2, 2, -2, 5),
#                relativeTimes = rep(0, K_true))
# 
# Ks <- 0:35
# data_type <- paste0("easy_ntaxa=", ntaxa, "K_true=", K_true)
# 
# XX <- simulate(phylo = tree,
#                process = "OU",
#                root.state = list(random = TRUE,
#                                  stationary.root = TRUE,
#                                  value.root = NA,
#                                  exp.root = beta_0,
#                                  var.root = gamma), 
#                variance = 2*alpha*gamma,
#                shifts = shifts, 
#                selection.strength = alpha, 
#                optimal.value = beta_0)
# 
# data = extract.simulate(XX, what="states", where="tips")
# 
# plot.process.actual(Y.state = extract.simulate(XX, what="states", where="tips"),
#                     Z.state = extract.simulate(XX, what="states", where="nodes"),
#                     phylo = tree, 
#                     paramsEstimate = list(shifts = shifts, optimal.value = beta_0))

## Easy data big tree
# set.seed(20141211)
# ntaxa <- 256
# lambda <- 0.1
# tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = lambda, mu = 0, age = 1, mrca = TRUE)[[1]]
# plot(tree, show.tip.label = FALSE); edgelabels();
# 
# beta_0 <- 0
# gamma <- 0.1
# alpha <- 3
# 
# K_true <- 5
# shifts <- list(edges = c(396, 314, 217, 80, 19),
#                values = c(2, -2, 2, -2, 5),
#                relativeTimes = rep(0, K_true))
# 
# Ks <- 0:35
# data_type <- paste0("easy_ntaxa=", ntaxa, "K_true=", K_true)
# 
# XX <- simulate(phylo = tree,
#                process = "OU",
#                root.state = list(random = TRUE,
#                                  stationary.root = TRUE,
#                                  value.root = NA,
#                                  exp.root = beta_0,
#                                  var.root = gamma), 
#                variance = 2*alpha*gamma,
#                shifts = shifts, 
#                selection.strength = alpha, 
#                optimal.value = beta_0)
# 
# data = extract.simulate(XX, what="states", where="tips")
# 
# plot.process.actual(Y.state = extract.simulate(XX, what="states", where="tips"),
#                     Z.state = extract.simulate(XX, what="states", where="nodes"),
#                     phylo = tree, 
#                     paramsEstimate = list(shifts = shifts, optimal.value = beta_0))

###############################################################################
## EM for several values of K
##############################################################################

estimationfunction <- function(K_t, alpha_known = FALSE, alpha = 0) {
  ## If an estimation fails, catch error with "try" and try again
  results_estim_EM <- estimateEM(phylo=tree, 
                                 Y_data = data, 
                                 tol = list(variance=10^(-4), 
                                            value.root=10^(-4), 
                                            exp.root=10^(-4), 
                                            var.root=10^(-4), 
                                            selection.strength=10^(-3)), 
                                 process = "OU", 
                                 method.variance = "simple", 
                                 method.init = "lasso",
                                 method.init.alpha = "estimation",
                                 Nbr_It_Max = 1000, 
                                 nbr_of_shifts = K_t, 
                                 alpha_known = alpha_known, ##
                                 known.selection.strength = alpha,
                                 min_params=list(variance = 10^(-4), 
                                                 value.root = -10^(4), 
                                                 exp.root = -10^(4), 
                                                 var.root = 10^(-4),
                                                 selection.strength = 10^(-4)),
                                 max_params=list(variance = 10^(4), 
                                                 value.root = 10^(4), 
                                                 exp.root = 10^(4), 
                                                 var.root = 10^(4),
                                                 selection.strength = 10^(4)),
                                 methods.segmentation = c("lasso", "best_single_move"))
  extract.edges <- function(x) {
    z <- unlist(lapply(x, function(y) y$shifts$edges))
    if (!is.null(z)) z <- matrix(z, nrow = K_t)      
    return(z)
  }
  selected.edges <- extract.edges(results_estim_EM$params_history)
  params <- results_estim_EM$params
  X <- NULL
  X$alpha_estim = params$selection.strength
  X$gamma_estim = params$root.state$var.root
  X$shifts_estim = params$shifts
  X$EM_steps <- attr(results_estim_EM, "Nbr_It")
  X$DV_estim <- attr(results_estim_EM, "Divergence")
  X$CV_estim <- (attr(results_estim_EM, "Nbr_It") != 1000) && !X$DV_estim
  X$Zhat <- results_estim_EM$ReconstructedNodesStates
  compute.quality <- function(i) {
    res <- 0 + (selected.edges == i)
    return(mean(colSums(res)))
  }
  if (K_t != 0){
    edge.quality <- unlist(lapply(params$shifts$edges, compute.quality))
    names(edge.quality) <- params$shifts$edges
  } else {
    edge.quality <- NA
  }
  X$edge.quality <- edge.quality
  X$start <- results_estim_EM$params_history["0"]
  X$beta_0_estim <- params$root.state$exp.root
  X$log_likelihood <- attr(params, "log_likelihood")[1]
  X$mahalanobis_distance_data_mean <- attr(params, "mahalanobis_distance_data_mean")
  X$number_new_shifts <- results_estim_EM$number_new_shifts
  X$mean_number_new_shifts <- mean(results_estim_EM$number_new_shifts)
  X$number_equivalent_solutions <- results_estim_EM$number_equivalent_solutions
  X$K_try <- K_t
  X$complexity <- choose(2*ntaxa-2-K_t, K_t)
  return(X)
}



## Parallelized estimations
cl <- makeCluster(Ncores)
registerDoParallel(cl)
estimations <- foreach(i = Ks, .packages = reqpckg) %dopar%
{
  estimationfunction(i)
}
stopCluster(cl)

save.image(paste0(PATH, data_type, "_estimation_K_in_", paste(Ks, collapse = "_"), ".RData"))

####################################################
## Arrange Data 
####################################################

## Arange data
dd <- do.call(rbind, estimations)
df <- apply(dd[ , colnames(dd) %in% c("alpha", "gamma", "K_try", "n",
                                              "grp", "alpha_estim",
                                              "gamma_estim", "EM_steps", 
                                              "DV_estim", "CV_estim",
                                              "difficulty", "beta_0_estim",
                                              "log_likelihood",
                                              "mahalanobis_distance_data_mean",
                                              "mean_number_new_shifts",
                                              "number_equivalent_solutions")],
            2, unlist)
df <- as.data.frame(df)

## Normalize mahalanobis distance (stationary root case)
df$least_squares <- df$mahalanobis_distance_data_mean * df$gamma_estim

## Model Complexity
#model_complexity <- sapply(Ks, function(z) extract.partitionsNumber(partitionsNumber(tree, z + 1)))
model_complexity <- sapply(Ks, function(z) choose(2*ntaxa-2-z, z))
naive_complexity <- sapply(Ks, function(z) choose(2*ntaxa-2, z))
df[["model_complexity"]] <- model_complexity
df[["naive_complexity"]] <- naive_complexity


save.image(paste0(PATH, data_type, "_estimation_K_in_", paste(Ks, collapse = "_"), ".RData"))

####################################################
## Likelihood plots
####################################################

load(paste0(PATH, data_type, "_estimation_K_in_", paste(Ks, collapse = "_"), ".RData"))

## First penalty
pen1 <- function(K_try, model_complexity, c){
  return(c*K_try + log(model_complexity))
}
c <- 3
K_select_1 <- which.max(df$log_likelihood - pen1(df$K_try, df$model_complexity, c))

## Second penalty
pen2 <- function(K_try, model_complexity, B){
  return((sqrt(K_try) + sqrt(2 * B * K_try + 2 * log(model_complexity)))^2)
}
B <- 0.1
K_select_2 <- which.min(df$maha_data_mean + pen2(df$K_try, df$model_complexity, B))

## Plot Likelihood
p <- ggplot(df, aes(x = K_try, y = log_likelihood))
p <- p + geom_point()
p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = "K",
              y = "Log Likelihood")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

## Plot Mahalanobis Distance
p <- ggplot(df, aes(x = K_try, y = least_squares))
p <- p + geom_point()
p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = "K",
              y = "Mahalanobis Distance")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

## Plot Likelihood penalty 1
p <- ggplot(df, aes(x = pen1(K_try, model_complexity, c), y = log_likelihood))
p <- p + geom_point()
#p <- p + geom_vline(xintercept = c*K_true + log(model_complexity[K_true]))
#p <- p + geom_vline(xintercept = c*6 + log(model_complexity[6]), linetype = 2)
p <- p + labs(x = "Penalty",
              y = "Log Likelihood")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

## Plot Likelihood penalty 2
p <- ggplot(df, aes(x = pen2(K_try, model_complexity, c), y = log_likelihood))
p <- p + geom_point()
#p <- p + geom_vline(xintercept = c*K_true + log(model_complexity[K_true]))
#p <- p + geom_vline(xintercept = c*6 + log(model_complexity[6]), linetype = 2)
p <- p + labs(x = "Penalty",
              y = "Log Likelihood")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

## Plot Mahalanobis penalty 1
p <- ggplot(df, aes(x = pen1(K_try, model_complexity, c), y = least_squares))
p <- p + geom_point()
#p <- p + geom_vline(xintercept = c*K_true + log(model_complexity[K_true]))
#p <- p + geom_vline(xintercept = c*6 + log(model_complexity[6]), linetype = 2)
p <- p + labs(x = "Penalty",
              y = "Log Likelihood")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

## Plot Mahalanobis penalty 2
p <- ggplot(df, aes(x = pen2(K_try, model_complexity, c), y = least_squares))
p <- p + geom_point()
#p <- p + geom_vline(xintercept = c*K_true + log(model_complexity[K_true]))
#p <- p + geom_vline(xintercept = c*6 + log(model_complexity[6]), linetype = 2)
p <- p + labs(x = "Penalty",
              y = "Log Likelihood")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

## Plot Likelihood penalty naive
p <- ggplot(df, aes(x = c*K_try + log(naive_complexity), y = log_likelihood))
p <- p + geom_point()
p <- p + geom_vline(xintercept = c*K_true + log(naive_complexity[K_true]))
p <- p + labs(x = "Naive Penalty",
              y = "Log Likelihood")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

## Naive complexity vs complexity
p <- ggplot(df, aes(x = log(naive_complexity), y = log(model_complexity)))
p <- p + geom_point()
p <- p + geom_vline(xintercept = log(naive_complexity[K_true]))
p <- p + labs(x = "Naive Penalty",
              y = "Tree Penalty")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Model Complexity with K
p <- ggplot(df , aes(x = K_try, y = model_complexity))
p <- p + geom_point()
p <- p + labs(x = "Number of Shifts",
              y = "Model Complexity")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Penalty with K_try
p <- ggplot(df , aes(x = K_try, y = c*K_try + log(model_complexity)))
p <- p + geom_point()
p <- p + labs(x = "Number of Shifts",
              y = "Model Complexity")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Penalty naive with K_try
p <- ggplot(df , aes(x = K_try, y = c*K_try + log(naive_complexity)))
p <- p + geom_point()
p <- p + labs(x = "Number of Shifts",
              y = "Model Complexity")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

## Plot Likelihood + penalty
p <- ggplot(df, aes(x = K_try, y = log_likelihood - c*K_try - log(model_complexity)))
p <- p + geom_point()
p <- p + geom_hline(data = df, aes(yintercept = max(log_likelihood - c*K_try - log(model_complexity))))
p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = "K",
              y = "Log Likelihood")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Estimations of alpha with K
p <- ggplot(subset(df, K_try<20) , aes(x = K_try, y = alpha_estim))
p <- p + geom_point()
#p <- p + geom_hline(yintercept = alpha)
p <- p + geom_vline(xintercept = K_true)
p <- p + geom_vline(xintercept = K_select, linetype = 2)
p <- p + labs(x = "Number of Shifts",
              y = "Estimated Selection Strength")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Estimations of gamma with K_try
p <- ggplot(df , aes(x = K_try, y = gamma_estim))
p <- p + geom_point()
#p <- p + geom_hline(yintercept = gamma)
p <- p + geom_vline(xintercept = K_true)
p <- p + geom_vline(xintercept = K_select, linetype = 2)
p <- p + labs(x = "Number of Shifts",
              y = "Estimated Root Variance")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Estimations of beta_0 with K_try
p <- ggplot(df , aes(x = K_try, y = beta_0_estim))
p <- p + geom_point()
#p <- p + geom_hline(yintercept = beta_0)
p <- p + geom_vline(xintercept = K_true)
p <- p + geom_vline(xintercept = K_select, linetype = 2)
p <- p + labs(x = "Number of Shifts",
              y = "Estimated Root Optimal Value")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Nbr of EM steps with K_try
p <- ggplot(df , aes(x = K_try, y = EM_steps))
p <- p + geom_point()
#p <- p + geom_hline(yintercept = 1000)
p <- p + geom_vline(xintercept = K_true)
p <- p + geom_vline(xintercept = K_select, linetype = 2)
p <- p + labs(x = "Number of Shifts",
              y = "Number of EM steps")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Nbr of equivalent solutions with K_try
p <- ggplot(subset(df, K_try<20) , aes(x = K_try, y = number_equivalent_solutions))
p <- p + geom_point()
p <- p + geom_vline(xintercept = K_true)
p <- p + geom_vline(xintercept = K_select, linetype = 2)
p <- p + labs(x = "Number of Shifts",
              y = "Number of Equivalent Solutions")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p
