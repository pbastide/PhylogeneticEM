##########
## Test Case : Chelonian carapace evolution.
## Used in :
## - Uyeda 2014
## - Eastman 2011 (link to data : https://github.com/eastman/auteur/tree/master/auteur/data)
## -> Jaffe 2011 (original data)

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

Ncores <- 3

###################################################
## Import or generate data
###################################################

## Import data set form package geiger
# library(geiger)
# data(chelonia)
# 
# tree <- chelonia$phy
# data <- chelonia$dat
# data <- data[match(tree$tip.label, names(data))]
# 
# plot(tree, show.tip.label = FALSE)

## Random data
set.seed(20141210)
ntaxa <- 64
lambda <- 0.1
tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = lambda, mu = 0, age = 1, mrca = TRUE)[[1]]
plot(tree, show.tip.label = FALSE)

beta_0 <- 0
gamma <- 1
alpha <- 3

K_true <- 3
shifts <- sample_shifts(tree, 18, K_true)

Ks <- 1:21
data_type <- paste0("random_ntaxa=", ntaxa, "lambda=", lambda)

XX <- simulate(phylo = tree,
               process = "OU",
               root.state = list(random = TRUE,
                                 stationary.root = TRUE,
                                 value.root = NA,
                                 exp.root = beta_0,
                                 var.root = gamma), 
               variance = 2*alpha*gamma,
               shifts = shifts, 
               selection.strength = alpha, 
               optimal.value = beta_0)

data = extract.simulate(XX, what="states", where="tips")

###############################################################################
## EM for several values of K
##############################################################################

estimationfunction <- function(K) {
  ## If an estimation fails, catch error with "try" and try again
  time <- system.time(
    results_estim_EM <- estimateEM(phylo = tree, 
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
                                   nbr_of_shifts = K, 
                                   alpha_known = FALSE,
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
                                   methods.segmentation = c("same_shifts",
                                                            "best_single_move",
                                                            "lasso"))
  )
  extract.edges <- function(x) {
    z <- unlist(lapply(x, function(y) y$shifts$edges))
    z <- matrix(z, nrow = K)      
    return(z)
  }
  selected.edges <- extract.edges(results_estim_EM$params_history)
  params <- results_estim_EM$params
  X <- NULL
  X$K <- K
  X$time <- time
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
  edge.quality <- unlist(lapply(params$shifts$edges, compute.quality))
  names(edge.quality) <- params$shifts$edges
  X$edge.quality <- edge.quality
  X$start <- results_estim_EM$params_history["0"]
  X$beta_0_estim <- params$root.state$exp.root
  X$log_likelihood <- attr(params, "log_likelihood")[1]
  X$number_new_shifts <- results_estim_EM$number_new_shifts
  X$mean_number_new_shifts <- mean(results_estim_EM$number_new_shifts)
  X$number_equivalent_solutions <- results_estim_EM$number_equivalent_solutions
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
## Exploitation for likelihood plot
####################################################

## Arange Data
dd <- do.call(rbind, estimations)
df <- apply(dd[ , colnames(dd) %in% c("alpha", "gamma", "K", "n",
                                              "grp", "alpha_estim",
                                              "gamma_estim", "EM_steps", 
                                              "DV_estim", "CV_estim",
                                              "difficulty", "beta_0_estim",
                                              "log_likelihood",
                                              "mean_number_new_shifts"
                                              ,"number_equivalent_solutions")],
            2, unlist)
df <- as.data.frame(df)

## Model Complexity
model_complexity <- sapply(Ks, function(z) extract.partitionsNumber(partitionsNumber(tree, z + 1)))
df[["model_complexity"]] <- model_complexity

## Plot Likelihood
p <- ggplot(df, aes(x = model_complexity, y = log_likelihood))
p <- p + geom_point()
p <- p + geom_vline(x_intercept = K_true)
p <- p + labs(x = "Model Complexity",
              y = "Log Likelihood")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Model Complexity with K
p <- ggplot(df , aes(x = K, y = model_complexity))
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

# Estimations of alpha with K
p <- ggplot(df , aes(x = K, y = alpha_estim))
p <- p + geom_point()
p <- p + geom_hline(yintercept = alpha)
p <- p + labs(x = "Number of Shifts",
              y = "Estimated Selection Strength")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Estimations of gamma with K
p <- ggplot(df , aes(x = K, y = gamma_estim))
p <- p + geom_point()
p <- p + geom_hline(yintercept = gamma)
p <- p + labs(x = "Number of Shifts",
              y = "Estimated Selection Strength")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Estimations of beta_0 with K
p <- ggplot(df , aes(x = K, y = beta_0_estim))
p <- p + geom_point()
p <- p + geom_hline(yintercept = beta_0)
p <- p + labs(x = "Number of Shifts",
              y = "Estimated Selection Strength")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Nbr of EM steps with K
p <- ggplot(df , aes(x = K, y = EM_steps))
p <- p + geom_point()
p <- p + geom_hline(yintercept = 1000)
p <- p + labs(x = "Number of Shifts",
              y = "Estimated Selection Strength")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

# Nbr of equivalent solutions with K
p <- ggplot(df , aes(x = K, y = number_equivalent_solutions))
p <- p + geom_point()
p <- p + labs(x = "Number of Shifts",
              y = "Estimated Selection Strength")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

save.image(paste0(PATH, data_type, "_estimation_K_in_", paste(Ks, collapse = "_"), ".RData"))
