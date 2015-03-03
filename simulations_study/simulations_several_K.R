#############################################
## Calibration of the penalty on simulations
#############################################

##############
## Parameters and functions
##############
reqpckg <- c("ape", "quadrupen", "robustbase")

require(doParallel)
require(foreach)
require(ape)
#require(glmnet) # For Lasso initialization
require(quadrupen) # For Lasso initialization
require(robustbase) # For robust fitting of alpha
library(TreeSim) # For simulation of the tree
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

## Set seed
set.seed(18051804)
savedatafile = "../Results/Simulation_Calibration/simulation_ou_on_tree_bayou_design"

## Set number of parallel cores
Ncores <- 3

## Fixed parameters for simulations
process <- "OU"
beta_0 <- 0
alpha_base <- 3
sigma_base <- 3
gamma_base <- sigma_base/(2*alpha_base)
K_base <- 9
sigma_delta_base <- 18
seg <- c("lasso", "best_single_move")

## alpha grid
alpha <- log(2)*1/c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 2, 10)

## gamma^2 grid (\gamma^2 = \sigma^2 / 2 \alpha)
gamma <- (1/(2*alpha_base))*c(0.1, 0.5, 1, 2, 5, 10, 25, 50) # enlevé : 3 (base)

## Number of shifts for simulation
K <- c(0, 1, 2)#5, 7, 8, 10, 11, 13, 20) # enlevé : 9 (base) et 50 (trop)

## Number of shifts for inference
K_try <- 0:5

## replication depth (number of replicates per )
n <- 1:200

## The combination of simulation parameters
simparams_alpha <- expand.grid(alpha, gamma_base, K_base, n, "alpha_var")
simparams_gamma <- expand.grid(alpha_base, gamma, K_base, n, "gamma_var")
simparams_K <- expand.grid(alpha_base, gamma_base, K, n, "K_var")
simparams_base <- expand.grid(alpha_base, gamma_base, K_base, n, "base")
simparams <- rbind(simparams_base,
                   simparams_alpha,
                   simparams_gamma,
                   simparams_K)
colnames(simparams) <- c("alpha", "gamma", "K", "n", "grp")

## Define date-stamp for file names
datestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
print(datestamp)


## simulation function 
# Return list of parameters + list of shifts + data at tips
datasetsim <- function(alpha, gamma, K, n, grp) {
  shifts <- sample_shifts(tree, sigma_delta_base, K)
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = beta_0,
                     var.root = gamma)
  params <-  list(variance = 2*alpha*gamma,
                  root.state = root.state,
                  shifts = shifts,
                  selection.strength = alpha, 
                  optimal.value = beta_0)
  XX <- simulate(phylo = tree,
                 process = process,
                 root.state = root.state, 
                 variance = 2*alpha*gamma,
                 shifts = shifts, 
                 selection.strength = alpha, 
                 optimal.value = beta_0)
  sim <- list(alpha = alpha,
              gamma = gamma, 
              K = K,
              n = n,
              K_try = K_try,
              grp = grp,
              shifts = shifts,
              Y_data = extract.simulate(XX, what="states", where="tips"),
              Z_data = extract.simulate(XX, what = "states", where = "nodes"),
              m_Y_data = extract.simulate(XX, what="expectations", where="tips"))
  sim$log_likelihood.true <- log_likelihood.OU(sim$Y_data, tree, params)
  return(sim)
}

sample_shifts <- function(tree, sigma_delta, K){
  if (K == 0) return(NULL)
  shifts_edges <- sample_shifts_edges(tree, K)
  shifts_values <- sample_shifts_values(sigma_delta, K)
  shifts <- list(edges = shifts_edges, 
                 values = shifts_values, 
                 relativeTimes = rep(0, K))
  return(shifts)
}

sample_shifts_values <- function(sigma_delta, K){
  return(rnorm(K, mean=0, sd=sqrt(sigma_delta)))
}

## estimation function
# Entrée : sortie de datasetsim + K_try number of shifts for EM
# Sortie : liste avec l'entrée, et les parametres estimés.
estimationfunction_alpha_known <- function(X, K_t) {
  ## If an estimation fails, catch error with "try" and try again
  results_estim_EM <- estimateEM(phylo=tree, 
                                 Y_data=X$Y_data, 
                                 tol = list(variance=10^(-4), 
                                            value.root=10^(-4), 
                                            exp.root=10^(-4), 
                                            var.root=10^(-4), 
                                            selection.strength=10^(-3)), 
                                 process = process, 
                                 method.variance = "simple", 
                                 method.init = "lasso",
                                 method.init.alpha = "estimation",
                                 Nbr_It_Max = 1000, 
                                 nbr_of_shifts = K_t, 
                                 alpha_known = TRUE, ##
                                 known.selection.strength = X$alpha,
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
                                 methods.segmentation = seg)
  extract.edges <- function(x) {
    z <- unlist(lapply(x, function(y) y$shifts$edges))
    if (!is.null(z)) z <- matrix(z, nrow = K_t)      
    return(z)
  }
  selected.edges <- extract.edges(results_estim_EM$params_history)
  params <- results_estim_EM$params
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
  X$number_new_shifts <- results_estim_EM$number_new_shifts
  X$mean_number_new_shifts <- mean(results_estim_EM$number_new_shifts)
  X$number_equivalent_solutions <- results_estim_EM$number_equivalent_solutions
  X$K_try <- K_t
  X$complexity <- choose(2*ntaxa-2-K_t, K_t)
  return(X)
}

estimations_several_K_alpha_known <- function(X){
  XX <- sapply(K_try, function(z) estimationfunction_alpha_known(X, z))
  return(XX)
}

##################################
## Simulations
##################################
savedatafile_sim <- paste0(savedatafile, "_simulated_tree")

## simulate tree
ntaxa <- 64
lambda <- 0.1
tree <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = lambda, mu = 0, age = 1, mrca = TRUE)[[1]]
plot(tree)

## Sequencial simulations (for reproductability)
simlist <- foreach(i = 1:nrow(simparams), .packages = reqpckg[1]) %do% {
  alpha <- simparams[i, "alpha"]
  gamma <- simparams[i, "gamma"]
  K <- simparams[i, "K"]
  n <- simparams[i, "n"]
  grp <- simparams[i, "grp"]
  sim <- datasetsim(alpha, gamma, K, n, grp)
  return(sim)
}

names(simlist) <- apply(simparams, 1, paste0, collapse = "_")

############
## Estimations (alpha known)
############

## Register parallel backend for computing
cl <- makeCluster(Ncores)
registerDoParallel(cl)


n.range <- 1:200
simulations2keep <- lapply(simlist, function(x) { x$n %in% n.range }, simplify = TRUE)
simlist <- simlist[simulations2keep]

## Parallelized estimations
time_alpha_known <- system.time(simestimations_alpha_known <- foreach(i = simlist, .packages = reqpckg) %dopar%
{
  estimations_several_K_alpha_known(i)
}
)

# Stop the cluster (parallel)
stopCluster(cl)
save.image(paste0(savedatafile_sim, "_alpha_known-", datestamp, ".RData"))
