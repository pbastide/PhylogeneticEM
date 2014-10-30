rm(list = ls())

WD <- "/home/bastide/Dropbox/These/Code/Phylogenetic-EM"
#WD <- "/Users/paulb/Dropbox/These/Code/Phylogenetic-EM" # (Mac)
setwd(WD)

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
source("simulate.R")
source("estimateEM.R")
source("init_EM.R")
source("E_step.R")
source("M_step.R")
source("shutoff.R")
source("generic_functions.R")
source("shifts_manipulations.R")
source("plot_functions.R")

## Set seed
#set.seed(1121983)
#set.seed(21031989)
set.seed(18051804)
savedatafile = "../Results/Simulation_Estimation_Bayou_Design_new_seg/simulation_ou_on_tree_bayou_design"

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
gamma <- (1/(2*alpha_base))*c(0.1, 1, 2, 3, 5, 10, 25, 50) # enlevé : 0.5 (base)

## snr (signal to noise ratio, computed as \sigma_delta^2 / \gamma^2)
#snr <- c(0.2, 0.5, 1, 2, 5)

## Number of shifts
K <- c(1, 5, 7, 8, 10, 11, 13, 20, 50) # enlevé : 9 (base)

## replication depth (number of replicates per )
n <- 1:2

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
  XX <- simulate(phylo = tree,
                 process = process,
                 root.state = root.state, 
                 variance = 2*alpha*gamma,
                 shifts = shifts, 
                 selection.strength = alpha, 
                 optimal.value = beta_0)
  return(list(alpha = alpha,
              gamma = gamma, 
              K = K,
              n = n,
              grp = grp,
              shifts = shifts,
              Y_data = extract.simulate(XX, what="states", where="tips"),
              Z_data = extract.simulate(XX, what = "states", where = "nodes")))
}

sample_shifts <- function(tree, sigma_delta, K){
  shifts_edges <- sample_shifts_edges(tree, K)
  shifts_values <- sample_shifts_values(sigma_delta, K)
  shifts <- list(edges = shifts_edges, 
                 values = shifts_values, 
                 relativeTimes = rep(0, K))
  return(shifts)
}

sample_shifts_edges <- function(tree, K){
  Nbranch <- nrow(tree$edge)
  ntaxa <- length(tree$tip.label)
  if (K > ntaxa) stop("A parcimonious repartition of the shifts cannot have more shifts than tips !")
  ## Generate shifts so that they are parcimonious
  Tr <- incidence.matrix(tree)
  Nbr_grps <- 0
  It <- 0
  while ((Nbr_grps < K+1 && It < 500)) {
    ## Generate K branches
    edges <- sample_edges_intervals(tree, K)
    ## Generate Groups
    groups <- rep(0, ntaxa); names(groups) <- tree$tip.label
    for (i in order(edges)) { # Do it in order
      ed <- edges[i]
      groups[tree$tip.label[Tr[,ed]]] <- i
    }
    Nbr_grps <- length(unique(groups))
    It <- It + 1
  }
  if (It==500) stop("I could not find a parcimonious repartition of the shifts after 500 iterations. Please consider taking a smaller number of shifts.")
  return(edges)
}

sample_edges_intervals <- function(tree, K){
  pas <- 1/K * 0:K
  node_heights <- node.depth.edgelength(tree)
  groups <- split(1:length(node_heights), findInterval(node_heights, pas))
  sh <- NULL; p <- 1
  for (k in 1:K) {
    grp <- groups[[paste(k)]]
    if (is.null(grp)){
      p <- p + 1
    } else {
      if (p <= length(grp)){
        if (length(grp) == 1) {
          sh <- c(sh, grp)
        } else {
          sh <- c(sh, sample(grp, p))
        }
        p <- 1
      } else {
        sh <- c(sh, grp)
        p <- p - length(grp) + 1
      }
    }
  }
  sh <- sapply(sh, function(x) which(tree$edge[,1]==x)[1])
  if (length(sh) < K) {
    p <- K-length(sh)
    sh <- c(sh, sample(tree$edge[-sh], p))
  }
  return(sh)
}

sample_shifts_values <- function(sigma_delta, K){
  return(rnorm(K, mean=0, sd=sqrt(sigma_delta)))
}

## estimation function
# Entrée : sortie de datasetsim
# Sortie : liste avec l'entrée, et les parametres estimés.
estimationfunction <- function(X) {
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
                                 nbr_of_shifts = X$K, 
                                 alpha_known = FALSE,
                                 min_params=list(variance = 10^(-3), 
                                                 value.root = -10^(3), 
                                                 exp.root = -10^(3), 
                                                 var.root = 10^(-3),
                                                 selection.strength = 10^(-3)),
                                 max_params=list(variance = 10^(3), 
                                                 value.root = 10^(3), 
                                                 exp.root = 10^(3), 
                                                 var.root = 10^(3),
                                                 selection.strength = 10^(3)),
                                 methods.segmentation = seg)
  extract.edges <- function(x) {
    z <- unlist(lapply(x, function(y) y$shifts$edges))
    z <- matrix(z, nrow = X$K)      
    return(z)
  }
  selected.edges <- extract.edges(results_estim_EM$params_history)
  params <- results_estim_EM$params
  X$alpha_estim = params$selection.strength
  X$gamma_estim = params$root.state$var.root
  X$shifts_estim = params$shifts
  X$EM_steps <- attr(results_estim_EM, "Nbr_It")
  X$DV_estim <- attr(results_estim_EM, "Divergence")
  X$CV_estim <- (attr(results_estim_EM, "Nbr_It") != 500) && !X$DV_estim
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
  return(X)
}

## estimation function
# Entrée : sortie de datasetsim
# Sortie : liste avec l'entrée, et les parametres estimés.
estimationfunction_alpha_known <- function(X) {
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
                                 nbr_of_shifts = X$K, 
                                 alpha_known = TRUE, ##
                                 known.selection.strength = X$alpha,
                                 min_params=list(variance = 10^(-3), 
                                                 value.root = -10^(3), 
                                                 exp.root = -10^(3), 
                                                 var.root = 10^(-3),
                                                 selection.strength = 10^(-3)),
                                 max_params=list(variance = 10^(3), 
                                                 value.root = 10^(3), 
                                                 exp.root = 10^(3), 
                                                 var.root = 10^(3),
                                                 selection.strength = 10^(3)),
                                 methods.segmentation = seg)
  extract.edges <- function(x) {
    z <- unlist(lapply(x, function(y) y$shifts$edges))
    z <- matrix(z, nrow = X$K)      
    return(z)
  }
  selected.edges <- extract.edges(results_estim_EM$params_history)
  params <- results_estim_EM$params
  X$alpha_estim = params$selection.strength
  X$gamma_estim = params$root.state$var.root
  X$shifts_estim = params$shifts
  X$EM_steps <- attr(results_estim_EM, "Nbr_It")
  X$DV_estim <- attr(results_estim_EM, "Divergence")
  X$CV_estim <- (attr(results_estim_EM, "Nbr_It") != 500) && !X$DV_estim
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
  return(X)
}

##################################
## Mammal Tree
##################################
savedatafile_mammal <- paste0(savedatafile, "_mammal_tree")
## import tree
load("../data/Several_Trees.RData")
tree <- Cetacea_Autocorrelated
tree <- reorder.phylo(tree, order="cladewise")

## normalize tree height
scale.tree <- function(phylo){
  if (!is.ultrametric(phylo)) stop("The tree is not ultrametric")
  ntaxa <- length(phylo$tip.label)
  height <- min(node.depth.edgelength(phylo)[1:ntaxa]) - .Machine$double.eps^0.5# take the min so that any error is above 1
  phylo$edge.length <- phylo$edge.length/height
  return(phylo)
}

tree <- scale.tree(tree)

## Sequencial simulations (for reproductibility)
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
## Alpha known

## Register parallel backend for computing
cl <- makeCluster(Ncores)
registerDoParallel(cl)

## Parallelized estimations
time_alpha_known <- system.time(simestimations_alpha_known <- foreach(i = simlist, .packages = reqpckg) %dopar%
{
  estimationfunction_alpha_known(i)
}
)

# Stop the cluster (parallel)
stopCluster(cl)
save.image(paste0(savedatafile_mammal, "_alpha_known-", datestamp, ".RData"))

############
## Alpha unknown

## Register parallel backend for computing
cl <- makeCluster(Ncores)
registerDoParallel(cl)

## Parallelized estimations
time <- system.time(simestimations <- foreach(i = simlist, .packages = reqpckg) %dopar%
{
  estimationfunction(i)
}
)

# Stop the cluster (parallel)
stopCluster(cl)
save.image(paste0(savedatafile_mammal, "-", datestamp, ".RData"))

rm(savedatafile_mammal, tree, scale.tree, simlist, cl, time_alpha_known, simestimations_alpha_known, time, simestimations)

##################################
## Simulated Tree
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
## Alpha known

## Register parallel backend for computing
cl <- makeCluster(Ncores)
registerDoParallel(cl)

## Parallelized estimations
time_alpha_known <- system.time(simestimations_alpha_known <- foreach(i = simlist, .packages = reqpckg) %dopar%
{
  estimationfunction_alpha_known(i)
}
)

# Stop the cluster (parallel)
stopCluster(cl)
save.image(paste0(savedatafile_sim, "_alpha_known-", datestamp, ".RData"))

############
## Alpha unknown

## Register parallel backend for computing
cl <- makeCluster(Ncores)
registerDoParallel(cl)

## Parallelized estimations
time <- system.time(simestimations <- foreach(i = simlist, .packages = reqpckg) %dopar%
{
  estimationfunction(i)
}
)

# Stop the cluster (parallel)
stopCluster(cl)
save.image(paste0(savedatafile_sim, "-", datestamp, ".RData"))
