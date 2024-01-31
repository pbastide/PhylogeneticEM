##############
## Parameters and functions
##############

#require(doParallel)
require(foreach)
require(ape)
#require(glmnet) # For Lasso initialization
#require(quadrupen) # For Lasso initialization
#require(robustbase) # For robust fitting of alpha
library(TreeSim) # For simulation of the tree
#require(ggplot2) # Plot
#require(reshape2) # Plot
#require(grid) # Plot
library(Matrix)

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

savedatafile = "../Results/Simulations_Several_K/several_K_simlist_uncertainties"

## Fixed parameters for simulations
process <- "OU"
beta_0 <- 0
alpha_base <- 3
sigma_base <- 3
gamma_base <- 0.5
K_base <- 5 # Change par rapport à Uyeda (au lieu de 9)
# sigma_delta_base <- 18
m1 <- 4; m2 <- - m1
s1 <- 1; s2 <- s1
seg <- c("lasso", "best_single_move")

# ## alpha grid
# alpha <- log(2)*1/c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 2, 10)
# 
# ## gamma^2 grid (\gamma^2 = \sigma^2 / 2 \alpha)
# gamma <- c(0.05, 0.1, 1, 2, 3, 5, 10, 25) # enlevé : 0.5 (base)
# 
# ## Number of shifts for simulation
# K <- c(0:4, 8, 11, 16) # enlevé : 5 (base)

## Noise for (gaussian) uncertainty
s_noise_base <- 0
s_noise <- c(0.01, 0.02, 0.05, 0.1, 0.5, 1, 5, 10)

## Number of taxa
ntaxa <- c(64, 128, 256)

# ## Number of shifts for inference (depends on the size of the tree)
# K_try <- 0:sqrt(nta)

## replication depth (number of replicates per )
n <- 1:200

## The combination of simulation parameters
# simparams_alpha <- expand.grid(alpha, gamma_base, K_base, ntaxa, n, "alpha_var")
# simparams_gamma <- expand.grid(alpha_base, gamma, K_base, ntaxa, n, "gamma_var")
# simparams_K <- expand.grid(alpha_base, gamma_base, K, ntaxa, n, "K_var")
simparams_noise <- expand.grid(alpha_base, gamma_base, K_base, s_noise, ntaxa, n, "s_noise_var")
simparams_base <- expand.grid(alpha_base, gamma_base, K_base, s_noise_base, ntaxa, n, "base")
simparams <- rbind(simparams_base,
#                    simparams_alpha,
#                    simparams_gamma,
#                    simparams_K,
                   simparams_noise)
colnames(simparams) <- c("alpha", "gamma", "K", "s_noise", "ntaxa", "n", "grp")

## Generation of trees
set.seed(03032015)
trees <- vector(mode = "list")
times_shared <- vector(mode = "list")
distances_phylo <- vector(mode = "list") ;
h_tree <- vector(mode = "list") ;
subtree.list <- NULL; T_tree <- NULL; K_try <- NULL
lambda <- 0.1
for (nta in ntaxa){
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
  K_try[[paste0(nta)]] <- 0:floor(sqrt(nta))
}


## Define date-stamp for file names
datestamp_data <- format(Sys.time(), "%Y-%m-%d")

## simulation function 
# Return list of parameters + list of shifts + data at tips
datasetsim <- function(alpha, gamma, K, s_noise, ntaxa, n, grp) {
  tree <- trees[[paste0(ntaxa)]]
  shifts <- sample_shifts_GMM(tree, m1, m2, s1, s2, K)
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
  attr(params, "ntaxa")  <- ntaxa
  attr(params, "p_dim")  <- 1
  XX <- simulate_internal(phylo = tree,
                 process = process,
                 root.state = root.state, 
                 variance = 2*alpha*gamma,
                 shifts = shifts, 
                 selection.strength = alpha, 
                 optimal.value = beta_0)
  ## Add noise
  Y_data <- extract_simulate_internal(XX, what="states", where="tips")
  Y_noise <- t(rnorm(n = ntaxa, sd = sqrt(s_noise)))
  ## Retuen results
  sim <- list(alpha = alpha,
              gamma = gamma, 
              K = K,
              s_noise = s_noise,
              n = n,
              ntaxa = ntaxa,
              grp = grp,
              shifts = shifts,
              Y_data = Y_data + Y_noise,
              Y_noise = Y_noise,
              Z_data = extract_simulate_internal(XX, what = "states", where = "nodes"),
              m_Y_data = extract_simulate_internal(XX, what="expectations", where="tips"))
  ## Compute true likelihood and difficulty of the problem
  # Moments
  Sigma <- compute_variance_covariance.scOU(times_shared = times_shared[[paste0(ntaxa)]], 
                                            distances_phylo = distances_phylo[[paste0(ntaxa)]],
                                            params_old = params)
  Sigma_YY <- extract_variance_covariance(Sigma, what="YY")
  Sigma_YY_chol <- chol(Sigma_YY)
  Sigma_YY_chol_inv <- backsolve(Sigma_YY_chol, diag(ncol(Sigma_YY_chol)))
  Sigma_YY_inv <- Sigma_YY_chol_inv %*% t(Sigma_YY_chol_inv)
  # Difficulty
  mu_0 <- (sum(Sigma_YY_inv))^(-1) * sum(Sigma_YY_inv%*%t(sim$m_Y_data))
  sim$difficulty <- as.vector((sim$m_Y_data - mu_0) %*% Sigma_YY_inv %*% t(sim$m_Y_data - mu_0))
  # Log likelihood
  sim$log_likelihood.true <- compute_log_likelihood.simple(phylo = tree,
                                                           Y_data = as.vector(sim$Y_data),
                                                           sim = XX,
                                                           Sigma = Sigma,
                                                           Sigma_YY_chol_inv = Sigma_YY_chol_inv)
  
  return(sim)
}

##################################
## Simulations
##################################
## Set seed
set.seed(18051804)

## Sequencial simulations (for reproductability)
simlist <- foreach(i = 1:nrow(simparams)) %do% {
  sim <- datasetsim(alpha = simparams[i, "alpha"],
                    gamma = simparams[i, "gamma"],
                    K = simparams[i, "K"],
                    s_noise = simparams[i, "s_noise"],
                    ntaxa = simparams[i, "ntaxa"],
                    n = simparams[i, "n"],
                    grp = simparams[i, "grp"])
  sim$it <- i
  return(sim)
}

names(simlist) <- apply(simparams, 1, paste0, collapse = "_")

save.image(paste0(savedatafile, "_", datestamp_data, ".RData"))