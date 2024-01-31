#############################################
## Calibration of the penalty on simulations
#############################################

##############
## Parameters and functions
##############

require(doParallel)
require(foreach)
require(ape)
#require(glmnet) # For Lasso initialization
require(quadrupen) # For Lasso initialization
require(robustbase) # For robust fitting of alpha
library(TreeSim) # For simulation of the tree
require(ggplot2) # Plot
require(reshape2) # Plot
require(grid) # Plot

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

savedatafile_sim = "../Results/Simulation_Calibration/simulation_ou_on_tree_bayou_design_simulated_tree"
datestamp <- "2015-02-09_13-49-35"

load(paste0(savedatafile_sim, "_alpha_known-", datestamp, ".RData"))

############
## Compute Distances form true mean
############
## Tree Structure
ntaxa <- length(tree$tip.label)
times_shared <- compute_times_ca(tree)
distances_phylo <- compute_dist_phy(tree)

## How good an estimation is ?
compute_distance <- function(simest){
  ## Parameters of the simulation
  root.state <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = beta_0,
                     var.root  = simest$gamma)
  params_simu <- list(variance = simest$gamma * 2 * simest$alpha,
                      root.state = root.state,
                      shifts = simest$shifts,
                      selection.strength = simest$alpha,
                      optimal.value = beta_0)
  ## Variance - Covariance of the tips (true parameters)
  Sigma <- compute_variance_covariance.OU(times_shared = times_shared, 
                                          distances_phylo = distances_phylo,
                                          params_old = params_simu)
  Sigma_YY <- extract_variance_covariance(Sigma, what="YY")
  Sigma_YY_inv <- solve(Sigma_YY)
  ## Mean at the tips (true parameters)
  XX <- simulate_internal(phylo = tree,
                 process = process,
                 root.state = root.state, 
                 variance = 2 * simest$alpha * simest$gamma,
                 shifts = simest$shifts, 
                 selection.strength = simest$alpha, 
                 optimal.value = beta_0)
  MU <- extract_simulate_internal(XX, what="expectations", where="tips")
  ## Mean at the tips (estimated parameters)
  root.state_estim <- list(random = TRUE,
                     stationary.root = TRUE,
                     value.root = NA,
                     exp.root = simest$beta_0_estim,
                     var.root  = simest$gamma_estim)
  XX <- simulate_internal(phylo = tree,
                 process = process,
                 root.state = root.state_estim, 
                 variance = 2 * simest$alpha_estim * simest$gamma_estim,
                 shifts = simest$shifts_estim, 
                 selection.strength = simest$alpha_estim, 
                 optimal.value = simest$beta_0_estim)
  MU_estim <- extract_simulate_internal(XX, what="expectations", where="tips")
  ## Distance
  dif <- t(MU - MU_estim) %*% Sigma_YY_inv %*% (MU - MU_estim)
  return(dif)
}

distances <- lapply(simestimations_alpha_known, function(z) sapply(z, compute_distance))
distances <- do.call(c, distances)

save.image(paste0(savedatafile_sim, "_alpha_known-", datestamp, "-with_distances.RData"))

#############
## Format the data
#############
## Simple entries
simples <- c("alpha", "gamma", "K", "n", "K_try",
             "grp", "alpha_estim",
             "gamma_estim", "EM_steps", 
             "DV_estim", "CV_estim", "beta_0_estim",
             "log_likelihood", "log_likelihood.true",
             "mean_number_new_shifts",
             "number_equivalent_solutions", "complexity")

dd <- lapply(simestimations_alpha_known, function(x) t(sapply(x, function(z) matrix(unlist(z[simples]), ncol=length(simples)))))
dd <- do.call(rbind, dd)
dd <- as.data.frame(dd)
colnames(dd) <- simples

## Inject distances
dd$distances <- distances

rm(simestimations_alpha_known)
save.image(paste0(savedatafile_sim, "_alpha_known-", datestamp, "-formated.RData"))

#############
## Compute penalties fo several constants
#############
csts <- seq(0, 10, 0.1)

fun <- function(z){
  ll <- as.list(z)
  ll$c <- csts
  return(expand.grid(ll))
}
dd <- ddply(dd, colnames(dd), fun)

dd$penalties <- dd$log_likelihood - dd$c * dd$K_try - log(dd$complexity)

#############
## Select K for each model / constant
#############
data_select <- ddply(dd,
                     .(alpha, gamma, K, grp, n, c),
                     summarize,
                     K_select = K_try[which.max(penalties)],
                     distance_select = distances[which.max(penalties)]
                     )
## Compute E[||m_Y - \tilde{m}_Y(c)||^2]
mean_distance_select <- ddply(data_select,
                              .(alpha, gamma, K, grp, c),
                              summarize,
                              mean_distance_select = mean(distance_select)
)

## Compute E[||m_Y - \hat{m}_Y^K||^2] for each K_try
data_inf_K <- ddply(dd,
                    .(alpha, gamma, K, grp, K_try, c),
                    summarize,
                    distance_K = mean(distances)
)

## Compute inf E[||m_Y - \hat{m}_Y^K||^2] on K_try for each configuration
inf_mean_K <- ddply(data_inf_K,
                    .(alpha, gamma, K, grp, c),
                    summarize,
                    inf_mean_K = min(distance_K)
)

###############
## Criteria : sup {E[||m_Y - \tilde{m}_Y(c)||^2] /  inf E[||m_Y - \hat{m}_Y^K||^2]}
###############
datasum <- merge(mean_distance_select, inf_mean_K)

criteria <- inf_mean_K <- ddply(datasum,
                                .(c),
                                summarize,
                                crit = max(mean_distance_select/inf_mean_K)
)

save.image(paste0(savedatafile_sim, "_alpha_known-", datestamp, "-formated-with-criteria.RData"))

############
## Plot
############

p <- ggplot(criteria, aes(x = c, y = crit))
p <- p + geom_point()
p <- p + labs(x = "c",
              y = "F(c)")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p
# ggsave(p, file = paste0(PATH, "gamma_rmse.pdf")
#        , width = 210, height = 148, units="mm"
# )