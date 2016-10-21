rm(list=ls())

library(ape)
library(plyr)
# library(quadrupen)
# library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(TreeSim)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

library(testthat)

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
sourceCpp("src/upward_downward.cpp")

## ntaxa fixed, number of traits change
set.seed(17920902)
ntaxa = 500
tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                                 age = 1, mrca = TRUE)[[1]]
tree <- reorder(tree, order = "postorder")
times_shared <- compute_times_ca(tree)
distances_phylo <- compute_dist_phy(tree)
T_tree <- incidence.matrix(tree)
U_tree <- incidence.matrix.full(tree)
h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
root_edge_length <- 0
if (!is.null(tree$root.edge)) root_edge_length <- tree$root.edge
F_moments <- compute_fixed_moments(times_shared + root_edge_length, ntaxa)

micro_time_old <- vector(length = 15)
micro_time_updown <- vector(length = 15)
for(p in 1:15){
  variance <- matrix(0.8, p, p) + diag(0.2, p, p)
  independent <- FALSE
  root.state <- list(random = FALSE,
                     value.root = rep(1, p),
                     exp.root = NA,
                     var.root = NA)
  shifts = list(edges = c(18, 32, 45, 109, 254, 398),
                values=cbind(c(4, -10, 3, 12, 32, -5),
                             c(-5, 5, 0, -1.1, 32.89, 16),
                             c(4, -10, 3, 12, 32, -5),
                             c(4, -10, 3, 12, 32, -5),
                             c(4, -10, 3, 12, 32, -5),
                             c(4, -10, 3, 12, 32, -5)),
                relativeTimes = 0)
  shifts$values <- apply(shifts$values, 1, function(z) rep(z, length.out = p))
  shifts$values <- matrix(shifts$values, nrow = p)
  paramsSimu <- list(variance = variance,
                     shifts = shifts,
                     root.state = root.state)
  attr(paramsSimu, "p_dim") <- p
  
  X1 <- simulate_internal(tree,
                 p = p,
                 root.state = root.state,
                 process = "BM",
                 variance = variance,
                 shifts = shifts)
  traits <- extract_simulate_internal(X1, where = "tips", what = "state")
  
  micro_time_old[p] <- print(microbenchmark(wrapper_E_step(phylo = tree,
                                                        times_shared = times_shared,
                                                        distances_phylo = distances_phylo,
                                                        process = "BM",
                                                        params_old = paramsSimu,
                                                        masque_data = masque_data,
                                                        F_moments = F_moments,
                                                        independent = independent,
                                                        Y_data_vec_known = as.vector(traits[!miss]),
                                                        miss = miss,
                                                        Y_data = traits,
                                                        U_tree = U_tree,
                                                        compute_E = compute_E.simple.nomissing.BM),
                                         times = 100))$median
  
  micro_time_updown[p] <- print(microbenchmark(wrapper_E_step(phylo = tree,
                                                           process = "BM",
                                                           params_old = paramsSimu,
                                                           independent = independent,
                                                           Y_data = traits,
                                                           compute_E = compute_E.upward_downward),
                                            times = 100))$median
}

library(ggplot2)
library(reshape2)
data <- data.frame("Factorized" = micro_time_old,
                   "Upward Downward" = micro_time_updown,
                   "Dimension" = 1:p)
data <- melt(data, id.vars = "Dimension", 
             variable.name = "Method", value.name = "time")
ggplot(data, aes(x = Dimension, y = time, color = Method)) + geom_point()

## p fixed, ntaxa change
set.seed(17920902)
p <- 4
variance <- matrix(0.8, p, p) + diag(0.2, p, p)
independent <- FALSE
root.state <- list(random = FALSE,
                   value.root = rep(1, p),
                   exp.root = NA,
                   var.root = NA)
shifts = list(edges = c(18, 32, 45, 95),
              values=cbind(c(4, -10, 3, 12),
                           c(-5, 5, 0, -1.1),
                           c(4, -10, 3, 12),
                           c(4, -10, 3, 12),
                           c(4, -10, 3, 12),
                           c(4, -10, 3, 12)),
              relativeTimes = 0)
shifts$values <- apply(shifts$values, 1, function(z) rep(z, length.out = p))
shifts$values <- matrix(shifts$values, nrow = p)
paramsSimu <- list(variance = variance,
                   shifts = shifts,
                   root.state = root.state)
attr(paramsSimu, "p_dim") <- p

ntaxa_list <- c(50, 100, 300, 500, 700, 1000)
micro_time_old <- vector(length = length(ntaxa_list))
micro_time_updown <- vector(length = length(ntaxa_list))
for(i in 1:length(ntaxa_list)){
  
  ntaxa <- ntaxa_list[i]
  tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                                   age = 1, mrca = TRUE)[[1]]
  tree <- reorder(tree, order = "postorder")
  times_shared <- compute_times_ca(tree)
  distances_phylo <- compute_dist_phy(tree)
  T_tree <- incidence.matrix(tree)
  U_tree <- incidence.matrix.full(tree)
  h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
  root_edge_length <- 0
  if (!is.null(tree$root.edge)) root_edge_length <- tree$root.edge
  F_moments <- compute_fixed_moments(times_shared + root_edge_length, ntaxa)
  
  X1 <- simulate_internal(tree,
                 p = p,
                 root.state = root.state,
                 process = "BM",
                 variance = variance,
                 shifts = shifts)
  traits <- extract_simulate_internal(X1, where = "tips", what = "state")
  
  micro_time_old[i] <- print(microbenchmark(wrapper_E_step(phylo = tree,
                                                           times_shared = times_shared,
                                                           distances_phylo = distances_phylo,
                                                           process = "BM",
                                                           params_old = paramsSimu,
                                                           masque_data = masque_data,
                                                           F_moments = F_moments,
                                                           independent = independent,
                                                           Y_data_vec_known = as.vector(traits[!miss]),
                                                           miss = miss,
                                                           Y_data = traits,
                                                           U_tree = U_tree,
                                                           compute_E = compute_E.simple.nomissing.BM),
                                            times = 100))$median
  
  micro_time_updown[i] <- print(microbenchmark(wrapper_E_step(phylo = tree,
                                                              process = "BM",
                                                              params_old = paramsSimu,
                                                              independent = independent,
                                                              Y_data = traits,
                                                              compute_E = compute_E.upward_downward),
                                               times = 100))$median
}

library(ggplot2)
library(reshape2)
data <- data.frame("Factorized" = micro_time_old,
                   "Upward Downward" = micro_time_updown,
                   "Dimension" = ntaxa_list)
data <- melt(data, id.vars = "Dimension", 
             variable.name = "Method", value.name = "time")
ggplot(data, aes(x = Dimension, y = time, color = Method)) + geom_point()

## Missing Data
set.seed(17920902)
p <- 4
variance <- matrix(0.8, p, p) + diag(0.2, p, p)
independent <- FALSE
root.state <- list(random = FALSE,
                   value.root = rep(1, p),
                   exp.root = NA,
                   var.root = NA)
shifts = list(edges = c(18, 32, 45, 95),
              values=cbind(c(4, -10, 3, 12),
                           c(-5, 5, 0, -1.1),
                           c(4, -10, 3, 12),
                           c(4, -10, 3, 12),
                           c(4, -10, 3, 12),
                           c(4, -10, 3, 12)),
              relativeTimes = 0)
shifts$values <- apply(shifts$values, 1, function(z) rep(z, length.out = p))
shifts$values <- matrix(shifts$values, nrow = p)
paramsSimu <- list(variance = variance,
                   shifts = shifts,
                   root.state = root.state)
attr(paramsSimu, "p_dim") <- p

ntaxa <- 500
tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = 0.1, mu = 0, 
                                 age = 1, mrca = TRUE)[[1]]
tree <- reorder(tree, order = "postorder")
times_shared <- compute_times_ca(tree)
distances_phylo <- compute_dist_phy(tree)
T_tree <- incidence.matrix(tree)
U_tree <- incidence.matrix.full(tree)
h_tree <- max(diag(as.matrix(times_shared))[1:ntaxa])
root_edge_length <- 0
if (!is.null(tree$root.edge)) root_edge_length <- tree$root.edge
F_moments <- compute_fixed_moments(times_shared + root_edge_length, ntaxa)

X1 <- simulate_internal(tree,
               p = p,
               root.state = root.state,
               process = "BM",
               variance = variance,
               shifts = shifts)

NA_per_list <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
micro_time_old <- vector(length = length(ntaxa_list))
micro_time_updown <- vector(length = length(ntaxa_list))
for(i in 1:length(NA_per_list)){
  
  NA_per <- NA_per_list[i]
  traits <- extract_simulate_internal(X1, where = "tips", what = "state")
  nMiss <- floor(ntaxa * p * NA_per)
  miss <- sample(1:(p * ntaxa), nMiss, replace = FALSE)
  chars <- (miss - 1) %% p + 1
  tips <- (miss - 1) %/% p + 1
  for (i in 1:nMiss){
    traits[chars[i], tips[i]] <- NA
  }
  miss <- as.vector(is.na(traits))
  masque_data <- rep(FALSE, (ntaxa + tree$Nnode) * p)
  masque_data[1:(p * ntaxa)] <- !miss
  
  micro_time_old[i] <- print(microbenchmark(wrapper_E_step(phylo = tree,
                                                           times_shared = times_shared,
                                                           distances_phylo = distances_phylo,
                                                           process = "BM",
                                                           params_old = paramsSimu,
                                                           masque_data = masque_data,
                                                           F_moments = F_moments,
                                                           independent = independent,
                                                           Y_data_vec_known = as.vector(traits[!miss]),
                                                           miss = miss,
                                                           Y_data = traits,
                                                           U_tree = U_tree,
                                                           compute_E = compute_E.simple),
                                            times = 100))$median
  
  micro_time_updown[i] <- print(microbenchmark(wrapper_E_step(phylo = tree,
                                                              process = "BM",
                                                              params_old = paramsSimu,
                                                              independent = independent,
                                                              Y_data = traits,
                                                              compute_E = compute_E.upward_downward),
                                               times = 100))$median
}

library(ggplot2)
library(reshape2)
data <- data.frame("Factorized" = micro_time_old,
                   "Upward Downward" = micro_time_updown,
                   "Dimension" = ntaxa_list)
data <- melt(data, id.vars = "Dimension", 
             variable.name = "Method", value.name = "time")
ggplot(data, aes(x = Dimension, y = time, color = Method)) + geom_point()
