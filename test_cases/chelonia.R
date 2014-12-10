##########
## Test Case : Chelonian carapace evolution.
## Used in :
## - Uyeda 2014
## - Eastman 2011 (link to data : https://github.com/eastman/auteur/tree/master/auteur/data)
## -> Jaffe 2011 (original data)

rm(list=ls())

PATH <- "../Results/Chelonia/"

library(ape)
library(glmnet) # For Lasso initialization
library(quadrupen)
#library(nlme) # For second derivative computation
#library(combinat) # For alpha prior robust estimation
library(robustbase) # For robust fitting of alpha
library(ggplot2) # Plot
library(reshape2) # Plot
library(grid) # Plot
library(TreeSim)
library(plyr)
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

## Import data set form package geiger
library(geiger)
data(chelonia)

tree <- chelonia$phy
data <- chelonia$dat
data <- data[match(tree$tip.label, names(data))]

plot(tree, show.tip.label = FALSE)

###############################################################################
## EM
##############################################################################
K <- 16
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

## Save History and process
name <- paste0("_K=", K)
history <- list_to_table.history(results_estim_EM$params_history)
write.csv2(history, paste0(PATH, "boite_noire", name, ".csv"))
plot.history.OU.stationnary(params_history = results_estim_EM$params_history,
                            tree = tree,
                            Y_data_ref = data,
                            PATH = PATH,
                            name = paste0("history_plot", name))
plot.process.actual(Y.state = data,
                    Z.state = results_estim_EM$ReconstructedNodesStates,
                    phylo = tree, 
                    paramsEstimate = results_estim_EM$params)

## Equivalent Solutions
times_shared <- compute_times_ca(tree)
t_tree <-  min(node.depth.edgelength(tree)[1:ntaxa])
T_tree <- incidence.matrix(tree)
ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                                                  selection.strength = results_estim_EM$params$selection.strength,
                                                  times_shared = times_shared)
T_tree_ac <- T_tree * ac_tree

eq_shifts_edges_sim <- equivalent_shifts_edges(tree, results_estim_EM$params$shifts$edges)
eq_shifts_values_sim <- equivalent_shifts_values(tree,
                                                 shifts = results_estim_EM$params$shifts,
                                                 beta_0 = results_estim_EM$params$optimal.value,
                                                 eq_shifts_edges = eq_shifts_edges_sim,
                                                 T_tree_ac = T_tree_ac)

plot_equivalent_shifts.actual(tree, eq_shifts_edges_sim, eq_shifts_values_sim, use.edge.length = FALSE, adj = 0)

## Save Image
save.image(paste0(PATH, "estimation", name, ".RData"))
rm(results_estim_EM, history, times_shared, t_tree, T_tree, ac_tree, T_tree_ac, eq_shifts_edges_sim, eq_shifts_values_sim)

###############################################################################
## OUwie
###############################################################################
library(OUwie)

## Import Data from Jaffe
vv <- read.csv("../data/Chelonia_habitats_vector.txt", header = FALSE)
chel_df <- NULL
chel_df$species <- as.character(vv[c(2:44, 178:222, 362:406, 546:590, 730:774, 914:916), ])
chel_df$accession_number <-  as.character(vv[c(46:88, 224:268, 408:452, 592:636, 776:820, 918:920), ])
tmp <- vv[c(90:132, 270:314, 454:498, 638:682, 822:866, 922:924), ]
chel_df$length <- as.numeric(levels(tmp))[tmp]
chel_df$habitat <-  factor(vv[c(134:176, 316:360, 500:544, 684:728, 868:912, 926:928), ])
chel_df <- as.data.frame(chel_df)
# Correspondance
chel_df <- chel_df[match(names(data), chel_df$species), ]

## Clustering Habitat
clusters <- chel_df$habitat
Nbr_ParSols <- extract.parsimonyNumber(parsimonyNumber(tree, clusters))
ParSols <- extract.enumerate_parsimony(enumerate_parsimony(tree, clusters))
Nbr_shifts_OUwie <- extract.parsimonyCost(parsimonyCost(tree, clusters), 227)

## format data for OUwie
treeWIE <- tree
treeWIE$node.label <- ParSols[1,(length(tree$tip.label)+1):dim(ParSols)[2]]
dataWIE <- chel_df[ ,c("species", "habitat", "length")]
dataWIE$length <- log(dataWIE$length)

res_OUwie <- OUwie(phy = treeWIE,
                   data = dataWIE,
                   model = "OUM")

save.image(paste0(PATH, "estimation_OUwie", ".RData"))
rm(vv, treeWIE, dataWIE, res_OUwie)

###############################################################################
## Summary
###############################################################################
load(paste0(PATH, "estimation", name, ".RData"))
load(paste0(PATH, "estimation_OUwie", ".RData"))

EM <- list(Nbr_shifts = length(results_estim_EM$params$shifts$edges),
           Nbr_regimes = length(results_estim_EM$params$shifts$edges) + 1,
           lnL = attr(results_estim_EM$params, "log_likelihood")[1],
           alpha = results_estim_EM$params$selection.strength,
           half_life = log(2)/results_estim_EM$params$selection.strength,
           sigma = results_estim_EM$params$variance)

OU_habitat <- list(Nbr_shifts = Nbr_shifts_OUwie,
                   Nbr_regimes = 4,
                   lnL = res_OUwie$loglik,
                   alpha = res_OUwie$solution["alpha", 1],
                   half_life = log(2)/res_OUwie$solution["alpha", 1],
                   sigma = res_OUwie$solution["sigma.sq", 1])

Summary <- cbind(EM, OU_habitat)

write.csv2(Summary, paste0(PATH, "summary", ".RData"))
save.image(paste0(PATH, "summary", ".RData"))
           