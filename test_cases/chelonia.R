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

###############################################################################
## Import DataSet with habitat clustering
###############################################################################
## Import data set form package geiger
library(geiger)
data(chelonia)

tree <- chelonia$phy
data <- chelonia$dat
data <- data[match(tree$tip.label, names(data))]

plot(tree, show.tip.label = FALSE)

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
clusters <- as.vector(chel_df$habitat)
Nbr_ParSols <- extract.parsimonyNumber(parsimonyNumber(tree, clusters))
ParSols <- extract.enumerate_parsimony(enumerate_parsimony(tree, clusters))
Nbr_shifts_OUwie <- extract.parsimonyCost(parsimonyCost(tree, clusters), 227)

## Plot clustering
colors_habitat <- chel_df$habitat
levels(colors_habitat) <- 1:length(levels(colors_habitat))

colors_regimes <- as.factor(ParSols[1, tree$edge[,2]])
levels(colors_regimes) <- 1:length(levels(colors_regimes))

plot(tree, show.tip.label = FALSE, edge.color = colors_regimes, x.lim = 300)
tiplabels(frame = "circle", col = colors_habitat, pch = 16, adj = 4)
legend("topright", legend = levels(chel_df$habitat), col = levels(colors_habitat), pch = 16)

###############################################################################
## EM
###############################################################################
# K <- 16
# name <- paste0("_K=", K)
# time <- system.time(
#   results_estim_EM <- estimateEM(phylo = tree, 
#                                  Y_data = data, 
#                                  tol = list(variance=10^(-4), 
#                                             value.root=10^(-4), 
#                                             exp.root=10^(-4), 
#                                             var.root=10^(-4), 
#                                             selection.strength=10^(-3)), 
#                                  process = "OU", 
#                                  method.variance = "simple", 
#                                  method.init = "lasso",
#                                  method.init.alpha = "estimation",
#                                  Nbr_It_Max = 1000, 
#                                  nbr_of_shifts = K, 
#                                  alpha_known = FALSE,
#                                  min_params=list(variance = 10^(-4), 
#                                                  value.root = -10^(4), 
#                                                  exp.root = -10^(4), 
#                                                  var.root = 10^(-4),
#                                                  selection.strength = 10^(-4)),
#                                  max_params=list(variance = 10^(4), 
#                                                  value.root = 10^(4), 
#                                                  exp.root = 10^(4), 
#                                                  var.root = 10^(4),
#                                                  selection.strength = 10^(4)),
#                                  methods.segmentation = c("same_shifts",
#                                                           "best_single_move",
#                                                           "lasso"))
# )
# 
# ## Save History and process
# history <- list_to_table.history(results_estim_EM$params_history)
# write.csv2(history, paste0(PATH, "boite_noire", name, ".csv"))
# plot.history.OU.stationnary(params_history = results_estim_EM$params_history,
#                             tree = tree,
#                             Y_data_ref = data,
#                             PATH = PATH,
#                             name = paste0("history_plot", name))
# plot.process.actual(Y.state = data,
#                     Z.state = results_estim_EM$ReconstructedNodesStates,
#                     phylo = tree, 
#                     paramsEstimate = results_estim_EM$params)

## Load EM Solutions (computed with script Likelihood_plot.R)
Ks <- 1:63
data_type <- "chelonia"
load(paste0("../Results/Likelihood_Plot/", data_type, "_estimation_K_in_", paste(Ks, collapse = "_"), ".RData"))

## Select one K
c <- 3
K_select <- which.max(df$log_likelihood - c*df$K - log(df$model_complexity))

## Equivalent Solutions
times_shared <- compute_times_ca(tree)
t_tree <-  min(node.depth.edgelength(tree)[1:ntaxa])
T_tree <- incidence.matrix(tree)
ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                      selection.strength = estimations[[K_select]]$alpha_estim,
                      times_shared = times_shared)
T_tree_ac <- T_tree * ac_tree

eq_shifts_edges_K_select <- equivalent_shifts_edges(tree, 
                                            estimations[[K_select]]$shifts$edges)
eq_shifts_values_K_select <- equivalent_shifts_values(tree,
                                  shifts = estimations[[K_select]]$shifts,
                                  beta_0 = estimations[[K_select]]$beta_0_estim,
                                  eq_shifts_edges = eq_shifts_edges_K_select,
                                  T_tree_ac = T_tree_ac)

plot_equivalent_shifts.actual(tree, eq_shifts_edges_K_select, eq_shifts_values_K_select, use.edge.length = FALSE, adj = 0)

## Solution compared with habitat
sol <- 1

colors_habitat <- chel_df$habitat
levels(colors_habitat) <- rainbow(length(levels(colors_habitat)), start = 0, v = 0.5)

regimes <- allocate_regimes_from_shifts(tree, eq_shifts_edges_K_select[, sol])
edges_regimes <- regimes[tree$edge[,2]]
tips_regimes <- regimes[1:ntaxa]

colors_regimes <- as.factor(edges_regimes)
levels(colors_regimes) <- 1:length(levels(colors_regimes))
tips_regimes <- as.factor(tips_regimes)
levels(tips_regimes) <- 1:length(levels(tips_regimes))

plot(tree, show.tip.label = FALSE, edge.color = as.vector(colors_regimes), x.lim = 370, y.lim = 220)
segments(220, 1:length(data), 220 + 10*(data), 1:length(data), col = as.vector(tips_regimes))
segments(212, 1:length(data), 212 + 5, col = as.vector(colors_habitat), lwd = 3)
segments(220, -3, 220 + 10, -3, lwd = 2)
text(250, -3, "1 cm", cex = 0.8)
#points(rep(215, length(tree$tip.label)), 1:length(tree$tip.label), pch=19, col=colors_habitat, cex=0.5)
#tiplabels(pch = 19, cex = abs(data)/mean(abs(data)), col = colors_habitat, adj = 4)
#tiplabels(frame = "circle", col = colors_habitat, pch = 16, adj = 4)
legend("topright", legend = levels(chel_df$habitat), col = levels(colors_habitat), pch = 16)

## Equivalent Solutions K_true
times_shared <- compute_times_ca(tree)
t_tree <-  min(node.depth.edgelength(tree)[1:ntaxa])
T_tree <- incidence.matrix(tree)
ac_tree <- incidence_matrix_actualization_factors(tree = tree, 
                                                  selection.strength = estimations[[K_true]]$alpha_estim,
                                                  times_shared = times_shared)
T_tree_ac <- T_tree * ac_tree

eq_shifts_edges_K_true <- equivalent_shifts_edges(tree, 
                                                    estimations[[K_true]]$shifts$edges)
eq_shifts_values_K_true <- equivalent_shifts_values(tree,
                                                      shifts = estimations[[K_true]]$shifts,
                                                      beta_0 = estimations[[K_true]]$beta_0_estim,
                                                      eq_shifts_edges = eq_shifts_edges_K_true,
                                                      T_tree_ac = T_tree_ac)

## Save Image
save.image(paste0(PATH, data_type, "_estimation_K_in_", paste(Ks, collapse = "_"), ".RData"))
# rm(results_estim_EM, history, times_shared, t_tree, T_tree, ac_tree, T_tree_ac, eq_shifts_edges_sim, eq_shifts_values_sim)

###############################################################################
## OUwie
###############################################################################
library(OUwie)

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

EM_select <- list(Nbr_shifts = length(estimations[[K_select]]$shifts$edges),
                  Nbr_regimes = length(estimations[[K_select]]$shifts$edges) + 1,
                  lnL = estimations[[K_select]]$log_likelihood,
                  alpha = estimations[[K_select]]$alpha_estim,
                  half_life = log(2)/estimations[[K_select]]$alpha_estim,
                  sigma = 2 * estimations[[K_select]]$alpha_estim * estimations[[K_select]]$gamma_estim)

EM_habitat <- list(Nbr_shifts = length(estimations[[K_true]]$shifts$edges),
                  Nbr_regimes = length(estimations[[K_true]]$shifts$edges) + 1,
                  lnL = estimations[[K_true]]$log_likelihood,
                  alpha = estimations[[K_true]]$alpha_estim,
                  half_life = log(2)/estimations[[K_true]]$alpha_estim,
                  sigma = 2 * estimations[[K_true]]$alpha_estim * estimations[[K_select]]$gamma_estim)

OU_habitat <- list(Nbr_shifts = Nbr_shifts_OUwie,
                   Nbr_regimes = 4,
                   lnL = res_OUwie$loglik,
                   alpha = res_OUwie$solution["alpha", 1],
                   half_life = log(2)/res_OUwie$solution["alpha", 1],
                   sigma = res_OUwie$solution["sigma.sq", 1])

Summary <- cbind(OU_habitat, EM_habitat, EM_select)

write.csv2(Summary, paste0(PATH, "summary", ".RData"))
save.image(paste0(PATH, "summary", ".RData"))

EM <- list(Nbr_shifts = length(results_estim_EM$params$shifts$edges),
           Nbr_regimes = length(results_estim_EM$params$shifts$edges) + 1,
           lnL = attr(results_estim_EM$params, "log_likelihood")[1],
           alpha = results_estim_EM$params$selection.strength,
           half_life = log(2)/results_estim_EM$params$selection.strength,
           sigma = results_estim_EM$params$variance)
           