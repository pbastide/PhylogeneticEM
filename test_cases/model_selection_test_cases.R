##########
## Test Case : Varying K in the inference.

rm(list=ls())

PATH <- "../Results/Model_Selection_Test_Cases/"

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
source("R/model_selection.R")

Ncores <- 3

###################################################
## Import or generate data
###################################################

# ## Import data set form package geiger
# library(geiger)
# data(chelonia)
# 
# tree <- chelonia$phy
# data <- chelonia$dat
# data <- data[match(tree$tip.label, names(data))]
# 
# plot(tree, show.tip.label = FALSE)
# 
# K_max <- 30
# data_type <- "chelonia"
# ntaxa <- length(data)
# K_true <- 16
# 
# alpha <- 0.038 # Pick a "reasonable" alpha
# # alpha <- 0.03204742 # K = 6 Baraud et al
# # alpha <- 0.04470039 # K = 8 Birgé Massart
# data_type <- paste0(data_type, "_alpha=", alpha)

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
# data_type <- paste0("random_ntaxa=", ntaxa, "_K_true=", K_true)
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

## Easy data small tree
set.seed(20141211)
ntaxa <- 64
lambda <- 0.1
tree <- TreeSim::sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = lambda, mu = 0, age = 1, mrca = TRUE)[[1]]
plot(tree, show.tip.label = FALSE); edgelabels();

beta_0 <- 0
gamma <- 0.1
alpha <- 3
alpha_try <- 100 # 0.1 1 2 4

K_true <- 5
shifts <- list(edges = c(61, 96, 87, 11, 30),
               values = c(2, -2, 2, -2, 5),
               relativeTimes = rep(0, K_true))

K_max <- 30
data_type <- paste0("easy_ntaxa=", ntaxa, "_K_true=", K_true, "_alpha_try=",alpha_try)

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

plot.data.process.actual(Y.state = extract.simulate(XX, what="states", where="tips"),
                    #Z.state = extract.simulate(XX, what="states", where="nodes"),
                    phylo = tree, 
                    params = list(shifts = shifts, optimal.value = beta_0),
                    automatic_colors = TRUE)

# ## Easy data big tree
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
# K_max <- 35
# data_type <- paste0("easy_ntaxa=", ntaxa, "_K_true=", K_true)
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

estimations <- estimateEM_several_K.OUsr(phylo = tree, 
                                         Y_data = data, 
                                         K_max = K_max,
                                         alpha_known = TRUE,
                                         alpha = alpha_try)

save.image(paste0(PATH, data_type, "_estimation_K_max=", K_max, ".RData"))

###############################################################################
## Birge massart
###############################################################################
load(paste0(PATH, data_type, "_estimation_K_max=", K_max, ".RData"))

df <- estimations$summary
B <- 0.1

## Plot - Least Squares against penalty
p <- ggplot(df, aes(x = penalty_BirgeMassart_shape1(K_try, complexity, B), y = -least_squares))
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

## Slope Heuristic
K_max <- 35
df_sub <- subset(df, K_try <= K_max)
pen_shape <- penalty_BirgeMassart_shape1(df_sub$K_try, df_sub$complexity, B)
data_capushe <- data.frame(names = df_sub$K_try, 
                           pen_shape = pen_shape,
                           complexity = df_sub$complexity,
                           contrast = df_sub$least_squares)
# Slope Heuristic
DDSE_results <- DDSE(data_capushe)
DDSE_results@model
DDSE_results@interval$interval
plot(DDSE_results, newwindow = FALSE)
# Dimension Jump
Djump_results <- Djump(data_capushe)
Djump_results
Djump_results@ModelHat$Kopt
plot(Djump_results, newwindow = FALSE)

## Plot Least Squares + penalty (comparisons)
pen_shape <- penalty_BirgeMassart_shape1(df$K_try, df$complexity, B)
crits <- data.frame(K_try = df$K_try, 
                    "0.2" = df$least_squares + 0.2 * pen_shape,
                    Slope = df$least_squares + 2*DDSE_results@interval$interval["max"]*pen_shape,
                    Djump = df$least_squares + Djump_results@ModelHat$Kopt*pen_shape)
crits_plot <- melt(crits, id.vars = "K_try", value.name = "criteria", variable.name = "calibration")
levels(crits_plot$calibration) <- c("0.2", "Slope Heuristic", "Dimension Jump")
crits_min <- data.frame(K_try = c(which.min(crits$X0.2), which.min(crits$Slope), which.min(crits$Djump)) - 1,
                        calibration = as.factor(c("0.2", "Slope Heuristic", "Dimension Jump")),
                        criteria = c(min(crits$X0.2), min(crits$Slope), min(crits$Djump)))

p <- ggplot(crits_plot, aes(x = K_try, y = criteria, color = calibration))
p <- p + geom_point()
p <- p + geom_point(data = crits_min, aes(x = K_try, y = criteria, size = 5))
p <- p + labs(x = "K",
              y = "Penalized Least Squares")
p <- p + scale_size(name = "", labels = "Min")
p <- p + scale_x_continuous(breaks = c(0, 5, 10, 13, 20, 30))
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

## Process Plots
# True process
plot.process.actual(Y.state = extract.simulate(XX, what="states", where="tips"),
                    Z.state = extract.simulate(XX, what="states", where="nodes"),
                    phylo = tree, 
                    paramsEstimate = list(shifts = shifts, optimal.value = beta_0))
# DDSE
Kt <- as.integer(DDSE_results@model) + 1
plot.process.actual(Y.state = extract.simulate(XX, what="states", where="tips"),
                    Z.state = estimations[[Kt]]$Zhat,
                    phylo = tree, 
                    paramsEstimate = list(shifts = estimations[[Kt]]$shifts, 
                                          optimal.value = estimations[[Kt]]$beta_0_estim))
# Djump
Kt <- as.integer(Djump_results@model) + 1
plot.process.actual(Y.state = extract.simulate(XX, what="states", where="tips"),
                    Z.state = estimations[[Kt]]$Zhat,
                    phylo = tree, 
                    paramsEstimate = list(shifts = estimations[[Kt]]$shifts, 
                                          optimal.value = estimations[[Kt]]$beta_0_estim))
# True K
Kt <- 5 + 1
plot.process.actual(Y.state = extract.simulate(XX, what="states", where="tips"),
                    Z.state = estimations[[Kt]]$Zhat,
                    phylo = tree, 
                    paramsEstimate = list(shifts = estimations[[Kt]]$shifts, 
                                          optimal.value = estimations[[Kt]]$beta_0_estim))

###############################################################################
## Baraud Giraud Huet
###############################################################################
df <- estimations$summary
D_max <- compute_K_max(ntaxa, kappa = 0.9)
pen_ll <- 1/2 * penalty_BaraudGiraudHuet_likelihood(df$K_try, df$complexity, ntaxa, C = 1.1)

## Plot Likelihood
p <- ggplot(df, aes(x = K_try, y = log_likelihood))
p <- p + geom_point()
# p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = "K",
              y = "Log Likelihood")
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

## Plot Likelihood + penalty
crit <- data.frame(K_try = df$K_try, 
                   crit_ll = -df$log_likelihood + pen_ll)
crit_min <- subset(crit, crit_ll == min(crit_ll))

p <- ggplot(crit, aes(x = K_try, y = crit_ll))
p <- p + geom_point()
p <- p + geom_point(data = crit_min, aes(x = K_try, y = crit_ll, size = 5))
# p <- p + geom_vline(xintercept = K_true)
p <- p + labs(x = "K",
              y = "Log Likelihood")
p <- p + scale_size(name = "", labels = "Min")
p <- p + scale_x_continuous(breaks = c(0, 5, 10, 20, 30))
p <- p + theme_bw()
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12)
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
p

