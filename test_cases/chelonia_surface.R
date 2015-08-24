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
## SURFACE
###############################################################################
library(surface)

## Import data set form package geiger
library(bayou)
data(chelonia)

tree <- chelonia$phy
dat <- chelonia$dat
dat <- dat[match(tree$tip.label, names(dat))]

tree <- nameNodes(tree)
dat <- as.data.frame(dat)
colnames(dat) <- "Carapace_length"
olist <- convertTreeData(tree, dat)
otree <- olist[[1]]
odata <- olist[[2]]

## Forward
time_fwd <- system.time(fwd <- surfaceForward(otree, odata,
                                              aic_threshold = 0, exclude = 0,
                                              verbose = FALSE, plotaic = FALSE))
k <- length(fwd)

fsum <- surfaceSummary(fwd)
names(fsum)
fsum$aics
fwd[[k]]

## Backward
time_bwd <- system.time(bwd <- surfaceBackward(otree, odata,
                                               starting_model = fwd[[k]],
                                               aic_threshold = 0,
                                               only_best = TRUE,
                                               verbose = FALSE, plotaic = FALSE))

bsum<-surfaceSummary(bwd)
kk<-length(bwd)

bsum$alpha
bsum$sigma_squared
bsum$theta

surfaceTreePlot(tree, bwd[[kk]], labelshifts = T)

surfaceColors <- function (tree, hansenfit) {
  fit <- hansenfit$fit[[1]]
  otree <- as(fit, "data.frame")
  otree <- data.frame(otree,
                      shifts = rep(NA, length(otree$nodes)))
  otree$shifts[match(names(hansenfit$savedshifts), otree$nodes)] <- 1:length(hansenfit$savedshifts)
  ntip <- (dim(otree)[1] + 1)/2
  nnode <- ntip - 1
  otree2 <- otree[match(c(tree$tip.label, tree$node.label), 
                        otree$labels), ]
  otree2 <- otree2[tree$edge[, 2], ]
  xx <- summary(factor(hansenfit$savedshifts))
  cols <- c("black", rainbow(length(xx) - 1))
  edgecols <- cols[as.numeric(factor(otree2[, 5]))]
  return(edgecols)
}

colorsSurface <- surfaceColors(tree, bwd[[kk]])

mapping <- function(tree, colors){
  nedges <- length(tree$edge.length)
  ## Produce the map
  maps <- vector("list", nedges)
  cols <- tree$edge.length
  names(cols) <- colors
  for (i in 1:nedges){
    maps[[i]] <- cols[i]
  }
  ## Produce the colors
  colors <- unique(colors)
  names(colors) <- colors
  return(list(maps = maps,
              colors = colors))
}

surfaceMapping <- mapping(tree, colorsSurface)

tree$maps <- surfaceMapping$maps
colors <- surfaceMapping$colors
plotSimmap(tree, colors = colors, fsize = 0, ftype = "off")

Q<-matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3)
rownames(Q)<-colnames(Q)<-letters[1:3]
treeb<-sim.history(pbtree(n=100,scale=1),Q)
cols<-setNames(c("blue","red","green"),letters[1:3])
# plot the mapping
plotSimmap(treeb,cols,ftype="i",fsize=0.7)

Surface_fwd <- list(Nbr_shifts = k,
                    Nbr_regimes = unname(fsum$n_regimes[2]),
                    lnL = fsum$lnls[k],
                    MlnL = NaN,
                    alpha = unname(fsum$alpha),
                    half_life = unname(fsum$phylhalflife),
                    sigma = unname(fsum$sigma_squared),
                    gamma = unname(fsum$sigma_squared / (2 * fsum$alpha)),
                    time = unname(time_fwd)[3])

Surface_bwd <- list(Nbr_shifts = k,
               Nbr_regimes = unname(bsum$n_regimes[2]),
               lnL = bsum$lnls[kk],
               MlnL = NaN,
               alpha = unname(bsum$alpha),
               half_life = unname(bsum$phylhalflife),
               sigma = unname(bsum$sigma_squared),
               gamma = unname(bsum$sigma_squared / (2 * bsum$alpha)),
               time = unname(time_fwd + time_bwd)[3])

save.image(paste0(PATH, "chelonia_surface", ".RData"))
save(Surface_fwd, Surface_bwd, colorsSurface, surfaceMapping, surfacebetas, savedshiftsSurface, file = paste0(PATH, "chelonia_surface_summary.RData"))

####################

surfaceBetas <- function (tree, hansenfit) {
  fit <- hansenfit$fit[[1]]
  otree <- as(fit, "data.frame")
  otree <- data.frame(otree,
                      shifts = rep(NA, length(otree$nodes)))
  otree$shifts[match(names(hansenfit$savedshifts), otree$nodes)] <- 1:length(hansenfit$savedshifts)
  ntip <- (dim(otree)[1] + 1)/2
  nnode <- ntip - 1
  otree2 <- otree[match(c(tree$tip.label, tree$node.label), 
                        otree$labels), ]
  #otree2 <- otree2[tree$edge[, 2], ]
  xx <- summary(factor(hansenfit$savedshifts))
  betas <- hansenfit$fit$Carapace_length@theta[[1]]
  betas <- betas[as.numeric(factor(otree2[, 5]))]
  return(betas)
}

surfacebetas <- surfaceBetas(tree, bwd[[kk]])
savedshiftsSurface <- bwd[[kk]]$savedshifts



