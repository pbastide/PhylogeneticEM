rm(list=ls())

WD <- "/Users/paulb/Dropbox/These/Code" # Dossier de travail (Mac)
#WD <- "/home/bastide/Dropbox/These/Code" # Dossier de travail (Ubuntu)
setwd(WD)


library(ape)
library(plyr)
library(microbenchmark)

source("Phylogenetic-EM/generic_functions.R")
source("Phylogenetic-EM/parcimonyNumber.R")

tree <- read.tree(text="(((T,T),C),C);")
plot(tree); tiplabels(); nodelabels()

clusters=c(1,2,3,3)

## Finds the correct number of parsimonious allocations
extract.parcimonyNumber(parcimonyNumber(tree,clusters))
microbenchmark(extract.parcimonyNumber(parcimonyNumber(tree,clusters)), times = 1000L)

extract.enumerate_parsimony(enumerate_parsimony(tree,clusters))
microbenchmark(extract.parcimonyNumber(parcimonyNumber(tree,clusters)), times = 1000L)

# Reconstruction (1,2,3,3,3,3,3) has also two shifts, and is missing !

# Problem with the algorithm at node 7 : (1,1,0) from the informations given by the
# tips, excluding solution that starts with state 3 (two shifts bellow it, but no 
# shift before)
