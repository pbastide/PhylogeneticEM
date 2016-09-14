# Package: pmc
# Version: 0.0-8
# Date: 2012-02-22
# Title: Phylogenetic Monte Carlo
# Author: Carl Boettiger <cboettig@gmail.com>
#   Maintainer: Carl Boettiger <cboettig@gmail.com>
#   Description: Monte Carlo based model choice for applied phylogenetics
# of continuous traits
# URL: https://github.com/cboettig/pmc
# BugReports: https://github.com/cboettig/pmc/issues
# License: CC0
# Authors@R: person("Carl", "Boettiger", role=c("aut", "cre"),
#                   email="cboettig@gmail.com")
# Type: Package
# Imports: ape, geiger, snowfall, TreeSim, methods, ouch, reshape2
# Depends: ggplot2, ouch
# Enhances: ouch, geiger, ape
# Collate: 'format_data.R' 'pmc.R' 'powercurve.R' 'wrappers_geiger.R'
# 'wrappers_laser.R' 'wrappers_ouch.R' 'wrappers_treepar.R'
# 'convert.R' 'treepalette.R' 'generics.R' 'data_documentation.R'
# 'maticce.R'
# Packaged: 2013-03-04 01:49:53 UTC; cboettig
# NeedsCompilation: no
# Repository: CRAN
# Date/Publication: 2013-03-04 07:19:19


# toggles between ouch and ape format trees
# @param ot a phylogenetic tree in ouch or ape format
# @param regimes if given in ape format, are appended to
# phylo format as phy$regimes.  If the ouchtree is a fitted
# hansen object, regimes will automatically be imported from it
# unless other regime choice is given here.
# @param safe mode for going from ape to ouch.  After converting, this
# writes to a temporary nexus file and reads the tree back in, because
# phylo format does not have a unique specification for a unique tree,
# resulting in all kinds of silly problems when developers haven't been careful.
# Defaults to true, and will clean up after itself.
# @return a phylogenetic tree in the opposite format
convert.pmc <- function(ot, regimes=NULL, safe=TRUE){
  if(!is.null(regimes))
    safe<-FALSE # cannot write to nexus if regimes are required
	if(is(ot, "ouchtree")){

		n <- ot@nnodes
		# find the tips vs internals
		anc <- as.integer(ot@ancestors[!is.na(ot@ancestors)])
		internal <- sort(anc)[seq(2,n-1, by=2)]
		tmp <- integer(n)
		tmp[internal] = 1
		tips <- which(tmp == 0)

		root <- which(is.na(ot@ancestors))
		internal <- internal[internal!=root] #remove root form internal list
		
		new_ids <- integer(n)
		new_ids[tips] <- 1:length(tips)
		new_ids[root] <- length(tips)+1
		new_ids[internal] <- (length(tips)+2):n
		
		new_ancestor <- new_ids[as.integer(ot@ancestors)]

		edge <- matrix(NA, n-1, 2)
		edge[,1] <- new_ancestor[!is.na(new_ancestor)]
		edge[,2] <- new_ids[!is.na(new_ancestor)]

		anc <- as.integer(ot@ancestors[!is.na(ot@ancestors)])
		lengths <- ot@times[!is.na(ot@ancestors)] - ot@times[anc]
	
		labels <- ot@nodelabels[tips]

		tree <- list(edge=edge, Nnode = (n-1)/2, tip.label = labels, edge.length= lengths )
		class(tree) <- "phylo"

    if(safe){
      write.nexus(file="tmp3612411.nex", tree)
      tree <- read.nexus("tmp3612411.nex")
      unlink("tmp.nex")
    }

		if (is(ot, "hansentree")) {
			regimes <- ot@regimes[[1]][-1]
			tree$regimes <- regimes
		}

    if(!is.null(regimes)){
      tree$regimes <- regimes[-1]
    }

	} else	if (is(ot, "phylo")){ 
			tree <- ape2ouch(ot)
		}
#	plot(tree)
	tree
}



