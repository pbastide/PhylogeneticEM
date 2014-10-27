# {Shifts Manipulations}
# Copyright (C) {2014}  {SR, MM, PB}
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

###############################################################################
## Here is a set of functions to handle both list-like and vector like 
## representations of the shifts (through the incedence matrix of the tree), and to
## make the correspondances between regimes and shifts.
## Dependencies : generic_functions.R
###############################################################################

#####################################################
## Incidences matrices and list-vector manipulations
#####################################################

##
#' @title Incidence matrix of a tree.
#'
#' @description
#' \code{incidence.matrix} computes the incidence matrix T of a tree : for a 
#' lineage i and a branch b, T[i,b]=1 if b is in the lineage i, and 0 
#' otherwise.
#'
#' @details
#' This function uses the general up tree recursion function \code{recursionUp}
#' with the initialization function \code{init.incidence.matrix} and update
#' function \code{update.incidence.matrix}.
#'
#' @param phylo Input tree.
#' 
#' @return Matrix of incidence.
#'
#'17/06/14 - Initial release
##
incidence.matrix <- function(phylo){
  ## Reorder and keep track
  phy <- reorder(phylo, order = "postorder")
  cor <- correspondanceEdges(edges=1:nrow(phy$edge), from=phylo, to=phy)
  ## Init and recurence
  T <- init.incidence.matrix(phy)
  T <- recursionUp(phy, T, update.incidence.matrix)
  ## Take for each node its parenting branch
  daughters <- phy$edge[,2]
  T <- T[daughters,]
  ## Return to original tree
  T <- T[cor, ]
  return(t(T))
}

##
#' @title Initialization for incidence matrix
#'
#' @description
#' \code{init.incidence.matrix} initialize the matrix updated in
#' \code{update.incidence.matrix} for the computation of the incidence matrix
#' in \code{incidence.matrix}.
#'
#' @details
#' The initialized matrix has ntaxa column and nNodes rows. Each node
#' represent its parental branch. A row corresponding to a tip i is initialized
#' to a vector of zeros, with only entry i equal to one. (Branch ending at 
#' tip i is only in the i^th lineage)
#'
#' @param phy Input tree.
#' 
#' @return Matrix with nNodes rows and ntaxa column.
#'
#'17/06/14 - Initial release
##
init.incidence.matrix <- function(phy){
  ntaxa <- length(phy$tip.label)
  T <- matrix(NA, nrow = 1 + nrow(phy$edge), ncol = ntaxa)
  for (i in 1:ntaxa){
    T[i,] <- 1:ntaxa == i
  }
  return(T)
}

##
#' @title Update function for incidence matrix
#'
#' @description
#' \code{update.incidence.matrix} updates the matrix initialized in 
#' \code{init.incidence.matrix} for the computation of the incidence matrix
#' in \code{incidence.matrix}.
#'
#' @details
#' A node belongs to all the lineages of its daughters.
#'
#' @param daughtersParams : rows of updated matrix corresponding to the
#' daughters of the current node.
#' 
#' @return Vector of length ntaxa, indicating to which lineages the branch
#' above the current node belongs to.
#'
#'17/06/14 - Initial release
##
update.incidence.matrix <- function(daughtersParams, ...){
  inc <- colSums(daughtersParams)
  return(inc>0)
}

##
#' @title Incidence matrix of a tree.
#'
#' @description
#' \code{incidence.matrix.full} computes the incidence matrix T of a tree : for a 
#' node i and a branch b, T[i,b]=1 if b is in the lineage i, and 0 
#' otherwise.
#'
#' @details
#' This function uses the general up tree recursion function \code{recursionUp}
#' with the initialization function \code{init.incidence.matrix.full} and update
#' function \code{update.incidence.matrix.full}.
#'
#' @param phylo Input tree.
#' 
#' @return Matrix of incidence, size ntaxa + nNodes
#'
#'06/10/14 - Initial release
##
incidence.matrix.full <- function(phylo){
  ## Reorder and keep track
  phy <- reorder(phylo, order = "postorder")
  cor <- correspondanceEdges(edges=1:nrow(phy$edge), from=phylo, to=phy)
  ## Init and recurence
  U <- init.incidence.matrix.full(phy)
  U <- recursionUp(phy, U, updateUp = update.incidence.matrix.full)
  ## Take for each node its parenting branch
  daughters <- phy$edge[,2]
  U <- U[daughters,]
  ## Return to original tree
  U <- U[cor, ]
  return(t(U))
}

##
#' @title Initialization for incidence matrix (full tree)
#'
#' @description
#' \code{init.incidence.matrix.full} initialize the matrix updated in
#' \code{update.incidence.matrix.full} for the computation of the incidence
#'  matrix of the full tree in \code{incidence.matrix.full}.
#'
#' @details
#' The initialized matrix is squared of size ntaxa + nNodes. Each node
#' represent its parental branch. A row corresponding to a tip i is initialized
#' to a vector of zeros, with only entry i equal to one. (Branch ending at 
#' tip i is only in the i^th lineage)
#'
#' @param phy Input tree.
#' 
#' @return Matrix of size ntaxa + nNodes.
#'
#'06/10/14 - Initial release
##
init.incidence.matrix.full <- function(phy){
  ntaxa <- length(phy$tip.label)
  nNodes <- phy$Nnode
  U <- matrix(NA, nrow = nNodes, ncol = ntaxa + nNodes)
  T_1 <- matrix(0, nrow = ntaxa, ncol = ntaxa + nNodes)
  U <- rbind(T_1, U)
  diag(U) <- rep(1, ntaxa + nNodes)
  return(U)
}

##
#' @title Update function for incidence matrix
#'
#' @description
#' \code{update.incidence.matrix.full} updates the matrix initialized in 
#' \code{init.incidence.matrix.full} for the computation of the incidence matrix
#' in \code{incidence.matrix.full}.
#'
#' @details
#' A node belongs to all the lineages of its daughters.
#'
#' @param daughtersParams : rows of updated matrix corresponding to the
#' daughters of the current node.
#' 
#' @return Vector of length ntaxa + nNodes, indicating to which lineages the 
#' branch above the current node belongs to.
#'
#'06/10/14 - Initial release
##
update.incidence.matrix.full <- function(daughtersParams, parent, ...){
  inc <- colSums(daughtersParams)
  inc[parent] <- 1
  return(inc > 0)
}

##
#' @title Compute the vector of shifts.
#'
#' @description
#' \code{shifts.list_to_vector} takes the list description of the shifts 
#' to give the vectorial representation of the shifts : the b th element of 
#' the vector has the value of the shift occuring on that branch b.
#'
#' @details
#'
#' @param phy Input tree.
#' @param shifts : list description of the shifts : shifts$edges, shifts$values
#' 
#' @return Vector of length nbranch.
#' 
#' @seealso \code{shifts.vector_to_list}
#' 
#'17/06/14 - Initial release
##
shifts.list_to_vector <- function(phy, shifts){
  delta <- rep(0, nrow(phy$edge))
  delta[shifts$edges] <- shifts$values
  return(delta)
}

##
#' @title Compute the list of shifts.
#'
#' @description
#' \code{shifts.vector_to_list} takes the vectorial description of the shifts 
#' to create the list description of the shifts.
#'
#' @details
#'
#' @param delta : vector description of the shift.
#' 
#' @return Vector of length nbranch.
#' 
#' @seealso \code{shifts.list_to_vector}
#' 
#'17/06/14 - Initial release
##
shifts.vector_to_list <- function(delta){
  edsh <- which(delta != 0)
  shifts <- list(edges = edsh, values = delta[edsh], relativeTimes = rep(0, length(edsh)))
  return(shifts)
}

##########################################
## Handle correspondance shifts - regimes
##########################################

##
#' @title Initialisation for the computation of the optimal values
#'
#' @description
#' \code{init.compute_betas} initialize the vector of optimal values at nodes and
#' tips, with the value at the root.
#'
#' @details
#' This function is used in function \code{compute_betas} and is designed to 
#' furnish function \code{update.compute_betas} with the right structure of data.
#'
#' @param phy Input tree.
#' @param optimal.value the optimal value at the root of the tree
#' 
#' @return Matrix of size (nNodes + ntaxa)x1 of NAs, with the optiaml value
#'  at the root.
#'
#'06/10/14 - Initial release
##
init.compute_betas <- function(phy, optimal.value, ...){
  ntaxa <- length(phy$tip.label)
  beta <- matrix(nrow = 1 + nrow(phy$edge), ncol = 1) # selection strength
  beta[ntaxa + 1,] <- optimal.value
  return(beta)
}

##
#' @title Update function ofr optimal value computation
#'
#' @description
#' \code{update.compute_betas} computes the optimal value at a daughter node, 
#' knowing the optimal value at the parent node and the vector of shifts.
#'
#' @details
#' This function is used in function \code{compute_betas} and is designed to 
#' furnish function \code{recursionDown} with the right structure of data.
#'
#' @param edgeNbr : Number of the edge considered
#' @param ancestral : Computed vector for the parental node
#' @param shifts position and values of the shifts 
#' 
#' @return Updated matrix of size (nNodes + ntaxa)x1.
#'
#'06/10/14 - Initial release
##
update.compute_betas <- function(edgeNbr, ancestral, shifts, ...){
  shiftsIndex <- which(shifts$edges == edgeNbr) #If no shifts = NULL, and sum = 0
  beta <- ancestral + sum(shifts$values[shiftsIndex])
  return(beta)
}

##
#' @title Computation of the optimal values at nodes and tips.
#'
#' @description
#' \code{compute_betas} computes the optimal values at the nodes and tips of the
#' tree, given the value at the root and the list of shifts occuring in the tree.
#'
#' @details
#' This function uses function \code{recursionDown} for recursion on the tree, 
#' with \code{init.compute_betas} for the initialization of the vector, and 
#' \code{update.compute_betas} for its actualization.
#'
#' @param phylo : imput tree
#' @param optimal.value the optimal value at the root of the tree
#' @param shifts position and values of the shifts 
#' 
#' @return Vector of size (ntaxa + nNodes) of the ptimal values at the tips
#' of the tree.
#'
#'06/10/14 - Initial release
##
compute_betas <- function(phylo, optimal.value, shifts){
  phy <- reorder(phylo, order = "cladewise")
  ## Trace edges
  shifts_ordered <- shifts
  shifts_ordered$edges <- correspondanceEdges(edges=shifts$edges,from=phylo,to=phy)
  betas <- init.compute_betas(phy, optimal.value)
  betas <- recursionDown(phy, betas, update.compute_betas, shifts_ordered)
  return(betas)
}

##
#' @title Initialisation for the allocation of shifts.
#'
#' @description
#' \code{init.allocate_regimes_from_shifts} initialize the vector of regimes 
#' at nodes and tips, with the regime (0) at the root.
#'
#' @details
#' This function is used in function \code{allocate_regimes_from_shifts} and is 
#' designed to furnish function \code{update.allocate_regimes_from_shifts} with 
#' the right structure of data.
#'
#' @param phy Input tree.
#' 
#' @return Matrix of size (nNodes + ntaxa)x1 of NAs, with the 0 at the root.
#'
#'10/10/14 - Initial release
##
init.allocate_regimes_from_shifts <- function(phy, ...){
  ntaxa <- length(phy$tip.label)
  regimes <- matrix(nrow = 1 + nrow(phy$edge), ncol = 1) 
  regimes[ntaxa + 1,] <- 0
  return(regimes)
}

##
#' @title Update function for regime allocation.
#'
#' @description
#' \code{update.allocate_regimes_from_shifts} computes the regime of a daughter 
#' node, knowing the regime at the parent node and the vector of shifts positions
#'
#' @details
#' This function is used in function \code{allocate_regimes_from_shifts} and is 
#' designed to furnish function \code{recursionDown} with the right structure of
#'  data.
#'
#' @param edgeNbr : Number of the edge considered
#' @param ancestral : regime of the parent node
#' @param shifts_edges positions on edges
#' 
#' @return regime of the daughter node.
#'
#'10/10/14 - Initial release
##
update.allocate_regimes_from_shifts <- function(edgeNbr, ancestral, shifts_edges, ...){
  shiftsIndex <- which(shifts_edges == edgeNbr) # If no shifts = integer(0)
  if (length(shiftsIndex) == 0){
    return(ancestral) # No shift : keep same regime
  } else {
    return(shiftsIndex) # Shift : new regime numbered by the position in the list
  }
}


##
#' @title Allocation of regimes to nodes.
#'
#' @description
#' \code{allocate_regimes_from_shifts} allocate a number (from 0 to the number 
#' of shifts) to each node, corresponding to its regime : all nodes below shift 
#' i are numbered by i.
#' 
#' @details
#' This function uses function \code{recursionDown} for recursion on the tree, 
#' with \code{init.allocate_regimes_from_shifts} for the initialization of the
#' vector, and \code{update.allocate_regimes_from_shifts} for its actualization.
#'
#' @param phylo : imput tree
#' @param shifts_edges edges were the shifts are
#' 
#' @return Vector of size (ntaxa + nNodes) of the regimes of each node and tip.
#'
#'10/10/14 - Initial release
##
allocate_regimes_from_shifts <- function(phylo, shifts_edges){
  phy <- reorder(phylo, order = "cladewise")
  ## Trace edges
  shifts_ordered <- correspondanceEdges(edges=shifts_edges, from=phylo, to=phy)
  regimes <- init.allocate_regimes_from_shifts(phy, shifts_ordered)
  regimes <- recursionDown(phy, regimes, update.allocate_regimes_from_shifts, shifts_ordered)
  return(regimes)
}


##
#' @title Allocation of shifts to edges
#'
#' @description
#' \code{allocate_shifts_from_regimes} returns the position of the shifts induced
#' by the allocation of the regimes. Only works in an "infinite site" model.
#' 
#' @details
#' This function uses function fun on each row of matrix of edges.
#'
#' @param phylo : imput tree
#' @param regimes : vector of size (ntaxa + nNodes) of the regimes of each node
#' and tip.
#' 
#' @return vector of edges numbers where the shifts are
#'
#'10/10/14 - Initial release
##
allocate_shifts_from_regimes <- function(phylo, regimes){
  fun <- function(edge){
    if (regimes[edge[1]] == regimes[edge[2]]) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
  shifts_edges <- apply(phylo$edge, 1, fun)
  return(which(shifts_edges))
}

##
#' @title Computation of shifts from the vector of optimal values
#'
#' @description
#' \code{compute_shifts_from_betas} computes the list of shifts corresponding
#' to the vector of optimal values on nodes.
#' 
#' @details
#' This function uses function fun on each row of matrix of edges.
#'
#' @param phylo : imput tree
#' @param betas : vector of size (ntaxa + nNodes) of the optimal values at each
#' node and tip.
#' 
#' @return list of shifts
#'
#'13/10/14 - Initial release
##
compute_shifts_from_betas <- function(phylo, betas){
  fun <- function(edge){
    if (betas[edge[1]] == betas[edge[2]]) {
      return(c(0, 0))
    } else {
      return(c(1, betas[edge[2]] - betas[edge[1]]))
    }
  }
  shifts_ed <- apply(phylo$edge, 1, fun)
  edges <- which(shifts_ed[1,] == 1)
  shifts <- list(edges = edges,
                 values = shifts_ed[2, edges],
                 relativeTimes = rep(0, length(edges)))
  return(shifts)
}