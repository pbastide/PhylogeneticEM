# {shutoff}
# Copyright (C) {2014} {SR, MM, PB}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

###############################################################################
## Here are the convergence monitoring functions.
## Dependencies : generic_functions.R
###############################################################################

#######################################################
## Shutoff
#######################################################
##
# shutoff.EM.BM (params_old,params,tol)
# PARAMETERS:
# @params_old (list) : old parameters
# @params (list) : new parameters
# @tol (list) : tolerance
# RETURNS:
# (bool) should we shutoff ?
# DEPENDENCIES:
# none
# PURPOSE:
# shutoff?
# NOTES:
# TO BE DEFINED
# REVISIONS:
# 22/05/14 - Initial release
##
shutoff.EM.BM <- function(params_old, params, tol) {
  if (params_old$root.state$random) {
    return(shutoff.EM.BM.randroot(params_old,params,tol))
  } else {
    return(shutoff.EM.BM.fixedroot(params_old,params,tol))
  }
}

shutoff.EM.BM.randroot <- function(params_old, params, tol){
  if (abs(params_old$variance-params$variance)<tol$variance &&
        abs(params_old$root.state$exp.root-params$root.state$exp.root)<tol$exp.root &&
        abs(params_old$root.state$var.root-params$root.state$var.root)<tol$var.root) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

shutoff.EM.BM.fixedroot <- function(params_old, params, tol){
  if (abs(params_old$variance-params$variance)<tol$variance &&
        abs(params_old$root.state$value.root-params$root.state$value.root)<tol$value.root) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

shutoff.EM.OU <- function(stationnary.root, shifts_at_nodes, alpha_known){
  if (stationnary.root && shifts_at_nodes && alpha_known) {
    return(shutoff.EM.OU.specialCase)
  } else if (stationnary.root && shifts_at_nodes) {
    return(shutoff.EM.OU.stationnary.root_AND_shifts_at_nodes)
  } else {
    stop("The EM algorithm for the OU is only defined (for the moment) for a stationnary root and shifts at nodes !")
  }
}

shutoff.EM.OU.specialCase <- function(params_old, params, tol){
  if (abs(params_old$variance-params$variance)<tol$variance &&
        abs(params_old$root.state$exp.root-params$root.state$exp.root)<tol$exp.root &&
        abs(params_old$root.state$var.root-params$root.state$var.root)<tol$var.root) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

shutoff.EM.OU.stationnary.root_AND_shifts_at_nodes <- function(params_old, params, tol){
  if (abs(params_old$variance-params$variance)<tol$variance &&
        abs(params_old$root.state$exp.root-params$root.state$exp.root)<tol$exp.root &&
        abs(params_old$root.state$var.root-params$root.state$var.root)<tol$var.root &&
        abs(params_old$selection.strength-params$selection.strength)<tol$selection.strength) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

######################################################
## Finite parameters ?
######################################################
##
#' @title Check wether parameters are finite.
#'
#' @description
#' \code{is.finite.params} checks whether calculated parameters in the EM are
#' finite or not.
#'
#' @details
#' This function is used to test the convergence of the algorithm.
#'
#' @param params list of parameters with the correct structure
#'
#' @return boolean
#'
#'10/06/14 - Initial release
##

is.finite.params.BM <- function(params) {
  if (params$root.state$random) {
    return(is.finite.params.BM.randroot(params))
  } else {
    return(is.finite.params.BM.fixedroot(params))
  }
}

is.finite.params.BM.randroot <- function(params) {
  if (is.finite(params$variance) &&
        is.finite(params$root.state$exp.root) &&
        is.finite(params$root.state$var.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.finite.params.BM.fixedroot <- function(params) {
  if ( is.finite(params$variance) &&
         is.finite(params$root.state$value.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.finite.params.OU <- function(stationnary.root, shifts_at_nodes, alpha_known){
  if (stationnary.root && shifts_at_nodes && alpha_known) {
    return(is.finite.params.OU.specialCase)
  } else if (stationnary.root && shifts_at_nodes) {
    return(is.finite.params.OU.stationnary.root_AND_shifts_at_nodes)
  } else {
    stop("The EM algorithm for the OU is only defined (for the moment) for a stationnary root and shifts at nodes !")
  }
}

is.finite.params.OU.specialCase <- function(params) {
  if (is.finite(params$variance) &&
        is.finite(params$root.state$exp.root) &&
        is.finite(params$root.state$var.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.finite.params.OU.stationnary.root_AND_shifts_at_nodes <- function(params) {
  if (is.finite(params$variance) &&
        is.finite(params$root.state$exp.root) &&
        is.finite(params$root.state$var.root) &&
        is.finite(params$selection.strength)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#################################################################
## Are parameters in ranges ?
#################################################################
##
#' @title Check wether parameters are in ranges.
#'
#' @description
#' \code{is.in.ranges.params} checks whether calculated parameters in the EM are
#' in the defined ranges.
#'
#' @details
#' This function is used to test the convergence of the algorithm.
#'
#' @param params list of parameters with the correct structure
#' @param min list of minimum values for the parameters
#' @param max list of maximum values for the parameters
#'
#' @return boolean
#'
#'16/07/14 - Initial release
##
is.in.ranges <- function(p, min, max){
  if (p < min || p > max){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.in.ranges.params.BM <- function(params, min, max) {
  if (params$root.state$random) {
    return(is.in.ranges.params.BM.randroot(params, min, max))
  } else {
    return(is.in.ranges.params.BM.fixedroot(params, min, max))
  }
}

is.in.ranges.params.BM.randroot <- function(params, min, max) {
  if (is.in.ranges(params$variance, min$variance, max$variance) &&
        is.in.ranges(params$root.state$exp.root, min$exp.root, max$exp.root) &&
        is.in.ranges(params$root.state$var.root, min$var.root, max$var.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.in.ranges.params.BM.fixedroot <- function(params, min, max) {
  if ( is.in.ranges(params$variance, min$variance, max$variance) &&
         is.in.ranges(params$root.state$value.root, min$value.root, max$value.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.in.ranges.params.OU <- function(stationnary.root, shifts_at_nodes, alpha_known){
  if (stationnary.root && shifts_at_nodes && alpha_known) {
    return(is.in.ranges.params.OU.specialCase)
  } else if (stationnary.root && shifts_at_nodes) {
    return(is.in.ranges.params.OU.stationnary.root_AND_shifts_at_nodes)
  } else {
    stop("The EM algorithm for the OU is only defined (for the moment) for a stationnary root and shifts at nodes !")
  }
}

is.in.ranges.params.OU.specialCase <- function(params, min, max) {
  if (is.in.ranges(params$variance, min$variance, max$variance) &&
        is.in.ranges(params$root.state$exp.root, min$exp.root, max$exp.root) &&
        is.in.ranges(params$root.state$var.root, min$var.root, max$var.root) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.in.ranges.params.OU.stationnary.root_AND_shifts_at_nodes <- function(params, min, max) {
  if (is.in.ranges(params$variance, min$variance, max$variance) &&
        is.in.ranges(params$root.state$exp.root, min$exp.root, max$exp.root) &&
        is.in.ranges(params$root.state$var.root, min$var.root, max$var.root) &&
        is.in.ranges(params$selection.strength, min$selection.strength, max$selection.strength)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}