# {Make Names}
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
## Here is a set of functions to plot processes
## Dependencies : generic_functions.R
## : simulate.R
###############################################################################
##
# make.name (Name, TreeType, Y.state, Z.state, phylo, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, ...)
# PARAMETERS:
# @process (string) : which process ?
# @paramsSimu (list) : parameters used for th simulation
# @paramsEstimate (list) : parameters found by an (optional) estimation from the data
# @estimate (bool) : whether the data is issued from an estimation or a direct simulation
# @...
# RETURNS:
# (string) standardized name for the data
# DEPENDENCIES:
# catch.ProcessParams, catch.TolParams
# PURPOSE:
# Generate a standardized name for the data
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
##
# make.name <- function(process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, params_algo_EM=NULL, ...) {
#   ## Choose Process
#   catch.ProcessParams <- switch(process,
#                                 BM = catch.ProcessParams.BM,
#                                 OU = catch.ProcessParams.OU)
#   catch.TolParams <- switch(process,
#                             BM = catch.TolParams.BM,
#                             OU = catch.TolParams.OU)
#   ## Define File name
#   RootState <- paste("_root-rand=",paramsSimu$root.state$random,"_val-root=",paramsSimu$root.state$value.root,"_exp-root=",paramsSimu$root.state$exp.root,"_var-root=",paramsSimu$root.state$var.root,sep="")
#   if (is.null(paramsEstimate$shifts$edges)){
#     ShiftsState <- "_no-shift"
#   } else {
#     ShiftsState <- paste("_shifts-edges=", paste(paramsSimu$shifts$edges,collapse="-"), "_shifts-val=", paste(paramsSimu$shifts$values, collapse="-"), "_shifts-T=", paste(paramsSimu$shifts$relativeTimes,collapse="-"), sep="")
#   }
#   ProcessParams <- catch.ProcessParams(paramsSimu)
#   ## Parametrers of the EM if relevant
#   if (estimate) {
#     TolParams <- "" # catch.TolParams(params_algo_EM)
#     EstimParams <- paste(TolParams, "_process-used=", params_algo_EM$process, "_met-variance=", params_algo_EM$method.variance, "_met-init=", params_algo_EM$method.init, "_nbr-shifts=", params_algo_EM$nbr_of_shifts, sep="")
#   } else {
#     EstimParams <- ""
#   }
#   return(paste(ProcessParams, RootState, ShiftsState, EstimParams, sep=""))
# }

##
# plot.process (Name, TreeType, Y.state, Z.state, phylo, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, ...)
# PARAMETERS:
# @Name (string) begening of the name
# @TreeType (string) : type pf tree
# @Y.state (vector) : value of the traits at the tips
# @Z.state (vector) : value of the traits at the internal nodes
# @process (string) : which process ?
# @paramsSimu (list) : parameters used for th simulation
# @paramsEstimate (list) : parameters found by an (optional) estimation from the data
# @estimate (bool) : whether the data is issued from an estimation or a direct simulation
# @directory (string) : where to store the plot
# @...
# RETURNS:
# (pdf) a pdf plot of the process
# DEPENDENCIES:
# catch.LegendProcess, make.name
# PURPOSE:
# Generate a pdf file with a plot of the process
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
# 06/06.14 - Add directory
##
# plot.process <- function(Name, TreeType, Y.state, Z.state, phylo, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, position.legend="bottomleft", directory, params_algo_EM=NULL) {
#   ## Define legend
#   LegendProcess <- catch.LegendProcess(process, paramsEstimate, estimate, params_algo_EM)
#   LegendRoot <- c(#paste("Random Root = ",round(paramsEstimate$root.state$random,2),sep=""),
#     #paste("Root Value (if not random) = ",round(paramsEstimate$root.state$value.root,2),sep=""),
#     paste("Root expectation = ",round(paramsEstimate$root.state$exp.root,2),sep=""),
#     paste("Root variance = ",round(paramsEstimate$root.state$var.root,2), sep=""))
#   ## Define File Name
#   FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, params_algo_EM)
#   ## Plot
#   FileName <- paste(Name, TreeType, FileName, sep="")
#   pdf(paste(directory, FileName, ".pdf",sep=""), height=10,width=20)
#   plot.process.actual <- function(Y.state = Y.state, 
#                                   Z.state = Z.state, 
#                                   phylo = phylo, 
#                                   process = process,
#                                   paramsEstimate = paramsEstimate,
#                                   estimate = estimate,
#                                   position.legend = position.legend,
#                                   params_algo_EM = params_algo_EM)
#     legend(paste(position.legend),legend=c(LegendProcess,LegendRoot), col="black", cex = 2)
#     dev.off()
# }

# plot.process.actual <- function(Y.state, Z.state, phylo, paramsEstimate, normalize = TRUE, adj = 1, shifts_bg = "chocolate4", root_bg = "chocolate4", quant.root = 0.25, ...){
#   ntaxa <- length(phylo$tip.label)
#   if (normalize){
#     norm <- mean(abs(Y.state))
#   } else {
#     norm <- 1
#   }
#   ## Plot
#   par(mar = c(0,0,0,0), omi = c(0,0,0,0))
#   # Take care of the root
#   phylo$root.edge <- quantile(phylo$edge.length, quant.root)
#   plot(phylo, show.tip.label = FALSE, root.edge = TRUE, ...)
#   tiplabels(pch = 19, cex = abs(Y.state)/norm, col = ifelse(Y.state >= 0, "orangered", "lightblue"))
#   nodelabels(pch = 19, cex = abs(Z.state)/norm, col = ifelse(Z.state >= 0, "orangered", "lightblue"))
#   my.labeller <- function(variable, value) {
#     value <- paste(variable, "==", as.character(value))
#     value <- lapply(value, function(x) parse(text = x))
#     return(value)
#   }
#   nodelabels(text = round(paramsEstimate$optimal.value, 2), node=ntaxa + 1, bg=root_bg, cex = 1, adj = adj)
#   if ( !is.null(paramsEstimate$shifts$edges) ) {
#     edgelabels(text=round(paramsEstimate$shifts$values,2), edge=paramsEstimate$shifts$edges, bg = shifts_bg, cex = 1)
#   }
# }

## add a parameter frac to deplace the position of the label on the edge
edgelabels_home <- function (text, edge, adj = c(0.5, 0.5), frame = "rect",
                             pch = NULL, thermo = NULL, pie = NULL,
                             piecol = NULL, col = "black", bg = "lightgreen",
                             horiz = FALSE, width = NULL, height = NULL, 
                             date = NULL,
                             beg = FALSE, ...) 
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (missing(edge)) {
    sel <- 1:dim(lastPP$edge)[1]
    subedge <- lastPP$edge
  }
  else {
    sel <- edge
    subedge <- lastPP$edge[sel, , drop = FALSE]
  }
  if (lastPP$type == "phylogram") {
    if (lastPP$direction %in% c("rightwards", "leftwards")) {
      XX <- (lastPP$xx[subedge[, 1]] + lastPP$xx[subedge[, 2]])/2
      if (beg) XX <- lastPP$xx[subedge[, 1]]
      YY <- lastPP$yy[subedge[, 2]]
    }
    else {
      XX <- lastPP$xx[subedge[, 2]]
      YY <- (lastPP$yy[subedge[, 1]] + lastPP$yy[subedge[, 2]])/2
      if (beg) YY <- lastPP$yy[subedge[, 1]]
    }
  }
  else {
    XX <- (lastPP$xx[subedge[, 1]] + lastPP$xx[subedge[, 2]])/2
    if (beg) XX <- lastPP$xx[subedge[, 1]]
    YY <- (lastPP$yy[subedge[, 1]] + lastPP$yy[subedge[, 2]])/2
    if (beg) YY <- lastPP$yy[subedge[, 1]]
  }
  if (!is.null(date)) 
    XX[] <- max(lastPP$xx) - date
  BOTHlabels(text, sel, XX, YY, adj, frame, pch, thermo, pie, 
             piecol, col, bg, horiz, width, height, ...)
}

##
#' @title Plot for class \code{PhyloEM}
#'
#' @description
#' This function takes an object of class \code{PhyloEM}, result of function
#' \code{\link{PhyloEM}}, and plots the result of the inference.
#'
#' @param x an object of class \code{PhyloEM}, result of function
#' \code{\link{PhyloEM}}.
#' @param traits a vector of integers giving the numbers of the trait to be plotted.
#' Default to 1:p (all the traits).
#' @param params (optional) some user-specified parameters.
#' Must be of class \code{\link{params_process}}. If left blank, they are extracted
#' using the \code{method.selection} argument (see below).
#' @param method.selection select the parameters to plot. One of "LINselect", "DDSE",
#' "Djump". Default to "LINselect". See
#' \code{\link{params_process.PhyloEM}}.
#' @param automatic_colors whether to color the edges automatically according to
#' their regimes. Default to TRUE. If FALSE, colors can be manually specified through
#' arguments \code{color_characters} and \code{colro_edges} (see below).
#' @param color_characters if \code{automatic_colors=FALSE}, a vector of colors for
#' the tips of the tree.
#' @param color_edges if \code{automatic_colors=FALSE}, a vector of colors for the
#' edges of the tree.
#' @param plot_ancestral_states whether to plot the ancestral traits inferred at the
#' internal nodes of the tree. Only available if only one trait is plotted. Default
#' to FALSE.
#' @param name_trait name of the trait scale bar for the ancestral states plotting.
#' Default to "Trait Value".
#' @param imposed_scale if \code{plot_ancestral_states=TRUE}, a vector specifying the
#' imposed scale for the ancestral states plotting. Useful to make comparisons.
#' Default to the plotted trait.
#' @param ancestral_cex if \code{plot_ancestral_states=TRUE}, the size of the
#' ancestral states on the tree. Default to 2.
#' @param ancestral_pch if \code{plot_ancestral_states=TRUE}, the symbol used of the 
#' ancestral states. Default to circles (\code{pch=19}).
#' @param value_in_box whether to plot the value of the shift in a box on the edges.
#' Only available when only one trait is plotted. Can be difficult to read on big
#' trees. The size of the text in the boxes is controlled by parameter.
#' Default to FALSE.
#' @param ancestral_as_shift whether to represent the ancestral value at the root
#' as an ancestral shift on the root edge. Default to FALSE.
#' \code{shifts_cex} (see below).
#' @param shifts_cex if \code{value_in_box=TRUE}, the size of the text in the boxes.
#' Default to 0.8.
#' @param shifts_bg if \code{value_in_box=TRUE}, the background color of the boxes.
#' @param root_bg if \code{value_in_box=TRUE} and \code{ancestral_as_shift=TRUE},
#' the background color of the ancestral box.
#' @param shifts_adj the adj parameter for the shifts position on the edges. Default
#' to 0 (beginning of the edge).
#' @param root_adj if \code{ancestral_as_shift=TRUE}, the adj parameter for the
#' ancestral value position on the root edge. Default to 1.
#' @param color_shifts_regimes whether to color each shift according to its regime
#' (default to the same color of the edge it's on). Default to FALSE.
#' @param regime_boxes whether to draw a box showing all the tips below a given.
#' The transparency of the border of the box is controlled by parameter
#' \code{alpha_border} (see below).
#' @param alpha_border if \code{regime_boxes=TRUE}, the alpha parameter of
#' the border of the box. Default to 70.
#' @param show.tip.label whether to show the tip labels. Default to FALSE.
#' @param label_cex if \code{show.tip.label=TRUE}, the size of the labels. Default
#' to 0.5.
#' @param label_font if \code{show.tip.label=TRUE}, the font of the labels (see \link{par}).
#' @param label_offset if \code{show.tip.label=TRUE}, the size of the offset between
#' the tree and the labels. Default to 0.
#' @param axis_cex cex for the label values of the plot. Default to 0.7.
#' @param axis_las las for the label values of the plot. Default to 0 (see \link{par}).
#' @param edge.width width of the edge. Default to 1.
#' @param margin_plot vector giving the margin to around the plot.
#' Default to \code{c(0, 0, 0, 0)}.
#' @param gray_scale if TRUE, the colors are replaced by a gray scale.
#' Default to FALSE.
#' @param root.edge a logical indicating whether to draw the root edge (defaults to TRUE)
#' @param ... further arguments to be passed to \code{\link{plot.phylo}}.
#' 
#' 
#' @return
#' NULL
#' 
#' @seealso \code{\link{params_process.PhyloEM}}, \code{\link{imputed_traits.PhyloEM}}
#' 
#' @export
#'

plot.PhyloEM <- function(x,
                         traits = 1:(x$p),
                         params = NULL,
                         method.selection = NULL,
                         automatic_colors = TRUE,
                         color_characters = "black",
                         color_edges = "black",
                         plot_ancestral_states = FALSE,
                         name_trait = "Trait Value",
                         imposed_scale,
                         ancestral_cex = 2,
                         ancestral_pch = 19,
                         value_in_box = FALSE,
                         ancestral_as_shift = FALSE,
                         shifts_cex = 0.6,
                         shifts_bg = "chocolate4",
                         root_bg = "chocolate4",
                         shifts_adj = 0,
                         root_adj = 1,
                         color_shifts_regimes = FALSE,
                         regime_boxes = FALSE,
                         alpha_border = 70,
                         show.tip.label = FALSE,
                         label_cex = 0.5,
                         label_font = 1,
                         label_offset = 0,
                         axis_cex = 0.7,
                         axis_las = 0,
                         edge.width = 1,
                         margin_plot = NULL,
                         gray_scale = FALSE,
                         root.edge = TRUE,
                         ...){
  ## Checking consistency
  if (plot_ancestral_states && length(traits) > 1) stop("Ancestral state plotting is only allowed for one single trait. Please select the trait you would like to plot with argument 'traits' (see documentation).")
  if (value_in_box && length(traits) > 1) stop("Showing the shifts values on the tree is only allowed for one single trait. Please select the trait you would like to plot with argument 'traits' (see documentation).")
  
  # ## Save curent par
  # .pardefault <- par(no.readonly = T)
  # on.exit(par(.pardefault), add = TRUE)
  
  ## parameters
  if (is.null(params)){
    params <- params_process.PhyloEM(x,
                                     method.selection = method.selection)
  } else {
    if (class(params) != "params_process") {
      stop("The user specified parameters must be of class 'params_process'.")
    }
  }
  # If on trait, select relevant quantities
  if (length(traits) == 1){
    if (length(as.vector(params$selection.strength)) == 0) params$selection.strength <- 0
    if (length(as.vector(params$selection.strength)) == 1){
      if (x$p == 1){
        dim(params$selection.strength) <- c(1,1)
      } else {
        params$selection.strength <- diag(rep(params$selection.strength, x$p))
      }
    }
    params <- split_params_independent(params)
    params <- params[[traits]]
    class(params) <- "params_process"
  }
  
  ## Ancestral and imputed traits
  reconstructed_traits <- imputed_traits.PhyloEM(x, trait = traits,
                                                 save_all = TRUE,
                                                 params = params,
                                                 method.selection = method.selection)
  ancestral_states <- imputed_traits.PhyloEM(x, trait = traits,
                                             where = "nodes",
                                             method.selection = method.selection,
                                             reconstructed_states = reconstructed_traits)
  Y_state <- imputed_traits.PhyloEM(x, trait = traits,
                                    where = "tips",
                                    method.selection = method.selection,
                                    reconstructed_states = reconstructed_traits)
  rownames(Y_state) <- rownames(x$Y_data)[traits]
  if (missing(imposed_scale)) imposed_scale <- Y_state
  
  ## Plotting
  plot.data.process.actual(Y.state = Y_state,
                           phylo = x$phylo,
                           params = params,
                           process = x$process,
                           miss = is.na(x$Y_data[traits, , drop = FALSE]),
                           imposed_scale = imposed_scale,
                           root_adj = root_adj,
                           shifts_adj = shifts_adj,
                           shifts_bg = shifts_bg,
                           root_bg = root_bg,
                           quant.root = 0.25,
                           color_characters = color_characters,
                           color_edges = color_edges,
                           edge.width = edge.width,
                           automatic_colors = automatic_colors,
                           regime_boxes = regime_boxes,
                           alpha_border = alpha_border,
                           value_in_box = value_in_box,
                           shifts_cex = shifts_cex,
                           axis_cex = axis_cex,
                           axis_las = axis_las,
                           margin_plot = margin_plot,
                           color_shifts_regimes = color_shifts_regimes,
                           # shifts_regimes = shifts_regimes,
                           plot_ancestral_states = plot_ancestral_states,
                           ancestral_states = ancestral_states,
                           name_trait = name_trait,
                           # imposed_scale.nodes = imposed_scale.nodes,
                           ancestral_cex = ancestral_cex,
                           ancestral_pch = ancestral_pch,
                           label_cex = label_cex,
                           label_font = label_font,
                           show.tip.label = show.tip.label,
                           # underscore = underscore,
                           # label.offset = label.offset,
                           ancestral_as_shift = ancestral_as_shift,
                           gray_scale = gray_scale,
                           root.edge = root.edge,
                           ...)
}

plot.data.process.actual <- function(Y.state, phylo, params,
                                     miss = is.na(Y.state),
                                     process = "BM",
                                     #norm = max(abs(Y.state)),
                                     imposed_scale = Y.state,
                                     root_adj = 1, shifts_adj = 0,
                                     shifts_bg = "chocolate4",
                                     root_bg = "chocolate4", quant.root = 0.25,
                                     color_characters = "black",
                                     color_edges = "black",
                                     edge.width = 1,
                                     automatic_colors = FALSE,
                                     regime_boxes = FALSE,
                                     alpha_border = 70,
                                     value_in_box = TRUE,
                                     shifts_cex = 1,
                                     axis_cex = 0.7,
                                     axis_las = 0,
                                     margin_plot = NULL,
                                     color_shifts_regimes = FALSE,
                                     # shifts_regimes = NULL,
                                     plot_ancestral_states = FALSE,
                                     ancestral_states = NULL,
                                     name_trait = "Trait Value",
                                     imposed_scale.nodes = ancestral_states,
                                     ancestral_cex = 2,
                                     ancestral_pch = 19,
                                     label_cex = 1,
                                     label_font = 1,
                                     show.tip.label = FALSE,
                                     underscore = FALSE,
                                     label.offset = 0,
                                     ancestral_as_shift = TRUE,
                                     gray_scale = FALSE,
                                     root.edge = TRUE,
                                     ...){
  # ## Save curent par
  # .pardefault <- par(no.readonly = T)
  # on.exit(par(.pardefault), add = TRUE)
  
  ntaxa <- length(phylo$tip.label)
  p_dim <- nrow(Y.state)
  if (is.null(p_dim)) p_dim <- 0
  #   if (normalize){
  #     norm <- max(abs(Y.state))
  #   } else {
  #     norm <- 1
  #   }
  Y.state <- Y.state # / norm
  unit <- 1 # / norm
  ## Root state
  if (process == "OU" || is.null(params$root.state)){
    root.val <- params$optimal.value
  } else {
    if (params$root.state$random){
      root.val <- params$root.state$exp.root
    } else {
      root.val <- params$root.state$value.root
    }
  }
  ## Automatic colors
  if (automatic_colors){
    nodes_regimes <- allocate_regimes_from_shifts(phylo,
                                                  params$shifts$edges)
    
    color_characters <- as.factor(nodes_regimes[1:ntaxa])
    if (!gray_scale){
      levels(color_characters) <- c("black", 
                                    rainbow(length(levels(color_characters)) - 1,
                                            start = 0, v = 0.5))
    } else {
      levels(color_characters) <- gray.colors(length(levels(color_characters)),
                                              start = 0, end = 0.8)
    }
    
    color_edges <- as.factor(nodes_regimes[phylo$edge[, 2]])
    if (!gray_scale){
      levels(color_edges) <- c("black", rainbow(length(levels(color_edges)) - 1,
                                                start = 0, v = 0.5))
    } else {
      levels(color_edges) <- gray.colors(length(levels(color_edges)),
                                         start = 0, end = 0.8)
    }
  }
  ## Plot ancestral states ?
  if (plot_ancestral_states){
    if (!requireNamespace("phytools", quietly = TRUE)) {
      stop("phytools is needed for plotting ancestral states. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("graphics", quietly = TRUE)) {
      stop("graphics is needed for plotting ancestral states. Please install it.",
           call. = FALSE)
    }
    # library(phytools)
    # library(graphics)
    if (is.null(ancestral_states)){
      warning("Plot option clash: the ancestral states could not be plotted (please provide values).")
    } else {
      imp.scale.nodes  <- range(c(imposed_scale, imposed_scale.nodes), na.rm = TRUE)
      map2color <- function(x, pal, limits = NULL) {
        if (is.null(limits)) {
          limits = range(x)
        }
        pal[findInterval(x,
                         seq(limits[1], limits[2], length.out = 1000 + 1),
                         all.inside = TRUE)]
      }
      pal <- rev(palette(rainbow(1001, start = 0, end = 0.7)))
      pal <- rev(palette(rainbow(1001, start = 0, end = 0.7)))
      col_ancestral <- map2color(ancestral_states, pal = pal, limits = imp.scale.nodes)
      # If plotting ancestral, colors of the tips values to match colors of the palette
      if (!is.null(Y.state)){
        color_characters_regimes <- color_characters
        color_characters <- map2color(as.vector(Y.state),
                                      pal, limits = imp.scale.nodes)
      }
    }
  }
  ## Plot
  if (!is.null(margin_plot)){
    opar <- par("mar", "omi")
    par(mar = margin_plot, omi = margin_plot)
    on.exit(par(opar))
  }
  # Take care of the root
  phylo$root.edge <- quantile(phylo$edge.length, quant.root)
  # Plot tree
  if (is.null(Y.state)){
    plot(phylo, show.tip.label = show.tip.label, root.edge = root.edge, 
         edge.color = as.vector(color_edges),
         edge.width = edge.width, ...)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  } else {
    h_p <- max(ape::node.depth.edgelength(phylo))
    if (show.tip.label){
      size_labels <- h_p / 4 
    } else {
      size_labels <- 0
    }
    x.lim.max <- h_p + p_dim * h_p / 3 + size_labels
    y.lim.min <- -ntaxa/10
    y.lim.max <- ntaxa + ntaxa/10
    plot(phylo, show.tip.label = FALSE, root.edge = root.edge, 
         x.lim = c(0, x.lim.max), 
         y.lim = c(y.lim.min, y.lim.max),
         edge.color = as.vector(color_edges),
         edge.width = edge.width, ...)
    if (show.tip.label){
      size_labels <- max(strwidth(phylo$tip.label, cex = label_cex))
    }
    # Plot data at tips
    ## length available for character plotting
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    pos_last_tip <- max(lastPP$xx)
    # label.offset <- 1/8 * (x.lim.max - pos_last_tip - size_labels)
    available_x <- x.lim.max - pos_last_tip - size_labels
    ell <- available_x / (p_dim +  (p_dim + 1) / 3)# lenght for the plot of one character
    offset <- ell / 3
    
    ## Plots characters
    for (t in 1:p_dim){
      imp.scale  <- c(min(0, min(imposed_scale[t, ], na.rm = TRUE)),
                      max(0, max(imposed_scale[t, ], na.rm = TRUE)))
      mult <- ell / (imp.scale[2] - imp.scale[1])
      Y.plot <- mult * Y.state[t, ]
      unit <- mult * unit
      minY <- min(Y.plot, na.rm = TRUE)
      maxY <- max(Y.plot, na.rm = TRUE)
      eccart_g <- -min(minY, 0) + offset
      # 0 bar
      segments(pos_last_tip + eccart_g, y.lim.min + ntaxa/15,
               pos_last_tip + eccart_g, y.lim.max - ntaxa/10,
               lty = 3)
    # text(pos_last_tip + eccart_g, y.lim.max - 2*ntaxa/30,
    #      "0", cex = lastPP$cex)
      if (!is.null(rownames(Y.state))){
        text(pos_last_tip + eccart_g, y.lim.max - 2*ntaxa/30,
             rownames(Y.state)[t], cex = axis_cex)
      }
      # unit length
      axis(1, at = pos_last_tip + eccart_g + range(Y.plot, na.rm = TRUE),
           labels = round(range(Y.state[t, ], na.rm = TRUE), digits = 2),
           pos = y.lim.min + ntaxa/15, 
           cex.axis = axis_cex,
           # padj = -0.5,
           las = axis_las)
      # segments(pos_last_tip + eccart_g, y.lim.min + ntaxa/15,
      #          pos_last_tip + eccart_g + unit, y.lim.min + ntaxa/15,
      #          lwd = 2)
      # text(pos_last_tip + eccart_g, y.lim.min + ntaxa/15,
      #      "Unit", cex = lastPP$cex,
      #      pos = 2)
      # characters
      segments(pos_last_tip + eccart_g, lastPP$yy[1:ntaxa][!miss[t, ]],
               pos_last_tip + eccart_g + Y.plot[!miss[t, ]], lastPP$yy[1:ntaxa][!miss[t, ]],
               col = as.vector(color_characters)[!miss[t, ]],
               lwd = edge.width)
      # missing ones as dotted
      if (any(miss[t, ])){
        segments(pos_last_tip + eccart_g, lastPP$yy[1:ntaxa][miss[t, ]],
                 pos_last_tip + eccart_g + Y.plot[miss[t, ]], lastPP$yy[1:ntaxa][miss[t, ]],
                 col = as.vector(color_characters)[miss[t, ]],
                 lwd = edge.width,
                 lty = 3) 
      }
    
      # report for next
      pos_last_tip <- pos_last_tip + ell + offset
    }
    ## Tip Labels
    if (show.tip.label){
      if (is.expression(phylo$tip.label)) underscore <- TRUE
      if (!underscore) phylo$tip.label <- gsub("_", " ", phylo$tip.label)
      x.lim.max.data <- max(pos_last_tip, na.rm = TRUE) + label.offset
      if (!exists("color_characters_regimes")) color_characters_regimes <- color_characters
      text(x.lim.max.data, lastPP$yy[1:ntaxa], phylo$tip.label, 
           cex = label_cex, pos = 4,
           col = as.vector(color_characters_regimes),
           font = label_font)
    }
  }
  ## Ancestral states
  if (plot_ancestral_states){
    nodelabels(pch = ancestral_pch, cex = ancestral_cex, col = col_ancestral)
    leg <- 0.5 * ape::node.depth.edgelength(phylo)[1]
    phytools::add.color.bar(leg, pal, title = name_trait,
                            lims = imp.scale.nodes,
                            digits = 2, prompt = FALSE,
                            lwd = 4, outline = TRUE,
                            x = 0,
                            y = 0.8 * par()$usr[3],
                            subtitle = "")
  }
  ## Plot beta_0
  if (value_in_box){ # Write value of shift in the box
    if (!is.null(root.val) && ancestral_as_shift){
      nodelabels(text = round(root.val, 1), 
                 node = ntaxa + 1,
                 bg = root_bg,
                 cex = shifts_cex,
                 adj = root_adj)
    }
    # Plot shifts
    if ( !is.null(params$shifts$edges) ) {
      edgelabels_home(text = round(params$shifts$values, 1), 
                      edge = params$shifts$edges, 
                      bg = shifts_bg,
                      cex = shifts_cex,
                      beg = TRUE,
                      adj = shifts_adj)
    }
  } else {
    if (color_shifts_regimes){ # Shift has one color for each regime
      nodes_regimes  <-  compute_betas_from_shifts(phylo, 
                                                   root.val,
                                                   params$shifts)
      color_edges <- as.factor(nodes_regimes[phylo$edge[, 2]])
      levels(color_edges) <- c("black", rainbow(length(levels(color_edges)) - 1,
                                                start = 0, v = 0.5))
      col_shifts <- as.vector(color_edges[params$shifts$edges])
      edgelabels_home(text = rep("", length(col_shifts)),
                      edge = params$shifts$edges, 
                      frame = "circle",
                      cex = shifts_cex,
                      bg = col_shifts,
                      beg = TRUE)
    } else if (p_dim == 1) { # Color code for shifts values
      values <- c(root.val, params$shifts$values)
      if (!gray_scale){
        col_shifts <- color_palette(values)
      } else {
        col_shifts <- gray(c(0.7, 0.3))[(values >= 0) + 1]
      }
      if (!is.null(root.val) && ancestral_as_shift){
        nodelabels(text = "", 
                   node = ntaxa + 1,
                   frame = "circle",
                   cex = shifts_cex,
                   bg = col_shifts[1])
      }
      col_shifts <- col_shifts[-1]
      # Plot shifts
      if ( !is.null(params$shifts$edges) ) {
        edgelabels_home(text = rep("", length(col_shifts)),
                        edge = params$shifts$edges, 
                        frame = "circle",
                        cex = shifts_cex,
                        bg = col_shifts,
                        beg = TRUE)
      }
    } else {
      # Plot shifts
      if ( !is.null(params$shifts$edges) ) {
        edgelabels_home(text = rep("", length(params$shifts$edges)),
                        edge = params$shifts$edges, 
                        frame = "circle",
                        cex = shifts_cex,
                        bg = "black",
                        beg = TRUE)
      }
    }
  }
  ## Boxes aroud regimes
  if (regime_boxes){
    nodes_regimes <- allocate_regimes_from_shifts(phylo,
                                                  params$shifts$edges)
    tips_regimes <- nodes_regimes[1:ntaxa]
    all_regimes <- 1:max(tips_regimes)
    for (reg in all_regimes){
      groupe <- which(tips_regimes == reg)
      prac <- getMRCA(phylo, groupe)
      if(is.null(prac)){
        prac_fa <- groupe # If only one tip in the group
        edge <- which(phylo$edge[,2]==prac_fa)
      } else {
        edge <- which(phylo$edge[,2]==prac)
        prac_fa <- phylo$edge[edge,1]
      }
      rect(lastPP$xx[prac_fa] + 0.5 * phylo$edge.length[edge],
           lastPP$yy[min(groupe)] - 0.5,
           lastPP$x.lim[2] + 2,
           lastPP$yy[max(groupe)] + 0.5,
           lwd = 2,
           border = paste0("#000000", alpha_border))
    }
  }
}

##
#' @title Plot Model Selection Criterion
#'
#' @description
#' This function takes an object of class \code{PhyloEM}, result of function
#' \code{\link{PhyloEM}}, and plots a model selection criterion.
#'
#' @param res an object of class \code{PhyloEM}, result of function
#' \code{\link{PhyloEM}}.
#' @param method.selection select the parameters to plot. One of "LINselect", "DDSE",
#' "Djump" or "likelihood" (for un-penalized likelihood). Default to "LINselect". See
#' \code{\link{params_process.PhyloEM}}.
#' @param add boolean: should the points be added to a urent plot (default to FALSE).
#' @param select.col the color of the point selected by the criterion. Default to "red".
#' @param ... further argument to be passed to base \code{\link{plot}}.
#' 
#' @return
#' NULL
#' 
#' @seealso \code{\link{params_process.PhyloEM}}, \code{\link{plot.PhyloEM}}, \code{\link{get_criterion}}
#' 
#' @export
#'
plot_criterion <- function(res, method.selection = NULL, add = FALSE, select.col = "red", ...) {
  m_sel <- get_method_selection(res, method.selection = method.selection)
  K_grid <- res[[m_sel[2]]]$results_summary$K_try
  if (m_sel[1] == "log_likelihood") {
    name_crit <- "log_likelihood"
  } else {
    name_crit <- paste0("crit_", m_sel[1])
  }
  Crit <- res[[m_sel[2]]]$results_summary[[name_crit]]
  if (add) {
    graphics::points(K_grid, Crit, ...)
  } else {
    plot(K_grid, Crit, xlab = "K", ylab = paste0("Criterion ", m_sel[3]), ...)
  }
  if (m_sel[4] != "no_max_min") {
    if (m_sel[4] == "min") {
      selected <- which.min(Crit)
    } else {
      selected <- which.max(Crit)
    }
    graphics::points(K_grid[selected], Crit[selected], col = select.col, ...)
  }
}

##
#' @title Get Model Selection Criterion
#'
#' @description
#' This function takes an object of class \code{PhyloEM}, result of function
#' \code{\link{PhyloEM}}, and return the values of the model selection criterion
#' for each value of K.
#'
#' @param res an object of class \code{PhyloEM}, result of function
#' \code{\link{PhyloEM}}.
#' @param method.selection select the parameters to plot. One of "LINselect", "DDSE",
#' "Djump" or "likelihood" (for un-penalized likelihood). Default to "LINselect". See
#' \code{\link{params_process.PhyloEM}}.
#' 
#' @return
#' A named vector with the values of the criterion for each number of shift K.
#' 
#' @seealso \code{\link{params_process.PhyloEM}}, \code{\link{plot.PhyloEM}}, \code{\link{plot_criterion}}
#' 
#' @export
#'
get_criterion <- function(res, method.selection = NULL) {
  m_sel <- get_method_selection(res, method.selection = method.selection)
  K_grid <- res[[m_sel[2]]]$results_summary$K_try
  if (m_sel[1] == "log_likelihood") {
    name_crit <- "log_likelihood"
  } else {
    name_crit <- paste0("crit_", m_sel[1])
  }
  Crit <- res[[m_sel[2]]]$results_summary[[name_crit]]
  names(Crit) <- K_grid
  return(Crit)
}

# save.process <- function(Name, TreeType, XX, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, directory, ...) {
#   ## Define File Name
#   FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, ...)
#   FileName <- paste(Name, TreeType, FileName, sep="")
#   save(XX, paramsSimu, paramsEstimate, file=paste(directory, FileName, ".RData",sep=""))
# }

# write.table.results <- function(Name, TreeType, res, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, directory, ...) {
#   ## Define File Name
#   FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, ...)
#   FileName <- paste(Name, TreeType, FileName, sep="")
#   write.table(res, paste(directory, FileName, ".csv",sep=""))
# }

# load.process <- function(Name, TreeType, XX, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, directory, ...) {
#   ## Define File Name
#   FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, ...)
#   FileName <- paste(Name, TreeType, FileName, sep="")
#   load(file=paste(directory, FileName, ".RData",sep=""))
# }

##
# catch.ProcessParams (paramsSimu)
# PARAMETERS:
# @paramsSimu (list) : parameters used for th simulation
# RETURNS:
# (string) string containing all the informations on the parameters used for the simulation
# DEPENDENCIES:
# none
# PURPOSE:
# generate adequate string
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
##
# catch.ProcessParams.BM <- function(paramsSimu){
#   return(paste("_process=BM","_var=",paramsSimu$variance,sep=""))
# }

# catch.ProcessParams.OU <- function(paramsSimu){
#   return(paste("_process=OU","_var=",paramsSimu$variance,"_opt-val=",paramsSimu$optimal.value,"_sel-strength=",paramsSimu$selection.strength,sep=""))
# }

##
# catch.LegendProcess (paramsEstimate)
# PARAMETERS:
# @paramsEstimate (list) : parameters estimated from the data
# RETURNS:
# (vector) strings containing all the informations on the parameters used for the simulation
# DEPENDENCIES:
# none
# PURPOSE:
# generate adequate vector of string for a legend
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
##
# catch.LegendProcess <- function(process, paramsEstimate, estimate, params_algo_EM=NULL) {
#   if (estimate) {
#     proc <- params_algo_EM$process
#   } else {
#     proc <- process
#   }
#   if (proc=="BM") {
#     return(catch.LegendProcess.BM(paramsEstimate))
#   } else if (proc=="OU") {
#     return(catch.LegendProcess.OU(paramsEstimate))
#   }
# }
# catch.LegendProcess.BM <- function(paramsEstimate){
#   return( c("Process : BM",
#             paste("Process Variance = ", round(paramsEstimate$variance,2), sep="")) )
# }
# 
# catch.LegendProcess.OU <- function(paramsEstimate){
#   return( c("Process : OU",
#             paste("Process Variance = ", round(paramsEstimate$variance,2),sep=""),
#             paste("Beta_0 = ", round(paramsEstimate$optimal.value,2), sep=""),
#             paste("Selection Strength = ", round(paramsEstimate$selection.strength,2), sep="")) )
# }

##
# catch.TolParams (params_algo_EM)
# PARAMETERS:
# @params_algo_EM (list) : parameters of the EM algorithm used
# RETURNS:
# (string) string containing all the informations on the parameters used for the tolerance of the parameters in the estimation
# DEPENDENCIES:
# none
# PURPOSE:
# generate adequate string
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
##
# catch.TolParams.BM <- function(params_algo_EM){
#   return(paste("_tol-variance=", params_algo_EM$tol$variance,
#                "_tol-exp-root=", params_algo_EM$tol$exp.root,
#                "_tol-var-root=", params_algo_EM$tol$var.root,sep=""))
# }
# 
# catch.TolParams.OU <- function(params_algo_EM){
#   return(paste("_tol-variance=", params_algo_EM$tol$variance,
#                "_tol-exp-root=", params_algo_EM$tol$exp.root,
#                "_tol-var-root=", params_algo_EM$tol$var.root,
#                "_tol-optim-value=", params_algo_EM$tol$optim.value,
#                "_tol-selection-strength=", params_algo_EM$tol$selection.strength,sep=""))
# }

#####################################################################
## Plot the history of the estimations
#####################################################################

# list_to_table.history <- function(params_history) {
#   ll <- unlist(sapply(params_history, function(x) attr(x, "log_likelihood")[1]))
# #  ll_bis <- unlist(sapply(params_history, function(x) attr(x, "log_likelihood_bis")[1]))
#   #method <- unlist(sapply(params_history, function(x) attr(x, "segmentation_algorithm_used")))
#   nbr_of_shifts <- length(params_history[['1']]$shifts$edges)
#   params_history[['0']] <- replaceInList(params_history[['0']], function(x) if(is.null(x))rep(0,nbr_of_shifts) else x)
#   params_history <- lapply(params_history, unlist)
#   history <- do.call(cbind, params_history)
#   history <- rbind(history, log_likelihood = c(ll))#, log_likelihood_bis = c(ll_bis, NA)), segmentation_algorithm = c(method, NA, NA))
#   return(history)
# }

# write.table.history <- function(history, params_algo_EM, PATH, ...) {
#   ## Define File Name
#   name <- paste(PATH, "history_parameters", "_init=", params_algo_EM$method.init, "_initalpha=", params_algo_EM$method.init.alpha, "_nbrofshifts=", params_algo_EM$nbr_of_shifts, ".csv", sep="")
#   write.csv2(history, name, ...)
# }

# plot.history.OU.stationary <- function(params_history, tree, params_ref, Y_data_ref, PATH, name, ref = "true"){
#   if (missing(params_ref)){
#     ref <- NULL
#     history <- list_to_table.history(params_history)
#   } else {
#     params_history[[ref]] <- params_ref
#     history <- list_to_table.history(params_history)
#     history[,ref]["log_likelihood"] <-log_likelihood.OU(Y_data_ref, 
#                                                         tree, 
#                                                         params_ref)
#   }
#   #params_simu  <-  unlist(paramsSimu)[names(history[,1])]
#   pdf(paste(PATH, name, ".pdf", sep=""), width = 12, height = 8)
#   ## Title of the page
#   #title <- paste("Initialization : ", params_algo_EM$method.init, "\n", "Alpha Initialization : ", params_algo_EM$method.init.alpha, sep="")
#   ## Plot
#   plot.history.OU.stationary.actual(history, ref)
#   dev.off()
# }
# 
# plot.history.OU.stationary.actual <- function (history, ref = "true", title = "Parameters estimations and log-likelihood of the model at each step of the EM algorithm.") {
#   ## Discriminate
#   if (!is.null(ref)){
#     params_simu <- history[, ref]
#     history <- history[, colnames(history) != ref]
#   }
#   ## Create grid
#   pushViewport(viewport(layout = grid.layout(2+1, 3, heights = unit(c(1,5,5), "null"))))
#   ## title
#   grid.text(title, vp = vplayout(1,1:3))
#   ## Continuous parameters
#   row <- c(2,2,2,3); col <- c(1,2,3,3)
#   params_to_plot <- c("selection.strength", "root.state.var.root", "optimal.value", "log_likelihood")
#   params_to_plot_legend <- c(expression(alpha), expression(gamma^2), expression(beta[0]), "log likelihood")
#   for (s in 1:4) {
#     # Choose right score
#     df <- history[params_to_plot[s],]
#     # Plot
#     p <- qplot(seq_along(df)-1, df, xlab="iterations", ylab=params_to_plot_legend[s])
#     if (!is.null(ref)) p <- p + geom_hline(yintercept = params_simu[params_to_plot[s]])
#     print(p, vp=vplayout(row[s],col[s]))
#   }
#   ## Plot with differents colours for each shift
#   # Values
#   df_val <- as.data.frame(t(history[grepl('shifts.values', names(history[,1])),, drop=F]))
#   df_val_long <- melt(df_val, variable.name = "shift", value.name="shift.value")
#   df_val_long$iterations <- rep(seq_along(df_val[,1])-1, dim(df_val)[2])
#   df_val_long$shift <- as.factor(rep(seq_along(df_val[1,]), each=dim(df_val)[1]))
#   if (!is.null(ref)) df_val_long$true <- rep(params_simu[grepl('shifts.values', names(history[,1]))], each=dim(df_val)[1])
#   p <- ggplot(data=df_val_long, aes(x=iterations, y=shift.value, colour = shift)) + geom_point()
#   if (!is.null(ref)) p <- p + geom_hline(aes(yintercept=true, colour=shift), data=df_val_long)
#   p <- p + ylab(expression(delta))
#   print(p, vp=vplayout(3,1))
#   # Edges
#   df_ed <- as.data.frame(t(history[grepl('shifts.edges', names(history[,1])),, drop=F]))
#   df_ed_long <- melt(df_ed, variable.name = "shift", value.name="shift.edge")
#   df_ed_long$iterations <- rep(seq_along(df_ed[,1])-1, dim(df_ed)[2])
#   df_ed_long$shift <- as.factor(rep(seq_along(df_ed[1,]), each=dim(df_ed)[1]))
#   if (!is.null(ref)) df_ed_long$true <- rep(params_simu[grepl('shifts.edges', names(history[,1]))], each=dim(df_ed)[1])
#   p <- ggplot(data=df_ed_long, aes(x=iterations, y=shift.edge, colour = shift)) + geom_point()
#   if (!is.null(ref)) p <- p + geom_hline(aes(yintercept=true, colour=shift), data=df_ed_long)
#   p <- p + ylab(expression(tau))
#   print(p, vp=vplayout(3,2))
# }

# vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
# dtoc <- function(x) gsub(".", ",", x, fixed=TRUE)

color_palette <- function (values) {
  indPos <- values > 0
  indNeg <- values < 0
  indNull <- values == 0
  nbrColPos <- sum(indPos)
  nbrColNeg <- sum(indNeg)
  nbrColNull <- sum(indNull)
  palettePos <- grDevices::colorRampPalette(c("orangered", "red"))(nbrColPos)
  paletteNeg <- grDevices::colorRampPalette(c("blue", "lightblue"))(nbrColNeg)
  col_shifts <- rep(NA, length(values))
  if (nbrColPos <= 1){
    col_shifts[indPos] <- palettePos
  } else {
    col_shifts[indPos] <- palettePos[cut(values[indPos], nbrColPos)]
  }
  if (nbrColNeg <= 1){
    col_shifts[indNeg] <- paletteNeg
  } else {
    col_shifts[indNeg] <- paletteNeg[cut(values[indNeg], nbrColNeg)]
  }
  col_shifts[indNull] <- "white"
  return(col_shifts)
}

# # Function to plot color bar
# # (http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette)
# color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
#   scale = (length(lut)-1)/(max-min)
#   
#   dev.new(width=1.75, height=5)
#   plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
#   axis(2, ticks, las=1)
#   for (i in 1:(length(lut)-1)) {
#     y = (i-1)/scale + min
#     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
#   }
# }