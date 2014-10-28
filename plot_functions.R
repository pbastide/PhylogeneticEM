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
# @estimate (bool) : wether the data is issued from an estimation or a direct simulation
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
make.name <- function(process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, params_algo_EM=NULL, ...) {
  ## Choose Process
  catch.ProcessParams <- switch(process,
                                BM = catch.ProcessParams.BM,
                                OU = catch.ProcessParams.OU)
  catch.TolParams <- switch(process,
                            BM = catch.TolParams.BM,
                            OU = catch.TolParams.OU)
  ## Define File name
  RootState <- paste("_root-rand=",paramsSimu$root.state$random,"_val-root=",paramsSimu$root.state$value.root,"_exp-root=",paramsSimu$root.state$exp.root,"_var-root=",paramsSimu$root.state$var.root,sep="")
  if (is.null(paramsEstimate$shifts$edges)){
    ShiftsState <- "_no-shift"
  } else {
    ShiftsState <- paste("_shifts-edges=", paste(paramsSimu$shifts$edges,collapse="-"), "_shifts-val=", paste(paramsSimu$shifts$values, collapse="-"), "_shifts-T=", paste(paramsSimu$shifts$relativeTimes,collapse="-"), sep="")
  }
  ProcessParams <- catch.ProcessParams(paramsSimu)
  ## Parametrers of the EM if relevent
  if (estimate) {
    TolParams <- "" # catch.TolParams(params_algo_EM)
    EstimParams <- paste(TolParams, "_process-used=", params_algo_EM$process, "_met-variance=", params_algo_EM$method.variance, "_met-init=", params_algo_EM$method.init, "_nbr-shifts=", params_algo_EM$nbr_of_shifts, sep="")
  } else {
    EstimParams <- ""
  }
  return(paste(ProcessParams, RootState, ShiftsState, EstimParams, sep=""))
}
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
# @estimate (bool) : wether the data is issued from an estimation or a direct simulation
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
plot.process <- function(Name, TreeType, Y.state, Z.state, phylo, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, position.legend="bottomleft", directory, params_algo_EM=NULL) {
  ## Define File Name
  FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, params_algo_EM)
  ## Define legend
  LegendProcess <- catch.LegendProcess(process, paramsEstimate, estimate, params_algo_EM)
  LegendRoot <- c(#paste("Random Root = ",round(paramsEstimate$root.state$random,2),sep=""),
    #paste("Root Value (if not random) = ",round(paramsEstimate$root.state$value.root,2),sep=""),
    paste("Root expectation = ",round(paramsEstimate$root.state$exp.root,2),sep=""),
    paste("Root variance = ",round(paramsEstimate$root.state$var.root,2), sep=""))
  ## Plot
  FileName <- paste(Name, TreeType, FileName, sep="")
  pdf(paste(directory, FileName, ".pdf",sep=""), height=10,width=20)
  plot(phylo, show.tip.label = FALSE)
  tiplabels(pch = 19, cex = abs(Y.state)/mean(abs(Y.state)), col = ifelse(Y.state >= 0, "orangered", "lightblue"))
  nodelabels(pch = 19, cex = abs(Z.state)/mean(abs(Z.state)), col = ifelse(Z.state >= 0, "orangered", "lightblue"))
  if ( !is.null(paramsEstimate$shifts$edges) ) {
    edgelabels(text=round(paramsEstimate$shifts$values,2), edge=paramsEstimate$shifts$edges, bg="chocolate4", cex=2.5)
  }
  legend(paste(position.legend),legend=c(LegendProcess,LegendRoot), col="black", cex = 2)
  dev.off()
}
save.process <- function(Name, TreeType, XX, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, directory, ...) {
  ## Define File Name
  FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, ...)
  FileName <- paste(Name, TreeType, FileName, sep="")
  save(XX, paramsSimu, paramsEstimate, file=paste(directory, FileName, ".RData",sep=""))
}
write.table.results <- function(Name, TreeType, res, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, directory, ...) {
  ## Define File Name
  FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, ...)
  FileName <- paste(Name, TreeType, FileName, sep="")
  write.table(res, paste(directory, FileName, ".csv",sep=""))
}
load.process <- function(Name, TreeType, XX, process = c("BM", "OU"), paramsSimu, paramsEstimate=paramsSimu, estimate=FALSE, directory, ...) {
  ## Define File Name
  FileName <- make.name(process, paramsSimu, paramsEstimate, estimate, ...)
  FileName <- paste(Name, TreeType, FileName, sep="")
  load(file=paste(directory, FileName, ".RData",sep=""))
}
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
catch.ProcessParams.BM <- function(paramsSimu){
  return(paste("_process=BM","_var=",paramsSimu$variance,sep=""))
}
catch.ProcessParams.OU <- function(paramsSimu){
  return(paste("_process=OU","_var=",paramsSimu$variance,"_opt-val=",paramsSimu$optimal.value,"_sel-strength=",paramsSimu$selection.strength,sep=""))
}
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
catch.LegendProcess <- function(process, paramsEstimate, estimate, params_algo_EM=NULL) {
  if (estimate) {
    proc <- params_algo_EM$process
  } else {
    proc <- process
  }
  if (proc=="BM") {
    return(catch.LegendProcess.BM(paramsEstimate))
  } else if (proc=="OU") {
    return(catch.LegendProcess.OU(paramsEstimate))
  }
}
catch.LegendProcess.BM <- function(paramsEstimate){
  return( c("Process : BM",
            paste("Process Variance = ", round(paramsEstimate$variance,2), sep="")) )
}
catch.LegendProcess.OU <- function(paramsEstimate){
  return( c("Process : OU",
            paste("Process Variance = ", round(paramsEstimate$variance,2),sep=""),
            paste("Beta_0 = ", round(paramsEstimate$optimal.value,2), sep=""),
            paste("Selection Strength = ", round(paramsEstimate$selection.strength,2), sep="")) )
}
##
# catch.TolParams (params_algo_EM)
# PARAMETERS:
# @params_algo_EM (list) : parameters of the EM algorithm used
# RETURNS:
# (string) string containing all the informations on the parameters used for the tolerence of the parameters in the estimation
# DEPENDENCIES:
# none
# PURPOSE:
# generate adequate string
# NOTES:
# none
# REVISIONS:
# 26/05/14 - Initial release
} ##
catch.TolParams.BM <- function(params_algo_EM){
  return(paste("_tol-variance=", params_algo_EM$tol$variance,
               "_tol-exp-root=", params_algo_EM$tol$exp.root,
               "_tol-var-root=", params_algo_EM$tol$var.root,sep=""))
}
catch.TolParams.OU <- function(params_algo_EM){
  return(paste("_tol-variance=", params_algo_EM$tol$variance,
               "_tol-exp-root=", params_algo_EM$tol$exp.root,
               "_tol-var-root=", params_algo_EM$tol$var.root,
               "_tol-optim-value=", params_algo_EM$tol$optim.value,
               "_tol-selection-strength=", params_algo_EM$tol$selection.strength,sep=""))
}
#####################################################################
## Plot the history of the estimations
#####################################################################
list_to_table.history <- function(params_history) {
  ll <- unlist(sapply(params_history, function(x) attr(x, "log_likelihood")[1]))
  ll_bis <- unlist(sapply(params_history, function(x) attr(x, "log_likelihood_bis")[1]))
  method <- unlist(sapply(params_history, function(x) attr(x, "segmentation_algorithm_used")))
  nbr_of_shifts <- length(params_history[['1']]$shifts$edges)
  params_history[['0']] <- replaceInList(params_history[['0']], function(x) if(is.null(x))rep(0,nbr_of_shifts) else x)
  params_history <- lapply(params_history, unlist)
  history <- do.call(cbind, params_history)
  history <- rbind(history, log_likelihood = c(ll, NA), log_likelihood_bis = c(ll_bis, NA), segmentation_algorithm = c(method, NA, NA))
  return(history)
}

write.table.history <- function(history, params_algo_EM, PATH, ...) {
  ## Define File Name
  name <- paste(PATH, "history_parameters", "_init=", params_algo_EM$method.init, "_initalpha=", params_algo_EM$method.init.alpha, "_nbrofshifts=", params_algo_EM$nbr_of_shifts, ".csv", sep="")
  write.csv2(history, name, ...)
}

plot.history.OU.stationnary <- function(params_history, paramsSimu, PATH, params_algo_EM, name){
  history <- list_to_table.history(params_history)
  params_simu  <-  unlist(paramsSimu)[names(history[,1])]
  pdf(paste(PATH, "history_parameters", "_init=", params_algo_EM$method.init, "_initalpha=", params_algo_EM$method.init.alpha, "_nbrofshifts=", params_algo_EM$nbr_of_shifts, ".pdf", sep=""), width = 12, height = 8)
  ## Create grid
  pushViewport(viewport(layout = grid.layout(2+1, 3, heights = unit(c(1,5,5), "null"))))
  ## Title of the page
  grid.text(paste("Initialization : ", params_algo_EM$method.init, "\n", "Alpha Initialization : ", params_algo_EM$method.init.alpha, sep=""), vp = vplayout(1,1:3))
  ## Continuous parameters
  row <- c(2,2,3,3); col <- c(1,2,1,2)
  params_to_plot <- c("variance", "selection.strength", "root.state.var.root", "optimal.value")
  params_to_plot_legend <- c(expression(sigma^2), expression(alpha), expression(gamma^2), expression(beta[0]))
  for (s in 1:4) {
    # Choose right score
    df <- history[params_to_plot[s],]
    # Plot
    p <- qplot(seq_along(df)-1, df, xlab="iterations", ylab=params_to_plot_legend[s])
    p <- p + geom_hline(yintercept=params_simu[params_to_plot[s]])
    print(p, vp=vplayout(row[s],col[s]))
  }
  ## Plot with differents colours for each shift
  # Values
  df_val <- as.data.frame(t(history[grepl('shifts.values', names(history[,1])),, drop=F]))
  df_val_long <- melt(df_val, variable.name = "shift", value.name="shift.value")
  df_val_long$iterations <- rep(seq_along(df_val[,1])-1, dim(df_val)[2])
  df_val_long$shift <- as.factor(rep(seq_along(df_val[1,]), each=dim(df_val)[1]))
  df_val_long$true <- rep(params_simu[grepl('shifts.values', names(history[,1]))], each=dim(df_val)[1])
  p <- ggplot(data=df_val_long, aes(x=iterations, y=shift.value, colour = shift)) + geom_point()
  p <- p + geom_hline(aes(yintercept=true, colour=shift), data=df_val_long)
  p <- p + ylab(expression(delta))
  print(p, vp=vplayout(2,3))
  # Edges
  df_ed <- as.data.frame(t(history[grepl('shifts.edges', names(history[,1])),, drop=F]))
  df_ed_long <- melt(df_ed, variable.name = "shift", value.name="shift.edge")
  df_ed_long$iterations <- rep(seq_along(df_ed[,1])-1, dim(df_ed)[2])
  df_ed_long$shift <- as.factor(rep(seq_along(df_ed[1,]), each=dim(df_ed)[1]))
  df_ed_long$true <- rep(params_simu[grepl('shifts.edges', names(history[,1]))], each=dim(df_ed)[1])
  p <- ggplot(data=df_ed_long, aes(x=iterations, y=shift.edge, colour = shift)) + geom_point()
  p <- p + geom_hline(aes(yintercept=true, colour=shift), data=df_ed_long)
  p <- p + ylab(expression(tau))
  print(p, vp=vplayout(3,3))
  dev.off()
}

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
dtoc <- function(x) gsub(".", ",", x, fixed=TRUE)