####################
## Parameters
####################
library(doParallel)
library(foreach)
library(ape)
library(quadrupen) # For Lasso initialization
library(robustbase) # For robust fitting of alpha
reqpckg <- c("ape", "quadrupen", "robustbase")

## Set number of parallel cores
Ncores <- 3

## Define date-stamp for file names
datestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
datestamp_day <- format(Sys.time(), "%Y-%m-%d")

## Load simulated data
datestamp_data <- "2015-03-17" #format(Sys.time(), "%Y-%m-%d")
savedatafile = "../Results/Simulations_Several_K/several_K_simlist"
saveresultfile <- "../Results/Simulations_Several_K/several_K_alpha_0"

load(paste0(savedatafile, "_", datestamp_data, ".RData"))

## Source functions
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

## These values should be erased by further allocations (generate_inference_files)
n.range <- n
inference.index <- 0

## Select data (according to the value of n)
n <- n

## Here n.range should be defined by generate_inference_files.R
simulations2keep <- sapply(simlist, function(x) { x$n %in% n.range }, simplify = TRUE)
simlist <- simlist[simulations2keep]
nbrSim <- length(simlist)

##########################################################
## Estimation Function
##########################################################

estimations_several_K_alpha_0 <- function(X){
  ## Inference function
  fun <- function(K_t){
    return(estimation_wrapper.OUsr(K_t, 
                                   phylo = trees[[paste0(X$ntaxa)]], 
                                   Y_data = X$Y_data, 
                                   times_shared = times_shared[[paste0(X$ntaxa)]], 
                                   distances_phylo = distances_phylo[[paste0(X$ntaxa)]], 
                                   T_tree = T_tree[[paste0(X$ntaxa)]],
                                   subtree.list = subtree.list[[paste0(X$ntaxa)]],
                                   h_tree = max(diag(times_shared[[paste0(X$ntaxa)]])[1:X$ntaxa]),
                                   alpha_known = FALSE,
                                   Nbr_It_Max = 0))
  }
  ## Apply function for all K_try
  XX <- lapply(K_try[[paste0(X$ntaxa)]], fun)
  names(XX) <- K_try[[paste0(X$ntaxa)]]
  ## Formate results
  dd <- do.call(rbind, XX)
  df <- do.call(rbind, dd[ , "summary"])
  df <- as.data.frame(df)
  df$alpha  <- X$alpha
  df$gamma  <- X$gamma
  df$K <- X$K
  df$n <- X$n
  df$ntaxa <- X$ntaxa
  df$grp <- X$grp
  df$log_likelihood_true <- X$log_likelihood.true[1]
  df$difficulty <- X$difficulty
  ## Results
  X$results_summary <- df
  X$params_estim <- dd[, "params"]
  X$params_init_estim <- dd[, "params_init"]
  X$alpha_0 <- dd[, "alpha_0"]
  X$Zhat <- dd[, "Zhat"]
  X$m_Y_estim <- dd[, "m_Y_estim"]
  X$edge.quality <- dd[, "edge.quality"]
  return(X)
}

cl <- makeCluster(Ncores)
registerDoParallel(cl)

## Parallelized estimations
time_alpha_0 <- system.time(
  simestimations_alpha_0 <- foreach(i = simlist[1:10], .packages = reqpckg) %dopar%
{
  estimations_several_K_alpha_0(i)
}
)
# Stop the cluster (parallel)
stopCluster(cl)

save.image(paste0(saveresultfile, "-", datestamp_day, "_", inference.index, ".RData"))

###################################################################
## Extract Results
###################################################################
saveresultfile <- "../Results/Simulations_Several_K/several_K_alpha_0"
datestamp_day <- "2015-03-31"
inference.index <- 0

load(paste0(saveresultfile, "-", datestamp_day, "_res_", inference.index, ".RData"))

extract_data_frame <- function(simests){
  dd <- do.call(rbind, simests)
  results_summary <- as.data.frame(do.call(rbind, dd[,"results_summary"]))
  return(results_summary)
}

result_summary_alpha_0 <- extract_data_frame(simestimations_alpha_0)

####################################################################
## Fails
####################################################################
result_summary_alpha_0 <- transform(result_summary_alpha_0,
                                    fail_regression = (result_summary_alpha_0$alpha_0_regression == 1))

mean(result_summary_alpha_0$fail_regression)

results_fails <- ddply(result_summary_alpha_0,
                       .(alpha, gamma, K, K_try, ntaxa, grp),
                       summarize,
                       n.fails = mean(fail_regression))

####################################################################
## Format for Plot
####################################################################

format_plot <- function(df){
  # Save true values
  df$alpha.true <- df$alpha
  df$gamma.true <- df$gamma
  df$K.true <- df$K
  
  # Transformations and names
  df[["ln(2)/alpha"]] <- log(2)/df$alpha
  df[["sigma^2 == 2*alpha*gamma^2"]] <- 2*df$alpha*df$gamma
  
  # Melt
  df_plot <- melt(df, 
                  measure.vars = c("ln(2)/alpha", "sigma^2 == 2*alpha*gamma^2", "K"),
                  value.name = "parameter.value")
  
  # Supress lines where group != variable
  grp_var <- c(alpha_var = "ln(2)/alpha", 
               gamma_var = "sigma^2 == 2*alpha*gamma^2", 
               K_var = "K")
  fun <- function(z){
    if (z["grp"] == "base") return(TRUE)
    return(unname(grp_var[z["grp"]] == z["variable"]))
  }
  masque <- apply(df_plot, 1, fun)
  df_plot <- df_plot[masque,]
  # Factor for boxplots
  df_plot$ntaxa <- as.factor(df_plot$ntaxa)
  df_plot$parameter.value <- round(df_plot$parameter.value, 2)
  df_plot$parameter.value <- as.factor(df_plot$parameter.value)
  return(df_plot)
}

result_summary_alpha_0_plot <- format_plot(result_summary_alpha_0)
results_fails_plot <- format_plot(results_fails)

## Melt for estimations of alpha_0
result_summary_alpha_0_plot_methods <- melt(result_summary_alpha_0_plot,
                                            measure.vars = c("alpha_0_regression", "alpha_0_median", "alpha_estim_init"),
                                            value.name = "alpha_0",
                                            variable.name = "method")

#####################################################################
## Plot
#####################################################################
## Utility functions
gtable_filter <- function (x, pattern, fixed = FALSE, trim = TRUE, complementary = FALSE) 
{
  matches <- grepl(pattern, x$layout$name, fixed = fixed)
  if (complementary) matches <- !matches
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  if (trim) 
    x <- gtable_trim(x)
  x
}

g_table_put_right_strips <- function(g, grid){
  # Find x label and delete it
  xlab <- gtable_filter(g, "xlab", trim = FALSE, complementary = FALSE)
  g <- gtable_filter(g, "xlab", trim = FALSE, complementary = TRUE)
  # delete unused strips
  g <- gtable_filter(g, "strip_t-[4-9]", trim = FALSE, complementary = TRUE)
  # deplace used strips
  matches <- grepl("strip_t", g$layout$name)
  g$layout[matches, "t"] <- xlab$layout[1,"t"]
  g$layout[matches, "b"] <- xlab$layout[1,"b"]
  return(g)
}

facet_wrap_labeller <- function(g, labels = NULL) {
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))
  
  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }
  g$grobs <- gg
  return(g)
}

## modified version of geom_boxplot
geom_boxplot_noOutliers <- function (mapping = NULL, data = NULL, stat = "boxplot",
                                     position = "dodge", outlier.colour = NULL,
                                     outlier.shape = NULL, outlier.size = NULL,
                                     notch = FALSE, notchwidth = .5, varwidth = FALSE,
                                     ...) {
  
  #outlier_defaults <- ggplot2:::Geom$find('point')$default_aes()
  
  #outlier.colour   <- outlier.colour %||% outlier_defaults$colour
  #outlier.shape    <- outlier.shape  %||% outlier_defaults$shape
  #outlier.size     <- outlier.size   %||% outlier_defaults$size
  
  GeomBoxplot_noOutliers$new(mapping = mapping, data = data, stat = stat,
                             position = position, outlier.colour = outlier.colour,
                             outlier.shape = outlier.shape, outlier.size = outlier.size, notch = notch,
                             notchwidth = notchwidth, varwidth = varwidth, ...)
}

GeomBoxplot_noOutliers <- proto(ggplot2:::Geom, {
  objname <- "boxplot_noOutliers"
  
  reparameterise <- function(., df, params) {
    df$width <- df$width %||%
      params$width %||% (resolution(df$x, FALSE) * 0.9)
    
    # if (!is.null(df$outliers)) {
    #    suppressWarnings({
    #      out_min <- vapply(df$outliers, min, numeric(1))
    #      out_max <- vapply(df$outliers, max, numeric(1))
    #    })
    #    
    #    df$ymin_final <- pmin(out_min, df$ymin)
    #    df$ymax_final <- pmax(out_max, df$ymax)
    #   }
    
    # if `varwidth` not requested or not available, don't use it
    if (is.null(params) || is.null(params$varwidth) || !params$varwidth || is.null(df$relvarwidth)) {
      df$xmin <- df$x - df$width / 2
      df$xmax <- df$x + df$width / 2
    } else {
      # make `relvarwidth` relative to the size of the largest group
      df$relvarwidth <- df$relvarwidth / max(df$relvarwidth)
      df$xmin <- df$x - df$relvarwidth * df$width / 2
      df$xmax <- df$x + df$relvarwidth * df$width / 2
    }
    df$width <- NULL
    if (!is.null(df$relvarwidth)) df$relvarwidth <- NULL
    
    df
  }
  
  draw <- function(., data, ..., fatten = 2, outlier.colour = NULL, outlier.shape = NULL, outlier.size = 2,
                   notch = FALSE, notchwidth = .5, varwidth = FALSE) {
    common <- data.frame(
      colour = data$colour,
      size = data$size,
      linetype = data$linetype,
      fill = alpha(data$fill, data$alpha),
      group = data$group,
      stringsAsFactors = FALSE
    )
    
    whiskers <- data.frame(
      x = data$x,
      xend = data$x,
      y = c(data$upper, data$lower),
      yend = c(data$ymax, data$ymin),
      alpha = NA,
      common)
    
    box <- data.frame(
      xmin = data$xmin,
      xmax = data$xmax,
      ymin = data$lower,
      y = data$middle,
      ymax = data$upper,
      ynotchlower = ifelse(notch, data$notchlower, NA),
      ynotchupper = ifelse(notch, data$notchupper, NA),
      notchwidth = notchwidth,
      alpha = data$alpha,
      common)
    
    #  if (!is.null(data$outliers) && length(data$outliers[[1]] >= 1)) {
    #    outliers <- data.frame(
    #      y = data$outliers[[1]],
    #      x = data$x[1],
    #      colour = outlier.colour %||% data$colour[1],
    #      shape = outlier.shape %||% data$shape[1],
    #      size = outlier.size %||% data$size[1],
    #      fill = NA,
    #      alpha = NA,
    #      stringsAsFactors = FALSE)
    #    outliers_grob <- GeomPoint$draw(outliers, ...)
    #  } else {
    outliers_grob <- NULL
    #  }
    
    ggname(.$my_name(), grobTree(
      outliers_grob,
      GeomSegment$draw(whiskers, ...),
      GeomCrossbar$draw(box, fatten = fatten, ...)
    ))
  }
  
  guide_geom <- function(.) "boxplot_noOutliers"
  draw_legend <- function(., data, ...)  {
    data <- aesdefaults(data, .$default_aes(), list(...))
    gp <- with(data, gpar(col=colour, fill=alpha(fill, alpha), lwd=size * .pt, lty = linetype))
    gTree(gp = gp, children = gList(
      linesGrob(0.5, c(0.1, 0.25)),
      linesGrob(0.5, c(0.75, 0.9)),
      rectGrob(height=0.5, width=0.75),
      linesGrob(c(0.125, 0.875), 0.5)
    ))
  }
  
  default_stat <- function(.) StatBoxplot
  default_pos <- function(.) PositionDodge
  default_aes <- function(.) aes(weight=1, colour="grey20", fill="white", size=0.5, alpha = NA, shape = 16, linetype = "solid")
  required_aes <- c("x", "lower", "upper", "middle", "ymin", "ymax")
  
})

##################################################################
## Plots
##################################################################
plot_alpha_0_boxplot <- function(estim, true = NULL, name, sub = "K_try", val){
  replace_name <- function(i){
    i <- as.integer(as.character(i))
    i <- formatC(i, width = 3, format = "d", flag = "0")
    i <- paste0("ntaxa == ", i)
    return(i)
  }
  data <- subset(result_summary_alpha_0_plot_methods,
                 result_summary_alpha_0_plot_methods[, sub] == val)
  data$ntaxa <- sapply(data$ntaxa, replace_name)
  data$ntaxa <- as.factor(data$ntaxa)
  data$parameter.value <- data$parameter.value
  if (!is.null(true)){
    minustrueontrue <- paste0("/", true, "- 1")
  } else {
    minustrue <- NULL
  }
  p <- ggplot(data,
              aes_string(x = "parameter.value",
                         y = paste0(estim, minustrueontrue),
                         color = "method"))
  p <- p + scale_colour_discrete(name = "Method",
                                 breaks = c("alpha_0_regression", "alpha_0_median", "alpha_estim_init"),
                                 labels = c("Robust Regression", "Median", "Mean of the two"))
  #p <- p + aes(group = interaction(ntaxa, alpha_known))
  p <- p + facet_wrap(ntaxa ~ variable,
                      #labeller = label_parsed,
                      scales = "free",
                      shrink = TRUE)
  #   p <- p + geom_point(alpha = 0.2,
  #                       position = position_jitter(width = 0.4, height = 0))
  p <- p + geom_boxplot_noOutliers(outlier.size = 0.5, alpha = 0.5)
  if (!is.null(true)) p <- p + geom_hline(yintercept = 0, alpha = 0.5)
  p <- p + labs(x = "",
                y = name)
  #p <- p + scale_x_continuous(trans = log_trans(10))
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
  p <- p + theme(axis.text = element_text(size = 12),
                 strip.text = element_text(size = 12),
                 strip.background = element_blank()
                 ##legend.position = c(0, 1),
                 ##legend.justification = c(0, 1)
  )
  g <- ggplotGrob(p)
  g <- g_table_put_right_strips(g)
  labels <- unname(sapply(levels(result_summary_alpha_0_plot_methods$variable),
                          function(z) parse(text = z)))
  g <- facet_wrap_labeller(g, labels = labels)
  grid.newpage()
  grid.draw(g)
}

plot_alpha_0_boxplot("alpha_0", "alpha.true", "alpha - alpha / alpha", "K_try", 0)


results_fails_plot$K_try <- as.factor(results_fails_plot$K_try)
p <- ggplot(results_fails_plot,
            aes(x = parameter.value,
                y = n.fails,
                color = K_try))
p <- p + aes(group = K_try)
p <- p + facet_wrap(ntaxa ~ variable, scales = "free")
p <- p + geom_line()
p <- p + labs(x = "",
              y = "Proportion of regression fails")
#p <- p + scale_x_continuous(trans = log_trans(10))
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
p <- p + theme(axis.text = element_text(size = 12),
               strip.text = element_text(size = 12),
               strip.background = element_blank()
               ##legend.position = c(0, 1),
               ##legend.justification = c(0, 1)
)
g <- ggplotGrob(p)
g <- g_table_put_right_strips(g)
labels <- unname(sapply(levels(results_fails_plot$variable),
                        function(z) parse(text = z)))
g <- facet_wrap_labeller(g, labels = labels)
grid.newpage()
grid.draw(g)
