context("Function merge rotation")

###############################################################################
test_that("Rotations", {
  testthat::skip_on_cran()
  data(monkeys)
  
  # Fit data
  res <- PhyloEM(Y_data = monkeys$dat,
                 phylo = monkeys$phy, 
                 process = "scOU", nbr_alpha = 2,
                 random.root = FALSE,
                 K_max = 10,
                 method.selection = "BGHmlraw",
                 parallel_alpha = F)
  
  # rotate data
  rot <- matrix(c(cos(pi/4), -sin(pi/4), sin(pi/4), cos(pi/4)),
                nrow= 2, ncol = 2)
  Yrot <- t(rot) %*% monkeys$dat
  
  # Fit rot
  res_rot <- PhyloEM(Y_data = Yrot,
                     phylo = monkeys$phy, 
                     process = "scOU", nbr_alpha = 2,
                     random.root = FALSE,
                     K_max = 10,
                     method.selection = "BGHmlraw",
                     parallel_alpha = F)
  
  # Find rotation
  expect_equal(find_rotation(res, res_rot), rot)
  
  # rotate params
  expect_equivalent(res$alpha_max$params_estim$`3`, rotate_params(res_rot$alpha_max$params_estim$`3`, rot))
  
  expect_equal(log_likelihood(params_process(res, K = 3), phylo = monkeys$phy, Y_data = monkeys$dat),
               log_likelihood(rotate_params(params_process(res_rot, K = 3), rot), phylo = monkeys$phy, Y_data = monkeys$dat))
  
  # another rotation
  rot2 <- matrix(c(cos(pi/3), -sin(pi/3), sin(pi/3), cos(pi/3)),
                 nrow= 2, ncol = 2)
  Yrot2 <- t(rot2) %*% monkeys$dat
  
  # Fit rot
  res_rot2 <- PhyloEM(Y_data = Yrot2,
                      phylo = monkeys$phy, 
                      process = "scOU", nbr_alpha = 2,
                      random.root = FALSE,
                      K_max = 10,
                      method.selection = "BGHmlraw",
                      parallel_alpha = F)
  
  # merge solution
  res_merge <- merge_rotations(res, res_rot, res_rot2)
  
  expect_equal(res_merge$alpha_max$results_summary$log_likelihood,
               apply(sapply(list(res, res_rot, res_rot2), function(x) x$alpha_max$results_summary[["log_likelihood"]]), 1, max))
  
  expect_equal(log_likelihood(params_process(res_merge, K = 3), phylo = monkeys$phy, Y_data = monkeys$dat),
               res_merge$alpha_max$results_summary$log_likelihood[4])

  ###############################################################################
  ## Not rotations
  
  ## transform data
  norot <- matrix(c(1, 0.5, 0.1, 2), 2)
  Ynorot <- t(norot) %*% monkeys$dat
  
  # Fit rot
  expect_warning(res_norot <- PhyloEM(Y_data = Ynorot,
                                      phylo = monkeys$phy, 
                                      process = "scOU", nbr_alpha = 2,
                                      random.root = FALSE,
                                      K_max = 10,
                                      method.selection = "BGHmlraw",
                                      parallel_alpha = F),
                 "The solution fund by lasso had non independent vectors. Had to modify this solution.")
  
  # Rotation ?
  expect_error(find_rotation(res, res_norot), "The datasets are not linked by a rotation.")
  
  ## Other completely unrelated data
  Ynorot2 <- matrix(rnorm(length(monkeys$dat)), dim(monkeys$dat))
  colnames(Ynorot2) <- colnames(Ynorot)
  
  # Fit rot
  res_norot2 <- PhyloEM(Y_data = Ynorot2,
                      phylo = monkeys$phy, 
                      process = "scOU", nbr_alpha = 2,
                      random.root = FALSE,
                      K_max = 10,
                      method.selection = "BGHmlraw",
                      parallel_alpha = F)
  
  # Rotation ?
  expect_error(find_rotation(res, res_norot2), "The datasets are not linearly mapped.")
  
  ###############################################################################
  ## Missing data
  
  dat <- monkeys$dat
  dat[, c(12, 23, 32)] <- NA
  
  # Fit data
  res <- PhyloEM(Y_data = dat,
                 phylo = monkeys$phy, 
                 process = "scOU", nbr_alpha = 2,
                 random.root = FALSE,
                 K_max = 10,
                 method.selection = "BGHmlraw",
                 parallel_alpha = F)
  
  # rotate data
  rot <- matrix(c(cos(pi/4), -sin(pi/4), sin(pi/4), cos(pi/4)),
                nrow= 2, ncol = 2)
  Yrot <- t(rot) %*% dat
  
  # Fit rot
  res_rot <- PhyloEM(Y_data = Yrot,
                     phylo = monkeys$phy, 
                     process = "scOU", nbr_alpha = 2,
                     random.root = FALSE,
                     K_max = 10,
                     method.selection = "BGHmlraw",
                     parallel_alpha = F)
  
  # Find rotation
  expect_equal(find_rotation(res, res_rot), rot)
  
  # rotate params
  expect_equivalent(res$alpha_max$params_estim$`3`,
                    rotate_params(res_rot$alpha_max$params_estim$`3`, rot),
                    1e-5)
  
  expect_equal(log_likelihood(params_process(res, K = 3), phylo = monkeys$phy, Y_data = monkeys$dat),
               log_likelihood(rotate_params(params_process(res_rot, K = 3), rot), phylo = monkeys$phy, Y_data = monkeys$dat))
  
  
  # another rotation
  rot2 <- matrix(c(cos(pi/3), -sin(pi/3), sin(pi/3), cos(pi/3)),
                 nrow= 2, ncol = 2)
  Yrot2 <- t(rot2) %*% dat
  
  # Fit rot
  res_rot2 <- PhyloEM(Y_data = Yrot2,
                      phylo = monkeys$phy, 
                      process = "scOU", nbr_alpha = 2,
                      random.root = FALSE,
                      K_max = 10,
                      method.selection = "BGHmlraw",
                      parallel_alpha = F)
  
  # merge solution
  res_merge <- merge_rotations(res, res_rot, res_rot2)
  
  expect_equal(res_merge$alpha_max$results_summary$log_likelihood,
               apply(sapply(list(res, res_rot, res_rot2), function(x) x$alpha_max$results_summary[["log_likelihood"]]), 1, max))
  
  expect_equal(log_likelihood(params_process(res_merge, K = 3), phylo = monkeys$phy, Y_data = dat),
               res_merge$alpha_max$results_summary$log_likelihood[4])
  
  
  ###############################################################################
  ## Missing data - Errors
  
  dat <- monkeys$dat
  dat[, c(12, 23, 32)] <- NA
  dat[1, 13] <- NA
  
  # Fit data
  res <- PhyloEM(Y_data = dat,
                 phylo = monkeys$phy, 
                 process = "scOU", nbr_alpha = 2,
                 random.root = FALSE,
                 K_max = 10,
                 method.selection = "BGHmlraw",
                 parallel_alpha = F)
  
  # rotate data
  rot <- matrix(c(cos(pi/4), -sin(pi/4), sin(pi/4), cos(pi/4)),
                nrow= 2, ncol = 2)
  Yrot <- t(rot) %*% dat
  
  # Fit rot
  res_rot <- PhyloEM(Y_data = Yrot,
                     phylo = monkeys$phy, 
                     process = "scOU", nbr_alpha = 2,
                     random.root = FALSE,
                     K_max = 10,
                     method.selection = "BGHmlraw",
                     parallel_alpha = F)
  
  # Find rotation
  expect_error(find_rotation(res, res_rot),
               "Rotations can only be applied to datasets that have entire species missing")
  expect_error(find_rotation(res_rot, res),
               "Rotations can only be applied to datasets that have entire species missing")
  
  # all missing plus another
  dat[2, 13] <- NA
  dat[, 14] <- NA
  
  res <- PhyloEM(Y_data = dat,
                 phylo = monkeys$phy, 
                 process = "scOU", nbr_alpha = 2,
                 random.root = FALSE,
                 K_max = 10,
                 method.selection = "BGHmlraw",
                 parallel_alpha = F)
  
  expect_error(find_rotation(res_rot, res),
               "The two datasets used in the analyses do not have the same missing data.")
  
})