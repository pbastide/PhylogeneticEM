context("Initialization Functions")

test_that("Default Init for BM", {
  ## Default default
  params <- init.EM.default.BM()
  default <- check_dimensions(1, params$root.state, params$shifts, params$variance, params$pheno_error)
  expect_that(default$shifts, equals(params$shifts))
  expect_that(default$root.state, equals(params$root.state))
  expect_that(default$variance, equals(params$variance))
  expect_that(default$pheno_error, equals(params$pheno_error))
  
  ## dim 1, user provided default
  params <- init.EM.default.BM(variance.init = 3.2,
                               pheno_error.init = 0.2,
                               random.init = TRUE,
                               value.root.init = 324.45,
                               exp.root.init = 1,
                               var.root.init = 4.5,
                               edges.init = c(12, 48),
                               values.init = c(-4, 2),
                               relativeTimes.init = 0,
                               nbr_of_shifts = 2)
  default <- check_dimensions(1, params$root.state, params$shifts, params$variance, params$pheno_error)
  expect_that(default$shifts, equals(params$shifts))
  expect_that(default$root.state, equals(params$root.state))
  expect_that(default$variance, equals(params$variance))
  expect_that(default$pheno_error, equals(params$pheno_error))
  
  ## dim p, user provided default
  p <- 4
  params <- init.EM.default.BM(variance.init = diag(3.3, p, p),
                               pheno_error = diag(0.9, p, p),
                               random.init = TRUE,
                               value.root.init = 324.45,
                               exp.root.init = 1:p,
                               var.root.init = diag(rep(4, p)),
                               edges.init = c(12, 48),
                               values.init = matrix(c(-4, 2), p, 2),
                               relativeTimes.init = 0,
                               nbr_of_shifts = 2,
                               Y_data = matrix(NA, p, 213))
  default <- check_dimensions(p, params$root.state, params$shifts, params$variance, params$pheno_error)
  expect_that(default$shifts, equals(params$shifts))
  expect_that(default$root.state, equals(params$root.state))
  expect_that(default$variance, equals(params$variance))
  expect_that(default$pheno_error, equals(params$pheno_error))
})

test_that("Default Init for OU", {
  ## Default default
  params <- init.EM.default.OU()
  default <- check_dimensions(1, params$root.state, params$shifts, params$variance,
                              params$pheno_error, params$selection.strength, params$optimal.value)
  expect_that(default$shifts, equals(params$shifts))
  expect_that(default$root.state, equals(params$root.state))
  expect_that(default$variance, equals(params$variance))
  expect_that(default$pheno_error, equals(params$pheno_error))
  expect_that(default$selection.strength, equals(params$selection.strength))
  expect_that(default$optimal.value, equals(params$optimal.value))
  
  ## dim 1, user provided default
  params <- init.EM.default.OU(variance.init = 3.2,
                               pheno_error.init = 0.8,
                               random.init = TRUE,
                               stationary.root.init = TRUE,
                               value.root.init = 324.45,
                               exp.root.init = 1,
                               var.root.init = 4.5,
                               optimal.value.init = 1,
                               selection.strength.init = 2,
                               edges.init = c(12, 48),
                               values.init = c(-4, 2),
                               relativeTimes.init = 0,
                               nbr_of_shifts = 2)
  default <- check_dimensions(1, params$root.state, params$shifts, params$variance, params$pheno_error,
                              params$selection.strength, params$optimal.value)
  expect_that(default$shifts, equals(params$shifts))
  expect_that(default$root.state, equals(params$root.state))
  expect_that(default$variance, equals(params$variance))
  expect_that(default$selection.strength, equals(params$selection.strength))
  expect_that(default$shifts, equals(params$shifts))
  expect_that(default$pheno_error, equals(params$pheno_error))
  
  ## dim p, user provided default
  p <- 4
  params <- init.EM.default.OU(variance.init = diag(3.3, p, p),
                               pheno_error = diag(0.9, p, p),
                               random.init = TRUE,
                               value.root.init = 324.45,
                               exp.root.init = 1:p,
                               var.root.init = diag(rep(4, p)),
                               optimal.value.init = 1:p,
                               selection.strength.init = diag(rep(2, p)),
                               edges.init = c(12, 48),
                               values.init = matrix(c(-4, 2), p, 2),
                               relativeTimes.init = 0,
                               nbr_of_shifts = 2,
                               Y_data = matrix(NA, p, 213))
  default <- check_dimensions(p, params$root.state, params$shifts, params$variance, params$pheno_error,
                              params$selection.strength, params$optimal.value)
  expect_that(default$shifts, equals(params$shifts))
  expect_that(default$root.state, equals(params$root.state))
  expect_that(default$variance, equals(params$variance))
  expect_that(default$selection.strength, equals(params$selection.strength))
  expect_that(default$shifts, equals(params$shifts))
  expect_that(default$pheno_error, equals(params$pheno_error))
})