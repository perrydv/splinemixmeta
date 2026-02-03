test_that("make_smm_smooth works (simple, univariate)", {
  sTerm <- mgcv::s(x, k = 5)
  x <- rnorm(30)
  ssm_smooth <- make_smm_smooth(sTerm, data = data.frame(x = x, n = length(x)), vnames = "x")
  expect_equal(length(ssm_smooth), 1)
  ssm_smooth <- ssm_smooth[[1]]
  expect_equal(dim(ssm_smooth$basisFxns), c(30, 3))
  expect_equal(dim(ssm_smooth$x_fixed), c(30, 1))
  # last arg FALSE
  ssm_smooth <- make_smm_smooth(sTerm, data = data.frame(x = x, n = length(x)), vnames = "x")
  expect_equal(length(ssm_smooth), 1)
  ssm_smooth <- ssm_smooth[[1]]
  expect_equal(dim(ssm_smooth$basisFxns), c(30, 3))
  # expect_true(is.null(ssm_smooth$x_fixed))
  # omit data
  ssm_smooth2 <- make_smm_smooth(sTerm, vnames = "x")
  expect_equal(length(ssm_smooth2), 1)
  ssm_smooth2 <- ssm_smooth2[[1]]
  expect_identical(ssm_smooth2, ssm_smooth)

  sTerm <- mgcv::s(x, k = 5, bs = "cr")
  ssm_smooth <- make_smm_smooth(sTerm, data = data.frame(x = x, n = length(x)), vnames = "x")
  expect_equal(length(ssm_smooth), 1)
  ssm_smooth <- ssm_smooth[[1]]
  expect_equal(dim(ssm_smooth$basisFxns), c(30, 3))
  expect_equal(dim(ssm_smooth$x_fixed), c(30, 1))
  ssm_smooth <- make_smm_smooth(sTerm, data = data.frame(x = x, n = length(x)), vnames = "x")
  expect_equal(length(ssm_smooth), 1)
  ssm_smooth <- ssm_smooth[[1]]
  expect_equal(dim(ssm_smooth$basisFxns), c(30, 3))
  # expect_true(is.null(ssm_smooth$x_fixed))

  # cyclic spline has no unpenalized part
  sTerm <- mgcv::s(x, k = 5, bs = "cc")
  ssm_smooth <- make_smm_smooth(sTerm, data = data.frame(x = x, n = length(x)), vnames = "x")
  expect_equal(length(ssm_smooth), 1)
  ssm_smooth <- ssm_smooth[[1]]
  expect_equal(dim(ssm_smooth$basisFxns), c(30, 3))
  expect_equal(dim(ssm_smooth$x_fixed), c(30, 0))
  ssm_smooth <- make_smm_smooth(sTerm, data = data.frame(x = x, n = length(x)), vnames = "x")
  expect_equal(length(ssm_smooth), 1)
  ssm_smooth <- ssm_smooth[[1]]
  expect_equal(dim(ssm_smooth$basisFxns), c(30, 3))
  # expect_true(is.null(ssm_smooth$x_fixed))
})

# Splitting by group should be done in mixmeta, not the s() term.
#
# test_that("make_smm_smooth works (univariate by factor)", {
#   sTerm <- mgcv::s(x, k = 5, by = group)
#   x <- rnorm(30)
#   group <- factor(rep(letters[1:3], each = 10))
#   data <- data.frame(x = x, n = length(x), group = group)
#   ssm_smooth <- make_smm_smooth(sTerm, data = data, vnames = "x", auto_fixed_effects = TRUE)
#   expect_equal(dim(ssm_smooth$basisFxns), c(30, 3))
#   expect_equal(dim(ssm_smooth$x_fixed), c(30, 1))
#
#   sTerm <- mgcv::s(x, k = 5, bs = "cr")
#   ssm_smooth <- make_smm_smooth(sTerm, data = data.frame(x = x, n = length(x)), vnames = "x", auto_fixed_effects = TRUE)
#   expect_equal(dim(ssm_smooth$basisFxns), c(30, 3))
#   expect_equal(dim(ssm_smooth$x_fixed), c(30, 1))
#
#   sTerm <- mgcv::s(x, k = 5, bs = "cc")
#   ssm_smooth <- make_smm_smooth(sTerm, data = data.frame(x = x, n = length(x)), vnames = "x", auto_fixed_effects = TRUE)
#   expect_equal(dim(ssm_smooth$basisFxns), c(30, 3))
#   expect_equal(dim(ssm_smooth$x_fixed), c(30, 0))
# })

test_that("make_smm_smooth works (simple, univariate)", {
  # For "tp", k is a single number across all dimensions
  set.seed(1)
  x <- rnorm(50)
  y <- rnorm(50)
  vnames <- c("x", "y")
  sTerm <- mgcv::s(x, y, k = 20, bs = "tp")
  data <- list(x = x, y = y, n = 50)
  ## error because n should not be included with multiple variables
  expect_error(ssm_smooth <- make_smm_smooth(sTerm, data = data, vnames = vnames))

  data <- data.frame(x=x, y = y)

  ssm_smooth <- make_smm_smooth(sTerm, data = data, vnames = vnames)
  expect_equal(length(ssm_smooth), 1)
  expect_equal(dim(ssm_smooth[[1]]$basisFxns), c(50, 17))
  expect_equal(dim(ssm_smooth[[1]]$x_fixed), c(50, 2))

  ssm_smooth <- make_smm_smooth(sTerm, vnames = vnames)
  expect_equal(length(ssm_smooth), 1)
  expect_equal(dim(ssm_smooth[[1]]$basisFxns), c(50, 17))
  expect_equal(dim(ssm_smooth[[1]]$x_fixed), c(50, 2))

  # mrf would take more work to set up
  # sTerm <- mgcv::s(x, y, k = -1, bs = "mrf")
  # ssm_smooth <- make_smm_smooth(sTerm, data = data, vnames = vnames, auto_fixed_effects = TRUE)

  # We can't support mgcv::te.
  # It evidently can't give diagonal penalties, but it does
  # give multiple random effects,
  # and we can't handle multiple bscov="prop" random effects
  # sTerm <- mgcv::te(x, y, bs = "cr")
  #test<-mgcv::smoothCon(sTerm, data = data,
  #                absorb.cons = TRUE,
  #                diagonal.penalty = TRUE)
  # ssm_smooth <- make_smm_smooth(sTerm, data = data, vnames = vnames, auto_fixed_effects = TRUE)

  # We also can't handle ti
  # sTerm <- mgcv::ti(x, y, bs = "cr")
  # test<-mgcv::smoothCon(sTerm, data = data,
  #                absorb.cons = TRUE,
  #                diagonal.penalty = TRUE)

})
