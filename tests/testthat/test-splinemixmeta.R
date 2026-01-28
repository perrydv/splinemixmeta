
test_that("splinemixmeta perfectly passes through to mixmeta if no smooth is specified", {
  sim_data <- sim_basic_meta()

  res_splinemix <- splinemixmeta(smooth = NULL, y ~ x, random = ~ 1 | group, data = sim_data, S = sim_data$S)
  res_mixmeta <- mixmeta::mixmeta(y ~ x, random = ~ 1 | group, data = sim_data, S = sim_data$S)
  expect_identical(res_splinemix$coefficients, res_mixmeta$coefficients)
  expect_identical(res_splinemix$Psi, res_mixmeta$Psi)

  # Without data argument, as if the columns of sim_data are in the environment
  attach(sim_data)
  on.exit(detach(sim_data))
  res_splinemix <- splinemixmeta(smooth = NULL, y ~ x, random = ~ 1 | group, S = sim_data$S)
  res_mixmeta <- mixmeta::mixmeta(y ~ x, random = ~ 1 | group, S = sim_data$S)
  expect_identical(res_splinemix$coefficients, res_mixmeta$coefficients)
  expect_identical(res_splinemix$Psi, res_mixmeta$Psi)
})

test_that("splinemixmeta works with a simple smooth in one x", {
  sim_data <- sim_basic_meta()
#debug(mixmeta::mixmeta)
  res_splinemixA <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y, data = sim_data, S = sim_data$S)
  res_splinemixB <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y ~ 1, data = sim_data, S = sim_data$S)
  res_splinemixC <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y ~ x, manual_fixed=TRUE, data = sim_data, S = sim_data$S)
  debug(splinemixmeta)
  res_splinemixD <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y ~ x, manual_fixed=TRUE, data = sim_data, se = sqrt(sim_data$S))
  expect_equal(AIC(res_splinemixA), AIC(res_splinemixB))
  expect_equal(AIC(res_splinemixA), AIC(res_splinemixC))
  expect_equal(AIC(res_splinemixA), AIC(res_splinemixD))
  predA <- predict_smm(res_splinemixA, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  predB <- predict_smm(res_splinemixB, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  predC <- predict_smm(res_splinemixC, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  predD <- predict_smm(res_splinemixD, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  expect_equal(predA[,'blup'], predB[,'blup'])
  expect_equal(predA[,'blup'], predC[,'blup'])
  expect_equal(predA[,'blup'], predD[,'blup'])
  expect_equal(predA[,'se'], predB[,'se'])
  expect_equal(predA[,'se'], predC[,'se'])
  expect_equal(predA[,'se'], predD[,'se'])

# Add group random effects
  res_splinemixA <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y, random = ~ 1 | group, data = sim_data, S = sim_data$S)
  res_splinemixB <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y ~ 1, random = ~ 1 | group, data = sim_data, S = sim_data$S)
  res_splinemixC <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y ~ x, manual_fixed=TRUE, random = ~ 1 | group, data = sim_data, S = sim_data$S)
  expect_equal(AIC(res_splinemixA), AIC(res_splinemixB))
  expect_equal(AIC(res_splinemixA), AIC(res_splinemixC))
  predA <- predict_smm(res_splinemixA, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  predB <- predict_smm(res_splinemixB, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  predC <- predict_smm(res_splinemixC, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  expect_equal(predA[,'blup'], predB[,'blup'])
  expect_equal(predA[,'blup'], predC[,'blup'])
  expect_equal(predA[,'se'], predB[,'se'])
  expect_equal(predA[,'se'], predC[,'se'])
})

# Original prototype that led to splinemixmeta
mixmetagam <- function(y, x, se, sTerm = s(x, bs = 'cr'),
                       mixed = TRUE, control=list()) {
  require(mixmeta)
  require(mgcv)
  len <- length(x)
  xf <- factor(x)
  x <- as.numeric(x)
  x_input <- x
  defaultControl <- list(
    smooth2random=TRUE,
    manualFixed=TRUE
  )
  for(ce in names(defaultControl)) {
    if(is.null(control[[ce]])) control[[ce]] <- defaultControl[[ce]]
  }

  # Make a factor to group the entire data vector
  all <- as.factor(rep("all", len))
  d <- data.frame(y = y, x = x, xf = xf, all = all)
  # Construct the pieces of a smooth term.
  sCdiag <- smoothCon(sTerm, data = d,
                      absorb.cons = TRUE,
                      diagonal.penalty = TRUE)

  # Obtain the basis function evaluations and spline penalty matrix.
  if(isTRUE(control[["smooth2random"]])) {
    RE <- smooth2random(sCdiag[[1]], c("x", "xf", "y", "all"), type = 2)
    #  inds_penalized <- which(RE$pen.ind != 0)
    #  inds_unpenalized <- which(RE$pen.ind == 0)
    basisFxns <- RE$rand$Xr
    Psi <- diag(ncol(basisFxns))
    if(!isTRUE(control[["manualFixed"]])) {
      x <- RE$Xf
    }
  } else {
    Smgcv <- sCdiag[[1]]$S[[1]]
    iUnpenalized <- which(diag(Smgcv)==0)
    numUnpenalized <- length(iUnpenalized)
    sDim <- dim(Smgcv)[1]
    sDim_fullrank <- sDim-numUnpenalized
    S_fullRank <- Smgcv[-iUnpenalized, -iUnpenalized]
    # Remove the unpenalized basis dimension from mgcv
    # and invert the precision matrix to be a covariance matrix
    # (although both are Identity because diagonal.penalty = TRUE).
    Psi <- solve(S_fullRank)
    # Obtain basis function evaluations at the x values and remove
    # the linear basis function, which is the last one.
    basisFxns <- PredictMat(sCdiag[[1]], data = d)[,-iUnpenalized]
  }
  # Fit the model in mixmeta.
  S <- se^2
  if(mixed) {
    random <- list(~ basisFxns - 1|all, ~1 | xf)
    bscov <- c("prop", "unstr")
  } else {
    random = list(~ basisFxns - 1|all)
    bscov <- "prop"
  }
  m1 <- mixmeta(y ~ x, data = d,  S = S,
                random = random,
                bscov = bscov,
                control = list(Psifix = Psi))

  if(mixed) {
    # Move xf random effects variance into S and refit (equivalent model)
    # so that BLUPs will include only the spline as random effects
    xf_sd <- m1$par['xf']
    S2 <- S + xf_sd*xf_sd
    m2 <- mixmeta(y ~ x, data = d, S = S2,
                  random = list(~ basisFxns - 1|all),
                  bscov = "prop",
                  control = list(Psifix = Psi))
    pred <- blup(m2, se = TRUE)
    pred_spline <- blup(m2, level = 1, se = TRUE)
  } else {
    m2 <- NULL
    pred <- blup(m1, se = TRUE)
    pred_spline <- blup(m1, level = 1, se = TRUE)
  }
  # Return list of various components of work done above.
  list(fit = m1,
       pred = pred,
       pred_spline = pred_spline,
       Psi = Psi,
       basisFxns = basisFxns,
       fit_for_blups = m2)
}


test_that("splinemixmeta gives correct answer for SFB30AO data", {
  SFB30A0 <- read.csv(file.path("fixtures", "SFB30AO.csv"))
  orig_res <- mixmetagam(y = SFB30A0$y,
                           x = SFB30A0$year,
                           se = SFB30A0$se,
                           sTerm = mgcv::s(x, k = 30, bs = "cr"),
                           mixed = TRUE)
  # S <- SFB30A0$se^2
  package_resA <- splinemixmeta(smooth = mgcv::s(year, k = 30, bs = "cr"),
                               y,
                               data = SFB30A0,
                               se = SFB30A0$se)
  package_resB <- splinemixmeta(smooth = mgcv::s(year, k = 30, bs = "cr"),
                               y ~ year,
                               manual_fixed = TRUE,
                               data = SFB30A0,
                               se = SFB30A0$se)
  expect_equal(AIC(orig_res$fit),
               AIC(package_resA))
  expect_equal(AIC(orig_res$fit),
               AIC(package_resB))
  expect_equal(as.numeric(orig_res$fit$coefficients),
               as.numeric(package_resB$coefficients))
  pred_orig <- orig_res$pred_spline
  pred_resA <- predict_smm(package_resA)
  pred_resB <- predict_smm(package_resB)

  expect_equal(pred_orig[,'blup'], pred_resA[,'blup'], tolerance= 0.01)
  expect_equal(pred_orig[,'blup'], pred_resB[,'blup'], tolerance = 0.01)
  expect_equal(pred_orig[,'se'], pred_resA[,'se'])
  expect_equal(pred_orig[,'se'], pred_resB[,'se'])
  })

