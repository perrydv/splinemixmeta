
test_that("splinemixmeta perfectly passes through to mixmeta if no smooth is specified", {
  sim_data <- sim_basic_meta()

  res_splinemix <- splinemixmeta(smooth = NULL, y ~ x, random = ~ 1 | group, data = sim_data, S = sim_data$S)
  expect_false(inherits(res_splinemix, "splinemixmeta"))
  res_mixmeta <- mixmeta::mixmeta(y ~ x, random = ~ 1 | group, data = sim_data, S = sim_data$S)
  expect_identical(res_splinemix$coefficients, res_mixmeta$coefficients)
  expect_identical(res_splinemix$Psi, res_mixmeta$Psi)

  # With se used as the way to provide S
  res_splinemix <- splinemixmeta(smooth = NULL, y ~ x, random = ~ 1 | group, data = sim_data, se = sqrt(sim_data$S))
  expect_identical(res_splinemix$coefficients, res_mixmeta$coefficients)
  expect_identical(res_splinemix$Psi, res_mixmeta$Psi)

  # Without data argument, as if the columns of sim_data are in the environment
  attach(sim_data)
  on.exit(detach(sim_data))
  res_splinemix <- splinemixmeta(smooth = NULL, y ~ x, random = ~ 1 | group, S = sim_data$S)
  res_mixmeta <- mixmeta::mixmeta(y ~ x, random = ~ 1 | group, S = sim_data$S)
  expect_identical(res_splinemix$coefficients, res_mixmeta$coefficients)
  expect_identical(res_splinemix$Psi, res_mixmeta$Psi)
  detach(sim_data)
  on.exit()
})

test_that("splinemixmeta works with a simple smooth in one x", {
  sim_data <- sim_basic_meta()
#debug(mixmeta::mixmeta)
  # y provided as a name, not a formula
  res_splinemixA <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y, data = sim_data, S = sim_data$S)
  # y ~ 1 provided as formula
  res_splinemixB <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y ~ 1, data = sim_data, S = sim_data$S)
  # provide linear x manually
  res_splinemixC <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y ~ x, manual_fixed=TRUE, data = sim_data, S = sim_data$S)
  # provide se instead of S
  res_splinemixD <- splinemixmeta(smooth = mgcv::s(x, k = 7, bs ="cr"), y ~ x, manual_fixed=TRUE, data = sim_data, se = sqrt(sim_data$S))

  expect_equal(AIC(res_splinemixA), AIC(res_splinemixB))
  expect_equal(AIC(res_splinemixA), AIC(res_splinemixC))
  expect_equal(AIC(res_splinemixA), AIC(res_splinemixD))
  predA <- predict(res_splinemixA, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  predB <- predict(res_splinemixB, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  predC <- predict(res_splinemixC, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  predD <- predict(res_splinemixD, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
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
  predA <- predict(res_splinemixA, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  predB <- predict(res_splinemixB, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  predC <- predict(res_splinemixC, include_smooths = TRUE, include_REs = FALSE, include_residuals = TRUE)
  expect_equal(predA[,'blup'], predB[,'blup'])
  expect_equal(predA[,'blup'], predC[,'blup'])
  expect_equal(predA[,'se'], predB[,'se'])
  expect_equal(predA[,'se'], predC[,'se'])
})

# test_that("splinemixmeta works with two smooths", {
#   set.seed(1)
  
#   sim_data <- sim_data_for_two_splines()

#   smm1 <- splinemixmeta(smooth = list(mgcv::s(x, k = 7, bs ="cr"),
#                                       mgcv::s(y, k = 7, bs ="cr")),
#                         z ~ 1,
#                         data = sim_data,
#                         S = sim_data$S)

#   expect_true(inherits(smm1, "splinemixmeta"))
#   pred <- predict(smm1)
#   res <- sim_data$z - pred[,'blup']

#   pred2 <- predict(smm1, type = "residual")
#   predfix <- predict(smm1, include_smooths=FALSE)
#   res2 <- sim_data$z - (predfix[,'blup'] + pred2[,'blup'])
#   expect_equal(res, res2)

#   expect_true(sd(res) < 1)
#   lmcheck <- lm(sim_data$z ~ I(pred[,'blup']))
#   expect_true(abs(coef(lmcheck)[2] - 1.0) < 0.03)
# })

# test_that("splinemixmeta works with bs='tp'", {
#   set.seed(1)
#   sim_data_for_two_splines
#   sim_data <- sim_data_for_two_splines()

#   smm1 <- splinemixmeta(smooth = mgcv::s(x, y, bs ="tp"),
#                         z ~ 1,
#                         data = sim_data,
#                         S = sim_data$S)

#   expect_true(inherits(smm1, "splinemixmeta"))
#   pred <- predict(smm1)
#   res <- sim_data$z - pred[,'blup']

#   pred2 <- predict(smm1, type = "residual")
#   predfix <- predict(smm1, include_smooths=FALSE)
#   res2 <- sim_data$z - (predfix[,'blup'] + pred2[,'blup'])
#   expect_equal(res, res2)

#   expect_true(sd(res) < 1)
#   lmcheck <- lm(sim_data$z ~ I(pred[,'blup']))
#   expect_true(abs(coef(lmcheck)[2] - 1.0) < 0.03)

# })

# # Original prototype that led to splinemixmeta
# # The test using this below is currently commented out
# # while we resolve a subtlety about blup calculations.
# mixmetagam <- function(y, x, se, sTerm = s(x, bs = 'cr'),
#                        mixed = TRUE, control=list()) {
#   require(mixmeta)
#   require(mgcv)
#   len <- length(x)
#   xf <- factor(x)
#   x <- as.numeric(x)
#   x_input <- x
#   defaultControl <- list(
#     smooth2random=TRUE,
#     manualFixed=TRUE
#   )
#   for(ce in names(defaultControl)) {
#     if(is.null(control[[ce]])) control[[ce]] <- defaultControl[[ce]]
#   }

#   # Make a factor to group the entire data vector
#   all <- as.factor(rep("all", len))
#   d <- data.frame(y = y, x = x, xf = xf, all = all)
#   # Construct the pieces of a smooth term.
#   sCdiag <- smoothCon(sTerm, data = d,
#                       absorb.cons = TRUE,
#                       diagonal.penalty = TRUE)

#   # Obtain the basis function evaluations and spline penalty matrix.
#   if(isTRUE(control[["smooth2random"]])) {
#     RE <- smooth2random(sCdiag[[1]], c("x", "xf", "y", "all"), type = 2)
#     #  inds_penalized <- which(RE$pen.ind != 0)
#     #  inds_unpenalized <- which(RE$pen.ind == 0)
#     basisFxns <- RE$rand$Xr
#     Psi <- diag(ncol(basisFxns))
#     if(!isTRUE(control[["manualFixed"]])) {
#       x <- RE$Xf
#     }
#   } else {
#     Smgcv <- sCdiag[[1]]$S[[1]]
#     iUnpenalized <- which(diag(Smgcv)==0)
#     numUnpenalized <- length(iUnpenalized)
#     sDim <- dim(Smgcv)[1]
#     sDim_fullrank <- sDim-numUnpenalized
#     S_fullRank <- Smgcv[-iUnpenalized, -iUnpenalized]
#     # Remove the unpenalized basis dimension from mgcv
#     # and invert the precision matrix to be a covariance matrix
#     # (although both are Identity because diagonal.penalty = TRUE).
#     Psi <- solve(S_fullRank)
#     # Obtain basis function evaluations at the x values and remove
#     # the linear basis function, which is the last one.
#     basisFxns <- PredictMat(sCdiag[[1]], data = d)[,-iUnpenalized]
#   }
#   # Fit the model in mixmeta.
#   S <- se^2
#   if(mixed) {
#     random <- list(~ basisFxns - 1|all, ~1 | xf)
#     bscov <- c("prop", "unstr")
#   } else {
#     random = list(~ basisFxns - 1|all)
#     bscov <- "prop"
#   }
#   m1 <- mixmeta(y ~ x, data = d,  S = S,
#                 random = random,
#                 bscov = bscov,
#                 control = list(Psifix = Psi))

#   if(mixed) {
#     # Move xf random effects variance into S and refit (equivalent model)
#     # so that BLUPs will include only the spline as random effects
#     xf_sd <- m1$par['xf']
#     S2 <- S + xf_sd*xf_sd
#     m2 <- mixmeta(y ~ x, data = d, S = S2,
#                   random = list(~ basisFxns - 1|all),
#                   bscov = "prop",
#                   control = list(Psifix = Psi))
#     pred <- blup(m2, se = TRUE)
#     pred_spline <- blup(m2, level = 1, se = TRUE)
#   } else {
#     m2 <- NULL
#     pred <- blup(m1, se = TRUE)
#     pred_spline <- blup(m1, level = 1, se = TRUE)
#   }
#   # Return list of various components of work done above.
#   list(fit = m1,
#        pred = pred,
#        pred_spline = pred_spline,
#        Psi = Psi,
#        basisFxns = basisFxns,
#        fit_for_blups = m2)
# }

# test_that("splinemixmeta gives correct answer for SFB30AO data", {
#   SFB30A0 <- read.csv(file.path("fixtures", "SFB30AO.csv"))
#   # SFB30A0$year <- as.numeric(SFB30A0$year)
#   orig_res <- mixmetagam(y = SFB30A0$y,
#                            x = SFB30A0$year,
#                            se = SFB30A0$se,
#                            sTerm = mgcv::s(x, k = 30, bs = "cr"),

#                            mixed = TRUE)
#   # S <- SFB30A0$se^2
#   package_resA <- splinemixmeta(smooth = mgcv::s(year, k = 30, bs = "cr"),
#                                y,
#                                data = SFB30A0,
#                                se = SFB30A0$se)
#   package_resB <- splinemixmeta(smooth = mgcv::s(year, k = 30, bs = "cr"),
#                                y ~ year,
#                                manual_fixed = TRUE,
#                                data = SFB30A0,
#                                se = SFB30A0$se)
#   expect_equal(AIC(orig_res$fit),
#                AIC(package_resA))
#   expect_equal(AIC(orig_res$fit),
#                AIC(package_resB))
#   expect_equal(as.numeric(orig_res$fit$coefficients),
#                as.numeric(package_resB$coefficients))
#   pred_orig <- orig_res$pred_spline
#   pred_resA <- predict(package_resA)
#   pred_resB <- predict(package_resB)

#   ###
#   # Here I will try to step through blup.mixmeta, since we are getting different blup results, especially in the std err's
#   # blups including spline but not residuals RE by setting level 1. This will have wrong SigmaInv.
#   # For each debugonce call, I step through until after the invUlist is set up, then record the "blup_info_<i>" object, then continue to the end.
#   debugonce(blup.mixmeta)
#   test1 <- mixmeta::blup(package_resB, vcov=TRUE, se=TRUE, level=1, type="outcome" )
#   .GlobalEnv$blup_info_1 <- list(invUlist = invUlist, groups = groups, ord = ord, Zlist = Zlist, Z = Z, ZPZlist = ZPZlist, Slist = Slist, reslist = reslist, Xlist = Xlist, object=object)
#   # blups including both spline and residuals RE
#   debugonce(blup.mixmeta)
#   test2 <- mixmeta::blup(package_resB, vcov=TRUE, se=TRUE, level=2, type="outcome" )
#   .GlobalEnv$blup_info_2 <- list(invUlist = invUlist, groups = groups, ord = ord, Zlist = Zlist, Z = Z, ZPZlist = ZPZlist, Slist = Slist, reslist = reslist, Xlist = Xlist, object=object)
#   # blups including spline but not residuals RE by moving residuals RE into S
#   debugonce(blup.mixmeta)
#   test3 <- mixmeta::blup(orig_res$fit_for_blups, vcov=TRUE, se=TRUE, level=1, type="outcome" )
#   .GlobalEnv$blup_info_3 <- list(invUlist = invUlist, groups = groups, ord = ord, Zlist = Zlist, Z = Z, ZPZlist = ZPZlist, Slist = Slist, reslist = reslist, Xlist = Xlist, object=object)

#   o1 <- order(blup_info_1$ord)
#   o2 <- order(blup_info_2$ord)
#   o3 <- order(blup_info_3$ord)
#   identical(o1, o2) # TRUE
#   identical(o1, o3) # FALSE
#   ## The blup results from the three calls above should be identical but are not.
#   ## Here I will show what is happening
#   maxdiff <- \(a,b) (a-b)[which.max(abs(a - b))]
#   plot(blup_info_1$Z[o1,], blup_info_2$Z[[1]][o2, ]); abline(0,1)
#   maxdiff(blup_info_1$Z[o1, ], blup_info_2$Z[[1]][o2,]) ## matches
#   maxdiff(blup_info_1$Z[o1, ], blup_info_3$Z[o3,]) ## matches, no that I have ordering fixed.

#     ## Sigma as it is constructed
#   Sigma_1 <- (blup_info_1$ZPZlist[[1]] + blup_info_1$Slist[[1]])[o1, o1]
#   Sigma_2 <- (blup_info_2$ZPZlist[[1]] + blup_info_2$Slist[[1]])[o2, o2]
#   Sigma_3 <- (blup_info_3$ZPZlist[[1]] + blup_info_3$Slist[[1]])[o3, o3]
#   maxdiff(Sigma_1, Sigma_2) # These should match but do not
#   maxdiff(Sigma_2, Sigma_3) # These should match and do
#   maxdiff(Sigma_1, Sigma_3) # These should match but do not
#   ## The inverse covariance of Y, Sigma^Inv, is tcrossprod(invUlist[[i]])
#   SigmaInv_1 <- tcrossprod(blup_info_1$invUlist[[1]])[o1, o1]
#   SigmaInv_2 <- tcrossprod(blup_info_2$invUlist[[1]])[o2, o2]
#   SigmaInv_3 <- tcrossprod(blup_info_3$invUlist[[1]])[o3, o3]
#   maxdiff(SigmaInv_1, SigmaInv_2) ## These should match but do not
#   maxdiff(SigmaInv_2, SigmaInv_3) ## These should match and do
#   maxdiff(SigmaInv_1, SigmaInv_3) ## These should match but do not
#   ## I have demonstrated that the test1 version is wrong because
#   ## the second level of random effects is omitted from Sigma.
#   ## Now I will check that ZPZlist[[1]] matches for 1 and 3 but not 2
#   ## and Slist[[1]] matches for 1 and 2 but not 3.
#   ## The reason ZPZ should match for 1 and 3 is that this is Z Psi Z
#   ## and both of these methods are including only the spline random effects in that
#   ## The reason Slist[[1]] should match for 1 and 2 is that for 3 I have moved
#   ## the residual RE into S.
#     ZPZ_1 <- blup_info_1$ZPZlist[[1]][o1, o1]
#   ZPZ_2 <- blup_info_2$ZPZlist[[1]][o2, o2]
#   ZPZ_3 <- blup_info_3$ZPZlist[[1]][o3, o3]
#   maxdiff(ZPZ_1, ZPZ_2) ## These should not match and do not
#   maxdiff(ZPZ_1, ZPZ_3) ## These should match and do
#   maxdiff(ZPZ_2, ZPZ_3) ## These should not match and do not
#   S_1 <- blup_info_1$Slist[[1]][o1, o1]
#   S_2 <- blup_info_2$Slist[[1]][o2, o2]
#   S_3 <- blup_info_3$Slist[[1]][o3, o3]
#   maxdiff(S_1, S_2) ## These should match and do
#   maxdiff(S_1, S_3) ## These should not match and don't. max diff is the same as for ZPZ mismatches above
#   maxdiff(S_2, S_3) ## These should not match and don't
#   ###

#   expect_equal(pred_orig[,'blup'], pred_resA[,'blup'], tolerance= 0.01)
#   expect_equal(pred_orig[,'blup'], pred_resB[,'blup'], tolerance = 0.01)
#   expect_equal(pred_orig[,'se'], pred_resA[,'se'])
#   expect_equal(pred_orig[,'se'], pred_resB[,'se'])
#   })

