# Simulate basic meta-analysis data
sim_basic_meta <- function(n_groups = 10,
                           n_per_group = 5,
                           sd_group = 0.5,
                           intercept = 2,
                           slope = 0.8,
                           S = rep(0.2^2, n_groups * n_per_group)) {
  n <- n_groups * n_per_group
  group <- rep(LETTERS[1:n_groups], each = n_per_group)
  x <- runif(n, 0, 1)
  u_group <- rnorm(n_groups, 0, sd_group)
  y <- intercept + slope * x + rep(u_group, each = n_per_group) + rnorm(n, 0, sqrt(S))
  data.frame(x = x, group = group, y = y, S = S)
}

sim_data_for_two_splines <- function(x = as.numeric(1:20),
                                     y = as.numeric(1:20),
                                     n = length(x) * length(y),
                                     intercept = 2,
                                     coeff_x = 0.8,
                                     coeff_y = -0.4,
                                     S = rep(0.2^2, n))
 {
  grid <- expand.grid(x = x, y = y)
  xg <- grid$x
  yg <- grid$y
  f_x <- sin(2 * pi * xg/12)
  f_y <- cos(2 * pi * yg/15)
  mu <- intercept + coeff_x * xg + coeff_y * yg + f_x + f_y
  z <- rnorm(n, mean = mu, sd = sqrt(S))
  data.frame(x = xg, y = yg, z = z, S = S)
 }

