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
