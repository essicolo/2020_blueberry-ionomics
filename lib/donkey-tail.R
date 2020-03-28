set.seed(245005) # random.org
max_iter <- 20 
radius <- c(0.5)
radius_factor <- 1.25 # for adaptative search
radius_limits <- c(0.5, 1.5)
n_rad <- 5000

n_obs <- nrow(df_npk)
opt_bal <- matrix(ncol = nrow(sbp_leaf), nrow = n_obs)
opt_yield <- matrix(ncol = 2, nrow = n_obs)
ait_dist <- rep(NA, n_obs)

for (n in 1:n_obs) {
  print(paste("row", n, "of", n_obs))

  misbal_observation <- df_npk[n, ]

  # the ionome is initialized
  ref_leaf <- misbal_observation %>%
    select(starts_with("Leaf")) %>%
    unlist()

  # predicted yield is initialized
  yield_init <- predict(m_fit, newdata = bake(npk_recipe, misbal_observation)) *
    sd(npk_train$yield) +
    mean(npk_train$yield)

  # the information complementary to to ionome is
  # extracted in a table, then replicated `n_rad` times to be binded further
  # on to the perturbed ionomes scanned around the best ionome
  misbal_observation_noleaf <- misbal_observation %>%
    select(-starts_with("Leaf"))
  misbal_observation_noleaf_stacked <- do.call(
    "rbind",
    replicate(n_rad, misbal_observation_noleaf, simplify = FALSE)
  )

  # the initial values are put into vectors and matrices
  ref_yield <- c(yield_init)
  ref_leaf <- matrix(ncol = nrow(sbp_leaf), nrow = max_iter)

  ref_leaf[1, ] <- misbal_observation %>%
    select(starts_with("Leaf")) %>%
    unlist()

  # fire the Markov chain
  for (i in 2:max_iter) {
    offset <- matrix(runif(ncol(ref_leaf) * n_rad, -1, 1),
                     ncol = ncol(ref_leaf),
                     nrow = n_rad)
    offset <- t(apply(offset, 1, function(x) radius[i-1] * x / sqrt(sum(x^2))))
    leaf_search <- t(apply(offset, 1, function(x) x + ref_leaf[i-1, ] ))

    maha_dist <- mahalanobis(leaf_search, bal_mean, bal_icov, inverted = TRUE)
    leaf_search <- data.frame(leaf_search) %>% filter(maha_dist < crit_dist)
    names(leaf_search) <- paste0("Leaf_", leaf_bal_def)

    df_search <- leaf_search %>%
      bind_cols(misbal_observation_noleaf_stacked %>% filter(maha_dist < crit_dist))

    if(nrow(df_search) == 0) { # if no points are generated in the hyper ellipsoid, keep the reference but increase radius
      ref_yield[i] <- ref_yield[i-1]
      ref_leaf[i, ] <- ref_leaf[i-1, ]
      # increase the radius
      radius[i] <- radius[i - 1] * radius_factor
      radius[i] <- ifelse(radius[i] > radius_limits[2], radius_limits[2], radius[i])
    } else {
      yield_stochastic <- predict(m_fit, newdata = bake(npk_recipe, df_search)) *
        sd(npk_train$yield) + mean(npk_train$yield)
      if(max(yield_stochastic) > ref_yield[i-1]) {
        ref_yield[i] <- max(yield_stochastic)
        ref_leaf[i, ] <- leaf_search[which.max(yield_stochastic), ] %>%
          unlist()
        # decrease the radius
        radius[i] <- radius[i - 1] / radius_factor
        radius[i] <- ifelse(radius[i] < radius_limits[1], radius_limits[1], radius[i])
      } else {
        ref_yield[i] <- ref_yield[i-1]
        ref_leaf[i, ] <- ref_leaf[i-1, ]
        # increase the radius
        radius[i] <- radius[i - 1] * radius_factor
        radius[i] <- ifelse(radius[i] > radius_limits[2], radius_limits[2], radius[i])
      }
    }
  }


  # we extract the last component of the chain and transform it
  # to a composition
  opt_bal[n, ] <- ref_leaf[max_iter, ]
  opt_yield[n, 1] <- yield_init
  opt_yield[n, 2] <- ref_yield[max_iter]
  ait_dist[n] <- sqrt(sum((ref_leaf[1, ] - ref_leaf[max_iter, ]) ^ 2))
}
