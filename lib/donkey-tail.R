set.seed(245005) # random.org

# The names of the managed features: here, the ionome
managed_features <- df_npk %>%
  select(starts_with("Leaf")) %>%
  names()

# statistics needed to compute the mahalanobis distance
contr_mean <- apply(train_baked %>% select(one_of(managed_features)), 2, mean) # should be zero
contr_sd <- apply(train_baked %>% select(one_of(managed_features)), 2, mean)
contr_icov <- solve(cov(train_baked %>% select(one_of(managed_features))))
crit_dist <- qchisq(p = 0.999, df = length(contr_mean))

# statistics to back transform to original scale
managed_mean <- mean_npk_train[names(mean_npk_train) %in% managed_features]
managed_sd <- sd_npk_train[names(sd_npk_train) %in% managed_features]

# markov chain parameters
max_iter <- 10
radius <- c(0.5)
radius_factor <- 1.25 # for adaptative search
radius_limits <- c(0.5, 1.5)
n_rad <- 5000

# initialisation of variables
n_obs <- nrow(df_npk)
opt_bal <- matrix(ncol = length(managed_features), nrow = n_obs)
opt_yield <- matrix(ncol = 2, nrow = n_obs)
ait_dist <- rep(NA, n_obs)


for (n_ in 1:n_obs) {
  print(paste("row", n_, "of", n_obs))

  observation <- df_npk[n_, ] #%>% select(-yield)

  # the features are initialized with max_iter rows and NAs where managed features 
  # will be optimized
  features_iterations <- observation %>%
    slice(rep(row_number(), max_iter)) %>%
    bake(npk_recipe, .) %>%
    select(-yield)
  for (col in managed_features) features_iterations[2:max_iter, col] <- NA
  
  iterations_yield <- c(step_back(predict(m_fit, newdata = bake(npk_recipe, observation))))
  
  iterations_managed <- features_iterations %>%
    select(managed_features) %>%
    as.matrix()
  
  # the information complementary to the managed features is
  # extracted in a table, then replicated `n_rad` times to be binded further
  # on to the managed featrures scanned around the best combination
  observation_conditional <- observation %>% 
    bake(npk_recipe, .) %>%
    select(-one_of(managed_features), -yield) %>%
    slice(rep(row_number(), n_rad))

  # fire the Markov chain
  for (i in 2:max_iter) {
    # print(paste(i, "/", max_iter))
    offset <- matrix(runif(length(managed_features) * n_rad, -1, 1), 
                     ncol = length(managed_features),
                     nrow = n_rad)
    offset <- t(apply(offset, 1, function(x) radius[i-1] * x / sqrt(sum(x^2))))
    offset <- offset * runif(length(offset), 0, 1)
    observation_search <- t(apply(offset, 1, function(x) x + iterations_managed[i-1, ]))

    # Compute the Mahalanobis distance
    maha_dist <- mahalanobis(observation_search, contr_mean, contr_icov, inverted = TRUE)
    # filter out search ionomes outside the Mahalanobis distance limit
    df_search <- data.frame(observation_search) %>%
      bind_cols(observation_conditional) %>%
      filter(maha_dist < crit_dist) %>%
      mutate(doseN = abs(doseN), # assure positive dosage
             doseP = abs(doseP),
             doseK = abs(doseK))

    if(nrow(df_search) == 0) { # if no points are generated in the hyper ellipsoid, keep the reference but increase radius
      print(paste("Iteration", i, "- all points are out of the hyperellipsoid."))
      iterations_yield[i] <- iterations_yield[i-1]
      iterations_managed[i, ] <- iterations_managed[i-1, ]
      # increase the radius
      radius[i] <- radius[i - 1] * radius_factor
      radius[i] <- ifelse(radius[i] > radius_limits[2], radius_limits[2], radius[i])
    } else {
      # Compute predicted yield
      yield_stochastic <- step_back(predict(m_fit, newdata = df_search))
      if(max(yield_stochastic) > iterations_yield[i-1]) {
        print(paste("Iteration", i, "- yield improved to", max(yield_stochastic)))
        iterations_yield[i] <- max(yield_stochastic)
        iterations_managed[i, ] <- df_search[which.max(yield_stochastic), managed_features] %>% unlist()
        # decrease the radius
        radius[i] <- radius[i - 1] / radius_factor
        radius[i] <- ifelse(radius[i] < radius_limits[1], radius_limits[1], radius[i])
      } else {
        print(paste("Iteration", i, "- no yield improvement."))
        iterations_yield[i] <- iterations_yield[i-1]
        iterations_managed[i, ] <- iterations_managed[i-1, ]
        # increase the radius
        radius[i] <- radius[i - 1] * radius_factor
        radius[i] <- ifelse(radius[i] > radius_limits[2], radius_limits[2], radius[i])
      }
    }
  }
  
  # backtransform managed features to their original scale
  iterations_managed_unsc <- apply(iterations_managed, 1, function(x) x * managed_sd + managed_mean) %>% t()
  
  # we extract the last component of the chain and transform it
  # to a composition
  opt_bal[n_, ] <- iterations_managed_unsc[max_iter, ]
  opt_yield[n_, 1] <- iterations_yield[1]
  opt_yield[n_, 2] <- iterations_yield[max_iter]
  ait_dist[n_] <- sqrt(sum((iterations_managed_unsc[1, ] - iterations_managed_unsc[max_iter, ]) ^ 2))
}
