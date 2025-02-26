# library(testthat)

# used to break weighted counts into weighted samples
break_sum <- function(summ, n) {
  if (n == 1) {
    dif <- c(summ)
    return(dif)
  }
  if (n == 0) {
    dif <- NA
    return(dif)
  }
  if (summ == 0) {
    dif <- NA
    return(dif)
  }
  random_numbers <- sort(c(0, runif(n, 0, summ)))
  dif <- diff(random_numbers)
  dif[length(dif)] <- summ - sum(dif[-length(dif)])
  return(dif)
}

# Scenario A: weighted, different bin_widths, different discrete weights, mixed case
test_that("Scenario A", {
  # create random observed (histogram) count data set
  set.seed(101)
  bin_number <- sample(10:20, 1)
  random_vector <- sample(c(1:100), bin_number, replace = TRUE)
  random_vector <- random_vector / 100 #=bin_width
  tics <- c(0, cumsum(random_vector))
  continous_domain <- c(0, max(tics))
  interval_widths <- diff(tics)
  mids_dist <- c(interval_widths / 2)
  mids <- tics[1:(bin_number)] + mids_dist
  n_discrete <- sample(1:4, 1, prob = c(0.25, 0.25, 0.25, 0.25))
  discrete_values <- runif(n_discrete, min = 0, max = max(tics))
  discrete_values <- sort(discrete_values)
  mids_and_discrete <- sort(c(mids, discrete_values))
  n <- bin_number + n_discrete
  indizes <- which(mids_and_discrete %in% discrete_values)
  weights_discrete <- runif(n_discrete, min = 0.1, max = 5)
  width_and_weights <- interval_widths
  for (i in seq_along(weights_discrete)) {
    width_and_weights <-
      append(width_and_weights, weights_discrete[i], after = indizes[i] - 1)
  }
  featureA <- rep(c('a', 'b'), 3)
  featureB <- rep(rep(c('f', 'g', 'h'), each = 2), 1)
  group_ID <- rep(1:6, each = 1)
  data <- data.frame(
    counts = rep(NA, 6),
    weighted_counts = rep(NA, 6),
    response =  rep(NA, 6),
    featureA = featureA,
    featureB = featureB,
    group_ID = group_ID,
    gam_weights = rep(NA, 6),
    gam_offsets = rep(NA, 6),
    Delta = rep(NA, 6)
  )
  data <- data[rep(row.names(data), each = n), ]
  data$response <- rep(mids_and_discrete, 6)
  data$Delta <- rep(width_and_weights, 6)
  data$counts <- sample(0:100, n * 6, replace = TRUE)
  data$weighted_counts <-
    data$counts + runif(n * 6, min = 0, max = 50)
  data$weighted_counts[data$counts == 0] <- 0
  data$gam_weights <-
    ifelse(data$counts != 0, data$weighted_counts / data$counts, 1)
  data$gam_offsets <- -log(data$gam_weights)
  data<-data%>%mutate(discrete=response%in%discrete_values)

  ############# break it into original Data set #########
  original <- data[rep(row.names(data), times = data$counts), ]

  sample_weighted <- c()
  for (i in 1:length(data$weighted_counts)) {
    if (data$weighted_counts[i] != 0) {
      vect <- break_sum(data$weighted_counts[i], data$counts[i])
      sample_weighted <- append(sample_weighted, vect)
    }
  }
  original$sample_weight <- na.omit(sample_weighted)
  colnames(original)[3] <- "mids"
  original$response <- 0
  for (i in 1:length(original$response)) {
    if (!original$mids[i] %in% discrete_values) {
      original$response[i] <-
        original$mids[i] + runif(1,
                                 min = -original$Delta[i] / 2,
                                 max = original$Delta[i] / 2)
    }
    else{
      original$response[i] <- original$mids[i]
    }
  }
  original <-
    original[, c("response", "featureA", "featureB", "sample_weight")]
  original = original[sample(1:nrow(original)), ]
  data_backtransformed <-
    data2counts(
      dta = original,
      var_vec = c(2:3),
      y = 1,
      sample_weights = 4,
      bin_width = interval_widths ,
      values_discrete = discrete_values,
      weights_discrete = weights_discrete ,
      domain_continuous = continous_domain
    )
  data$combi <- paste0(data$featureA, data$featureB)
  data_backtransformed$combi <-
    paste0(data_backtransformed$featureA,
           data_backtransformed$featureB)
  data_check <- data %>% arrange(combi)
  data_bt_check <- data_backtransformed %>% arrange(combi)

  expect_equal(data_check$counts, data_bt_check$counts)
  expect_equal(round(data_check$weighted_counts, 5),
               round(data_bt_check$weighted_counts, 5))
  expect_equal(round(data_check$response, 5),
               round(data_bt_check$response, 5))
  expect_equal(data_check$combi, data_bt_check$combi)
  expect_equal(round(data_check$gam_weights, 5),
               round(data_bt_check$gam_weights, 5))
  expect_equal(round(data_check$gam_offsets, 5),
               round(data_bt_check$gam_offsets, 5))
  expect_equal(round(data_check$Delta, 5), round(data_bt_check$Delta, 5))
  expect_equal(data_check$discrete, data_bt_check$discrete)

})

# Scenario B: weighted, equidistant bins, different discrete weights, mixed case
test_that("Scenario B", {
  set.seed(101)
  # create random observed (histogram) count data set
  bin_number <- sample(10:20, 1)# n-1=number of bins
  max_continous <- sample(10:100, 1) / 10
  continous_domain <- c(0, max_continous)
  interval_width <- rep(max_continous / bin_number, bin_number)
  mids_dist <- c(interval_width / 2)
  tics <- c(0, cumsum(interval_width))
  mids <- tics[1:(bin_number)] + mids_dist
  n_discrete <- sample(1:4, 1, prob = c(0.25, 0.25, 0.25, 0.25))
  discrete_values <- runif(n_discrete, min = 0, max = max(tics))
  discrete_values <- sort(discrete_values)
  mids_and_discrete <- sort(c(mids, discrete_values))
  n <- bin_number + n_discrete
  indizes <- which(mids_and_discrete %in% discrete_values)
  weights_discrete <- runif(n_discrete, min = 0.1, max = 5)
  width_and_weights <- interval_width
  for (i in seq_along(weights_discrete)) {
    width_and_weights <-
      append(width_and_weights, weights_discrete[i], after = indizes[i] - 1)
  }
  featureA <- rep(c('a', 'b'), 3)
  featureB <- rep(rep(c('f', 'g', 'h'), each = 2), 1)
  group_ID <- rep(1:6, each = 1)
  data <- data.frame(
    counts = rep(NA, 6),
    weighted_counts = rep(NA, 6),
    response =  rep(NA, 6),
    featureA = featureA,
    featureB = featureB,
    group_ID = group_ID,
    gam_weights = rep(NA, 6),
    gam_offsets = rep(NA, 6),
    Delta = rep(NA, 6)
  )
  data <- data[rep(row.names(data), each = n), ]
  data$response <- rep(mids_and_discrete, 6)
  data$Delta <- rep(width_and_weights, 6)
  data$counts <- sample(0:100, n * 6, replace = TRUE)
  data$weighted_counts <-
    data$counts + runif(n * 6, min = 0, max = 50)
  data$weighted_counts[data$counts == 0] <- 0
  data$gam_weights <-
    ifelse(data$counts != 0, data$weighted_counts / data$counts, 1)
  data$gam_offsets <- -log(data$gam_weights)
  data<-data%>%mutate(discrete=response%in%discrete_values)
  ############# break it into original Data set #########
  original <- data[rep(row.names(data), times = data$counts), ]
  sample_weighted <- c()
  for (i in 1:length(data$weighted_counts)) {
    if (data$weighted_counts[i] != 0) {
      vect <- break_sum(data$weighted_counts[i], data$counts[i])
      sample_weighted <- append(sample_weighted, vect)
    }
  }
  original$sample_weight <- na.omit(sample_weighted)
  colnames(original)[3] <- "mids"
  original$response <- 0
  for (i in 1:length(original$response)) {
    if (!original$mids[i] %in% discrete_values) {
      original$response[i] <-
        original$mids[i] + runif(1,
                                 min = -original$Delta[i] / 2,
                                 max = original$Delta[i] / 2)
    }
    else{
      original$response[i] <- original$mids[i]
    }
  }
  original <-
    original[, c("response", "featureA", "featureB", "sample_weight")]
  original = original[sample(1:nrow(original)), ]
  data_backtransformed <-
    data2counts(
      original,
      var_vec = c(2:3),
      y = 1,
      sample_weights = 4,
      bin_number = bin_number,
      values_discrete = discrete_values,
      weights_discrete = weights_discrete ,
      domain_continuous = continous_domain
    )
  data$combi <- paste0(data$featureA, data$featureB)
  data_backtransformed$combi <-
    paste0(data_backtransformed$featureA,
           data_backtransformed$featureB)
  data_check <- data %>% arrange(combi)
  data_bt_check <- data_backtransformed %>% arrange(combi)
  expect_equal(data_check$counts, data_bt_check$counts)
  expect_equal(round(data_check$weighted_counts, 5),
               round(data_bt_check$weighted_counts, 5))
  expect_equal(round(data_check$response, 5),
               round(data_bt_check$response, 5))
  expect_equal(data_check$combi, data_bt_check$combi)
  expect_equal(round(data_check$gam_weights, 5),
               round(data_bt_check$gam_weights, 5))
  expect_equal(round(data_check$gam_offsets, 5),
               round(data_bt_check$gam_offsets, 5))
  expect_equal(round(data_check$Delta, 5), round(data_bt_check$Delta, 5))
  expect_equal(data_check$discrete, data_bt_check$discrete)

})

# Scenario C: weighted, equidistant bins, only continuous
test_that("Scenario C", {
  set.seed(101)
  # create random observed (histogram) count data set
  bin_number <- sample(10:20, 1)
  max_continous <- sample(10:100, 1) / 10
  continous_domain <- c(0, max_continous)
  interval_width <- rep(max_continous / bin_number, bin_number)
  mids_dist <- c(interval_width / 2)
  tics <- c(0, cumsum(interval_width))
  mids <- tics[1:(bin_number)] + mids_dist
  n <- length(mids)
  width_and_weights <- interval_width
  featureA <- rep(c('a', 'b'), 3)
  featureB <- rep(rep(c('f', 'g', 'h'), each = 2), 1)
  group_ID <- rep(1:6, each = 1)
  data <- data.frame(
    counts = rep(NA, 6),
    weighted_counts = rep(NA, 6),
    response =  rep(NA, 6),
    featureA = featureA,
    featureB = featureB,
    group_ID = group_ID,
    gam_weights = rep(NA, 6),
    gam_offsets = rep(NA, 6),
    Delta = rep(NA, 6)
  )
  data <- data[rep(row.names(data), each = n), ]
  data$response <- rep(mids, 6)
  data$Delta <- rep(width_and_weights, 6)
  data$counts <- sample(0:100, n * 6, replace = TRUE)
  data$weighted_counts <-
    data$counts + runif(n * 6, min = 0, max = 50)
  data$weighted_counts[data$counts == 0] <- 0
  data$gam_weights <-
    ifelse(data$counts != 0, data$weighted_counts / data$counts, 1)
  data$gam_offsets <- -log(data$gam_weights)
  data$discrete<-FALSE
  ############# break it into original Data set #########
  original <- data[rep(row.names(data), times = data$counts), ]
  sample_weighted <- c()
  for (i in 1:length(data$weighted_counts)) {
    if (data$weighted_counts[i] != 0) {
      vect <- break_sum(data$weighted_counts[i], data$counts[i])
      sample_weighted <- append(sample_weighted, vect)
    }
  }
  original$sample_weight <- na.omit(sample_weighted)
  colnames(original)[3] <- "mids"
  original$response <- 0
  for (i in 1:length(original$response)) {
    original$response[i] <-
      original$mids[i] + runif(1,
                               min = -original$Delta[i] / 2,
                               max = original$Delta[i] / 2)
  }
  original <-
    original[, c("response", "featureA", "featureB", "sample_weight")]
  original = original[sample(1:nrow(original)), ]
  data_backtransformed_2 <-
    data2counts(
      original,
      var_vec = c(2:3),
      y = 1,
      sample_weights = 4,
      bin_number = bin_number,
      values_discrete = FALSE,
      weights_discrete = FALSE ,
      domain_continuous = continous_domain
    )
  data_backtransformed_3 <-
    data2counts(
      original,
      var_vec = c(2:3),
      y = 1,
      sample_weights = 4,
      bin_width = interval_width,
      values_discrete = FALSE,
      weights_discrete = FALSE ,
      domain_continuous = continous_domain
    )
  data_backtransformed_4 <-
    data2counts(
      original,
      var_vec = c(2:3),
      y = 1,
      sample_weights = 4,
      bin_number = bin_number,
      bin_width = interval_width,
      values_discrete = FALSE,
      weights_discrete = FALSE ,
      domain_continuous = continous_domain
    )

  expect_equal(data_backtransformed_2, data_backtransformed_3)
  expect_equal(data_backtransformed_2, data_backtransformed_4)
  data$combi <- paste0(data$featureA, data$featureB)
  data_backtransformed_2$combi <-
    paste0(data_backtransformed_2$featureA,
           data_backtransformed_2$featureB)
  data_check <- data %>% arrange(combi)
  data_bt_check <- data_backtransformed_2 %>% arrange(combi)
  expect_equal(data_check$counts, data_bt_check$counts)
  expect_equal(round(data_check$weighted_counts, 5),
               round(data_bt_check$weighted_counts, 5))
  expect_equal(round(data_check$response, 5),
               round(data_bt_check$response, 5))
  expect_equal(data_check$combi, data_bt_check$combi)
  expect_equal(round(data_check$gam_weights, 5),
               round(data_bt_check$gam_weights, 5))
  expect_equal(round(data_check$gam_offsets, 5),
               round(data_bt_check$gam_offsets, 5))
  expect_equal(round(data_check$Delta, 5), round(data_bt_check$Delta, 5))
  expect_equal(data_check$discrete, data_bt_check$discrete)
})

# Scenario D: weighted, equidistant bins, different discrete weights, only discrete
test_that("Scenario D", {
  set.seed(101)
  # create random observed (histogram) count data set
  break_sum <- function(summ, n) {
    if (n == 1) {
      dif <- c(summ)
      return(dif)
    }
    if (n == 0) {
      dif <- NA
      return(dif)
    }
    if (summ == 0) {
      dif <- NA
      return(dif)
    }
    random_numbers <- sort(c(0, runif(n, 0, summ)))
    dif <- diff(random_numbers)
    dif[length(dif)] <- summ - sum(dif[-length(dif)])
    return(dif)
  }
  n_discrete <- sample(1:10, 1)
  discrete_values <- runif(n_discrete)
  discrete_values <- sort(discrete_values)
  n <- n_discrete
  weights_discrete <- runif(n_discrete, min = 0.1, max = 5)
  featureA <- rep(c('a', 'b'), 3)
  featureB <- rep(rep(c('f', 'g', 'h'), each = 2), 1)
  group_ID <- rep(1:6, each = 1)
  data <- data.frame(
    counts = rep(NA, 6),
    weighted_counts = rep(NA, 6),
    response =  rep(NA, 6),
    featureA = featureA,
    featureB = featureB,
    group_ID = group_ID,
    gam_weights = rep(NA, 6),
    gam_offsets = rep(NA, 6),
    Delta = rep(NA, 6)
  )
  data <- data[rep(row.names(data), each = n), ]
  data$response <- rep(discrete_values, 6)
  data$Delta <- rep(weights_discrete, 6)
  data$counts <- sample(0:100, n * 6, replace = TRUE)
  data$weighted_counts <-
    data$counts + runif(n * 6, min = 0, max = 50)
  data$weighted_counts[data$counts == 0] <- 0
  data$gam_weights <-
    ifelse(data$counts != 0, data$weighted_counts / data$counts, 1)
  data$gam_offsets <- -log(data$gam_weights)
  data<-data%>%mutate(discrete=response%in%discrete_values)
  ############# break it into original Data set #########
  original <- data[rep(row.names(data), times = data$counts), ]
  sample_weighted <- c()
  for (i in 1:length(data$weighted_counts)) {
    if (data$weighted_counts[i] != 0) {
      vect <- break_sum(data$weighted_counts[i], data$counts[i])
      sample_weighted <- append(sample_weighted, vect)
    }
  }
  original$sample_weight <- na.omit(sample_weighted)
  original <-
    original[, c("response", "featureA", "featureB", "sample_weight")]
  original = original[sample(1:nrow(original)), ]
  data_backtransformed_5 <-
    data2counts(
      original,
      var_vec = c(2:3),
      y = 1,
      sample_weights = 4,
      values_discrete = discrete_values,
      weights_discrete = weights_discrete ,
      domain_continuous = FALSE
    )
  data$combi <- paste0(data$featureA, data$featureB)
  data_backtransformed_5$combi <-
    paste0(data_backtransformed_5$featureA,
           data_backtransformed_5$featureB)
  data_check <- data %>% arrange(combi)
  data_bt_check <- data_backtransformed_5 %>% arrange(combi)
  expect_equal(data_check$counts, data_bt_check$counts)
  expect_equal(round(data_check$weighted_counts, 5),
               round(data_bt_check$weighted_counts, 5))
  expect_equal(round(data_check$response, 5),
               round(data_bt_check$response, 5))
  expect_equal(data_check$combi, data_bt_check$combi)
  expect_equal(round(data_check$gam_weights, 5),
               round(data_bt_check$gam_weights, 5))
  expect_equal(round(data_check$gam_offsets, 5),
               round(data_bt_check$gam_offsets, 5))
  expect_equal(round(data_check$Delta, 5), round(data_bt_check$Delta, 5))
  expect_equal(data_check$discrete, data_bt_check$discrete)
})

# Scenario E: unweighted, different bin_widths, different discrete weights, mixed case
test_that("Scenario E", {
  set.seed(101)
  # create random observed (histogram) count data set
  bin_number <- sample(10:20, 1)# n-1=number of bins
  random_vector <- sample(c(1:100), bin_number, replace = TRUE)
  random_vector <- random_vector / 100 #=bin_width
  tics <- c(0, cumsum(random_vector))
  continous_domain <- c(0, max(tics))
  interval_widths <- diff(tics)
  mids_dist <- c(interval_widths / 2)
  mids <- tics[1:(bin_number)] + mids_dist
  n_discrete <- sample(1:4, 1, prob = c(0.25, 0.25, 0.25, 0.25))
  discrete_values <- runif(n_discrete, min = 0, max = max(tics))
  discrete_values <- sort(discrete_values)
  mids_and_discrete <- sort(c(mids, discrete_values))
  n <- bin_number + n_discrete
  indizes <- which(mids_and_discrete %in% discrete_values)
  weights_discrete <- runif(n_discrete, min = 0.1, max = 5)
  width_and_weights <- interval_widths
  for (i in seq_along(weights_discrete)) {
    width_and_weights <-
      append(width_and_weights, weights_discrete[i], after = indizes[i] - 1)
  }
  featureA <- rep(c('a', 'b'), 3)
  featureB <- rep(rep(c('f', 'g', 'h'), each = 2), 1)
  group_ID <- rep(1:6, each = 1)
  data <- data.frame(
    counts = rep(NA, 6),
    response =  rep(NA, 6),
    featureA = featureA,
    featureB = featureB,
    group_ID = group_ID,
    gam_weights = rep(NA, 6),
    gam_offsets = rep(NA, 6),
    Delta = rep(NA, 6)
  )
  data <- data[rep(row.names(data), each = n), ]
  data$response <- rep(mids_and_discrete, 6)
  data$Delta <- rep(width_and_weights, 6)
  data$counts <- sample(0:100, n * 6, replace = TRUE)
  data$gam_weights <- 1
  data$gam_offsets <- 0
  data<-data%>%mutate(discrete=response%in%discrete_values)
  ############# break it into original Data set #########
  original <- data[rep(row.names(data), times = data$counts), ]
  break_sum <- function(summ, n) {
    if (n == 1) {
      dif <- c(summ)
      return(dif)
    }
    if (n == 0) {
      dif <- NA
      return(dif)
    }
    if (summ == 0) {
      dif <- NA
      return(dif)
    }
    random_numbers <- sort(c(0, runif(n, 0, summ)))
    dif <- diff(random_numbers)
    dif[length(dif)] <- summ - sum(dif[-length(dif)])
    return(dif)
  }
  colnames(original)[2] <- "mids"
  original$response <- 0
  for (i in 1:length(original$response)) {
    if (!original$mids[i] %in% discrete_values) {
      original$response[i] <-
        original$mids[i] + runif(1,
                                 min = -original$Delta[i] / 2,
                                 max = original$Delta[i] / 2)
    }
    else{
      original$response[i] <- original$mids[i]
    }
  }
  original <- original[, c("response", "featureA", "featureB")]
  original = original[sample(1:nrow(original)), ]
  data_backtransformed_6 <-
    data2counts(
      original,
      var_vec = c(2:3),
      y = 1,
      bin_width = interval_widths,
      values_discrete = discrete_values,
      weights_discrete = weights_discrete ,
      domain_continuous = continous_domain
    )
  data$combi <- paste0(data$featureA, data$featureB)
  data_backtransformed_6$combi <-
    paste0(data_backtransformed_6$featureA,
           data_backtransformed_6$featureB)
  data_check <- data %>% arrange(combi)
  data_bt_check <- data_backtransformed_6 %>% arrange(combi)
  expect_equal(data_check$counts, data_bt_check$counts)
  expect_equal(round(data_check$response, 5),
               round(data_bt_check$response, 5))
  expect_equal(data_check$combi, data_bt_check$combi)
  expect_equal(round(data_check$gam_weights, 5),
               round(data_bt_check$gam_weights, 5))
  expect_equal(round(data_check$gam_offsets, 5),
               round(data_bt_check$gam_offsets, 5))
  expect_equal(round(data_check$Delta, 5), round(data_bt_check$Delta, 5))
  expect_equal(data_check$discrete, data_bt_check$discrete)
})

# Scenario F: weighted, different bin_widths, only one discrete value + discrete weights, mixed case
test_that("Scenario F", {
  set.seed(101)
  # create random observed (histogram) count data set
  bin_number <- sample(10:20, 1)# n-1=number of bins
  random_vector <- sample(c(1:100), bin_number, replace = TRUE)
  random_vector <- random_vector / 100 #=bin_width
  tics <- c(0, cumsum(random_vector))
  continous_domain <- c(0, max(tics))
  interval_widths <- diff(tics)
  mids_dist <- c(interval_widths / 2)
  mids <- tics[1:(bin_number)] + mids_dist
  n_discrete <- 1
  discrete_values <- runif(n_discrete, min = 0, max = max(tics))
  discrete_values <- sort(discrete_values)
  mids_and_discrete <- sort(c(mids, discrete_values))
  n <- bin_number + n_discrete
  indizes <- which(mids_and_discrete %in% discrete_values)
  weights_discrete <- runif(n_discrete, min = 0.1, max = 5)
  width_and_weights <- interval_widths
  for (i in seq_along(weights_discrete)) {
    width_and_weights <-
      append(width_and_weights, weights_discrete[i], after = indizes[i] - 1)
  }
  featureA <- rep(c('a', 'b'), 3)
  featureB <- rep(rep(c('f', 'g', 'h'), each = 2), 1)
  group_ID <- rep(1:6, each = 1)
  data <- data.frame(
    counts = rep(NA, 6),
    weighted_counts = rep(NA, 6),
    response =  rep(NA, 6),
    featureA = featureA,
    featureB = featureB,
    group_ID = group_ID,
    gam_weights = rep(NA, 6),
    gam_offsets = rep(NA, 6),
    Delta = rep(NA, 6)
  )
  data <- data[rep(row.names(data), each = n), ]
  data$response <- rep(mids_and_discrete, 6)
  data$Delta <- rep(width_and_weights, 6)
  data$counts <- sample(0:100, n * 6, replace = TRUE)
  data$weighted_counts <-
    data$counts + runif(n * 6, min = 0, max = 50)
  data$weighted_counts[data$counts == 0] <- 0
  data$gam_weights <-
    ifelse(data$counts != 0, data$weighted_counts / data$counts, 1)
  data$gam_offsets <- -log(data$gam_weights)
  data<-data%>%mutate(discrete=response%in%discrete_values)
  ############# break it into original Data set #########
  original <- data[rep(row.names(data), times = data$counts), ]
  sample_weighted <- c()
  for (i in 1:length(data$weighted_counts)) {
    if (data$weighted_counts[i] != 0) {
      vect <- break_sum(data$weighted_counts[i], data$counts[i])
      sample_weighted <- append(sample_weighted, vect)
    }
  }
  original$sample_weight <- na.omit(sample_weighted)
  colnames(original)[3] <- "mids"
  original$response <- 0
  for (i in 1:length(original$response)) {
    if (!original$mids[i] %in% discrete_values) {
      original$response[i] <-
        original$mids[i] + runif(1,
                                 min = -original$Delta[i] / 2,
                                 max = original$Delta[i] / 2)
    }
    else{
      original$response[i] <- original$mids[i]
    }
  }
  original <-
    original[, c("response", "featureA", "featureB", "sample_weight")]
  original = original[sample(1:nrow(original)), ]
  data_backtransformed <-
    data2counts(
      original,
      var_vec = c(2:3),
      y = 1,
      sample_weights = 4,
      bin_width = interval_widths ,
      values_discrete = discrete_values,
      weights_discrete = weights_discrete ,
      domain_continuous = continous_domain
    )
  data$combi <- paste0(data$featureA, data$featureB)
  data_backtransformed$combi <-
    paste0(data_backtransformed$featureA,
           data_backtransformed$featureB)
  data_check <- data %>% arrange(combi)
  data_bt_check <- data_backtransformed %>% arrange(combi)
  expect_equal(data_check$counts, data_bt_check$counts)
  expect_equal(round(data_check$weighted_counts, 5),
               round(data_bt_check$weighted_counts, 5))
  expect_equal(round(data_check$response, 5),
               round(data_bt_check$response, 5))
  expect_equal(data_check$combi, data_bt_check$combi)
  expect_equal(round(data_check$gam_weights, 5),
               round(data_bt_check$gam_weights, 5))
  expect_equal(round(data_check$gam_offsets, 5),
               round(data_bt_check$gam_offsets, 5))
  expect_equal(round(data_check$Delta, 5), round(data_bt_check$Delta, 5))
  expect_equal(data_check$discrete, data_bt_check$discrete)
})

# Scenario G: mids have to be shifted
test_that("Szenatio G", {
  # create small observed (histogram) count data set where one mid (0.25) has to
  # be shifted due to a discrete value of 0.25
  domain_continuous <- c(0, 1)
  values_discrete <- c(0, 0.25, 1)
  bin_number <- 2
  n = 5
  data <-
    data.frame(
      counts = NA,
      response = NA,
      featureA = c("a", "a", "b", "b"),
      featureB = c("c", "d", "c", "d"),
      group_ID = c(1:4),
      gam_weights = 1,
      gam_offsets = 0,
      Delta = NA
    )
  data <- data[rep(row.names(data), each = n), ]
  # include shift in data$response for mid=0.25
  data$response <- rep(c(0, 0.25, 0.25 + 0.0001 * 0.25, 0.75, 1), 4)
  data$Delta <- rep(c(1, 1, 0.5, 0.5, 1), 4)
  data$counts <- c(1:20)
  data<-data%>%mutate(discrete=response%in%values_discrete)
  original <- data[rep(row.names(data), times = data$counts), ]
  colnames(original)[2] <- "mids"
  original$response <- 0
  for (i in 1:length(original$response)) {
    if (!original$mids[i] %in% values_discrete) {
      original$response[i] <-
        original$mids[i] + runif(1,
                                 min = -original$Delta[i] / 2,
                                 max = original$Delta[i] / 2)
    }
    else{
      original$response[i] <- original$mids[i]
    }
  }
  original <- original[, c("response", "featureA", "featureB")]
  original = original[sample(1:nrow(original)), ]
  data_backtransformed <-
    data2counts(
      original,
      var_vec = c(2:3),
      y = 1,
      bin_width = 0.5 ,
      values_discrete = c(0, 0.25, 1)
    )
  data$combi <- paste0(data$featureA, data$featureB)
  data_backtransformed$combi <-
    paste0(data_backtransformed$featureA,
           data_backtransformed$featureB)
  data_check <- data %>% arrange(combi)
  data_bt_check <- data_backtransformed %>% arrange(combi)
  expect_equal(data_check$counts, data_bt_check$counts)
  expect_equal(round(data_check$response, 5),
               round(data_bt_check$response, 5))
  expect_equal(data_check$combi, data_bt_check$combi)
  expect_equal(round(data_check$gam_weights, 5),
               round(data_bt_check$gam_weights, 5))
  expect_equal(round(data_check$gam_offsets, 5),
               round(data_bt_check$gam_offsets, 5))
  expect_equal(round(data_check$Delta, 5), round(data_bt_check$Delta, 5))
  expect_equal(data_check$discrete, data_bt_check$discrete)
  expect_warning(
    data2counts(
      original,
      var_vec = c(2:3),
      y = 1,
      bin_width = 0.5 ,
      values_discrete = c(0, 0.25, 1)
    )
  )
})

# testing errors
test_that("Errors", {
  set.seed(101)
  dt <- data.frame(a = c(2, 3),
                   b = c(2, 3),
                   c = FALSE)
  # invalid type of data
  expect_error(
    data2counts(
      "hhha",
      var_vec = c(2:3),
      y = 1,
      bin_width = 0.01 ,
      values_discrete = discrete_values,
      domain_continuous = continous_domain
    )
  )
  # invalid type of y
  expect_error(data2counts(
    dt,
    var_vec = c(1, 2),
    y = FALSE,
    bin_width = 0.01
  ))
  # invalid type of var_vec
  expect_error(data2counts(
    dt,
    var_vec = FALSE,
    y = 3,
    bin_width = 0.01
  ))
  # invalid type of sample_weights
  expect_error(data2counts(
    dt,
    var_vec = c(1, 2),
    y = 3,
    sample_weights = FALSE,
    bin_width = 0.01
  ))
  # invalid type of bin_width
  expect_error(data2counts(
    dt,
    var_vec = c(1, 2),
    y = 3,
    bin_width = FALSE
  ))
  # invalid type of bin_number
  expect_error(data2counts(
    dt,
    var_vec = c(1, 2),
    y = 3,
    bin_number = FALSE
  ))
  # invalid type of values_discrete
  expect_error(data2counts(
    dt,
    var_vec = c(1, 2),
    y = 3,
    bin_number = 10,
    values_discrete = "hhh"
  ))
  # invalid type of domain_continuous
  expect_error(
    data2counts(
      dt,
      var_vec = c(1, 2),
      y = 3,
      bin_number = 10,
      domain_continuous = "hhh"
    )
  )
  # invalid type of Weights_discrete
  expect_error(data2counts(
    dt,
    var_vec = c(1, 2),
    y = 3,
    bin_number = 10,
    weights_discrete =  "hhh"
  ))
  # invalid column name in var_vec
  expect_error(data2counts(
    dt,
    var_vec = c("nichtDrin"),
    y = 3,
    bin_number = 10
  ))
  # invalid column name in y
  expect_error(data2counts(
    dt,
    var_vec = 1,
    y = c("nichtDrin"),
    bin_number = 10
  ))
  # invalid column name in sample_weights
  expect_error(data2counts(
    dt,
    var_vec = 1,
    y = 2,
    sample_weights = c("nichtDrin"),
    bin_number = 10
  ))
  # domain_continuous and values_discrete both FALSE
  expect_error(
    data2counts(
      dt,
      var_vec = 1,
      y = 2,
      bin_number = 10,
      values_discrete = FALSE,
      domain_continuous = FALSE
    )
  )
  # uncompatible bin_number and bin_width
  expect_error(data2counts(
    dt,
    var_vec = 1,
    y = 2,
    bin_number = 10,
    bin_width = c(22, 2)
  ))
  # uncompatible bin_number and bin_width
  expect_error(data2counts(
    dt,
    var_vec = 1,
    y = 2,
    bin_number = 10,
    bin_width = 22
  ))
  # column index of var_vec out of range
  expect_error(data2counts(
    dt,
    var_vec = 5,
    y = 3,
    bin_number = 10
  ))
  # column index of y out of range
  expect_error(data2counts(
    dt,
    var_vec = 1,
    y = 5,
    bin_number = 10
  ))
  # column index of sample_weights of range
  expect_error(data2counts(
    dt,
    var_vec = 1,
    y = 2,
    sample_weights = 5,
    bin_number = 10
  ))
  # unequal number of values_discrete and weights_discrete
  expect_error(
    data2counts(
      dt,
      var_vec = 1,
      y = 2,
      values_discrete = c(0, 1),
      weights_discrete = c(1, 2, 3)
    )
  )
  # values_discrete out of continuous_domain range
  expect_error(data2counts(
    dt,
    var_vec = 1,
    y = 2,
    values_discrete = c(0, 10)
  ))
  # values of response out of continuous_domain range (data frame)
  expect_error(data2counts(
    dt,
    var_vec = 1,
    y = 2,
   domain_continuous = c(0, 1)
  ))
  # bin_width incompatible with continuous_domain
  expect_error(data2counts(
    dt,
    var_vec = 1,
    y = 2,
    bin_width  = c(19, 10, 10, 10)
  ))
  # non-positive weights_discrete
  expect_error(data2counts(
    dt,
    var_vec = 1,
    y = 2,
    weights_discrete = -20,
    domain_continuous = c(0,10)
  ))
  # values of response out of continuous_domain range (data table)
  expect_error(data2counts(
    as.data.table(dt),
    var_vec = 1,
    y = 2,
    domain_continuous = c(0,1)
  ))
})
