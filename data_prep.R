# define function which bins all features and replaces the actual values with bin centers
bin_features <- function(data, features, n_bins) {
  data_copy <- copy(data)
  for (feature in features) {
    feature_col <- data_copy[,..feature][[1]]
    max_value <- max(feature_col)
    min_value <- min(feature_col)
    step_size <- (max_value - min_value) / n_bins
    grid <- seq(from = min_value, to = max_value, by = step_size)
    feature_hist <- hist(feature_col,
                         breaks = grid, right = FALSE, # changed to FALSE to comply with binning of response
                         plot = FALSE)
    list_of_mids <- feature_hist$mids
    binned_feature_col <- list_of_mids[findInterval(feature_col,
                                                    grid,
                                                    left.open = FALSE, # changed to FALSE to comply with binning of response
                                                    rightmost.closed = TRUE)]
    data_copy[,(feature):=binned_feature_col]
  }

  data_copy
}

# I define a function, which constructs the dataset based on the names of the regressors
construct_dataset <- function(data_grouped, response = list("share"),
                              step_size, weighted = FALSE) {

  # join response and regressors since response is viewed as regressor in poisson density regression
  regressors <- colnames(data_grouped)[2:length(colnames(dta_grouped))]
  print(regressors)
  # construct hists for each unique covariable combination and targets
  dta_targets <- construct_hist(data_grouped, step_size, weighted = weighted)
  # construct regression features
  dta_feat <- construct_features(data_grouped, step_size, regressors)
  # construct complete dataset
  dta_est <- cbind(dta_targets, dta_feat)
  # return dataset
  dta_est
}

# # hist function
# construct_hist <- function(data_grouped, step_size) {
#   # define grid for hist
#   grid_hist <- seq(from = 0, to = 1, by = step_size)
#
#   # construct a list of histograms, where each group in data_grouped gets its own hist
#   num_groups <- max(data_grouped$group_id)
#
#   list_of_mids <- list() # initialize lists
#   list_of_counts <- list()
#
#   for (i in 1:num_groups) {
#     data_group <- data_grouped[group_id == i, share]
#     counts_zero <- sum(data_group == 0)
#     counts_one <- sum(data_group == 1)
#     data_hist <- data_group[data_group > 0 & data_group < 1]
#     if (length(data_hist) == 0) {
#       counts_hist <- rep_len(0, length(grid_hist) - 1)
#       hist <- hist(0.5,
#                    breaks = grid_hist, right = FALSE, # In the paper, we defined the binning intervals to be left-closed
#                    plot = FALSE)
#     }
#     else {
#       hist <- hist(data_hist, breaks = grid_hist, plot = FALSE, right = FALSE) # In the paper, we defined the binning intervals to be left-closed
#       counts_hist <- hist$counts
#     }
#     list_of_mids[[i]] <- c(0, hist$mids, 1)
#     counts <- c(counts_zero, counts_hist, counts_one)
#     list_of_counts[[i]] <- counts
#   }
#
#   dta_targets <- data.table::data.table(counts = unlist(list_of_counts, recursive = FALSE), # unlist to get data set
#                                         share = unlist(list_of_mids, recursive = FALSE))
#   dta_targets
# }

# weighted hist function
construct_hist <- function(data_grouped, step_size, weighted = FALSE) {
  if (!weighted) {
    data_grouped$weighting_factor <- 1
  }
  # define grid for hist
  grid_hist <- seq(from = 0, to = 1, by = step_size)

  # construct a list of histograms, where each group in data_grouped gets its own hist
  num_groups <- max(data_grouped$group_id)

  list_of_mids <- list() # initialize lists
  list_of_counts <- list()

  for (i in 1:num_groups) {
    data_group <- data_grouped[group_id == i, c("share", "weighting_factor")]
    counts_zero <- sum(data_group[share == 0, weighting_factor])
    counts_one <- sum(data_group[share == 1, weighting_factor])
    data_hist <- data_group[data_group$share > 0 & data_group$share < 1]
    if (length(data_hist) == 0) {
      counts_hist <- rep_len(0, length(grid_hist) - 1)
      hist <- hist(0.5,
                   breaks = grid_hist, right = FALSE, # In the paper, we defined the binning intervals to be left-closed
                   plot = FALSE)
    }
    else {
      hist <- weights::wtd.hist(data_hist$share, breaks = grid_hist, plot = FALSE,
                                weight = data_hist$weighting_factor, right = FALSE) # In the paper, we defined the binning intervals to be left-closed
      # hist <- wtd.hist(data_hist$share, breaks = grid_hist, plot = FALSE,
      #                           weight = data_hist$weighting_factor, right = FALSE) # In the paper, we defined the binning intervals to be left-closed
      counts_hist <- hist$counts
    }
    list_of_mids[[i]] <- c(0, hist$mids, 1)
    counts <- c(counts_zero, counts_hist, counts_one)
    list_of_counts[[i]] <- counts
  }

  dta_targets <- data.table::data.table(counts = unlist(list_of_counts, recursive = FALSE), # unlist to get data set
                                        share = unlist(list_of_mids, recursive = FALSE))
  dta_targets
}

# features function
construct_features <- function(data_grouped, step_size, regressors) {

  num_groups <- max(data_grouped$group_id)

  # replicate each covariable combination nbins + 3 times
  nbins = length(seq(from = 0, to = 1, by = step_size)) - 1
  index_vec <- as.vector(sapply(1:num_groups, rep, times = nbins + 2,
                                simplify = "vector"))

  dta_feat <- data_grouped[, head(.SD, 1), # each hist is associated with one covariable combination
                           by = group_id][index_vec][, ..regressors] # -> all counts share same regressor values

  # return features
  dta_feat

}
