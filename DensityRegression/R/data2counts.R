#' From observations to a vector of observed (histogram) counts
#'
#' \code{data2counts} prepares data containing the individual response observations
#' \eqn{y_i} appropriately to be used in \code{\link[mgcv]{gam}} Poisson models
#' by combining all observations of the same conditional distribution (i.e., all
#' observations sharing identical values in all covariates) into a vector of counts
#' via a histogram on \eqn{I\setminus D} and counts on \eqn{D} where \eqn{I} is
#' the interval of the continuous domain and \eqn{D} the set of discrete values.
#' the discrete component of the underlying mixed (continuous/discrete)
#' Bayes Hilbert space \eqn{B^2(\mu) = B^2(\mathcal{Y}, \mathcal{A}, \mu)}{B^2(\mu) = B^2(Y, A, \mu)}.
#' We briefly summarize the approach below in the details. Please see Section 2.3
#' in Maier et al. (2025b) for comprehensive description.
#'
#' @encoding UTF-8
#' @importFrom dplyr "%>%" arrange mutate summarise group_by ungroup cur_group_id across
#' @import data.table
#' @importFrom Rdpack reprompt
#' @importFrom weights wtd.hist
#'
#' @param var_vec Vector of variables of \code{dta} on which the covariate
#' combinations are based. The vector can either contain the variable names as
#' strings or the column positions of the respective variables in \code{dta}.
#' @param y Variable in \code{dta} containing the response observations
#' \eqn{y_i}. Either the variable name can be given as string or the column
#' position of the variable in \code{dta} as integer. If missing (\code{NULL}),
#' if there is a unique column of \code{dta} not specified in \code{var_vec},
#' \code{sample_weights}, \code{counts}, and \code{weighted_counts} (see below),
#' this unique column is used.
#'
#' @inheritParams densreg
#'
#' @return The function returns an object of the class \code{histogram_count_data},
#' which is a \code{\link[data.table]{data.table}} with columns:
#' \itemize{
#' \item \code{counts} - For each by \code{var_vec} defined covariate combination
#' the observed (histogram) counts in the first column.
#' \item \code{weighted_counts} - If \code{sample_weights} is not \code{NULL}:
#' Weighted (histogram) counts for the respective bin/discrete value incorporating
#' the sample weights given by \code{sample_weights}.
#' \item name of variable given to \code{y} - Marks the mid of the respective
#' histogram bin (for values in \eqn{I\setminus D}) or the discrete value. If a
#' mid corresponds to a discrete value, the mid is shifted to the right by
#' \eqn{0.0001} times the minimal distance to the next interval limit OR discrete
#' value so that no mid is exactly corresponding to a discrete value. A warning
#' message is generated in this case.
#' \item names of all variable columns which where specified by \code{var_vec} -
#' These columns contain the values of the respective variables.
#' \item \code{group_id} - ID of each covariate combination.
#' \item \code{gam_weights} - Vector to be passed to argument \code{weights} in
#' \code{\link[mgcv]{gam}} when fitting the Poisson Model, if \code{dta} contains
#' sample weights, see Appendix C of Maier et al. (2023).
#' \item \code{gam_offset} - Negative logarithm of \code{gam_weights} to be used
#' as offset to the predictor of the Poisson Model, if \code{dta} contains sample
#' weights, see Appendix C of Maier et al. (2023).
#' \item \code{Delta} - Width of the histogram bin or weight of the Dirac measure
#' for a discrete value defined by \code{weights_discrete}. The Poisson model uses
#' \code{offset(log(Delta))} to add the necessary additive term in the predictor
#' that includes binwidths/dirac weights into the estimation.
#' \item \code{discrete} - Logical value indicating whether the respective \code{y}
#' is a discrete value in \eqn{D}.
#' }
#' Note that a \code{plot}-method for objects of class \code{histogram_count_data}
#' is available via \code{DensityRegression:::histogram_count_data}, however, it is
#' not exported, since it is not tested/documented appropriately, yet.
#'
#' @author Lea Runge, Eva-Maria Maier
#'
#' @seealso \code{\link{densreg}}
#'
#' @examples
#' set.seed(101)
#'
#' # create data where 0 and 1 are the discrete observations, values
#' # equal 2 are replaced below by drawing from a beta distribution
#'
#' dta <- data.frame(obs_density = sample(0:2, 100, replace = TRUE,
#'                   prob = c(0.15, 0.1, 0.75)),
#'                   covariate1 = sample(c("a", "b"), 100, replace = TRUE),
#'                   covariate2 = sample(c("c", "d"), 100, replace = TRUE),
#'                   sample_weights = runif (100, 0, 2))
#' dta[which(dta$obs_density == 2), ]$obs_density <- rbeta(length(which(dta$obs_density == 2)),
#'                                                         shape1 = 3, shape2 = 3)
#'
#' # Create histogram count dataset for dta using 10 equidistant
#' # bins and default values for continuous domain, discrete
#' # values and discrete weights while considering a mixed case
#' # of continuous and discrete domains. The following function calls are
#' # equivalent:
#'
#' data2counts(dta, var_vec = c("covariate1", "covariate2"), y = "obs_density",
#'            sample_weights = "sample_weights", bin_number = 10)
#' data2counts(dta, var_vec = c(2, 3), y = 1, sample_weights = 4, bin_width = 0.1)
#'
#' # Use the vector bin_width to define non-equidistant bins and
#' # specify with values_discrete and weights_discrete discrete values
#' # and weights besides the default (0,1) and weight 1:
#'
#' data2counts(dta, var_vec = c(2, 3), y = 1, sample_weights = 4,
#'            bin_width = c(0.1, 0.5, 0.4), values_discrete = c(0, 1),
#'            weights_discrete = c(0.5, 2))
#'
#' # The use of "values_discrete=FALSE" refers to histogram data in a
#' # purely continous setting (note that now the observations at 0 and 1 are
#' # counted towards the outer bins):
#'
#' data2counts(dta, var_vec = c(2, 3), y = 1, sample_weights = 4, bin_width = 0.1,
#'            values_discrete = FALSE)
#'
#' # filter data set for only observations valued in discrete domain
#'
#' dta_discrete <- dta[which(dta$obs_density %in% c(0, 1)), ]
#'
#' # The use of "domain_continuous=FALSE" refers to histogram data in a
#' # purely discrete setting:
#'
#' data2counts(dta_discrete, var_vec = c(2, 3), y = 1, sample_weights = 4,
#'            bin_width = 0.1, values_discrete = c(0, 1),
#'            weights_discrete = c(0.5, 2), domain_continuous  = FALSE)
#'
#' # use the optional argument counts to "data2counts" already preprocessed
#' # data and select only one of two variables for the grouping
#' dta_discrete <- unique(dta_discrete[, 1:3])
#' dta_discrete$counts <- sample(0:10, nrow(dta_discrete), replace = TRUE)
#'
#' data2counts(dta_discrete, var_vec = "covariate1", y = "obs_density",
#'            counts = "counts", bin_width = 0.1, values_discrete = c(0, 1),
#'            weights_discrete = c(0.5, 2), domain_continuous  = FALSE)
#'
#' @export
#'
#' @references
#' Maier, E.-M., Fottner, A., Stoecker, A., Okhrin, Y., & Greven, S. (2025b):
#' Conditional density regression for individual-level data.
#' arXiv preprint arXiv:XXXX.XXXXX.

data2counts <- function(dta, var_vec, y = NULL, sample_weights = NULL, counts = NULL,
                       weighted_counts = NULL, bin_width = NULL,
                       bin_number = NULL, values_discrete = c(0, 1),
                       weights_discrete = 1, domain_continuous = c(0, 1)) { # already_formatted=FALSE
  if (is.null(counts)) {
    # check for invalid arguments
    checking_data2counts_1(dta, var_vec, sample_weights, bin_width, bin_number,
             values_discrete, weights_discrete, domain_continuous)
    # convert to standardised format and convert given indice vectors into column name vecors
    if ("response" %in% colnames(dta)) {
      names(dta)[names(dta) == 'response'] <- 'response_'
    }
    if (!is.data.table(dta)) {
      dta <- as.data.table(dta)
    }
    if (all(sapply(var_vec, is.numeric))) {
      var_vec <- colnames(dta)[var_vec]
    }
    if (is.numeric(sample_weights)) {
      sample_weights <- colnames(dta)[sample_weights]
    }
    if (is.numeric(counts)) {
      counts <- colnames(dta)[counts]
    }
    if (is.numeric(weighted_counts)) {
      weighted_counts <- colnames(dta)[weighted_counts]
    }
    if (is.null(y)) {
      y <- setdiff(names(dta), c(var_vec, sample_weights, counts, weighted_counts))
    }
    stopifnot("y must be one variable contained in dta. Can only be NULL, if uniquely determined by sample_weights and covariates used to specify effects." = length(y)==  1)
    checking_data2counts_2(dta, y, domain_continuous)
    if (is.numeric(y)) {
      y <- colnames(dta)[y]
    }
    if (isFALSE(weights_discrete)) {
      weights_discrete <- NULL
    }
    if (isFALSE(values_discrete)) {
      values_discrete <- NULL
    }
    if (is.null(bin_width) && is.null(bin_number)) {
      bin_number <- 100
    }
    if (!isFALSE(domain_continuous) &
        is.null(bin_width) && !is.null(bin_number)) {
      bin_width <-
        (domain_continuous[2] - domain_continuous[1]) / bin_number
    }
    if (!isFALSE(domain_continuous)) {
      if (is.null(bin_number) && length(bin_width) == 1) {
        bin_number <-
          length(seq(
            from = domain_continuous[1],
            to = domain_continuous[2],
            by = bin_width
          )) - 1
      }
      if (is.null(bin_number) && length(bin_width) > 1) {
        bin_number <- length(bin_width)
      }
    }
    # group data by covariate combinations
    keys <- paste(var_vec, collapse = ",")
    selection <- c(y, var_vec, "group_id", sample_weights)
    dta_grouped <- dta[, group_id := .GRP,
                       keyby = keys][, ..selection]
    # standardise name of sample weights
    if (!is.null(sample_weights)) {
      colnames(dta_grouped)[length(colnames(dta_grouped))] <-
        "weighting_factor"
    }
    colnames(dta_grouped)[1] <- "response"
    # construct data appropiately for the poisson model
    # join response and regressors since response is viewed as regressor in poisson
    # density regression
    if (is.null(sample_weights)) {
      dta_grouped$weighting_factor <- 1
    }
    regressors <-
      colnames(dta_grouped)[2:(length(colnames(dta_grouped)) - 1)]
    # if unweighted: use 1 as sample weights
    dta_grouped$unweighting_factor <- 1
    # define grid for hist
    if (!isFALSE(domain_continuous)) {
      if (length(bin_width) == 1) {
        grid_hist <-
          seq(from = domain_continuous[1],
              to = domain_continuous[2],
              by = bin_width)
      } else{
        grid_hist <- c(domain_continuous[1], cumsum(bin_width))
      }
      interval_width <- diff(grid_hist)
    }
    shifted <- FALSE
    # initialize histogram structure
    num_groups <- max(dta_grouped$group_id)
    list_of_mids <- list() # initialize lists
    list_of_counts_weighted <- list()
    list_of_counts_unweighted <- list()
    # create hist-data for all groups
    for (i in 1:num_groups) {
      dta_temp <-
        dta_grouped[group_id == i, c("response", "weighting_factor", "unweighting_factor")]
      counts_discrete_unweighted <- rep(0, length(values_discrete))
      counts_discrete_weighted <- rep(0, length(values_discrete))
      if (!is.null(values_discrete)) {
        # counts for discrete values
        for (j in 1:length(values_discrete)) {
          counts_discrete_weighted[j] <-
            sum(dta_temp[response == values_discrete[j], weighting_factor])
          counts_discrete_unweighted[j] <-
            sum(dta_temp[response == values_discrete[j], unweighting_factor])
        }
      }
      dta_hist <- dta_temp[!dta_temp$response %in% values_discrete]
      if (!isFALSE(domain_continuous)) {
        # if only discrete values considered: hist equals discrete count data
        hist_weighted <-
          weights::wtd.hist(
            dta_hist$response,
            breaks = grid_hist,
            plot = FALSE,
            weight = dta_hist$weighting_factor,
            right = FALSE
          )
        counts_hist_weighted <- hist_weighted$counts
        hist_unweighted <-
          weights::wtd.hist(
            dta_hist$response,
            breaks = grid_hist,
            plot = FALSE,
            weight = dta_hist$unweighting_factor,
            right = FALSE
          )
        counts_hist_unweighted <- hist_unweighted$counts
        # if any hsitogram mid is equal to a discrete value, the mid is shifted
        # slightly (0.0001*min(bin_width/2, distance to next discrete value))
        if (any(hist_weighted$mids %in% values_discrete) &
            all(hist_weighted$mids[hist_weighted$mids %in% values_discrete] < hist_weighted$breaks[-1][hist_weighted$mids %in% values_discrete])) {
          dist2next_tic_or_discrete <-
            rep(Inf, length(hist_weighted$mids))
          # find minimal distance to next interval border or discrete value
          for (j in 1:length(hist_weighted$mids)) {
            closest_values <-
              values_discrete[values_discrete > hist_weighted$mids[j]]
            next_discrete <- min(stats::na.omit(closest_values))
            dist_discrete <- next_discrete - hist_weighted$mids[j]
            dist_tic <- interval_width[j] / 2
            dist2next_tic_or_discrete[j] <-
              min(dist_discrete, dist_tic)
          }
          hist_weighted$mids[hist_weighted$mids %in% values_discrete] <-
            hist_weighted$mids[hist_weighted$mids %in% values_discrete] + 0.0001 * dist2next_tic_or_discrete[hist_weighted$mids %in% values_discrete]
          shifted <- TRUE
        }
        list_of_mids[[i]] <- c(values_discrete, hist_weighted$mids)
        counts_weighted <-
          c(counts_discrete_weighted, counts_hist_weighted)
        list_of_counts_weighted[[i]] <- counts_weighted
        counts_unweighted <-
          c(counts_discrete_unweighted, counts_hist_unweighted)
        list_of_counts_unweighted[[i]] <- counts_unweighted
      } else {
        list_of_mids[[i]] <- c(values_discrete)
        counts_weighted <- c(counts_discrete_weighted)
        list_of_counts_weighted[[i]] <- counts_weighted
        counts_unweighted <- c(counts_discrete_unweighted)
        list_of_counts_unweighted[[i]] <- counts_unweighted
      }
    }
    # construct hists from counts and mids vectors for each unique covariate
    # combination and targets
    dta_targets <-
      data.table::data.table(
        counts = unlist(list_of_counts_unweighted, recursive = FALSE),
        weighted_counts = unlist(list_of_counts_weighted, recursive = FALSE),
        # unlist to get dta set
        response = unlist(list_of_mids, recursive = FALSE)
      )
    # construct regression features dta_grouped
    if (!isFALSE(domain_continuous)) {
      index_vec <-
        as.vector(sapply(
          1:num_groups,
          rep,
          times = bin_number + length(values_discrete),
          simplify = "vector"
        ))
    } else{
      index_vec <-
        as.vector(sapply(
          1:num_groups,
          rep,
          times = length(values_discrete),
          simplify = "vector"
        ))
    }
    dta_feat <-
      dta_grouped[, utils::head(.SD, 1), # each hist is associated with one covariate combination
                  by = group_id][index_vec][, ..regressors]
    # construct complete dtaset
    dta_est <- cbind(dta_targets, dta_feat)
    # We have to specify the argument weights in gam and use an offset in the
    # formula to include the weights properly; For each bin, the weight to
    # pass to gam is the weighted count of observations in this bin, divided by
    # the usual count (i.e., the number) of observations in this bin; the offset
    # that has to be used is -log of the respective weight. See Paper of Maier et al.(2023).
    dta_est[, gam_weights := ifelse(counts != 0, weighted_counts / counts, 1)]
    # If counts == 0, the value of weight can actually be arbitrary, but its log
    # has to exist
    dta_est[, gam_offsets := -log(gam_weights)]
    #}
    # if no sample weights are given, the column of weighted counts (equal to the
    # count coulumn in this case) is removed from the dataset
    if (is.null(sample_weights)) {
      dta_est <- dta_est[, -2]
    }
    # Delta equals bin width/ Dirac weight of discrete value
    if (length(weights_discrete) == 1) {
      #same dirac weight for all discrete values
      dta_est$Delta <- NA
      dta_est$Delta[!dta_est$response %in% values_discrete] <-
        rep(bin_width, max(dta_est$group_id))
      dta_est$Delta[dta_est$response %in% values_discrete] <-
        weights_discrete
    }
    if (length(weights_discrete) > 1) {
      #different dirac weights
      discretes <-
        data.frame(values = values_discrete,
                   weights = weights_discrete,
                   count = 0)
      dta_est$Delta <- NA
      dta_est$Delta[!dta_est$response %in% values_discrete] <-
        rep(bin_width, max(dta_est$group_id))
      for (j in 1:length(discretes$values)) {
        # Delta= bin width for all bins (=non-discrete values)
        dta_est$Delta[is.na(dta_est$Delta) &
                        dta_est$response == discretes$values[j]] <-
          discretes$weights[j]
      }
    }
    if (is.null(weights_discrete)) {
      if (length(bin_width) == 1) {
        dta_est[, Delta := rep(bin_width, max(dta_est$group_id) * bin_number)]
      }
      else{
        dta_est[, Delta := rep(bin_width, max(dta_est$group_id))]
      }
    }
    # sort the data by group_id and within the group by the response variable
    dta_est <- dta_est %>% arrange(group_id, response)
    # add column which indicates if a bin is in continuous or discrete domain
    if (!isFALSE(values_discrete)) {
      dta_est <- dta_est%>%mutate(discrete=response%in%values_discrete)
    }
    # reconversion of changed variable names
    if (!is.null(sample_weights)) {
      colnames(dta_est)[3] <- y
    } else {
      colnames(dta_est)[2] <- y
    }

    if ("response_" %in% colnames(dta_est)) {
      names(dta_est)[names(dta_est) == 'response_'] <- 'response'
    }
    if (isTRUE(shifted)) {
      warning(
        "Some bin mids have been shifted minimally as they correspond to the discrete values.",
        call. = FALSE
      )
    }
  } else { # if !is-null(counts)
    if (!is.null(sample_weights)) {
      warning("Since counts are specified, argument sample_weights is ignored (weighted_counts is not).")
    }
    # if (is.numeric(counts)) {
    #   counts <- colnames(dta)[counts]
    # }
    # if (!is.null(weighted_counts)) {
    #   if (is.numeric(weighted_counts)) {
    #     weighted_counts <- colnames(dta)[weighted_counts]
    #   }
    # }
    if (!is.data.table(dta)) {
      dta <- as.data.table(dta)
    }

    # Delta
    if (!isFALSE(domain_continuous)) {
      cont_values <- unlist(unique(dta[, ..y]))
      if (!isFALSE(values_discrete)) {
        cont_values <- setdiff(cont_values, values_discrete)
      }

      if (is.null(bin_width) & is.null(bin_number)) {
        cont_values_sorted <- sort(cont_values)
        breaks <- c(domain_continuous[1],
                     sapply(seq_along(cont_values_sorted)[-1],
                            function(i) cont_values_sorted[i - 1] +
                              (cont_values_sorted[i] - cont_values_sorted[i - 1]) / 2),
                     domain_continuous[2])
        Delta <- diff(breaks) # [order(cont_values)]
        # breaks <- c(domain_continuous[1], sort(cont_valuessort(, domain_continuous[2])
        # diffs <- diff(breaks)
        # Delta <- c(diffs[1] + 0.5 * diffs[2])
        # for (i in 2:(length(cont_values)-1)) {
        #   Delta <- append(Delta, 0.5 * diffs[i] + 0.5 * diffs[i+1])
        # }
        # Delta <- append(Delta, diffs[i+1] * 0.5 + diffs[i+2])
        if (!isFALSE(values_discrete)) {
          ordered_values <- order(c(cont_values, values_discrete))
          if (length(weights_discrete)==1) {
            weights_discrete <- rep(weights_discrete,length(values_discrete))
          }
          Delta <- c(Delta, weights_discrete)[ordered_values]
        }
      } else {
        if (length(bin_width) == 1) {
          Delta <- rep(bin_width, length(cont_values))
        } else {
          Delta <- bin_width
        }
        if (!"weighted_counts" %in% colnames(dta)) {
          # For histogram construction, counts weight the respective observations,
          # i.e., we "abuse" the argument sample_weights for construction
          dta_tmp <- Recall(dta, y=y, var_vec=var_vec, bin_number=bin_number,
                            bin_width= bin_width, # already_formatted = FALSE,
                            sample_weights = counts, domain_continuous = domain_continuous,
                            values_discrete = values_discrete)
          # The real counts of interest are now in the column "weighted_counts",
          # while the first column, "counts", is abundant, and "gam_weights"
          # and "gam_offsets" should be the neutral values (without influence
          # on the estimation), 1 and 0
          dta_est <- dta_tmp[,-1]
          dta_est$gam_weights <- 1
          dta_est$gam_offsets <- 0
          colnames(dta_est)[1] <- "counts"
          return(dta_est)
        } else {
          dta_tmp_1 <- Recall(dta, y=y, var_vec=var_vec, bin_number=bin_number,
                              bin_width= bin_width, sample_weights = counts,
                              domain_continuous = domain_continuous, values_discrete = values_discrete)
          dta_tmp_2 <- Recall(dta, y=y, var_vec=var_vec, bin_number=bin_number,
                              bin_width= bin_width, sample_weights = weighted_counts,
                              domain_continuous = domain_continuous, values_discrete = values_discrete)
          dta_est <- dta_tmp_2
          dta_est[,1] <- dta_tmp_1[,2]
          dta_est$gam_weights <- ifelse(dta_est$counts != 0, dta_est$weighted_counts / dta_est$counts, 1)
          dta_est$gam_offsets <-  -log(dta_est$gam_weights)
          return(dta_est)
        }

        # if (!isFALSE(values_discrete))
        # { ordered_values <- order(c(cont_values, values_discrete))
        # if (length(weights_discrete)==1) {
        #   weights_discrete <- rep(weights_discrete,length(values_discrete))
        # }
        # Delta <- c(Delta, weights_discrete)[ordered_values]
        # }
      }

      if (!isFALSE(values_discrete)) {
        ordered_values <- order(c(cont_values, values_discrete))
        discrete <- c(rep(FALSE, length(cont_values)), rep(TRUE,length(values_discrete)))[ordered_values]
      } else {
        discrete <- c(rep(FALSE, length(cont_values)))
      }
    } else {
      if (length(weights_discrete)==1) {
        weights_discrete <- rep(weights_discrete,length(values_discrete))
      }
      Delta <- weights_discrete
      discrete <- c(rep(TRUE, length(values_discrete)))
    }

    if (is.numeric(y)) {
      y <- colnames(dta)[y]
    }
    if (is.numeric(var_vec)) {
      var_vec <- colnames(dta)[var_vec]
    }

    if (!"weighted_counts"%in% colnames(dta)) {
      dta_est <- cbind(dta$counts,dta[,..y], dta[,..var_vec])
      colnames(dta_est)[1] <- "counts"
      keys <- paste(var_vec, collapse = ",")
      selection <- c(y, var_vec, "group_id")
      dta_est <- dta_est %>%
        group_by(across(all_of(c(y,var_vec)))) %>%
        summarise(counts = sum(counts), .groups = 'drop')

      dta_est <- dta_est %>%
        group_by(across(all_of(c(var_vec)))) %>%
        mutate(group_id = cur_group_id()) %>%
        ungroup()
      dta_est <- dta_est[order(dta_est$group_id),]


      dta_est$Delta <- rep(Delta, max(dta_est$group_id))
      dta_est$discrete <- rep(discrete, max(dta_est$group_id))

      dta_est$gam_weights <- 1
      dta_est$gam_offsets <- 0
      l <- length(var_vec)
      dta_est <- dta_est[,c(l+2,1,2:(l+1),l+3,l+6,l+7,l+4,l+5)]
    } else {
      dta_est <- cbind(dta$counts,dta$weighted_counts,dta[,..y], dta[,..var_vec])
      colnames(dta_est)[1] <- "counts"
      colnames(dta_est)[2] <- "weighted_counts"
      keys <- paste(var_vec, collapse = ",")
      selection <- c(y, var_vec, "group_id")
      dta_est <- dta_est %>%
        group_by(across(all_of(c(y,var_vec)))) %>%
        summarise(counts = sum(counts),weighted_counts = sum(weighted_counts), .groups = 'drop')

      dta_est <- dta_est %>%
        group_by(across(all_of(c(var_vec)))) %>%
        mutate(group_id = cur_group_id()) %>%
        ungroup()
      dta_est <- dta_est[order(dta_est$group_id),]


      dta_est$Delta <- rep(Delta, max(dta_est$group_id))
      dta_est$discrete <- rep(discrete, max(dta_est$group_id))

      dta_est$gam_weights <- ifelse(dta_est$counts != 0, dta_est$weighted_counts / dta_est$counts, 1)
      dta_est$gam_offsets <-  -log(dta_est$gam_weights)
      l <- length(var_vec)
      dta_est <- dta_est[,c(l+2,l+3,1,2:(l+1),l+4,l+7,l+8,l+5,l+6)]
    }
  }
  dta_est[] <- dta_est
  attr(dta_est, "class") <- c("histogram_count_data", class(dta_est))
  if (length(attr(dta_est, "class"))==4) {
    attr(dta_est, "class") <- attr(dta_est, "class")[c(1,4)]
  }
  return(dta_est)
}
################################################################################

#' Check for validity of parameters of \code{\link{data2counts}}
#'
#' These functions check validity of parameters for \code{\link{data2counts}}.
#' Returns an error if any parameters requirement is violated.
#'
#' @noRd
#'
#' @encoding UTF-8
#' @importFrom dplyr "%>%" arrange
#' @import data.table
#' @importFrom Rdpack reprompt
#'
#' @inheritParams data2counts

checking_data2counts_1 <- function(dta, var_vec, # y,
                        sample_weights, bin_width, bin_number,
                        values_discrete, weights_discrete, domain_continuous) {

  if (!is.data.table(dta) & !is.data.frame(dta)) {
    stop("Invalid type of argument 'dta'! data.frame or data.table is required")
  }
  # if (!is.null(y)) {
  #   if (!is.numeric(y) & !is.character(y)) {
  #     stop("Invalid type of argument 'y'! character or numeric is required")
  #   }
  # }
  if (!is.numeric(var_vec) & !is.character(var_vec)) {
    stop("Invalid type of argument 'var_vec'! Vector of type character or numeric is required")
  }
  if (!is.numeric(sample_weights) &
      !is.character(sample_weights) & !is.null(sample_weights)) {
    stop("Invalid type of argument 'sample_weights'! character or numeric is required")
  }
  if (!is.numeric(bin_width) & !is.null(bin_width)) {
    stop("Invalid type of argument 'bin_width'! numeric is required")
  }
  if (!is.numeric(bin_number) & !is.null(bin_number)) {
    stop("Invalid type of argument 'bin_number'! numeric is required")
  }
  if ((
    !is.numeric(values_discrete) &
    !is.null(values_discrete) & !isFALSE(values_discrete)
  )) {
    stop("Invalid type of argument 'values_discrete'! numeric is required")
  }
  if ((
    !is.numeric(weights_discrete) &
    !is.null(weights_discrete) & !isFALSE(weights_discrete)
  )) {
    stop("Invalid type of argument 'weights_discrete'! numeric is required")
  }
  if ((
    !is.numeric(domain_continuous) &
    !is.null(domain_continuous) & !isFALSE(domain_continuous)
  )) {
    stop("Invalid type of argument 'domain_continuous'! numeric is required")
  }
  if (!is.numeric(var_vec) & !all(var_vec %in% colnames(dta))) {
    stop("Invalid variable names!")
  }
  # if (!is.numeric(y) & !all(y %in% colnames(dta))) {
  #   stop("Invalid name for variable of density definition!")
  # }
  if (!is.numeric(sample_weights) &
      !all(sample_weights %in% colnames(dta))) {
    stop("Invalid name for sample weights variable!")
  }
  if (isFALSE(values_discrete) & isFALSE(domain_continuous)) {
    stop("Discrete and continuous components are both empty!")
  }
  if (!is.null(bin_width) & !is.null(bin_number)) {
    if (length(bin_width) != bin_number & length(bin_width) > 1) {
      stop("Number of bins and bin width not compatible!")
    }
  }
  if (!isFALSE(domain_continuous)) {
    if (!is.null(bin_width) &
        !is.null(bin_number) &
        length(bin_width) == 1 &
        all((domain_continuous[2] - domain_continuous[1]) / bin_number != bin_width)) {
      stop("Number of bins and bin width not compatible!")
    }
  }
  if (is.numeric(var_vec) & max(var_vec) > length(colnames(dta))) {
    stop("Invalid range of given variables!")
  }
  if (!is.null(sample_weights)) {
    if (is.numeric(sample_weights) &
        max(sample_weights) > length(colnames(dta))) {
      stop("Invalid range of given sample weight variable!")
    }
  }
  # if (is.numeric(y) & max(y) > length(colnames(dta))) {
  #   stop("Invalid range of given density variable!")
  # }
  if (!isFALSE(values_discrete) &
      !isFALSE(weights_discrete) &
      length(weights_discrete) != 1 &
      length(values_discrete) != length(weights_discrete)) {
    stop("Unequal number of discrete values and discrete weights!")
  }
  if (!isFALSE(domain_continuous) &
      length(bin_width) > 1 &
      sum(bin_width) != domain_continuous[2]) {
    stop("bin width vector is not compatible with continuous domain!")
  }
  if (!isFALSE(domain_continuous) & !isFALSE(values_discrete)) {
    if (max(values_discrete) > max(domain_continuous) |
        min(values_discrete) < min(domain_continuous)) {
      stop("Discrete values out of range of the continuous domain!")
    }
  }
  # if (!is.data.table(dta)&!isFALSE(domain_continuous)&is.data.frame(dta)) {
  #   if (max(dta[,y]) > max(domain_continuous) |
  #       min(dta[,y]) < min(domain_continuous)) {
  #     stop("Values of response out of range of the continuous domain!")
  #   }
  # }
  # if (is.data.table(dta)&!isFALSE(domain_continuous)&is.data.frame(dta)) {
  #   if (max(dta[,..y]) > max(domain_continuous) |
  #       min(dta[,..y]) < min(domain_continuous)) {
  #     stop("Values of response out of range of the continuous domain!")
  #   }
  # }
  if (!isFALSE(weights_discrete)) {
    if (!all(weights_discrete > 0)) {
      stop("Weights for discrete values have to be positive!")
    }
  }
}

#' @describeIn checking_data2counts_1 Further parameter checks for \code{\link{data2counts}}
checking_data2counts_2 <- function(dta, y, domain_continuous) {
  if (!is.numeric(y) & !is.character(y)) {
    stop("Invalid type of argument 'y'! character or numeric is required")
  }

  if (!is.numeric(y) & !all(y %in% colnames(dta))) {
    stop("Invalid name for variable of density definition!")
  }

  if (is.numeric(y) & max(y) > length(colnames(dta))) {
    stop("Invalid range of given density variable!")
  }

  if (!is.data.table(dta)&!isFALSE(domain_continuous)&is.data.frame(dta)) {
    if (max(dta[,y]) > max(domain_continuous) |
        min(dta[,y]) < min(domain_continuous)) {
      stop("Values of response out of range of the continuous domain!")
    }
  }
  if (is.data.table(dta)&!isFALSE(domain_continuous)&is.data.frame(dta)) {
    if (max(dta[,..y]) > max(domain_continuous) |
        min(dta[,..y]) < min(domain_continuous)) {
      stop("Values of response out of range of the continuous domain!")
    }
  }
}
