################################################################################
######################## Functions to get model effects ########################
################################################################################

# get_single_effects computes the model effect matrix from a specific effect
# (group specific interceps) from the predicted terms of a model estimated with mgcv::gam
# Arguments:
# pred_terms: object created via predict(..., type = "terms")
# positions: column positions of the effects of interest in pred_terms
# G: number of bins
# Return:
# list of vectors named with the respective effect


get_single_effects <- function(pred_terms, positions, G = 400) {
  estimated_effects = list()
  for (i in c(1:length(positions))) {
    effect <- matrix(pred_terms[, positions[i]], nrow = G)
    zero_cols <- list()
    for (j in c(1:ncol(effect))) {
      if (all(effect[, j] == 0)) {
        zero_cols <- append(zero_cols, j)
      }
    }
    if (!is.null(unlist(zero_cols))) {
      effect <- effect[, -unlist(zero_cols)]
    }
    effect <- effect[, 1]
    estimated_effects <- append(estimated_effects, list(effect))
  }
  names(estimated_effects) <- colnames(pred_terms)[positions]
  names(estimated_effects)[1] <- "intercept"
  return(estimated_effects)
}


# get_single_effects computes the model effect matrix from a specific effect
# (group specific interceps) from the predicted terms of a model estimated with mgcv::gam
# Arguments:
# pred_terms: object created via predict(..., type = "terms")
# positions: column positions of the effects of interest in pred_terms
# G: number of bins
# cols_in_origData: vector column numbers in origData of the included covariates,
#                   for interaction effects the last element of the vector is "inter"
#                   for functional vary. coeff. the last element of the vector is "cont"
# origData: histogram data on which the model is based
# Return:
# list of vectors named with the respective effect

get_smooth_effects <- function(pred_terms, positions, G = 400, cols_in_origData,
                               origData) {
    estimated_effects <- list()
    if (is.list(cols_in_origData)) {
      j <- 1
      for (i in c(1:length(cols_in_origData))) {
        if (length(cols_in_origData[[i]]) == 1) {
          effect <- matrix(pred_terms[, positions[j]], nrow = G)
          groupID_index <-
            unique(origData %>% dplyr::select(group_id, cols_in_origData[[i]]))
          groupID_index$ind <- c(1:ncol(effect))
          uniq_ind <- unique(groupID_index %>% dplyr::select(2, ind))
          uniq_ind <- uniq_ind[!duplicated(uniq_ind[, 1]), ]
          uniq_ind <- uniq_ind[order(uniq_ind[, 1]), ]
          effect <- effect[, uniq_ind$ind]
          estimated_effects <- append(estimated_effects, list(effect))
          j <- j + 1
        }
        if (length(cols_in_origData[[i]]) == 2) {
          n_levels <- nrow(unique(origData %>% dplyr::select(cols_in_origData[[i]][2])))
          for (k in (2:n_levels)) {
            x <- cols_in_origData[[i]][2]
            lev <- origData %>% dplyr::select(all_of(x)) %>% sapply(levels)
            effect <- matrix(pred_terms[, positions[j]], nrow = G)
            groupID_index <-
              unique(origData %>% dplyr::select(group_id, cols_in_origData[[i]][1], cols_in_origData[[i]][2]))
            groupID_index$ind <- c(1:ncol(effect))
            groupID_index <- groupID_index %>% dplyr::filter(.[[3]] == lev[k])
            uniq_ind <- unique(groupID_index %>% dplyr::select(2, ind))
            uniq_ind <- uniq_ind[!duplicated(uniq_ind[, 1]), ]
            uniq_ind <- uniq_ind[order(uniq_ind[, 1]), ]
            effect <- effect[, uniq_ind$ind]
            estimated_effects <- append(estimated_effects, list(effect))
            j <- j + 1
          }
        }
        if (length(cols_in_origData[[i]]) >= 3) {
          if (cols_in_origData[[i]][3] == "cont")
          {
            n_levels <-
              nrow(unique(origData %>% dplyr::select(
                as.numeric(cols_in_origData[[i]][2])
              )))
            x <- as.numeric(cols_in_origData[[i]][2])
            lev <- origData %>% dplyr::select(all_of(x)) %>% sapply(unique)
            effect <- matrix(pred_terms[, positions[j]], nrow = G)
            to_select <- c()
            for (k in (1:n_levels)) {
              groupID_index <-
                unique(origData %>% dplyr::select(
                  group_id,
                  as.numeric(cols_in_origData[[i]][1]),
                  as.numeric(cols_in_origData[[i]][2])
                ))
              groupID_index$ind <- c(1:ncol(effect))
              groupID_index <- groupID_index %>% dplyr::filter(.[[3]] == lev[k])
              uniq_ind <- unique(groupID_index %>% dplyr::select(2, ind))
              uniq_ind <- uniq_ind[!duplicated(uniq_ind[, 1]), ]
              uniq_ind <- uniq_ind[order(uniq_ind[, 1]), ]
              to_select <- append(to_select, uniq_ind$ind)
            }
            effect <- effect[, to_select]
            estimated_effects <-
              append(estimated_effects, list(effect))
            j <- j + 1
          } else{
            if(cols_in_origData[[i]][length(cols_in_origData[[i]])]=="inter"){
            effect <- matrix(pred_terms[, positions[j]], nrow = G)
            to_select <- c()
            inds <-
              c(as.numeric(cols_in_origData[[i]][c(1:(length(cols_in_origData[[i]]) -
                                                        1))]))
            uniq <-
              origData[!duplicated(origData[, ..inds]),]
            for (covCol in rev(inds)) {
              uniq <- uniq[order(uniq[, ..covCol]), ]
            }
            uniq_groups <- unique(uniq$group_id)
            to_select <- append(to_select, uniq_groups)
            effect <- effect[, to_select]
            estimated_effects <-
              append(estimated_effects, list(effect))
            j <- j + 1
            }
            else{
              n_levels <- nrow(unique(origData %>% dplyr::select(as.numeric(cols_in_origData[[i]][length(cols_in_origData[[i]])-1]))))
              for (k in (2:n_levels)) {

                x <- as.numeric(cols_in_origData[[i]][length(cols_in_origData[[i]])-1])
                lev <- origData %>% dplyr::select(all_of(x)) %>% sapply(levels)
                effect <- matrix(pred_terms[, positions[j]], nrow = G)
                groupID_index <-
                  unique(origData %>% dplyr::select(group_id, as.numeric(cols_in_origData[[i]][1:(length(cols_in_origData[[i]])-1)])))
                groupID_index$ind <- c(1:ncol(effect))
                groupID_index <- groupID_index %>% dplyr::filter(.[[length(cols_in_origData[[i]])]] == lev[k])
                to_select <- c()
                for (covCol in rev(2:(length(cols_in_origData[[i]]) -
                                  1))) {
                  groupID_index <- groupID_index[order(groupID_index[, ..covCol]), ]
                }
                uniq_groups <- unique(groupID_index$group_id)
                to_select <- append(to_select, uniq_groups)
                effect <- effect[, to_select]
                estimated_effects <-
                  append(estimated_effects, list(effect))
                j <- j + 1

              }
            }
          }
        }
      }
    } else {
      for (i in c(1:length(positions))) {
        effect <- matrix(pred_terms[, positions[i]], nrow = G)
        groupID_index <-
          unique(origData %>% dplyr::select(group_id, cols_in_origData[i]))
        groupID_index$ind <- c(1:ncol(effect))
        uniq_ind <- unique(groupID_index %>% dplyr::select(2, ind))
        uniq_ind <- uniq_ind[!duplicated(uniq_ind[, 1]), ]
        uniq_ind <- uniq_ind[order(uniq_ind[, 1]), ]
        effect <- effect[, uniq_ind$ind]
        estimated_effects <- append(estimated_effects, list(effect))
      }
    }
    names(estimated_effects) <- colnames(pred_terms)[positions]
    return(estimated_effects)
  }



# get_estimated_model_effects calculates the model effects (Bayes- and clr-level).
# Arguments:
# model: gam model object
# pred_terms: optional if already predicted
# G: number of bins for continuous component
# positions_singles: column positions of the group-specific intercepts of interest in pred_terms
# positions_smooths: column positions of the smooth effects of interest in pred_terms
# smooth_cols: vector column numbers in origData of the in the smooth effects included covariates,
#                   for interaction effects the last element of the vector is "inter"
#                   for functional vary. coeff. the last element of the vector is "cont"
# origData: histogram data
# positions_random_effects: indices of random effects in pred_terms, if non included: NA
# domain_continuous
# values_discrete=FALSE
# weights_discrete=FALSE
# Return:
# List of estimated effects, named after the respective terms


get_estimated_model_effects <-
  function(model,
           pred_terms = NULL,
           G = 400,
           positions_singles,
           positions_smooths,
           smooth_cols,
           origData,
           positions_random_effects = -1,
           domain_continuous,
           values_discrete = FALSE,
           weights_discrete = FALSE) {
    # Package to fit functional regression models
    if (is.null(pred_terms)) {
      pred_terms <- stats::predict(model, type = "terms")
    }
    if (isFALSE(values_discrete)) {
      values_discrete <- NULL
    }
    if (isFALSE(domain_continuous)) {
      G <- 0
    }
    num_obs <- length(values_discrete) + G
    stopifnot(
      "Given number of bins for continuous component is not compatible with given model" = nrow(pred_terms) %% (G +
                                                                                                                  length(values_discrete)) == 0
    )
    if (!is.null(positions_singles)) {
      singles <-
        get_single_effects(pred_terms, positions_singles, G + length(values_discrete))
    } else{
      singles <- NULL
    }
    if (!is.null(positions_smooths)) {
      smooths <-
        get_smooth_effects(
          pred_terms,
          positions_smooths,
          G + length(values_discrete),
          smooth_cols,
          origData
        )
    }
    else{
      smooths <- NULL
    }
    estimated_effects <- c(singles, smooths)
    if (!is.na(positions_random_effects))
    {
      if (positions_random_effects == -1) {
        names(estimated_effects)[length(estimated_effects)] <-
          paste("random", names(estimated_effects)[length(estimated_effects)], sep =
                  "_")
      }
      else{
        positions_random_effects - 1
        names(estimated_effects)[positions_random_effects] <-
          paste("random", names(estimated_effects)[positions_random_effects], sep =
                  "_")
      }
    }
    #estimated_effects_PLUS_intercept<-estimated_effects
    #for (i in c(2:length(positions_singles))){
    # estimated_effects_PLUS_intercept[[i]]<-estimated_effects_PLUS_intercept[[i]]+estimated_effects_PLUS_intercept[[1]]
    #}
    #if (!is.null(positions_smooths)){
    #   for (i in c((length(positions_singles)+1):length(estimated_effects))){
    #intc<-replicate(ncol(estimated_effects_PLUS_intercept[[i]]),estimated_effects_PLUS_intercept[[1]])
    #estimated_effects_PLUS_intercept[[i]]<-estimated_effects_PLUS_intercept[[i]]+estimated_effects_PLUS_intercept[[1]]
    # } }
    pdf_estimated_effects <- estimated_effects
    for (i in c(1:length(pdf_estimated_effects))) {
      if (is.null(values_discrete)) {
        if (!is.null(ncol(pdf_estimated_effects[[i]]))) {
          pdf_estimated_effects[[i]] <-
            apply(
              pdf_estimated_effects[[i]],
              MARGIN = 2,
              FUN = clr,
              inverse = TRUE,
              w = origData$Delta[1:num_obs]
            )

        } else{
          pdf_estimated_effects[[i]] <-
            clr(pdf_estimated_effects[[i]],
                inverse = TRUE,
                w = origData$Delta[1:num_obs])

        }
      }
      if (!is.null(values_discrete) & !isFALSE(domain_continuous)) {
        interval_width <-
          (domain_continuous[2] - domain_continuous[1]) / G
        t <- seq(
          from = domain_continuous[1] + 1 / 2 * interval_width,
          to = domain_continuous[2] - 1 / 2 * interval_width,
          by = interval_width
        )
        t <- sort(c(t, values_discrete))
        discrete_ind <- match(values_discrete, t)
        if (!is.null(ncol(pdf_estimated_effects[[i]]))) {
          pdf_estimated_effects[[i]] <-
            apply(
              pdf_estimated_effects[[i]],
              MARGIN = 2,
              FUN = clr,
              inverse = TRUE,
              w = origData$Delta[1:num_obs]
            )
        } else{
          pdf_estimated_effects[[i]] <-
            clr(pdf_estimated_effects[[i]],
                inverse = TRUE,
                w = origData$Delta[1:num_obs])
        }
      }
      if (!is.null(values_discrete) & isFALSE(domain_continuous))
      {
        if (!is.null(ncol(pdf_estimated_effects[[i]]))) {
          pdf_estimated_effects[[i]] <-
            apply(
              pdf_estimated_effects[[i]],
              MARGIN = 2,
              FUN = clr,
              inverse = TRUE,
              w = origData$Delta[1:num_obs]
            )

        } else{
          pdf_estimated_effects[[i]] <-
            clr(pdf_estimated_effects[[i]],
                inverse = TRUE,
                w = origData$Delta[1:num_obs])

        }
      }
    }
    # pdf_estimated_effects_PLUS_intercept <-estimated_effects_PLUS_intercept
    # for(i in c(1:length(pdf_estimated_effects_PLUS_intercept))){
    #   if (is.null(values_discrete)) {
    #     if(!is.null(ncol(pdf_estimated_effects[[i]]))){
    #       pdf_estimated_effects_PLUS_intercept[[i]]<-apply(pdf_estimated_effects_PLUS_intercept[[i]],MARGIN=2,FUN=clr,inverse = TRUE, w=origData$Delta[1:num_obs])
    #     }else{
    #       pdf_estimated_effects_PLUS_intercept[[i]]<-clr(pdf_estimated_effects_PLUS_intercept[[i]],inverse = TRUE, w=origData$Delta[1:num_obs])
    #     }}
    #   if (!is.null(values_discrete)&!isFALSE(domain_continuous)){
    #     interval_width <-
    #       (domain_continuous[2] - domain_continuous[1]) / G
    #     t <- seq(
    #       from = domain_continuous[1] + 1 / 2 * interval_width,
    #       to = domain_continuous[2] - 1 / 2 * interval_width,
    #       by = interval_width
    #     )
    #     t <- sort(c(t, values_discrete))
    #     discrete_ind <- match(values_discrete, t)
    #     if(!is.null(ncol(pdf_estimated_effects[[i]]))){
    #       pdf_estimated_effects_PLUS_intercept[[i]]<-apply(pdf_estimated_effects_PLUS_intercept[[i]],MARGIN=2,FUN=clr,inverse = TRUE, w=origData$Delta[1:num_obs])
    #     }else{
    #       pdf_estimated_effects_PLUS_intercept[[i]]<-clr(pdf_estimated_effects_PLUS_intercept[[i]],inverse = TRUE, w=origData$Delta[1:num_obs])
    #     }
    #   }
    #   if(!is.null(values_discrete)&isFALSE(domain_continuous))
    #   { if(!is.null(ncol(pdf_estimated_effects[[i]]))){
    #     pdf_estimated_effects_PLUS_intercept[[i]]<-apply(pdf_estimated_effects_PLUS_intercept[[i]],MARGIN=2,FUN=clr,inverse = TRUE, w=origData$Delta[1:num_obs])
    #   }else{
    #     pdf_estimated_effects_PLUS_intercept[[i]]<-clr(pdf_estimated_effects_PLUS_intercept[[i]],inverse = TRUE, w=origData$Delta[1:num_obs])
    #   }
    #   }
    # }
    return(
      list(
        G = G ,
        estimated_effects = estimated_effects,
        #estimated_effects_PLUS_intercept=estimated_effects_PLUS_intercept,
        pdf_estimated_effects = pdf_estimated_effects
      )
    )
    # pdf_estimated_effects_PLUS_intercept=pdf_estimated_effects_PLUS_intercept))
  }



#' PLOT FOR CONTINUOUS MODELS
#'
#' plot_densities is used to plot response densities and predictions (Bayes- and
#' clr-level) for our model.
#'
#' @noRd
#'
#' @encoding UTF-8
#'
#' @param dens_matrix matrix with each column containing the density values of one covariable combination
#' @param pdf boolean value indicating if the plot shall be saved (as pdf)
#' @param name, width, height, path indicating the name (with the respective path),
#'   width and height of the file, if pdf = TRUE;
#'   if NULL, name is created automatically out of the name of the matrix given to
#'   dens_matrix
#' @param color, ylim, lty, main, ylab specify the color gradient (each year corresponds
#'   to one color), ylim, the line type, the main, and the y-axis limits of the plot
#' @param G position of the grid points, where the densities are evaluated
#' domain_continuous
#' @param breaks_x, breaks_y
#' @param xlab, ylab names for x and y axis
#' @param single boolean if all curves should be plotted in one plot
#' @param axes, axis_labels if single=FALSE, should every plot have its own x-axis and label
#' @param breaks_x, breaks_y breaks for the axes
#' @param legend_names names for the different levels in the legend
#' @param legend_title legend title
#' @param monochrome boolean if all curves should be plotted in the same color if single=TRUE
#' @param ncol number of columns if single=FALSE for the plot grid
#' @param ... Further arguments passed to matplot/plot (e.g. las, etc.)

plot_densities <- function(dens_matrix,
                           pdf = FALSE,
                           name = NULL,
                           width = NA,
                           height = NA,
                           unit = c("in", "cm", "mm", "px"),
                           path = "./Images/",
                           color = NULL,
                           ylim = NULL,
                           xlim = NULL,
                           main = NULL,
                           space = 0.3,
                           breaks_x = waiver(),
                           breaks_y = waiver(),
                           ylab = NULL,
                           G = 100,
                           domain = c(0, 1),
                           single = FALSE,
                           axes = "all_x",
                           axis_labels = "all_x",
                           xlab = "s",
                           legend_names = NULL,
                           legend_title = NULL,
                           monochrome = FALSE,
                           ncol = 3,
                           ...) {
  if (pdf && is.null(name)) {
    stop("File name has to be given if pdf is TRUE.")
  }
  interval_width <- (domain[2] - domain[1]) / G
  t <- seq(
    from = domain[1] + 1 / 2 * interval_width,
    to = domain[2] - 1 / 2 * interval_width,
    by = interval_width
  )
  n <- max(ncol(dens_matrix), 1)
  melted <- reshape2::melt(dens_matrix)
  if (ncol(melted) == 1) {
    melted$variable <- 1
  }
  if (ncol(melted) > 2) {
    melted <- reshape2::melt(dens_matrix)[, c(2, 3)]
    colnames(melted)[1] <- "variable"
  }
  melted$t <- rep(t, n)
  melted$variable <- as.factor(melted$variable)
  if (!is.null(legend_names)) {
    levels(melted$variable) <- legend_names
  }
  if (monochrome) {
    mcol <- color
    if (is.null(mcol)) {
      mcol <- "black"
    }
    color <- scale_color_manual(values = (rep(mcol, n)))
  }
  if (is.null(xlim)) {
    xlim <- domain
  }
  if (is.null(ylim)) {
    ylim <- c(min(dens_matrix, na.rm = TRUE),
              max(dens_matrix, na.rm = TRUE))
  }
  if (single) {
    plt <-
      ggplot(data = melted, aes(
        x = t,
        y = value,
        group = variable,
        color = variable
      )) +
      geom_line() +
      color +
      ylab(ylab) + xlab(xlab) +
      guides(color = guide_legend(title = legend_title)) + ggtitle(main) +
      theme(plot.title = element_text(size = 14, hjust = 0.5), ...) +
      scale_x_continuous(breaks = breaks_x, limits = xlim) +
      scale_y_continuous(breaks = breaks_y, limits = ylim)
  } else{
    plt <-
      ggplot(data = melted, aes(
        x = t,
        y = value,
        group = variable,
        color = variable
      )) +
      geom_line() +
      theme(legend.position = "none") +

      theme(
        panel.spacing = unit(space, "lines"),
        strip.text.x = element_text(size = 8),
        plot.title = element_text(size = 14, hjust = 0.5)
      ) +
      color +
      facet_wrap( ~ variable,
                  axes = axes,
                  axis.labels = axis_labels,
                  ncol = ncol) +
      ylab(ylab) + xlab(xlab) + ggtitle(main) +
      scale_x_continuous(breaks = breaks_x, limits = xlim) +
      scale_y_continuous(breaks = breaks_y, limits = ylim)
  }
  if (pdf) {
    ggsave(
      name,
      plot = plt,
      path = path,
      create.dir = TRUE,
      device = "pdf",
      height = height ,
      width = width,
      units = unit
    )
  }
  return(plt)
}




# PLOT FOR MIXED MODELS
# plot_densities_mixed is used to plot response densities and predictions (Bayes- and
# clr-level) for our model.
# Arguments:
# dens_matrix: matrix with each column containing the density values of one covariable combination
# pdf: boolean value indicating if the plot shall be saved (as pdf)
# name, width, height, path: indicating the name (with the respective path),
#   width and height of the file, if pdf = TRUE;
#   if NULL, name is created automatically out of the name of the matrix given to
#   dens_matrix
# color, ylim, lty, main, ylab: specify the color gradient (each year corresponds
#   to one color), ylim, the line type, the main, and the y-axis limits of the plot
# G: position of the grid points, where the densities are evaluated
# domain_continuous
# values_discrete
# breaks_x, breaks_y
# xlab, ylab: names for x and y axis
# single: boolean if all curves should be plotted in one plot
# axes, axis_labels: if single=FALSE, should every plot have its own x-axis and label
# breaks_x, breaks_y: breaks for the axes
# legend_names: names for the different levels in the legend
# legend_title: legend title
# monochrome: boolean if all curves should be plotted in the same color if single=TRUE
# ncol: number of columns if single=FALSE for the plot grid
# ...: Further arguments passed to matplot/plot (e.g. las, etc.)

plot_densities_mixed <-
  function(dens_matrix,
           pdf = FALSE,
           name = NULL,
           width = NA,
           height = NA,
           unit = c("in", "cm", "mm", "px"),
           path = "./Images/",
           color = NULL,
           ylim = NULL,
           xlim = NULL,
           main = NULL,
           space = 0.3,
           breaks_x = waiver(),
           breaks_y = waiver(),
           ylab = NULL,
           G = 100,
           domain = c(0, 1),
           single = FALSE,
           axes = "all_x",
           axis_labels = "all_x",
           xlab = "s",
           legend_names = NULL,
           legend_title = NULL,
           monochrome = FALSE,
           ncol = 3,
           values_discrete = c(0, 1),
           ...) {
    if (pdf && is.null(name)) {
      stop("File name has to be given if pdf is TRUE.")
    }
    interval_width <- (domain[2] - domain[1]) / G
    t <- seq(
      from = domain[1] + 1 / 2 * interval_width,
      to = domain[2] - 1 / 2 * interval_width,
      by = interval_width
    )
    n <- max(ncol(dens_matrix), 1)
    t <- sort(c(t, values_discrete))
    melted <- reshape2::melt(dens_matrix)
    if (ncol(melted) == 1) {
      melted$variable <- 1
    }
    if (ncol(melted) > 2) {
      melted <- reshape2::melt(dens_matrix)[, c(2, 3)]
      colnames(melted)[1] <- "variable"
    }
    melted$t <- rep(t, n)
    discrete_ind <- which(values_discrete == melted$t)
    melted$variable <- as.factor(melted$variable)
    if (!is.null(legend_names)) {
      levels(melted$variable) <- legend_names
    }
    if (monochrome) {
      mcol <- color
      if (is.null(mcol)) {
        mcol <- "black"
      }
      color <- scale_color_manual(values = (rep(mcol, n)))
    }
    if (is.null(xlim)) {
      xlim <- domain
    }
    if (is.null(ylim)) {
      ylim <- c(min(dens_matrix, na.rm = TRUE),
                max(dens_matrix, na.rm = TRUE))
    }
    if (single) {
      ##so machen
      plt <- ggplot() +
        geom_line(data = melted[-discrete_ind, ], aes(
          x = t,
          y = value,
          group = variable,
          color = variable
        )) + geom_point(data = melted[discrete_ind, ], aes(
          x = t,
          y = value,
          group = variable,
          color = variable
        )) +
        color +
        ylab(ylab) + xlab(xlab) +
        guides(color = guide_legend(title = legend_title)) + ggtitle(main) +
        theme(plot.title = element_text(size = 14, hjust = 0.5), ...) +
        scale_x_continuous(breaks = breaks_x, limits = xlim) +
        scale_y_continuous(breaks = breaks_y, limits = ylim)
    } else{
      plt <-
        ggplot() + geom_point(data = melted[discrete_ind, ], aes(
          x = t,
          y = value,
          group = variable,
          color = variable
        )) +
        geom_line(data = melted[-discrete_ind, ], aes(
          x = t,
          y = value,
          group = variable,
          color = variable
        )) +
        theme(legend.position = "none") +
        theme(
          panel.spacing = unit(space, "lines"),
          strip.text.x = element_text(size = 8),
          plot.title = element_text(size = 14, hjust = 0.5)
        ) +
        color +
        facet_wrap( ~ variable,
                    axes = axes,
                    axis.labels = axis_labels,
                    ncol = ncol) +
        ylab(ylab) + xlab(xlab) + ggtitle(main) +
        scale_x_continuous(breaks = breaks_x, limits = xlim) +
        scale_y_continuous(breaks = breaks_y, limits = ylim)
    }
    if (pdf) {
      ggsave(
        name,
        plot = plt,
        path = path,
        create.dir = TRUE,
        device = "pdf",
        height = height ,
        width = width,
        units = unit
      )
    }
    return(plt)
  }






# PLOT FOR DISCRETE MODELS
# plot_densities_discrete is used to plot response densities and predictions (Bayes- and
# clr-level) for our model.
# Arguments:
# dens_matrix: matrix with each column containing the density values of one covariable combination
# pdf: boolean value indicating if the plot shall be saved (as pdf)
# name, width, height, path: indicating the name (with the respective path),
#   width and height of the file, if pdf = TRUE;
#   if NULL, name is created automatically out of the name of the matrix given to
#   dens_matrix
# color, ylim, lty, main, ylab: specify the color gradient (each year corresponds
#   to one color), ylim, the line type, the main, and the y-axis limits of the plot
# G: position of the grid points, where the densities are evaluated
# domain_continuous
# values_discrete
# breaks_x, breaks_y
# xlab, ylab: names for x and y axis
# single: boolean if all curves should be plotted in one plot
# axes, axis_labels: if single=FALSE, should every plot have its own x-axis and label
# breaks_x, breaks_y: breaks for the axes
# legend_names: names for the different levels in the legend
# legend_title: legend title
# monochrome: boolean if all curves should be plotted in the same color if single=TRUE
# ncol: number of columns if single=FALSE for the plot grid
# ...: Further arguments passed to matplot/plot (e.g. las, etc.)

plot_densities_discrete <-
  function(dens_matrix,
           pdf = FALSE,
           name = NULL,
           width = NA,
           height = NA,
           unit = c("in", "cm", "mm", "px"),
           path = "./Images/",
           color = NULL,
           ylim = NULL,
           xlim = NULL,
           main = NULL,
           space = 0.3,
           breaks_x = waiver(),
           breaks_y = waiver(),
           ylab = NULL,
           G = 100,
           domain = c(0, 1),
           single = FALSE,
           axes = "all_x",
           axis_labels = "all_x",
           xlab = "s",
           legend_names = NULL,
           legend_title = NULL,
           monochrome = FALSE,
           ncol = 3,
           values_discrete = c(0, 1),
           ...) {
    if (pdf && is.null(name)) {
      stop("File name has to be given if pdf is TRUE.")
    }
    t <- values_discrete
    n <- max(ncol(dens_matrix), 1)
    melted <- reshape2::melt(dens_matrix)
    if (ncol(melted) == 1) {
      melted$variable <- 1
    }
    if (ncol(melted) > 2) {
      melted <- reshape2::melt(dens_matrix)[, c(2, 3)]
      colnames(melted)[1] <- "variable"
    }
    melted$t <- rep(t, n)
    melted$variable <- as.factor(melted$variable)
    if (!is.null(legend_names)) {
      levels(melted$variable) <- legend_names
    }
    if (monochrome) {
      mcol <- color
      if (is.null(mcol)) {
        mcol <- "black"
      }
      color <- scale_color_manual(values = (rep(mcol, n)))
    }
    if (is.null(xlim)) {
      xlim <- c(min(values_discrete), max(values_discrete))
    }
    if (is.null(ylim)) {
      ylim <- c(min(dens_matrix, na.rm = TRUE),
                max(dens_matrix, na.rm = TRUE))
    }
    if (single) {
      plt <-
        ggplot(data = melted, aes(
          x = t,
          y = value,
          group = variable,
          color = variable
        )) +
        geom_point() +
        color +
        ylab(ylab) + xlab(xlab) +
        guides(color = guide_legend(title = legend_title)) + ggtitle(main) +
        theme(plot.title = element_text(size = 14, hjust = 0.5), ...) +
        scale_x_continuous(breaks = breaks_x, limits = xlim) +
        scale_y_continuous(breaks = breaks_y, limits = ylim)

    } else{
      plt <-
        ggplot(data = melted, aes(
          x = t,
          y = value,
          group = variable,
          color = variable
        )) +
        geom_point() +
        theme(legend.position = "none") +
        theme(
          panel.spacing = unit(space, "lines"),
          strip.text.x = element_text(size = 8),
          plot.title = element_text(size = 14, hjust = 0.5)
        ) +
        color +
        facet_wrap( ~ variable,
                    axes = axes,
                    axis.labels = axis_labels,
                    ncol = ncol) +
        ylab(ylab) + xlab(xlab) + ggtitle(main) +
        scale_x_continuous(breaks = breaks_x, limits = xlim) +
        scale_y_continuous(breaks = breaks_y, limits = ylim)
    }
    if (pdf) {
      ggsave(
        name,
        plot = plt,
        path = path,
        create.dir = TRUE,
        device = "pdf",
        height = height ,
        width = width,
        units = unit
      )
    }
    return(plt)

  }


# PLOT HISTO AND DENSITY FOR CONTINUOUS MODELS
# Plot underlying observed histogram and estimated densities for one covariate combination
# Arguments:
# dta: histogram data
# dens: estimated density
# main: Main title for plot
# automatic_main: if TRUE, an automatic main containing the group_id and covariate values is generated
# G: number of bins
# domain: domain continuous

plot_histo_and_dens <- function(dta,
                                dens,
                                main = NULL,
                                automatic_main = TRUE,
                                G,
                                domain,
                                ...)
{
  interval_width <- (domain[2] - domain[1]) / G
  t <- seq(
    from = domain[1] + 1 / 2 * interval_width,
    to = domain[2] - 1 / 2 * interval_width,
    by = interval_width
  )
  t <- t - domain[1]
  hist <-
    plot_single_histo(dta,
                      main = main,
                      automatic_main = automatic_main,
                      case = "continuous")
  graphics::lines(t, dens, col = "red")
}


# PLOT HISTO AND DENSITY FOR MIXED MODELS
# Plot underlying observed histogram and estimated densities for one covariate combination
# Arguments:
# dta: histogram data
# dens: estimated density
# main: Main title for plot
# automatic_main: if TRUE, an automatic main containing the group_id and covariate values is generated
# G: number of bins
# domain: domain continuous
# discretes: values discrete
plot_histo_and_dens_mixed <- function(dta,
                                      dens,
                                      main = NULL,
                                      automatic_main = TRUE,
                                      G,
                                      domain,
                                      discretes = c(0, 1),
                                      ...)
{
  if ("weighted_counts" %in% colnames(dta)) {
    sum_weights <- sum(dta$weighted_counts)
    dta$freq <- dta$weighted_counts / sum_weights
    dta$height <- dta$freq / dta$Delta
    interval_width <- (domain[2] - domain[1]) / G
    t <- seq(
      from = domain[1] + 1 / 2 * interval_width,
      to = domain[2] - 1 / 2 * interval_width,
      by = interval_width
    )
    t <- sort(c(t, discretes))
    discrete_ind <- match(discretes, t)
    t <- t - domain[1]
    hist <-
      plot_single_histo(
        dta,
        main = main,
        automatic_main = automatic_main,
        ylim = c(0, max(dta$height, dens)),
        xaxt = "n"
      )
    graphics::lines(t[-discrete_ind], dens[-discrete_ind], col = "blue")
    graphics::points(
      t[discrete_ind],
      y = dens[discrete_ind],
      type = "p",
      pch = 20,
      col = "blue"
    )
  } else{
    sum_weights <- sum(dta$counts)
    dta$freq <- dta$counts / sum_weights
    dta$height <- dta$freq / dta$Delta
    interval_width <- (domain[2] - domain[1]) / G
    t <- seq(
      from = domain[1] + 1 / 2 * interval_width,
      to = domain[2] - 1 / 2 * interval_width,
      by = interval_width
    )
    t <- sort(c(t, discretes))
    discrete_ind <- match(discretes, t)
    t <- t - domain[1]
    hist <-
      plot_single_histo(
        dta,
        main = main,
        automatic_main = automatic_main,
        ylim = c(0, max(dta$height, dens)),
        xaxt = "n"
      )
    graphics::lines(t[-discrete_ind], dens[-discrete_ind], col = "blue")
    graphics::points(
      t[discrete_ind],
      y = dens[discrete_ind],
      type = "p",
      pch = 20,
      col = "blue"
    )
  }
}


# PLOT HISTO AND DENSITY FOR DISCRETE MODELS
# Plot underlying observed histogram and estimated densities for one covariate combination
# Arguments:
# dta: histogram data
# dens: estimated density
# main: Main title for plot
# automatic_main: if TRUE, an automatic main containing the group_id and covariate values is generated
# discretes: values discrete
plot_histo_and_dens_discrete <- function(dta,
                                         dens,
                                         main = NULL,
                                         automatic_main = TRUE,
                                         discretes = c(0, 1),
                                         ...)
{
  if ("weighted_counts" %in% colnames(dta)) {
    sum_weights <- sum(dta$weighted_counts)
    dta$freq <- dta$weighted_counts / sum_weights
    dta$height <- dta$freq / dta$Delta
    t <- sort(discretes)
    #t<-t-domain[1]
    hist <-
      plot_single_histo(
        dta,
        main = main,
        automatic_main = automatic_main,
        ylim = c(0, max(dta$height, dens)),
        xaxt = "n",
        case = "discrete"
      )
    graphics::points(
      t,
      y = dens,
      type = "p",
      pch = 20,
      col = "blue"
    )
  } else{
    sum_weights <- sum(dta$counts)
    dta$freq <- dta$counts / sum_weights
    dta$height <- dta$freq / dta$Delta
    t <- sort(discretes)
    #t<-t-domain[1]
    hist <-
      plot_single_histo(
        dta,
        main = main,
        automatic_main = automatic_main,
        ylim = c(0, max(dta$height, dens)),
        xaxt = "n",
        case = "discrete"
      )
    graphics::points(
      t,
      y = dens,
      type = "p",
      pch = 20,
      col = "blue"
    )

  }

}

## INTERACTIVE PLOT FOR CONTINUOUS MODELS
# interactive version of plot_histo_and_dens for all contained groups in dta
interactive_histo_and_dens <- function(dta, dens, G, domain, ...)
{
  manipulate(
    plot_histo_and_dens(
      dta %>% dplyr::filter(group_id == k),
      dens[, k],
      automatic_main = FALSE,
      main = paste("Group-ID", k),
      G,
      domain
    ),
    k = slider(
      min = min(dta$group_id),
      max = max(dta$group_id)
    )
  )

}

## INTERACTIVE PLOT FOR MIXED MODELS
# interactive version of plot_histo_and_dens_mixed for all contained groups in dta
interactive_histo_and_dens_mixed <-
  function(dta, dens, G, domain, discretes, ...)
  {
    manipulate(
      plot_histo_and_dens_mixed(
        dta %>% dplyr::filter(group_id == k),
        dens[, k],
        automatic_main = FALSE,
        main = paste("Group-ID", k),
        G,
        domain,
        discretes
      ),
      k = slider(
        min = min(dta$group_id),
        max = max(dta$group_id)
      )
    )

  }

## INTERACTIVE PLOT FOR DISCRETE MODELS
# interactive version of plot_histo_and_dens_discrete for all contained groups in dta
interactive_histo_and_dens_discrete <-
  function(dta, dens, discretes, ...)
  {
    manipulate(
      plot_histo_and_dens_discrete(
        dta %>% dplyr::filter(group_id == k),
        dens[, k],
        automatic_main = FALSE,
        main = paste("Group-ID", k),
        discretes = discretes
      ),
      k = slider(
        min = min(dta$group_id),
        max = max(dta$group_id)
      )
    )

  }
