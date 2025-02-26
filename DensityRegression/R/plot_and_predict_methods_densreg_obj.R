#' Plot function for conditional density regression models
#'
#' \code{plot.densreg_obj} is the default plot method for data of the class \code{densreg_obj}.
#' @encoding UTF-8
#' @param x \code{densreg_obj}-object, i.e. the output of the \code{densreg}-function.
#' @param type "histo" or "effects": If \code{type = "histo"}, the underlying
#' histogram and the estimated conditional density is plotted for each unique
#' covariate combination. If \code{type = "effects"} the different partial effects
#' are plotted.
#' @param interactive Indicates for \code{type = "effects"} whether the plots are
#' displayed as an interactive plot (\code{TRUE}) or if all plots are plotted consecutively.
#' @param level Indicates for \code{type = "effects"} with "clr" that the clr-transformed
#' effects in \eqn{L^2_0(\mu)} and with "pdf" that the partial effects in
#' \eqn{B^2(\mu)} are plotted.
#' @param display_all Indicates for \code{type = "effects"} and \code{predict = NULL}
#' whether all curves are plotted for effects depending on continuous covariates,
#' i.e., for each observed value of the covariate one curve is plotted,
#' (\code{display_all = TRUE}) or if only three curves (min, med, max) are
#' shown (\code{display_all = FALSE}).
#' If the effects are predicted for new data, i.e., \code{predict} is not
#' \code{NULL}, \code{display_all = FALSE} causes that only the predicted effects
#' for the first, middle and last row are displayed.
#' @param pick_sites Vector of integers from 1 to 6. Indicates for \code{type = "effects"}
#' if only certain effects should be plotted. Number corresponds to the order of
#' the different effect types in the argument of \code{densreg}. If missing
#' (\code{NULL}) all included effects are plotted.
#' @param predict If plots for certain covariate values are desired, \code{predict}
#' has to be a data frame with new data in form of a data table with columns named
#' as the relevant covariates. In each row, the user can specify a value of the
#' respective covariate.
#' @param terms If predict is not \code{NULL}. Vector of term names or the indices
#' of the terms (starting with 1 = intercept) specifying which terms should be
#' predicted and plotted.
#' @param ... Further arguments passed to \code{\link[graphics]{matplot}}/\code{\link[graphics]{plot.default}}
#'
#' @return Plot(s) as specified.
#' @examples
#'
#' \donttest{# please run the examples of densreg to estimate the needed models
#'
#' example("densreg")
#'
#' #' # create new data for predict
#'
#' nd <- data.frame(covariate1 = c("a", "b", "c", "a"), covariate4 = c(0.4, 0.5, 0.1, 0.3),
#'                  covariate2 = c("d", "d", "c", "d"), covariate3 = c(1, 0, 0.2, 2),
#'                  covariate5 = c(0.2, 0.4, 1, 2))
#'
#'
#' # plot mixed model (default settings: histogram, not interactive, density-level, display all)
#'
#' plot(m_mixed)
#'
#' # plot partial effects on clr-level of the continuous model in an interactive plot,
#' # do not show all groups
#'
#' plot(m_cont, type = "effects", interactive = TRUE, level = "clr", display_all = FALSE)
#'
#' # show only first plot
#'
#' plot(m_cont, type = "effects", interactive = FALSE, pick_sites = 1, level = "clr",
#'      display_all = FALSE)
#'
#' # plot partial effects on density-level estimated for new data based on the mixed model
#'
#' plot(m_mixed, type = "effects", level = "pdf", display_all = TRUE, predict = nd)
#'
#' # estimate and plot only the intercept (first term)
#'
#' plot(m_mixed, type = "effects", level = "pdf", display_all = TRUE, predict = nd, terms = 1)
#'
#' #' # estimate and plot only second term
#'
#' plot(m_mixed, type = "effects", level = "pdf", display_all = TRUE, predict = nd, terms = 2)
#'
#' # customize some plot parameters, other plot parameters can also be specified by the user
#'
#' plot(m_mixed, type = "effects", display_all = FALSE, level = "clr", main = "Your title",
#'      legend.position = "none", xlab = "your xlab")
#' }
#'
#' @noRd

plot.densreg_obj <- function(x, type = "histo", interactive = FALSE,
                              level = "pdf", display_all = TRUE, pick_sites = FALSE,
                              predict = NULL, terms = NULL, ...) {
  obj <- x
  if (type == "histo") {
    if (any(obj$count_data$discrete==  TRUE) &
        any(obj$count_data$discrete==  FALSE)) {
      interactive_histo_and_dens_mixed(
        obj$count_data,
        obj$f_hat,
        domain = obj$params[[1]],
        G = obj$params[[3]],
        discretes = obj$params[[2]]
      )
    }
    if (all(obj$count_data$discrete==  FALSE)) {
      interactive_histo_and_dens(
        obj$count_data,
        obj$f_hat,
        domain = obj$params[[1]],
        G = obj$params[[3]],
        ...
      )
    }
    if (all(obj$count_data$discrete==  TRUE)) {
      interactive_histo_and_dens_discrete(obj$count_data,
                                          obj$f_hat, discretes = obj$params[[2]], ...)
    }
  }
  if (type == "effects") {
    if (!is.null(predict)) {
      plot_list <- list()
      G <- obj$params[[3]]
      domain <- obj$params[[1]]
      if (level == "pdf") {
        effects <-
          stats::predict(
            obj,
            type = "terms",
            which = terms,
            new_data = predict,
            level = level
          )
        j <- 1
        for (eff in effects) {

          sp <- "[(), :]"

          parts <- unlist(strsplit(names(effects)[j], sp))
          included <- sapply(colnames(predict), function(substr) {
            any(grepl(substr, parts, ignore.case = TRUE))
          })
          parts <- colnames(predict)[included]



          #parts <- parts[parts%in%(colnames(predict))]
          terms_title <- paste(parts, collapse = ", ")
          legend_labels <- paste(as.data.frame(unlist(t((predict%>%select(parts))))))
          #legend_labels <- gsub('["c(")]', '', legend_labels)
          if ( names(effects)[j] == "intercept") {
            eff <- eff[, 1]
            single <- TRUE
            terms_title <- "intercept"
            legend.position = "none"
            legend_labels <- "intercept"
          } else {
            if (display_all == FALSE) {
              if (length(legend_labels) == 1) {
                indx <- c(1)}
              if (length(legend_labels) == 2) {
                indx <- c(1, length(legend_labels))}
              if (length(legend_labels) > 2) {
                indx <- c(1, round(length(legend_labels)/2), length(legend_labels))}

            } else {
              indx <- 1:length(legend_labels)
            }
            single <- TRUE
            legend.position = NULL
            eff <- eff[, indx]
            legend_labels[indx]
          }



          if (all(obj$count_data$discrete==  FALSE)) {

            p <-
              plot_densities(
                eff,
                G = G,
                domain = domain,
                #single = TRUE,
                main = names(effects)[j],
                legend_title = terms_title,
                ylab = "density",
                single = single,
                legend.position = legend.position,
                legend_names = legend_labels
              )
            plot_list <- append(plot_list, list(p))
            j <- j + 1}

          if (all(obj$count_data$discrete==  TRUE)) {
            p <-
              plot_densities_discrete(
                eff,
                G = G,
                values_discrete = unlist(obj$params[2]),
                # single = TRUE,
                main = names(effects)[j],
                legend_title = terms_title,
                ylab = "density",
                single = single,
                legend.position = legend.position,
                legend_names = legend_labels
              )
            plot_list <- append(plot_list, list(p))
            j <- j + 1
          }
          if (any(obj$count_data$discrete==  FALSE) &
              any(obj$count_data$discrete==  TRUE)) {
            p <-
              plot_densities_mixed(
                eff,
                G = G,
                domain = domain,
                values_discrete = unlist(obj$params[2]),
                # single = TRUE,
                main = names(effects)[j],
                legend_title = terms_title,
                ylab = "density",
                single = single,
                legend.position = legend.position,
                legend_names = legend_labels
              )
            plot_list <- append(plot_list, list(p))
            j <- j + 1
          }

        }
      }
      else{
        effects <-
          stats::predict(
            obj,
            type = "terms",
            which = terms,
            new_data = predict
          )
        j <- 1
        for (eff in effects) {
          if (all(obj$count_data$discrete==  FALSE)) {
            lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
            p <-
              plot_densities(
                eff,
                G = G,
                domain = domain,
                single = TRUE,
                main = names(effects)[j],
                legend_title = terms_title,
                legend_names = "",
                ylab = "clr(density)"
              )
            plot_list <- append(plot_list, list(p))
            j <- j + 1
          }
          if (all(obj$count_data$discrete==  TRUE)) {
            p <-
              plot_densities_discrete(
                eff,
                G = G,
                values_discrete = unlist(obj$params[2]),
                single = TRUE,
                main = names(effects)[j],
                legend_title = terms_title,
                ylab = "clr(density)"
              )
            plot_list <- append(plot_list, list(p))
            j <- j + 1
          }
          if (any(obj$count_data$discrete==  FALSE) &
              any(obj$count_data$discrete==  TRUE)) {
            p <-
              plot_densities_mixed(
                eff,
                G = G,
                domain = domain,
                values_discrete = unlist(obj$params[2]),
                single = TRUE,
                main = names(effects)[j],
                legend_title = terms_title,
                ylab = "clr(density)"
              )
            plot_list <- append(plot_list, list(p))
            j <- j + 1
          }

        }
      }

      if (!interactive) {
        return(plot_list)
      } else{
        (manipulate(plot_list[[k]], k = slider(
          min = 1, max = length(plot_list), step = 1
        )))
      }


    } else{
      if (all(obj$count_data$discrete==  FALSE)) {
        if (level==  "pdf") {
          if (display_all) {
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities(
                obj$effects[["pdf_estimated_effects"]][["intercept"]],
                G = obj$params[[3]],
                legend_names = "intercept",
                single = TRUE,
                ylab = expression(hat(beta)[0]),
                legend.position = "none",
                domain = obj$params[[1]],
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities(
                  cbind(
                    data.frame(ref = rep(1 / (
                      domain[2] - domain[1]
                    ), G)),
                    data.frame(obj$effects[["pdf_estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0("Group-specific intercept for ", group_spec),
                  legend_title = paste0(group_spec),
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                num_lev <- length(lev)
                p <-
                  plot_densities(
                    obj$effects[["pdf_estimated_effects"]][[j]],
                    G = G,
                    single = TRUE,
                    legend_names = round(lev, digits = 3),
                    ylab = paste0("Smooth effect of ", smooth_effects[[1]]),
                    legend_title = paste0(smooth_effects[[1]]),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev_by <- unique(obj$count_data[[smooth_effects[[5]]]])

                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  lev <-
                    obj$count_data %>% select(smooth_effects[[5]], smooth_effects[[1]])
                  lev <- lev %>% filter(.[[1]]==  lev_by[i + 1])
                  lev <- sort(unname(unlist(unique(
                    lev[, 2]
                  ))))
                  num_lev <- length(lev)
                  p <-
                    plot_densities(
                      obj$effects[["pdf_estimated_effects"]][[j]],
                      G = G,
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0(
                        "Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1]
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              p <-
                plot_densities(
                  obj$effects[["pdf_estimated_effects"]][[j]],
                  G = G,
                  single = TRUE,
                  legend_names = round(lev, digits = 3),
                  ylab = paste0("Linear effect of ", linear),
                  legend_title = paste0(linear)
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities(
                  obj$effects[["pdf_estimated_effects"]][[j]],
                  G = G,
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]]
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }
            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities(
                  obj$effects[["pdf_estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", ")
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          } else{
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities(
                obj$effects[["pdf_estimated_effects"]][["intercept"]],
                G = obj$params[[3]],
                legend_names = "intercept",
                single = TRUE,
                ylab = expression(hat(beta)[0]),
                legend.position = "none",
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities(
                  cbind(
                    data.frame(ref = rep(1 / (
                      domain[2] - domain[1]
                    ), G)),
                    data.frame(obj$effects[["pdf_estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0("Group-specific intercept for ", group_spec),
                  legend_title = paste0(group_spec),
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                lev_used <- c(min(lev), stats::median(lev), max(lev))
                used_cols <- c(1, stats::median(c(1:length(
                  lev
                ))), length(lev))
                p <-
                  plot_densities(
                    obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                    G = G,
                    single = TRUE,
                    legend_names = paste0(
                      c("min ", "med ", "max "),
                      round(lev_used, digits = 3)
                    ),
                    ylab = paste0("Smooth effect of ", smooth_effects[[1]]),
                    legend_title = paste0(smooth_effects[[1]]),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev_by <- sort(unique(obj$count_data[[smooth_effects[[5]]]]))

                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  lev <-
                    obj$count_data %>% select(smooth_effects[[5]], smooth_effects[[1]])
                  lev <- lev %>% filter(.[[1]]==  lev_by[i + 1])
                  lev <- sort(unname(unlist(unique(
                    lev[, 2]
                  ))))
                  num_lev <- length(lev)
                  lev_used <- c(min(lev), stats::median(lev), max(lev))
                  used_cols <- c(1, stats::median(c(1:length(
                    lev
                  ))), length(lev))

                  p <-
                    plot_densities(
                      obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                      G = G,
                      single = TRUE,
                      legend_names = paste0(
                        c("min ", "med ", "max "),
                        round(lev_used, digits = 3)
                      ),
                      ylab = paste0(
                        "Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1]
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              lev_used <- c(min(lev), stats::median(lev), max(lev))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              p <-
                plot_densities(
                  obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                  G = G,
                  single = TRUE,
                  legend_names = paste0(
                    c("min ", "med ", "max "),
                    round(lev_used, digits = 3)
                  ),
                  ylab = paste0("Linear effect of ", linear),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              lev_used <- lev[used_cols]

              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities(
                  obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                  G = G,
                  single = TRUE,
                  legend_names =  lev_used,
                  ylab = paste0(
                    "varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]]
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }
            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities(
                  obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", ")
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )

              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          }


        }
        if (level==  "clr") {
          if (display_all) {
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities(
                obj$effects[["estimated_effects"]][["intercept"]],
                G = obj$params[[3]],
                legend_names = "intercept",
                single = TRUE,
                ylab = expression(paste("clr(", hat(beta)[0], ")")),
                legend.position = "none",
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities(
                  cbind(
                    data.frame(ref = rep(0, G)),
                    data.frame(obj$effects[["estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "clr(Group-specific intercept for ",
                    group_spec,
                    ")"
                  ),
                  legend_title = paste0(group_spec),
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                num_lev <- length(lev)
                p <-
                  plot_densities(
                    obj$effects[["estimated_effects"]][[j]],
                    G = G,
                    single = TRUE,
                    legend_names = round(lev, digits = 3),
                    ylab = paste0("clr(Smooth effect of ", smooth_effects[[1]]),
                    legend_title = paste0(smooth_effects[[1]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                lev_by <-
                  unique(obj$count_data[[smooth_effects[[5]]]])
                num_lev <- length(lev)
                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  p <-
                    plot_densities(
                      obj$effects[["estimated_effects"]][[j]],
                      G = G,
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0(
                        "clr(Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1],
                        ")"
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              p <-
                plot_densities(
                  obj$effects[["estimated_effects"]][[j]],
                  G = G,
                  single = TRUE,
                  legend_names = round(lev, digits = 3),
                  ylab = paste0("clr(Linear effect of ", linear, ")"),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities(
                  obj$effects[["estimated_effects"]][[j]],
                  G = G,
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "clr(varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]],
                    ")"
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities(
                  obj$effects[["estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "clr(Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", "),
                    ")"
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          }
          else{
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities(
                obj$effects[["estimated_effects"]][["intercept"]],
                G = obj$params[[3]],
                legend_names = "clr(intercept)",
                single = TRUE,
                ylab = expression(hat(beta)[0]),
                legend.position = "none",
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities(
                  cbind(
                    data.frame(ref = rep(0, G)),
                    data.frame(obj$effects[["estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "clr(Group-specific intercept for ",
                    group_spec,
                    ")"
                  ),
                  legend_title = paste0(group_spec),
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                lev_used <- c(min(lev), stats::median(lev), max(lev))
                used_cols <- c(1, stats::median(c(1:length(
                  lev
                ))), length(lev))
                p <-
                  plot_densities(
                    obj$effects[["estimated_effects"]][[j]][, used_cols],
                    G = G,
                    single = TRUE,
                    legend_names = paste0(
                      c("min ", "med ", "max "),
                      round(lev_used, digits = 3)
                    ),
                    ylab = paste0(
                      "clr(Smooth effect of ",
                      smooth_effects[[1]],
                      ")"
                    ),
                    legend_title = paste0(smooth_effects[[1]]),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev_by <- sort(unique(obj$count_data[[smooth_effects[[5]]]]))

                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  lev <-
                    obj$count_data %>% select(smooth_effects[[5]], smooth_effects[[1]])
                  lev <- lev %>% filter(.[[1]]==  lev_by[i + 1])
                  lev <- sort(unname(unlist(unique(
                    lev[, 2]
                  ))))
                  num_lev <- length(lev)
                  lev_used <- c(min(lev), stats::median(lev), max(lev))
                  used_cols <- c(1, stats::median(c(1:length(
                    lev
                  ))), length(lev))

                  p <-
                    plot_densities(
                      obj$effects[["estimated_effects"]][[j]][, used_cols],
                      G = G,
                      single = TRUE,
                      legend_names = paste0(
                        c("min ", "med ", "max "),
                        round(lev_used, digits = 3)
                      ),
                      ylab = paste0(
                        "clr(Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1],
                        ")"
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              lev_used <- c(min(lev), stats::median(lev), max(lev))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              p <-
                plot_densities(
                  obj$effects[["estimated_effects"]][[j]][, used_cols],
                  G = G,
                  single = TRUE,
                  legend_names = paste0(
                    c("min ", "med ", "max "),
                    round(lev_used, digits = 3)
                  ),
                  ylab = paste0("clr(Linear effect of ", linear, ")"),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              lev_used <- lev[used_cols]

              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities(
                  obj$effects[["estimated_effects"]][[j]][, used_cols],
                  G = G,
                  single = TRUE,
                  legend_names =  lev_used,
                  ylab = paste0(
                    "clr(varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]],
                    ")"
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities(
                  obj$effects[["estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "clr(Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", "),
                    ")"
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )

              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          }



        }
      }

      if (all(obj$count_data$discrete==  TRUE)) {
        if (level==  "pdf") {
          if (display_all) {
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities_discrete(
                obj$effects[["pdf_estimated_effects"]][["intercept"]],
                G = obj$params[[3]],
                legend_names = "intercept",
                single = TRUE,
                ylab = expression(hat(beta)[0]),
                legend.position = "none",
                values_discrete = obj$params[[2]],
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities_discrete(
                  cbind(
                    data.frame(ref = rep(
                      1 / length(obj$params[[2]]), length(obj$params[[2]])
                    )),
                    data.frame(obj$effects[["pdf_estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0("Group-specific intercept for ", group_spec),
                  legend_title = paste0(group_spec),
                  values_discrete = obj$params[[2]],
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                num_lev <- length(lev)
                p <-
                  plot_densities_discrete(
                    obj$effects[["pdf_estimated_effects"]][[j]],
                    G = G,
                    single = TRUE,
                    legend_names = round(lev, digits = 3),
                    ylab = paste0("Smooth effect of ", smooth_effects[[1]]),
                    legend_title = paste0(smooth_effects[[1]]),
                    values_discrete = obj$params[[2]],
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev_by <- unique(obj$count_data[[smooth_effects[[5]]]])

                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  lev <-
                    obj$count_data %>% select(smooth_effects[[5]], smooth_effects[[1]])
                  lev <- lev %>% filter(.[[1]]==  lev_by[i + 1])
                  lev <- sort(unname(unlist(unique(
                    lev[, 2]
                  ))))
                  num_lev <- length(lev)
                  p <-
                    plot_densities_discrete(
                      obj$effects[["pdf_estimated_effects"]][[j]],
                      G = G,
                      values_discrete = obj$params[[2]],
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0(
                        "Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1]
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              p <-
                plot_densities_discrete(
                  obj$effects[["pdf_estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = round(lev, digits = 3),
                  ylab = paste0("Linear effect of ", linear),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities_discrete(
                  obj$effects[["pdf_estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]]
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }
            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities_discrete(
                  obj$effects[["pdf_estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", ")
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )

              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          } else{
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities_discrete(
                obj$effects[["pdf_estimated_effects"]][["intercept"]],
                values_discrete = obj$params[[2]],
                G = obj$params[[3]],
                legend_names = "intercept",
                single = TRUE,
                ylab = expression(hat(beta)[0]),
                legend.position = "none",
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities_discrete(
                  cbind(
                    data.frame(ref = rep(
                      1 / length(obj$params[[2]]), length(obj$params[[2]])
                    )),
                    data.frame(obj$effects[["pdf_estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0("Group-specific intercept for ", group_spec),
                  legend_title = paste0(group_spec),
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                lev_used <- c(min(lev), stats::median(lev), max(lev))
                used_cols <- c(1, stats::median(c(1:length(
                  lev
                ))), length(lev))
                p <-
                  plot_densities_discrete(
                    obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                    G = G,
                    values_discrete = obj$params[[2]],
                    single = TRUE,
                    legend_names = paste0(
                      c("min ", "med ", "max "),
                      round(lev_used, digits = 3)
                    ),
                    ylab = paste0("Smooth effect of ", smooth_effects[[1]]),
                    legend_title = paste0(smooth_effects[[1]]),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev_by <- sort(unique(obj$count_data[[smooth_effects[[5]]]]))

                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  lev <-
                    obj$count_data %>% select(smooth_effects[[5]], smooth_effects[[1]])
                  lev <- lev %>% filter(.[[1]]==  lev_by[i + 1])
                  lev <- sort(unname(unlist(unique(
                    lev[, 2]
                  ))))
                  num_lev <- length(lev)
                  lev_used <- c(min(lev), stats::median(lev), max(lev))
                  used_cols <- c(1, stats::median(c(1:length(
                    lev
                  ))), length(lev))

                  p <-
                    plot_densities_discrete(
                      obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                      G = G,
                      values_discrete = obj$params[[2]],
                      single = TRUE,
                      legend_names = paste0(
                        c("min ", "med ", "max "),
                        round(lev_used, digits = 3)
                      ),
                      ylab = paste0(
                        "Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1]
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              lev_used <- c(min(lev), stats::median(lev), max(lev))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              p <-
                plot_densities_discrete(
                  obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = paste0(
                    c("min ", "med ", "max "),
                    round(lev_used, digits = 3)
                  ),
                  ylab = paste0("Linear effect of ", linear),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              lev_used <- lev[used_cols]

              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities_discrete(
                  obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev_used,
                  ylab = paste0(
                    "varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]]
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }
            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities_discrete(
                  obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", ")
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )

              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          }


        }
        if (level==  "clr") {
          if (display_all) {
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities_discrete(
                obj$effects[["estimated_effects"]][["intercept"]],
                values_discrete = obj$params[[2]],
                G = obj$params[[3]],
                legend_names = "intercept",
                single = TRUE,
                ylab = expression(paste("clr(", hat(beta)[0], ")")),
                legend.position = "none",
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities_discrete(
                  cbind(
                    data.frame(ref = rep(0, length(
                      obj$params[[2]]
                    ))),
                    data.frame(obj$effects[["estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "clr(Group-specific intercept for ",
                    group_spec,
                    ")"
                  ),
                  legend_title = paste0(group_spec),
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                num_lev <- length(lev)
                p <-
                  plot_densities_discrete(
                    obj$effects[["estimated_effects"]][[j]],
                    G = G,
                    values_discrete = obj$params[[2]],
                    single = TRUE,
                    legend_names = round(lev, digits = 3),
                    ylab = paste0("clr(Smooth effect of ", smooth_effects[[1]]),
                    legend_title = paste0(smooth_effects[[1]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                lev_by <-
                  unique(obj$count_data[[smooth_effects[[5]]]])
                num_lev <- length(lev)
                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  p <-
                    plot_densities_discrete(
                      obj$effects[["estimated_effects"]][[j]],
                      G = G,
                      values_discrete = obj$params[[2]],
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0(
                        "clr(Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1],
                        ")"
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              p <-
                plot_densities_discrete(
                  obj$effects[["estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = round(lev, digits = 3),
                  ylab = paste0("clr(Linear effect of ", linear, ")"),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities_discrete(
                  obj$effects[["estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "clr(varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]],
                    ")"
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities_discrete(
                  obj$effects[["estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "clr(Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", "),
                    ")"
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )

              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }
          }
          else{
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities_discrete(
                obj$effects[["estimated_effects"]][["intercept"]],
                G = obj$params[[3]],
                values_discrete = obj$params[[2]],
                legend_names = "clr(intercept)",
                single = TRUE,
                ylab = expression(hat(beta)[0]),
                legend.position = "none",
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities_discrete(
                  cbind(
                    data.frame(ref = rep(0, length(
                      obj$params[[2]]
                    ))),
                    data.frame(obj$effects[["estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "clr(Group-specific intercept for ",
                    group_spec,
                    ")"
                  ),
                  legend_title = paste0(group_spec),
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                lev_used <- c(min(lev), stats::median(lev), max(lev))
                used_cols <- c(1, stats::median(c(1:length(
                  lev
                ))), length(lev))
                p <-
                  plot_densities_discrete(
                    obj$effects[["estimated_effects"]][[j]][, used_cols],
                    G = G,
                    values_discrete = obj$params[[2]],
                    single = TRUE,
                    legend_names = paste0(
                      c("min ", "med ", "max "),
                      round(lev_used, digits = 3)
                    ),
                    ylab = paste0(
                      "clr(Smooth effect of ",
                      smooth_effects[[1]],
                      ")"
                    ),
                    legend_title = paste0(smooth_effects[[1]]),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev_by <- sort(unique(obj$count_data[[smooth_effects[[5]]]]))

                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  lev <-
                    obj$count_data %>% select(smooth_effects[[5]], smooth_effects[[1]])
                  lev <- lev %>% filter(.[[1]]==  lev_by[i + 1])
                  lev <- sort(unname(unlist(unique(
                    lev[, 2]
                  ))))
                  num_lev <- length(lev)
                  lev_used <- c(min(lev), stats::median(lev), max(lev))
                  used_cols <- c(1, stats::median(c(1:length(
                    lev
                  ))), length(lev))

                  p <-
                    plot_densities_discrete(
                      obj$effects[["estimated_effects"]][[j]][, used_cols],
                      G = G,
                      values_discrete = obj$params[[2]],
                      single = TRUE,
                      legend_names = paste0(
                        c("min ", "med ", "max "),
                        round(lev_used, digits = 3)
                      ),
                      ylab = paste0(
                        "clr(Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1],
                        ")"
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              lev_used <- c(min(lev), stats::median(lev), max(lev))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              p <-
                plot_densities_discrete(
                  obj$effects[["estimated_effects"]][[j]][, used_cols],
                  G = G,
                  single = TRUE,
                  values_discrete = obj$params[[2]],
                  legend_names = paste0(
                    c("min ", "med ", "max "),
                    round(lev_used, digits = 3)
                  ),
                  ylab = paste0("clr(Linear effect of ", linear, ")"),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              lev_used <- lev[used_cols]

              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities_discrete(
                  obj$effects[["estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev_used,
                  ylab = paste0(
                    "clr(varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]],
                    ")"
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities_mixed(
                  obj$effects[["estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "clr(Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", "),
                    ", "
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )

              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          }



        }
      }
      if (any(obj$count_data$discrete==  TRUE) &&
          any(obj$count_data$discrete==  FALSE)) {
        if (level==  "pdf") {
          if (display_all) {
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities_mixed(
                obj$effects[["pdf_estimated_effects"]][["intercept"]],
                G = obj$params[[3]],
                legend_names = "intercept",
                single = TRUE,
                ylab = expression(hat(beta)[0]),
                legend.position = "none",
                values_discrete = obj$params[[2]],
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities_mixed(
                  cbind(
                    data.frame(ref = rep(
                      1 / length(obj$params[[2]]), length(obj$params[[2]])
                    )),
                    data.frame(obj$effects[["pdf_estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0("Group-specific intercept for ", group_spec),
                  legend_title = paste0(group_spec),
                  values_discrete = obj$params[[2]],
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                num_lev <- length(lev)
                p <-
                  plot_densities_mixed(
                    obj$effects[["pdf_estimated_effects"]][[j]],
                    G = G,
                    single = TRUE,
                    legend_names = round(lev, digits = 3),
                    ylab = paste0("Smooth effect of ", smooth_effects[[1]]),
                    legend_title = paste0(smooth_effects[[1]]),
                    values_discrete = obj$params[[2]],
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev_by <- unique(obj$count_data[[smooth_effects[[5]]]])

                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  lev <-
                    obj$count_data %>% select(smooth_effects[[5]], smooth_effects[[1]])
                  lev <- lev %>% filter(.[[1]]==  lev_by[i + 1])
                  lev <- sort(unname(unlist(unique(
                    lev[, 2]
                  ))))
                  num_lev <- length(lev)
                  p <-
                    plot_densities_mixed(
                      obj$effects[["pdf_estimated_effects"]][[j]],
                      G = G,
                      domain = obj$params[[1]],
                      values_discrete = obj$params[[2]],
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0(
                        "Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1]
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              p <-
                plot_densities_mixed(
                  obj$effects[["pdf_estimated_effects"]][[j]],
                  G = G,
                  domain = obj$params[[1]],
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = round(lev, digits = 3),
                  ylab = paste0("Linear effect of ", linear),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities_mixed(
                  obj$effects[["pdf_estimated_effects"]][[j]],
                  G = G,
                  domain = obj$params[[1]],
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]]
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }
            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities_mixed(
                  obj$effects[["pdf_estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", ")
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          } else{
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities_mixed(
                obj$effects[["pdf_estimated_effects"]][["intercept"]],
                values_discrete = obj$params[[2]],
                G = obj$params[[3]],
                legend_names = "intercept",
                single = TRUE,
                ylab = expression(hat(beta)[0]),
                legend.position = "none",
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities_mixed(
                  cbind(
                    data.frame(ref = rep(
                      1 / length(obj$params[[2]]), length(obj$params[[2]])
                    )),
                    data.frame(obj$effects[["pdf_estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0("Group-specific intercept for ", group_spec),
                  legend_title = paste0(group_spec),
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                lev_used <- c(min(lev), stats::median(lev), max(lev))
                used_cols <- c(1, stats::median(c(1:length(
                  lev
                ))), length(lev))
                p <-
                  plot_densities_mixed(
                    obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                    G = G,
                    values_discrete = obj$params[[2]],
                    single = TRUE,
                    legend_names = paste0(
                      c("min ", "med ", "max "),
                      round(lev_used, digits = 3)
                    ),
                    ylab = paste0("Smooth effect of ", smooth_effects[[1]]),
                    legend_title = paste0(smooth_effects[[1]]),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev_by <- sort(unique(obj$count_data[[smooth_effects[[5]]]]))

                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  lev <-
                    obj$count_data %>% select(smooth_effects[[5]], smooth_effects[[1]])
                  lev <- lev %>% filter(.[[1]]==  lev_by[i + 1])
                  lev <- sort(unname(unlist(unique(
                    lev[, 2]
                  ))))
                  num_lev <- length(lev)
                  lev_used <- c(min(lev), stats::median(lev), max(lev))
                  used_cols <- c(1, stats::median(c(1:length(
                    lev
                  ))), length(lev))

                  p <-
                    plot_densities_mixed(
                      obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                      G = G,
                      values_discrete = obj$params[[2]],
                      single = TRUE,
                      legend_names = paste0(
                        c("min ", "med ", "max "),
                        round(lev_used, digits = 3)
                      ),
                      ylab = paste0(
                        "Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1]
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              lev_used <- c(min(lev), stats::median(lev), max(lev))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              p <-
                plot_densities_mixed(
                  obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = paste0(
                    c("min ", "med ", "max "),
                    round(lev_used, digits = 3)
                  ),
                  ylab = paste0("Linear effect of ", linear),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              lev_used <- lev[used_cols]

              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities_mixed(
                  obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev_used,
                  ylab = paste0(
                    "varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]]
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }
            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities_mixed(
                  obj$effects[["pdf_estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", ")
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )

              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          }


        }
        if (level==  "clr") {
          if (display_all) {
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities_mixed(
                obj$effects[["estimated_effects"]][["intercept"]],
                values_discrete = obj$params[[2]],
                G = obj$params[[3]],
                legend_names = "intercept",
                single = TRUE,
                ylab = expression(paste("clr(", hat(beta)[0], ")")),
                legend.position = "none",
                ...
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities_mixed(
                  cbind(
                    data.frame(ref = rep(0, length(
                      obj$params[[2]]
                    ))),
                    data.frame(obj$effects[["estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "clr(Group-specific intercept for ",
                    group_spec,
                    ")"
                  ),
                  legend_title = paste0(group_spec),
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                num_lev <- length(lev)
                p <-
                  plot_densities_mixed(
                    obj$effects[["estimated_effects"]][[j]],
                    G = G,
                    values_discrete = obj$params[[2]],
                    single = TRUE,
                    legend_names = round(lev, digits = 3),
                    ylab = paste0("clr(Smooth effect of ", smooth_effects[[1]]),
                    legend_title = paste0(smooth_effects[[1]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                lev_by <-
                  unique(obj$count_data[[smooth_effects[[5]]]])
                num_lev <- length(lev)
                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  p <-
                    plot_densities_mixed(
                      obj$effects[["estimated_effects"]][[j]],
                      G = G,
                      values_discrete = obj$params[[2]],
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0(
                        "clr(Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1],
                        ")"
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              p <-
                plot_densities_mixed(
                  obj$effects[["estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = round(lev, digits = 3),
                  ylab = paste0("clr(Linear effect of ", linear, ")"),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities_mixed(
                  obj$effects[["estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "clr(varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]],
                    ")"
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities_mixed(
                  obj$effects[["estimated_effects"]][[j]],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "clr(Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", "),
                    ")"
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          }
          else{
            G <- obj$params[[3]]
            domain <- obj$params[[1]]
            plot_intercept <-
              plot_densities_mixed(
                obj$effects[["estimated_effects"]][["intercept"]],
                G = obj$params[[3]],
                values_discrete = obj$params[[2]],
                legend_names = "clr(intercept)",
                single = TRUE,
                ylab = expression(hat(beta)[0]),
                legend.position = "none"
              )
            plot_list <- list(list(plot_intercept))
            j <- 1
            for (group_spec in obj$predicted_effects[[1]]) {
              lev <- sort(unique(obj$count_data[[group_spec]]))
              num_lev <- length(lev)
              ind_eff <- c((j + 1):(j + num_lev - 1))

              p <-
                plot_densities_mixed(
                  cbind(
                    data.frame(ref = rep(0, length(
                      obj$params[[2]]
                    ))),
                    data.frame(obj$effects[["estimated_effects"]])[, ind_eff]
                  ),
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names = lev,
                  ylab = paste0(
                    "clr(Group-specific intercept for ",
                    group_spec,
                    ")"
                  ),
                  legend_title = paste0(group_spec),
                  ...
                )
              j <- j + num_lev - 1
              plot_list <- append(plot_list, list(p))
            }
            j <- j + 1
            for (smooth_effects in obj$predicted_effects[[2]]) {
              if (is.null(smooth_effects[[5]])) {
                lev <- sort(unique(obj$count_data[[smooth_effects[[1]]]]))
                lev_used <- c(min(lev), stats::median(lev), max(lev))
                used_cols <- c(1, stats::median(c(1:length(
                  lev
                ))), length(lev))
                p <-
                  plot_densities_mixed(
                    obj$effects[["estimated_effects"]][[j]][, used_cols],
                    G = G,
                    values_discrete = obj$params[[2]],
                    single = TRUE,
                    legend_names = paste0(
                      c("min ", "med ", "max "),
                      round(lev_used, digits = 3)
                    ),
                    ylab = paste0(
                      "clr(Smooth effect of ",
                      smooth_effects[[1]],
                      ")"
                    ),
                    legend_title = paste0(smooth_effects[[1]]),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              else{
                lev_by <- sort(unique(obj$count_data[[smooth_effects[[5]]]]))

                num_lev_by <- length(lev_by)
                for (i in c(1:(num_lev_by - 1))) {
                  lev <-
                    obj$count_data %>% select(smooth_effects[[5]], smooth_effects[[1]])
                  lev <- lev %>% filter(.[[1]]==  lev_by[i + 1])
                  lev <- sort(unname(unlist(unique(
                    lev[, 2]
                  ))))
                  num_lev <- length(lev)
                  lev_used <- c(min(lev), stats::median(lev), max(lev))
                  used_cols <- c(1, stats::median(c(1:length(
                    lev
                  ))), length(lev))

                  p <-
                    plot_densities_mixed(
                      obj$effects[["estimated_effects"]][[j]][, used_cols],
                      G = G,
                      values_discrete = obj$params[[2]],
                      single = TRUE,
                      legend_names = paste0(
                        c("min ", "med ", "max "),
                        round(lev_used, digits = 3)
                      ),
                      ylab = paste0(
                        "clr(Smooth effect of ",
                        smooth_effects[[1]],
                        " given ",
                        smooth_effects[[5]],
                        " = ",
                        lev_by[i + 1],
                        ")"
                      ),
                      legend_title = paste0(smooth_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }

              }
            }

            for (linear in obj$predicted_effects[[3]]) {
              lev <- sort(unique(obj$count_data[[linear]]))
              num_lev <- length(lev)
              lev_used <- c(min(lev), stats::median(lev), max(lev))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              p <-
                plot_densities_mixed(
                  obj$effects[["estimated_effects"]][[j]][, used_cols],
                  G = G,
                  single = TRUE,
                  values_discrete = obj$params[[2]],
                  legend_names = paste0(
                    c("min ", "med ", "max "),
                    round(lev_used, digits = 3)
                  ),
                  ylab = paste0("clr(Linear effect of ", linear, ")"),
                  legend_title = paste0(linear),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (varying_coef in obj$predicted_effects[[4]]) {
              lev_base <-
                paste0(
                  "(",
                  round(obj$count_data[[varying_coef[[1]]]], digits = 3),
                  ", ",
                  round(obj$count_data[[varying_coef[[2]]]], digits = 3),
                  ")"
                )
              lev <- sort(unique(lev_base))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))
              lev_by <- unique(obj$count_data[[varying_coef[[2]]]])
              lev_used <- lev[used_cols]

              num_lev <- length(lev)
              num_lev_by <- length(lev_by)
              p <-
                plot_densities_mixed(
                  obj$effects[["estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev_used,
                  ylab = paste0(
                    "clr(varying effect of ",
                    varying_coef[[1]],
                    " given ",
                    varying_coef[[2]],
                    ")"
                  ),
                  legend_title = paste0("(", varying_coef[[1]], ", ", varying_coef[[2]], ")"),
                  ...
                )
              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }

            for (flexi in obj$predicted_effects[[5]]) {
              lev_base <- round(obj$count_data %>% select(flexi[[1]]), digits = 3)

              lev <- unite(lev_base,
                           col = 'com',
                           flexi[[1]],
                           sep = ', ')

              lev <- sort(unique(lev$com))
              used_cols <-
                c(1, stats::median(c(1:length(lev))), length(lev))

              lev_used <- lev[used_cols]

              num_lev <- length(lev)

              p <-
                plot_densities_mixed(
                  obj$effects[["estimated_effects"]][[j]][, used_cols],
                  G = G,
                  values_discrete = obj$params[[2]],
                  single = TRUE,
                  legend_names =  lev,
                  ylab = paste0(
                    "clr(Smooth interaction of ",
                    paste0(unlist(flexi[[1]]), collapse = ", "),
                    ")"
                  ),
                  legend_title = paste0("(", paste0(unlist(
                    flexi[[1]]
                  ), collapse = ", "), ")"),
                  ...
                )

              j <- j + 1
              plot_list <- append(plot_list, list(p))
            }



          }



        }
      }




      if (!isFALSE(pick_sites)) {
        plot_list <- plot_list[pick_sites]
      }

      if (!interactive) {
        return(plot_list)
      } else{
        (manipulate(plot_list[[k]], k = slider(
          min = 1, max = length(plot_list), step = 1
        )))
      }
    }
  }

}





#' Predict function for conditonal density regression models
#'
#' \code{predict.densreg_obj} is the default predict method for data of the class \code{densreg_obj}.

#' @param object \code{densreg_obj}-object, i.e. the output of the \code{densreg}-function.
#' @param new_data New data in form of a data table with columns named as the relevant covariates for the terms which should be predicted. In each row, the user can specify a value of the respective covariate. If not specyfied (\code{NULL}) the values of the original data set (see object$count_data) are used for the predictions.
#' @param which Only terms (or their index number) named in this array will be predicted. Covariates only needed for these terms have to be given in \code{new_data}. If \code{which} and \code{exclude} are given, \code{exclude} will be used.
#' @param exclude Terms (or their index number) named in this array will not be predicted. Covariates only needed for terms which are excluded do not have to be given in \code{new_data}. If \code{which} and \code{exclude} are given, \code{exclude} will be used.
#' @param type "terms", "pdf" or "clr": If \code{type = "terms"} each component of the linear predictor is returned seperately, if
#' @param level If \code{type = "terms"}: "pdf" results in predicted terms on pdf-level, "clr" in terms on clr-level.  \code{type = "pdf"}  or \code{type = "clr"} specifies the prediction of \eqn{\hat f} or \eqn{clr(\hat f)} of the respective covariate values.
#'
#' @return A list of matrices (if type = "terms") with one matrix for each term, different columns for every predicted covariate value.
#'A matrix with columns for the different covariate combinations if \code{type = "pdf"} or \code{ = "clr"} containing the estimated \eqn{\hat f} or \eqn{clr(\hat f)}.
#' @examples
#'
#' \donttest{# please run the examples of densreg to estimate the needed models
#'
#' example("densreg")
#'
#' # see names of the models' terms
#' sapply(m_mixed$model$smooth, "[[",  "label")
#' sapply(m_dis$model$smooth, "[[",  "label")
#' sapply(m_cont$model$smooth, "[[",  "label")
#'
#' # create new data for predict
#'
#' nd <- data.frame(covariate1 = c("a", "b", "c", "a"), covariate4 = c(0.4, 0.5, 0.1, 0.3), covariate2 = c("d", "d", "c", "d"), covariate3 = c(1, 0, 0.2, 2), covariate5 = c(0.2, 0.4, 1, 2))
#'
#' # predict mixed model, all terms on pdf-level without new data
#' p1 <- stats::predict(m_mixed, type =  "terms", level = "pdf")
#'
#' # predict  partial effects at clr-level of the mixed model for new data
#' p2 <- stats::predict(m_mixed, type =  "terms",  new_data = nd)
#'
#' # predict only the second partial effect at density-level of the mixed model for new data
#' p3 <- stats::predict(m_mixed, type =  "terms", which =  "ti(obs_density):covariate1b", new_data = nd, level = "pdf")
#'
#' # predict partial effects of the mixed model for new data and at density-level without the second term
#' p4 <- stats::predict(m_mixed, type =  "terms", exclude =  "ti(obs_density):covariate1b", new_data = nd, level = "pdf")
#'
#' # predict f_hat on density-level for new data
#' p5 <- stats::predict(m_mixed, type =  "pdf", new_data = nd)
#'
#' # predict clr(f_hat) for new data
#' p6 <- stats::predict(m_mixed, type =  "clr", new_data = nd)
#'
#' }
#'
#' @noRd
predict.densreg_obj <-
  function(object,
           new_data = NULL,
           which = NULL,
           exclude = NULL,
           type = "terms",
           level = "clr") {
    if (!isFALSE(object$params[[2]]) & !isFALSE(object$params[[1]])) {
      n_bins <- object$params$G + length(object$params[[2]])
    }
    if (isFALSE(object$params[[2]]) & !isFALSE(object$params[[1]])) {
      n_bins <- object$params$G
    }
    if (!isFALSE(object$params[[2]]) & isFALSE(object$params[[1]])) {
      n_bins <- length(object$params[[2]])
    }
    if (is.null(exclude) & is.null(which) | !type==  "terms") {
      needed <- as.list(attr(object$model$terms, "variables"))[-c(1:5)]
    }
    if (!is.null(exclude) & !is.null(which)) {
      warning(
        "The parameters exclude and which are both specified. exclude is primarily used and may overwrite information in which"
      )
    }
    all_terms <- sapply(object$model$smooth, "[[",  "label")
    if (!is.null(which) | !is.null(exclude)) {
      if (is.numeric(which)) {
        which <- all_terms[which]
      }
      if (is.numeric(exclude)) {
        exclude <- all_terms[exclude]
      }
    }
    if (is.null(exclude)) {
      if (is.null(which)) {
        exclude <- NULL
      } else{
        exclude <- all_terms[-match(which, all_terms)]
      }
    }
    if (is.null(which) & !is.null(exclude)) {
      which <- all_terms[-match(exclude, all_terms)]
    }

    if (!is.null(which) & type==  "terms") {
      needed <- c()
      for (trm in which) {
        for (cov in colnames(object$count_data)) {
          if (grepl(cov, trm)) {
            needed <- append(needed, cov)
          }
        }

      }
      needed <- unique(needed)[-1]
    }
    if (!is.null(new_data)) {
      if (!all(needed %in% colnames(new_data))) {
        stop("Not all needed covariates are included in new data!")
      }

      nd <- object$count_data[rep(c(1:n_bins), nrow(new_data))]
      relevant <- match(colnames(new_data), colnames(nd))
      j <- 1
      for (index in relevant) {
        nd[, index] <- rep(new_data[, j], each = n_bins)
        j <- j + 1
      }
    }
    else{
      nd <- object$count_data
    }
    if (type==  "terms") {
      prediction <-
        predict.gam(
          object$model,
          newdata = nd,
          exclude = exclude,
          type = "terms",
          se.fit = FALSE
        )[, -1]
      pred_list <- list()
      prediction <- data.frame(prediction)
      for (c in 1:ncol(prediction)) {
        pr <- matrix(prediction[, c], nrow = n_bins)
        if (level==  "pdf") {
          pr <- apply(
            pr,
            MARGIN = 2,
            FUN = clr,
            inverse = TRUE,
            w = nd$Delta[1:n_bins]
          )
        }
        pred_list[[length(pred_list) + 1]] <- pr

      }
      if (is.null(exclude)) {
        names(pred_list) <- c("intercept", all_terms[-1])
      } else{
        names(pred_list) <- c(which)
      }
      return(pred_list)
    } else{
      X <- stats::predict(object$model, type = "lpmatrix", newdata = nd)
      # remove intercepts per covariate combination (not of interest for estimated densities)
      intercepts <- which(grepl("group_id", colnames(X)))
      X <- X[, -intercepts]
      theta_hat <- object$model$coefficients[-intercepts]

      f_hat_clr <- X %*% theta_hat
      f_hat_clr <-
        matrix(f_hat_clr,
               nrow = n_bins)
      if (!is.null(ncol(f_hat_clr))) {
        f_hat <-
          apply(
            f_hat_clr,
            MARGIN = 2,
            FUN = clr,
            inverse = TRUE,
            w = nd$Delta[1:n_bins]
          )

      } else{
        f_hat <- clr(f_hat_clr, inverse = TRUE, w = nd$Delta[1:n_bins])

      }
      if (type==  "pdf") {
        return(f_hat)
      }
      if (type==  "clr") {
        return(f_hat_clr)
      }
    }
  }
