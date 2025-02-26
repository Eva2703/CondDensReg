#' A plot function for observed (histogram) count data of one covariate combination
#'
#' \code{\link{plot_single_histo}} can display the resulting
#' histogram or bar plot of a subset of data which was prepared by
#' \code{\link{data2counts}} which contains a vector of (weighted)
#' counts via a histogram on \eqn{I\setminus D} and (weighted) counts on \eqn{D}
#' where \eqn{I} is the interval of the continuous domain and \eqn{D} the set of
#' discrete values. It is meant to take only one group, i.e., one of the different
#' covariate combinations, into account. If no sample weights were included in the
#' data in the \code{\link{data2counts}}-step, the unweighted counts
#' are used instead of weighted counts.
#'
#' @encoding UTF-8
#' @importFrom dplyr "%>%" arrange mutate filter select distinct
#' @import data.table
#' @importFrom Rdpack reprompt
#' @param dta \code{data.table}-object which is an object of the sub-class
#' \code{histogram_count_data}, i.e., the output of the
#' \code{\link{data2counts}}-function. Subset in such a way, that
#' only one covariate combination, i.e. one \code{group_id} is contained.
#' @param case Plotting case as a \code{character}-object: can be specified as
#' "all", "continuous" or "discrete". Determines the components of the data which
#' will be plotted: "all" leads to a plot of the complete data, "continuous" only
#' plots the bins displaying data of the continuous domain \eqn{I} and "discrete"
#' shows only bins displaying the counts of the respective discrete values in \eqn{D}.
#' Please note that the plotting case is not the same as the specification of the
#' case in \code{\link{data2counts}} where the bins and discrete values
#' are defined. If the case is not "all" and \code{hist} is TRUE, the sum of the
#' areas under the bins and the discrete bars is not one. If missing (\code{NULL}),
#' the default is "all".
#' @param abs \code{logical} which determines if the resulting plot shows the absolute,
#' (weighted) counts (\code{abs=TRUE}) or the relative frequencies of the weighed
#' counts (\code{abs=FALSE}). If missing (\code{NULL}), the default is \code{FALSE}.
#' @param main Optional \code{character} passing a user-defined main title for the
#' resulting plot. If missing (\code{NULL}), the main will be generated automatically,
#' see also \code{automatic_main}.
#' @param automatic_main \code{logical} determining the automatic generation of a
#' main title for the plot which contains the number of bins, the plot case defined
#' by \code{case}, the values of the respective covariate combination and the group_id.
#' If \code{TRUE} and a user-defined \code{main} is passed, the automatic main is
#' overwriting the passed main. If missing (\code{NULL}), the default is \code{TRUE}.
#' @param ylim An interval (i.e., a \code{numeric} vector of length 2) specifying
#' the limits of the plot's y-axis. If missing (\code{NULL}), the limits are set
#' to 0 and the maximal height of a bar or histogram bin.
#' @param xlab Optional \code{character} passing a user-defined title of the x-axis.
#' If missing (\code{NULL}), the axis title is set to be the name of the density
#' variable of \code{dta}.
#' @param ylab Optional \code{character} passing a user-defined title of the y-axis.
#' If missing (\code{NULL}), the axis title contains the information if absolute
#' values or frequencies are shown in a histogram (for the continuous components)
#' or a bar plot.
#' @param col_discrete Name or hexadecimal RGB triplet of a color contained in
#' base R. Determines the color of the bars of the discrete values. If missing
#' (\code{NULL}), the default color is "red". It is recommended to use a color
#' contrasting the black outlines of the bins of the continuous domain especially
#' if one discrete value lies exactly at a bin border.
#' @param hist \code{logical} which determines if the resulting plot is a histogram
#' or a simple bar plot of the data. If \code{case="all"} and \code{hist=TRUE},
#' the areas under the bins and the height of discrete bars times the respective
#' weight of the discrete dirac measure sum up to one. Otherwise the height of
#' the bins is directly referring to the frequency or the (weighted) counts of
#' the respective bin without taking the bin width or the dirac measure into account.
#' If missing (\code{NULL}), the default is \code{TRUE}.
#' @param ... Other parameters which can be passed to the respective plotting function
#' specifying plotting parameters.
#'
#' @return The returned object is a plot object specified as explained above as
#' a histogram or bar plot of frequencies or absolute counts. The display of
#' mixed case count (histogram) data is as considered in the work of Maier et
#' al. (2025b).
#'
#' @examples
#' set.seed(101)
#'
#' # create data where 0 and 1 are the discrete observations, values
#' # equal 2 are replaced below by drawing from a beta distribution
#'
#' dta <-data.frame(
#' obs_density = sample(0:2, 100, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
#' covariate1 = sample(c("a", "b"), 100, replace = TRUE),
#' covariate2 = sample(c("c", "d"), 100, replace = TRUE),
#' sample_weights = runif(100, 0, 2)
#' )
#' dta[which(dta$obs_density == 2), ]$obs_density <- rbeta(length(which(dta$obs_density == 2)),
#'                                                         shape1 = 3, shape2 = 3)
#'
#' # Create histogram count dataset for dta using 10 equidistant
#' # bins and default values for continuous domain, discrete
#' # values and discrete weights while considering a mixed case
#' # of continuous and discrete domains.
#'
#' histo_data <- data2counts(
#' dta,
#' var_vec = c("covariate1", "covariate2"),
#' y = "obs_density",
#' sample_weights = "sample_weights",
#' bin_number = 10
#' )
#'
#' # Plot a single histogram of the weighted frequencies for
#' # group_id 1. Show all observations.
#'
#' plot_single_histo(histo_data%>%filter(group_id==1))
#'
#' # Restrict the plot to only discrete observations, use as
#' # bar color blue and show the absolute weighted counts
#' # instead of the frequencies:
#'
#' plot_single_histo(
#' histo_data %>% filter(group_id == 1),
#' case = "discrete",
#' abs = TRUE,
#' col_discrete = "blue"
#' )
#'
#' # Restrict now the plot to only continuous observations
#' # and plot a bar plot without the histogram property and
#' # with user-defined main, axis titles and axis limits
#'
#' plot_single_histo(
#'   histo_data %>% filter(group_id == 1),
#'   case = "continuous",
#'   hist = FALSE,
#'   main ="Title",
#'   automatic_main = FALSE,
#'   xlab = "x-axis",
#'   ylab ="y-axis",
#'   ylim = c(0,1)
#'   )
#'
#' @noRd
#
#' @references
#' Maier, E.-M., Fottner, A., Stoecker, A., Okhrin, Y., & Greven, S. (2025b):
#' Conditional density regression for individual-level data.
#' arXiv preprint arXiv:XXXX.XXXXX.
#'
plot_single_histo <- function(dta, case = "all", abs = FALSE, main = NULL,
                      automatic_main = TRUE, ylim = NULL, xlab = NULL,
                      ylab = NULL, col_discrete = "red", hist = TRUE, ...) {
  # library(data.table)
  # library(stringr)
  # initialize important variables
  n_bins <- sum((!dta$discrete))
  x_values_discrete <- NULL
  # if data2counts was not weighted, dta has no column "weighted_counts"
  # find names and values of response variable and covariate columns for both cases
  if ("weighted_counts" %in% colnames(dta)) {
    values_discrete <- distinct(dta[dta$discrete, 3])
    response_name <- colnames(dta)[3]
    var_ind <- 4:(ncol(dta) - 5)
    var_names <- colnames(dta)[c(4:(ncol(dta) - 5))]
    var_values <- dta[1, ..var_ind]
  }
  else{
    response_name <- colnames(dta)[2]
    dta$weighted_counts <- dta$counts
    var_ind <- c(3:(ncol(dta) - 6))
    var_names <- colnames(dta)[c(3:(ncol(dta) - 6))]
    var_values <- dta[1, ..var_ind]
  }
  # check if passed parameters are valid
  if (!is.null(case) &
      !case %in% c("all", "continuous", "discrete")) {
    stop("Invalid case parameter")
  }
  if (!is.null(abs) & !abs %in% c(TRUE, FALSE)) {
    stop("Invalid abs parameter")
  }
  if (!is.null(hist) & !hist %in% c(TRUE, FALSE)) {
    stop("Invalid hist parameter")
  }
  if (!is.null(main) & !is.character(main)) {
    stop("Invalid main parameter")
  }
  if (!is.null(xlab) & !is.character(xlab)) {
    stop("Invalid xlab parameter")
  }
  if (!is.null(ylab) & !is.character(ylab)) {
    stop("Invalid ylab parameter")
  }
  if (!is.null(ylim) & !(is.numeric(ylim) & length(ylim) == 2)) {
    stop("Invalid ylim parameter")
  }
  if (!is.null(main) & automatic_main) {
    warning("automatic_main = TRUE overwrites argument passed to main.")
  }
  if (length(unique(dta$group_id))>1) {
    warning("More than one covariate combination.")
  }
  # if not specified use name of y as x-axis label
  if (is.null(xlab)) {
    xlab <- response_name
  }
  # only discrete values in the data set -> create one pseudo bin with height zero
  if (nrow(dta%>%filter(discrete==FALSE))==0){
    pseudo_col<-as.data.frame(t(rep(0,ncol(dta))))
    n_bins<-1
    colnames(pseudo_col)<-colnames(dta)
    pseudo_col$group_id<-dta$group_id[1]
    pseudo_col[response_name]<-min(dta[,..response_name])+1/2*(max(dta[,..response_name])-min(dta[,..response_name]))
    pseudo_col$Delta<-1/2*(max(dta[,..response_name])-min(dta[,..response_name]))
    dta<-rbind(dta,pseudo_col)
    boarders<-c(min(dta[,..response_name]), max(dta[,..response_name]))
  }
  else{
    # define (n_bins+1) boarders of the bins for continuous components
    boarders <-
      append(distinct(dta[discrete == FALSE, ..response_name] - 1 / 2 * dta[discrete ==FALSE, Delta]),
             dta[discrete == FALSE, ..response_name][n_bins] + 1 / 2 * dta[discrete == FALSE, Delta][n_bins])
    boarders <- c(boarders[[1]], boarders[2])
  }
  # only continuous values in the data set -> create one discrete value with height zero
  if (nrow(dta%>%filter(discrete==TRUE))==0){
    pseudo_col<-as.data.frame(t(rep(0,ncol(dta))))
    colnames(pseudo_col)<-colnames(dta)
    pseudo_col$group_id<-dta$group_id[1]
    pseudo_col[response_name]<-min(dta[,..response_name])+1/2*(max(dta[,..response_name])-min(dta[,..response_name]))
    pseudo_col$Delta<-1/2*(max(dta[,..response_name])-min(dta[,..response_name]))
    pseudo_col$discrete<-TRUE
    dta<-rbind(dta,pseudo_col)
  }

  # sum_weights as normalizing factor if (weighted) counts are transformed to frequencies
  sum_weights <- sum(dta$weighted_counts)
  # isFALSE(abs) <=> plot frequencies instead of absolute counts
  if (isFALSE(abs)) {
    dta$freq <- dta$weighted_counts / sum_weights
    # frequency case for case="all"
    if (case == "all") {
      # generate automatic main
      if (automatic_main) {
        main <-
          paste0("All observations, number of bins: ", n_bins, "\n")
        x <- ""
        for (i in 1:length(var_values)) {
          x <- paste0(x, var_names[i], ": ", var_values[[i]], ", ")
        }
        x <- stringr::str_sub(x, end = -3)
        main <-
          paste0(main, x, "\n", "group_id: ", dta[1, group_id])
      }
      # dta$height needed for histogram plot
      ## Delta= bin_width (for continuous values) or = weight of Dirac measure (for discrete values)
      dta$height <- dta$freq / dta$Delta
      # histograms only for continuous components
      plot_hist <- dta %>% filter(discrete == FALSE)
      # histogram case:
      if (hist) {
        # create vectors of the heights for discrete values
        x_values_discrete <-
          c(unname(dta[discrete == TRUE, ..response_name]))[[1]]
        y_values_discrete <-
          c(unname(dta %>% filter(discrete == TRUE) %>% select(height)))[[1]]
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(plot_hist$height, y_values_discrete))
          if (is.na(ylim[2])){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "frequency (histogram)"
        }
        # barplot with histogram properties (height*Delta=frequency)
        graphics::barplot(plot_hist$height, main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim,
                ...)
        # add bars of discrete values
        graphics::segments(x0 = x_values_discrete-min(plot_hist[,..response_name])+1/2*plot_hist$Delta[1],
                 y0 = ylim[1],
                 x1 = x_values_discrete-min(plot_hist[,..response_name])+1/2*plot_hist$Delta[nrow(plot_hist)],
                 y1 = y_values_discrete,col = col_discrete)
        tcks<-round(unlist(plot_hist[,..response_name])-min(plot_hist[,..response_name])+1/2*plot_hist$Delta,digits=2)
        lbs<-round(unlist(plot_hist[,..response_name]),digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
      # barplot case (height=value of the frequency for this bin):
      if (!hist) {
        # create vectors of the frequencies for discrete values
        x_values_discrete <-
          c(unname(dta[discrete == TRUE, ..response_name]))[[1]]
        y_values_discrete <-
          c(unname(dta %>% filter(discrete == TRUE) %>% select(freq)))[[1]]
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(plot_hist$freq, y_values_discrete))
          if (is.na(ylim[2])){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "frequency (bar plot)"
        }
        # bar plot of frequencies for continuous domain
        graphics::barplot(plot_hist$freq, main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim, ...)
        # add bars of discrete values
        graphics::segments(x0 = x_values_discrete-min(plot_hist[,..response_name])+1/2*plot_hist$Delta[1],
                 y0 = ylim[1],
                 x1 = x_values_discrete-min(plot_hist[,..response_name])+1/2*plot_hist$Delta[nrow(plot_hist)],
                 y1 = y_values_discrete,col = col_discrete)
        tcks<-round(unlist(plot_hist[,..response_name])-min(plot_hist[,..response_name])+1/2*plot_hist$Delta,digits=2)
        lbs<-round(unlist(plot_hist[,..response_name]),digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
    }
    # only continuous domain is displayed
    if (case == "continuous") {
      # generate automatic main
      if (automatic_main) {
        main <-
          paste0("Continuous observations, number of bins: ", n_bins, "\n")
        x <- ""
        for (i in 1:length(var_values)) {
          x <- paste0(x, var_names[i], ": ", var_values[[i]], ", ")
        }
        x <- stringr::str_sub(x, end = -3)
        main <- paste0(main, x, "\n", "group_id: ", dta[1, group_id])
      }
      # dta$height needed for histogram plot
      ## Delta= bin_width (for continuous values) or = weight of Dirac measure (for discrete values)
      dta$height <- dta$freq / dta$Delta
      # histograms only for continuous components
      plot_hist <- dta %>% filter(discrete == FALSE)
      # histogram case:
      if (hist) {
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(plot_hist$height))
          if (is.na(ylim[2])){
            ylim<-c(0,1)
          }
          if (ylim[1] == ylim[2]){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "frequency (histogram)"
        }
        # barplot with histogram properties (height*Delta=frequency)
        graphics::barplot(plot_hist$height,main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim,
                ...)

        tcks<-round(unlist(plot_hist[,..response_name])-min(plot_hist[,..response_name])+1/2*plot_hist$Delta,digits=2)
        lbs<-round(unlist(plot_hist[,..response_name]),digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
      # barplot case (height=value of the frequency for this bin):
      if (!hist) {
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(plot_hist$freq))#
          if (is.na(ylim[2])){
            ylim<-c(0,1)
          }
          if (ylim[1] == ylim[2]){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "frequency (bar plot)"
        }
        # bar plot of frequencies for continuous domain
        graphics::barplot(plot_hist$freq, main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim, ...)

        tcks<-round(unlist(plot_hist[,..response_name])-min(plot_hist[,..response_name])+1/2*plot_hist$Delta,digits=2)
        lbs<-round(unlist(plot_hist[,..response_name]),digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
    }
    # only discrete values are displayed
    if (case == "discrete") {
      # generate automatic main
      if (automatic_main) {
        main <-
          paste0("Discrete observations, number of bins: ", n_bins, "\n")
        x <- ""
        for (i in 1:length(var_values)) {
          x <- paste0(x, var_names[i], ": ", var_values[[i]], ", ")
        }
        x <- stringr::str_sub(x, end = -3)
        main <-
          paste0(main, x, "\n", "group_id: ", dta[1, group_id])
      }
      # dta$height needed for histogram plot
      ## Delta= bin_width (for continuous values) or = weight of Dirac measure (for discrete values)
      dta$height <- dta$freq / dta$Delta
      # plot_hist$height for continuous domain set to zero (invisible in the plot)
      plot_hist <- dta %>% filter(discrete == FALSE)
      plot_hist$height <- 0
      plot_hist$freq <- 0
      # histogram case
      if (hist) {
        # create vectors of the heights for discrete values
        x_values_discrete <-
          c(unname(dta[discrete == TRUE, ..response_name]))[[1]]
        y_values_discrete <-
          c(unname(dta %>% filter(discrete == TRUE) %>% select(height)))[[1]]
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(y_values_discrete))
          if (is.na(ylim[2])){
            ylim<-c(0,1)
          }
          if (ylim[1] == ylim[2]){
            ylim<-c(0,1)
          }

        }
        if (is.null(ylab)) {
          ylab <- "frequency (histogram)"
        }
        # empty histogram of continous domain
        graphics::barplot(plot_hist$freq, main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim,
                xlim=c(min(x_values_discrete),max(x_values_discrete)),
                ...)
        # add bars of discrete values
        graphics::segments(x0 = x_values_discrete, y0 = ylim[1], x1 = x_values_discrete,
                 y1 = y_values_discrete,col = col_discrete)
        tcks<-round(x_values_discrete,digits=2)
        lbs<-round(x_values_discrete,digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
      # barplot case (height=value of the frequency for this bin):
      if (!hist) {
        # create vectors of the frequencies for discrete values
        x_values_discrete <-
          c(unname(dta[discrete == TRUE, ..response_name]))[[1]]
        y_values_discrete <-
          c(unname(dta %>% filter(discrete == TRUE) %>% select(freq)))[[1]]
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(y_values_discrete))
          if (is.na(ylim[2])){
            ylim<-c(0,1)
          }
          if (ylim[1] == ylim[2]){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "frequency (bar plot)"
        }
        # empty barplot of continuous domain
        graphics::barplot(plot_hist$freq, main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim, xlim=c(min(x_values_discrete),max(x_values_discrete)),
                ...)
        # add bars with discrete values frequencies
        graphics::segments(x0 = x_values_discrete, y0 = ylim[1], x1 = x_values_discrete,
                 y1 = y_values_discrete,col = col_discrete)
        tcks<-round(x_values_discrete,digits=2)
        lbs<-round(x_values_discrete,digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
    }
  }
  ### !isFALSE(abs) <=> display absolute (weighted) counts instead of frequencies
  if (!isFALSE(abs)) {
    # absolute count case for case="all"
    if (case == "all") {
      # generate automatic main
      if (automatic_main) {
        main <-
          paste0("All observations, number of bins: ", n_bins, "\n")
        x <- ""
        for (i in 1:length(var_values)) {
          x <- paste0(x, var_names[i], ": ", var_values[[i]], ", ")
        }
        x <- stringr::str_sub(x, end = -3)
        main <-
          paste0(main, x, "\n", "group_id: ", dta[1, group_id])
      }
      # dta$height needed for histogram plot
      ## Delta= bin_width (for continuous values) or = weight of Dirac measure (for discrete values)
      dta$height <- dta$weighted_counts / dta$Delta
      # histograms only for continuous components
      plot_hist <- dta %>% filter(discrete == FALSE)
      # histogram case:
      if (hist) {
        # create vectors of the heights for discrete values
        x_values_discrete <-
          c(unname(dta[discrete == TRUE, ..response_name]))[[1]]
        y_values_discrete <-
          c(unname(dta %>% filter(discrete == TRUE) %>% select(height)))[[1]]
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(plot_hist$height, y_values_discrete))
          if (ylim[1] == ylim[2]){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "counts (histogram)"
        }
        # barplot with histogram properties (height*Delta=count)
        graphics::barplot(plot_hist$height, main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim,
                ...)
        # add bars of discrete values
        graphics::segments(x0 = x_values_discrete-min(plot_hist[,..response_name])+1/2*plot_hist$Delta[1], y0 = ylim[1], x1 = x_values_discrete-min(plot_hist[,..response_name])+1/2*plot_hist$Delta[nrow(plot_hist)],
                 y1 = y_values_discrete,col = col_discrete)
        tcks<-round(unlist(plot_hist[,..response_name])-min(plot_hist[,..response_name])+1/2*plot_hist$Delta,digits=2)
        lbs<-round(unlist(plot_hist[,..response_name]),digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
      # barplot case (height=value of the frequency for this bin):
      if (!hist) {
        # create vectors of the counts for discrete values
        x_values_discrete <-
          c(unname(dta[discrete == TRUE, ..response_name]))[[1]]
        y_values_discrete <-
          c(unname(dta %>% filter(discrete == TRUE) %>% select(weighted_counts)))[[1]]
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(plot_hist$weighted_counts, y_values_discrete))
          if (ylim[1] == ylim[2]){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "counts (bar plot)"
        }
        # bar plot of counts for continuous domain
        graphics::barplot(plot_hist$weighted_counts,main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim,
                ...)
        # add bars of discrete values
        graphics::segments(x0 = x_values_discrete-min(plot_hist[,..response_name])+1/2*plot_hist$Delta[1], y0 = ylim[1], x1 = x_values_discrete-min(plot_hist[,..response_name])+1/2*plot_hist$Delta[nrow(plot_hist)],
                 y1 = y_values_discrete,col = col_discrete)
        tcks<-round(unlist(plot_hist[,..response_name])-min(plot_hist[,..response_name])+1/2*plot_hist$Delta,digits=2)
        lbs<-round(unlist(plot_hist[,..response_name]),digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
    }
    # only continuous domain is displayed
    if (case == "continuous") {
      # generate automatic main
      if (automatic_main) {
        main <-
          paste0("Continuous observations, number of bins: ", n_bins, "\n")
        x <- ""
        for (i in 1:length(var_values)) {
          x <- paste0(x, var_names[i], ": ", var_values[[i]], ", ")
        }
        x <- stringr::str_sub(x, end = -3)
        main <-
          paste0(main, x, "\n", "group_id: ", dta[1, group_id])
      }
      # dta$height needed for histogram plot
      ## Delta= bin_width (for continuous values) or = weight of Dirac measure (for discrete values)
      dta$height <- dta$weighted_counts / dta$Delta
      # histograms only for continuous components
      plot_hist <- dta %>% filter(discrete == FALSE)
      # histogram case
      if (hist) {
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(plot_hist$height))
          if (ylim[1] == ylim[2]){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "counts (histogram)"
        }
        # barplot with histogram properties (height*Delta=count)
        graphics::barplot(plot_hist$height, main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim,
                ...)
        tcks<-round(unlist(plot_hist[,..response_name])-min(plot_hist[,..response_name])+1/2*plot_hist$Delta,digits=2)
        lbs<-round(unlist(plot_hist[,..response_name]),digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
      # barplot case (height=count for this bin):
      if (!hist) {
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(plot_hist$weighted_counts))
          if (ylim[1] == ylim[2]){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "counts (bar plot)"
        }
        # bar plot of counts for continuous domain
        graphics::barplot(plot_hist$weighted_counts,main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim,
                ... )
        tcks<-round(unlist(plot_hist[,..response_name])-min(plot_hist[,..response_name])+1/2*plot_hist$Delta,digits=2)
        lbs<-round(unlist(plot_hist[,..response_name]),digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
    }
    # only discrete values are displayed
    if (case == "discrete") {
      # generate automatic main
      if (automatic_main) {
        main <-
          paste0("Discrete observations, number of bins: ", n_bins, "\n")
        x <- ""
        for (i in 1:length(var_values)) {
          x <- paste0(x, var_names[i], ": ", var_values[[i]], ", ")
        }
        x <- stringr::str_sub(x, end = -3)
        main <-
          paste0(main, x, "\n", "group_id: ", dta[1, group_id])
      }
      # dta$height needed for histogram plot
      ## Delta= bin_width (for continuous values) or = weight of Dirac measure (for discrete values)
      dta$height <- dta$weighted_counts / dta$Delta
      # plot_hist$height for continuous domain set to zero (invisible in the plot)
      plot_hist <- dta %>% filter(discrete == FALSE)
      plot_hist$weighted_counts <- 0
      plot_hist$height <- plot_hist$weighted_counts / plot_hist$Delta
      # histogram case
      if (hist) {
        # create vectors of the heights for discrete values
        x_values_discrete <-
          c(unname(dta[discrete == TRUE, ..response_name]))[[1]]
        y_values_discrete <-
          c(unname(dta %>% filter(discrete == TRUE) %>% select(height)))[[1]]
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(y_values_discrete))
          if (ylim[1] == ylim[2]){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "counts (histogram)"
        }
        # empty histogram of continous domain
        graphics::barplot(plot_hist$freq, main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim, xlim=c(min(x_values_discrete),max(x_values_discrete)),
                ...)
        # add bars of discrete values
        graphics::segments(x0 = x_values_discrete, y0 = ylim[1], x1 = x_values_discrete,
                 y1 = y_values_discrete,col = col_discrete)

        tcks<-round(x_values_discrete,digits=2)
        lbs<-round(x_values_discrete,digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
      # barplot case (height= count for this bin):
      if (!hist) {
        # create vectors of the counts for discrete values
        x_values_discrete <-
          c(unname(dta[discrete == TRUE, ..response_name]))[[1]]
        y_values_discrete <-
          c(unname(dta %>% filter(discrete == TRUE) %>% select(weighted_counts)))[[1]]
        if (is.null(ylim)) {
          ylim <- c(0, 1.2*max(y_values_discrete))
          if (ylim[1] == ylim[2]){
            ylim<-c(0,1)
          }
        }
        if (is.null(ylab)) {
          ylab <- "counts (bar plot)"
        }
        # empty histogram of continous domain
        graphics::barplot(plot_hist$freq, main = main, xlab = xlab, ylab = ylab,
                width = plot_hist$Delta,space = 0,ylim = ylim, xlim=c(min(x_values_discrete),max(x_values_discrete)),
                ...)
        # add bars of discrete values
        graphics::segments(x0 = x_values_discrete, y0 = ylim[1], x1 = x_values_discrete,
                 y1 = y_values_discrete,col = col_discrete)
        tcks<-round(x_values_discrete,digits=2)
        lbs<-round(x_values_discrete,digits=2)
        graphics::axis(1, at=tcks,labels=lbs )
      }
    }
  }
}


#' Saving plots of observed (histogram) count data
#'
#' \code{save_plots} offers the possibility to create plots for selected groups
#' with \code{\link{plot_single_histo}} and then save them locally as pdf.
#'
#' @encoding UTF-8
#'
#' @importFrom dplyr "%>%" arrange mutate filter select distinct
#' @import data.table
#' @importFrom Rdpack reprompt
#' @import png
#' @import cowplot
#' @import grid
#' @param dta \code{data.table}-object which is also an object of the sub-class
#' \code{histogram_count_data}, i.e., the output of the \code{\link{data2counts}}-function.
#' Subset in such a way, that only one covariate combination, i.e. one \code{group_id}
#' is contained.
#' @param selected_groups Numeric vector of \code{group_id}s for which the respective
#' plots will be created and saved. If missing (\code{NULL}) plots of all groups
#' within the range of the minimum and maximum \code{group_id} will be saved.
#' @param single \code{logical} indicating if all plots are saved as single PDFs
#' if \code{single=TRUE} or jointly in one document with four plots per page if
#' \code{single=FALSE}. If missing (\code{NULL}), the default is \code{TRUE}.
#' @param path Directory path where resulting plots will be saved. If not specified
#' by the user then the storage directory is the subdirectory "Images" of the current
#' working directory. If this does not yet exist, it will be created.
#' @param case Plotting case as a \code{character}-object: can be specified as
#' "all", "continuous" or "discrete". Determines the components of the data which
#' will be plotted: "all" leads to a plot of the complete data, "continuous" only
#' plots the bins displaying data of the continuous domain \eqn{I} and "discrete"
#' shows only bins displaying the counts of the respective discrete values in \eqn{D}.
#' Please note that the plotting case is not the same as the specification of the
#' case in \code{\link{data2counts}} where the bins and discrete values are defined.
#' If the case is not "all" and \code{hist} is TRUE, the sum of the areas under the
#' bins and the discrete bars is not one. If missing (\code{NULL}), the default is "all".
#' @param abs \code{logical} which determines if the resulting plot shows the absolute,
#' (weighted) counts (\code{abs=TRUE}) or the relative frequencies of the weighed
#' counts (\code{abs=FALSE}). If missing (\code{NULL}), the default is \code{FALSE}.
#' @param main Optional \code{character} passing a user-defined main title for the
#' resulting plot. If missing (\code{NULL}), the main will be generated automatically,
#' see also \code{automatic_main}.
#' @param automatic_main \code{logical} determining the automatic generation of a
#' main title for the plot which contains the number of bins, the plot case defined
#' by \code{case}, the values of the respective covariate combination and the group_id.
#' If \code{TRUE} and a user-defined \code{main} is passed, the automatic main is
#' overwriting the passed main. If missing (\code{NULL}), the default is \code{TRUE}.
#' @param ylim An interval (i.e., a \code{numeric} vector of length 2) specifying
#' the limits of the plot's y-axis. If missing (\code{NULL}), the limits are set
#' to 0 and the maximal height of a bar or histogram bin.
#' @param xlab Optional \code{character} passing a user-defined title of the x-axis.
#' If missing (\code{NULL}), the axis title is set to be the name of the density
#' variable of \code{dta}.
#' @param ylab Optional \code{character} passing a user-defined title of the y-axis.
#' If missing (\code{NULL}), the axis title contains the information if absolute
#' values or frequencies are shown in a histogram (for the continuous components)
#' or a bar plot.
#' @param col_discrete Name or hexadecimal RGB triplet of a color contained in
#' base R. Determines the color of the bars of the discrete values. If missing
#' (\code{NULL}), the default color is "red". It is recommended to use a color
#' contrasting the black outlines of the bins of the contunuous domain especially
#' if one discrete value lies exactly at a bin border.
#' @param file_names Vector of file names can be passed by the user. If \code{single=TRUE}
#' the number of file names has to be the same as the number of selected groups.
#' If \code{single=FALSE} only one file name can be passed. If missing (\code{NULL}),
#' the names of the single plots will be their \code{group_id} with the prefix "hist_".
#' The name of the file if \code{single=FALSE} is "all.pdf". A warning occurs if
#' \code{file_names_stem} is not \code{NULL} because only one of the two parameters
#' \code{file_names_stem} and \code{file_names} can be specified.
#' @param file_names_stem \code{character} which replaces the prefix "hist_" in
#' the automatically generated file name if \code{file_names=NULL}. If missing
#' (\code{NULL}), the names of the single plots will be their \code{group_id}
#' with the prefix "hist_". The name of the file if \code{single=FALSE} is "all.pdf".
#' A warning occurs if \code{file_names} is not \code{NULL} because only one of
#' the two parameters \code{file_names_stem} and \code{file_names} can be specified.
#' @param plot_width Value for the plot width in the case of \code{single=TRUE}.
#' For \code{single=FALSE} the plot width is not adjustable. If missing (\code{NULL}),
#' the default value if \code{single=TRUE} is as for the \code{pdf()}-function 7 inches.
#' @param plot_height Value for the plot height in the case of \code{single=TRUE}.
#' For \code{single=FALSE} the plot height is not adjustable. If missing (\code{NULL}),
#' the default value if \code{single=TRUE} is as for the \code{pdf()}-function 7 inches.
#' @param hist \code{logical} which determines if the resulting plot is a histogram
#' or a simple bar plot of the data. If \code{case="all"} and \code{hist=TRUE},
#' the areas under the bins and the height of discrete bars times the respective
#' weight of the discrete dirac measure sum up to one. Otherwise the height of the
#' bins is directly referring to the frequency or the (weighted) counts of the
#' respective bin without taking the bin width or the dirac measure into account.
#' If missing (\code{NULL}), the default is \code{TRUE}.
#' @param ... Other parameters which can be passed to the respective plotting
#' function specifying plotting parameters
#'
#' @return Saved plots in PDF-format.
#'
#' @noRd
#'
#' @examples
#' set.seed(101)
#'
#' # create data where 0 and 1 are the discrete observations, values
#' # equal 2 are replaced below by drawing from a beta distribution
#'
#' dta <-data.frame(
#' obs_density = sample(0:2, 100, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
#' covariate1 = sample(c("a", "b"), 100, replace = TRUE),
#' covariate2 = sample(c("c", "d"), 100, replace = TRUE),
#' sample_weights = runif(100, 0, 2)
#' )
#' dta[which(dta$obs_density == 2), ]$obs_density <- rbeta(length(which(dta$obs_density == 2)),
#'                                                         shape1 = 3, shape2 = 3)
#'
#' # Create histogram count dataset for dta using 10 equidistant
#' # bins and default values for continuous domain, discrete
#' # values and discrete weights while considering a mixed case
#' # of continuous and discrete domains.
#'
#' histo_data <- data2counts(
#' dta,
#' var_vec = c("covariate1", "covariate2"),
#' y = "obs_density",
#' sample_weights = "sample_weights",
#' bin_number = 10
#' )
#'
#' # Plot a single histogram of the weighted frequencies.
#' # Show all observations and save all plots in single PDFs.
#'
#' save_plots(histo_data)
#'
#' # Restrict the plot to only discrete observations
#' # of the groups 1-3, use as bar color blue and
#' # show the absolute weighted counts instead of
#' # the frequencies. Name the three resulting PDF
#' # with user-defined names.
#'
#' save_plots(
#' histo_data,
#' selected_groups= c(1,2,3),
#' case = "discrete",
#' abs = TRUE,
#' col_discrete = "blue",
#' file_names= c("one","two","three")
#' )
#'
#' # Restrict now the plot to only continuous observations
#' # and plot a bar plot without the histogram property and
#' # with user-defined main, axis titles and axis limits.
#' # Save all plots in one PDF.
#'
#' save_plots(
#'   histo_data,
#'   case = "continuous",
#'   hist = FALSE,
#'   main ="Title",
#'   automatic_main = FALSE,
#'   xlab = "x-axis",
#'   ylab ="y-axis",
#'   ylim = c(0,1),
#'   single=FALSE
#'   )
#'
save_plots <- function(dta, selected_groups = unique(dta$group_id), single = TRUE,
                       path = "./", case = "all", main = NULL, automatic_main = TRUE,
                       ylim = NULL, xlab = NULL, ylab = NULL, col_discrete = "red",
                       file_names = NULL, file_names_stem = NULL,
                       plot_width = ifelse(single, 7, 480), #default values for width and height in pdf()/png()
                       plot_height = ifelse(single, 7, 480), abs =FALSE, hist =TRUE, ...)  {
  ## check if passed parameters are valid
  # path
  if (!is.null(path) & !is.character(path)){
    stop("Invalid path parameter")
  }
  # selected groups format
  if (!is.null(selected_groups) & !(is.numeric(selected_groups))){
    stop("Invalid selected_groups parameter")
  }
  # selected groups range
  if (!is.null(selected_groups) & !(all(selected_groups%in%dta$group_id)) ){
    stop("Invalid range of the selected_groups parameter")
  }
  # single
  if (!is.null(single) & !single %in% c(TRUE, FALSE)){
    stop("Invalid single parameter")
  }
  # filenames
  if (!is.null(file_names) & !is.character(file_names)){
    stop("Invalid file_names parameter")
  }
  # filenamesstem
  if (!is.null(file_names_stem) & !is.character(file_names_stem)){
    stop("Invalid file_names_stem parameter")
  }
  # names+stem conflict
  if (!is.null(file_names_stem) & !is.null(file_names)){
    stop("Only one of the two parameters file_names_stem and file_names can be specified.")
  }
  # plot width valid
  if (!is.null( plot_width) & !(is.numeric( plot_width)&length( plot_width)==1)){
    stop("Invalid  plot_width parameter")
  }
  # plot height valid
  if (!is.null( plot_height) & !(is.numeric( plot_height)&length( plot_height)==1)){
    stop("Invalid  plot_height parameter")
  }
  # if path is not existing yet, create it
  if (!file.exists(path)) {
    dir.create(file.path(path))
  }
  # number of plots = number of selected groups
  n_plots <- length(selected_groups)
  # if single is TRUE, all plots are saved separately in pdfs, respective plot height and width used
  if (single) {
    for (k in selected_groups) {
      # create automatic file name if nothing is specified by the user: "hist_'k'.pdf"
      if (is.null(file_names_stem) & is.null(file_names)) {
        grDevices::pdf(paste0(path, "hist_", k, ".pdf"),
            width = plot_width,
            height = plot_height)
      }
      # create automatic file name if nothing the file name stem is specified by the user:
      # use stem instead of "hist_": "'stem''k'.pdf"
      if (!is.null(file_names_stem) & is.null(file_names)) {
        grDevices::pdf(
          paste0(path, file_names_stem, k, ".pdf"),
          width = plot_width,
          height = plot_height
        )
      }
      # if filenames are passed: use them
      if (!is.null(file_names)) {
        i<- which (k == selected_groups)
        grDevices::pdf(paste0(path, file_names[i], ".pdf"),
            width = plot_width,
            height = plot_height)
      }
      # select for the current group the subdata set
      dta_sg <- dta %>% filter(group_id == k)
      # create with plot_single_histo respective plot
      plot_single_histo(dta_sg, case, abs, main, automatic_main, ylim, xlab, ylab,
                col_discrete, hist, ...)
      grDevices::dev.off()
    }
  }
  # if single is FALSE, save all plots in one pdf file, 4 per page
  else{
    # if no filen ame is passed, use "all.pdf"
    if (is.null(file_names)) {
      grDevices::pdf(paste0(path, "all", ".pdf"))
    }
    # or use passed file name
    else{
      grDevices::pdf(paste0(path, file_names, ".pdf"))
    }
    # 2x2 grid per PDF-page
    graphics::par(mfrow=c(2,2))
    for (k in selected_groups) {
      # select for the current group the subdata set
      dta_sg <- dta %>% filter(group_id == k)
      # create with plot_single_histo respective plot
      plot_single_histo(dta_sg, case, abs, main, automatic_main, ylim, xlab, ylab,
                                  col_discrete, hist, ...)
    }
    grDevices::dev.off()
  }
}


#' Interactive application for plots for all covariate combinations
#'
#' \code{interactive_plots} creates for all groups, i.e., all different covariate combinations considered in the data, a plot based on  \code{\link{plot_single_histo}} and integrates them in an interactive application where the user can select a \code{group_id} with a slider. All specified parameters will be passed to \code{\link{plot_single_histo}}. \code{interactive_plots} is also used as default plotting method for data of the sub-class \code{histogram_count_data}, i.e., the output of the \code{\link{data2counts}}-function.
#' @encoding UTF-8
#' @importFrom dplyr "%>%" arrange mutate filter select distinct
#' @import data.table
#' @import stringr
#' @import manipulate
#' @importFrom Rdpack reprompt
#' @param dta \code{data.table}-object which is also an object of the sub-class \code{histogram_count_data}, i.e., the output of the \code{\link{data2counts}}-function with multiple covariate combinations, i.e. \code{group_id}s, contained.
#' @param selected_groups Numeric vector of \code{group_id}s for which the respective plots will be created and included in the interactive plot. If missing (\code{NULL}) plots of all groups within the range of the minimum and maximum \code{group_id} will be included.
#' @param case Plotting case as a \code{character}-object: can be specified as "all", "continuous" or "discrete". Determines the components of the data which will be plotted: "all" leads to a plot of the complete data, "continuous" only plots the bins displaying data of the continuous domain \eqn{I} and "discrete" shows only bins displaying the counts of the respective discrete values in \eqn{D}. Please note that the plotting case is not the same as the specification of the case in \code{\link{data2counts}} where the bins and discrete values are defined. If the case is not "all" and \code{hist} is TRUE, the sum of the areas under the bins and the discrete bars is not one. If missing (\code{NULL}), the default is "all".
#' @param abs \code{logical} which determines if the resulting plot shows the absolute, (weighted) counts (\code{abs=TRUE}) or the relative frequencies of the weighed counts (\code{abs=FALSE}). If missing (\code{NULL}), the default is \code{FALSE}.
#' @param main Optional \code{character} passing a user-defined main title for the resulting plot. If missing (\code{NULL}), the main will be generated automatically, see also \code{automatic_main}.
#' @param automatic_main \code{logical} determining the automatic generation of a main title for the plot which contains the number of bins, the plot case defined by \code{case}, the values of the respective covariate combination and the group_id. If \code{TRUE} and a user-defined \code{main} is passed, the automatic main is overwriting the passed main. If missing (\code{NULL}), the default is \code{TRUE}.
#' @param ylim An interval (i.e., a \code{numeric} vector of length 2) specifying the limits of the plot's y-axis. If missing (\code{NULL}), the limits are set to 0 and the maximal height of a bar or histogram bin.
#' @param xlab Optional \code{character} passing a user-defined title of the x-axis. If missing (\code{NULL}), the axis title is set to be the name of the density variable of \code{dta}.
#' @param ylab Optional \code{character} passing a user-defined title of the y-axis. If missing (\code{NULL}), the axis title contains the information if absolute values or frequencies are shown in a histogram (for the continuous components) or a bar plot.
#' @param col_discrete Name or hexadecimal RGB triplet of a color contained in base R. Determines the color of the bars of the discrete values. If missing (\code{NULL}), the default color is "red". It is recommended to use a color contrasting the black outlines of the bins of the contunuous domain especially if one discrete value lies exactly at a bin border.
#' @param hist \code{logical} which determines if the resulting plot is a histogram or a simple bar plot of the data. If \code{case="all"} and \code{hist=TRUE}, the areas under the bins and the height of discrete bars times the respective weight of the discrete dirac measure sum up to one. Otherwise the height of the bins is directly referring to the frequency or the (weighted) counts of the respective bin without taking the bin width or the dirac measure into account. If missing (\code{NULL}), the default is \code{TRUE}.
#' @param ... Other parameters which can be passed to the respective plotting function specifying plotting parameters.
#' @return The returned object is an interactive plot object specified as explained above as a histogram or bar plot of frequencies or absolute counts. The display of mixed case count (histogram) data is as considered in the work of Maier et al. (2025b). If the slider contains group IDs which are not available in the data, an empty plot and a warning occure.
#'
#'@examples
#' set.seed(101)
#'
#' # create data where 0 and 1 are the discrete observations, values
#' # equal 2 are replaced below by drawing from a beta distribution
#'
#' dta <-data.frame(
#' obs_density = sample(0:2, 100, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
#' covariate1 = sample(c("a", "b"), 100, replace = TRUE),
#' covariate2 = sample(c("c", "d"), 100, replace = TRUE),
#' sample_weights = runif(100, 0, 2)
#' )
#' dta[which(dta$obs_density == 2), ]$obs_density <- rbeta(length(which(dta$obs_density == 2)), shape1 = 3, shape2 = 3)
#'
#' # Create histogram count dataset for dta using 10 equidistant
#' # bins and default values for continuous domain, discrete
#' # values and discrete weights while considering a mixed case
#' # of continuous and discrete domains.
#'
#' histo_data <- data2counts(
#' dta,
#' var_vec = c("covariate1", "covariate2"),
#' y = "obs_density",
#' sample_weights = "sample_weights",
#' bin_number = 10
#' )
#'
#' # Plot a single histogram of the weighted frequencies.
#' # Show all observations.
#'
#' interactive_plots(histo_data)
#'
#' # Restrict the plot to only discrete observations
#' # of the groups 1-3, use as bar color blue and
#' # show the absolute weighted counts instead of
#' # the frequencies:
#'
#' interactive_plots(
#' histo_data,
#' selected_groups= c(1,2,3),
#' case = "discrete",
#' abs = TRUE,
#' col_discrete = "blue"
#' )
#'
#' # Restrict now the plot to only continuous observations
#' # and plot a bar plot without the histogram property and
#' # with user-defined main, axis titles and axis limits:
#'
#' interactive_plots(
#'   histo_data,
#'   case = "continuous",
#'   hist = FALSE,
#'   main ="Title",
#'   automatic_main = FALSE,
#'   xlab = "x-axis",
#'   ylab ="y-axis",
#'   ylim = c(0,1)
#'   )
#'
#' @noRd
#'
#' @references
#' Maier, E.-M., Fottner, A., Stoecker, A., Okhrin, Y., & Greven, S. (2025b):
#' Conditional density regression for individual-level data.
#' arXiv preprint arXiv:XXXX.XXXXX.
#'
interactive_plots <- function(dta, selected_groups = unique(dta$group_id), case = "all", abs =FALSE, main = NULL,
                              automatic_main = TRUE, ylim = NULL, xlab = NULL,
                              ylab =NULL, col_discrete = "red", hist =TRUE, ...){
  # selected groups format
  if (!is.null(selected_groups) & !(is.numeric(selected_groups))){
    stop("Invalid selected_groups parameter")
  }
  # selected groups range
  if (!is.null(selected_groups) & !(all(selected_groups%in%dta$group_id)) ){
    stop("Invalid range of the selected_groups parameter")
  }
  dta_selected <- dta%>%filter(group_id %in% selected_groups)
  # create manipulate of plots which were created with plot_single_histo for all included groups with passed plot parameters
  manipulate(plot_single_histo(dta_selected %>% filter(group_id == k), case, abs, main,
                       automatic_main , ylim, xlab, ylab, col_discrete, hist, ... ),
           k = slider( min = min(dta_selected$group_id), max = max(dta_selected$group_id) ) )
}



#' Plot function for the observed (histogram) count data
#'
#' \code{plot.histogram_count_data} is the default plot method for data of the class \code{histogram_count_data}. It is using the function \code{interactive_plots} with its default parameters, i.e. it is creating an interactive plot of the frequencies of the respective (weighted) counts with histograms on \eqn{I\setminus D} and of the counts on \eqn{D} with simple bars where \eqn{I} is the interval of the continuous domain and \eqn{D} the set of discrete values.. In the plot the areas under the bins and the height of discrete bars times the respective weight of the discrete dirac measure sum up to one.
#' @encoding UTF-8
#' @param x \code{data.table}-object which is an object of the sub-class \code{histogram_count_data}, i.e., the output of the \code{\link{data2counts}}-function with multiple covariate combinations, i.e. \code{group_id}s, contained.
#' @param ... Further parameters passed to \code{interactive_plots}.
#' @return Interactive plot displaying all histograms for the contained groups where the respective \code{group_id} is selectable with a slider. The display of mixed case count (histogram) data is as considered in the work of Maier et al. (2025b).
#'
#' @noRd
#'
#' @references
#' Maier, E.-M., Fottner, A., Stoecker, A., Okhrin, Y., & Greven, S. (2025b):
#' Conditional density regression for individual-level data.
#' arXiv preprint arXiv:XXXX.XXXXX.
#'
#' @examples
#' set.seed(101)
#'
#' # create data where 0 and 1 are the discrete observations, values
#' # equal 2 are replaced below by drawing from a beta distribution
#'
#' dta <-data.frame(
#' obs_density = sample(0:2, 100, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
#' covariate1 = sample(c("a", "b"), 100, replace = TRUE),
#' covariate2 = sample(c("c", "d"), 100, replace = TRUE),
#' sample_weights = runif(100, 0, 2)
#' )
#' dta[which(dta$obs_density == 2), ]$obs_density <- rbeta(length(which(dta$obs_density == 2)), shape1 = 3, shape2 = 3)
#'
#' # Create histogram count dataset for dta using 10 equidistant
#' # bins and default values for continuous domain, discrete
#' # values and discrete weights while considering a mixed case
#' # of continuous and discrete domains.
#'
#' histo_data <- data2counts(
#' dta,
#' var_vec = c("covariate1", "covariate2"),
#' y = "obs_density",
#' sample_weights = "sample_weights",
#' bin_number = 10
#' )
#'
#' # Plot histo_data:
#' plot(histo_data)

plot.histogram_count_data <- function(x, ...) {
  tmp <- as.data.table(x)
  interactive_plots(tmp)
}
