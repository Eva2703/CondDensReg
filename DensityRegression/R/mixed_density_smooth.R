#' P-Splines with integrate-to-zero constraint for density regression
#'
#' @description
#' Specifying \code{bs = "md"} in the tensor product smooth \code{\link[mgcv]{ti}}
#' (setting \code{mc = FALSE}) constructs a smoother to be used in \code{mgcv}'s
#' \code{\link[mgcv]{gam}} with \code{family = poisson()} for the case that the
#' response variable of interest is a density in a mixed Bayes Hilbert space
#' \eqn{B^2(\mu) = B^2(\mathcal{Y}, \mathcal{A}, \mu)}{B^2(\mu) = B^2(Y, A, \mu)}
#' (including continuous or discrete ones as special cases) as in Maier et al. (2025b).
#' The subset of the domain \eqn{\mathcal{Y}}{Y} corresponding to the continuous
#' part of the densities is denoted with \eqn{I}, the one corresponding to the
#' discrete part with \eqn{\mathcal{D}}{D}.
#' Corresponding covariate observations are the discrete values for the discrete
#' component and the midpoints of the bins underlying the histograms approximating
#' the densities for the continuous component. Response observations are the
#' counts of the discrete observations and the histograms (the count data can
#' easily be constructed from individual data with \code{\link{data2counts}}).
#' Specification of basis and penalization details for the continuous component
#' is as for \code{bs = "ps"}, see \code{mgcv}'s
#' \code{\link[mgcv]{smooth.construct.ps.smooth.spec}}, which our implementation
#' is oriented on.
#' Specification of the domain and for the discrete component is possible via
#' the argument \code{xt} in \code{\link[mgcv]{ti}} (see below). % \cr % Linebreak
#'
#' Note that the implementation includes the integrate-to-zero constraint of
#' \eqn{L^2_0}, i.e., further centering is not reasonable! This means, that the
#' constructor should not be used with \code{\link[mgcv]{s}} or
#' \code{\link[mgcv]{te}}, which always include centering (sum-to-zero) constraints
#' for all marginals. Instead, it is supposed to be used with \code{\link[mgcv]{ti}},
#' setting \code{mc = FALSE} for the corresponding component.
#' % (It would be desirable to include a warning, however, object does not seem
#' % to contain information on which kind of function it was called in.)
#'
#' We recommend to use \code{link{densreg}} to specify density regression models
#' (based on \code{smooth.construct.md.smooth.spec}), instead of via \code{\link[mgcv]{gam}}
#' directly, since the is quite cumbersome and specification has to be done with
#' extreme care to obtain a reasonable model.
#'
#' @details
#' The basis and penalty are constructed from a P-spline basis (continuous component)
#' respectively indicator functions with optional difference penalty (discrete
#' component) transformed to \eqn{L^2_0} as described in Appendix D of Maier et al.
#' (2025a) and embedded to the mixed Bayes Hilbert space as described in Section
#' 2.2 of Maier et al. (2025b).
#' The argument \code{xt} in \code{\link[mgcv]{ti}} is used for further specification
#' regarding the underlying Bayes Hilbert space. \code{xt} has to be a list of
#' the same length as the vector \code{bs} specifying the marginal bases. For all
#' \code{"md"} type marginal bases, the corresponding \code{xt}-list element is
#' again a list with the following elements:
#' % ti(..., xt = list(list(values_discrete = NULL,
#' %                        weights_discrete = NULL,
#' %                        domain_continuous = NULL,
#' %                        penalty_discrete = NULL)))
#' \itemize{
#' \item
#' \code{values_discrete}: Vector of values in \eqn{\mathcal{D}}{D} (the subset of
#' the domain corresponding to the discrete part of the densities). Defaults to
#' missing (\code{NULL}) in which case it is set to \code{c(0, 1)}. If set to
#' \code{FALSE}, the discrete component is considered to be empty, i.e., the
#' Lebesgue measure is used as reference measure (continuous special case).
#' \item
#' \code{weights_discrete}: Vector of weights for the Dirac measures corresponding
#' to \code{values_discrete}. If missing (\code{NULL}) it is set to 1 in all
#' components as default. Can be a scalar for equal weights for all discrete values
#' or a vector with specific weights for each corresponding discrete value.
#' \item
#' \code{domain_continuous}: An interval (i.e., a vector of length 2) specifying
#' \eqn{I}. If missing (\code{NULL}) it is set to \code{c(0, 1)} as default. If
#' set to \code{FALSE}, the continuous component is considered to be empty, i.e.,
#' a weighted sum of dirac measures is used as reference measure (discrete special
#' case).
#' \item
#' \code{penalty_discrete}: integer (or \code{NULL}) giving the order of
#' differences to be used for the penalty of the discrete component, with \code{0}
#' corresponding to the identity matrix as penalty matrix (analogously to \code{m[2]}
#' for the continuous component); Note that the order of differences must be smaller
#' than the number of values in the discrete component, i.e., the length of
#' \code{values_discrete}; if missing (\code{NULL}), in the mixed case, it is
#' set to a zero matrix (corresponding to no penalization), while in the discrete
#' case, it is set to a diagonal matrix (as in a ridge penalty). In the discrete
#' case, please set the corresponding argument of \code{sp} in \code{\link[mgcv]{ti}}
#' to \code{0} to estimate an unpenalized model.
#' % For marginal bases of type \code{"ps"}/\code{"bs"}, the corresponding xt-list
#' % element has to be set to \code{NULL}.
#' }
#'
#' @importFrom MASS Null
#' @importFrom splines splineDesign
#'
#' @param object a smooth specification object, usually generated by a term
#' \code{ti(x, bs="md", ...)}
#' @inheritParams mgcv::smooth.construct.ps.smooth.spec
#'
#' @return An object of class \code{"mdspline.smooth"}. See
#' \code{\link[mgcv]{smooth.construct}}, for the elements that this object will contain.
#'
#' @author Eva-Maria Maier, Alexander Fottner
#'
#' @references
#' Maier, E.-M., Fottner, A., Stoecker, A., Okhrin, Y., & Greven, S. (2025b):
#' Conditional density regression for individual-level data.
#' arXiv preprint arXiv:XXXX.XXXXX.
#'
#' Maier, E.-M., Stoecker, A., Fitzenberger, B., Greven, S. (2025a):
#' Additive Density-on-Scalar Regression in Bayes Hilbert Spaces with an Application to Gender Economics.
#' Annals of Applied Statistics, 19(1), ???-???.
#'
#' @examples
#' \donttest{
#' ### create data
#' set.seed(101)
#'
#' N <- 100
#' dta <- data.frame(covariate1 = sample(c("a", "b", "c"), N, replace = TRUE),
#'                   covariate2 = rnorm(n = N))
#' dta$obs_density <- sapply(seq_len(N),
#'                           function(i) {
#'                             a_0 <- ifelse(dta$covariate1[i] == "a", 0.1,
#'                                           ifelse(dta$covariate1[i] == "b", 0.2, 0.3))
#'                             p_0 <- a_0 * sin(dta$covariate2[i]) + a_0 + 0.05
#'                             a_1 <- ifelse(dta$covariate1[i] == "a", 0.25,
#'                                           ifelse(dta$covariate1[i] == "b", 0.15, 0.05))
#'                             p_1 <- a_1 * cos(dta$covariate2[i]) + a_1 + 0.1
#'                             sample(0:2, 1, prob = c(p_0, p_1, 1 - p_0 - p_1))
#'                           })
#'
#' dta$covariate1 <- ordered(dta$covariate1)
#'
#' # data for the mixed case
#' dta_mixed <- dta
#' ind_cont <- which(dta_mixed$obs_density == 2)
#' dta_mixed[ind_cont, ]$obs_density <-
#'   sapply(seq_along(ind_cont), function(i)
#'     rbeta(1, shape1 = 1 + exp(dta_mixed$covariate2[i]),
#'           shape2 = 1 + as.numeric(dta_mixed$covariate1[i])))
#' n_bins <- 20
#' dta_mixed <- data2counts(dta = dta_mixed, var_vec = c("covariate1", "covariate2"),
#'                         bin_number = n_bins, values_discrete = c(0, 1),
#'                         domain_continuous = c(0, 1))
#'
#' # data for the continuous case
#' dta_cont <- dta_mixed[which(!dta_mixed$discrete), ]
#'
#' # data for the discrete case
#' dta_dis <- data2counts(dta = dta, c("covariate1", "covariate2"),
#'                       values_discrete = c(0, 1, 2), domain_continuous = FALSE)
#'
#' ### fit model in the mixed case
#' # All marginals in density-direction are specified via the mixed density basis "md",
#' # which per default uses a mixed Bayes Hilbert space with continuous domain (0, 1)
#' # and discrete values at 0 and 1. Since this is exactly our scenario, we don't
#' # actually need to specify xt here.
#'
#' m_mixed <- gam(counts ~
#'                  # no scalar global intercept (we add a density-intercept instead)
#'                  - 1 +
#'                  # intercept (corresponding to reference covariate1 = "a")
#'                  ti(obs_density, bs = "md", m = list(c(2, 2)), k = 8,
#'                     mc = FALSE, np = FALSE) +
#'                  # group specific intercept for leveles "b" and "c" of covariate2
#'                  ti(obs_density, bs = "md", m = list(c(2, 2)), k = 8,
#'                     mc = FALSE, np = FALSE, by = covariate1) +
#'                  # smooth effect of covariate2 (modeled via P-splines)
#'                  ti(covariate2, obs_density, bs = c("ps", "md"),
#'                     m = list(c(2, 2), c(2, 2)), k = c(8, 8), mc = c(TRUE, FALSE),
#'                     np = FALSE) +
#'                  # intercepts per covariate combination (modeling absolute
#'                  # counts for poisson regression)
#'                  as.factor(group_id) +
#'                  # offsets: log(Delta) accounting for bin width, gam_offsets
#'                  # for weighted observations
#'                  offset(log(Delta)),
#'                family = poisson(), method = "REML", data = dta_mixed)
#'
#' ### fit model in the continuous case
#' # Continuous domain (0, 1) without discrete values
#' xt_c <- list(values_discrete = FALSE, domain_continuous = c(0, 1))
#'
#' m_cont <- gam(counts ~
#'                 # no scalar global intercept (we add a density-intercept instead)
#'                 - 1 +
#'                 # intercept (corresponding to reference covariate1 = "a")
#'                 ti(obs_density, bs = "md", m = list(c(2, 2)), k = 8,
#'                    mc = FALSE, np = FALSE, xt = list(xt_c)) +
#'                 # group specific intercept for leveles "b" and "c" of covariate2
#'                 ti(obs_density, bs = "md", m = list(c(2, 2)), k = 8,
#'                    mc = FALSE, np = FALSE, xt = list(xt_c), by = covariate1) +
#'                 # smooth effect of covariate2 (modeled via P-splines)
#'                 ti(covariate2, obs_density, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
#'                    k = c(8, 8), mc = c(TRUE, FALSE), np = FALSE,
#'                    xt = list(NULL, xt_c)) +
#'                 # intercepts per covariate combination (modeling absolute
#'                 # counts for poisson regression)
#'                 as.factor(group_id) +
#'                 # offsets: log(Delta) accounting for bin width, gam_offsets
#'                 # for weighted observations
#'                 offset(log(Delta)),
#'               family = poisson(), method = "REML", data = dta_cont)
#'
#' ### fit model in the discrete case
#' # Continuous domain empty, discrete values 0, 1, 2
#' xt_d <- list(domain_continuous = FALSE, values_discrete = c(0, 1, 2))
#'
#' m_dis <- gam(counts ~
#'                # no scalar global intercept (we add a density-intercept instead)
#'                - 1 +
#'                # intercept (corresponding to reference covariate1 = "a")
#'                ti(obs_density, bs = "md", m = list(c(2, 2)), k = 8,
#'                   mc = FALSE, np = FALSE, xt = list(xt_d), sp = 0) +
#'                # group specific intercept for leveles "b" and "c" of covariate2
#'                ti(obs_density, bs = "md", m = list(c(2, 2)), k = 8,
#'                   mc = FALSE, np = FALSE, xt = list(xt_d), sp = 0, by = covariate1) +
#'                # smooth effect of covariate2 (modeled via P-splines)
#'                ti(covariate2, obs_density, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
#'                   k = c(8, 8), mc = c(TRUE, FALSE), np = FALSE,
#'                   xt = list(NULL, xt_d), sp = c(-1, 0)) +
#'                # intercepts per covariate combination (modeling absolute
#'                # counts for poisson regression)
#'                as.factor(group_id) +
#'                # offsets: log(Delta) accounting for bin width, gam_offsets
#'                # for weighted observations
#'                offset(log(Delta)),
#'              family = poisson(), method = "REML", data = dta_dis)
#' }
#'
#' @export
smooth.construct.md.smooth.spec <- function (object, data, knots) {
  if (is.list(object$xt[[1]])) {
    if (length(object$xt) == 1) {
      object$xt <- object$xt[[1]]
    } else {
      stop("xt is not specified correctly.")
    }
  }
  discrete <- ifelse(isFALSE(object$xt$values_discrete), FALSE, TRUE)
  continuous <- ifelse(isFALSE(object$xt$domain_continuous), FALSE, TRUE)
  if (!discrete & !continuous) {
    stop("Discrete and continuous components are both empty!")
  }
  x <- data[[object$term]]
  # getting specifications for continuous component
  if (continuous) {
    if (length(object$p.order) == 1)
      m <- rep(object$p.order, 2)
    else m <- object$p.order  # m[1] - basis order, m[2] - penalty order
    m[is.na(m)] <- 2 ## default
    object$p.order <- m
    if (object$bs.dim < 0)
      object$bs.dim <- max(10, m[1] + 1) ## default
    nk <- object$bs.dim - m[1] # basis dimension - order of spline -> number of interior knots for continuous component
    if (nk <= 0)
      stop("Basis dimension too small for b-spline order")
    if (length(object$term) != 1)
      stop("Basis only handles 1D smooths")
    cont_dim <- object$bs.dim
  }

  # initializing discrete component
  if (is.null(object$xt$values_discrete)) {
    object$xt$values_discrete <- c(0, 1)
  }
  t_discrete <- object$xt$values_discrete
  if (is.null(object$xt$weights_discrete)) {
    object$xt$weights_discrete <- rep(1, length(t_discrete))
  }
  w_discrete <- object$xt$weights_discrete
  if (length(t_discrete) != length(w_discrete)) {
    stop("Lengths of values_discrete and weights_discrete have to be the same.")
  }

  # initializing continuous component
  if (!discrete) {
    cont_positions <- seq_along(x)
    x_cont <- x
    if (is.null(object$xt$domain_continuous)) {
      object$xt$domain_continuous <- range(x)
    }
    discrete_dim <- 0
    design_discrete <- Z_discrete <- D_discrete <- NULL
  } else if (!continuous) { # discrete case
    if (any(!(x %in% t_discrete))) {
      stop("Given domain (which is discrete) does not include all observations!")
    }
    discrete_dim <- length(t_discrete)
    x_discrete <- x
    design_cont <- Z_cont <- D_cont <- NULL
    object$bs.dim <- discrete_dim
    cont_dim <- 0
    m <- c(0, 0)
    # cont_positions <- NULL
    # x_cont <- NULL
  } else {
    cont_positions <- which(!(x %in% t_discrete))
    x_cont <- x[cont_positions]
    if (is.null(object$xt$domain_continuous)) {
      object$xt$domain_continuous <- range(t_discrete, x)
    }
  }

  if (continuous) {
    if (object$xt$domain_continuous[1] > min(x_cont) ||
        object$xt$domain_continuous[2] < max(x_cont)) {
      stop("Given domain does not include data corresponding to continuous component.")
    }
  }

  # in mixed case: add discrete value t_{D+1} and weight corresponding to the continuous component
  if (continuous & discrete) {
    t_discrete[length(t_discrete) + 1] <- range(t_discrete, x)[2] + 1
    w_discrete[length(w_discrete) + 1] <- diff(object$xt$domain_continuous)
    discrete_dim <- length(t_discrete)
    # set all continuous values to t_{D+1} for discrete basis
    x_discrete <- x
    x_discrete[cont_positions] <- max(t_discrete)
    object$bs.dim <- object$bs.dim + discrete_dim # combined dimension (before implementing constraints!)
  }

  if (length(unique(x)) < object$bs.dim)
    warning("Basis dimension is larger than number of (histogram/discrete) counts.")

  if (continuous) {
    k <- sort(knots[[object$term]])
    if (is.null(k)) {
      xl <- object$xt$domain_continuous[1]
      xu <- object$xt$domain_continuous[2]
    } else if (length(k) == 2) {
      xl <- min(k)
      xu <- max(k)
      if (xl > min(x_cont) || xu < max(x_cont))
        stop("Knot range does not include data corresponding to continuous component.")
    }
    if (is.null(k) || length(k) == 2) {
      # nk = basis dimension - order of spline
      xr <- xu - xl
      # xl <- xl - xr * 0.001 ### E: this shifts the knots slightly outwards, i.e., not (exactly) the supplied knots (for length(k) == 2) are used
      # xu <- xu + xr * 0.001 ### E: I commented it out, since for me it's unclear why not to use the supplied knots
      dx <- (xu - xl)/(nk - 1)
      k <- seq(xl - dx * (m[1] + 1), xu + dx * (m[1] + 1),
               length = nk + 2 * m[1] + 2)
    } else {
      if (length(k) != nk + 2 * m[1] + 2)
        stop(paste("There should be ", nk + 2 * m[1] + 2, " supplied knots"))
    }
    ord <- m[1] + 2
    # if (k[ord] != object$xt$domain_continuous[1] ||
    #     k[length(k) - (ord - 1)] != object$xt$domain_continuous[2]) {
    if (isFALSE(all.equal(k[c(ord, length(k) - (ord - 1))],
                          object$xt$domain_continuous))) {
      warning("Knots do not match domain of continuous component.")
    }
    if (is.null(object$deriv)) {
      object$deriv <- 0
    } else if (object$deriv != 0) {
      warning("The mixed density smoother is not intended to be used for derivatives. Reasonable behavior is only guaranteed for deriv = 0.")
    }
    # construct design matrices from transformed B-splines integrating to zero (see
    # Appendix B of Maier et al., 2025a, based on Wood, 2017, Section 1.8.1) and
    # apply embedding to combine them to one design matrix of a mixed basis (Maier
    # et al., 2025a, Proposition A.4; Maier et al., 2022, Section 2.2)
    C_cont <- sapply(1:(length(k) - ord), function(j) stats::integrate(function(x)
      splines::splineDesign(knots = k, x, ord = ord, derivs = object$deriv)[, j],
      lower = k[ord], upper = k[length(k) - (ord - 1)])$value)
    Z_cont <- MASS::Null(C_cont)
    design_cont <- matrix(0, nrow = length(x), ncol = ncol(Z_cont))
    design_cont[cont_positions, ] <- splines::splineDesign(knots = k, x_cont, ord = ord,
                                                           derivs = object$deriv) %*% Z_cont

    if (is.null(object$mono))
      object$mono <- 0
    if (object$mono != 0) {
      stop("SCOP splines are not supported yet!")
    } else {
      # construct penalty matrices including necessary transformation (see Appendix
      # B of Maier et al., 2021) and combine them to one penalty (block) matrix
      # for the mixed design matrix
      if (m[2] > 0) {
        D_cont <- diff(diag(cont_dim), differences = m[2])
      } else {
        D_cont <- diag(cont_dim)
      }
      S_cont <- t(Z_cont) %*% crossprod(D_cont) %*% Z_cont
    }
  }

  if (discrete) {
    k_discrete <- sapply(seq_len(length(t_discrete) + 1),
                         function(j) mean(c(min(t_discrete) - 1, t_discrete, max(t_discrete) + 1)[j:(j+1)]))
    object$knots_discrete <- k_discrete
    object$values_discrete <- t_discrete
    # integrals of basis functions are equal to weights
    C_discrete <- object$weights_discrete <- w_discrete
    Z_discrete <- MASS::Null(C_discrete)
    design_discrete <- splines::splineDesign(k_discrete, x_discrete, 1)  %*% Z_discrete

    if (is.null(object$xt$penalty_discrete)) {
      if (continuous) {
        # default in mixed case: discrete component unpenalized
        D_discrete <- matrix(0, nrow = discrete_dim, ncol = discrete_dim)
      } else {
        # default in discrete case: Ridge penalty; Unpenalized unfortunately not
        # possible, since gam() still computes an eigenvalue decomposition later,
        # which yields an error for a 0-matrix. Thus an unpenalized model has to
        # be specified by setting the sp-argument of ti() to 0 (in the corresponding
        # component)
        D_discrete <- diag(discrete_dim)
      }
    } else if (object$xt$penalty_discrete > 0) {
      if (object$xt$penalty_discrete >= discrete_dim) {
        warning(paste0("Order of differences for discrete penalty must be smaller than number of discrete values. Maximal possible value ",
                       discrete_dim - 1, " is used for penalty_discrete instead of specified value ",
                       object$xt$penalty_discrete, "."))
        object$xt$penalty_discrete <- discrete_dim - 1
      }
      D_discrete <- diff(diag(discrete_dim), differences = object$xt$penalty_discrete)
    } else if (object$xt$penalty_discrete == 0) {
      D_discrete <- diag(discrete_dim)
    } else {
      warning("penalty_discrete has to be a non-negative integer. No penalty for discrete component is used.")
      D_discrete <- matrix(0, nrow = discrete_dim, ncol = discrete_dim)
    }
    S_discrete <- t(Z_discrete) %*% crossprod(D_discrete) %*% Z_discrete
  }

  object$X <- cbind(design_cont, design_discrete) # combine both design matrices to get final design matrix
  object$Z <- list(cont = Z_cont, discrete = Z_discrete)

  if (!discrete) {
    if (!is.null(k)) {
      if (sum(colSums(object$X) == 0) > 0)
        warning("There is *no* information about some basis coefficients")
    }
  }

  object$D <- list(continuous = D_cont, discrete = D_discrete)

  if (!discrete) {
    S <- S_cont
    discrete_dim <- 0
    object$xt$penalty_discrete <- NULL
  } else if (!continuous) {
    S <- S_discrete
  } else {
    # combine penalties in a block matrix
    S <- cbind(rbind(S_cont, matrix(0, ncol = ncol(S_cont), nrow = nrow(S_discrete))),
               rbind(matrix(0, nrow = nrow(S_cont), ncol = ncol(S_discrete)), S_discrete))
  }

  object$S <- list(S)
  object$rank <- Matrix::rankMatrix(S)[1]
  if (length(object$S) != 1) {
    warning("Length of penalties does not match length of marginals (which should be one)")
  }
  object$null.space.dim <- ifelse(is.null(object$xt$penalty_discrete), m[2],
                                  m[2] + object$xt$penalty_discrete)
  if (continuous) {
    object$knots <- k
    object$m <- m
  }
  object$type <- ifelse(!continuous, "discrete",
                        ifelse(!discrete, "continuous", "mixed"))
  class(object) <- "mdspline.smooth"
  return(object)
}

#' Predict matrix method function for mixed density smooth
#'
#' \code{\link[mgcv]{Predict.matrix}} method function for smooth class
#' \code{mdspline.smooth} to enable prediction from a model fitted with
#' \code{mgcv}'s \code{\link[mgcv]{gam}}.
#'
#' The Predict matrix function is not normally called directly, but is rather
#' used internally by \code{mgcv}'s \code{\link[mgcv]{predict.gam}} etc. to
#' predict from a fitted \code{\link[mgcv]{gam}} model. See \code{\link[mgcv]{Predict.matrix}}
#' for more details, or \code{\link{smooth.construct.md.smooth.spec}} for details
#' on the mixed density smooth \code{"md"}.
#'
#' @importFrom splines splineDesign
#'
#' @param object a smooth specification object, usually generated by a term
#' \code{ti(x,bs="md",...)}
#' @inheritParams mgcv::Predict.matrix.cr.smooth
#'
#' @return A matrix mapping the coefficients for the smooth term to its values at
#' the supplied data values.
#'
#' @author Eva-Maria Maier, Alexander Fottner
#'
#' @references
#' Maier, E.-M., Fottner, A., Stoecker, A., Okhrin, Y., & Greven, S. (2025b):
#' Conditional density regression for individual-level data.
#' arXiv preprint arXiv:XXXX.XXXXX.
#'
#' @examples
#' \donttest{
#' ### create data
#' set.seed(101)
#'
#' N <- 100
#' dta <- data.frame(covariate1 = sample(c("a", "b", "c"), N, replace = TRUE),
#'                   covariate2 = rnorm(n = N))
#' dta$obs_density <- sapply(seq_len(N),
#'                           function(i) {
#'                             a_0 <- ifelse(dta$covariate1[i] == "a", 0.1,
#'                                           ifelse(dta$covariate1[i] == "b", 0.2, 0.3))
#'                             p_0 <- a_0 * sin(dta$covariate2[i]) + a_0 + 0.05
#'                             a_1 <- ifelse(dta$covariate1[i] == "a", 0.25,
#'                                           ifelse(dta$covariate1[i] == "b", 0.15, 0.05))
#'                             p_1 <- a_1 * cos(dta$covariate2[i]) + a_1 + 0.1
#'                             sample(0:2, 1, prob = c(p_0, p_1, 1 - p_0 - p_1))
#'                           })
#'
#' dta$covariate1 <- ordered(dta$covariate1)
#'
#' dta_mixed <- dta
#' ind_cont <- which(dta_mixed$obs_density == 2)
#' dta_mixed[ind_cont, ]$obs_density <-
#'   sapply(seq_along(ind_cont), function(i)
#'     rbeta(1, shape1 = 1 + exp(dta_mixed$covariate2[i]),
#'           shape2 = 1 + as.numeric(dta_mixed$covariate1[i])))
#' n_bins <- 20
#' dta_mixed <- data2counts(dta = dta_mixed, var_vec = c("covariate1", "covariate2"),
#'                         bin_number = n_bins, values_discrete = c(0, 1),
#'                         domain_continuous = c(0, 1))
#'
#' ### fit model
#' m_mixed <- gam(counts ~
#'                  # no scalar global intercept (we add a density-intercept instead)
#'                  - 1 +
#'                  # intercept (corresponding to reference covariate1 = "a")
#'                  ti(obs_density, bs = "md", m = list(c(2, 2)), k = 8,
#'                     mc = FALSE, np = FALSE) +
#'                  # group specific intercept for leveles "b" and "c" of covariate2
#'                  ti(obs_density, bs = "md", m = list(c(2, 2)), k = 8,
#'                     mc = FALSE, np = FALSE, by = covariate1) +
#'                  # smooth effect of covariate2 (modeled via P-splines)
#'                  ti(covariate2, obs_density, bs = c("ps", "md"),
#'                     m = list(c(2, 2), c(2, 2)), k = c(8, 8), mc = c(TRUE, FALSE),
#'                     np = FALSE) +
#'                  # intercepts per covariate combination (modeling absolute
#'                  # counts for poisson regression)
#'                  as.factor(group_id) +
#'                  # offsets: log(Delta) accounting for bin width, gam_offsets
#'                  # for weighted observations
#'                  offset(log(Delta)),
#'                family = poisson(), method = "REML", data = dta_mixed)
#'
#' ### Functions using S3 Predict.matrix-method for class 'mdspline.smooth' are, e.g.,
#' ### model.matrix and predict
#'
#' # Using model.matrix (and internally Predict.matrix.mdspline.smooth) to extract
#' # design matrix and construct estimated densities thereof
#' X <- model.matrix(m_mixed)
#' # The design matrix includes intercepts per covariate combination modeling the
#' # absolute counts for the Poisson model, which are not of interest on density-level.
#' # Thus, we remove them.
#' intercepts <- which(grepl("group_id", colnames(X)))
#' X <- X[, -intercepts]
#' theta_hat <- m_mixed$coefficients[-intercepts]
#' # compute estimated conditional clr-transformed densities
#' f_hat_clr <- matrix(c(X %*% theta_hat), nrow = length(unique(dta_mixed$obs_density)))
#' f_hat <- apply(f_hat_clr, 2,
#'                FDboost::clr, w = c(1, rep(1/n_bins, n_bins), 1), inverse = TRUE)
#' t <- unique(dta_mixed$obs_density)
#' matplot(t[2:(length(t) - 1)], f_hat[2:(length(t) - 1), ], type = "l",
#'         col = rainbow(ncol(f_hat)), lty = rep(1:5, ceiling(ncol(f_hat) / 5)),
#'         ylim = range(f_hat))
#' matpoints(t[c(1, length(t))], f_hat[c(1, length(t)), ],
#'           col = rainbow(ncol(f_hat)), pch = 1)
#'
#' # Using predict (and internally Predict.matrix.mdspline.smooth) to extract
#' # estimated partial effects
#' pred_terms <- predict(m_mixed, type = "terms")
#' }
#'
#' @export
Predict.matrix.mdspline.smooth <- function (object, data) {
  m <- object$m[1] + 2 # object$m[1] + 1
  x <- data[[object$term]]
  t_discrete <- object$values_discrete

  if (object$type != "discrete") {
    ll <- object$xt$domain_continuous[1] # object$knots[m]
    ul <- object$xt$domain_continuous[2] # object$knots[length(object$knots) - m + 1]
  }
  # m <- m + 1
  if (object$type == "continuous") {
    cont_positions <- seq_along(x)
    x_cont <- x
    ind <- list(cont = (x_cont <= ul & x_cont >= ll))
    design_discrete <- NULL
  }
  if (object$type == "mixed") {
    cont_positions <- which(!(x %in% t_discrete[-length(t_discrete)])) # last knot is artificial, corresponding to continuous component
    x_cont <- x[cont_positions]
    x_discrete <- x
    x_discrete[cont_positions] <- max(t_discrete)
    ind <- list(cont = (x_cont <= ul & x_cont >= ll), discrete = (x %in% t_discrete[-length(t_discrete)]))
  }
  if (object$type == "discrete") {
    x_discrete <- x
    ind <- list(discrete = (x %in% t_discrete))
    design_cont <- NULL
  }

  if (is.null(object$deriv)) {
    object$deriv <- 0
  } else if (object$deriv != 0) {
    warning("The mixed density smoother is not intended to be used for derivatives. Reasonable behavior is only guaranteed for deriv = 0.")
  }
  if (sum(sapply(ind, sum)) == length(x)) {
    if (object$type != "discrete") {
      k <- object$knots
      Z_cont <- object$Z$cont
      design_cont <- matrix(0, nrow = length(x), ncol = ncol(Z_cont))
      design_cont[cont_positions, ] <- splines::splineDesign(knots = k, x_cont, ord = m,
                                                             derivs = object$deriv) %*% Z_cont
    }
    if (object$type != "continuous") {
      k_discrete <- object$knots_discrete
      Z_discrete <- object$Z$discrete
      design_discrete <- splines::splineDesign(k_discrete, x_discrete, 1)  %*% Z_discrete
    }
    X <- cbind(design_cont, design_discrete)
  } else {
    stop("Supplied data is not in support of underlying density!")
  }
  if (object$type != "discrete") {
    if (object$mono == 0){
      return(X)
    } else {
      stop("SCOP splines are not supported yet!")
    }
  } else {
    return(X)
  }
}
