#' Clr and inverse clr transformation
#'
#' \code{clr} computes the (inverse) clr transformation of a vector \code{f}
#' of (clr transformed) density evaluations with respect to integration
#' weights \code{w} corresponding to a Bayes Hilbert space
#' \eqn{B^2(\mu) = B^2(\mathcal{Y}, \mathcal{A}, \mu)}{B^2(\mu) = B^2(Y, A, \mu)}.
#' This is a copy of \code{FDboost::clr} (with slightly adjusted help page).
#'
#' @param f A vector containing the function values (evaluated on a grid) of the
#' function \eqn{f} to transform. If \code{inverse = FALSE}, \code{f} must be a density,
#' i.e., all entries must be positive and usually \code{f} integrates to one.
#' If \code{inverse = TRUE}, \code{f} should integrate to zero, see Details.
#' For \code{\link{densreg}} objects \code{dr}, one column of \code{dr$f_hat} or
#' \code{dr$f_hat_clr} (for inverse clr) are appropriate choices for \code{f}
#' (or use \code{apply(dr$f_hat, MARGIN = 2, FUN = clr, w = dr$domain_data$Delta)}
#' to clr transform all columns of \code{dr$f_hat} at once).
#' @param w A vector of length one or of the same length as \code{f} containing
#' positive integration weights. If \code{w} has length one, this
#' weight is used for all function values. The integral of \eqn{f} is approximated
#' via \eqn{\int_{\mathcal{Y}} f \, \mathrm{d}\mu \approx
#' \sum_{j=1}^m}{\sum_{j=1}^m} \code{w}\eqn{_j} \code{f}\eqn{_j},
#' where \eqn{m} equals the length of \code{f}.
#' For \code{\link{densreg}} objects \code{dr}, the vector \code{dr$domain_data$Delta}
#' corresponds the appropriate weights to be used for \code{w}, when applying
#' \code{clr} to columns of \code{dr$f_hat} or \code{dr$f_hat_clr}
#' (for inverse clr).
#' @param inverse Logical, indicating if clr (\code{inverse = FALSE}; default)
#' or inverse clr transformation (\code{inverse = TRUE}) shall be computed.
#'
#' @details The clr transformation maps a density \eqn{f} from \eqn{B^2(\mu)} to
#' \eqn{L^2_0(\mu) := \{ f \in L^2(\mu) ~|~ \int_{\mathcal{Y}} f \, \mathrm{d}\mu = 0\}}{L^2_0(\mu) := {f \in L^2(\mu) | \int_Y f d\mu = 0}}
#' via
#' \deqn{\mathrm{clr}(f) := \log f - \frac{1}{\mu (\mathcal{Y})} \int_{\mathcal{Y}} \log f \, \mathrm{d}\mu.}{clr(f) := log f - 1/\mu(Y) * \int_Y log f d\mu.}
#' The inverse clr transformation maps a function \eqn{f} from
#' \eqn{L^2_0(\mu)} to \eqn{B^2(\mu)} via
#' \deqn{\mathrm{clr}^{-1}(f) := \frac{\exp f}{\int_{\mathcal{Y}} \exp f \, \mathrm{d}\mu}.}{clr^{-1}(f) := (exp f) / (\int_Y \exp f d\mu).}
#' Note that in contrast to Maier et al. (2025b), this definition of the inverse
#' clr transformation includes normalization, yielding the respective probability
#' density function (representative of the equivalence class of proportional
#' functions in \eqn{B^2(\mu)}).
#'
#' The (inverse) clr transformation depends not only on \eqn{f}, but also on the
#' underlying measure space \eqn{\left( \mathcal{Y}, \mathcal{A}, \mu\right)}{(Y, A, \mu)},
#' which determines the integral. In \code{clr} this is specified via the
#' integration weights \code{w}. E.g., for a discrete set \eqn{\mathcal{Y}}{Y}
#' with \eqn{\mathcal{A} = \mathcal{P}(\mathcal{Y})}{A = P(Y)} the power set of
#' \eqn{\mathcal{Y}}{Y} and \eqn{\mu = \sum_{t \in Y} \delta_t} the sum of dirac
#' measures at \eqn{t \in \mathcal{Y}}{t \in Y}, the default \code{w = 1} is
#' the correct choice. In this case, integrals are indeed computed exactly, not
#' only approximately.
#' For an interval \eqn{\mathcal{Y} = [a, b]}{Y = [a, b]}
#' with \eqn{\mathcal{A} = \mathcal{B}}{A = B} the Borel \eqn{\sigma}-algebra
#' restricted to \eqn{\mathcal{Y}}{Y} and \eqn{\mu = \lambda} the Lebesgue measure,
#' the choice of \code{w} depends on the grid on which the function was evaluated:
#' \code{w}\eqn{_j} must correspond to the length of the subinterval of \eqn{[a, b]}, which
#' \code{f}\eqn{_j} represents.
#' E.g., for a grid with equidistant distance \eqn{d}, where the boundary grid
#' values are \eqn{a + \frac{d}{2}}{a + d/2} and \eqn{b - \frac{d}{2}}{b - d/2}
#' (i.e., the grid points are centers of intervals of size \eqn{d}),
#' equal weights \eqn{d} should be chosen for \code{w}.
#'
#' @return A vector of the same length as \code{f} containing the (inverse) clr
#' transformation of \code{f}.
#'
#' @author Eva-Maria Maier
#'
#' @references
#' \code{FDboost::clr}
#'
#' Maier, E.-M., Fottner, A., Stoecker, A., Okhrin, Y., & Greven, S. (2025b):
#' Conditional density regression for individual-level data.
#' arXiv preprint arXiv:XXXX.XXXXX.
#'
#' @examples
#' ### Continuous case (Y = [0, 1] with Lebesgue measure):
#' # evaluate density of a Beta distribution on an equidistant grid
#' g <- seq(from = 0.005, to = 0.995, by = 0.01)
#' f <- dbeta(g, 2, 5)
#' # compute clr transformation with distance of two grid points as integration weight
#' f_clr <- clr(f, w = 0.01)
#' # visualize result
#' plot(g, f_clr , type = "l")
#' abline(h = 0, col = "grey")
#' # compute inverse clr transformation (w as above)
#' f_clr_inv <- clr(f_clr, w = 0.01, inverse = TRUE)
#' # visualize result
#' plot(g, f, type = "l")
#' lines(g, f_clr_inv, lty = 2, col = "red")
#'
#' @export
clr <- function(f, w = 1, inverse = FALSE) {
  stopifnot("inverse must be TRUE or FALSE." = inverse %in% c(TRUE, FALSE))
  stopifnot("f must be numeric." = is.numeric(f))
  stopifnot("f contains missing values." = !(anyNA(f)))
  n <- length(f)
  if (length(w) == 1) {
    w <- rep(w, n)
  }
  stopifnot("w must be contain positive weights of length 1 or the same length as f." =
              length(w) == n & is.numeric(w) & all(w > 0))
  int_f <- sum(f * w)
  if (!inverse) {
    stopifnot("As a density, f must be positive." = all(f > 0))
    if (!isTRUE(all.equal(int_f, 1, tolerance = 0.01))) {
      warning(paste0("f is not a probability density with respect to w. Its integral is ",
                     round(int_f, 2), "."))
    }
    return(log(f) - 1 / sum(w) * sum(log(f) * w))
  } else {
    if (!isTRUE(all.equal(int_f, 0, tolerance = 0.01))) {
      warning(paste0("f does not integrate to zero with respect to w. Its integral is ",
                     round(int_f, 2), "."))
    }
    return(exp(f) / sum(exp(f) * w))
  }
}
