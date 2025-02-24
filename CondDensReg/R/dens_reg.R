#' Conditional density regression for individual-level data
#'
#' This function implements the approach by Maier et al. (2025b) for fitting
#' structured additive regression models with densities in a mixed (continuous/
#' discrete) Bayes Hilbert space
#' \eqn{B^2(\mu) = B^2(\mathcal{Y}, \mathcal{A}, \mu)}{B^2(\mu) = B^2(Y, A, \mu)}
#' as response given scalar covariates \eqn{x_i},
#' based on observations \eqn{y_i} from the conditional distributions given
#' \eqn{x_i}, via a (penalized) maximum likelihood approach. The (penalized)
#' log-likelihood function is approximated via the (penalized) log-likelihood of
#' an appropriate poisson regression model, which (after constructing the count
#' data appropriately) is then fitted using \code{mgcv}'s \code{\link[mgcv]{gam}}
#' with a new smooth for mixed densities in \eqn{L^2_0(\mu)}.
#' We briefly summarize the approach below in the details.
#' % to enable proper description of the function.
#' Please see Maier et al. (2025b) for comprehensive description.
#'
#' The function \code{dens_reg} estimates the densities \eqn{f_{x_i}} of
#' conditional distributions of random variables \eqn{Y_i | x_i} from independent
#' observations \eqn{(y_i, x_i)}, where \eqn{y_i} are realizations of \eqn{Y_i | x_i},
#' and \eqn{x_i} is a vector of covariate observations, \eqn{i = 1, ..., N}. The
#' densities are considered as elements of a mixed Bayes Hilbert space
#' \eqn{B^2(\mu) = B^2(\mathcal{Y}, \mathcal{A}, \mu)}{B^2(\mu) = B^2(Y, A, \mu)}
#' (including continuous and discrete ones as special cases). The subset of the
#' domain \eqn{\mathcal{Y}}{Y} corresponding to the continuous part of the densities
#' is denoted with \eqn{I}, the one corresponding to the discrete part with
#' \eqn{\mathcal{D}}{D}. The densities are modeled via a structured additive
#' regression model
#' \deqn{f_{x_i} = \bigoplus_{j = 1}^J h_j (x_i)}
#' with partial effects \eqn{h_j (x_i) \in B^2(\mu)} depending on no, one, or
#' several covariates \eqn{x_i}. Each partial effect is represented using a
#' tensor product basis, consisting of an appropriate vector of basis functions
#' \eqn{b_{\mathcal{Y}}{b_Y}}{b_Y} in \eqn{B^2(\mu)} over the domain of
#' \eqn{B^2(\mu)}, and a vector of basis functions \eqn{b_{\mathcal{X}, j}}{b_{X, j}}
#' over the respective covariates.
#' To obtain basis functions \eqn{b_{\mathcal{Y}}}{b_Y} consider the orthogonal
#' decomposition of the mixed Bayes Hilbert space \eqn{B^2(\mu)} into a discrete
#' Bayes Hilbert space \eqn{B^2(\delta^\bullet)} and a continuous Bayes Hilbert
#' space \eqn{B^2(\lambda)} developed in Section 3.4 of Maier et al. (2025a).
#' Note that for the domain of the discrete Bayes Hilbert space, the discrete part
#' \eqn{\mathcal{D}}{D} has to be extended by an additional arbitrary discrete value
#' (for the discrete summary of the continuous part).
#' We construct basis functions in both these spaces by transforming appropriate
#' basis functions in the corresponding (unconstrained) \eqn{L^2} spaces (more
#' precisely, indicator functions with optional difference penalty for the discrete
#' and a P-spline basis for the continuous Bayes Hilbert space), to the respective
#' \eqn{L^2_0} spaces, i.e., the subspaces of \eqn{L^2} containing only functions
#' that integrate to zero.
#' % \eqn{L^2_0(\mu) : =  \{ f \in L^2(\mu) ~|~ \int_{\mathcal{T}} f \, \mathrm{d}\mu = 0\}}{L^2_0(\mu) : =  {f \in L^2(\mu) | \int_T f d\mu = 0}}
#' See appendix D of Maier et al. (2025a) for details on this transformation.
#' Applying the inverse centered log-ratio (clr) transformation (which is an
#' isometric isomorphism, allowing to consider densities in some \eqn{B^2} space
#' equivalently in the respective \eqn{L^2_0} space) yields one basis in
#' \eqn{B^2(\delta^\bullet)} and one in \eqn{B^2(\lambda)}. The basis in the actually
#' considered mixed Bayes Hilbert space \eqn{B^2(\mu)} is then obtained by
#' applying the respective embeddings to these two bases (compare Section 3.4 of
#' Maier et al. (2025a) and Section 2.2 of Maier et al. (2025b)).
#' The choice of \eqn{b_{\mathcal{X}, j}}{b_{X, j}} determines
#' the type of the partial effect, e.g., linear or smooth for a continuous covariate.
#' For smooth effects, the same marginal bases as in \code{\link[mgcv]{gam}}
#' can be used. In particular, penalization is also possible.
#' The resulting log-likelihood is then approximated by the log-likelihood of a
#' corresponding multinomial model, which can equivalently be estimated by a
#' Poisson model. The data for these models is obtained from the original
#' observations \eqn{y_1, ..., y_N} by combining all observations of the same
#' conditional distribution into a vector of counts via a histogram on the
#' continuous part of the domain of the densities and counts on the discrete
#' part of the domain. For details, see Maier et al. (2025b).
#' \code{dens_reg} constructs the respective count data from the individual
#' observations and estimates the corresponding Poisson model. Furthermore, the
#' resulting densities on pdf- and clr-level as well as the estimated partial
#' effects on both levels are calculated.
#'
#' @encoding UTF-8
#'
#' @import mgcv
#' @import ggplot2
#' @import data.table
#' @import tidyr
#' @importFrom dplyr "%>%" arrange mutate
#' @importFrom FDboost clr
#' @importFrom Rdpack reprompt
#'
#' @param dta Data set of type \code{\link[base]{data.frame}} or
#' \code{\link[data.table]{data.table}} containing the observations \eqn{(y_i, x_i)}
#' of response and covariates as well as optional sample weights (compare
#' \code{sample_weights}) for each observation in the rows (\eqn{i = 1, ..., N}).
#' @param y Variable in \code{dta} containing the response observations
#' \eqn{y_i}. Either the variable name can be given as string or the column
#' position of the variable in \code{dta} as integer. If missing (\code{NULL}),
#' if there is a unique column of \code{dta} not specified in \code{sample_weights},
#' \code{counts}, and \code{weighted_counts} and not used as covariate during
#' effect specification (see below), this unique column is used.
#' @param sample_weights (Optional) variable in \code{dta} which contains a sample
#' weight for each observation. Either the variable name can be given as string
#' or the column position of the variable in \code{dta} as integer. If missing
#' (\code{NULL}), no sample weights are included per default (i.e., all
#' observations have the same weight 1).
#' @param counts (Optional) variable in \code{dta} which contains a count for
#' each observation (for cases, where the available data contains counts instead
#' of individual observations, which is common in particular for discrete data).
#' Either the variable name can be given as string or the column position of the
#' variable in \code{dta} as integer. If \code{bin_number} and \code{bin_width}
#' are both \code{NULL}, midpoints of unique observations \eqn{y_i} are used as
#' boundaries of histogram bins (if \eqn{I} is not empty, i.e., if there is a
#' continuous component), which are then used to compute the bin widths (which are
#' later required es offset in the regression model).
#' If argument \code{counts} is missing (\code{NULL}; default), counts are
#' constructed from individual observations (which is equivalent to each observation
#' counted once). Note that if \code{counts} and \code{sample_weights} are both
#' not \code{NULL}, the latter is ignored (also indicated via a warning message).
#' Please use \code{weighted_counts} (additionally to \code{counts}) to include
#' possible weighted data.
#' @param weighted_counts (Optional) variable in \code{dta} which contains a
#' weighted count for each observation (compare Appendix D of Maier et al. (2025b)).
#' In this case, also absolute counts have to be specified via \code{counts}.
#' Otherwise, \code{weighted_counts} is ignored. Either the variable name can be
#' given as string or the column position of the variable in \code{dta} as integer.
#' If missing (\code{NULL}), counts are constructed from individual observations
#' (which is equivalent to all observations counted once)
#' % @param already_formatted A logical indicating if the data in \code{dta} is
#' % already formatted as count data. If \code{already_formatted = TRUE}, the data
#' % have to have a column named \code{"counts"}, an additional column with the
#' % name \code{"weighted_counts"} is optional.
#' % % The relevant variables used for further aggregation are submitted via \code{var_vec}, the relevant column of the observed density via \code{y}.
#' % The bin width is computed automatically based on the observed continuous
#' % values unless an integer or a vector is submitted by the user via
#' % \code{bin_width} and/or \code{bin_number}. In these cases, the binning is
#' % based on the given bin number or width.
#' .
#' @param values_discrete Vector of values in \eqn{\mathcal{D}}{D} (the subset of
#' the domain corresponding to the discrete part of the densities). Defaults to
#' missing (\code{NULL}) in which case it is set to \code{c(0, 1)}. If set to
#' \code{FALSE}, the discrete component is considered to be empty, i.e., the
#' Lebesgue measure is used as reference measure (continuous special case).
#' @param weights_discrete Vector of weights for the Dirac measures corresponding
#' to \code{values_discrete}. If missing (\code{NULL}) it is set to 1 in all
#' components as default. Can be a scalar for equal weights for all discrete values
#' or a vector with specific weights for each corresponding discrete value.
#' @param domain_continuous An interval (i.e., a vector of length 2) specifying
#' \eqn{I} (the subset of the domain corresponding to the continuous part of the
#' densities). If missing (\code{NULL}) it is set to \code{c(0, 1)} as default. If
#' set to \code{FALSE}, the continuous component is considered to be empty, i.e.,
#' a weighted sum of dirac measures is used as reference measure (discrete special
#' case).
#' @param bin_number Number of equidistant histogram bins partitioning \eqn{I}.
#' Alternative to \code{bin_width}. If neither parameter is specified,
#' \code{bin_number = 100} is used as default. If \code{bin_number} and \code{bin_width} are
#' both given and the two values are not compatible, an error is returned.
#' @param bin_width Width of histogram bins partitioning \eqn{I}. Can be one
#' scalar value (specifying an equidistant bin width), or a vector containing the
#' width of each bin. The combined length of the specified bins must match the
#' length of the continuous part of the domain. Alternative to \code{bin_number}.
#' If \code{bin_number} and \code{bin_width} are both given and the two values
#' are not compatible, an error is returned.
#' @param m_continuous Vector of two integers specifying the order of the B-spline
#' basis over \eqn{I} and the order of the difference penalty (like the argument
#' \code{m} for P-Spline smooth terms \code{bs = "ps"} in \code{\link[mgcv]{ti}},
#' etc.) for basis functions in \eqn{L^2(\lambda)}, before transformation to
#' \eqn{L^2_0(\lambda)} (compare details). If missing it is set to cubic splines
#' with second order difference penalty, i.e., \code{m_continuous = c(2, 2)},
#' as default.
#' @param k_continuous Integer specifying the number of B-spline basis functions
#' in \eqn{L^2(\lambda)} (like the argument \code{k} for P-Spline smooth terms
#' \code{bs = "ps"} in \code{\link[mgcv]{ti}}, etc.), before transformation
#' to \eqn{L^2_0(\lambda)} (compare details). Note that the transformation reduces
#' the basis number by 1. The basis \eqn{b_{\mathcal{Y}}{b_Y}}{b_Y} in the mixed
#' Bayes Hilbert space \eqn{B^2(\mu)} will thus have
#' \code{k_continuous - 1 + length(values_discrete)} elements. See also
#' \code{\link[mgcv]{choose.k}} for more information on choosing this parameter.
#' If missing (\code{NULL}) it is set to 10.
#' @param sp_y Integer or vector specifying the smoothing parameter for the marginal
#' penalty matrix for \eqn{b_{\mathcal{Y}}{b_Y}}{b_Y} (anisotropic penalty; like
#' the argument \code{sp} in \code{\link[mgcv]{ti}}, etc.). If a vector is submitted,
#' its length must be the number of partial effects (including the intercept), with
#' the \eqn{j}-th entry specifying the smoothing parameter for the \eqn{j}-th
#' partial effect. The order of parameters within this vector corresponds to the
#' intercept followed by the partial effects in the order as specified below
#' (\code{group_specific_intercepts}, \code{linear_effects}, \code{smooth_effects},
#' \code{varying_coefficients}, and \code{smooth_interactions}). Upenalized
#' estimation is accomplished by setting \code{sp_y} to zero. If missing (\code{NULL})
#' or a negative value is supplied, the parameter will be estimated by
#' \code{\link[mgcv]{gam}} via the estimation method specified in \code{method}
#' (with REML-estimation per default). Any positive or zero value is treated as
#' fixed (see \code{\link[mgcv]{gam}}).
#' @param method String characterizing the smoothing parameter estimation method
#' as in \code{\link[mgcv]{gam}}. Defaults to \code{"REML"} (REML-estimation).
#' @param penalty_discrete Integer or \code{NULL} giving the order of differences
#' to be used for the penalty of the discrete component, with 0 corresponding to
#' the identity matrix as penalty matrix (analogously to m[2] for the continuous
#' component); Note that the order of differences must be smaller than the number
#' of values in the discrete component, i.e., the length of \code{values_discrete}.
#' % In the discrete case, please set the corresponding argument of sp in ti() to 0 to estimate an unpenalized model.
#' If missing (\code{NULL}), the discrete component is estimated unpenalized (in
#' the mixed case, by setting the discrete component of the penalty matrix to zero
#' - with the continuous component non-zero -, in the discrete case by setting
#' the smoothing parameter to zero (with a diagonal penalty matrix) since due to
#' technical reasons it is not possible to set the whole penalty matrix to zero).
#' @param group_specific_intercepts Vector of the form \code{c("x_a", "x_b", ...)}
#' of names of categorical covariates, i.e., \code{\link[base]{factor}} variables
#' contained in \code{dta} for which group-specific intercepts \eqn{\beta_{x_a},
#' \beta_{x_b}, ...} (one per category of the respective covariate) shall be
#' estimated. For \code{\link[base]{ordered}} factors, the first level is used as reference category
#' (see \code{\link[mgcv]{gam.models}} for more details).
#' If missing (\code{NULL}), no group-specific intercept is included.
#' @param linear_effects Vector of the form \code{c("x_a", "x_b", ...)} of names
#' of numeric covariates contained in \code{dta} for which linear effects
#' \eqn{x_a \odot \beta_{x_a}, x_b \odot \beta_{x_b}, ...} % {x_a * \beta_{x_a}, x_b * \beta_{x_b}, ...}
#' shall be estimated. If missing (\code{NULL}), no linear effect is included.
#' @param smooth_effects List of (named) lists of the form
#' \code{list(list(cov = "x_a", bs = bs_a, m = m_a, k = k_a, mc = mc_a, by = "x_b"), ...)}.
#' If the lists are unnamed, the names are assigned in the order of the given
#' elements. The list is filled with \code{NULL} if it contains less than 6 elements.
#' Each list is adding one (group-specific) smooth effect of the form \eqn{g_{x_b}(x_a)} to the
#' model with:
#' \itemize{
#' \item \code{cov}: Name of a numeric covariate contained in \code{dta}.
#' \item \code{bs}: Character string specifying the type for the marginal
#' basis in covariate direction. See \code{\link[mgcv]{smooth.terms}} for details and
#' full list. If not specified or \code{NULL}, a P-spline basis \code{"ps"}
#' is used.
#' \item \code{m}: Vector of two integers or \code{NULL} giving the order of the
#' marginal spline basis \eqn{b_j} in covariate direction and the order of its
#' penalty (as in \code{\link[mgcv]{ti}}). If not specified or \code{NULL},
#' \code{c(2, 2)} is used.
#' \item \code{k}: Integer or \code{NULL} giving the dimension of the marginal
#' basis \eqn{b_j} (as in \code{\link[mgcv]{ti}}). See \code{\link[mgcv]{choose.k}}
#' for more information. If not specified or \code{NULL}, \code{10} is used.
#' \item \code{mc}: Logical indicating if the marginal in covariate direction
#' should have centering constraints applied. By default all marginals are
#' constrained, i.e., \code{mc = TRUE}.
#' \item \code{by}: Optional, name of a categorical covariate, i.e.,
#' \code{\link[base]{factor}} variable contained in \code{dta} specifying whether
#' the smooth effect should be modeled specifically for each level of the
#' \code{by}-covariate (group-specific). For \code{\link[base]{ordered}} factors,
#' the first level is used as reference category (see \code{\link[mgcv]{gam.models}}
#' for more details). If missing or \code{NULL}, the smooth effect is not depending
#' on the level of an additional covariate.
#' }
#' If missing (\code{NULL}), no (group-specific) smooth effect is included.
#' @param varying_coefficients  List of lists of the form
#' \code{list(list(cov = "x_a", by = "x_b", bs = bs_a, m = m_a, k = k_a, mc = mc_a), ...)}.
#' If the lists are unnamed, the names are assigned in the order of the given
#' elements. The list is filled with \code{NULL} if it contains less than 6 elements.
#' Each list is adding one varying coefficient of the form \eqn{x_b \odot g(x_a)}
#' to the model with:
#' \itemize{
#' \item \code{cov}: Name of a numeric covariate contained in \code{dta}.
#' \item \code{by}: Name of a numeric covariate contained in \code{dta}.
#' \item \code{bs}: Character string specifying the type for the marginal
#' basis in covariate direction. See \code{\link[mgcv]{smooth.terms}} for details and
#' full list. If not specified or \code{NULL}, a P-spline basis \code{"ps"}
#' is used.
#' \item \code{m}: Vector of two integers or \code{NULL} giving the order of the
#' marginal spline basis \eqn{b_j} in covariate direction and the order of its
#' penalty (as in \code{\link[mgcv]{ti}}). If not specified or \code{NULL},
#' \code{c(2, 2)} is used
#' \item \code{k}: Integer or \code{NULL} giving the dimension of the marginal
#' basis \eqn{b_j} (as in \code{\link[mgcv]{ti}}). See \code{\link[mgcv]{choose.k}}
#' for more information. If not specified or \code{NULL}, \code{10} is used.
#' \item \code{mc}: Logical indicating if the marginal in the direction of the
#' first covariate should have centering constraints applied. By default all
#' marginals are constrained, i.e., \code{mc = TRUE}.
#' }
#' If missing (\code{NULL}), no varying coeffecient is included.
#' @param smooth_interactions  List of lists of the form
#' \code{list(list(covs = c("x_a", "x_b", ...), bs = c(bs_a, bs_b, ...),
#' m = list(m_a, m_b, ...), k = c(k_a, k_b, ...), mc = c(mc_a, mc_b, ...), by = "x_by"),
#' list(...), ...)}. If the lists are unnamed, the names are assigned in the order
#' of the given elements. The list is filled with vectors/lists of \code{NULL}
#' if it contains less than 6 elements. Each list is specifying a (group-specific)
#' smooth interaction effect between at least two continuous covariates of the
#' form \eqn{g_{x_by}(x_a, x_b, ...)} in the model with:
#' \itemize{
#' \item \code{covs}: Vector of names of numeric covariates contained in \code{dta}.
#' \item \code{bs}: Vector of character strings specifying the types for each
#' marginal basis of each covariate. See \code{\link[mgcv]{smooth.terms}} for
#' details and full list. If not specified or \code{NULL}, a P-spline basis \code{"ps"}
#' is used.
#' \item\code{m}: List containing vectors of two integers or \code{NULL} giving
#' the order of the marginal spline basis in direction of the respective covariate
#' in \code{covs} for the smooth effect and the order of its penalty (as in
#' \code{\link[mgcv]{ti}}). If not specified or \code{NULL}, a list of
#' \code{c(2, 2)} is used.
#' \item \code{k}: Vector of integers or \code{NULL} giving the
#' dimension of the marginal basis in direction of the respective covariate
#' in \code{covs} for the smooth effect. See \code{\link[mgcv]{choose.k}} for
#' more information. If not specified or \code{NULL}, a vector specifying each
#' parameter as \code{10} is used.
#' \item \code{mc}: Logical vector indicating if the marginals in covariate direction
#' should have centering constraints applied. By default all marginals are
#' constrained, i.e., \code{TRUE}.
#' \item \code{by}: Optional, name of a categorical covariate, i.e.,
#' \code{\link[base]{factor}} variable contained in \code{dta} specifying whether
#' the smooth interaction effect should be modeled specifically for each level of the
#' \code{by}-covariate (group-specific). For \code{\link[base]{ordered}} factors,
#' the first level is used as reference category (see \code{\link[mgcv]{gam.models}}
#' for more details). If missing or \code{NULL}, the smooth interaction effect
#' is not depending on the level of an additional covariate.
#' }
#' If missing (\code{NULL}), no (group-specific) smooth interaction is included.
#' @param effects Indicates if estimated partial effects should be returned
#' (\code{TRUE}; default) or not (\code{FALSE}).
#' @param ...  further arguments for passing on to \code{\link[mgcv]{gam}}.
#'
#'
#' @return The function returns an object of the class \code{dens_reg_obj}, which
#' is a \code{\link[base]{list}} with elements:
#' \itemize{
#' \item \code{count_data}: \code{\link[data.table]{data.table}}-object containing
#' the count data obtained by using \code{\link{preprocess}} for the given data
#' \code{dta} and covariates, which is also an object of the sub-class
#' \code{histogram_count_data}. See \code{?preprocess} for more information.
#' \item \code{model}: \code{\link[mgcv]{gam}}-object of the estimated model.
#' \item \code{model_matrix}: Model matrix of the model.
#' \item \code{theta_hat}: Estimated coeffecient vector \eqn{\hat{\theta}}.
#' \item \code{f_hat_clr}: Estimated conditional densities on clr-level
#' \eqn{clr(\hat f)} for every covariate combination in the order of the respective
#' \code{group_id} in \code{count_data}.
#' \item \code{f_hat}: Estimated conditional densities on clr-level \eqn{clr(\hat{f})}
#' for every covariate combination in the order of the respective group_id in
#' \code{count_data}.
#' \item \code{effects}: Only contained, if \code{effects = TRUE}; List of lists.
#' Each list gives one estimated partial effect. Both clr- and density-level are
#' included.
#' \item \code{params}: List of \code{domain_continuous, values_discrete} and
#' \code{bin_number} as given to the function.
#' \item \code{predicted_effects}: List of lists (\code{group_specific_intercepts,
#' smooth effects, linear_effects, varying_coefficient, smooth_interactions})
#' collecting the specification of all partial effects as given to the function
#' in the respective parameters.
#' \item \code{ID_covCombi}: Data frame which gives an overview over the assignment
#' of the unique covariate combinations to the group IDs.
#' }
#' Note that \code{plot}- and \code{predict}-methods for objects of class
#' \code{dens_reg_obj} are available via \code{CondDensReg:::plot.dens_reg_obj}
#' and \code{CondDensReg:::predict.dens_reg_obj}, however, they are not exported,
#' since they are not tested/documented appropriately, yet.
#'
#' @author Lea Runge, Eva-Maria Maier
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
#' ### Note that the following simulated data are only to illustrate
#' ### function usage and do not possess significant covariate effects
#'
#' # for further information on the parameters of the preprocessing step see ?preprocess
#'
#' # create data for the mixed case
#' set.seed(101)
#'
#' dta <- data.frame(obs_density = sample(0:2, 150, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
#'                   covariate1 = sample(c("a", "b", "c"), 150, replace = TRUE),
#'                   covariate2 = sample(c("c", "d"), 150, replace = TRUE),
#'                   covariate3 = rep(rnorm(n = 15), 10),
#'                   covariate4 = rep(rnorm(n = 10), 15),
#'                   covariate5 = rep(rnorm(n = 10), 15), sample_weights = runif (150, 0, 2))
#' dta[which(dta$obs_density == 2), ]$obs_density <- rbeta(length(which(dta$obs_density==  2)),
#'                                                        shape1 = 3, shape2 = 3)
#' dta$covariate1 <- ordered(dta$covariate1)
#' dta$covariate2 <- ordered(dta$covariate2)
#'
#' # create discrete data
#'
#' dta_dis <- data.frame(obs_density = sample(0:2, 150, replace = TRUE, prob = c(0.25, 0.45, 0.3)),
#'                       covariate1 = sample(c("a", "b", "c"), 150, replace = TRUE),
#'                       covariate2 = sample(c("c", "d"), 150, replace = TRUE),
#'                       covariate3 = rep(rnorm(n = 15), 10),
#'                       covariate4 = rep(rnorm(n = 10), 15),
#'                       covariate5 = rep(rnorm(n = 10), 15),
#'                       sample_weights = runif (150, 0, 2))
#' dta_dis$covariate1 <- ordered(dta_dis$covariate1)
#' dta_dis$covariate2 <- ordered(dta_dis$covariate2)
#'
#' # examples for different partial effects
#'
#' ## group specific intercepts
#' group_specific_intercepts <- c("covariate1", "covariate2")
#' ## linear effects
#' linear_effects <- c("covariate4")
#' ## smooth effects
#' smooth_effects <- list(list(cov = "covariate3", bs = "ps", m = c(2, 2), k = 4),
#'                        list(cov = "covariate3", bs = "ps", m = c(2, 2), k = 4,
#'                             mc = FALSE, by = "covariate1"))
#' ## varying coefficient
#' varying_coef <- list(list(cov = "covariate3", by = "covariate4", bs = "ps", m = c(2, 2) , k = 4))
#' ## smooth interaction
#' smooth_inter <- list(list(covs = c("covariate3", "covariate4", "covariate5"),
#'                           bs = c("ps", "ps", "ps"),
#'                           m = list( c(2, 2), c(2, 2), c(2, 2)), k = c( 4, 4, 5),
#'                           mc = c(TRUE, FALSE, TRUE), by = NULL))
#'
#' # fit models (warning: calculation may take a few minutes)
#'
#' ## fit model for the mixed case with group specific intercepts and linear effects
#' ### use fixed smoothing parameters in density direction and calculate also the partial effects
#'
#' m_mixed <- dens_reg(dta = dta, y = 1, m_continuous = c(2, 2),
#'   k_continuous = 4, group_specific_intercepts = group_specific_intercepts,
#'   linear_effects = linear_effects, effects = TRUE, sp_y = c(1, 3, 5, 0.5))
#'
#' ## fit model for the discrete case with smooth effects and smooth interaction
#' ### do not calculate effects
#'
#' m_dis <- dens_reg(
#'   dta = dta_dis, y = 1, values_discrete = c(0, 1, 2),
#'   weights_discrete = c(1, 1, 1), domain_continuous = FALSE, m_continuous = c(2, 2),
#'   k_continuous = 4, group_specific_intercepts = group_specific_intercepts,
#'   smooth_effects = smooth_effects, smooth_interactions = smooth_inter, effects = FALSE)
#'
#' # fit model for the continuous case with a functional varying coeffecient
#'
#' m_cont <- dens_reg(dta = dta[which(!(dta$obs_density %in% c(0, 1))), ],
#'   y = 1, values_discrete = FALSE, m_continuous = c(2, 2),
#'   k_continuous = 12, varying_coefficients = varying_coef, effects = TRUE)
#' }
#'
#' @export


dens_reg <- function(dta, y = NULL, sample_weights = NULL, counts = NULL,
                     weighted_counts = NULL, values_discrete = c(0, 1),
                     weights_discrete = 1, domain_continuous = c(0, 1),
                     bin_number = NULL, bin_width = NULL, m_continuous = c(2, 2),
                     k_continuous = 10, sp_y = NULL, method = "REML",
                     penalty_discrete = NULL, group_specific_intercepts = NULL,
                     # list: c("var_name1", ...)
                     ##  ti(y, bs = "md", m = m_continuous, k = k_continuous, mc = FALSE, np = FALSE, by = var_name1)
                     ## beta_k
                     linear_effects = NULL,
                     # list: c("var_name1", ...)
                     ## ti(y, bs = "md", m = m_continuous, k = k_continuous, mc = FALSE, np = FALSE, by = var_name1)
                     ## x*beta
                     # linear_interaction,
                     ### not working      # list of lists: list(list("var_nameA_1", "var_nameB_1"), ...)
                     ## ti(y, bs = "md", m = m_continuous, k = k_continuous, mc = FALSE, np = FALSE)
                     ## x1*(x2*beta)
                     smooth_effects = NULL,
                     # list of lists: list(list("var_name1", "basis1", m1, k1, "by1"), list(...), ...) #if by = empty - >  by = NULL)
                     ## ti(var_name1, y, bs = c(basis1, "md"), m = list(m1, m_continuous), k = c(k1, k_continuous), mc = c(TRUE, FALSE), np = FALSE, by = by)
                     ## g(x) or g_k(x)
                     varying_coefficients = NULL,
                     # list of lists: list(list("var_nameA_1", "var_nameB_1", basis1, m1, k1), ...)
                     ## ti(var_nameA_1, y, bs = c(basis1, "md"), m = list(m1, m_continuous), k = c(k1, k_continuous), mc = c(TRUE, FALSE), np = FALSE, by = var_nameB_1)
                     ## x1*g(x2)
                     smooth_interactions = NULL,
                     # list of lists: list(list(c("var_nameA_1", "var_nameB_1", ...), c("basis1A", "basis1B", ...), c(m1A, m1B, ...), c(k1A, k1B, ...)), list(...), ...)
                     ## ti(var_nameA_1, var_nameB_1, y, bs = c(basis1A, basis1B, "md"), m = list(m1A, m1B, m_continuous), k = c(k1a, k1B, k_continuous), mc = c(TRUE, TRUE, FALSE), np = FALSE)
                     ## g(x1, x2)
                     effects = TRUE, ...) {
  if (isFALSE(values_discrete)) {
    weights_discrete <- NULL
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
  var_vec <- intersect(names(dta),
                       unlist(c(group_specific_intercepts, linear_effects,
                                smooth_effects, varying_coefficients,
                                smooth_interactions)))
  if (is.null(y)) {
    y <- setdiff(names(dta), c(var_vec, sample_weights, counts, weighted_counts))
  }
  stopifnot("y must be one variable contained in dta. Can only be NULL, if uniquely determined by sample_weights and covariates used to specify effects." = length(y)==  1)
  dta_est <- preprocess(dta = dta, var_vec = var_vec, y = y,
                        sample_weights = sample_weights, counts = counts,
                        weighted_counts = weighted_counts, bin_width = bin_width,
                        bin_number = bin_number, values_discrete = values_discrete,
                        weights_discrete = weights_discrete, domain_continuous = domain_continuous)
  cov_combi_id <- unique(subset(dta_est, TRUE,
                                names(dta_est) %in% c(var_vec, sample_weights)))

  checking_dens_reg_1(m_continuous, k_continuous, sp_y, penalty_discrete,
                      effects, dta)

  if (is.numeric(y)) {
    y <- colnames(dta)[y]
  }
  if (is.numeric(var_vec)) {
    var_vec <- colnames(dta)[var_vec]
  }
  if (!is.null(weights_discrete)) {
    if (length(weights_discrete)==  1) {
      weights_discrete <- rep(weights_discrete, length(values_discrete))
    }
  }
  if (is.null(sp_y)) {
    sp_y_vec <- -1
    sp_y_ <- "NULL"
  } else {

    if (length(sp_y) == 1) {
      sp_y_ <- sp_y
      sp_y_vec <- sp_y
    }
  }

  if (is.null(penalty_discrete)) {
    if (isFALSE(domain_continuous)) {
      fx <- "FALSE"
      penalty_discrete <- "NULL"
      sp_y_ <- 0 # in discrete case, unpenalized estimation has to be specified via zero smoothing parameter
    }  else{
      fx <- "FALSE"
      penalty_discrete <- "NULL"
    }
  } else{
    fx = "FALSE"
  }
  j <- 1
  if (length(sp_y) > 1) {
    stopifnot("If sp_y is supplied as vector, its length must equal the number of partial effect (including the intercept)." = j <= length(sp_y))
    sp_y_ <- sp_y[j]
    sp_y_vec <- sp_y[j]
  }

  f_main <- paste0("counts~ti(", y, ", bs = 'md', mc = FALSE, np = FALSE, m = list(",
                   deparse(m_continuous), "), k = ", k_continuous, ", sp = ", sp_y_,
                   ", fx = ", fx, ", xt = list(list(values_discrete =  ", deparse(values_discrete),
                   ",  domain_continuous = ", deparse(domain_continuous),
                   ", weights_discrete = ", deparse(weights_discrete),
                   ", penalty_discrete = ", penalty_discrete, ")))")
  f_intercepts <- ""
  f_flexibles <- ""
  f_linear <- ""
  f_function_var_coef <- ""
  f_flex_interact <- ""
  j <- j+1
  for (effect in group_specific_intercepts) {
    if (length(sp_y) > 1) {
      stopifnot("If sp_y is supplied as vector, its length must equal the number of partial effect (including the intercept)." = j <= length(sp_y))
      sp_y_ <- sp_y[j]
      sp_y_vec <- sp_y[j]
    }
    f_intercepts <- paste0(f_intercepts, " + ti(", y, ", bs = 'md', m = list(",
                           deparse(m_continuous), "), k = ", k_continuous,
                           ", mc = FALSE, np = FALSE, by = ", effect, ", sp = ",
                           sp_y_, ", fx = ", fx, ", xt = list(list(values_discrete = ",
                           deparse(values_discrete), ", domain_continuous = ",
                           deparse(domain_continuous), ", weights_discrete = ",
                           deparse(weights_discrete), ", penalty_discrete =  ",
                           penalty_discrete, ")))")
    j <- j+1
  }
  l <- 1
  for (effect in smooth_effects) {
    if (!is.null(names(effect))) {
      if (!all(names(effect) %in% c("cov", "bs", "m", "k", "mc", "by"))) {
        stop("Names in list specifying smooth effect are incorrect!")
      }
      if (is.null(effect$bs)) {
       effect$bs <- "ps"
      }
      if (is.null(effect$m)) {
        effect$m <- c(2, 2)
      }
      if (is.null(effect$k)) {
        effect$k <- 10
      }
      if (is.null(effect$mc)) {
        effect$mc <- TRUE
      }
      if (is.null(effect$by)) {
        effect$by <- NULL
      }
      effect <- list(cov = effect$cov, bs = effect$bs, m = effect$m, k = effect$k,
                     mc = effect$mc, by = effect$by)
      smooth_effects[[l]] <- effect
    } else {
      effect <- append(effect, rep(list(NULL), 6-length(effect)))
      names(effect) <- c("cov", "bs", "m", "k", "mc", "by")
      if (is.null(effect$bs)) {
        effect$bs <- "ps"
      }
      if (is.null(effect$m)) {
        effect$m <- c(2, 2)
      }
      if (is.null(effect$k)) {
        effect$k <- 10
      }
      if (is.null(effect$mc)) {

        effect$mc <- TRUE
      }
      if (is.null(effect$by)) {
        effect$by <- NULL
      }
      smooth_effects[[l]] <- effect
    }

    if (length(sp_y) > 1) {
      stopifnot("If sp_y is supplied as vector, its length must equal the number of partial effect (including the intercept)." = j <= length(sp_y))
      sp_y_ <- sp_y[j]
      sp_y_vec <- sp_y[j]
    }
    if (is.null(effect[5][[1]])) {
      mc <- TRUE
    }
    else{
      mc <- effect[5][[1]]
    }
    l <- l+1
    if (is.null(effect[6][[1]])) {
      f_flexibles <- paste0(f_flexibles, "+ ti(", effect[1], ", ", y, ", bs = c(",
                            '"', effect[2][[1]], '"', ", ", '"md")', ", m = list(",
                            deparse(effect[3][[1]]), ", ", deparse(m_continuous),
                            "), k = c(", effect[4], ", ", k_continuous,
                            "), mc = c(", mc, ", FALSE), np = FALSE, sp = array(c(-1, ",
                            sp_y_vec, " ))", ", fx = ", fx, ", xt = list(NULL, list(values_discrete = ",
                            deparse(values_discrete), ",  domain_continuous = ",
                            deparse(domain_continuous), ", weights_discrete = ",
                            deparse(weights_discrete), ", penalty_discrete =  ",
                            penalty_discrete, ")))")
    } else {
      f_flexibles <- paste0(f_flexibles, "+ ti(", effect[1], ", ", y, ", bs = c(",
                            '"', effect[2][[1]], '"', ", ", '"md")', ", m = list(",
                            deparse(effect[3][[1]]), ", ", deparse(m_continuous),
                            "), k  = c( ", effect[4], ", ", k_continuous, ")",
                            ", mc = c(", mc, ", FALSE), np = FALSE, by = ",
                            effect[6], ", sp = array(c(-1, ", sp_y_vec,
                            " )), fx = ", fx, ", xt = list(NULL, list(values_discrete = ",
                            deparse(values_discrete), ",  domain_continuous = ",
                            deparse(domain_continuous), ", weights_discrete = ",
                            deparse(weights_discrete), ", penalty_discrete =  ",
                            penalty_discrete, ")))")
    }
    j <- j+1
  }
  for (effect in linear_effects) {
    if (length(sp_y) > 1) {
      stopifnot("If sp_y is supplied as vector, its length must equal the number of partial effect (including the intercept)." = j <= length(sp_y))
      sp_y_ <- sp_y[j]
      sp_y_vec <- sp_y[j]
    }
    f_linear <- paste0(f_linear, "+ ti(", y, ", bs = ", '"md"', ", m = list(",
                       deparse(m_continuous), "), k = ", k_continuous,
                       ", mc = FALSE, np = FALSE, by = ", effect, ", sp = ",
                       sp_y_, ", fx = ", fx, ", xt = list(list(values_discrete = ",
                       deparse(values_discrete), ", domain_continuous = ",
                       deparse(domain_continuous), ", weights_discrete = ",
                       deparse(weights_discrete), ", penalty_discrete =  ",
                       penalty_discrete, ")))")
    j <- j+1
  }
  l <- 1
  for (effect in varying_coefficients) {
    if (length(sp_y) > 1) {
      stopifnot("If sp_y is supplied as vector, its length must equal the number of partial effect (including the intercept)." = j <= length(sp_y))
      sp_y_ <- sp_y[j]
      sp_y_vec <- sp_y[j]
    }

    if (!is.null(names(effect))) {
      if (!all(names(effect)%in%c("cov", "by", "bs", "m", "k", "mc"))) {
        stop("Names of varying coefficient incorrect!")
      }
      if (is.null(effect$bs)) {
        effect$bs <- "ps"
      }
      if (is.null(effect$m)) {
        effect$m <- c(2, 2)
      }
      if (is.null(effect$k)) {
        effect$k <- 10
      }
      if (is.null(effect$mc)) {

        effect$mc <- TRUE
      }

      effect <- list(cov = effect$cov, by = effect$by, bs = effect$bs, m = effect$m, k = effect$k, mc = effect$mc)
      varying_coefficients[[l]] <- effect
    } else {
      effect <- append(effect, rep(list(NULL), 6-length(effect)))
      names(effect) <- c("cov", "by", "bs", "m", "k", "mc")
      if (is.null(effect$bs)) {
        effect$bs <- "ps"
      }
      if (is.null(effect$m)) {
        effect$m <- c(2, 2)
      }
      if (is.null(effect$k)) {
        effect$k <- 10
      }
      if (is.null(effect$mc)) {

        effect$mc <- TRUE
      }

      varying_coefficients[[l]] <- effect
    }
    l <- l+1

    if (is.null(effect[6][[1]])) {
      mc <- TRUE
    }
    else{
      mc <- effect[6][[1]]
    }
    f_function_var_coef <- paste0(f_function_var_coef, "+ ti(", effect[1], ", ",
                                  y, ", bs = c(", '"', effect[3][[1]], '"', ", ",
                                  '"md")', ", m = list(", deparse(effect[4][[1]]),
                                  ", ", deparse(m_continuous), "), k  = c( ",
                                  effect[5], ", ", k_continuous, ")",
                                  ", mc = c(", mc, ", FALSE), np = FALSE, by = ",
                                  effect[2], ", sp = array(c(-1, ", sp_y_vec,
                                  " )), fx = ", fx, ", xt = list(NULL, list(values_discrete = ",
                                  deparse(values_discrete), ",  domain_continuous = ",
                                  deparse(domain_continuous), ", weights_discrete = ",
                                  deparse(weights_discrete), ", penalty_discrete =  ",
                                  penalty_discrete, ")))")
    j <- j+1
  }
  l <- 1
  for (effect in smooth_interactions) {
    if (length(sp_y) > 1) {
      stopifnot("If sp_y is supplied as vector, its length must equal the number of partial effect (including the intercept)." = j <= length(sp_y))
      sp_y_ <- sp_y[j]
      sp_y_vec <- sp_y[j]
    }
    if (!is.null(names(effect))) {
      if (!all(names(effect)%in%c("covs", "bs", "m", "k", "mc", "by"))) {
        stop("Names of smooth interaction effect incorrect!")
      }
      n_c <- length(effect$covs)
      if (is.null(effect$bs)) {
        effect$bs <- rep("ps", n_c)
      }
      if (is.null(effect$m)) {
        effect$m <- rep(list(c(2, 2)), n_c)
      }
      if (is.null(effect$k)) {
        effect$k <- rep(10, n_c)
      }
      if (is.null(effect$mc)) {

        effect$mc <- rep(TRUE, n_c)
      }
      if (is.null(effect$by)) {
        effect$by <- NULL
      }
      effect <- list(covs = effect$covs, bs = effect$bs, m = effect$m, k = effect$k, mc = effect$mc, by = effect$by)
      smooth_interactions[[l]] <- effect
    } else {
      effect <- append(effect, rep(list(NULL), 6-length(effect)))
      names(effect) <- c("covs", "bs", "m", "k", "mc", "by")
      n_c <- length(effect$covs)
      if (is.null(effect$bs)) {
        effect$bs <- rep("ps", n_c)
      }
      if (is.null(effect$m)) {
        effect$m <- rep(list(c(2, 2)), n_c)
      }
      if (is.null(effect$k)) {
        effect$k <- rep(10, n_c)
      }
      if (is.null(effect$mc)) {

        effect$mc <- rep(TRUE, n_c)
      }
      if (is.null(effect$by)) {
        effect$by <- NULL
      }
      smooth_interactions[[l]] <- effect
    }

    l <- l+1
    if (is.null(effect[[5]])) {
      mc <- paste0(rep("TRUE", length(effect[[1]])), collapse = ", ")
    }
    else{
      mc <- paste0(effect[[5]], collapse = ", ")
    }
    if (is.null(effect[6][[1]])) {
      f_flex_interact <- paste0(f_flex_interact, "+ ti(",
                                paste0(effect[[1]], collapse = ", "), ", ", y,
                                paste0(", bs = c('", paste0(effect[[2]], collapse = "', '"), "', 'md')") ,
                                ", m = list(", paste0(effect[[3]], collapse = ", "),
                                ", ", deparse(m_continuous), "), k  = c( ",
                                paste(effect[[4]], collapse = ", "), ", ", k_continuous,
                                ")", ", mc = c(", mc, ", FALSE), np = FALSE, sp = array(c(",
                                paste0(rep("-1", length(effect[[1]])), collapse = ", "),
                                ", ", sp_y_vec, " )), fx = ", fx, ", xt = list(",
                                paste0(rep("NULL", length(effect[[1]])), collapse = ", "),
                                ", list(values_discrete = ", deparse(values_discrete),
                                ",  domain_continuous = ", deparse(domain_continuous),
                                ", weights_discrete = ", deparse(weights_discrete),
                                ", penalty_discrete = ", penalty_discrete, ")))")
      j <- j+1}
    else{
      f_flex_interact <- paste0(f_flex_interact, "+ ti(",
                                paste0(effect[[1]], collapse = ", "), ", ", y,
                                paste0(", bs = c('", paste0(effect[[2]], collapse = "', '"), "', 'md')") ,
                                ", m = list(", paste0(effect[[3]], collapse = ", "),
                                ", ", deparse(m_continuous), "), k  = c( ",
                                paste(effect[[4]], collapse = ", "), ", ",
                                k_continuous, ")", ", mc = c(", mc,
                                ", FALSE), np = FALSE, sp = array(c(",
                                paste0(rep("-1", length(effect[[1]])), collapse = ", "),
                                ", ", sp_y_vec, " )), fx = ", fx, ", by = ",
                                effect[6], ", xt = list(",
                                paste0(rep("NULL", length(effect[[1]])), collapse = ", "),
                                ", list(values_discrete = ", deparse(values_discrete),
                                ",  domain_continuous = ", deparse(domain_continuous),
                                ", weights_discrete = ", deparse(weights_discrete),
                                ", penalty_discrete = ", penalty_discrete, ")))")
      j <- j+1
    }
  }
  checking_dens_reg_2(group_specific_intercepts, linear_effects, smooth_effects,
                      varying_coefficients, smooth_interactions, dta)

  f <- paste0(f_main, f_intercepts, f_flexibles, f_linear, f_function_var_coef,
              f_flex_interact, "+ as.factor(group_id) -1 + offset(log(Delta) + gam_offsets)")
  if (!is.null(sample_weights)) {
    m <-
      do.call(
        "gam",
        list(
          formula = stats::as.formula((f)),
          data = as.name("dta_est"),
          method = method,
          family = stats::poisson(),
          weights = dta_est$gam_weights,
          ...
        )
      )
  }
  else{
    m <-
      do.call(
        "gam",
        list(
          formula = stats::as.formula((f)),
          data = as.name("dta_est"),
          method = method,
          family = stats::poisson(),
          ...
        )
      )
  }
  X <- stats::model.matrix(m)
  # remove intercepts per covariate combination (not of interest for estimated densities)
  intercepts <- which(grepl("group_id", colnames(X)))
  X <- X[, -intercepts]
  theta_hat <- m$coefficients[-intercepts]
  f_hat_clr <- X %*% theta_hat
  if (is.numeric(y)) {
    densi <- colnames(dta)[y]
  } else{
    densi <- y
  }
  dta_est <- as.data.table(dta_est)
  obs_density <- unique(dta_est[, ..densi])
  f_hat_clr <-
    matrix(f_hat_clr,
           nrow = nrow(obs_density),
           ncol = length(intercepts))
  if (isFALSE(values_discrete)) {
    if (!is.null(ncol(f_hat_clr))) {
      f_hat <-
        apply(
          f_hat_clr,
          MARGIN = 2,
          FUN = clr,
          inverse = TRUE,
          w = dta_est$Delta[1:nrow(obs_density)]
        )
    } else{
      f_hat <-
        clr(f_hat_clr, inverse = TRUE, w = dta_est$Delta[1:nrow(obs_density)])
    }
  }
  if (colnames(dta_est)[2] == "weighted_counts") {
    t <- unique(dta_est[, 3])}
  else{
    t <- unique(dta_est[, 2])
  }
  if (!isFALSE(values_discrete) & !isFALSE(domain_continuous)) {
    # interval_width <-
    #   (domain_continuous[2] - domain_continuous[1]) / bin_number
    # t <- seq(
    #   from = domain_continuous[1] + 1 / 2 * interval_width,
    #   to = domain_continuous[2] - 1 / 2 * interval_width,
    #   by = interval_width
    # )

    if (is.null(bin_number)) {
      bin_number <- nrow(t)-length(values_discrete)
    }
    discrete_ind <- match(values_discrete, t)
    if (!is.null(ncol(f_hat_clr))) {
      f_hat <-
        apply(
          f_hat_clr,
          MARGIN = 2,
          FUN = clr,
          inverse = TRUE,
          w = dta_est$Delta[1:nrow(obs_density)]
        )
    } else{
      f_hat <-
        clr(f_hat_clr, inverse = TRUE, w = dta_est$Delta[1:nrow(obs_density)])
    }
  }
  if (!isFALSE(values_discrete) & isFALSE(domain_continuous))
  {

    if (!is.null(ncol(f_hat_clr))) {
      f_hat <-
        apply(
          f_hat_clr,
          MARGIN = 2,
          FUN = clr,
          inverse = TRUE,
          w = dta_est$Delta[1:nrow(obs_density)]
        )

    } else{
      f_hat <-
        clr(f_hat_clr, inverse = TRUE, w = dta_est$Delta[1:nrow(obs_density)])
    }
  }
  if (isFALSE(values_discrete) & !isFALSE(domain_continuous))
  {    if (is.null(bin_number)) {
    bin_number <- nrow(t)
  }
  }
  if (!effects) {
    result <-    list(
      count_data = dta_est,
      model = m,
      model_matrix = X,
      theta_hat = theta_hat,
      f_hat_clr = f_hat_clr,
      f_hat = f_hat,
      params = list(domain_continuous = domain_continuous,
                    values_discrete = values_discrete, G = bin_number),
      predicted_effects = list(
        group_specific_intercepts = group_specific_intercepts,
        smooth_effects = smooth_effects,
        linear_effects = linear_effects,
        varying_coefficients = varying_coefficients,
        smooth_interactions = smooth_interactions
      ),
      ID_covCombi = cov_combi_id
    )
    attr(result, "class") <- c("dens_reg_obj", class(result))
    return(result)
  } else {
    pred_terms = stats::predict(m, type = "terms")
    G <- nrow(obs_density)
    n_singles <- length(group_specific_intercepts)
    levels_singles <- c()
    dta <- as.data.table(dta)
    ## find
    for (single in group_specific_intercepts) {
      n <- length(unlist(unique(dta[, ..single])))
      levels_singles <- append(levels_singles, n)
    }
    positions_singles <-
      c(2:(2 + sum(levels_singles) - length(levels_singles)))
    if ((max(positions_singles)) < length(colnames(pred_terms))) {
      positions_smooths <-
        c((max(positions_singles) + 1):length(colnames(pred_terms)))
    } else{
      positions_smooths <- c()
    }
    smooth_cols <- c()
    for (effect in smooth_effects) {
      #
      if (is.null(effect[6][[1]])) {
        ind <- match(effect[1][[1]], colnames(dta_est))
        smooth_cols <- append(smooth_cols, ind)
      } else{
        ind <- match(c(effect[1][[1]], effect[6][[1]]), colnames(dta_est))
        smooth_cols <- append(smooth_cols, list(ind))
      }
    }
    for (effect in linear_effects) {
      ind <- match(c(effect), colnames(dta_est))
      smooth_cols <- append(smooth_cols, ind)
    }
    for (effect in varying_coefficients) {
      ind <- match(c(effect[1][[1]], effect[2][[1]]), colnames(dta_est))
      ind <- c(ind, "cont")
      smooth_cols <- append(smooth_cols, list(ind))
    }
    for (effect in smooth_interactions) {
      if (is.null(effect[6][[1]])) {
        number_covs <- length(effect[[1]])
        ind <- match(effect[[1]], colnames(dta_est))
        ind <- c(ind, "inter")
        smooth_cols <- append(smooth_cols, list(ind))
      } else {
        ind <- match(c(effect[[1]], effect[6][[1]]), colnames(dta_est))
        ind <- c(ind, "inter_by")
        smooth_cols <- append(smooth_cols, list(ind))
      }
    }
    effects <-
      get_estimated_model_effects(
        m,
        pred_terms = pred_terms,
        G = bin_number,
        positions_singles = positions_singles,
        positions_smooths = positions_smooths,
        smooth_cols = smooth_cols,
        origData = dta_est,
        positions_random_effects = NA,
        domain_continuous = domain_continuous,
        values_discrete = values_discrete,
        weights_discrete = weights_discrete
      )
    result <-
      list(
        count_data = dta_est,
        model = m,
        model_matrix = X,
        theta_hat = theta_hat,
        f_hat_clr = f_hat_clr,
        f_hat = f_hat,
        effects = effects,
        params = list(domain_continuous = domain_continuous,
                      values_discrete = values_discrete, G = bin_number),
        predicted_effects = list(
          group_specific_intercepts = group_specific_intercepts,
          smooth_effects = smooth_effects,
          linear_effects = linear_effects,
          varying_coefficients = varying_coefficients,
          smooth_interactions = smooth_interactions
        ),
        ID_covCombi = cov_combi_id
      )
    attr(result, "class") <- c("dens_reg_obj", class(result))
    return(result)
  }
}

#' Check for validity of parameters of \code{\link{dens_reg}}
#'
#' These functions check validity of parameters for \code{\link{dens_reg}}.
#' Returns an error if any parameters requirement is violated.
#' @encoding UTF-8
#' @import data.table
#' @noRd
#'
#' @inheritParams dens_reg
checking_dens_reg_1 <- function(m_continuous, k_continuous, sp_y, penalty_discrete,
                                effects, dta) {

    # check gam-parameters
    if (!(length(m_continuous) == 2) &&
        !all(check_natural_number_vec(m_continuous))) {
      stop("m_continuous must be a vector of 2 natural numbers.")
    }

    if (!check_natural_number(k_continuous)) {
      stop("k_continuous must be a natural number.")
    }

    if (!is.numeric(sp_y)&!is.null(sp_y)) {
      stop("sp_y must be numeric or NULL.")
    }

    if (!check_integer_or_null(penalty_discrete)) {
      stop("penalty_discrete must be an integer or NULL.")
    }

    # check effects
    if (!is.logical(effects) || length(effects) !=  1) {
      stop("effects must be TRUE or FALSE.")
    }

    return(TRUE)
}


#' Checks used for parameter checks of \code{\link{dens_reg}}
#'
#' These functions check if an object is of a specific type.
#'
#' @encoding UTF-8
#' @noRd
#'
#' @param x An object to be checked
check_natural_number <- function(x) {
  return(is.numeric(x) && length(x) == 1 && x == floor(x) && x > 0)
}

#' @describeIn check_natural_number check if x is vector of natural numbers
check_natural_number_vec <- function(x) {
  return(is.numeric(x) && all(x == floor(x)) && all(x > 0))
}

#' @describeIn check_natural_number check if x is integer (or NULL)
check_integer_or_null <- function(x) {
  return(is.null(x) ||
           (is.numeric(x) && length(x) == 1 && x == floor(x)))
}

#' @describeIn checking_dens_reg_1 Check for validity of parameters of \code{\link{dens_reg}}
checking_dens_reg_2 <- function(group_specific_intercepts, linear_effects,
                                smooth_effects, varying_coefficients,
                                smooth_interactions, dta) {

  # check format of group_specific_intercepts
  if (!is.null(group_specific_intercepts)) {
    if (!all(sapply(group_specific_intercepts, is.character))) {
      stop("group_specific_intercepts must be NULL or characters.")
    }

    if (!all(group_specific_intercepts %in% colnames(dta))) {
      stop("All covariates in group_specific intercepts must be variables contained in dta.")
    }

    for (col in group_specific_intercepts) {
      if (!is.factor(dta[[col]])) {
        stop(paste("Variable", col, "in dta must be of class 'factor'."))
      }
    }
  }

  # check linear_effects
  if (!is.null(linear_effects)) {
    if (!all(sapply(linear_effects, is.character))) {
      stop("linear_effects must be NULL or characters.")
    }

    if (!all(linear_effects %in% colnames(dta))) {
      stop("All covariates in linear_effects must be variables contained in dta.")
    }

    for (col in linear_effects) {
      if (!is.numeric(dta[[col]])) {
        stop(paste("Variable", col, "in dta must be of class 'numeric'."))
      }
    }
  }

  # check function for structure of smooth_effects
  check_smooth_effects_structure <- function(effects) {
    if (!is.list(effects)) {
      stop("smooth_effects must be a list of lists.")
    }

    for (effect in effects) {
      # if (length(effect) < 5 || length(effect) > 6) {
      #   stop("Each smooth_effects list must have 4 or 5 elements.")
      # }

      if (!is.character(effect[[1]]) || !is.numeric(dta[[effect[[1]]]])) {
        stop("(First) Element cov of each smooth_effect must be a variable contained in dta and the corresponding variable in dta must be numeric.")
      }

      if (!is.character(effect[[2]])) {
        stop(
          "(Second) Element bs of each smooth_effect must be a character which specifies a basis (see bs in ?mgcv::ti)."
        )
      }

      if (!check_natural_number_vec(effect[[3]]) || !check_natural_number(effect[[4]])) {
        stop(
          "(Third and fourth) Elements m and k of each smooth_effect must be natural numbers specifying m and k for the covariate direction."
        )
      }
      if (!is.logical(effect[[5]])) {
        stop(
          paste(
            "(Fifth) Element mc of each smooth_effect must be a variable contained in dta and the corresponding variable must be numeric."
          )
        )
      }

      if (length(effect) ==  6 &&
          !is.null(effect[[6]]) &&
          (!is.character(effect[[6]]) ||
           !is.factor(dta[[effect[[6]]]]))) {
        stop(
          paste(
            "(Sixth) Element by of each smooth_effect must be a variable contained in dta and the corresponding variable must be factor."
          )
        )
      }
    }
  }

  # check smooth_effects
  if (!is.null(smooth_effects)) {
    check_smooth_effects_structure(smooth_effects)
  }

  # check function for structure of varying_coefficients
  check_varying_coefficients_structure <-
    function(effects) {
      if (!is.list(effects)) {
        stop("varying_coefficients must be a list of lists.")
      }

      for (effect in effects) {
        if (length(effect) !=  6) {
          stop("Each varying_coefficients list must have 6 elements.")
        }

        if (!is.character(effect[[1]]) ||
            !is.numeric(dta[[effect[[1]]]])) {
          stop(
            paste(
              "(First) Element cov of each varying_coefficient must be a variable contained in dta and the corresponding variable must be numeric."
            )
          )
        }

        if (!is.character(effect[[2]]) || !is.numeric(dta[[effect[[2]]]])) {
          stop(
            paste(
              "(Second) Element by of each varying_coefficient must be a variable contained in dta and the corresponding variable must be numeric."
            )
          )
        }

        if (!is.character(effect[[3]])) {
          stop(
            "(Third) Element bs of each varying_coefficient must be a character which specifies a basis (see bs in ?mgcv::ti)."
          )
        }

        if (!check_natural_number_vec(effect[[4]]) ||
            !check_natural_number(effect[[5]])) {
          stop(
            "(Fourth and fifth) Elements m and k of each varying_coefficient must be natural numbers specifying m and k for the covariate direction."
          )
        }

        if (!is.logical(effect[[6]])) {
          stop(
            paste(
              "(Sixth) Element mc of each smooth_effect must be a variable contained in dta and the corresponding variable must be numeric."
            )
          )
        }
      }
    }

  # check varying_coefficients
  if (!is.null(varying_coefficients)) {
    check_varying_coefficients_structure(varying_coefficients)
  }

  # check function for structure of smooth_interactions
  check_smooth_interactions_structure <- function(effects) {
    if (!is.list(effects)) {
      stop("smooth_interactions must be a list of lists.")
    }

    for (effect in effects) {
      # if (!(
      #   length(effect[[1]]) == length(effect[[2]]) &&
      #   length(effect[[3]]) == length(effect[[4]]) &&
      #   length(effect[[1]]) == length(effect[[4]])
      # )) {
      # l <- lengths(effect)
      # l <- l[which(l > 0)]
      if (!diff(range((lengths(effect[which(lengths(effect) > 0)])))) == 0) {
        stop("Unequal number of elements in the lists for a smooth interaction effect.")
      }
      if (!is.character(effect[[1]]) || !is.numeric(unlist(dta[, effect[[1]]]))) {
        stop("(First) Element covs of each smooth_interaction must be a character vector of variable names contained in dta and the corresponding variables must be numeric.")
      }

      if (!is.character(effect[[2]])) {
        stop("(Second) Element bs of each smooth_interaction must be a character vector, each specifying a basis (see bs in ?mgcv::ti).")
      }

      if (!check_natural_number_vec(unlist(effect[[3]])) || !is.list(effect[[3]]) || any(lengths(effect[[3]]) != 2)) {
        stop(
          "(Third) Elements m of each smooth_interaction must be a list of vectors of two integers specifying m for all covariate marginals."
        )
      }

      if (!check_natural_number_vec(unlist(effect[[4]]))) {
        stop(
          "(Fourth) Element k of each smooth_interaction must be a vector of natural numbers specifying k for all covariate marginals."
        )
      }

      if (!is.logical(effect[[5]])) {
        stop(
          paste(
            "(Fifth) Element mc of each smooth_interaction must be a logical vector"
          )
        )
      }

      if (length(effect) ==  6 && !is.null(effect[[6]]) &&
          (!is.character(effect[[6]]) || !is.factor(dta[[effect[[6]]]]))) {
        stop(
          paste(
            "(Sixth) Element by of each smooth_interaction must be a variable contained in dta and the corresponding variable must be factor."
          )
        )
      }
    }
  }

  # check smooth_interactions
  if (!is.null(smooth_interactions)) {
    check_smooth_interactions_structure(smooth_interactions)
  }
  return(TRUE)
}



