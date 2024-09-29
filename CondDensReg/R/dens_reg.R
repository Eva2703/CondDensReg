#' Conditional density regression for individual-level data
#'
#' This function implements the approach by Maier et al. (2023) for fitting
#' structured additive regression models with densities in a mixed (continuous/
#' discrete) Bayes Hilbert space
#' \eqn{B^2(\mu) = B^2(\mathcal{Y}, \mathcal{A}, \mu)}{B^2(\mu) = B^2(Y, A, \mu)}
#' as response given scalar covariates \eqn{x_i},
#' based on observations \eqn{y_i} from the conditional distributions given
#' \eqn{x_i}, via a (penalized) maximum likelihood approach. The (penalized)
#' log-likelihood function is approximated via the (penalized) log-likelihood of
#' an appropriate poisson regression model, which is then fitted using
#' \code{mgcv}'s \code{\link[mgcv]{gam}} with a new smooth for mixed densities.
#' We briefly summarize the approach below in the details.
#' % to enable proper description of the function.
#' Please see Maier et al. (2023) for comprehensive discussion.
#'
#' The goal of \code{dens_reg} is to estimate the densities \eqn{f_{x_i}} of
#' conditional distributions of random variables \eqn{Y_i | x_i} from independent
#' observations \eqn{(y_i, x_i)}, where \eqn{y_i} are realizations of \eqn{Y_i | x_i},
#' and \eqn{x_i} is a vector of covariate observations, \eqn{i = 1,...,N}. The
#' densities are considered as elements of a mixed Bayes Hilbert space
#' \eqn{B^2(\mu) = B^2(\mathcal{Y}, \mathcal{A}, \mu)}{B^2(\mu) = B^2(Y, A, \mu)}
#' (including continuous and discrete ones as special cases). The subset of the
#' domain \eqn{\mathcal{Y}}{Y} corresponding to the continuous part of the densities
#' is denoted with \eqn{I}, the one corresponding to the discrete part with
#' \eqn{\mathcal{D}}{D}. The densities are modeled via a structured additive
#' regression model
#' \deqn{f_{x_i} = \bigoplus_{j=1}^J h_j (x_i)}
#' with partial effects \eqn{h_j (x_i) \in B^2(\mu)} depending on no, one, or
#' several covariates \eqn{x_i}. Each partial effect is represented using a
#' tensor product basis, consisting of an appropriate vector of basis functions
#' \eqn{b_{\mathcal{Y}}{b_Y}}{b_Y} in \eqn{B^2(\mu)} over the domain of
#' \eqn{B^2(\mu)}, and a vector of basis functions \eqn{b_{\mathcal{X}, j}}{b_{X, j}}
#' over the respective covariates. The basis functions \eqn{b_{\mathcal{Y}}}{b_Y}
#' are obtained by transforming a P-spline basis from \eqn{L^2(\mu)}
#' to
#' \eqn{L^2_0(\mu) := \{ f \in L^2(\mu) ~|~ \int_{\mathcal{T}} f \, \mathrm{d}\mu = 0\}}{L^2_0(\mu) := {f \in L^2(\mu) | \int_T f d\mu = 0}}
#' and then applying the inverse centered log-ratio (clr) transformation (the
#' clr transformation is an isometric isomorphism, allowing to consider densities
#' in \eqn{B^2(\mu)} equivalently in \eqn{L^2_0(\mu)}). See appendix D of Maier
#' et al. (2021) for details on the transformation from \eqn{L^2(\mu)} to
#' \eqn{L^2_0(\mu)}. The choice of \eqn{b_{\mathcal{X}, j}}{b_{X, j}} determines
#' the type of the partial effect (e.g., linear or smooth for a continuous covariate).
#' Penalization is also possible.
#' The resulting log-likelihood is then approximated by the log-likelihood of a
#' corresponding multinomial model, which can equivalently be estimated by a
#' poisson model. The data for these models is obtained from the original
#' observations \eqn{y_1,...,y_N} by combining all observations of the same
#' conditional distribution into a vector of counts via a histogram on the
#' continuous part of the domain of the densities and counts on the discrete
#' part of the domain. For details, see Maier et al. (2023).
#' \code{dens_reg} performs the transformation into appropriate count data and
#' estimates the corresponding poisson model. Furthermore, the resulting densities on
#' pdf- and clr-level as well as the estimated partial effects on both levels are
#' calculated.
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
#' of response and covariates as well as optional sampling weights for each
#' observation in the rows (\eqn{i = 1,...,N}).
#' @param var_vec Vector of variables in \code{dta} on which the covariate
#' combinations are based. (All covariates or a subset in case some covariates
#' are a refinement of others.) The vector can either contain the variable names
#' as strings or the column positions of the respective variables in \code{dta}.
#' @param density_var Variable in \code{dta} which corresponds to the response.
#' Variable name as string or numeric column position of the variable in \code{dta}.
#' @param sample_weights Variable in \code{dta} which contains a sample weight
#' for each observation. Variable name as string or numeric column position of
#' the variable in \code{dta}. If missing (\code{NULL}), no sample weights are
#' included per default.
#' @param bin_width Width of histogram bins partitioning \eqn{I} (the subset of
#' the domain corresponding to the continuous part of the densities). Can be one
#' scalar value (specifying an equidistant bin width), or a vector containing the
#' width of each bin. The combined length of the specified bins must match the
#' length of the continuous part of the domain. Alternative to \code{bin_number}.
#' @param bin_number Number of equidistant histogram bins partitioning \eqn{I}.
#' Alternative to \code{bin_width}. If \code{bin_number} and \code{bin_width} are
#' both given and the two values are not compatible, an error is returned. If
#' neither parameter is specified, \code{bin_number=100} is used as default.
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
#' \eqn{I}. If missing (\code{NULL}) it is set to \code{c(0, 1)} as default. If
#' set to \code{FALSE}, the continuous component is considered to be empty, i.e.,
#' a weighted sum of dirac measures is used as reference measure (discrete special
#' case).
#' @param m_density_var Vector of two integers specifying the order of the B-spline
#' basis over \eqn{I} and the order of the difference penalty (like the argument
#' \code{m} for P-Spline smooth terms, i.e., \code{bs="ps"} in \code{\link[mgcv]{ti}},
#' etc.) for basis functions in \eqn{L^2(\lambda)} (before transformation to \eqn{L^2_0(\lambda)},
#' compare details). If missing it is set to cubic splines with second
#' order difference penalty, i.e., \code{m_density_var=c(2,2)}, as default.
#' @param k_density_var Integer specifying the number of B-spline basis functions
#' over \eqn{I}, i.e., the length of the vector \eqn{b_{\mathcal{Y}}{b_Y}}{b_Y}
#' (like the argument \code{k} for P-Spline smooth terms, i.e., \code{bs="ps"} in
#' \code{\link[mgcv]{ti}}, etc.) for basis functions in \eqn{L^2(\lambda)}
#' (before transformation to \eqn{L^2_0(\lambda)}, compare details). See
#' \code{\link[mgcv]{choose.k}} for more information.
#' If missing (\code{NULL}) it is set to 10.
#' @param sp_density_var Integer or vector specifying the smoothing parameter for the marginal
#' penalty matrix for \eqn{b_{\mathcal{Y}}{b_Y}}{b_Y} (anisotropic penalty; like
#' the argument \code{sp} in \code{\link[mgcv]{ti}}, etc.). If a vector is submitted, all smoothing parameters
#' in density direction are treated as fixed but with different values for the different partial effects. The
#' order of parameters within this vector is the same as the order of partial effects as arguments of \code{dens_reg} with the firtst parameter corresponding to the intercept.
#' It must be of the same length as the number of partial effects plus one.
#' In the discrete case of only discrete values and a zero \code{penalty_discrete}, the smoothing
#' parameter \code{sp_density_var} is set to zero in order to obtain an unpenalized
#' model. If missing (\code{NULL}) or a negative value is supplied, the parameter
#' will be estimated by \code{\link[mgcv]{gam}}. Any positive or zero value is
#' treated as fixed. (see \code{mgcv})
#' @param penalty_discrete Integer or \code{NULL} giving the order of differences
#' to be used for the penalty of the discrete component, with 0 corresponding to
#' the identity matrix as penalty matrix (analogously to m[2] for the continuous
#' component); Note that the order of differences must be smaller than the number
#' of values in the discrete component, i.e., the length of values_discrete; if
#' missing (\code{NULL}), in the mixed case, it is set to a zero matrix (corresponding
#' to no penalization), while in the discrete case, it is set to a diagonal matrix
#' (as in a ridge penalty). In the discrete case, please set the corresponding
#' argument of sp in ti() to 0 to estimate an unpenalized model.
#' @param group_specific_intercepts Vector of the form \code{c("cov_1","cov_2",...)}
#' of names of categorical variables included by \code{var_vec} for which
#' covariate-specific intercepts \eqn{\beta_{cov_1},\ \beta_{cov_2},\ ...} are
#' included in the model. If mising (\code{NULL}), no covariate-specific intercept
#' is included.
#' @param linear_effects Vector of the form \code{c("cov_1","cov_2",...)} of names
#' of numeric variables included by \code{var_vec} for which linear effects
#' \eqn{\beta_{cov_1}\cdot cov_1,\ \beta_{cov_2}\cdot cov_2,\ ...} are included
#' in the model. If mising (\code{NULL}), no linear is included.
#' @param flexible_effects List of lists of the form
#' \code{list(list("cov_1","basis_1", m_1,k_1,mc1,"by_1"),list(...),...)}. Each list
#' is adding one flexible effect of the form \eqn{f_{by}(cov)} to the model with:
#' \itemize{
#' \item \code{"cov"}: Name of a numeric variable included in \code{var_vec}.
#' \item \code{"basis"}: Character string specifying the type for the marginal
#' basis in covariate direction. See \code{mgcv::smooth.terms} for details and
#' full list.
#' \item \code{m}: Vector of two integers or \code{NULL} giving the order of the
#' marginal spline basis \eqn{b_j} in covariate direction and the order of its
#' penalty. (see \code{mgcv})
#' \item \code{k}: Integer or \code{NULL} giving the dimension of the marginal
#' basis \eqn{b_j}. See \code{mgcv::choose.k} for more information.
#'  \item \code{mc}: Logical indicating if the marginal in covariate direction should have centering constraints applied. By default all marginals are constrained, i.e., \code{mc=TRUE}.
#'  \item \code{"by"}: Optional, name of a categorical variable included in
#' \code{var_vec} specifying whether the flexible effect should be modeled
#' specifically for each level of the by-covariate. If missing or \code{NULL},
#' the flexible effect is not depending on the level of an additional covariate.
#' } If mising (\code{NULL}), no flexible effect is included.
#' @param varying_coefficients  List of lists of the form
#' \code{list(list("covA_1","covB_1","basis1",m1,k1, mc1),...)}. Each list is adding
#' one varying coefficient of the form \eqn{ cov_A*f(cov_B)} to the
#' model with:
#' \itemize{
#' \item \code{"covA"}: Name of a numeric variable included in \code{var_vec}.
#'  \item \code{"covB"}: Name of a numeric variable included in \code{var_vec}.
#' \item \code{"basis"}: Character string specifying the type for the marginal
#' basis in covariate direction. See \code{mgcv::smooth.terms} for details and
#' full list.
#' \item \code{m}: Vector of two integers or \code{NULL} giving the order of the
#' marginal spline basis \eqn{b_j} in covariate direction for the smooth effect
#' and the order of its penalty. (see \code{mgcv})
#' \item \code{k}: Integer or \code{NULL} giving the dimension of the marginal
#' basis \eqn{b_j} for the smooth effect. See \code{mgcv::choose.k} for more
#' information.
#' }If mising (\code{NULL}), no varying coeffecient is included.
#' \item \code{mc}: Logical indicating if the marginal in the direction of the first covariate should have centering constraints applied. By default all marginals are constrained, i.e., \code{mc=TRUE}.
#' @param flexible_interaction  List of lists of the form
#' \code{list(list(c("covA_1","covB_1",...),c("basisA_1","basisB_1",...), list(mA_1, mB_1,...), c(kA_1, kB_1,...), c(mcA_1, mcB_1,...), "by"),list(...),...)}.
#' Each list is specifying flexible interaction effect between at least two
#' continuous covariates of the form \eqn{f(cov_A,cov_B,...)} in the model with:
#' \itemize{
#' \item \code{c("covA_1","covB_1", ...)}: Vector of names of numeric variables
#' included in \code{var_vec}.
#' \item \code{c("basisA", "basisB", ...}:Vector of character strings specifying
#' the type for the marginal bases in covariate direction. See
#' \code{mgcv::smooth.terms} for details and full list.
#' \item\code{list(mA,mB,...)}: List containing vectors of two integers or
#' \code{NULL} giving the order of the marginal spline basis in direction of
#' \eqn{cov_A} for the smooth effect and the order of its penalty. (see \code{mgcv})
#' \item \code{c(kA, kB, ...)}: Vector of integers or \code{NULL} giving the
#' dimension of the marginal basis in direction of \eqn{cov_A} for the smooth
#' effect. See \code{mgcv::choose.k} for more information.
#'  \item \code{mc}: Logical indicating if the marginal in covariate direction should have centering constraints applied. By default all marginals are constrained, i.e., \code{mc=TRUE}.
#'   \item \code{"by"}: Optional, name of a categorical variable included in
#' \code{var_vec} specifying whether the flexible interaction effect should be modeled
#' specifically for each level of the by-covariate. If \code{NULL},
#' the flexible interaction effect is not depending on the level of an additional covariate.
#' }
#' @param effects Indicates if partial effects should be returned (\code{TRUE})
#' or not (\code{FALSE}).
#' @param ...  further arguments for passing on to \code{gam}
#'
#'
#' @return The returned object is a \code{list}-object which is also an object
#' of the sub-class \code{dens_reg_obj} with elements:
#' \itemize{
#' \item \code{model}: \code{mgcv::gam}-object of the estimated model.
#' \item \code{model_matrix}: Model matrix of the model.
#' \item \code{theta_hat}: Estimated coeffecient vector \eqn{\hat\theta}.
#' \item \code{f_hat_clr}: Estimated conditional densities on clr-level
#' \eqn{clr(\hat f)} for every covariate combination in the order of the respective
#' group_id in \code{histo_data}.
#' \item \code{f_hat}: Estimated conditional densities on clr-level \eqn{clr(\hat f)}
#' for every covariate combination in the order of the respective group_id in
#' \code{histo_data}.
#' \item \code{histo_data}: \code{data.table}-object containing the histogram data
#' resulting from the \code{preprocess}-function with the original data \code{dta}
#' which is also an object of the sub-class \code{histogram_count_data}. See
#' \code{?preprocess} for more information.
#' \item \code{effects}: If \code{effects=TRUE}; List of lists. Each list gives
#' one estimated partial effect. Both clr- and density-level are included.
#' \item \code{params}: List of \code{domain_continuous, values_discrete} and
#' \code{bin_number}.
#' \item \code{predicted_effects}: List of lists (\code{group_specific_intercepts,
#' flexible effects, linear_effects, varying_coefficient, flexible_interaction})
#' collecting all modeled partial effects as specified in the function's parameters.
#' \item \code{ID_covCombi}: Data frame which gives an overview over the assignment
#' of the unique covariate combinations to the group IDs.
#' }
#' @examples
#' \donttest{# for further information on the parameters of the preprocessing step
#' # see ?preprocess
#'
#' library(dplyr)
#'
#' # create data (mixed)
#'
#'dta <- data.frame(obs_density = sample(0:2, 150, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
#'                   covariate1 = sample(c("a", "b", "c"), 150, replace = TRUE),
#'                   covariate2 = sample(c("c", "d"), 150, replace = TRUE),
#'                   covariate3 = rep(rnorm(n = 15), 10),
#'                   covariate4 = rep(rnorm(n = 10), 15),
#'                   covariate5=rep(rnorm(n = 10), 15), sample_weights = runif(150, 0, 2))
#'dta[which(dta$obs_density == 2),]$obs_density <-rbeta(length(which(dta$obs_density == 2)), shape1 = 3, shape2 = 3)
#'dta$covariate1 <- ordered(dta$covariate1)
#'dta$covariate2 <- ordered(dta$covariate2)
#'
#'# create discrete data
#'
#'dta_dis <- data.frame(obs_density = sample(0:2, 150, replace = TRUE, prob = c(0.25, 0.45,0.3)),
#'                       covariate1 = sample(c("a", "b", "c"), 150, replace = TRUE),
#'                       covariate2 = sample(c("c", "d"), 150, replace = TRUE),
#'                       covariate3 = rep(rnorm(n = 15), 10),
#'                       covariate4 = rep(rnorm(n = 10), 15),
#'                       covariate5=rep(rnorm(n = 10), 15),
#'                       sample_weights = runif(150, 0, 2))
#'dta_dis$covariate1 <- ordered(dta_dis$covariate1)
#'dta_dis$covariate2 <- ordered(dta_dis$covariate2)
#'
#'# examples for different partial effects
#'
#' ## group specific intercepts
#' group_specific_intercepts <- c("covariate1", "covariate2")
#' ## linear effects
#' linear_effects <- c("covariate4")
#'## flexible effects
#' flexible_effects <-list(list("covariate3", "ps", c(2, 2), 4, NULL,NULL), list("covariate3", "ps", c(2, 2), 4,FALSE, "covariate1"))
#' ## varying coefficient
#' fvc <-list(list("covariate3", "covariate4", "ps", c(2, 2) , 4, TRUE))
#'## flexible interaction
#' flex_inter <-list(list(c("covariate3", "covariate4","covariate5"), c("ps", "ps","ps"),list( c(2, 2), c(2, 2), c(2, 2)),c( 4, 4,5), list(TRUE, FALSE, TRUE), NULL))
#'
#'# fit models (warning: calculation may take a few minutes)
#'
#' ## fit model for the mixed case with group specific intercepts and linear effects
#' ### use fixed smoothing parameters in density direction and calculate also the partial effects
#'
#' m_mixed <- dens_reg(dta = dta, var_vec = c(2:6), density_var = 1, m_density_var = c(2, 2),
#'   k_density_var = 4, group_specific_intercepts = group_specific_intercepts,
#'   linear_effects = linear_effects, effects = TRUE, sp_density_var=c(1,3,5,0.5))
#'
#' ## fit model for the discrete case with flexible effects and flexible interaction
#' ### do not calculate effects
#'
#' m_dis <- dens_reg(
#'   dta = dta_dis, var_vec = c(2:6), density_var = 1, values_discrete = c(0, 1,2),
#'   weights_discrete = c(1,1,1), domain_continuous = FALSE, m_density_var = c(2, 2),
#'   k_density_var = 4, group_specific_intercepts = group_specific_intercepts,
#'   flexible_effects = flexible_effects, flexible_interaction = flex_inter, effects = FALSE)
#'
#' # fit model for the continuous case with a functional varying coeffecient
#'
#' m_cont <- dens_reg(dta = dta%>%filter(obs_density!=0&obs_density!=1),
#'   var_vec = c(2:6), density_var = 1, values_discrete = FALSE, m_density_var = c(2, 2),
#'   k_density_var = 12, varying_coefficients = fvc, effects = TRUE)
#'
#'
#' }
#' @export
#' @references Maier, E. M., Fottner, A., Stöcker, A., Okhrin, Y., & Greven, S. (2023). Conditional density regression for individual-level data.


dens_reg <- function(dta,
                     var_vec,
                     density_var,
                     sample_weights = NULL,
                     bin_width = NULL,
                     bin_number = NULL,
                     values_discrete = c(0, 1),
                     weights_discrete = 1,
                     domain_continuous = c(0, 1),
                     m_density_var= c(2, 2),
                     k_density_var=10,
                     sp_density_var = NULL,
                     penalty_discrete = NULL,
                     group_specific_intercepts = NULL,
                     # list: c("var_name1",...)
                     ##  ti(density_var, bs = "md", m = m_density_var, k = k_densitdy var, mc = FALSE, np = FALSE, by = var_name1)
                     ## beta_k
                     linear_effects = NULL,
                     # list: c("var_name1",...)
                     ## ti(density_var, bs = "md", m = m_density_var, k = k_densitdy var, mc = FALSE, np = FALSE, by=var_name1)
                     ## x*beta
                     # linear_interaction,
                     ### not working      # list of lists: list(list("var_nameA_1","var_nameB_1"),...)
                     ## ??ti(density_var, bs = "md", m = m_density_var, k = k_densitdy var, mc = FALSE, np = FALSE)
                     ## x1*(x2*beta)
                     flexible_effects = NULL,
                     # list of lists: list(list("var_name1","basis1", m1,k1,"by1"),list(...),...) #if by=empty -> by=NULL)
                     ## ti(var_name1, density_var, bs = c(basis1, "md"), m = list(m1, m_density_var), k = c(k1, k_densityVar), mc = c(TRUE, FALSE), np = FALSE, by=by)
                     ## g(x) or g_k(x)
                     varying_coefficients = NULL,
                     # list of lists: list(list("var_nameA_1","var_nameB_1",basis1,m1,k1),...)
                     ## ti(var_nameA_1, density_var, bs = c(basis1, "md"), m = list(m1, m_density_var), k = c(k1, k_densityVar), mc = c(TRUE, FALSE), np = FALSE, by=var_nameB_1)
                     ## x1*g(x2)
                     flexible_interaction = NULL,
                     # list of lists: list(list(c("var_nameA_1","var_nameB_1",...),c("basis1A","basis1B",...), c(m1A, m1B,...), c(k1A,k1B,...)),list(...),...)
                     ## ti(var_nameA_1,var_nameB_1, density_var, bs = c(basis1A, basis1B, "md"), m = list(m1A,m1B, m_density_var), k = c(k1a,k1B, k_densityVar), mc = c(TRUE, TRUE, FALSE), np = FALSE)
                     ## g(x1,x2)
                     effects = FALSE,
                     ...)
{
  if (isFALSE(values_discrete)) {
    weights_discrete <- NULL
  }


  dta_est <- preprocess(
    dta,
    var_vec,
    density_var,
    sample_weights,
    bin_width,
    bin_number ,
    values_discrete ,
    weights_discrete ,
    domain_continuous
  )
  cov_combi_id <-
    unique(dta_est %>% select((length(
      colnames(dta_est)
    ) - 4 - length(var_vec)):(length(
      colnames(dta_est)
    ) - 4)))
  checking_dens_reg(
    m_density_var,
    k_density_var,
    sp_density_var,
    penalty_discrete,
    group_specific_intercepts,
    linear_effects,
    flexible_effects,
    varying_coefficients,
    flexible_interaction,
    effects,
    dta
  )
  if (is.numeric(density_var)) {
    density_var <- colnames(dta)[density_var]
  }
  if (is.numeric(var_vec)) {
    var_vec <- colnames(dta)[var_vec]
  }
  if (!is.null(weights_discrete)) {
    if (length(weights_discrete) == 1) {
      weights_discrete <- rep(weights_discrete, length(values_discrete))
    }
  }
  if (is.null(sp_density_var)) {
    sp_density_var_vec <- -1
    sp_density_var_<-"NULL"
  }else{

  if (length(sp_density_var)==1){
    sp_density_var_<-sp_density_var
    sp_density_var_vec<-sp_density_var
  }
  }



  if (is.null(penalty_discrete)) {
    if (isFALSE(domain_continuous)) {
      fx <- "FALSE"
      penalty_discrete <- "NULL"
      sp_density_var_ <- 0
    }
    else{
      fx <- "FALSE"
      penalty_discrete <- "NULL"
    }
  } else{
    fx = "FALSE"
  }
  ### E: Habe ich auskommentiert und lade es außerhalb der Funktion
  # source("mixed_density_smooth.R")
  # source("plain_fcts.R")
  j<-1
  if (length(sp_density_var)>1){
    sp_density_var_<-sp_density_var[j]
    sp_density_var_vec<-sp_density_var[j]
  }
  f_main <-
    paste0(
      "counts~ti(",
      density_var,
      ", bs=",
      '"md"' ,
      ",mc=FALSE, np=FALSE, m = list(",
      deparse(m_density_var),
      "), k = ",
      k_density_var,
      ",sp=",
      sp_density_var_,
      ",fx=",
      fx,
      ",xt=list(list(values_discrete=",
      deparse(values_discrete),
      ",  domain_continuous=",
      deparse(domain_continuous),
      ", weights_discrete=",
      deparse(weights_discrete),
      ",penalty_discrete= ",
      penalty_discrete,

      ")))"
    )
  f_intercepts <- ""
  f_flexibles <- ""
  f_linear <- ""
  f_function_var_coef <- ""
  f_flex_interact <- ""
  j<-j+1
  for (effect in group_specific_intercepts) {
    if (length(sp_density_var)>1){
      sp_density_var_<-sp_density_var[j]
      sp_density_var_vec<-sp_density_var[j]
    }
    f_intercepts <-
      paste0(
        f_intercepts,
        "+ ti(",
        density_var,
        ", bs=",
        '"md"' ,
        ",m = list(",
        deparse(m_density_var),
        "), k = ",
        k_density_var,
        ",mc = FALSE, np = FALSE, by=",
        effect,
        ",sp=",
        sp_density_var_,
        ",fx=",
        fx,
        ",xt=list(list(values_discrete=",
        deparse(values_discrete),
        ",  domain_continuous=",
        deparse(domain_continuous),
        ", weights_discrete=",
        deparse(weights_discrete),
        ",penalty_discrete= ",
        penalty_discrete,

        ")))"
      )
    j<-j+1
  }
  for (effect in flexible_effects) {
    if (length(sp_density_var)>1){
      sp_density_var_<-sp_density_var[j]
      sp_density_var_vec<-sp_density_var[j]
    }
    if (is.null(effect[5][[1]])){
      mc<-TRUE
    }
    else{
      mc<-effect[5][[1]]
    }
    if (is.null(effect[6][[1]])) {
      f_flexibles <-
        paste0(
          f_flexibles,
          "+ ti(",
          effect[1],
          ",",
          density_var,
          ", bs=c(",
          '"',
          effect[2][[1]],
          '"',
          ",",
          '"md")' ,
          ",m = list(",
          deparse(effect[3][[1]]),
          ",",
          deparse(m_density_var),
          "), k =c( ",
          effect[4],
          ",",
          k_density_var,
          "),mc = c(",mc,", FALSE), np = FALSE, sp=array(c(-1,",
          sp_density_var_vec,
          " ))",
          ",fx=",
          fx,
          ",xt=list(NULL,list(values_discrete=",
          deparse(values_discrete),
          ",  domain_continuous=",
          deparse(domain_continuous),
          ", weights_discrete=",
          deparse(weights_discrete),
          ",penalty_discrete= ",
          penalty_discrete,
          ")))"
        )
    }
    else{
      f_flexibles <-
        paste0(
          f_flexibles,
          "+ ti(",
          effect[1],
          ",",
          density_var,
          ", bs=c(",
          '"',
          effect[2][[1]],
          '"',
          ",",
          '"md")' ,
          ",m = list(",
          deparse(effect[3][[1]]),
          ",",
          deparse(m_density_var),
          "), k =c( ",
          effect[4],
          ",",
          k_density_var,
          ")",
          ",mc = c(",mc,", FALSE), np = FALSE, by=",
          effect[6],
          ",sp=array(c(-1,",
          sp_density_var_vec,
          " )),fx=",
          fx,
          ",xt=list(NULL,list(values_discrete=",
          deparse(values_discrete),
          ",  domain_continuous=",
          deparse(domain_continuous),
          ", weights_discrete=",
          deparse(weights_discrete),
          ",penalty_discrete= ",
          penalty_discrete,

          ")))"
        )
    }
    j<-j+1
  }
  for (effect in linear_effects) {
    if (length(sp_density_var)>1){
      sp_density_var_<-sp_density_var[j]
      sp_density_var_vec<-sp_density_var[j]
    }
    f_linear <-
      paste0(
        f_linear,
        "+ ti(",
        density_var,
        ", bs=",
        '"md"' ,
        ",m = list(",
        deparse(m_density_var),
        "), k = ",
        k_density_var,
        ",mc = FALSE, np = FALSE, by=",
        effect,
        ",sp=",
        sp_density_var_,
        ",fx=",
        fx,
        ",xt=list(list(values_discrete=",
        deparse(values_discrete),
        ",  domain_continuous=",
        deparse(domain_continuous),
        ", weights_discrete=",
        deparse(weights_discrete),
        ",penalty_discrete= ",
        penalty_discrete,

        ")))"
      )
  j<-j+1
    }
  for (effect in varying_coefficients) {
    if (length(sp_density_var)>1){
      sp_density_var_<-sp_density_var[j]
      sp_density_var_vec<-sp_density_var[j]
    }
    if (is.null(effect[5][[1]])){
      mc<-TRUE
    }
    else{
      mc<-effect[5][[1]]
    }
    f_function_var_coef <-
      paste0(
        f_function_var_coef,
        "+ ti(",
        effect[1],
        ",",
        density_var,
        ", bs=c(",
        '"',
        effect[3][[1]],
        '"',
        ",",
        '"md")' ,
        ",m = list(",
        deparse(effect[4][[1]]),
        ",",
        deparse(m_density_var),
        "), k =c( ",
        effect[5],
        ",",
        k_density_var,
        ")",
        ",mc = c(",mc,", FALSE), np = FALSE, by=",
        effect[2],
        ",sp=array(c(-1,",
        sp_density_var_vec,
        " )),fx=",
        fx,
        ",xt=list(NULL,list(values_discrete=",
        deparse(values_discrete),
        ",  domain_continuous=",
        deparse(domain_continuous),
        ", weights_discrete=",
        deparse(weights_discrete),
        ",penalty_discrete= ",
        penalty_discrete,

        ")))"
      )
    j<-j+1
    }
  for (effect in flexible_interaction) {
    if (length(sp_density_var)>1){
      sp_density_var_<-sp_density_var[j]
      sp_density_var_vec<-sp_density_var[j]
    }
    if (is.null(effect[[5]])){
      mc<-paste0(rep("TRUE", length(effect[[1]])), collapse = ",")
    }
    else{
      mc<-paste0(effect[[5]], collapse=",")
    }
    if (is.null(effect[6][[1]])) {
    f_flex_interact <-
      paste0(
        f_flex_interact,
        "+ ti(",
        paste0(effect[[1]], collapse = ","),
        ",",
        density_var
        ,
        paste0(",bs=c('", paste0(effect[[2]], collapse = "','"), "','md')") ,
        ",m = list(",
        paste0(effect[[3]], collapse = ","),
        ",",
        deparse(m_density_var),
        "), k =c( ",
        paste(effect[[4]], collapse = ","),
        ",",
        k_density_var,
        ")",
        ",mc = c(",
        mc,
        ",FALSE), np = FALSE,sp=array(c(",
        paste0(rep("-1", length(effect[[1]])), collapse = ","),
        ",",
        sp_density_var_vec,
        " )),fx=",
        fx,
        ",xt=list(",
        paste0(rep("NULL", length(effect[[1]])), collapse = ","),
        ",list(values_discrete=",
        deparse(values_discrete),
        ",  domain_continuous=",
        deparse(domain_continuous),
        ", weights_discrete=",
        deparse(weights_discrete),
        ",penalty_discrete=",
        penalty_discrete,
        ")))"
      )
    j<-j+1}
    else{
      f_flex_interact <-
        paste0(
          f_flex_interact,
          "+ ti(",
          paste0(effect[[1]], collapse = ","),
          ",",
          density_var
          ,
          paste0(",bs=c('", paste0(effect[[2]], collapse = "','"), "','md')") ,
          ",m = list(",
          paste0(effect[[3]], collapse = ","),
          ",",
          deparse(m_density_var),
          "), k =c( ",
          paste(effect[[4]], collapse = ","),
          ",",
          k_density_var,
          ")",
          ",mc = c(",
          mc,
          ",FALSE), np = FALSE,sp=array(c(",
          paste0(rep("-1", length(effect[[1]])), collapse = ","),
          ",",
          sp_density_var_vec,
          " )),fx=",
          fx,
          ", by=",
          effect[6],",xt=list(",
          paste0(rep("NULL", length(effect[[1]])), collapse = ","),
          ",list(values_discrete=",
          deparse(values_discrete),
          ",  domain_continuous=",
          deparse(domain_continuous),
          ", weights_discrete=",
          deparse(weights_discrete),
          ",penalty_discrete=",
          penalty_discrete,
          ")))"
        )
      j<-j+1

    }
  }
  f <-
    paste0(
      f_main,
      f_intercepts,
      f_flexibles,
      f_linear,
      f_function_var_coef,
      f_flex_interact,
      "+ as.factor(group_id) -1 + offset(log(Delta) + gam_offsets)"
    )
  if (!is.null(sample_weights)) {
    m <-
      do.call(
        "gam",
        list(
          formula = as.formula((f)),
          data = as.name("dta_est"),
          method = "REML",
          family = poisson(),
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
          formula = as.formula((f)),
          data = as.name("dta_est"),
          method = "REML",
          family = poisson(),
          ...
        )
      )
  }
  X <- model.matrix(m)
  # remove intercepts per covariate combination (not of interest for estimated densities)
  intercepts <- which(grepl("group_id", colnames(X)))
  X <- X[,-intercepts]
  theta_hat <- m$coefficients[-intercepts]
  f_hat_clr <- X %*% theta_hat
  if (is.numeric(density_var)) {
    densi <- colnames(dta)[density_var]
  } else{
    densi <- density_var
  }
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
  if(colnames(dta_est)[2]=="weighted_counts"){
    t <- unique(dta_est[,3])}
  else{
    t<-unique(dta_est[,2])
  }
  if (!isFALSE(values_discrete) & !isFALSE(domain_continuous)) {
    # interval_width <-
    #   (domain_continuous[2] - domain_continuous[1]) / bin_number
    # t <- seq(
    #   from = domain_continuous[1] + 1 / 2 * interval_width,
    #   to = domain_continuous[2] - 1 / 2 * interval_width,
    #   by = interval_width
    # )

    if (is.null(bin_number)){
      bin_number<-nrow(t)-length(values_discrete)
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
  {    if (is.null(bin_number)){
    bin_number<-nrow(t)
  }
  }
  if (!effects) {
    result <-    list(
      model = m,
      model_matrix = X,
      theta_hat = theta_hat,
      f_hat_clr = f_hat_clr,
      f_hat = f_hat,
      histo_data = dta_est,
      params = list(domain_continuous, values_discrete, G = bin_number),
      predicted_effects = list(
        group_specific_intercepts = group_specific_intercepts,
        flexible_effects = flexible_effects,
        linear_effects = linear_effects,
        varying_coefficients = varying_coefficients,
        flexible_interaction = flexible_interaction
      ),
      ID_covCombi = cov_combi_id
    )
    attr(result, "class") <- c("dens_reg_obj", class(result))
    return(result)
  } else{
    pred_terms = predict(m, type = "terms")
    G <- nrow(obs_density)
    n_singles <- length(group_specific_intercepts)
    levels_singles <- c()
    dta<-as.data.table(dta)
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
    }
    else{
      positions_smooths <- c()
    }
    smooth_cols <- c()
    for (effect in flexible_effects) {
      #
      if (is.null(effect[5][[1]])) {
        ind <- match(effect[1][[1]], colnames(dta_est))
        smooth_cols <- append(smooth_cols, ind)
      } else{
        ind <- match(c(effect[1][[1]], effect[5][[1]]), colnames(dta_est))
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
    for (effect in flexible_interaction) {
      number_covs <- length(effect[[1]])
      ind <- match(effect[[1]], colnames(dta_est))
      ind <- c(ind, "inter")
      smooth_cols <- append(smooth_cols, list(ind))
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
        model = m,
        model_matrix = X,
        theta_hat = theta_hat,
        f_hat_clr = f_hat_clr,
        f_hat = f_hat,
        histo_data = dta_est,
        effects = effects,
        params = list(domain_continuous, values_discrete, G = bin_number),
        predicted_effects = list(
          group_specific_intercepts = group_specific_intercepts,
          flexible_effects = flexible_effects,
          linear_effects = linear_effects,
          varying_coefficients = varying_coefficients,
          flexible_interaction = flexible_interaction
        ),
        ID_covCombi = cov_combi_id
      )
    attr(result, "class") <- c("dens_reg_obj", class(result))
    return(result)
  }
}



## Check arguments of dens_reg (preprocess checks remaining argumets)
checking_dens_reg <-
  function(m_density_var,
           k_density_var,
           sp_density_var,
           penalty_discrete,
           group_specific_intercepts,
           linear_effects,
           flexible_effects,
           varying_coefficients,
           flexible_interaction,
           effects,
           dta) {
    # check if x is natural number
    check_natural_number <- function(x) {
      return(is.numeric(x) && length(x) == 1 && x == floor(x) && x > 0)
    }

    check_natural_number_vec <- function(x) {
      return(is.numeric(x) && all(x == floor(x)) && all(x > 0))
    }


    # check if x is integer (or NULL)
    check_integer_or_null <- function(x) {
      return(is.null(x) ||
               (is.numeric(x) && length(x) == 1 && x == floor(x)))
    }

    # check gam-parameters
    if (!(length(m_density_var) == 2) &&
        !all(check_natural_number_vec(m_density_var))) {
      stop("m_density_var must be a vector of 2 natural numbers.")
    }

    if (!check_natural_number(k_density_var)) {
      stop("k_density_var must be a natural number.")
    }

    if (!is.numeric(sp_density_var)&!is.null(sp_density_var)) {
      stop("sp_density_var must be numeric or NULL.")
    }

    if (!check_integer_or_null(penalty_discrete)) {
      stop("penalty_discrete must be an integer or NULL.")
    }

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

    # check function for structure of flexible_effects
    check_flexible_effects_structure <- function(effects) {
      if (!is.list(effects)) {
        stop("flexible_effects must be a list of lists.")
      }

      for (effect in effects) {
        if (length(effect) < 5 || length(effect) > 6) {
          stop("Each flexible_effects list must have 4 or 5 elements.")
        }

        if (!is.character(effect[[1]]) ||
            !is.numeric(dta[[effect[[1]]]])) {
          stop(
            paste(
              "First element of each flexible_effect must be a variable contained in dta and the corresponding variable in dta must be numeric."
            )
          )
        }

        if (!is.character(effect[[2]])) {
          stop(
            "Second element of each flexible_effect must be a character which specifies a basis (see bs in ?mgcv::ti)."
          )
        }

        if (!check_natural_number_vec(effect[[3]]) ||
            !check_natural_number(effect[[4]])) {
          stop(
            "The third and fourth elements of each flexible_effect must be natural numbers specifying m and k for the covariate direction."
          )
        }

        if (length(effect) == 5 &&
            !is.null(effect[[5]]) &&
            (!is.character(effect[[5]]) ||
             !is.factor(dta[[effect[[5]]]]))) {
          stop(
            paste(
              "The fifth element of each flexible_effect must be a variable contained in dta and the corresponding variable must be factor."
            )
          )
        }
      }
    }

    # check flexible_effects
    if (!is.null(flexible_effects)) {
      check_flexible_effects_structure(flexible_effects)
    }

    # check function for structure of varying_coefficients
    check_varying_coefficients_structure <-
      function(effects) {
        if (!is.list(effects)) {
          stop("varying_coefficients must be a list of lists.")
        }

        for (effect in effects) {
          if (length(effect) != 6) {
            stop("Each varying_coefficients list must have 6 elements.")
          }

          if (!is.character(effect[[1]]) ||
              !is.numeric(dta[[effect[[1]]]])) {
            stop(
              paste(
                "First element of each varying_coefficient must be a variable contained in dta and the corresponding variable must be numeric."
              )
            )
          }

          if (!is.character(effect[[2]]) ||
              !is.numeric(dta[[effect[[2]]]])) {
            stop(
              paste(
                "Second element of each varying_coefficient must be a variable contained in dta and the corresponding variable must be numeric."
              )
            )
          }

          if (!is.character(effect[[3]])) {
            stop(
              "Third element of each varying_coefficient must be a character which specifies a basis (see bs in ?mgcv::ti)."
            )
          }

          if (!check_natural_number_vec(effect[[4]]) ||
              !check_natural_number(effect[[5]])) {
            stop(
              "The fourth and fifth elements of each varying_coefficient must be natural numbers specifying m and k for the covariate direction."
            )
          }
        }
      }

    # check varying_coefficients
    if (!is.null(varying_coefficients)) {
      check_varying_coefficients_structure(varying_coefficients)
    }

    # check function for structure of flexible_interaction
    check_flexible_interaction_structure <- function(effects) {
      if (!is.list(effects)) {
        stop("flexible_interaction must be a list of lists.")
      }

      for (effect in effects) {
        if (!(
          length(effect[[1]]) == length(effect[[2]]) &&
          length(effect[[3]]) == length(effect[[4]]) &&
          length(effect[[1]]) == length(effect[[4]])
        )) {
          stop("Unequal number of elements in the lists for a flexible interaction effect.")
        }
        if (!is.character(effect[[1]]) ||
            !is.numeric(unlist(dta[, effect[[1]]]))) {
          stop(
            paste(
              "All elements in the first list of each flexible_interaction must be variables contained in dta  and the corresponding variable must be numeric."
            )
          )
        }

        if (!is.character(effect[[2]])) {
          stop(
            "The second list of each flexible_interaction must consist of characters which specify a respective basis (see bs in ?mgcv::ti)."
          )
        }

        if (!check_natural_number_vec(unlist(effect[[3]])) ||
            !is.list(effect[[3]]) ||
            !check_natural_number_vec(unlist(effect[[4]]))) {
          stop(
            "All elements of the third and fourth list of each flexible_interaction must be natural numbers specifying m and k for all covariate directions."
          )
        }
      }
    }

    # check flexible_interaction
    if (!is.null(flexible_interaction)) {
      check_flexible_interaction_structure(flexible_interaction)
    }

    # check effects
    if (!is.logical(effects) || length(effects) != 1) {
      stop("effects must be TRUE or FALSE.")
    }

    return(TRUE)
  }




#' Plot function for conditional density regression models
#'
#' \code{plot.dens_reg_obj} is the default plot method for data of the class \code{dens_reg_obj}.
#' @encoding UTF-8
#' @param obj \code{dens_reg_obj}-object, i.e. the output of the \code{dens_reg}-function.
#' @param type "histo" or "effects": If \code{type="histo"}, the underlying histogram and the estimated conditional density is plotted for each unique covariate combination. If \code{type="effects"} the different partial effects are plotted.
#' @param interactive Indicates for \code{type="effects"} wether the plots are displayed as an interactive plot (\code{TRUE}) or if all plots are plotted consecutively.
#' @param level Indicates for \code{type="effects"} with "clr" that the clr-transformed effects in \eqn{L^2_0(\mu)} and with "pdf" that the partial effects in \eqn{B^2(\mu)} are plotted.
#' @param display_all Indicates for \code{type="effects"} and \code{predict=NULL} wether all curves are plotted for effects depending on continuous covariates, i.e., for each observed value of the covariate one curve is plotted, (\code{display_all=TRUE}) or if only three curves (min, med, max) are shown (\code{display_all=FALSE}).
#' If the effects are predicted for new data, i.e., \code{predict} is not \code{NULL}, \code{display_all=FALSE} causes that only the predicted effects for the first, middle and last row are displayed.
#' @param pick_sites Vector of integers from 1 to 6. Indicates for \code{type="effects"} if only certain effects should be plotted. Number corresponds to the order of the different effect types in the argument of \code{dens_reg}. If missing (\code{NULL}) all included effects are plotted.
#' @param predict If plots for certain covariate values are desired, \code{predict} has to be a data frame with new data in form of a data table with columns named as the relevant covariates. In each row, the user can specify a value of the respective covariate.
#' @param terms If predict is not \code{NULL}. Vector of term names or the indices of the terms (starting with 1=intercept) specifying which terms should be predicted and plotted.
#' @return Plot(s) as specified.
#' @examples
#'
#' \donttest{# please run the examples of dens_reg to estimate the needed models
#'
#' example("dens_reg")
#'
#' #' # create new data for predict
#'
#' nd<-data.frame(covariate1=c("a","b","c","a"),covariate4=c(0.4,0.5,0.1,0.3), covariate2=c("d","d","c","d"),covariate3=c(1,0,0.2,2),covariate5=c(0.2,0.4,1,2))
#'
#'
#' # plot mixed model (default settings: histogram, not interactive, density-level, display all)
#'
#' plot(m_mixed)
#'
#' # plot partial effects on clr-level of the continuous model in an interactive plot, do not show all groups
#'
#' plot(m_cont, type="effects", interactive=TRUE, level="clr", display_all=FALSE)
#'
#' # show only first plot
#'
#' plot(m_cont, type="effects", interactive=FALSE, pick_sites=1, level="clr", display_all=FALSE)
#'
#' # plot partial effects on density-level estimated for new data based on the mixed model
#'
#' plot(m_mixed, type="effects", level="pdf", display_all=TRUE,predict=nd)
#'
#' # estimate and plot only the intercept (first term)
#'
#' plot(m_mixed, type="effects", level="pdf", display_all=TRUE,predict=nd, terms=1)
#'
#' #' # estimate and plot only second term
#'
#' plot(m_mixed, type="effects", level="pdf", display_all=TRUE,predict=nd, terms=2)
#' }
#'
#' @export
#'
plot.dens_reg_obj <-
  function(obj,
           type = "histo",
           interactive = FALSE,
           level = "pdf",
           display_all = TRUE,
           pick_sites = FALSE,
           predict = NULL,
           terms = NULL,
           ...) {
    #source("plain_fcts.R")
    if (type == "histo") {
      if (any(obj$histo_data$discrete == TRUE) &
          any(obj$histo_data$discrete == FALSE)) {
        interactive_histo_and_dens_mixed(
          obj$histo_data,
          obj$f_hat,
          domain = obj$params[[1]],
          G = obj$params[[3]],
          discretes = obj$params[[2]]
        )
      }
      if (all(obj$histo_data$discrete == FALSE)) {
        interactive_histo_and_dens(
          obj$histo_data,
          obj$f_hat,
          domain = obj$params[[1]],
          G = obj$params[[3]],
          ...
        )
      }
      if (all(obj$histo_data$discrete == TRUE)) {
        interactive_histo_and_dens_discrete(obj$histo_data,
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
            predict(
              obj,
              type = "terms",
              which = terms,
              new_data = predict,
              level = level
            )
          j <- 1
          for (eff in effects) {

            sp <- "[(),:]"

            parts <- unlist(strsplit(names(effects)[j], sp))
            included<-sapply(colnames(predict), function(substr) {
              any(grepl(substr,parts, ignore.case = TRUE))
            })
            parts<-colnames(predict)[included]



            #parts<-parts[parts%in%(colnames(predict))]
            terms_title<-paste(parts,collapse=", ")
            legend_labels<-paste(as.data.frame(unlist(t((predict%>%select(parts))))))
            #legend_labels<-gsub('["c(")]', '', legend_labels)
            if( names(effects)[j]=="intercept"){
              eff<-eff[,1]
              single<-TRUE
              terms_title<-"intercept"
              legend.position="none"
              legend_labels<-"intercept"
            }else{
              if(display_all==FALSE){
                if(length(legend_labels)==1){
                indx<-c(1)}
                if(length(legend_labels)==2){
                  indx<-c(1,length(legend_labels))}
                if(length(legend_labels)>2){
                  indx<-c(1,round(length(legend_labels)/2))}

              }else{
                indx<-1:length(legend_labels)
              }
              single<-TRUE
              legend.position=NULL
              eff<-eff[,indx]
              legend_labels[indx]
            }



            if (all(obj$histo_data$discrete == FALSE)) {

              p <-
                plot_densities(
                  eff,
                  G = G,
                  domain = domain,
                  #single = TRUE,
                  main = names(effects)[j],
                  legend_title = terms_title,
                  ylab = "density",
                  single=single,
                  legend.position=legend.position,
                  legend_names = legend_labels
                )
              plot_list <- append(plot_list, list(p))
              j <- j + 1}

            if (all(obj$histo_data$discrete == TRUE)) {
              p <-
                plot_densities_discrete(
                  eff,
                  G = G,
                  values_discrete = unlist(obj$params[2]),
                 # single = TRUE,
                  main = names(effects)[j],
                  legend_title = terms_title,
                  ylab = "density",
                  single=single,
                 legend.position=legend.position,
                 legend_names = legend_labels
                )
              plot_list <- append(plot_list, list(p))
              j <- j + 1
            }
            if (any(obj$histo_data$discrete == FALSE) &
                any(obj$histo_data$discrete == TRUE)) {
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
                  single=single,
                 legend.position=legend.position,
                 legend_names = legend_labels
                )
              plot_list <- append(plot_list, list(p))
              j <- j + 1
            }

          }
        }
        else{
          effects <-
            predict(
              obj,
              type = "terms",
              which = terms,
              new_data = predict
            )
          j <- 1
          for (eff in effects) {
            if (all(obj$histo_data$discrete == FALSE)) {
              lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
              p <-
                plot_densities(
                  eff,
                  G = G,
                  domain = domain,
                  single = TRUE,
                  main = names(effects)[j],
                  legend_title = terms_title,
                  legend_names="",
                  ylab = "clr(density)"
                )
              plot_list <- append(plot_list, list(p))
              j <- j + 1
            }
            if (all(obj$histo_data$discrete == TRUE)) {
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
            if (any(obj$histo_data$discrete == FALSE) &
                any(obj$histo_data$discrete == TRUE)) {
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
            min = 1, max = length(plot_list)
          )))
        }


      } else{
        if (all(obj$histo_data$discrete == FALSE)) {
          if (level == "pdf") {
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  num_lev <- length(lev)
                  p <-
                    plot_densities(
                      obj$effects[["pdf_estimated_effects"]][[j]],
                      G = G,
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0("Flexible effect of ", flexible_effects[[1]]),
                      legend_title = paste0(flexible_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev_by <- unique(obj$histo_data[[flexible_effects[[5]]]])

                  num_lev_by <- length(lev_by)
                  for (i in c(1:(num_lev_by - 1))) {
                    lev <-
                      obj$histo_data %>% select(flexible_effects[[5]], flexible_effects[[1]])
                    lev <- lev %>% filter(.[[1]] == lev_by[i + 1])
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
                          "Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1]
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]]
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "Flexible interaction of ",
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  lev_used <- c(min(lev), median(lev), max(lev))
                  used_cols <- c(1, median(c(1:length(
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
                      ylab = paste0("Flexible effect of ", flexible_effects[[1]]),
                      legend_title = paste0(flexible_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev_by <- sort(unique(obj$histo_data[[flexible_effects[[5]]]]))

                  num_lev_by <- length(lev_by)
                  for (i in c(1:(num_lev_by - 1))) {
                    lev <-
                      obj$histo_data %>% select(flexible_effects[[5]], flexible_effects[[1]])
                    lev <- lev %>% filter(.[[1]] == lev_by[i + 1])
                    lev <- sort(unname(unlist(unique(
                      lev[, 2]
                    ))))
                    num_lev <- length(lev)
                    lev_used <- c(min(lev), median(lev), max(lev))
                    used_cols <- c(1, median(c(1:length(
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
                          "Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1]
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
                num_lev <- length(lev)
                lev_used <- c(min(lev), median(lev), max(lev))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]]
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "Flexible interaction of ",
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
          if (level == "clr") {
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  num_lev <- length(lev)
                  p <-
                    plot_densities(
                      obj$effects[["estimated_effects"]][[j]],
                      G = G,
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0("clr(Flexible effect of ", flexible_effects[[1]]),
                      legend_title = paste0(flexible_effects[[1]], ")"),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  lev_by <-
                    unique(obj$histo_data[[flexible_effects[[5]]]])
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
                          "clr(Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1],
                          ")"
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]],
                      ")"
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }

              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "clr(Flexible interaction of ",
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  lev_used <- c(min(lev), median(lev), max(lev))
                  used_cols <- c(1, median(c(1:length(
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
                        "clr(Flexible effect of ",
                        flexible_effects[[1]],
                        ")"
                      ),
                      legend_title = paste0(flexible_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev_by <- sort(unique(obj$histo_data[[flexible_effects[[5]]]]))

                  num_lev_by <- length(lev_by)
                  for (i in c(1:(num_lev_by - 1))) {
                    lev <-
                      obj$histo_data %>% select(flexible_effects[[5]], flexible_effects[[1]])
                    lev <- lev %>% filter(.[[1]] == lev_by[i + 1])
                    lev <- sort(unname(unlist(unique(
                      lev[, 2]
                    ))))
                    num_lev <- length(lev)
                    lev_used <- c(min(lev), median(lev), max(lev))
                    used_cols <- c(1, median(c(1:length(
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
                          "clr(Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1],
                          ")"
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
                num_lev <- length(lev)
                lev_used <- c(min(lev), median(lev), max(lev))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]],
                      ")"
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }

              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "clr(Flexible interaction of ",
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

        if (all(obj$histo_data$discrete == TRUE)) {
          if (level == "pdf") {
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  num_lev <- length(lev)
                  p <-
                    plot_densities_discrete(
                      obj$effects[["pdf_estimated_effects"]][[j]],
                      G = G,
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0("Flexible effect of ", flexible_effects[[1]]),
                      legend_title = paste0(flexible_effects[[1]]),
                      values_discrete = obj$params[[2]],
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev_by <- unique(obj$histo_data[[flexible_effects[[5]]]])

                  num_lev_by <- length(lev_by)
                  for (i in c(1:(num_lev_by - 1))) {
                    lev <-
                      obj$histo_data %>% select(flexible_effects[[5]], flexible_effects[[1]])
                    lev <- lev %>% filter(.[[1]] == lev_by[i + 1])
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
                          "Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1]
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]]
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "Flexible interaction of ",
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  lev_used <- c(min(lev), median(lev), max(lev))
                  used_cols <- c(1, median(c(1:length(
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
                      ylab = paste0("Flexible effect of ", flexible_effects[[1]]),
                      legend_title = paste0(flexible_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev_by <- sort(unique(obj$histo_data[[flexible_effects[[5]]]]))

                  num_lev_by <- length(lev_by)
                  for (i in c(1:(num_lev_by - 1))) {
                    lev <-
                      obj$histo_data %>% select(flexible_effects[[5]], flexible_effects[[1]])
                    lev <- lev %>% filter(.[[1]] == lev_by[i + 1])
                    lev <- sort(unname(unlist(unique(
                      lev[, 2]
                    ))))
                    num_lev <- length(lev)
                    lev_used <- c(min(lev), median(lev), max(lev))
                    used_cols <- c(1, median(c(1:length(
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
                          "Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1]
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
                num_lev <- length(lev)
                lev_used <- c(min(lev), median(lev), max(lev))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]]
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "Flexible interaction of ",
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
          if (level == "clr") {
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  num_lev <- length(lev)
                  p <-
                    plot_densities_discrete(
                      obj$effects[["estimated_effects"]][[j]],
                      G = G,
                      values_discrete = obj$params[[2]],
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0("clr(Flexible effect of ", flexible_effects[[1]]),
                      legend_title = paste0(flexible_effects[[1]], ")"),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  lev_by <-
                    unique(obj$histo_data[[flexible_effects[[5]]]])
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
                          "clr(Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1],
                          ")"
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]],
                      ")"
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }

              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "clr(Flexible interaction of ",
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  lev_used <- c(min(lev), median(lev), max(lev))
                  used_cols <- c(1, median(c(1:length(
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
                        "clr(Flexible effect of ",
                        flexible_effects[[1]],
                        ")"
                      ),
                      legend_title = paste0(flexible_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev_by <- sort(unique(obj$histo_data[[flexible_effects[[5]]]]))

                  num_lev_by <- length(lev_by)
                  for (i in c(1:(num_lev_by - 1))) {
                    lev <-
                      obj$histo_data %>% select(flexible_effects[[5]], flexible_effects[[1]])
                    lev <- lev %>% filter(.[[1]] == lev_by[i + 1])
                    lev <- sort(unname(unlist(unique(
                      lev[, 2]
                    ))))
                    num_lev <- length(lev)
                    lev_used <- c(min(lev), median(lev), max(lev))
                    used_cols <- c(1, median(c(1:length(
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
                          "clr(Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1],
                          ")"
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
                num_lev <- length(lev)
                lev_used <- c(min(lev), median(lev), max(lev))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]],
                      ")"
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }

              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "clr(Flexible interaction of ",
                      paste0(unlist(flexi[[1]]), collapse = ", "),
                      ","
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
        if (any(obj$histo_data$discrete == TRUE) &&
            any(obj$histo_data$discrete == FALSE)) {
          if (level == "pdf") {
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  num_lev <- length(lev)
                  p <-
                    plot_densities_mixed(
                      obj$effects[["pdf_estimated_effects"]][[j]],
                      G = G,
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0("Flexible effect of ", flexible_effects[[1]]),
                      legend_title = paste0(flexible_effects[[1]]),
                      values_discrete = obj$params[[2]],
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev_by <- unique(obj$histo_data[[flexible_effects[[5]]]])

                  num_lev_by <- length(lev_by)
                  for (i in c(1:(num_lev_by - 1))) {
                    lev <-
                      obj$histo_data %>% select(flexible_effects[[5]], flexible_effects[[1]])
                    lev <- lev %>% filter(.[[1]] == lev_by[i + 1])
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
                          "Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1]
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]]
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "Flexible interaction of ",
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  lev_used <- c(min(lev), median(lev), max(lev))
                  used_cols <- c(1, median(c(1:length(
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
                      ylab = paste0("Flexible effect of ", flexible_effects[[1]]),
                      legend_title = paste0(flexible_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev_by <- sort(unique(obj$histo_data[[flexible_effects[[5]]]]))

                  num_lev_by <- length(lev_by)
                  for (i in c(1:(num_lev_by - 1))) {
                    lev <-
                      obj$histo_data %>% select(flexible_effects[[5]], flexible_effects[[1]])
                    lev <- lev %>% filter(.[[1]] == lev_by[i + 1])
                    lev <- sort(unname(unlist(unique(
                      lev[, 2]
                    ))))
                    num_lev <- length(lev)
                    lev_used <- c(min(lev), median(lev), max(lev))
                    used_cols <- c(1, median(c(1:length(
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
                          "Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1]
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
                num_lev <- length(lev)
                lev_used <- c(min(lev), median(lev), max(lev))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]]
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }
              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "Flexible interaction of ",
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
          if (level == "clr") {
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  num_lev <- length(lev)
                  p <-
                    plot_densities_mixed(
                      obj$effects[["estimated_effects"]][[j]],
                      G = G,
                      values_discrete = obj$params[[2]],
                      single = TRUE,
                      legend_names = round(lev, digits = 3),
                      ylab = paste0("clr(Flexible effect of ", flexible_effects[[1]]),
                      legend_title = paste0(flexible_effects[[1]], ")"),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  lev_by <-
                    unique(obj$histo_data[[flexible_effects[[5]]]])
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
                          "clr(Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1],
                          ")"
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]],
                      ")"
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }

              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "clr(Flexible interaction of ",
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
                lev <- sort(unique(obj$histo_data[[group_spec]]))
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
              for (flexible_effects in obj$predicted_effects[[2]]) {
                if (is.null(flexible_effects[[5]])) {
                  lev <- sort(unique(obj$histo_data[[flexible_effects[[1]]]]))
                  lev_used <- c(min(lev), median(lev), max(lev))
                  used_cols <- c(1, median(c(1:length(
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
                        "clr(Flexible effect of ",
                        flexible_effects[[1]],
                        ")"
                      ),
                      legend_title = paste0(flexible_effects[[1]]),
                      ...
                    )
                  j <- j + 1
                  plot_list <- append(plot_list, list(p))
                }
                else{
                  lev_by <- sort(unique(obj$histo_data[[flexible_effects[[5]]]]))

                  num_lev_by <- length(lev_by)
                  for (i in c(1:(num_lev_by - 1))) {
                    lev <-
                      obj$histo_data %>% select(flexible_effects[[5]], flexible_effects[[1]])
                    lev <- lev %>% filter(.[[1]] == lev_by[i + 1])
                    lev <- sort(unname(unlist(unique(
                      lev[, 2]
                    ))))
                    num_lev <- length(lev)
                    lev_used <- c(min(lev), median(lev), max(lev))
                    used_cols <- c(1, median(c(1:length(
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
                          "clr(Flexible effect of ",
                          flexible_effects[[1]],
                          " given ",
                          flexible_effects[[5]],
                          "=",
                          lev_by[i + 1],
                          ")"
                        ),
                        legend_title = paste0(flexible_effects[[1]]),
                        ...
                      )
                    j <- j + 1
                    plot_list <- append(plot_list, list(p))
                  }

                }
              }

              for (linear in obj$predicted_effects[[3]]) {
                lev <- sort(unique(obj$histo_data[[linear]]))
                num_lev <- length(lev)
                lev_used <- c(min(lev), median(lev), max(lev))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
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

              for (fvc in obj$predicted_effects[[4]]) {
                lev_base <-
                  paste0(
                    "(",
                    round(obj$histo_data[[fvc[[1]]]], digits = 3),
                    ",",
                    round(obj$histo_data[[fvc[[2]]]], digits = 3),
                    ")"
                  )
                lev <- sort(unique(lev_base))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))
                lev_by <- unique(obj$histo_data[[fvc[[2]]]])
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
                      fvc[[1]],
                      " given ",
                      fvc[[2]],
                      ")"
                    ),
                    legend_title = paste0("(", fvc[[1]], ",", fvc[[2]], ")"),
                    ...
                  )
                j <- j + 1
                plot_list <- append(plot_list, list(p))
              }

              for (flexi in obj$predicted_effects[[5]]) {
                lev_base <- round(obj$histo_data %>% select(flexi[[1]]), digits = 3)

                lev <- unite(lev_base,
                             col = 'com',
                             flexi[[1]],
                             sep = ',')

                lev <- sort(unique(lev$com))
                used_cols <-
                  c(1, median(c(1:length(lev))), length(lev))

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
                      "clr(Flexible interaction of ",
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
            min = 1, max = length(plot_list)
          )))
        }
      }
    }

  }





#' Predict function for conditonal density regression models
#'
#' \code{predict.dens_reg_obj} is the default predict method for data of the class \code{dens_reg_obj}.

#' @param obj \code{dens_reg_obj}-object, i.e. the output of the \code{dens_reg}-function.
#' @param new_data New data in form of a data table with columns named as the relevant covariates for the terms which should be predicted. In each row, the user can specify a value of the respective covariate. If not specyfied (\code{NULL}) the values of the original data set (see obj$histo_data) are used for the predictions.
#' @param which Only terms (or their index number) named in this array will be predicted. Covariates only needed for these terms have to be given in \code{new_data}. If \code{which} and \code{exclude} are given, \code{exclude} will be used.
#' @param exclude Terms (or their index number) named in this array will not be predicted. Covariates only needed for terms which are excluded do not have to be given in \code{new_data}. If \code{which} and \code{exclude} are given, \code{exclude} will be used.
#' @param type "terms", "pdf" or "clr": If \code{type="terms"} each component of the linear predictor is returned seperately, if
#' @param level If \code{type="terms"}: "pdf" results in predicted terms on pdf-level, "clr" in terms on clr-level.  \code{type="pdf"}  or \code{type="clr"} specifies the prediction of \eqn{\hat f} or \eqn{clr(\hat f)} of the respective covariate values.
#'
#' @return A list of matrices (if type="terms") with one matrix for each term, different columns for every predicted covariate value.
#'A matrix with columns for the different covariate combinations if \code{type="pdf"} or \code{="clr"} containing the estimated \eqn{\hat f} or \eqn{clr(\hat f)}.
#' @examples
#'
#' \donttest{# please run the examples of dens_reg to estimate the needed models
#'
#' example("dens_reg")
#'
#' # see names of the models' terms
#' sapply(m_mixed$model$smooth, "[[",  "label")
#' sapply(m_dis$model$smooth, "[[",  "label")
#' sapply(m_cont$model$smooth, "[[",  "label")
#'
#' # create new data for predict
#'
#' nd<-data.frame(covariate1=c("a","b","c","a"),covariate4=c(0.4,0.5,0.1,0.3), covariate2=c("d","d","c","d"),covariate3=c(1,0,0.2,2),covariate5=c(0.2,0.4,1,2))
#'
#' # predict mixed model, all terms on pdf-level without new data
#' p1<-predict(m_mixed, type= "terms",level="pdf")
#'
#' # predict  partial effects at clr-level of the mixed model for new data
#' p2<-predict(m_mixed, type= "terms",  new_data=nd)
#'
#' # predict only the second partial effect at density-level of the mixed model for new data
#' p3<-predict(m_mixed, type= "terms", which= "ti(obs_density):covariate1b", new_data=nd, level="pdf")
#'
#' # predict partial effects of the mixed model for new data and at density-level without the second term
#' p4<-predict(m_mixed, type= "terms", exclude= "ti(obs_density):covariate1b", new_data=nd, level="pdf")
#'
#' # predict f_hat on density-level for new data
#' p5<-predict(m_mixed, type= "pdf", new_data=nd)
#'
#' # predict clr(f_hat) for new data
#' p6<-predict(m_mixed, type= "clr", new_data=nd)
#'
#' }
#' @export
predict.dens_reg_obj <-
  function(obj,
           new_data = NULL,
           which = NULL,
           exclude = NULL,
           type = "terms",
           level = "clr") {
    if (!isFALSE(obj$params[[2]]) & !isFALSE(obj$params[[1]])) {
      n_bins <- obj$params$G + length(obj$params[[2]])
    }
    if (isFALSE(obj$params[[2]]) & !isFALSE(obj$params[[1]])) {
      n_bins <- obj$params$G
    }
    if (!isFALSE(obj$params[[2]]) & isFALSE(obj$params[[1]])) {
      n_bins <- length(obj$params[[2]])
    }
    if (is.null(exclude) & is.null(which) | !type == "terms") {
      needed <- as.list(attr(obj$model$terms, "variables"))[-c(1:5)]
    }
    if (!is.null(exclude) & !is.null(which)) {
      warning(
        "The parameters exclude and which are both specified. exclude is primarily used and may overwrite information in which"
      )
    }
    all_terms <- sapply(obj$model$smooth, "[[",  "label")
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

    if (!is.null(which) & type == "terms") {
      needed <- c()
      for (trm in which) {
        for (cov in colnames(obj$histo_data)) {
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

      nd <- obj$histo_data[rep(c(1:n_bins), nrow(new_data))]
      relevant <- match(colnames(new_data), colnames(nd))
      j <- 1
      for (index in relevant) {
        nd[, index] <- rep(new_data[, j], each = n_bins)
        j <- j + 1
      }
    }
    else{
      nd <- obj$histo_data
    }
    if (type == "terms") {
      prediction <-
        predict.gam(
          obj$model,
          newdata = nd,
          exclude = exclude,
          type = "terms",
          se.fit = FALSE
        )[, -1]
      pred_list <- list()
      prediction <- data.frame(prediction)
      for (c in 1:ncol(prediction)) {
        pr <- matrix(prediction[, c], nrow = n_bins)
        if (level == "pdf") {
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
      X <- predict(obj$model, type = "lpmatrix", newdata = nd)
      # remove intercepts per covariate combination (not of interest for estimated densities)
      intercepts <- which(grepl("group_id", colnames(X)))
      X <- X[,-intercepts]
      theta_hat <- obj$model$coefficients[-intercepts]

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
      if (type == "pdf") {
        return(f_hat)
      }
      if (type == "clr") {
        return(f_hat_clr)
      }
    }
  }
