% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DensityRegression-package.R
\docType{package}
\name{DensityRegression-package}
\alias{DensityRegression}
\alias{DensityRegression-package}
\alias{DensityRegression_package}
\alias{package-DensityRegression}
\title{DensityRegression: Conditional density regression for individual-level data}
\description{
Structured additive regression models for mixed (discrete/continuous)
densities as response with scalar covariates, fitted via Poisson regression
models based on \code{\link[mgcv]{gam}}, where the continuous part of the density
is approximated via histograms.
}
\details{
This package implements the approach of Maier et al. (2025b) for fitting regression
models with densities in mixed (continuous/discrete) Bayes Hilbert spaces (including
continuous and discrete ones as special cases) as response given scalar covariates, based on observed
samples of the conditional distributions via a (penalized) maximum likelihood
approach. The main function (which estimates the models) is \code{\link{densreg}}.
The (penalized) log-likelihood function is approximated via the
(penalized) log-likelihood of an appropriate Poisson regression model, which -
after constructing the count data appropriately via \code{\link{data2counts}}
- is then fitted using \code{mgcv}'s \code{\link[mgcv]{gam}} with a new smooth
term \code{bs="md"} (see \code{\link{smooth.construct.md.smooth.spec}})
for mixed densities.
% The main fitting function is \code{\link{densreg}}. Besides the (individual-
% level) data, it takes a variety of arguments related to specification of
% the underlying mixed Bayes Hilbert space, histogram construction for its
% continuous part, and partial effects (including penalties).
% We also provide \code{plot} and \code{predict} methods for the resulting fitted object.
For details see Maier et al. (2025b).
}
\references{
Maier, E.-M., Fottner, A., Stoecker, A., Okhrin, Y., & Greven, S. (2025b):
Conditional density regression for individual-level data.
arXiv preprint arXiv:XXXX.XXXXX.
}
\seealso{
\code{\link{densreg}} for the main fitting function.
}
\author{
Eva-Maria Maier, Lea Runge, and Alexander Fottner
}
