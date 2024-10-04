library(devtools)
devtools::document("C:/Users/learu/CondDensReg/new/CondDensReg")
devtools::install("C:/Users/learu/CondDensReg/new/CondDensReg")
load_all("C:/Users/learu/CondDensReg/new/CondDensReg")






library(CondDensReg)
?preprocess
?dens_reg
?plot.dens_reg_obj
?predict.dens_reg_obj

example("preprocess")
example("dens_reg")
example("predict.dens_reg_obj")
example("plot.dens_reg_obj")


vignette("CondDensReg")
# for further information on the parameters of the preprocessing step
# see ?preprocess

library(dplyr)

# create data (mixed)

dta <- data.frame(obs_density = sample(0:2, 150, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
                  covariate1 = sample(c("a", "b", "c"), 150, replace = TRUE),
                  covariate2 = sample(c("c", "d"), 150, replace = TRUE),
                  covariate3 = rep(rnorm(n = 15), 10),
                  covariate4 = rep(rnorm(n = 10), 15),
                  covariate5=rep(rnorm(n = 10), 15), sample_weights = runif(150, 0, 2))
dta[which(dta$obs_density == 2),]$obs_density <-rbeta(length(which(dta$obs_density == 2)), shape1 = 3, shape2 = 3)
dta$covariate1 <- ordered(dta$covariate1)
dta$covariate2 <- ordered(dta$covariate2)

# create discrete data

dta_dis <- data.frame(obs_density = sample(0:2, 150, replace = TRUE, prob = c(0.25, 0.45,0.3)),
                      covariate1 = sample(c("a", "b", "c"), 150, replace = TRUE),
                      covariate2 = sample(c("c", "d"), 150, replace = TRUE),
                      covariate3 = rep(rnorm(n = 15), 10),
                      covariate4 = rep(rnorm(n = 10), 15),
                      covariate5=rep(rnorm(n = 10), 15),
                      sample_weights = runif(150, 0, 2))
dta_dis$covariate1 <- ordered(dta_dis$covariate1)
dta_dis$covariate2 <- ordered(dta_dis$covariate2)

# examples for different partial effects

## group specific intercepts
group_specific_intercepts <- c("covariate1", "covariate2")
## linear effects
linear_effects <- c("covariate4")
## flexible effects
flexible_effects <-list(list(cov="covariate3", bs="ps", m=c(2, 2), k=4), list(cov="covariate3", bs="ps", m=c(2, 2), k=4, mc=FALSE, by="covariate1"))
## varying coefficient
fvc <-list(list(cov="covariate3", by="covariate4", bs="ps", m=c(2, 2) , k=4))
## flexible interaction
flex_inter <-list(list(covs=c("covariate3", "covariate4","covariate5"), bs=c("ps", "ps","ps"), m=list( c(2, 2), c(2, 2), c(2, 2)), k=c( 4, 4,5), mc=list(TRUE, FALSE, TRUE), by=NULL))

# fit models (warning: calculation may take a few minutes)

## fit model for the mixed case with group specific intercepts and linear effects
### use fixed smoothing parameters in density direction and calculate also the partial effects

m_mixed <- dens_reg(dta = dta, var_vec = c(2:6), density_var = 1, m_density_var = c(2, 2),
                    k_density_var = 4, group_specific_intercepts = group_specific_intercepts,
                    linear_effects = linear_effects, effects = TRUE, sp_density_var=c(1,3,5,0.5))

