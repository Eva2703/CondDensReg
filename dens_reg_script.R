setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

devtools::document("./CondDensReg/R")
# devtools::build_manual("./CondDensReg/")
# devtools::check_man("./CondDensReg/")
# roxygen2::roxygenize("./CondDensReg")
# setwd("./CondDensReg/")
# usethis::use_description()
# usethis::use_roxygen_md()
devtools::load_all("./CondDensReg/R")
?dens_reg
?preprocess

# devtools::load_all("./preprocessing/R")

?smooth.construct.md.smooth.spec

# source("./CondDensReg/R/mixed_density_smooth.R")
# source("./CondDensReg/R/plain_fcts.R")

#  in der Beschreibung/Details wären noch mehr math. Background gut (?)
?dens_reg
?predict.dens_reg_obj

# create data (mixed)
set.seed(1335)
dta <- data.frame(
  obs_density = sample(0:2, 150, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
  covariate1 = sample(c("a", "b", "c"), 150, replace = TRUE),
  covariate2 = sample(c("c", "d"), 150, replace = TRUE),
  covariate3 = rep(rnorm(n = 15), 10),
  covariate4 = rep(rnorm(n = 10), 15),
  covariate5=rep(rnorm(n = 10), 15),
  sample_weights = runif(150, 0, 2)
)
dta[which(dta$obs_density == 2), ]$obs_density  <-
  rbeta(length(which(dta$obs_density == 2)), shape1 = 3, shape2 = 3)
dta$covariate1 <- ordered(dta$covariate1)
dta$covariate2 <- ordered(dta$covariate2)

# create discrete data

dta_dis <- data.frame(
  obs_density = sample(0:2, 150, replace = TRUE, prob = c(0.25, 0.45, 0.3)),
  covariate1 = sample(c("a", "b", "c"), 150, replace = TRUE),
  covariate2 = sample(c("c", "d"), 150, replace = TRUE),
  covariate3 = rep(rnorm(n = 15), 10),
  covariate4 = rep(rnorm(n = 10), 15),
  covariate5=rep(rnorm(n = 10), 15),
  sample_weights = runif(150, 0, 2)
)
dta_dis$covariate1 <- ordered(dta_dis$covariate1)
dta_dis$covariate2 <- ordered(dta_dis$covariate2)


# examples for different partial effects

smooth_effects <- list(list(m = c(2, 2), bs = "ps", cov = "covariate3", k = 4)
                       # , list(cov = "covariate3", bs = "ps", m = c(2, 2), k = 4, by = "covariate1")
                       )
smooth_inter <- list(list(covs = c("covariate3", "covariate4", "covariate5"),
                          bs = c("ps", "ps", "ps"),
                          m = list(c(2, 2), c(2, 2), c(2, 2)),
                          k = c(4, 4, 5)))
varying_coef <- list(list(cov = "covariate3", by = "covariate4", bs = "ps",
                          m = c(2, 2), k = 4))
linear_effects <- c("covariate4")
group_specific_intercepts <- c("covariate1", "covariate2")


### fit models ### (Achtung: kann je nach Auswahl d. Effekte ein paar Minuten dauern)

# devtools::document("./CondDensReg/R")
devtools::load_all("./CondDensReg/R")
debug(dens_reg)
debug(checking_dens_reg_2)

# devtools::load_all("C:/Users/learu/CondDensReg/new/CondDensReg")
# fit model mixed
m_mixed <- dens_reg(
  dta = dta,
  # var_vec = c(2:6),
  y = 1,
  sample_weights = NULL,
  bin_width = NULL,
  bin_number = 100,
  values_discrete = c(0, 1),
  weights_discrete = c(1, 1),
  domain_continuous = c(0, 1),
  m_continuous = c(2, 2),
  k_continuous = 4,
  group_specific_intercepts = group_specific_intercepts,
  # smooth_effects = smooth_effects,
  # smooth_interactions =  smooth_inter,
  # linear_effects = linear_effects,
  # varying_coefficients = varying_coef,
  effects = TRUE,
  # sp_y = c(1, 3, 5),
  penalty_discrete = NULL
)

# all.equal(m_mixed, m_mixed_)

# default  = histo+dens
plot(m_mixed)
## plot effects with different specifications:
# pdf
plot(m_mixed, type="effects", level="clr", interactive=TRUE)
# only min, med, max
plot(m_mixed, type="effects", display_all = FALSE)
# select specific sides
plot(m_mixed, type="effects", pick_sites=c(2, 3))

# customize plot parameters
plot(m_mixed, type="effects", display_all = FALSE, level="clr", main="test", legend.position="none", xlab="blub")


# fit model discrete
m_dis <- dens_reg(
  dta = dta_dis,
  var_vec = c(2:6),
  y = 1,
  sample_weights = NULL,
  bin_width = NULL,
  bin_number = 100,
  values_discrete = c(0, 1, 2),
  weights_discrete = c(1, 1, 1),
  domain_continuous = FALSE,
  m_continuous = c(2, 2),
  k_continuous = 4,
  group_specific_intercepts = group_specific_intercepts,
  smooth_effects = NULL,
  flexible_interaction = NULL,
  linear_effects = NULL,
  functional_varying_coefficients = NULL,
  effects = TRUE
)



# fit model continuous

m_cont <- dens_reg(
  dta = dta%>%filter(obs_density!=0&obs_density!=1),
  var_vec = c(2:6),
  y = 1,
  sample_weights = NULL,
  bin_width = NULL,
  bin_number = 100,
  values_discrete = FALSE,
  weights_discrete = c(1, 1),
  domain_continuous = c(0, 1),
  m_continuous = c(2, 2),
  k_continuous = 12,
  group_specific_intercepts = NULL,
  smooth_effects = smooth_effects,
  flexible_interaction = NULL,
  linear_effects = NULL,
  functional_varying_coefficients = NULL,
  effects = TRUE
)


### PLOT ####

# default  = histo+dens
plot(m_mixed)
plot(m_cont)
plot(m_dis)
## plot effects with different specifications:
# pdf
plot(m_mixed, type="effects")
# clr
plot(m_cont, type="effects", level="clr")
# interactive
plot(m_dis, type="effects", interactive=TRUE)
# only min, med, max
plot(m_mixed, type="effects", display_all = FALSE)
# select specific sides
plot(m_mixed, type="effects", pick_sites=c(2, 3))

# customize plot parameters
plot(m_mixed, type="effects", display_all = FALSE, level="clr", main="test", legend.position="none", xlab="blub")


### predict ###
# terms names
all_terms <- sapply(m_mixed$model$smooth, "[[",  "label")

# create newdata for predict
nd <- data.frame(covariate1=c("a", "b", "c", "a"), covariate4=c(0.4, 0.5, 0.1, 0.3), covariate2=c("d", "d", "c", "d"), covariate3=c(1, 0, 0.2, 2), covariate5=c(0.2, 0.4, 1, 2))

# predict mixed model, all terms on pdf-level without newdata
p1 <- predict(m_mixed, type= "terms", level="pdf")
# new data and clr-level
p2 <- predict(m_mixed, type= "terms",  new_data=nd)
# only second term
p3 <- predict(m_mixed, type= "terms", which=all_terms[2], new_data=nd, level="pdf")
# without second term
p4 <- predict(m_mixed, type= "terms", exclude=all_terms[2], new_data=nd, level="pdf")
# predict f_hat for new data
p5 <- predict(m_mixed, type= "pdf", new_data=nd)
# predict clr(f_hat) for new data
p6 <- predict(m_mixed, type= "clr", new_data=nd)
# try to predict f_hat with not all needed covariates given in new_data => ERROR
p7 <- predict(m_mixed, type= "pdf", new_data=nd%>%select(covariate1))



### plot predicted effects ####
devtools::load_all("C:/Users/learu/CondDensReg/new/CondDensReg")
nd <- rbind(nd, nd)
all_terms <- sapply(m_cont$model$smooth, "[[",  "label")
plot(m_cont, type="effects", level="pdf", display_all =  FALSE)

# only intercept and interactive effect

plot(m_cont, type="effects", level="pdf", predict = nd, display_all = FALSE)

debug(plot.dens_reg_obj)



nd


## display all =FALSE mit predict
## rundung +align, leerzeichen hinter komma
## sowie bei ohne predict
##








devtools::load_all("C:/Users/learu/CondDensReg/new/CondDensReg")
debug(preprocess)
dta_weighted <- data_age_birth
dta_weighted$weighted_counts <- dta_weighted$counts+9
dta_pre <- preprocess(dta=dta_weighted,
                    var_vec = c("year", "marital_status"),
                    y = "age",
                    values_discrete = FALSE,
                    domain_continuous = c(15, 50),
                    already_formatted = TRUE, bin_width = 5)
