load("data_age_birth.RData")
load_all("C:/Users/learu/CondDensReg/new/CondDensReg")
document("C:/Users/learu/CondDensReg/new/CondDensReg")
library(dplyr)

debug(preprocess)

dta<-data_age_birth
dta$weighted_counts<-dta$counts*2
dta_est<-preprocess(data_age_birth, var_vec = c(1,3),
                    density_var = 2,
                    domain_continuous =FALSE,
                    values_discrete = 15:49+0.5,
                    bin_width = 15:49,
                    weights_discrete =2,
                    already_formatted = TRUE)


# mixed case nur vektor bin_width
?preprocess


dta <-data.frame(
  obs_density = sample(0:2, 100, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
  covariate1 = sample(c("a", "b"), 100, replace = TRUE),
  covariate2 = sample(c("c", "d"), 100, replace = TRUE),
  sample_weights = runif(100, 0, 2)
)
dta[which(dta$obs_density == 2), ]$obs_density <- rbeta(length(which(dta$obs_density == 2)), shape1 = 3, shape2 = 3)
p<- preprocess(
  dta,
  var_vec = c("covariate1", "covariate2"),
  density_var = "obs_density",
  sample_weights = "sample_weights",
  bin_number = 10
)

a<-preprocess(p, var_vec="covariate1", density_var = "obs_density", already_formatted = TRUE)
