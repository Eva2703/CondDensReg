load("data_age_birth.RData")
load_all("C:/Users/learu/CondDensReg/new/CondDensReg")
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
