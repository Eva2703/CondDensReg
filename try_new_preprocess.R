load("data_age_birth.RData")
load_all("C:/Users/learu/CondDensReg/new/CondDensReg")
library(dplyr)

debug(preprocess)
dta_est<-preprocess(data_age_birth, var_vec = c(1,3),
                    density_var = 2,
                    domain_continuous = c(15,50),
                    values_discrete = FALSE,
                   # bin_width = 1,
                    already_formatted = TRUE)
