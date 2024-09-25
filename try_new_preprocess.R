load("data_age_birth.RData")
load_all("C:/Users/learu/CondDensReg/new/CondDensReg")
library(dplyr)

debug(preprocess)
dta_est<-preprocess(data_age_birth, var_vec = c(1,3),
                    density_var = 2,
                    domain_continuous =15:50,
                    values_discrete = 16.5,
                    bin_width = 15:49,
                    weights_discrete =c(20),
                    already_formatted = TRUE)


# mixed case nur vektor bin_width
