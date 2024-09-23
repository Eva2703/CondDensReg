library(devtools)
devtools::document("C:/Users/learu/CondDensReg/new/CondDensReg")
devtools::install("C:/Users/learu/CondDensReg/new/CondDensReg")
library(CondDensReg)
?preprocess
?dens_reg
?plot.dens_reg_obj
?predict.dens_reg_obj

example("preprocess")
example("dens_reg")
example("predict.dens_reg_obj")
example("plot.dens_reg_obj")
