################################################################################
################################ Initialization ################################
################################################################################
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(current_working_dir, "/../.."))

# Package to deal with big datasets
library(data.table)
# Package to estimate GAMs (like our Poisson model)
library(mgcv)
library(devtools)

# load self-written smoother for mixed reference measure
load_all("Z:/densREG/CondDensReg")
# load functions used to prepare data appropriately to be used in Poisson model
source("Z:/densREG/data_prep.R")

# load data, add missing variables and transform characters to ordered factors
# (necessary for reference coding in gam(); lowest factor level is used as
# reference category)
# load data, add missing variables and transform characters to ordered factors
# (necessary for reference coding in gam(); lowest factor level is used as
# reference category)
dta <- as.data.table(readRDS("../income_share_data.rds"))
dta$West_East <- factor(ifelse(dta$region %in% c("northeast", "east"), "East", "West"),
                        levels = c("West", "East"), ordered = TRUE)
# dta$region <- factor(dta$region) # I don't know how to include region such that it is centered around West/East -> Need to talk to Almond
dta$c_age <- factor(dta$child_group, levels = c(3, 2, 1), ordered = TRUE)

# Create new factor variable for interaction of West_East and c_age
dta$West_East_c_age <- paste0(dta$West_East, "_", dta$c_age)
dta$West_East_c_age <- factor(dta$West_East_c_age, levels = c(paste0("West_", 3:1), paste0("East_", 3:1)),
                              ordered = TRUE)

# Normalize SOEP weight (compare to ?gam -> weights);
dta$weighting_factor_orig <- dta$weighting_factor
dta$weighting_factor <- dta$weighting_factor_orig / mean(dta$weighting_factor_orig)
# fivenum(dta$weighting_factor)


# observations with weight = 0 would have no influence, but we can't take their
# logarithm -> remove them from data set
dta <- dta[weighting_factor > 0]

# prepare data appropriately to be used in Poisson model
dta_grouped <- dta[, group_id := .GRP,
                   keyby = .(West_East, #region,
                             c_age,
                             West_East_c_age,
                             syear)][, .(share,
                                         West_East,
                                         # region,
                                         c_age,
                                         West_East_c_age,
                                         syear,
                                         weighting_factor,
                                         group_id)]

# Construct data set containing (unweighted) counts and weighted counts
step_size <- 0.01
n_bins <- 1/ step_size

dta_est <- construct_dataset(dta_grouped, response = list("share"),
                             step_size = step_size)
dta_est_weight <- construct_dataset(dta_grouped, response = list("share"),
                                    step_size = step_size, weighted = TRUE)
dta_est[, weighted_counts := dta_est_weight$counts]

# dta_est[weighted_counts == 0 & counts != 0]

# We have to specify the argument weights in gam and use an offset in the formula
# to include the SOEP-weights properly; For each bin, the weight to pass to gam
# is the weighted count of observations in this bin, divided by the usual count
# (i.e., the number) of observations in this bin; the offset that has to be used
# is -log of the respective weight
dta_est[, gam_weights := ifelse(counts != 0, weighted_counts / counts, 1)] # If counts == 0, the value of weight can actually be arbitrary, but its log has to exist
dta_est[, gam_offsets := - log(gam_weights)]

# Check whether weights are approximately normalized: mean(dta_est$gam_weights)

# # Sanity check: number of covariate combinations
# max(dta_est$group_id)

# weights for discrete points are 1, equal equidistant stepsize for all covariate
# combinations (histograms) in continuous component
dta_est[, Delta := ifelse(dta_est$share %in% 0:1, 1, step_size)]

# In East Germany, observations start only from 1991 (1984 in West). However, when
# constructing a (penalized) B-spline basis for the interaction of syear with East,
# the whole range of syear is used, starting from 1984. Thus, for East Germany,
# we need a modified syear-variable, which contains the values of syear for East
# Germany and some value larger than 1991 for West Germany. (Note that since
# West Germany is the reference, there is no effect estimated at all; Furthermore,
# for interaction with categorical covariates, each category has a separate (usually
# identical, which is however not reasonable in our case) design matrix)
dta_est$syear_ <- ifelse(dta_est$West_East == "East", dta_est$syear, mean(dta_est$syear))

# For the interaction of syear with West_East and c_age, we need two separate
# categorical variables due to the different range of years in East. Note that
# for West, we no effect should be estimated for the reference West_3 (no minor
# children), but also no effect should be estimated for East Germany, which gets
# a separate variable. Thus, we combine both to West_3_East and use it as
# reference, i.e., no effect is estimated for both.
dta_est$West_c_age <- factor(ifelse(dta_est$West_East_c_age %in% c("West_3", paste0("East_", 1:3)),
                                    "West_3_East", as.character(dta_est$West_East_c_age)),
                             levels = c("West_3_East", paste0("West_", 2:1)), ordered = TRUE)
dta_est$East_c_age <- factor(ifelse(grepl("East", dta_est$West_East_c_age),
                                    as.character(dta_est$West_East_c_age), "West"),
                             levels = c("West", paste0("East_", 3:1)), ordered = TRUE)



##########################################################################################
### change individual-level data such that syear_, east_c_age and west_c_age are included
dta_<-dta
dta_$syear_ <- ifelse(dta_$West_East == "East", dta_$syear, mean(dta_est$syear))
dta_$West_c_age <- factor(ifelse(dta_$West_East_c_age %in% c("West_3", paste0("East_", 1:3)),
                                    "West_3_East", as.character(dta_$West_East_c_age)),
                             levels = c("West_3_East", paste0("West_", 2:1)), ordered = TRUE)
dta_$East_c_age <- factor(ifelse(grepl("East", dta_$West_East_c_age),
                                    as.character(dta_$West_East_c_age), "West"),
                             levels = c("West", paste0("East_", 3:1)), ordered = TRUE)


# create preprocessed data set to compare dta_est to with preprocess
dta_vgl<-preprocess(dta_, var_vec = c(6,7,8,3,11,12,13), density_var = 1, sample_weights = 4)

dta_vgl_<-dta_vgl[,c(1,3,4,5,6,7,11,2,12,13,14,8,9,10)]

# check if both preprocessed data sets are identical
sum(dta_est[,-7]!=dta_vgl_)




## construct knots

knots_east <- smooth.construct(ti(syear, share, bs = c("ps", "md"),
                                  m = list(c(2,2), c(2,2)),
                                  k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE),
                               data = dta_est, knots = NULL)$margin[[1]]$knots

which(knots_east < 1991)

knots_east <- knots_east[-1]



# Model "by hand"
model_soep <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                              mc = FALSE, np = FALSE) # default: no penalty for discrete component
                  + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                       mc = FALSE, np = FALSE, by = West_East)
                  + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                       mc = FALSE, np = FALSE, by = c_age)
                  + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                       mc = FALSE, np = FALSE, by = West_East_c_age)
                  + ti(syear, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE)#
                  + ti(syear_, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(7, 12), mc = c(TRUE, FALSE), np = FALSE, by = West_East)#
                  + ti(syear, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE, by = c_age)#
                  + ti(syear, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE, by = West_c_age)#
                  + ti(syear_, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(7, 12), mc = c(TRUE, FALSE), np = FALSE, by = East_c_age)
                  + as.factor(group_id) - 1 # no scalar intercept, but one intercept per covariate combination
                  + offset(log(Delta) + gam_offsets),
                  data = dta_est, weights = dta_est$gam_weights,
                  knots = list(syear_ = knots_east),
                  method = "REML", family = poisson())



saveRDS(model_soep, "model_soep_standard.rds")


## model with dens_reg
load_all("Z:/densREG/CondDensReg")

model_new<- dens_reg(dta=dta_, var_vec = c(6,7,8,3,11,12,13), density_var = 1, sample_weights = 4, 
                     m_density_var = c(2, 2),
                     k_density_var = 12, 
                     group_specific_intercepts = c("West_East"),#,"c_age","West_East_c_age"),
                     flexible_effects = list(
                    list("syear","ps",c(2,2),8)),
                     #   list("syear_","ps",c(2,2),7,"West_East"),
                     #   list("syear","ps",c(2,2),8,"c_age"),
                     #   list("syear","ps",c(2,2),8,"West_c_age"),
                     #   list("syear_","ps",c(2,2),7,"East_c_age"),)
                     effects=FALSE,
                     knots = list(syear_ = knots_east))








"counts~ti(share, bs=\"md\",mc=FALSE, np=FALSE, m = list(c(2, 2)), k = 12,sp=NULL,fx=TRUE,xt=list(list(values_discrete=c(0, 1),  domain_continuous=c(0, 1), weights_discrete=c(1, 1),penalty_discrete= 1)))+ ti(share, bs=\"md\",m = list(c(2, 2)), k = 12,mc = FALSE, np = FALSE, by=West_East,sp=NULL,fx=TRUE,xt=list(list(values_discrete=c(0, 1),  domain_continuous=c(0, 1), weights_discrete=c(1, 1),penalty_discrete= 1)))+ ti(share, bs=\"md\",m = list(c(2, 2)), k = 12,mc = FALSE, np = FALSE, by=c_age,sp=NULL,fx=TRUE,xt=list(list(values_discrete=c(0, 1),  domain_continuous=c(0, 1), weights_discrete=c(1, 1),penalty_discrete= 1)))+ ti(share, bs=\"md\",m = list(c(2, 2)), k = 12,mc = FALSE, np = FALSE, by=West_East_c_age,sp=NULL,fx=TRUE,xt=list(list(values_discrete=c(0, 1),  domain_continuous=c(0, 1), weights_discrete=c(1, 1),penalty_discrete= 1)))+ ti(syear,share, bs=c(\"ps\",\"md\"),m = list(c(2, 2),c(2, 2)), k =c( 8,12),mc = c(TRUE, FALSE), np = FALSE, sp=array(c(-1,-1 )),fx=TRUE,xt=list(NULL,list(values_discrete=c(0, 1),  domain_continuous=c(0, 1), weights_discrete=c(1, 1),penalty_discrete= 1)))+ as.factor(group_id) -1 + offset(log(Delta) + gam_offsets)"




model_soep <- gam(counts ~ ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                              mc = FALSE, np = FALSE,sp=NULL) # default: no penalty for discrete component

                         + ti(share, bs = "md", m = list(c(2, 2)), k = 12,
                       mc = FALSE, np = FALSE, by = West_East_c_age)
                  + ti(syear, share, bs = c("ps", "md"), m = list(c(2,2), c(2,2)),
                       k = c(8, 12), mc = c(TRUE, FALSE), np = FALSE, sp=array(c(-1,-1)))#

                  + as.factor(group_id) - 1 # no scalar intercept, but one intercept per covariate combination
                  + offset(log(Delta) + gam_offsets),
                  data = dta_est, weights = dta_est$gam_weights,
                  knots = list(syear_ = knots_east),
                  method = "REML", family = poisson())

