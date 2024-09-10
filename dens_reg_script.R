
devtools::load_all("C:/Users/learu/CondDensReg/new/CondDensReg")

#  in der Beschreibung/Details wären noch mehr math. Background gut (?)
?dens_reg
?predict.dens_reg_obj
# create data (mixed)
dta <- data.frame(
  obs_density = sample(0:2, 150, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
  covariate1 = sample(c("a", "b", "c"), 150, replace = TRUE),
  covariate2 = sample(c("c", "d"), 150, replace = TRUE),
  covariate3 = rep(rnorm(n = 15), 10),
  covariate4 = rep(rnorm(n = 10), 15),
  covariate5=rep(rnorm(n = 10), 15),
  sample_weights = runif(150, 0, 2)
)
dta[which(dta$obs_density == 2),]$obs_density <-
  rbeta(length(which(dta$obs_density == 2)), shape1 = 3, shape2 = 3)
dta$covariate1 <- ordered(dta$covariate1)
dta$covariate2 <- ordered(dta$covariate2)

# create discrete data

dta_dis <- data.frame(
  obs_density = sample(0:2, 150, replace = TRUE, prob = c(0.25, 0.45,0.3)),
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

flexible_effects <-
  list(list("covariate3", "ps", c(2, 2), 4, NULL),
      list("covariate3", "ps", c(2, 2), 4, "covariate1"))
flex_inter <-
  list(list(c("covariate3", "covariate4","covariate5"), c("ps", "ps","ps"),list( c(2, 2), c(2, 2), c(2, 2)),c( 4, 4,5)))
fvc <-
  list(list("covariate3", "covariate4", "ps", c(2, 2), 4))
linear_effects <- c("covariate4")
group_specific_intercepts <- c("covariate1", "covariate2")


### fit models ### (Achtung: kann je nach Auswahl d. Effekte ein paar Minuten dauern)


devtools::load_all("C:/Users/learu/CondDensReg/new/CondDensReg")
# fit model mixed
m_mixed <- dens_reg(
  dta = dta,
  var_vec = c(2:6),
  density_var = 1,
  sample_weights = NULL,
  bin_width = NULL,
  bin_number = 100,
  values_discrete = c(0, 1),
  weights_discrete = c(1,1),
  domain_continuous = c(0,1),
  m_density_var = c(2, 2),
  k_density_var = 4,
  group_specific_intercepts = group_specific_intercepts,
  flexible_effects = NULL,
  flexible_interaction =  NULL,
  linear_effects = linear_effects,
  functional_varying_coefficients = NULL,
  effects = TRUE,
  sp_density_var=NULL,
  penalty_discrete = NULL
)


# fit model discrete
m_dis <- dens_reg(
  dta = dta_dis,
  var_vec = c(2:6),
  density_var = 1,
  sample_weights = NULL,
  bin_width = NULL,
  bin_number = 100,
  values_discrete = c(0, 1,2),
  weights_discrete = c(1,1,1),
  domain_continuous = FALSE,
  m_density_var = c(2, 2),
  k_density_var = 4,
  group_specific_intercepts = group_specific_intercepts,
  flexible_effects = NULL,
  flexible_interaction = NULL,
  linear_effects = NULL,
  functional_varying_coefficients = NULL,
  effects = TRUE,

  penalty_discrete = NULL
)



# fit model continuous

m_cont <- dens_reg(
  dta = dta%>%filter(obs_density!=0&obs_density!=1),
  var_vec = c(2:6),
  density_var = 1,
  sample_weights = NULL,
  bin_width = NULL,
  bin_number = 100,
  values_discrete = FALSE,
  weights_discrete = c(1,1),
  domain_continuous = c(0,1),
  m_density_var = c(2, 2),
  k_density_var = 12,
  group_specific_intercepts = group_specific_intercepts,
  flexible_effects = NULL,
  flexible_interaction = flex_inter,
  linear_effects = NULL,
  functional_varying_coefficients = NULL,
  effects = TRUE,
  penalty_discrete = NULL
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
plot(m_mixed, type="effects", pick_sites=c(2,3))

# customize plot parameters
plot(m_mixed, type="effects", display_all = FALSE,level="clr", main="test",legend.position="none", xlab="blub")


### predict ###
# terms names
all_terms <- sapply(m_mixed$model$smooth, "[[",  "label")

# create newdata for predict
nd<-data.frame(covariate1=c("a","b","c"),covariate4=c(0.4,0.5,0.1), covariate2=c("d","d","c"),covariate3=c(1,0,0.2),covariate5=c(0.2,0.4,1))

# predict mixed model, all terms on pdf-level without newdata
p1<-predict(m_mixed, type= "terms",level="pdf")
# new data and clr-level
p2<-predict(m_mixed, type= "terms",  new_data=nd)
# only second term
p3<-predict(m_mixed, type= "terms", which=all_terms[2], new_data=nd, level="pdf")
# without second term
p4<-predict(m_mixed, type= "terms", exclude=all_terms[2], new_data=nd, level="pdf")
# predict f_hat for new data
p5<-predict(m_mixed, type= "pdf", new_data=nd)
# predict clr(f_hat) for new data
p6<-predict(m_mixed, type= "clr", new_data=nd)
# try to predict f_hat with not all needed covariates given in new_data => ERROR
p7<-predict(m_mixed, type= "pdf", new_data=nd%>%select(covariate1))



### plot predicted effects ####

plot(m_cont, type="effects", level="pdf", predict = nd)
# only intercept and interactive effect

plot(m_cont, type="effects", level="pdf", predict = nd, terms=c(1,5))







