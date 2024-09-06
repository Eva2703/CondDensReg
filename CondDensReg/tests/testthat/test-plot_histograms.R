library(testthat)


dta <-data.frame(
  obs_density = sample(0:2, 100, replace = TRUE, prob = c(0.15, 0.1, 0.75)),
  covariate1 = sample(c("a", "b"), 100, replace = TRUE),
  covariate2 = sample(c("c", "d"), 100, replace = TRUE),
  sample_weights = runif(100, 0, 2)
)
dta[which(dta$obs_density == 2), ]$obs_density <- rbeta(length(which(dta$obs_density == 2)), shape1 = 3, shape2 = 3)
histo_data<-preprocess(
  dta,
  var_vec = c(2, 3),
  density_var = 1,
  sample_weights = 4,
  bin_width = 0.1
)
histo_data_unweighted<-preprocess(
  dta,
  var_vec = c(2, 3),
  density_var = 1,
  bin_width = 0.1)

histo_conti<-histo_data%>%filter(discrete==FALSE)
histo_disc<-histo_data%>%filter(discrete==TRUE)
# save_plots
test_that("save_plots", {
 
 expect_no_condition(save_plots(histo_data))
  expect_no_condition(save_plots(histo_data, single=FALSE,abs=TRUE,hist=FALSE))
  expect_no_condition(save_plots(histo_data, single=FALSE,abs=TRUE,hist=FALSE, file_names = "all__"))
  expect_no_condition(save_plots(histo_data, file_names = c("1","3", "4","x")))
  expect_no_condition(save_plots(histo_data, plot_width = 100, plot_height=100, file_names_stem = "changed_format"))
  ## warning: creates new directories
  #dir_name<-paste("./",paste(sample(9, replace = TRUE, 6),collapse=""),"/", sep = "")
  #expect_no_condition(save_plots(histo_data, path=dir_name))
  expect_error(save_plots(histo_data, path=2))
  expect_error(save_plots(histo_data, selected_groups = "a"))
  expect_error(save_plots(histo_data, selected_groups = c(100)))
  expect_error(save_plots(histo_data, single = "a"))
  expect_error(save_plots(histo_data, file_names = 1))
  expect_error(save_plots(histo_data, file_names_stem  = 1))
  expect_error(save_plots(histo_data, file_names  = "1", file_names_stem = "1"))
  expect_error(save_plots(histo_data, plot_width = "a"))
  expect_error(save_plots(histo_data, plot_height = "a"))
  })


# interactiive_plots and plot
test_that("plot and interactive_plots", {
  
  expect_no_condition(plot(histo_data))
  expect_no_condition(plot(histo_data%>%filter(discrete==FALSE)))
  expect_error(interactive_plots(histo_data, selected_groups = "a"))
  expect_error(interactive_plots(histo_data, selected_groups = c(100)))
    })

histo_data2<-histo_data%>%filter(group_id==1)

histo_data_unweighted2<-histo_data_unweighted%>%filter(group_id==1)
test_that("plot_single_histo",{
  expect_no_condition(plot_single_histo(histo_data2))
  expect_no_condition(plot_single_histo(histo_data_unweighted2))
  expect_no_condition(plot_single_histo(histo_data2, case="discrete"))
  expect_no_condition(plot_single_histo(histo_data2, case="continuous"))
  expect_no_condition(plot_single_histo(histo_data2, hist=FALSE))
  expect_no_condition(plot_single_histo(histo_data2, case="discrete", hist=FALSE))
  expect_no_condition(plot_single_histo(histo_data2, case="continuous",hist=FALSE))
  expect_no_condition(plot_single_histo(histo_data2, hist=FALSE, abs=TRUE))
  expect_no_condition(plot_single_histo(histo_data2, case="discrete", hist=FALSE, abs=TRUE))
  expect_no_condition(plot_single_histo(histo_data2, case="continuous",hist=FALSE, abs=TRUE))
  expect_no_condition(plot_single_histo(histo_data2, abs=TRUE))
  expect_no_condition(plot_single_histo(histo_data2, case="discrete", abs=TRUE))
  expect_no_condition(plot_single_histo(histo_data2, case="continuous", abs=TRUE))
  
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4)))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), case="continuous"))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), case="discrete"))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), hist=FALSE))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), case="continuous", hist=FALSE))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), case="discrete", hist=FALSE))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), abs=TRUE))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), abs=TRUE, case="continuous"))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), case="discrete", abs=TRUE))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), hist=FALSE, abs=TRUE))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), case="continuous", abs=TRUE, hist=FALSE))
  expect_no_condition(plot_single_histo(histo_disc%>%filter(group_id==4), case="discrete", hist=FALSE, abs=TRUE))
  
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4)))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), case="continuous"))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), case="discrete"))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), hist=FALSE))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), case="continuous", hist=FALSE))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), case="discrete", hist=FALSE))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), abs=TRUE))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), abs=TRUE, case="continuous"))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), case="discrete", abs=TRUE))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), hist=FALSE, abs=TRUE))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), case="continuous", abs=TRUE, hist=FALSE))
  expect_no_condition(plot_single_histo(histo_conti%>%filter(group_id==4), case="discrete", hist=FALSE, abs=TRUE))
  
  expect_error(plot_single_histo(histo_data2, case="test"))
  expect_error(plot_single_histo(histo_data2, abs="test"))
  expect_error(plot_single_histo(histo_data2, hist="test"))
  expect_error(plot_single_histo(histo_data2, xlab=1))
  expect_error(plot_single_histo(histo_data2, ylab=2))
  expect_error(plot_single_histo(histo_data2, ylim="test"))
  expect_warning(plot_single_histo(histo_data2, main="test", automatic_main = TRUE))
  expect_error(plot_single_histo(histo_data))
  expect_error(plot_single_histo(histo_data2, main=histo_data))
  
  histo_zero<-histo
  histo_zero$weighted_counts<-0
  expect_no_condition(interactive_plots(histo_zero))
  expect_no_condition(interactive_plots(histo_zero, case="discrete"))
  expect_no_condition(interactive_plots(histo_zero, case="continuous"))
  expect_no_condition(interactive_plots(histo_zero, hist=FALSE))
  expect_no_condition(interactive_plots(histo_zero, case="discrete", hist=FALSE))
  expect_no_condition(interactive_plots(histo_zero, case="continuous",hist=FALSE))
  expect_no_condition(interactive_plots(histo_zero, hist=FALSE, abs=TRUE))
  expect_no_condition(interactive_plots(histo_zero, case="discrete", hist=FALSE, abs=TRUE))
  expect_no_condition(interactive_plots(histo_zero, case="continuous",hist=FALSE, abs=TRUE))
  expect_no_condition(interactive_plots(histo_zero, abs=TRUE))
  expect_no_condition(interactive_plots(histo_zero, case="discrete", abs=TRUE))
  expect_no_condition(interactive_plots(histo_zero, case="continuous", abs=TRUE))
  

  })



