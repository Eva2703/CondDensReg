library(dplyr)
data<-read.csv("raw_age_birth.csv",sep=";", header=FALSE)
data<-data[-c(1:10),]
data<-data[-c(1217:1220),]
colnames(data)<-c("year", "age", "married_male", "unmarried_male", "married_female", "unmarried_female")
unique(data$age)
data$age<-rep(c("<15",15:49,">=50", "unknown"),32)
data$year<-as.numeric(data$year)

data<-data%>%filter(age!="unknown"&age!="<15"&age!=">=50")
data$age<-as.numeric(data$age)
data_long<-data[rep(1:1120,each=4),c(1,2)]
data_long$marital_status<-factor(rep(c("married","unmarried"), 2240))
data_long$sex<-factor(rep(rep(c("male","female"),each=2),1120))

data_long$counts<-as.vector(rbind(data$married_male,data$unmarried_male,data$married_female,data$unmarried_female))

data_long$counts<-as.numeric(data_long$counts)
data_long$counts[is.na(data_long$counts)] <- 0
data_long$age<-data_long$age+0.5

data_long<-na.omit(data_long)
data_age_birth<-data_long
save(data_age_birth, file="data_age_birth.RData")
