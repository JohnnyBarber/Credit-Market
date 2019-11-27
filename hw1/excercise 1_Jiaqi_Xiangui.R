library(data.table)
library(readxl)

setwd("C:/Users/Jiaqi Li/Desktop/class materials/quarter 4/Credit Markets/hw1")

data = read_excel("sample_logit.xlsx")
data = as.data.table(data)
summary(data)
default = data[default_f1 == 1,]
non_default = data[default_f1 == 0,]

names(data)
for(i in names(data)){
  assign(i,lm(default_f1 ~ get(i), data = data))
}

correlation_data = data[,-c(1,12)]

cor(correlation_data)[,"default_f1"]

logit = glm(default_f1~., data = correlation_data, family = "binomial")
summary(logit)
