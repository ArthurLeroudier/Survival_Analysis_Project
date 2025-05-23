###IMPORTS

library(KMsurv)
library(survival)
library(ggplot2)
library(survminer)
library(data.table)
library(joineRML)
library(tidyr)


#setup
set.seed(0)
data(pbc2)


###PART A

##Missing Values
colSums(is.na(pbc2))
missing_row <- data.frame(id= pbc2$id, missing = rowSums(is.na(pbc2)))
missing_row <-missing_row[missing_row$missing>0,]


##Constant variables
dif_values <- data.table(pbc2)[, lapply(.SD, function(x) length(unique(x))), by=id]
#^for each individual we check how many different values for each variable
check_constants <- colSums(dif_values[,!'id'])
check_constants #if the value is 312, the variable is constant for a given individual
#constant variables are: years, status, drug, age, sex, status2
pbc2_constants <- unique(pbc2[,c('years', 'status', 'drug', 'age', 'sex', 'status2')])
hist(pbc2$years)
barplot(prop.table(table(pbc2_constants$status)))
barplot(prop.table(table(pbc2_constants$drug)))
hist(pbc2$age)
barplot(prop.table(table(pbc2_constants$sex)))
barplot(prop.table(table(pbc2_constants$status2)))

## Quantitative variables
#Mean curve
plot(pbc2$year, pbc2$serBilir)
abline(lm(pbc2$serBilir ~ pbc2$year), col = "red")
plot(pbc2$year, pbc2$serChol)
abline(lm(pbc2$serChol ~ pbc2$year), col = "red")
plot(pbc2$year, pbc2$alkaline)
abline(lm(pbc2$alkaline ~ pbc2$year), col = "red")
plot(pbc2$year, pbc2$SGOT)
abline(lm(pbc2$SGOT ~ pbc2$year), col = "red")
plot(pbc2$year, pbc2$platelets)
abline(lm(pbc2$platelets ~ pbc2$year), col = "red")
plot(pbc2$year, pbc2$prothrombin)
abline(lm(pbc2$prothrombin ~ pbc2$year), col = "red")

ggplot(pbc2[c('id','serBilir','year')], aes(year, serBilir, group = id, color = id)) + 
  geom_point() + 
  geom_line() +
  theme(legend.position = "none")
ggplot(pbc2[c('id','serChol','year')], aes(year, serChol, group = id, color = id)) + 
  geom_point() + 
  geom_line() +
  theme(legend.position = "none")
ggplot(pbc2[c('id','alkaline','year')], aes(year, alkaline, group = id, color = id)) + 
  geom_point() + 
  geom_line() +
  theme(legend.position = "none")
ggplot(pbc2[c('id','SGOT','year')], aes(year, SGOT, group = id, color = id)) + 
  geom_point() + 
  geom_line() +
  theme(legend.position = "none")
ggplot(pbc2[c('id','platelets','year')], aes(year, platelets, group = id, color = id)) + 
  geom_point() + 
  geom_line() +
  theme(legend.position = "none")
ggplot(pbc2[c('id','prothrombin','year')], aes(year, prothrombin, group = id, color = id)) + 
  geom_point() + 
  geom_line() +
  theme(legend.position = "none")
