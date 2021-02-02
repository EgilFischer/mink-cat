#############################################################
#        Transmission of SARS-CoV2 on mink farms            #
#           from mink to cat                                #
#           between cats                                    #
#                                                           #
#       Data Analists: Egil A.J. Fischer / Anna van Aart    #
#                                                           #
#############################################################

## load packages / user defined functions ####
library(readxl)
source("FinalSizeFastImplementation.R")

## load data  ####
full.data = read_xlsx("C:/Surfdrive/Projecten/Covid/CatsAndMink/mink-cat/Dataset hond-katFV20210128EF20210202.xlsx",
                      sheet = "NC_KAT 11-Aug-2020",
                      col_types = c("numeric","text","numeric","date",rep("guess",15)))
exposure.data = read_xlsx("C:/Surfdrive/Projecten/Covid/CatsAndMink/mink-cat/Dataset hond-katFV20210128EF20210202.xlsx",
                          sheet = "exposure time")



## remove dog data ###
cat.data = full.data[full.data$katspec != "dog",]

## calculate exposure time
exposure.data.average <-  aggregate(exposure.data$culling,
                                  by =list(exposure.data$`code loc`),
                                  FUN = mean,id.vals =);
colnames(exposure.data.average) <-  c("code.loc","culling.date");
exposure.data.average$diagnose <-  aggregate(exposure.data$diagnose,
                                  by =list(exposure.data$`code loc`),
                                  FUN = mean)[,2]
exposure.data.average$firstsymp1 <-  aggregate(exposure.data$firstsymp1,
                                           by =list(exposure.data$`code loc`),
                                           FUN = mean)[,2]
exposure.data.average$firstsymp2 <-  aggregate(exposure.data$firstsymp2,
                                           by =list(exposure.data$`code loc`),
                                           FUN = mean)[,2]
exposure.data.average$firstsympmean <- aggregate(exposure.data$firstsympmean,
                                                 by =list(exposure.data$`code loc`),
                                                 FUN = mean)[,2]



#set the exposure time in days
cat.data$exposure <- sapply(c(1:length(cat.data$katnum)), 
                            FUN = function(x){difftime(cat.data[x,]$katdatum,
                                                        exposure.data.average[exposure.data.average$code.loc == cat.data[x,]$codeloc,]$firstsympmean,unit = "days")})  






## Mink to cat ####
#Assumption is that the exposure to a infected farm causes cats to be infected with a constant infection rate
#The probability p of a cat being infected given a certain exposure time t =
#p = 1 - (1- E^(beta*t))
# to estimate this value we use a GLM with a complementary log-log link function
fit <- glm(SARS2_ELISA ~ 1,offset = log(1 * exposure),family = binomial(link = "cloglog"), data = cat.data)
summary(fit)

point.est <- exp(summary(fit)$coefficients[1])
point.est
exp(sum(summary(fit)$coefficients[1:2]))
exp(-diff(summary(fit)$coefficients[1:2]))

# results in
1-exp(-point.est * 31)
1-exp(-point.est * 41)


#use farm location as covariate
fit <- glm(SARS2_ELISA ~ codeloc,offset = log(1 * exposure),family = binomial(link = "cloglog"), data = cat.data)
summary(fit)
drop1(fit)

exp(summary(fit)$coefficients[1])
exp(sum(summary(fit)$coefficients[1:2]))
exp(-diff(summary(fit)$coefficients[1:2]))

#use farm location as covariate
fit <- glm(SARS2_ELISA ~ as.numeric(1-cat.data$codeloc %in%c(1,4)),offset = log(1 * exposure),family = binomial(link = "cloglog"), data = cat.data)
summary(fit)
drop1(fit)

exp(summary(fit)$coefficients[1])
exp(sum(summary(fit)$coefficients[1:2]))
exp(-diff(summary(fit)$coefficients[1:2]))

#select farms NB1 and NB4
fit <- glm(SARS2_ELISA ~ codeloc,offset = log(1 * exposure),family = binomial(link = "cloglog"), data = cat.data[cat.data$codeloc == 1||cat.data$codeloc == 4],)
summary(fit)
drop1(fit)

exp(summary(fit)$coefficients[1])
exp(sum(summary(fit)$coefficients[1:2]))
exp(-diff(summary(fit)$coefficients[1:2]))

## R0
# estimate R with final size distribution for each farm 
fs.data <- cbind(aggregate(SARS2_ELISA ~codeloc,sum,data = cat.data),
      samples = aggregate(katserum ~ codeloc,length,data = cat.data)[1:8,2])
fs.data<- fs.data[fs.data$SARS2_ELISA>0&fs.data$samples>1,]#only those where at least one infectious and more than 1 sample is present
fs.data

x.cats <- fs.data$SARS2_ELISA-1 #cases 
s0.cats <- fs.data$samples -1   #susceptible cats at start is samples minus 1
i0.cats<- rep(1,length(fs.data$codeloc)) #initial infectious is 1
FinalSize(x.cats,s0.cats,i0.cats)
