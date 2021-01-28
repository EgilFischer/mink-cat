#############################################################
#        Transmission of SARS-CoV2 on mink farms            #
#           from mink to cat                                #
#           between cats                                    #
#                                                           #
#       Data Analists: Egil A.J. Fischer / Anna van Aart    #
#                                                           #
#############################################################

## load packages ####
library(readxl)

## load data  ####
cat.data = read_xlsx("C:/Surfdrive/Projecten/Covid/CatsAndMink/mink-cat/Dataset hond-kat.xlsx")
cat.data$exposure = 35
## remove dog data ###
cat.data = cat.data[cat.data$katspec != "dog",]

## Mink to cat ####
#Assumption is that the exposure to a infected farm causes cats to be infected with a constant infection rate
#The probability p of a cat being infected given a certain exposure time t =
#p = 1 - (1- E^(beta*t))
# to estimate this value we use a GLM with a complementary log-log link function
fit <- glm(SARS2_ELISA ~ 1,offset = log(1 * exposure),family = binomial(link = "cloglog"), data = cat.data)
summary(fit)

exp(summary(fit)$coefficients[1])
exp(sum(summary(fit)$coefficients[1:2]))
exp(-diff(summary(fit)$coefficients[1:2]))

#use farm location as covariate
fit <- glm(SARS2_ELISA ~ codeloc,offset = log(1 * exposure),family = binomial(link = "cloglog"), data = cat.data)
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

## Calculate R0
## calculation assumes that R0 is above 1 and that the final size has been reached in a deterministic model
R0 = function(z){-log(1-z)/z}

#estimate overall prevalence
R0(mean(cat.data$SARS2_ELISA, na.rm = T))

# estimate prevalence and R0 per farm
prevalencesR0 = aggregate(SARS2_ELISA ~ codeloc,FUN = mean, data = cat.data)
prevalencesR0$R0= sapply(X = prevalencesR0$SARS2_ELISA,FUN = R0)
prevalencesR0

# estimate R with final size distribution for each farm 



