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
fs.data <- cbind(aggregate(SARS2_ELISA ~codeloc,sum,data = cat.data),
      samples = aggregate(katserum ~codeloc,sum,data = cat.data)[1:8,2])
fs.data<- fs.data[fs.data$SARS2_ELISA>0&fs.data$samples>1,]#only those where at least one infectious and more than 1 sample is present
fs.data

x.cats <- fs.data$SARS2_ELISA-1 #cases 
s0.cats <- fs.data$samples -1   #susceptible cats at start is samples minus 1
i0.cats<- rep(1,length(fs.data$codeloc)) #initial infectious is 1
FinalSize(x.cats,s0.cats,i0.cats,max.val = 200)
FinalSize(x.cats[1],s0.cats[1],i0.cats[1],max.val = 1)
FinalSize(x.cats[2],s0.cats[2],i0.cats[2],max.val = 10)
