#############################################################
#        SARS-CoV2 and companion animals on mink farms      #
#           - Risk factors                                  #
#           - Transmission                                  #
#            from mink to cat                               #
#            between cats (R0 final size)                   #
#                                                           #
#       Authors:                                            #
#       Egil A.J. Fischer (e.a.j.fischer@uu.nl)             #
#       Anna E. van Aart (a.e.vanaart@students.uu.nl)       #
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
adult.cat.data = cat.data[cat.data$katspec != "kitten",]

### Possible risk factors SARS-CoV-2 infection blood sampled adult feral cats ##
attach(adult.cat.data)
names(adult.cat.data)
summary(adult.cat.data)

#number of records
length(adult.cat.data$katnum)

#positive cats
cat.data$sarscov2pos <- katpcr_o + SARS2_ELISA
table(sarscov2pos)

#number of cats with a valid serological or pcr test
length(sarscov2pos[!is.na(sarscov2pos)])


# positive (%) , negative (%)
table(SARS2_ELISA, katgesla)
prop.table(table(SARS2_ELISA, katgesla))
table(SARS2_ELISA, katleeft)
table(SARS2_ELISA, katpet)
table(SARS2_ELISA, katlact)
table(SARS2_ELISA, katprgnt)


fisher.test(SARS2_ELISA, katgesla) #sex
t.test(katleeft ~ SARS2_ELISA, var.equal=TRUE) #age
fisher.test(SARS2_ELISA, katpet) #pet or stray
fisher.test(SARS2_ELISA, katlact) #lactating
fisher.test(SARS2_ELISA, katprgnt) #pregnant

# risk factors pregnancy, lactation and sex -> OR + 95% CL was calculated with 2x2-table (manually)
detach(adult.cat.data)


#############################################################################
#            Estimation of transmission                                     #
#############################################################################
## calculate exposure time
exposure.data.average <-  aggregate(as.Date(exposure.data$culling, format = "%Y%m%D"),
                                    by =list(exposure.data$`code loc`),
                                    FUN = mean,id.vals =) ;
colnames(exposure.data.average) <-  c("code.loc","culling.date");
exposure.data.average$diagnose <-  as.Date(aggregate(exposure.data$diagnose,
                                             by =list(exposure.data$`code loc`),
                                             FUN = mean)[,2], format = "%Y%m%D")  
exposure.data.average$firstsymp1 <-  as.Date(aggregate(exposure.data$firstsymp1,
                                               by =list(exposure.data$`code loc`),
                                               FUN = mean)[,2], format = "%Y%m%D") 
exposure.data.average$firstsymp2 <-  as.Date(aggregate(exposure.data$firstsymp2,
                                               by =list(exposure.data$`code loc`),
                                               FUN = mean)[,2], format = "%Y%m%D") 
exposure.data.average$firstsympmean <- as.Date(aggregate(exposure.data$firstsympmean,
                                                 by =list(exposure.data$`code loc`),
                                                 FUN = mean)[,2], format = "%Y%m%D") 

#set the exposure time in days
cat.data$katdatum <- as.Date(cat.data$katdatum , format = "%Y%m%D") 
cat.data$exposure <- sapply(c(1:length(cat.data$katnum)), 
                            FUN = function(x){difftime(min(cat.data[x,]$katdatum,exposure.data.average[exposure.data.average$code.loc == cat.data[x,]$codeloc,]$culling.date),
                                                       exposure.data.average[exposure.data.average$code.loc == cat.data[x,]$codeloc,]$firstsympmean,unit = "days")})  
table(data.frame(cat.data$codeloc,cat.data$exposure))

adult.cat.data <- cat.data[cat.data$katspec != "kitten",]
## Mink to cat ####
#Assumption is that the exposure to a infected farm causes cats to be infected with a constant infection rate
#The probability p of a cat being infected given a certain exposure time t =
#p = 1 - (1- E^(-beta*t))
# to estimate this value we use a GLM with a complementary log-log link function
fit <- glm(SARS2_ELISA ~ 1,offset = log(1 * exposure),family = binomial(link = "cloglog"), data = cat.data)
summary(fit)

point.est <- exp(summary(fit)$coefficients[1])
point.est
exp(-diff(summary(fit)$coefficients[1:2]))
exp(sum(summary(fit)$coefficients[1:2]))


# results in
1-exp(-point.est * mean(cat.data$exposure))
1-exp(-point.est * mean(cat.data[cat.data$codeloc ==1, ]$exposure))
1-exp(-point.est * mean(cat.data[cat.data$codeloc ==4, ]$exposure))
sort(1-exp(-point.est * unique(cat.data$exposure)))
summary(1-exp(-point.est * unique(cat.data$exposure)))

#Expected number of infected cats per farm
exp.num.inf<- cbind(aggregate(1- exp(-point.est*cat.data$exposure),
          mean,
          by = list(NB = cat.data$codeloc)),
          aggregate(1- exp(-point.est*cat.data$exposure),
                    length,
                    by = list(NB = cat.data$codeloc)))
colnames(exp.num.inf)<-c("NB","E(prev)","NB","tested")
exp.num.inf$Obs <- aggregate(cat.data$SARS2_ELISA,
                             sum,
                             by = list(NB = cat.data$codeloc),
                             na.rm = T)$x          


exp.num.inf$Expected <- exp.num.inf$tested * exp.num.inf$`E(prev)`
exp.num.inf$pObs <- (mapply(pbinom,exp.num.inf$Obs,exp.num.inf$tested,exp.num.inf$`E(prev)`, list(lower.tail = T)))

#chisq.test(x = as.table(rbind(exp.num.inf$tested,exp.num.inf$Obs)))
prop.test(n = exp.num.inf$tested,
          x = exp.num.inf$Obs, 
          p = exp.num.inf$`E(prev)`)


## R0
# estimate R with final size distribution for each farm 
fs.data <- cbind(aggregate(SARS2_ELISA ~codeloc,sum,data = adult.cat.data),
      samples = aggregate(katserum ~ codeloc,length,data = adult.cat.data)[1:8,2])
fs.data<- fs.data[fs.data$SARS2_ELISA>0&fs.data$samples>1,]#only those where at least one infectious and more than 1 sample is present
fs.data

x.cats <- fs.data$SARS2_ELISA-1 #cases 
s0.cats <- fs.data$samples -1   #susceptible cats at start is samples minus 1
i0.cats<- rep(1,length(fs.data$codeloc)) #initial infectious is 1
FinalSize(x.cats,s0.cats,i0.cats)

#exclude farm 6
FinalSize(x.cats[1:2],s0.cats[1:2],i0.cats[1:2])

