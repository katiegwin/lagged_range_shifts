################################################################################
################################################################################
###                                                                          ###
###                            ANALYSIS SCRIPT 1:                            ###
###                 RECRUITMENT WITH RANGE POSITION MODELLING                ###
###                   HAS OPTIMAL RECRUITMENT PEAK SHIFTED?                  ###
###                                                                          ###
################################################################################
################################################################################


# Created Nov 10, 2022
# Last Updated Dec 2, 2024
# Authors: Katie Goodwin

# Analyses for Goodwin et al. 2024 Lagged climate-driven range shifts at species’ leading, but not trailing, range edges revealed by multispecies seed addition experiment. Ecography. DOI: 10.1111/ecog.07331
# testing whether optimal seedling recruitment has shifted from the adult range centre to cooler range regions. (see Fig. 3)
# tests the lagged response hypothesis at the optima

###########################################################################################
#### OUTPUTS FROM THIS SCRIPT ------------------------------------------------
###########################################################################################
load("outputs/models/optimarecruit_biascorrpredict.Rdata") # predicted recruitment with temp range position with bias-corrected CIs
load("outputs/models/optimalzi_biascorrpredict.Rdata") # predicted recruitment with temp range position with bias-corrected CIs (zero inflaiton part) 
load("outputs/models/ppt_biascorrpredict.Rdata") # predicted recruitment with ppt range position with bias-corrected CIs

# these data frames were then inputed into scripts/outputscripts/pretty_plots.R to generate manuscript figs.

###########################################################################################
#### INPUTS AND PACKAGES ------------------------------------------------
###########################################################################################
library(glmmTMB)
library(tidyverse)
library(ggplot2)
library(ggeffects)
library(DHARMa)
library(GGally)

# count census data 
load('data/recruit_census.RData') 
str(census)

###########################################################################################
##### MISC. INFO ---------------------------------------
###########################################################################################

ggpairs(census, columns = c("temprange_pos", "pptrange_pos", "canopycont")) # not correlated

##### SOME RANDOM NUMBERS FOR MANUSCRIPT
# how many species-rep combos were zero? (25 species in 180 plots)
censuszero <- subset(census, new_all == "0") # 4107 zeroes
4107 / 4500
# 91.3 % of species-rep combos were zero!

# how many plots did each species recruit in
plotsgerm <- subset(census, new_all > 0) # subset data set of non-zero recruitment
plotscount <- plotsgerm %>% # count of number of plots species had non-zero recruitment
  group_by(species) %>%
  summarize(count = n())



###########################################################################################
##### RECRUITMENT COUNT RECRUITMENT OPTIMA MODEL---------------------------
###########################################################################################

# new data frame for species that germinated in at least 8 plots to be included in this analysis
census8 <- census %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ABILAS" | species == "ANEOCC" |
           species == "ERIPER" | species == "PICENG" |
           species == "RUBURS" | species == "SORSIT" | species == "TELGRA" |
           species == "TOLMEN" | species == "VACDEL" | species == "VACPAR" |
           species == "MAHNER" | species == "ERILAN" | species == "MAHAQU") # remove species with indiv in less than 8 plots

# remove values where new_all = NA (i.e., was more recruitment in control than seed addition plot)
census8a <- subset(census8, new_all != "NA")
census8 <- census8a


### mopt -optimal recruitment model
mopt <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2) 
              + (1|site/rep) + (1|family), 
              ziformula = ~ poly(temprange_pos, 2) + poly(pptrange_pos, 2) 
              + (1|site/rep) + (1|family),
              data = census8, 
              family = nbinom2(link = "log")) 
summary(mopt)


#### MODEL DIAGNOSTICS ---
# residual plots using DHARMa
m1_simres <- simulateResiduals(mopt)
testUniformity(m1_simres) 
testOutliers(m1_simres) 
testQuantiles(m1_simres)
plotResiduals(m1_simres, census8$temprange_pos)
sd(resid(mopt))
plot(census8$temprange_pos, m1_simres$scaledResiduals, 
     ylab="Simulated Residuals", xlab="Thermal Range Position") 


##### BOOTSTRAPPING  -----
# i.e., do 95% confidence intervals of peak of curve include zero?
# separate code file to generate bootstraps as this loop takes >24 hours and was run on the zoology cluster
# (see tempbootstrap_cond_cluster.R, tempbootstrap_zi_cluster.R, pptbootstrap_cluster.R in scripts/analysis/bootstraps/ folder for where these all came from)

### import bootstrap runs from cluster server
# topt, toptzi, and popt record what the optima was for each bootstrap to generate confidence intervals around the recruitment optimum
# tpred, tpredzi, and ppred are the predicted recruitment counts for each bootstraps to generate 95% confidence intervals for Figures.

# temperature - conditional part of model
load('outputs/bootstraps/temp-cond/temprange_cond_bootstraps08DEC22.RData')
load('outputs/bootstraps/temp-cond/temprange_cond_optimaCI08DEC22.RData')

topt <- optima
tpred <- pred

# temperature - zero inflation part of model 
load('outputs/bootstraps/temp-zi/temprange_zi_bootstraps01DEC22.RData')

tpredzi <- pred

# ppt - conditional part of model
load('outputs/bootstraps/ppt/pptrange_cond_bootstraps01DEC22.RData')
load('outputs/bootstraps/ppt/pptrange_cond_optimaCI01DEC22.RData')

popt <- optima
ppred <- pred

# fix classes for optima dataframes
topt$opt_temprange <- as.numeric(as.character(topt$opt_temprange))
topt$opt_recruits <- as.numeric(as.character(topt$opt_recruits))
str(topt)

popt$opt_pptrange <- as.numeric(as.character(popt$opt_pptrange))
popt$opt_recruits <- as.numeric(as.character(popt$opt_recruits))
str(popt)


##### CI FOR PEAK OF RECRUITMENT CURVE -----

# temp - cond - explore distribution of recruitment optimum from bootstraps
hist(topt$opt_temprange, breaks = 50)
quantile(topt$opt_temprange, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
quantile(topt$opt_temprange, probs = c(0.9394), na.rm = TRUE) # 93.94 % of bootstraps have optimal thermal recruitment below zero in cold part of range

# ppt
hist(popt$opt_pptrange, breaks = 30)
quantile(popt$opt_pptrange, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
quantile(popt$opt_pptrange, probs = c(0.0554), na.rm = TRUE, names = TRUE, digits = 10) # 5.54 % of bootstraps have optimal thermal recruitment below zero
100 - 5.54 # proportion of bootstrap optima in wet part of range
summary(popt$opt_pptrange) # mean = 0.128, median = 0.120


##### model predictions ----
#  temp conditional part
predmopt <- expand.grid(temprange_pos = seq(-2.6,1.6, by = .01), pptrange_pos = 0, 
                    site = "NA", rep = "NA", family = "NA") # note: no interaction terms in models so setting these RE to means and holding ppt range position fix doesn't affect location of recruitment optimum
predmopt$pred <- predict(mopt, newdata = predmopt, type = "conditional", se = F,  allow.new.levels = T) 
predmopt$temprange_pos[which.max(predmopt$pred)] # peak of curve is -0.26

# ppt
predmoptp <- expand.grid(pptrange_pos = seq(-1.6, 2.3, by = .01), temprange_pos = 0, 
                      site = "NA", rep = "NA", family = "NA") 
predmoptp$pred <- predict(mopt, newdata = predmoptp, type = "conditional", se = F,  allow.new.levels = T) 
predmoptp$pptrange_pos[which.max(predmoptp$pred)] # peak of curve is 0.17

# temp zero inflation part
predmoptz <- expand.grid(temprange_pos = seq(-2.6, 1.6, by = .01), pptrange_pos = 0,
                       site = "NA", rep = "NA", family = "NA")
predmoptz$pred <- predict(mopt, newdata = predmoptz, type = "zprob", se = F,  allow.new.levels = T) 


##### GENERATE BIAS-CORRECTED 95% CONFIDENCE INTERVALS ----
# code adapted from: doi: 10.3389/fpsyg.2022.810258, Sheth and Angert 2018

### temperature conditional part (i.e., creating dataframe for what will be figure 3)
tboot<- tpred[sapply(tpred, is.numeric)] # remove character columns that are all NAs (i.e., bootstrap run didnt converge) 

# add in actual predicted values
tpredictopt <- predmopt[, c(1,6)]
tpred1 <- cbind(tpredictopt,tboot, by = "temprange_pos")
tpred <- tpred1[, -c(3)]

# proportion of bootstraps less than predicted and converted to z score 
head(tpred) # check first and last column that converged for next line of code
tpredz <- tpred %>% 
  rowwise() %>% 
  mutate(z = qnorm(sum(c_across(V10:V10003) < pred) / length(c_across(V10:V10003)))) 

# calculate new probabilities at which to place bias corrected lower and upper CIs
tpredictprob <-tpredz %>%
  mutate(lower.prob = pnorm(2*z-1.96),# bias corrected lower quantile to calculate lower CI
         upper.prob = pnorm(2*z+1.96)) # bias corrected upper quantile to calculate upper CI

# bias corrected CIs
tpredictbias <- tpredictprob %>% 
  rowwise() %>% 
  transmute(temprange_pos = temprange_pos,
            pred = pred,
            lower.ci = quantile((c_across(V10:V10003)), probs = lower.prob),
            upper.ci = quantile((c_across(V10:V10003)), probs = upper.prob))

# save prediction with bias corrected 95% CI
save(tpredictbias, file = "outputs/models/optimarecruit_biascorrpredict.Rdata")

### temperature zero inflation part (i.e., creating dataframe for what will be Fig. S8)
tbootzi <- tpredzi[sapply(tpredzi, is.numeric)] # remove character columns that are all NAs (i.e., boostrap runs that didnt converge)

# add in actual predicted values
tpredictzi <- predmoptz[, c(1,6)]
tpredz1 <- cbind(tpredictzi,tbootzi, by = "temprange_pos")
tpredz <- tpredz1[, -c(3)]

head(tbootzi) # check first and last column that converged for next line of code
# proportion of bootstraps less than predicted and converted to z score
tpredziz <- tpredz %>% 
  rowwise() %>% 
  mutate(z = qnorm(sum(c_across(V10:V20004) < pred) / length(c_across(V10:V20004)))) 
str(tpredziz$z)

# calculate new probabilities at which to place bias corrected lower and upper CIs
tpredictziprob <-tpredziz %>%
  mutate(lower.prob = pnorm(2*z-1.96),# bias corrected lower quantile to calculate
         upper.prob = pnorm(2*z+1.96)) # bias corrected upper quantile to calculate
str(tpredictziprob$lower.prob)
str(tpredictziprob$upper.prob)

# bias corrected CIs
tpredictzibias <- tpredictziprob %>% 
  rowwise() %>% 
  transmute(temprange_pos = temprange_pos,
            pred = pred,
            lower.ci = quantile((c_across(V10:V20004)), probs = lower.prob),
            upper.ci = quantile((c_across(V10:V20004)), probs = upper.prob))

# save prediction with bias corrected 95% CI
save(tpredictzibias, file = "outputs/models/optimalzi_biascorrpredict.Rdata")

### precipitation (i.e., creating dataframe that will become Fig. S6)
pboot<- ppred[sapply(ppred, is.numeric)] # remove character columns that are all NAs (i.e., bootstrap didnt converge) 

# add in actual predicted values
ppredictopt <- predmoptp[, c(1,6)]
ppred1 <- cbind(ppredictopt,pboot, by = "temprange_pos")
ppred <- ppred1[, -c(3)]

head(ppred) # check range of columns that converged for next line of code
# proportion of bootstraps less than predicted and converted to z score
ppredz <- ppred %>% 
  rowwise() %>% 
  mutate(z = qnorm(sum(c_across(V9:V10002) < pred) / length(c_across(V9:V10002))))
str(ppredz$z)

# calculate new probabilities at which to place bias corrected lower and upper CIs
ppredictprob <-ppredz %>%
  mutate(lower.prob = pnorm(2*z-1.96),# bias corrected lower quantile to calculate
         upper.prob = pnorm(2*z+1.96)) # bias corrected upper quantile to calculate
str(ppredictprob$lower.prob)
str(ppredictprob$upper.prob)

# bias corrected CIs
ppredictbias <- ppredictprob %>% 
  rowwise() %>% 
  transmute(pptrange_pos = pptrange_pos,
            pred = pred,
            lower.ci = quantile((c_across(V9:V20005)), probs = lower.prob),
            upper.ci = quantile((c_across(V9:V20005)), probs = upper.prob))

### save prediction with bias corrected 95% CI
save(ppredictbias, file = "outputs/models/ppt_biascorrpredict.Rdata")

##### SENSITIVITY TEST - ADDRESSING UNCERTAINTY IN OCCURRENCE RECORDS ----
# some of these species dont have many adult occurrence records and thus have higher uncertainty in adult climatic ranges. does removing those change results? 
# remove species with < 50 records ANEOCC, VACDEL, TELGRA, TOLMEN
# CONCLUSION: peak of recruitment curve is still in cooler part of thermal range, doesnt change results but full model doesnt run with this small a sample size so not using this for main analyses just to build confidence in results
census.highgbif <- census8 %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ABILAS" | 
           species == "ERIPER" | species == "PICENG" |
           species == "RUBURS" | species == "SORSIT" | 
           species == "VACPAR" |
           species == "MAHNER" | species == "ERILAN" | species == "MAHAQU") 

# full model doesnt converge, removed rep and ppt range position for it to run
mopttest <- glmmTMB(new_all ~  poly(temprange_pos, 2) 
                    + (1|site) + (1|family), 
                    ziformula = ~ poly(temprange_pos, 2) 
                    + (1|site) + (1|family),
                    data = census.highgbif, 
                    family = nbinom2(link = "log")) 
summary(mopttest)


mopttestpred <- as.data.frame(ggpredict(mopttest, terms = c("temprange_pos [all]"), type = "fe"))

#figure S4A
opttest <- ggplot(mopttestpred) +
  geom_point(data = census.highgbif, aes(x = temprange_pos, y = log(new_all + 1) ), col ="black", alpha = 0.4) +
  geom_line(aes(x = x, y = log(predicted + 1)), col = "black") +
  geom_ribbon(aes(x = x, ymin = log(conf.low + 1), ymax = log(conf.high +1)), fill = "black", alpha=0.3) +
  scale_x_reverse() +
  ylab("ln (Seedling Count + 1)") + xlab("Thermal Range Position") + 
  theme_classic()
opttest
# saved plot of this in recruitment_edge_analysis.R 



