################################################################################
################################################################################
###                                                                          ###
###                            ANALYSIS SCRIPT 3:                            ###
###                 RECRUITMENT WITH RANGE POSITION MODELLING                ###
###               how do survival and growth vary across range?              ###
###                                                                          ###
################################################################################
################################################################################


# Created Nov 28, 2022
# Last Updated Dec. 4, 2024
# Authors: Katie Goodwin

# Analyses for Goodwin et al. 2024 Lagged climate-driven range shifts at species’ leading, but not trailing, range edges revealed by multispecies seed addition experiment. Ecography. DOI: 10.1111/ecog.07331
# testing whether range patterns in survival and growth alter initial recruitment results

###########################################################################################
#### PACKAGES AND DATA ------------------------------------------------
###########################################################################################
library(glmmTMB)
library(tidyverse)
library(ggplot2)
library(GGally)
library(ggeffects)
library(DHARMa)
library(cowplot)

## height data
load('data/year1size.RData') # size data frame combining both years with mean height and leaf count for each plot

# Survival data
load('data/recruit_surv1yr.RData') # data set of surviving one year after recruitment
load('data/recruit_surv2yrs.RData') # data set of surviving two years after recruitment 
load('data/recruit_survto2022.RData') # data set of surviving to 2022 (up to 5 years)


################################################################################
##### HEIGHT 1st year with range position -----
################################################################################

# remove one very large outlier (testOutliers() found outliers with full dataset, removing this obs fixes the issue)
# summary of outputs the same either way removing outlier makes no difference to inferences.
# for now removing this outlier since the model diagnostics look better. 
year1sum <- subset(year1sum, stht < 4)


#### model - including region does not change results
mheight <- glmmTMB(stht ~ temprange_pos  + pptrange_pos + 
                     (1|site/rep),  
                   data = year1sum, family = gaussian) 
summary(mheight)

# quadratic
mheightq <- glmmTMB(stht ~ poly(temprange_pos, 2)  + poly(pptrange_pos, 2) + 
                     (1|site/rep),  
                   data = year1sum, family = gaussian) 
summary(mheightq)

AIC(mheight, mheightq) # linear relationship is better fit

#### model diagnostics
# residual plots using DHARMa
mheight_simres <- simulateResiduals(fittedModel = mheight, plot = T)

testUniformity(mheight_simres) 
testOutliers(mheight_simres)
testQuantiles(mheight_simres)

##### BOOTSTRAPPING  ---------------
# i.e., to create confidence intervals for prediction line
# separate code file as this loop takes >24 hours and is run on the zoology cluster
# (see heightbootstrap_cluster.R)

###### import bootstrap runs from cluster server 

# temperature - conditional part of model
load('outputs/bootstraps/growth/height_pred08DEC22.RData')

hpred <- pred
summary(year1sum$temprange_pos)

##### model predictions (not bootstrap but for solid part of prediction line) ----
predmheight <- expand.grid(temprange_pos = seq(-2.4,1.3, by = .01), pptrange_pos = 0, 
                           site = "NA", rep = "NA") # for predicted values. will take pop means for random effects. Since no interaction terms value of canopy or ppt shouldnt matter
predmheight$pred <- predict(mheight, newdata = predmheight, type = "response", se = F,  allow.new.levels = T) # fitted values for model

#### bias corrected bootstrap predictions -----
hboot<- hpred[sapply(hpred, is.numeric)] # remove character columns that are all NAs (i.e., bootstrap didnt converge) 

# add in actual predicted values
heightpredict <- predmheight[, c(1,5)]
hpred1 <- cbind(heightpredict,hboot, by = "temprange_pos")
hpred <- hpred1[, -c(3)]

head(hpred) # check first and last column that converged for next line of code
hpredz <- hpred %>% 
  rowwise() %>% 
  mutate(z = qnorm(sum(c_across(V7:V5660) < pred) / length(c_across(V7:V5660)))) # proportion of bootstraps less than predicted and converted to z score
str(hpredz$z)

# calculate new probabilities at which to place bias corrected lower and upper CIs
hpredictprob <-hpredz %>%
  mutate(lower.prob = pnorm(2*z-1.96),# bias corrected lower quantile to calculate
         upper.prob = pnorm(2*z+1.96)) # bias corrected upper quantile to calculate
str(hpredictprob$lower.prob)
str(hpredictprob$upper.prob)

# bias corrected CIs
hpredictbias <- hpredictprob %>% 
  rowwise() %>% 
  transmute(temprange_pos = temprange_pos,
            pred = pred,
            lower.ci = quantile((c_across(V7:V5660)), probs = lower.prob),
            upper.ci = quantile((c_across(V7:V5660)), probs = upper.prob))

heightpredictbias <- hpredictbias

### save prediction with bias corrected 95% CI
save(heightpredictbias, file = "outputs/models/height_biascorrpredict.Rdata")

plotheight <- ggplot(heightpredictbias) +
  geom_point(data = year1sum, aes(x = temprange_pos, y = stht ), alpha = 0.6, col ="black") +
  geom_line(data = predmheight, aes(x = temprange_pos, y = pred)) +
  geom_ribbon(aes(x = temprange_pos, ymin = lower.ci, ymax = upper.ci), alpha=0.4) +
  scale_x_reverse() +
  ylab("Standardized height") + xlab("Thermal range position") + theme_classic() +
  annotate("text", x = 1.1, y=3.3, label = "P < 0.01")


################################################################################
##### ### SURVIVAL -----
################################################################################

# #### model of survival for first year ----

# need to add in weights of the denominator for the trials. which is new18 for 2018 cohort and new19 for 2019 cohort. This makes new variable "weight" which is the denominator (i.e., original number of recruits)
surv1.18 <- subset(surv1, yr1 == "prop18surv19")
surv1.18$weight <- surv1.18$new18
surv1.19 <- subset(surv1, yr1 == "prop19surv20")
surv1.19$weight <- surv1.19$new19
surv1 <- rbind(surv1.18, surv1.19)
summary(surv1$weight)

# model
msurv <- glmmTMB(surv1yr ~  temprange_pos  + pptrange_pos + yr1 +
                   (1|site/rep),  
                 data = surv1, family = binomial(link = "logit"), weights = weight) 
summary(msurv)

# model - quadratic
msurvq <- glmmTMB(surv1yr ~  poly(temprange_pos, 2)  + poly(pptrange_pos, 2) + yr1 +
                   (1|site/rep),  
                 data = surv1, family = binomial(link = "logit"), weights = weight) 
summary(msurvq)


AIC(msurv, msurvq) # quadratic relationship has lower AIC. using that

#### model diagnostics
# residual plots using DHARMa
# see https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#general-remarks-on-interperting-residual-patterns-and-tests
msurv_simres <- simulateResiduals(fittedModel = msurvq, plot = T)

testUniformity(msurv_simres) 
testOutliers(msurv_simres)
testQuantiles(msurv_simres) 

### cannot run bootstraps because always give binomial non-integer warning (which is okay)

msurv_temp <- as.data.frame(ggpredict(msurvq, terms = c("temprange_pos [all]", "pptrange_pos [0]")))
plotsurv1 <- ggplot(msurv_temp) +
  geom_line(aes(x = x, y = predicted)) + #slope
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.4) + 
  geom_point(data = surv1,                      # adding the raw data 
             aes(x = temprange_pos, y = surv1yr), alpha = 0.6) + 
  scale_x_reverse() +
  ylab("Probability surviving to year 2") + xlab("Thermal range position")+ theme_classic() +
  annotate("text", x = -1.9, y=.95, label = "Thermal: P = 0.96") +
  annotate("text", x = -1.9, y=.9, label = "Thermal²: P = 0.69")

#### model of survival to third year (from start of experiment  to year 3) -----

# model
msurv2 <- glmmTMB(prop18surv20a ~  temprange_pos  + pptrange_pos +
                    (1|site/rep),  
                  data = surv2, family = binomial(link = "logit"), weights = new18) 
summary(msurv2)

# relationship could be quadratic
msurv2q <- glmmTMB(prop18surv20a ~  poly(temprange_pos, 2)  + poly(pptrange_pos, 2) +
                     (1|site/rep),  
                   data = surv2, family = binomial(link = "logit"), weights = new18) 
summary(msurv2q)
AIC(msurv2, msurv2q) # quadratic is best

#### model diagnostics
# residual plots using DHARMa
# see https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#general-remarks-on-interperting-residual-patterns-and-tests
msurv_simres <- simulateResiduals(fittedModel = msurv2q, plot = T)

testUniformity(msurv_simres) 
testOutliers(msurv_simres) 
testQuantiles(msurv_simres)  
### CURRENTLY cannot run bootstraps because always give binomial non-integer warning (which is okay)

msurv_temp <- as.data.frame(ggpredict(msurv2q, terms = c("temprange_pos [all]", "pptrange_pos [0]")))
plotsurv2 <- ggplot(msurv_temp) +
  geom_line(aes(x = x, y = predicted)) + #slope
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.4) + 
  geom_point(data = surv2,                      # adding the raw data 
             aes(x = temprange_pos, y = prop18surv20a), alpha = 0.4) + 
  scale_x_reverse() +
  ylab("Probability surviving to year 3") + xlab("Thermal range position")+ theme_classic() +
  annotate("text", x = -1.3, y=0.95, label = "Thermal: P = 0.77") +
  annotate("text", x = -1.3, y=0.9, label = "Thermal²: P = <0.01") +
  geom_vline(xintercept = 0, linetype = 2)



#### model of survival to 2022 ----
# # sample size is low due to uncertainty from no survey in 2021.

# model
msurv22 <- glmmTMB(propxxsurv22 ~  temprange_pos  + pptrange_pos +
                     (1|site/rep),  
                   data = surv22, family = binomial(link = "logit"), weights = new_all) 
summary(msurv22)

# model - relationship could be quadratic. 
msurv22q <- glmmTMB(propxxsurv22 ~  poly(temprange_pos, 2)  + pptrange_pos +
                      (1|site/rep),  
                    data = surv22, family = binomial(link = "logit"), weights = new_all) 
summary(msurv22q)

AIC(msurv22, msurv22q) # no difference in AIC keeping linear

#### model diagnostics
# residual plots using DHARMa
# see https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#general-remarks-on-interperting-residual-patterns-and-tests
msurv_simres <- simulateResiduals(fittedModel = msurv22, plot = T)

testUniformity(msurv_simres) # looks good!
testOutliers(msurv_simres) # looks good! (NOTE if one max obs of 4.13 there is an issue here)
testQuantiles(msurv_simres)  # looks good

### CURRENTLY cannot run bootstraps because always give binomial non-integer warning (which is okay)

msurv_temp <- as.data.frame(ggpredict(msurv22, terms = c("temprange_pos [all]", "pptrange_pos [0]")))
plotsurv22 <- ggplot(msurv_temp) +
  geom_line(aes(x = x, y = predicted)) + #slope
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.4) + 
  geom_point(data = surv22,                      # adding the raw data 
             aes(x = temprange_pos, y = propxxsurv22), alpha = 0.4) + 
  scale_x_reverse() +
  ylab("Probability surviving to 2022") + xlab("Thermal range position")+ theme_classic() +
  annotate("text", x = 1, y=.9, label = "P = 0.98")



######## PANEL OF PLOTS -----
# plotheight - A
# plotsurv1 - B
# plotsurv2 - C
# plotsurv22 - D

##### panel plot 
pdf(file = "outputs/final/heightsurvpred.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 7)
ggdraw() +
  draw_plot(plotheight, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(plotsurv1, x = 0.5, y = .5, width = .5, height = .5) +
  draw_plot(plotsurv2, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(plotsurv22, x = 0.5, y =0, width = .5, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))
dev.off()



