################################################################################
################################################################################
###                                                                          ###
###                            ANALYSIS SCRIPT 3:                            ###
###                 RECRUITMENT WITH RANGE POSITION MODELLING                ###
###                      OPTIMA WITH CHANGING CANOPY COVER                   ###
###                                                                          ###
################################################################################
################################################################################

# Created Dec. 9, 2022
# Last Updated Dec. 2, 2024
# Authors: Katie Goodwin

# Analyses for Goodwin et al. 2024 Lagged climate-driven range shifts at species’ leading, but not trailing, range edges revealed by multispecies seed addition experiment. Ecography. DOI: 10.1111/ecog.07331
# testing whether optimal seedling recruitment is different with changing canopy cover 
# also code here to confirm there was in fact a microclimatic buffering effect due to canopy cover (Fig. S3)
# tests the microclimatic buffering hypothesis at the optima

###########################################################################################
#### OUTPUTS FROM THIS SCRIPT ------------------------------------------------
###########################################################################################
load("outputs/models/canopy_optima.Rdata") # data to make canopy optima histogram
load("outputs/models/canopy_biascorrpredict.Rdata") # predicted recruitment with temp range position for different canopy covers with bias-corrected CIs

# these data frames were then inputed into scripts/outputscripts/pretty_plots.R to generate manuscript figs.


###########################################################################################
#### PACKAGES AND DATA ------------------------------------------------
###########################################################################################
library(glmmTMB)
library(tidyverse)
library(ggplot2)
library(GGally)
library(ggeffects)
library(DHARMa)

# count census data 
load('data/recruit_census.RData') 
load('data/microclimate.RData') # environmental data for each replicate plot (from Chardon et al. 2024 Ecography paper. The paired sister paper to this one on microhabitat. You should check it out! https://doi.org/10.1111/ecog.07144 )

str(census)
str(census.rep)





################################################################################
##### Does canopy cover actually have a microclimatic buffering effect? - using Nathalies TOMST loggers -----
################################################################################
# see Chardon et al 2024 https://doi.org/10.1111/ecog.07144 for detailed exploration of microscale recruitment patterns with these data
# I use these data here just to confirm that was a microclimatic buffering affect due to canopy cover

#### PREP DATA
## add canopy cover open or closed site 
repsum <- census %>%
  group_by(rep, canopy) %>%
  summarise()

rep.env <- cbind(census.rep, repsum)
rep.env <- rep.env[, -c(7)]

# summarize canopy cover to site level because that is the level the microclimatic data are at
rep.envsum <- rep.env %>%
  group_by(site, tmoist_avg, tmoist_range, t1_avg, t1_range, t3_avg, t3_range) %>%
  summarize(canopymean = mean(canopycont))

# convert canopy cover (currently % sky but want to change to % canopy)
rep.envsum$canopycont1 <- 100 - rep.envsum$canopymean

hist(rep.envsum$canopycont1)


### PLOTS OF DIFFERENT MICROCLIMATE VARIABLES 
# average  Soil moisture 
# Conclusion: soil moisture doesnt change with canopy cover
a <- ggplot(rep.envsum, aes(x = canopycont1, y = tmoist_avg)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', col = "black") +
  theme_classic() +
  xlab("% Canopy Cover") + ylab("Average Soil Moisture") +
  annotate("text", x = 15, y=.45, label = "P = 0.80", col = "black")
a

ma <- lm(canopycont1 ~ tmoist_avg, data = rep.envsum)
summary(ma)

# Soil moisture range
# Conclusion: soil moisture range doesnt change with canopy cover
b <- ggplot(rep.envsum, aes(x = canopycont1, y = tmoist_range)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', col = "black") +
  theme_classic() +
  xlab("% Canopy Cover") + ylab("Average Soil Moisture Range") +
  annotate("text", x = 15, y=.4, label = "P = 0.84", col = "black")
b

mb <- lm(canopycont1 ~ tmoist_range, data = rep.envsum)
summary(mb)

# Soil temp 
# Conclusion: soil temp decreases with canopy cover
c <- ggplot(rep.envsum, aes(x = canopycont1, y = t1_avg)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', col = "black") +
  theme_classic() +
  xlab("% Canopy Cover") + ylab("Average Soil Temperature (°C)") +
  annotate("text", x = 15, y=20, label = "P = 0.07", col = "black")
c

mc <- lm(canopycont1 ~ t1_avg, data = rep.envsum)
summary(mc)

# Soil temp range 
# no change with canopy cover
d <- ggplot(rep.envsum, aes(x = canopycont1, y = t1_range)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', col = "black") +
  theme_classic() +
  xlab("% Canopy Cover") + ylab("Average Soil Temperature Range (°C)") +
  annotate("text", x = 15, y=25, label = "P = 0.75", col = "black")
d

  md <- lm(canopycont1 ~ t1_range, data = rep.envsum)
  summary(md)


# Above ground temp
# Conclusion: decreases with canopy cover

e <- ggplot(rep.envsum, aes(x = canopycont1, y = t3_avg)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', col = "black") +
  theme_classic() +
  xlab("% Canopy Cover") + ylab("Average Plant Height Temperature (°C)") +
  annotate("text", x = 15, y=19.5, label = "P = 0.08", col = "black")
e

me <- lm(canopycont1 ~ t3_avg, data = rep.envsum)
summary(me)


# Above ground temp range 
# Conclusion: decreases with canopy cover
f <- ggplot(rep.envsum, aes(x = canopycont1, y = t3_range)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', col = "black") +
  theme_classic() +
  xlab("% Canopy Cover") + ylab("Average Plant Height Temperature Range (°C)")+
  annotate("text", x = 15, y=52, label = "P = < 0.01", col = "black")
f

mf <- lm(canopycont1 ~ t3_range, data = rep.envsum)
summary(mf)

## supp plot for MS showing microclimatic effects of canopycover - Fig. S3
pdf(file = "outputs/final/canopymicroclim.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 8)
ggdraw() +
  draw_plot(e, x = 0, y = 0.5, width = .33, height = .5) +
  draw_plot(c, x = 0.33, y = 0.5, width = .33, height = .5) +
  draw_plot(a, x = 0.67, y = 0.5, width = .33, height = .5) +
  draw_plot(f, x = 0, y = 0, width = .33, height = .5) +
  draw_plot(d, x = 0.33, y = 0, width = .33, height = .5) +
  draw_plot(b, x = .67, y = 0, width = .33, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), size = 15,
                  x = c(0, 0.33, 0.67, 0, 0.33, 0.67), y = c(1, 1, 1, 0.5, 0.5, 0.5))
dev.off()





################################################################################
### testing microrefugia hypothesis ----------
################################################################################
##### PREP DATA SET -----
# remove species that recruited in less than 5 plots - smallest number of species to be able to include interaction terms
censussuc <- census %>% 
  group_by(species) %>% 
  filter(species != "CARSPE" & species != "MAIDIL" & species !="MAIRAC" & species != "PINCON") 

# species that recruited in at least 5 plots
census5 <- censussuc %>% #already subseted to species that recruited at least once
  group_by(species) %>% 
  filter(species != "CARSTI" & species != "SAMCER" & species !="RUBSPE" & species != "PINPON" & species != "SAMRAC" &
           species != "CARSPE" & species != "MAIDIL" & species !="MAIRAC" & species != "PINCON") 

# remove NAs
census5a <- subset(census5, new_all != "NA")
census5 <- census5a


#### mcan: recruitment with range position with an interaction with canopy cover
# For model to converge removed pptrange_pos and added in region instead to account for some regional ppt effects
mcan <- glmmTMB(new_all ~  poly(temprange_pos, 2)  + canopycont + temprange_pos:canopycont + region
              + (1|site/rep) + (1|family), 
              ziformula = ~ poly(temprange_pos, 2)  + canopycont
              + (1|site/rep) + (1|family),
              data = census5, 
              family = nbinom2(link = "log")) 
summary(mcan) # no interaction effect between temprange_pos and canopy cover 

# # save model output
# save(mcan, file = "outputs/range_position/models/canopy_model.Rdata")
# load("outputs/range_position/models/canopy_model.Rdata")

#### what is the mean open and mean closed canopy cover
open <- subset(census, canopy == "open")
summary(open$canopycont) #35.979 %

closed <- subset(census, canopy == "closed")
summary(closed$canopycont) #12.034 %

# model diagnostics
# residual plots using DHARMa
# see https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#general-remarks-on-interperting-residual-patterns-and-tests
mcan_simres <- simulateResiduals(fittedModel = mcan, plot = T)

testUniformity(mcan_simres)
testOutliers(mcan_simres) 
testQuantiles(mcan_simres) 

##### BOOTSTRAPPING  ---------------
# i.e., do 95% confidence intervals of open and closed peak overlap?
# separate code file as this loop takes >24 hours and is run on the zoology cluster
# (see canopybootstrap_cluster.R)

###### import bootstrap runs from cluster server (note ~ half didnt converge so ran twice to get at least 10000 runs)

# bootstraps (file too big for git hub adding manually from Katies computer. not in repo)). The bootstrapped datasets can be generated by running the canopybootstrap_cluster.R script

# load("~/OneDrive - UBC/PhD/Cascades Legacy Collab/canopy_bootstraps/canopy/canopy_cond_optimaCI06OCT23.RData")
# load("~/OneDrive - UBC/PhD/Cascades Legacy Collab/canopy_bootstraps/canopy/canopy_cond_bootstraps06OCT23.RData")

canopt <- optima
canpred <- pred

# fix classes
canopt$opt_12 <- as.numeric(as.character(canopt$opt_12)) # 12% sky or 88% canopy cover i.e., mean of closed canopy sites
canopt$opt_36 <- as.numeric(as.character(canopt$opt_36)) # 36% sky or 64 % canopy cover (i.e., mean of open canopy sites)
canopt$opt_50 <- as.numeric(as.character(canopt$opt_50)) # 50% canopy cover
canopt$opt_75 <- as.numeric(as.character(canopt$opt_75)) # 25 % canopy cover
str(canopt)

##### CI FOR PEAK OF RECRUITMENT CURVE -----

# temp - cond - optima gets cooler with larger spread with increasingly open canopies. Less data for more open canopies hence the larger spread (note that canopycont is % sky. Any plot mentioning this is reversed for the MS so in the MS it is percent canopy cover)
hist(canopt$opt_12, breaks = 40)
hist(canopt$opt_36, breaks = 40)
hist(canopt$opt_50, breaks = 40)
hist(canopt$opt_75, breaks = 40)

quantile(canopt$opt_12, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE) 
quantile(canopt$opt_36, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE) 
quantile(canopt$opt_50, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE) 
quantile(canopt$opt_75, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE) 

# plotting both on one histogram (just 12 and 36)
toptplot1 <- canopt[, c(1,2)]
toptplot <- toptplot1 %>%
  pivot_longer(c(opt_12, opt_36), names_to = "canopy", values_to = "opt_temprange")
toptplot$canopy <- as.factor(toptplot$canopy)
toptplot$canopy <- plyr::revalue(toptplot$canopy, c(opt_12 = "mean closed", opt_36 = "mean open"))

# save optima dataframe to make pretty plot elsewhere
toptcanopy <- toptplot
save(toptcanopy, file = "outputs/models/canopy_optima.Rdata")


##### BIAS-CORRECTED BOOTSTRAP PREDICTIONS ----
# model predictions
predmcan <- expand.grid(temprange_pos = seq(-2.6, 1.6, by = .01), canopycont = c(12,36, 50, 75), region = "RP",
                         site = "NA", rep = "NA", family = "NA" ) 
predmcan$pred <- predict(mcan, newdata = predmcan, type = "cond", se = F,  allow.new.levels = T) # fitted values for model

canboot<- canpred[sapply(canpred, is.numeric)] # remove character columns that are all NAs (i.e., bootstrap didnt converge) 

# add in actual predicted values
canpredict <- predmcan[, c(1,7)]
canpred1 <- cbind(canpredict,canboot)
canpred <- canpred1[, -c(3)]

head(canpred) # check first and last column that converged for next line of code
# proportion of bootstraps less than predicted and converted to z score - long running time
canpredz <- canpred %>% 
  rowwise() %>% 
  mutate(z = qnorm(sum(c_across(V7:V20005) < pred) / length(c_across(V7:V20005))))
str(canpredz$z)

# calculate new probabilities at which to place bias corrected lower and upper CIs
canpredictprob <-canpredz %>%
  mutate(lower.prob = pnorm(2*z-1.96),# bias corrected lower quantile to calculate
         upper.prob = pnorm(2*z+1.96)) # bias corrected upper quantile to calculate
str(canpredictprob$lower.prob)
str(canpredictprob$upper.prob)

# bias corrected CIs
canpredictbias <- canpredictprob %>% 
  rowwise() %>% 
  transmute(temprange_pos = temprange_pos,
            pred = pred,
            canopycont = canopycont,
            lower.ci = quantile((c_across(V7:V20005)), probs = lower.prob),
            upper.ci = quantile((c_across(V7:V20005)), probs = upper.prob)) # long running time

### save prediction with bias corrected 95% CI
save(canpredictbias, file = "outputs/models/canopy_biascorrpredict.Rdata")
