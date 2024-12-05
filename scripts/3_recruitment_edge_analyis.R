################################################################################
################################################################################
###                                                                          ###
###                            ANALYSIS SCRIPT 2:                            ###
###                 RECRUITMENT WITH RANGE POSITION MODELLING                ###
###                      RECRUITMENT AT RANGE EDGES                          ###
###                                                                          ###
################################################################################
################################################################################


# Created Nov 10, 2022
# Last Updated Dec. 2, 2024
# Authors: Katie Goodwin

# Analyses for Goodwin et al. 2024 Lagged climate-driven range shifts at species’ leading, but not trailing, range edges revealed by multispecies seed addition experiment. Ecography. DOI: 10.1111/ecog.07331
# testing recruitment within and beyond range edges (Fig. 5 for temperature and Fig. S7 for precipitation)
# includes tests for both the lagged responses hypothesis and microclimatic buffering hypothesis at the edges

###########################################################################################
#### OUTPUTS FROM THIS SCRIPT ------------------------------------------------
###########################################################################################
load("outputs/models/coldedge_biascorrpredict.Rdata") # predicted cold edge recruitment with bias-corrected CIs
load("outputs/models/hotedge_biascorrpredict.Rdata") # predicted warm edge recruitment with bias-corrected CIs
load("outputs/models/wetedge_biascorrpredict.Rdata") # predicted wet edge recruitment with bias-corrected CIs
load("outputs/models/dryedge_biascorrpredict.Rdata") # predicted dry edge recruitment with bias-corrected CIs

# "outputs/final/highgbiftest.pdf" = Fig. S4 - sensitivity tests for uncertainty in occurrence records
# these data frames were then inputed into scripts/outputscripts/pretty_plots.R to generate manuscript figs.

###########################################################################################
#### PACKAGES AND DATA ------------------------------------------------
###########################################################################################
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(ggeffects)


# seedling count census data 
load('data/recruit_census.RData') 
str(census)


###########################################################################################
#### RECRUITMENT AT AND BEYOND COLD EDGE ----
###########################################################################################
##### GENERATE DATASET ----
# of species that successfully recruited and had sites beyond cold edge with binned thermal range position of cold half range of beyond cold edge

# subset to species that successfully recruited
speciesplots <- census %>%
  group_by(species) %>%
  summarize(plotsgerm = sum(new_all, na.rm = TRUE))
# CARSPE, MAIDIL, MAIRAC, PINCON have zero recruitment in any plots
censussuc <- census %>% 
  group_by(species) %>% 
  filter(species != "CARSPE" & species != "MAIDIL" & species !="MAIRAC" & species != "PINCON") 

# subset to sites in cool half range or beyond cold edge (i.e., temprange_pos <0)
leading1 <- subset(censussuc, temprange_pos < 0)

# subset to species that successfully recruited and had sites beyond leading edges
beyond <- subset(leading1, temprange_pos < -1.1) # a little smaller than -1 to properly be beyond range
unique(beyond$species) 
# ABIGRA PICSIT PINPON MAHNER RUBURS SAMCER VACPAR TOLMEN CARSTI ERILAN MAHAQU RUBSPE TELGRA (13 species)
# ^^^ species that recruited successfully at least once and had sites beyond cold range edge

leading <- leading1 %>% 
  group_by(species) %>% 
  filter(species == "ABIGRA" | species == "PICSIT" | species == "PINPON" |
           species == "MAHNER" | species == "RUBURS" | species == "SAMCER" |
           species == "VACPAR" | species == "TOLMEN" | species == "CARSTI" |
           species == "ERILAN" | species == "MAHAQU" | species == "RUBSPE" | 
           species == "TELGRA")
unique(leading$species) # all these ^ species are in here! 

# creating thermal range position bins: within = cold half of range, beyond = beyond cold edge
hist(leading$temprange_pos, breaks = 50)
within <- subset(leading, temprange_pos > -1)
beyond <- subset(leading, temprange_pos < -1)


within$temprange_fac <- "within"
beyond$temprange_fac <- "beyond"

leading <- rbind(within, beyond)
leading$temprange_fac <- as.factor(leading$temprange_fac)

# remove values where new_all = NA (i.e., was more recruitment in control than seed addition plot)
leadinga <- subset(leading, new_all != "NA")
leading <- leadinga

# save data frame
save(leading, file = "outputs/models/coolhypodata.RData")

##### mcold: model comparing recruitment within and beyond cold edge ----
# model only runs with a single zero-inflation parameter applied to all observations (ziformula = ~1)
# but diagnostics look good and this is acceptable to do just cant give inferences of variables driving prob. of recruitment
# reference: https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf
mcold <- glmmTMB(new_all ~  temprange_fac + poly(pptrange_pos, 2)
                 + (1|site/rep) + (1|family), 
                 ziformula = ~ 1,
                 data = leading, 
                 family = nbinom2(link = "log")) 
summary(mcold)

# save model output 
# save(mcold, file = "outputs/models/coldedge_model.Rdata")


### MODEL DIAGNOSTICS 

# residual plots using DHARMa
m1_simres <- simulateResiduals(mcold)
testUniformity(m1_simres) 
testOutliers(m1_simres)
testQuantiles(m1_simres)

##### BOOTSTRAPPING  ---------------
# i.e., to create confidence intervals for prediction line
# separate code file as this loop takes >24 hours and is run on the zoology cluster

###### import bootstrap runs from cluster server 
load('outputs/bootstraps/coldedge/coldedge_pred17JAN24.RData')
cpred <- pred

##### MODEL PREDICTIONS ----
predmcold <- expand.grid(temprange_fac = factor(levels(leading$temprange_fac)), pptrange_pos = 0, 
                           site = "NA", rep = "NA", family = "NA") # for predicted values. will take pop means for random effects. Since no interaction terms value of ppt shouldnt matter
predmcold$pred <- predict(mcold, newdata = predmcold, type = "conditional", se = F,  allow.new.levels = T) # fitted values for model

##### BIAS CORRECTED BOOTSTRAP CIS -----
# code adapted from doi: 10.3389/fpsyg.2022.810258, Sheth and Angert 2018 

cboot <- cpred[sapply(cpred, is.numeric)] # remove character columns that are all NAs (i.e., bootstrap didnt converge) 

# add in actual predicted values
cpredict1 <- predmcold[, c(1,6)]
cpred <- cbind(cpredict1, cboot)

head(cpred) # check first and last column that converged for next line of code
# proportion of bootstraps less than predicted and converted to z score
cpredz <- cpred %>% 
  rowwise() %>% 
  mutate(z = qnorm(sum(c_across(V9:V20002) < pred) / length(c_across(V9:V20002)))) # z
cpredz$z

# calculate new probabilities at which to place bias corrected lower and upper CIs
cpredictprob <- cpredz %>%
  mutate(lower.prob = pnorm(2*z-1.96),# bias corrected lower quantile to calculate
         upper.prob = pnorm(2*z+1.96)) # bias corrected upper quantile to calculate

# bias corrected CIs
cpredictbias <- cpredictprob %>% 
  rowwise() %>% 
  transmute(temprange_fac = temprange_fac,
            pred = pred,
            lower.ci = quantile((c_across(V9:V20002)), probs = lower.prob),
            upper.ci = quantile((c_across(V9:V20002)), probs = upper.prob))

### save prediction with bias corrected 95% CI
save(cpredictbias, file = "outputs/models/coldedge_biascorrpredict.Rdata")

##### SENSITIVITY TEST - ADDRESSING UNCERTAINTY IN OCCURRENCE RECORDS ----
## some species did have uncertainty in their thermal range b/c few occurrence records. same conclusion excluding those species?
# CONCLUSION: doesn't affect inferences at all. Still successful, but reduced recruitment beyond the cold edge

# subset to species with >50 occurrence records (removing) ABIGRA, PICSIT, RUBSPE, TELGRA, TOLMEN
leading.highgbif <- leading %>% 
  group_by(species) %>% 
  filter(species == "PINPON" |
           species == "MAHNER" | species == "RUBURS" | species == "SAMCER" |
           species == "VACPAR" | species == "CARSTI" |
           species == "ERILAN" | species == "MAHAQU" )

# full model doesnt converge. removed rep
mcoldtest <- glmmTMB(new_all ~  temprange_fac + poly(pptrange_pos, 2) 
                     + (1|site) + (1|family), 
                     ziformula = ~ 1,
                     data = leading.highgbif, 
                     family = nbinom2(link = "log")) 
summary(mcoldtest) # now there is significantly reduced recruitment beyond cold edge, but still recruitment occurred. Same conclusion as full model

mcold_temp <- ggpredict(mcoldtest, terms = c("temprange_fac [all]", "pptrange_pos[0]"), type = "fe")

# save ggplot of this
coldtest <- ggplot(mcold_temp) +
  geom_violin(data = leading.highgbif, aes(x = temprange_fac, y = log(new_all+1)), alpha = 0.6, col = "gray25") +
  geom_point(aes(x = x, y = log(predicted +1)), size = 1, col = "black") +
  geom_linerange(aes(x = x, ymin = log(conf.low + 1), ymax = log(conf.high + 1)), size = .5, col = "black") +
  ylab("ln(Seedling Count + 1)") + xlab("Thermal range position") + 
  scale_x_discrete(labels = c('Beyond cold edge', 'Cool half range')) +
  theme_classic() +
  theme(legend.position = "none")
coldtest


##### COLD EDGE RECRUITMENT - INTERACTION WITH CANOPY COVER? -----
# for testing the microclimatic buffering hypothesis
mcoldcan <- glmmTMB(new_all ~  temprange_fac + poly(pptrange_pos, 2) + temprange_fac*canopycont 
                    + (1|site/rep) + (1|family), 
                    ziformula = ~ 1,
                    data = leading, 
                    family = nbinom2(link = "log")) 
summary(mcoldcan)
# canopy cover interaction p 0.86



###########################################################################################
#### RECRUITMENT AT AND BEYOND WET EDGES ----
###########################################################################################
##### GENERATE DATASET ----
# dataset of species that successfully recruited and had sites beyond wet edge

# subset to species that successfully recruited
speciesplots <- census %>%
  group_by(species) %>%
  summarize(plotsgerm = sum(new_all, na.rm = TRUE))
# CARSPE, MAIDIL, MAIRAC, PINCON have zero recruitment in any plots
censussuc <- census %>% 
  group_by(species) %>% 
  filter(species != "CARSPE" & species != "MAIDIL" & species !="MAIRAC" & species != "PINCON") 

# subset to sites in wet half of range 
wetedge1 <- subset(censussuc, pptrange_pos > 0)

# subset to species that successfully recruited and had sites beyond wet edges
wetbeyond <- subset(wetedge1, pptrange_pos > 1.1) 
unique(wetbeyond$species) 
# ABIGRA CARSTI PINPON SAMCER ERILAN PICENG ABILAS SAMRAC (8 species)
# ^^^ recruited successfully at least once and had sites beyond range

wetedge <- wetedge1 %>% 
  group_by(species) %>% 
  filter(species == "ABIGRA" | species == "CARSTI" | species == "PINPON" |
           species == "SAMCER" | species == "ERILAN" | species == "PICENG" |
           species == "ABILAS" | species == "SAMRAC")
unique(wetedge$species) # all these ^ species are in here! 


### create bins of within and beyond wet edge as factor
withinwet <- subset(wetedge, pptrange_pos < 1)
beyondwet <- subset(wetedge, pptrange_pos > 1) 

# combining
withinwet$pptrange_fac <- "within"
beyondwet$pptrange_fac <- "beyond"

leadingwet <- rbind(withinwet, beyondwet)
leadingwet$pptrange_fac <- as.factor(leadingwet$pptrange_fac)

# remove values where new_all = NA (i.e., was more recruitment in control than seed addition plot)
leadingweta <- subset(leadingwet, new_all != "NA")
leadingwet <- leadingweta

# save
save(leadingwet, file = "outputs/models/wethypodata.RData")


##### mwet: model comparing recruitment within and beyond wet edge ------

# removing rep due to model convergence problems - false convergence warning 
mwet <- glmmTMB(new_all ~  pptrange_fac + poly(temprange_pos, 2)  
                 + (1|site) + (1|family), 
                 ziformula = ~ 1,
                 data = leadingwet, 
                 family = nbinom2(link = "log")) 
summary(mwet) # no difference within and beyond wet edge

# # save model
# save(mwet, file = "outputs/range_position/models/wetmodel.RData")

### MODEL DIAGNOSTICS
# residual plots using DHARMa
m1_simres <- simulateResiduals(mwet)
testUniformity(m1_simres) # looks good
testOutliers(m1_simres) # looks good
testQuantiles(m1_simres) # looks good

##### BOOTSTRAPPING  ---------------
# i.e., to create confidence intervals for prediction line
# separate code file as this loop takes >24 hours and is run on the zoology cluster

###### import bootstrap runs from cluster server 
load('outputs/bootstraps/wetedge/wetedge_pred17JAN24.RData')
wetpred <- pred

##### MODEL PREDICTIONS  ----
predmwet <- expand.grid(pptrange_fac = factor(levels(leadingwet$pptrange_fac)), temprange_pos = 0, 
                        site = "NA", rep = "NA", family = "NA") # for predicted values. will take pop means for random effects. Since no interaction terms shouldnt matter
predmwet$pred <- predict(mwet, newdata = predmwet, type = "conditional", se = F,  allow.new.levels = T) # fitted values for model

##### BIAS CORRECTED BOOTSTRAP CIS -----
wetboot<- wetpred[sapply(wetpred, is.numeric)] # remove character columns that are all NAs (i.e., bootstrap didnt converge) 

# add in actual predicted values
wetpredict1 <- predmwet[, c(1,6)]
wetpred <- cbind(wetpredict1, wetboot)

head(wetpred) # check first and last column that converged for next line of code
# proportion of bootstraps less than predicted and converted to z score
wetpredz <- wetpred %>% 
  rowwise() %>% 
  mutate(z = qnorm(sum(c_across(V9:V10004) < pred) / length(c_across(V9:V10004))))
wetpredz$z

# calculate new probabilities at which to place bias corrected lower and upper CIs
wetpredictprob <- wetpredz %>%
  mutate(lower.prob = pnorm(2*z-1.96),# bias corrected lower quantile to calculate
         upper.prob = pnorm(2*z+1.96)) # bias corrected upper quantile to calculate
wetpredictprob$lower.prob
wetpredictprob$upper.prob

# bias corrected CIs
wetpredictbias <- wetpredictprob %>% 
  rowwise() %>% 
  transmute(pptrange_fac = pptrange_fac,
            pred = pred,
            lower.ci = quantile((c_across(V9:V10004)), probs = lower.prob),
            upper.ci = quantile((c_across(V9:V10004)), probs = upper.prob))

### save prediction with bias corrected 95% CI
save(wetpredictbias, file = "outputs/models/wetedge_biascorrpredict.Rdata")

###########################################################################################
#### RECRUITMENT AT AND BEYOND WARM EDGES ----
###########################################################################################
##### GENERATE DATASET ----
# subset to species that successfully recruited
speciesplots <- census %>%
  group_by(species) %>%
  summarize(plotsgerm = sum(new_all, na.rm = TRUE))
# CARSPE, MAIDIL, MAIRAC, PINCON have zero recruitment in any plots
censussuc <- census %>% 
  group_by(species) %>% 
  filter(species != "CARSPE" & species != "MAIDIL" & species !="MAIRAC" & species != "PINCON") 

# subset to sites in warm half of range (i.e., temprange_pos >0)
trailing1 <- subset(censussuc, temprange_pos > 0)

# subset to species that successfully recruited and had sites beyond presumed seedling warm limit
# this value is based off of the recruitment optimum shift where we assume the seedling warm limit would shift the same amount (see manuscript for description 1 - 0.26 = 0.74)
hot <- subset(trailing1, temprange_pos > 0.74) 
unique(hot$species) 
# SORSIT ABILAS PICENG VACDEL ERIPER ANEOCC (6 species)
# ^^^ recruited successfully at least once and had sites that went beyond at least half of the trailing range edge

trailing <- trailing1 %>% 
  group_by(species) %>% 
  filter(species == "SORSIT" | species == "ABILAS" | species == "PICENG" |
           species == "VACDEL" | species == "ERIPER" | species == "ANEOCC" )
unique(trailing$species) # all these ^ species are in here! 

# bin to within and beyond this presumed seedling warm edge
newbeyond <- subset(trailing, temprange_pos > 0.74)
newwithin <- subset(trailing, temprange_pos < 0.74) 

# combining
newbeyond$temprange_fac <- "newbeyond"
newwithin$temprange_fac <- "newwithin"

trailing <- rbind(newwithin, newbeyond)
trailing$temprange_fac <- as.factor(trailing$temprange_fac)

# remove values where new_all = NA (i.e., was more recruitment in control than seed addition plot)
trailinga <- subset(trailing, new_all != "NA")
trailing <- trailinga

# save data frame
save(trailing, file = "outputs/models/warmhypodata.RData")


##### mhot:  comparing recruitment within and beyond presumed seedling warm edge  -----
mhot <- glmmTMB(new_all ~  temprange_fac + poly(pptrange_pos, 2) 
                + (1|site/rep) + (1|family), 
                ziformula = ~ 1,
                data = trailing, 
                family = nbinom2(link = "log")) 
summary(mhot)

# # save model
# save(mhot, file = "outputs/range_position/models/warmedge_model.Rdata")


### MODEL DIAGNOSTICS 
# residual plots using DHARMa
m1_simres <- simulateResiduals(mhot)
testUniformity(m1_simres)
testOutliers(m1_simres)
testQuantiles(m1_simres)

##### BOOTSTRAPPING  ---------------
# i.e., to create confidence intervals for prediction line
# separate code file as this loop takes >24 hours and is run on the zoology cluster

###### import bootstrap runs from cluster server 
load('outputs/bootstraps/warmedge/warmedge_pred17JAN24.RData')
hotpred <- pred

##### MODEL PREDICTIONS ----
predmhot <- expand.grid(temprange_fac = factor(levels(trailing$temprange_fac)), pptrange_pos = 0, 
                         site = "NA", rep = "NA", family = "NA") 
predmhot$pred <- predict(mhot, newdata = predmhot, type = "conditional", se = F,  allow.new.levels = T) # fitted values for model

###### BIAS CORRECTED BOOTSTRAP CIS -----

hotboot<- hotpred[sapply(hotpred, is.numeric)] # remove character columns that are all NAs (i.e., bootstrap didnt converge) 

# add in actual predicted values
hotpredict1 <- predmhot[, c(1,6)]
hotpred <- cbind(hotpredict1, hotboot)

head(hotpred) # check range of columns that converged for next line of code
# proportion of bootstraps less than predicted and converted to z score
hotpredz <- hotpred %>% 
  rowwise() %>% 
  mutate(z = qnorm(sum(c_across(V9:V5500) < pred) / length(c_across(V9:V5500))))
hotpredz$z

# calculate new probabilities at which to place bias corrected lower and upper CIs
hotpredictprob <- hotpredz %>%
  mutate(lower.prob = pnorm(2*z-1.96),# bias corrected lower quantile to calculate
         upper.prob = pnorm(2*z+1.96)) # bias corrected upper quantile to calculate
hotpredictprob$lower.prob
hotpredictprob$upper.prob

# bias corrected CIs
hotpredictbias <- hotpredictprob %>% 
  rowwise() %>% 
  transmute(temprange_fac = temprange_fac,
            pred = pred,
            lower.ci = quantile((c_across(V9:V5500)), probs = lower.prob),
            upper.ci = quantile((c_across(V9:V5500)), probs = upper.prob))

### save prediction with bias corrected 95% CI
save(hotpredictbias, file = "outputs/models/hotedge_biascorrpredict.Rdata")

##### SENSITIVITY TEST - ADDRESSING UNCERTAINTY IN OCCURRENCE RECORDS -----
## some species did have uncertainty in their thermal range b/c low adult occurrence records. same conclusion excluding those species?
# CONCLUSION: doesn't affect inferences at all. Still not significant difference in recruitment within or beyond the range

# subset to species with >50 adult occurrence records (removing) ANEOCC, VACDEL
trailing.highgbif <- trailing %>% 
  group_by(species) %>% 
  filter(species == "SORSIT" | species == "ABILAS" | species == "PICENG" |
           species == "ERIPER")

mhottest <- glmmTMB(new_all ~  temprange_fac + poly(pptrange_pos, 2) 
                    + (1|site/rep) + (1|family), 
                    ziformula = ~ 1,
                    data = trailing.highgbif, 
                    family = nbinom2(link = "log")) 
summary(mhottest)

mhot_temp <- ggpredict(mhottest, terms = c("temprange_fac [all]", "pptrange_pos[0]"), type = "fe")

# save ggplot of this
hottest <- ggplot(mhot_temp) +
  geom_violin(data = trailing.highgbif, aes(x = temprange_fac, y = log(new_all+1)), alpha = 0.6, col = "gray25") +
  geom_point(aes(x = x, y = log(predicted +1)), size = 1, col = "black") +
  geom_linerange(aes(x = x, ymin = log(conf.low + 1), ymax = log(conf.high + 1)), size = .5, col = "black") +
  ylab("ln(Seedling Count + 1)") + xlab("Thermal range position") + 
  scale_x_discrete(labels = c('Presumed unsuitably warm', 'Warm half range')) +
  theme_classic() +
  theme(legend.position = "none")
hottest

# combine all the high adult occurrence tests (note optima plot is in recruitment_optima_analysis.R and saved as object opttest)
library(cowplot)
pdf(file = "outputs/final/highgbiftest.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 4)
ggdraw() +
  draw_plot(opttest, x = 0, y = 0, width = .33, height = 1) +
  draw_plot(coldtest, x = .33, y = 0, width = .33, height = 1) +
  draw_plot(hottest, x = .66, y = 0, width = .33, height = 1) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0.005, 0.305, 0.635), y = c(1, 1, 1))
dev.off()


##### WARM EDGE RECRUITMENT - INTERACTION WITH CANOPY COVER? ----
mhotcan <- glmmTMB(new_all ~  temprange_fac + poly(pptrange_pos, 2) + canopycont + canopycont*temprange_fac
                   + (1|site/rep) + (1|family), 
                   ziformula = ~ 1,
                   data = trailing, 
                   family = nbinom2(link = "log")) 
summary(mhotcan) # canopy cover * temprange pos interaction term p = 0.64
# canopy cover is not important

###########################################################################################
#### RECRUITMENT AT AND BEYOND DRY EDGES ----
###########################################################################################
##### GENERATE DATASET ----
# subset to species that successfully recruited
speciesplots <- census %>%
  group_by(species) %>%
  summarize(plotsgerm = sum(new_all, na.rm = TRUE))
# CARSPE, MAIDIL, MAIRAC, PINCON have zero recruitment in any plots
censussuc <- census %>% 
  group_by(species) %>% 
  filter(species != "CARSPE" & species != "MAIDIL" & species !="MAIRAC" & species != "PINCON") 

# subset to sites in dry half of range 
dryedge1 <- subset(censussuc, pptrange_pos < 0)

# subset to species that successfully recruited and had sites beyond dry edges
drybeyond <- subset(dryedge1, pptrange_pos < -1.1) # a little smaller than -1 to properly be beyond range
unique(drybeyond$species) 
# MAHNER RUBURS SORSIT VACPAR TOLMEN VACDEL TELGRA (7 species)
# ^^^ recruited successfully at least once and had sites beyond range

dryedge <- dryedge1 %>% 
  group_by(species) %>% 
  filter(species == "MAHNER" | species == "RUBURS" | species == "SORSIT" |
           species == "VACPAR" | species == "TOLMEN" | species == "VACDEL" |
           species == "TELGRA")
unique(dryedge$species) # all these ^ species are in here! 

# within and beyond range as factor. 0.88 is the presumed seedling dry limit based on the optima ppt shift to the wet part of the range 
withindry <- subset(dryedge, pptrange_pos > - 0.88)
beyonddry <- subset(dryedge, pptrange_pos < -0.88) 

withindry$pptrange_fac <- "within"
beyonddry$pptrange_fac <- "beyond"

trailingdry <- rbind(withindry, beyonddry)
trailingdry$pptrange_fac <- as.factor(trailingdry$pptrange_fac)

# remove values where new_all = NA (i.e., was more recruitment in control than seed addition plot)
trailingdrya <- subset(trailingdry, new_all != "NA")
trailingdry <- trailingdrya

 # save as dataframe
 save(trailingdry, file = "outputs/models/dryhypodata.RData")

##### mdry: comparing recruitment within and beyond presumed seedling dry edge ----------------
# remove rep for model to converge
mdry <- glmmTMB(new_all ~  pptrange_fac + poly(temprange_pos, 2)  
                + (1|site) + (1|family), 
                ziformula = ~ 1,
                data = trailingdry, 
                family = nbinom2(link = "log")) 
summary(mdry)

# # save model
# save(mdry, file = "outputs/range_position/models/drymodel.RData")
# load("outputs/range_position/models/drymodel.RData")


### MODEL DIAGNOSTICS
# residual plots using DHARMa
m1_simres <- simulateResiduals(mdry)
testUniformity(m1_simres)
testOutliers(m1_simres) 
testQuantiles(m1_simres)

##### BOOTSTRAPPING  ------
# i.e., to create confidence intervals for prediction line
# separate code file as this loop takes >24 hours and is run on the zoology cluster

###### import bootstrap runs from cluster server 
load('outputs/bootstraps/dryedge/dryedge_pred17JAN24.RData')
drypred <- pred

##### MODEL PREDICTIONS ----
predmdry <- expand.grid(pptrange_fac = factor(levels(trailingdry$pptrange_fac)), temprange_pos = 0, 
                        site = "NA", rep = "NA", family = "NA") 
predmdry$pred <- predict(mdry, newdata = predmdry, type = "conditional", se = F,  allow.new.levels = T) # fitted values for model

##### BIAS CORRECTED BOOTSTRAP CIS -----
dryboot<- drypred[sapply(drypred, is.numeric)] # remove character columns that are all NAs (i.e., bootstrap didnt converge) 

# add in actual predicted values
drypredict1 <- predmdry[, c(1,6)]
drypred <- cbind(drypredict1, dryboot)

head(drypred) # check first and last column that converged for next line of code
# proportion of bootstraps less than predicted and converted to z score
drypredz <- drypred %>% 
  rowwise() %>% 
  mutate(z = qnorm(sum(c_across(V9:V12005) < pred) / length(c_across(V9:V12005))))
drypredz$z

# calculate new probabilities at which to place bias corrected lower and upper CIs
drypredictprob <- drypredz %>%
  mutate(lower.prob = pnorm(2*z-1.96),# bias corrected lower quantile to calculate
         upper.prob = pnorm(2*z+1.96)) # bias corrected upper quantile to calculate
drypredictprob$lower.prob
drypredictprob$upper.prob

# bias corrected CIs
drypredictbias <- drypredictprob %>% 
  rowwise() %>% 
  transmute(pptrange_fac = pptrange_fac,
            pred = pred,
            lower.ci = quantile((c_across(V9:V12005)), probs = lower.prob),
            upper.ci = quantile((c_across(V9:V12005)), probs = upper.prob))

### save prediction with bias corrected 95% CI
save(drypredictbias, file = "outputs/models/dryedge_biascorrpredict.Rdata")





