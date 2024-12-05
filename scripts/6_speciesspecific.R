################################################################################
################################################################################
###                                                                          ###
###                            ANALYSIS SCRIPT 5:                            ###
###                 RECRUITMENT WITH RANGE POSITION MODELLING                ###
###                   SPECIES SPECIFIC MODELS AND PLOTS                      ###
###                                                                          ###
################################################################################
################################################################################


# Created Dec. 14, 2022
# Last Updated Dec. 4, 2024
# Authors: Katie Goodwin

# Analyses for Goodwin et al. 2024 Lagged climate-driven range shifts at species’ leading, but not trailing, range edges revealed by multispecies seed addition experiment. Ecography. DOI: 10.1111/ecog.07331
# modelling species specific recruitment patterns

# note that these models have low sample size and wide confidence intervals. Do not recommend using for species specific inferences due to low power. Rather these species-specific models were used to demonstrate the wide variablity in recruitment patterns across climatic ranges at the species level that underlies the community level patterns we observed in the other scripts.

###########################################################################################
#### PACKAGES AND DATA ------------------------------------------------
###########################################################################################
library(glmmTMB)
library(tidyverse)
library(ggplot2)
library(GGally)
library(ggeffects)
library(DHARMa)



# data
load('data/recruit_census.RData') 


###########################################################################################
##### RECRUITMENT COUNT - SPECIES IN AT LEAST 8 PLOTS -----------------------------------------------------
###########################################################################################
# data frame of species that recruited in at least 8 plots
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

###########################################################################################
##### INDIVIDUAL SPECIES MODELS -----------------------------------------------------
###########################################################################################
# start at species that recruited in most plots and go until model no longer runs
# note: models get more simplistic as we go for models to converge

### m.SORSIT -supports shift in optima -------
sorsit <- subset(census8, species == "SORSIT")

m.sorsit <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2) 
              + (1|site) , 
              ziformula = ~ 1,
              data = sorsit, 
              family = nbinom2(link = "log")) 
summary(m.sorsit) 

# diagnostics
msorsit_simres <- simulateResiduals(m.sorsit)
plot(msorsit_simres)
testUniformity(msorsit_simres) 
testOutliers(msorsit_simres) 
testQuantiles(msorsit_simres) 

### m.VACPAR - supports shift in optima -------
vacpar <- subset(census8, species == "VACPAR")

m.vacpar <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2) 
                    + (1|site) , 
                    ziformula = ~ 1,
                    data = vacpar, 
                    family = nbinom2(link = "log")) 
summary(m.vacpar) 

# diagnostics
mvacpar_simres <- simulateResiduals(m.vacpar)
plot(mvacpar_simres)
testUniformity(mvacpar_simres) 
testOutliers(mvacpar_simres) 
testQuantiles(mvacpar_simres) 


### m.RUBURS - neither refute or support -------
ruburs <- subset(census8, species == "RUBURS")

m.ruburs <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2) 
                    + (1|site) , 
                    ziformula = ~ 1,
                    data = ruburs, 
                    family = nbinom2(link = "log")) 
summary(m.ruburs) 

# diagnostics
mruburs_simres <- simulateResiduals(m.ruburs)
plot(mruburs_simres)
testUniformity(mruburs_simres) 
testOutliers(mruburs_simres) 
testQuantiles(mruburs_simres) 


### m.tolmen - neither support/refute based on range coverage (dont capture peak) -------
tolmen <- subset(census8, species == "TOLMEN")

m.tolmen <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2) 
                    + (1|site) , 
                    ziformula = ~ 1,
                    data = tolmen, 
                    family = nbinom2(link = "log")) 
summary(m.tolmen) 

# diagnostics
mtolmen_simres <- simulateResiduals(m.tolmen)
plot(mtolmen_simres)
testUniformity(mtolmen_simres) 
testOutliers(mtolmen_simres) 
testQuantiles(mtolmen_simres) 

### m.luplat - supports shift in optima -------
luplat <- subset(census8, species == "LUPLAT")

m.luplat <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2) 
                    + (1|site) , 
                    ziformula = ~ 1,
                    data = luplat, 
                    family = nbinom2(link = "log")) 
summary(m.luplat) 

# diagnostics
mluplat_simres <- simulateResiduals(m.luplat)
plot(mluplat_simres)
testUniformity(mluplat_simres) 
testOutliers(mluplat_simres) 
testQuantiles(mluplat_simres) 


### m.telgra - neither support/refute due to range converage (dont capture peak)-------
telgra <- subset(census, species == "TELGRA")

m.telgra <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2) 
                    + (1|site) , 
                    ziformula = ~ 1,
                    data = telgra, 
                    family = nbinom2(link = "log")) 
summary(m.telgra) 

# diagnostics
mtelgra_simres <- simulateResiduals(m.telgra)
plot(mtelgra_simres)
testUniformity(mtelgra_simres) 
testOutliers(mtelgra_simres) 
testQuantiles(mtelgra_simres) 


### m.vacdel NO SUPPORT -recruitment increases towards warm edge --------
vacdel <- subset(census, species == "VACDEL")

m.vacdel <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2) 
                    + (1|site) , 
                    ziformula = ~ 1,
                    data = vacdel, 
                    family = nbinom2(link = "log")) 
summary(m.vacdel) 

# diagnostics
mvacdel_simres <- simulateResiduals(m.vacdel)
plot(mvacdel_simres)
testUniformity(mvacdel_simres) 
testOutliers(mvacdel_simres) 
testQuantiles(mvacdel_simres) 


### m.eriper --------
eriper <- subset(census, species == "ERIPER")

# removed random effect of site, not great but doesnt run otherwise

m.eriper <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2),
                    zi = ~1,
                    data = eriper, 
                    family = nbinom2(link = "log")) 
summary(m.eriper) 


# diagnostics
meriper_simres <- simulateResiduals(m.eriper)
plot(meriper_simres)
testUniformity(meriper_simres) 
testOutliers(meriper_simres) 
testQuantiles(meriper_simres) 


### m.piceng --------
piceng <- subset(census, species == "PICENG")

# removed random effect of site and zero inflation, not great but doesnt run otherwise 

m.piceng <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2),
                    zi = ~1,
                data = piceng, 
                family = nbinom2(link = "log")) 
summary(m.piceng) 


# diagnostics
mpiceng_simres <- simulateResiduals(m.piceng)
plot(mpiceng_simres)
testUniformity(mpiceng_simres) 
testOutliers(mpiceng_simres) 
testQuantiles(mpiceng_simres) 

              
### m.aneocc --------
aneocc <- subset(census, species == "ANEOCC")


m.aneocc <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2),
                    zi = ~1,
                    data = aneocc, 
                    family = nbinom2(link = "log")) 
summary(m.aneocc) 


# diagnostics
mabilas_simres <- simulateResiduals(m.abilas)
plot(mabilas_simres)
testUniformity(mabilas_simres) 
testOutliers(mabilas_simres) 
testQuantiles(mabilas_simres) 


### m.abilas --------
abilas <- subset(census, species == "ABILAS")


m.abilas <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2),
                zi = ~1,
                data = abilas, 
                family = nbinom2(link = "log")) 
summary(m.abilas) 


# diagnostics
mabilas_simres <- simulateResiduals(m.abilas)
plot(mabilas_simres)
testUniformity(mabilas_simres) 
testOutliers(mabilas_simres) 
testQuantiles(mabilas_simres) 

### m.mahner --------
mahner <- subset(census, species == "MAHNER")

# removed random effect of site and zero inflation, not great but doesnt run otherwise 

m.mahner <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2),
                data = mahner, 
                family = nbinom2(link = "log")) 
summary(m.mahner) 


# diagnostics
mmahner_simres <- simulateResiduals(m.mahner)
plot(mmahner_simres)
testUniformity(mmahner_simres) 
testOutliers(mmahner_simres) 
testQuantiles(mmahner_simres) 


### m.erilan --------
erilan <- subset(census, species == "ERILAN")

# removed random effect of site and zero inflation, not great but doesnt run otherwise 

m.erilan <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2),
                data = erilan, 
                family = nbinom2(link = "log")) 
summary(m.erilan) 


# diagnostics
merilan_simres <- simulateResiduals(m.erilan)
plot(merilan_simres)
testUniformity(merilan_simres) 
testOutliers(merilan_simres) 
testQuantiles(merilan_simres) 


### m.mahaqu --------
mahaqu <- subset(census, species == "MAHAQU")

# removed random effect of site and zero inflation, not great but doesnt run otherwise 

m.mahaqu <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2),
                    zi = ~1,
                data = mahaqu, 
                family = nbinom2(link = "log")) 
summary(m.mahaqu) 


# diagnostics
mmahaqu_simres <- simulateResiduals(m.mahaqu)
plot(mmahaqu_simres)
testUniformity(mmahaqu_simres) 
testOutliers(mmahaqu_simres) 
testQuantiles(mmahaqu_simres) 


