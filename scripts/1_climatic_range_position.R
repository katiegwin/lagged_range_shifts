################################################################################
################################################################################
###                                                                          ###
###                         CLIMATIC RANGE FROM GBIF                         ###
###                             CL seed addition                             ###
###                                                                          ###
################################################################################
################################################################################

# Created Apr 14, 2022
# Last Updated Dec. 5 2024
# Authors: Katie Goodwin 

# Quantify regional climatic range for focal species from Cascades
# Legacy Seed Addition Experiment using GBIF occurrence data and Climate NA data
# GBIF data downloaded on Jan. 31, 2022. Download DOI: https://doi.org/10.15468/dd.wjqxtc 
# Climate NA data downloaded on Oct. 17, 2022 Wang et al 2016

# Quantifies climatic range for each species for a defined region and specifies the relative
# climatic range position for each focal species at each seed addition site


#### INPUTS ----
load('data/occur_biome_clim_Normal_1981_2010Y.RData') # tidied gbif occurrence records within the biome with climate NA data (full untidied version is way tooo large for github so starting here) 
head(occbio)

siteclim <- read.csv('data/CLseed_sitesclim_Normal_1981_2010Y.csv') # climate data for seed addition sites
head(siteclim)

species <- read.csv('data/species_desc.csv') # species description data

load("data/recruit_census.RData")

library(ggplot2)
library(tidyverse)

##############################################################################################
###### MAT and MAP supp plot #############################################
##############################################################################################
pdf(file = "outputs/final/correlations.pdf",   # The directory you want to save the file in
   width = 4, # The width of the plot in inches
    height = 3)
ggplot() + 
  geom_point(data = occbio, aes(x = MAT, y = MAP), alpha = 0.3) +
  geom_point(data = siteclim, aes(x = MAT, y = MAP), colour = "#00CC66", alpha = 0.8) +
  annotate("text", x = -1, y=6800, label = "r = 0.17", col = "black")+ 
  annotate("text", x = -1, y=6300, label = "r = 0.23", col = "#00CC66")+ 
  theme_classic()
dev.off()

##############################################################################################
##### CALCULATE CLIMATIC RANGE FOR FOCAL SPECIES #############################################
##############################################################################################

# TEMPERATURE (MAT)
tempbio <- occbio %>%
  group_by(species) %>% #for each species
  summarize(n = n(), #sample size
            lower_temp = quantile(MAT, probs = 0.025), #2.5% quantile
            upper_temp = quantile(MAT, probs = 0.975)) #97.5% quantile

# PPT (MAP)
pptbio <- occbio %>%
  group_by(species) %>% #for each species
  summarize(n = n(), #sample size
            lower_ppt = quantile(MAP, probs = 0.025), #2.5% quantile
            upper_ppt = quantile(MAP, probs = 0.975)) #97.5% quantile

rangebio <- full_join(tempbio, pptbio, by = c("species", "n"))



##############################################################################################
##### CALCULATE CLIMATIC RANGE CENTRE AND EXTENT FOR EACH SPECIES -  ################
##############################################################################################

### calculate climatic range centre and extent (i.e., upper limit - lower limit) for each focal species 
# TEMPERATURE
tempbio <- tempbio %>%
  mutate(centre_temp = ((upper_temp + lower_temp) / 2)) %>%
  mutate(extent_temp = (upper_temp - lower_temp))

# PRECIPITATION
pptbio <- pptbio %>%
  mutate(centre_ppt = ((upper_ppt + lower_ppt) / 2)) %>%
  mutate(extent_ppt = (upper_ppt - lower_ppt))

### combine to one data frame of climatic range
rangebio <- full_join(tempbio, pptbio, by = c("species", "n"))

#rename to be same as census data
rangebio <- rangebio %>% 
  rename(full_species = species) 

# add in 6 letter species codes
rangebio1 <- full_join(rangebio, species)

# reorder columns
rangebio <- rangebio1[, c(11, 3, 4, 5, 6, 7, 8, 9, 10)]

# add species climatic range data to census data
census2 <- full_join(census, rangebio, by = "species") 
mism <- anti_join(census, rangebio, by = "species") # no mismatches
mism1 <- anti_join(rangebio, census, by = "species") # no mismatches

summary(census2)
NAs <- census2[is.na(census2$lower_temp),] #no NAs for newly added climate data :)

censusbio <- census2


##############################################################################################
##### STANDARDIZE CLIMATIC RANGE POSITION FOR EACH SITE FOR EACH SPECIES - BIOME #####################
##############################################################################################
# where within a species climatic range the site is, standardized to compare between species. 
# Formula is:
# (site temp - midpoint temp of species range) / (extent of range i.e., upper - lower limit) * 0.5)
# determines how many standardized units below or above midpoint of range it is
# upper temperature or ppt limit will have a value of +1
# range midpoint will have a value of 0
# lower temperature or ppt limit will have a value of -1
# values >1 mean the site is warmer or wetter than the climatic range of the species
# values <1 mean the site is colder or drier than the climatic range of the species
censusbio <- mutate(censusbio, temprange_pos =  ((MAT - centre_temp) / (extent_temp * .5)))
censusbio <- mutate(censusbio, pptrange_pos = ((MAP - centre_ppt) / (extent_ppt * .5)))


###########################################################################################
##### DATA FRAME JUST OF RANGE POSITION FOR EACH SPECIES AT EACH SITE ---------
######################################################################################

rangepos <- (data = censusbio) %>%
  distinct(region, site, lat, long, elev, canopy, MAT, MAP, species, temprange_pos, pptrange_pos)

# save csv
write.csv(rangepos, "data/siterangepositions.csv")

