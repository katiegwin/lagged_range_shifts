################################################################################
################################################################################
###                                                                          ###
###                 RECRUITMENT WITH RANGE POSITION MODELLING                ###
###               CODE FOR CLUSTER BOOTSTRAP CANOPY - CONDITIONAL
###                                                                          ###
################################################################################
################################################################################


# Created Dec. 12, 2022
# Last Updated Dec. 12, 2022
# Authors: Katie Goodwin


# code for bootstrap run copied from canopy_analysis.R for cluster - conditional part


###########################################################################################
#### PACKAGES AND DATA ------------------------------------------------
###########################################################################################
#install.packages("glmmTMB", lib="~/R_Libs/")

#install.packages("glmmTMB", dependencies = T, repos="http://cran.r-project.org", lib="~/R_Libs/")
#install.packages("tidyverse", repos="http://cran.r-project.org", lib="~/R_Libs/")

getRversion() # double check right version of r for glmmTMB (should be 4.2-1)

library(glmmTMB, lib="~/R/x86_64-pc-linux-gnu-library/4.2/") # cluster run only
#library(glmmTMB, lib="~/R_Libs/")
library(tidyverse)

load('canopy/census5_for_cluster.RData') # cluster run only

###########################################################################################
##### MODEL  -----------------------------------------------------
###########################################################################################

mcan <- glmmTMB(new_all ~  poly(temprange_pos, 2)  + canopycont + temprange_pos:canopycont + region
                + (1|site/rep) + (1|family), 
                ziformula = ~ poly(temprange_pos, 2)  + canopycont 
                + (1|site/rep) + (1|family),
                data = census5, 
                family = nbinom2(link = "log")) 
summary(mcan)


##### BOOTSTRAPPING - DO OPEN AND CLOSED CANOPY OPTIMA DIFFER? ---------------
# i.e., do 95% confidence intervals of peak of curve include zero?

nreps <- 20000 # number of reps - note that lots of these reps will not converge

# empty data frames: 
optima <- data.frame(opt_12 = rep(NA, nreps), opt_36 = rep(NA, nreps), opt_50 = rep(NA, nreps),
                     opt_75 = rep(NA, nreps)) # for optima
pred <- expand.grid(temprange_pos = seq(-2.6,1.6, by = .01), canopycont = c(12, 36, 50, 75), 
                    site = "NA", rep = "NA", family = "NA", region = "RP") # for predicted values. will take pop means for random effects. Since no interaction terms value of region shouldnt matter
warnings <- data.frame(warning = rep(NA, nreps)) # for bootstrap warnings


# the looooop - optimal thermal range position for conditional part of model
for (i in 1:nreps) {
  save(pred, file = 'canopy/canopy_cond_bootstraps06OCT23.RData') # intermittent save if case crashes
  save(optima, file = 'canopy/canopy_cond_optimaCI06OCT23.RData') #intermittent save if case crashes
  save(warnings, file = 'canopy/canopy_cond_bootstrapwarnings06OCT23.RData') # intermittent save if crashes
  print(paste0("Running model ", i)) # keep track of progress
  newmodel <- quietly(safely((refit)))(mcan, simulate(mcan)[[1]]) # bootstrap run that saves warnings (quietly) AND errors (safely) - warnings are convergence issues and errors are fatal (e.g., entire simulated dataset are zeroes)
  if (is.null(newmodel$result$error) == TRUE) { # if there were no fatal errors...
    if(length(newmodel$warnings) == 0) { # and there were no warnings...
      print(paste0("model ", i, " converged!! yay!"))
      pred[ , i + 6] <- predict(newmodel$result$result, newdata = pred, type = "conditional", se = F,  allow.new.levels = T) # save fitted values from bootstrap run as new column
      can12 <- subset(pred, canopycont == "12") # closed canopy mean
      optima[i, 1 ] <- can12[which.max(can12[[i + 6]]), c(1)]
      can36 <- subset(pred, canopycont == "36") # open canopy mean
      optima[i, 2 ] <- can36[which.max(can36[[i + 6]]), c(1)]
      can50 <- subset(pred, canopycont == "50") # open 50% open canopy
      optima[i, 3 ] <- can50[which.max(can50[[i + 6]]), c(1)]
      can75 <- subset(pred, canopycont == "75") # open 75% open canopy
      optima[i, 4 ] <- can75[which.max(can75[[i + 6]]), c(1)]
    } else{ # if there were warnings everything is NA
      print(paste0("model ", i," had warnings"))
      warnings[i,] <- "warning" #records no outputs b/c warning 
      pred[ , i+6] <- "NA"
      optima[i, ] <- "NA"
    } 
  } else { # if there was a fatal error everything is NA
    print(paste0("model ", i, " had a fatal error"))
    warnings[i,] <- "error" # records no outputs b/c error
    pred[ , i+6] <- "NA"
    optima[i, ] <- "NA"
  }
} 


