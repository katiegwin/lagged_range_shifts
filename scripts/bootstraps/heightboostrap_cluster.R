################################################################################
################################################################################
###                                                                          ###
###                 RECRUITMENT WITH RANGE POSITION MODELLING                ###
###               CODE FOR CLUSTER BOOTSTRAP GROWTH MODEL
###                                                                          ###
################################################################################
################################################################################


# Created Dec 5, 2022
# Last Updated Dec 5, 2022
# Authors: Katie Goodwin


# code for bootstrap run testing growth model


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

load('growthsurv/year1sum_for_cluster.RData') # cluster run only

###########################################################################################
##### MODEL  -----------------------------------------------------
###########################################################################################

#### model - including region does not change results
mheight <- glmmTMB(stht ~ temprange_pos  + pptrange_pos + 
                     (1|site/rep),  
                   data = year1sum, family = gaussian)

summary(mheight)

##### BOOTSTRAPPING - Confidence limits for cold edge hypothesis ---------------
# i.e., do 95% confidence intervals of peak of curve include zero?

nreps <- 20000 # number of reps - note that lots of these reps will not converge

# empty data frames: 
pred <- expand.grid(temprange_pos = seq(-2.4,1.3, by = .01), pptrange_pos = 0, canopy = "open", 
                    site = "NA", rep = "NA", family = "NA") # for predicted values. will take pop means for random effects. Since no interaction terms value of canopy or ppt shouldnt matter
warnings <- data.frame(warning = rep(NA, nreps)) # for bootstrap warnings


# the looooop
for (i in 1:nreps) {
  save(pred, file = 'growthsurv/height_pred08DEC22.RData') #intermittent save if case crashes
  save(warnings, file = 'growthsurv/height_warnings08DEC22.RData') # intermittent save if crashes
  print(paste0("Running model ", i)) # keep track of progress
  newmodel <- quietly(safely((refit)))(mheight, simulate(mheight)[[1]]) # bootstrap run that saves warnings (quietly) AND errors (safely) - warnings are convergence issues and errors are fatal (e.g., entire simulated dataset are zeroes)
  if (is.null(newmodel$result$error) == TRUE) { # if there were no fatal errors...
    if(length(newmodel$warnings) == 0) { # and there were no warnings...
      print(paste0("model ", i, " converged!! yay!"))
      pred[ , i + 6] <- predict(newmodel$result$result, newdata = pred, type = "response", se = F,  allow.new.levels = T) # save fitted values from bootstrap run as new column
    } else{ # if there were warnings everything is NA
      print(paste0("model ", i," had warnings"))
      warnings[i,] <- "warning" #records no outputs b/c warning 
      pred[ , i+6] <- "NA"
    } 
  } else { # if there was a fatal error everything is NA
    print(paste0("model ", i, " had a fatal error"))
    warnings[i,] <- "error" # records no outputs b/c error
    pred[ , i+6] <- "NA"
  }
} 
