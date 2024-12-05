################################################################################
################################################################################
###                                                                          ###
###                 RECRUITMENT WITH RANGE POSITION MODELLING                ###
###               CODE FOR CLUSTER BOOTSTRAP COLD EDGE
###                                                                          ###
################################################################################
################################################################################


# Created Dec 5, 2022
# Last Updated Jan. 17, 2024
# Authors: Katie Goodwin


# code for bootstrap run testing cold edge prediction


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

load('edge/leading_forcluster.RData') # cluster run only

###########################################################################################
##### MODEL  -----------------------------------------------------
###########################################################################################

mcold <- glmmTMB(new_all ~  temprange_fac + poly(pptrange_pos, 2) 
                 + (1|site/rep) + (1|family), 
                 ziformula = ~ 1,
                 data = leading, 
                 family = nbinom2(link = "log")) 

summary(mcold)

##### BOOTSTRAPPING - Confidence limits for cold edge hypothesis ---------------
# i.e., do 95% confidence intervals of peak of curve include zero?

nreps <- 20000 # number of reps - note that lots of these reps will not converge

# empty data frames: 
pred <- expand.grid(temprange_fac = levels(leading$temprange_fac), pptrange_pos = 0, canopy = "open", 
                    site = "NA", rep = "NA", family = "NA") # for predicted values. will take pop means for random effects. Since no interaction terms value of canopy or ppt shouldnt matter
warnings <- data.frame(warning = rep(NA, nreps)) # for bootstrap warnings


# the looooop
for (i in 1:nreps) {
  save(pred, file = 'edge/coldedge_pred17JAN24.RData') #intermittent save if case crashes
  save(warnings, file = 'edge/coldedge_warnings17JAN24.RData') # intermittent save if crashes
  print(paste0("Running model ", i)) # keep track of progress
  newmodel <- quietly(safely((refit)))(mcold, simulate(mcold)[[1]]) # bootstrap run that saves warnings (quietly) AND errors (safely) - warnings are convergence issues and errors are fatal (e.g., entire simulated dataset are zeroes)
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
