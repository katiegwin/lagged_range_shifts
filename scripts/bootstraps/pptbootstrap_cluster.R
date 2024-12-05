################################################################################
################################################################################
###                                                                          ###
###                 RECRUITMENT WITH RANGE POSITION MODELLING                ###
###                       CODE FOR CLUSTER BOOTSTRAP  PPT                       ###
###                                                                          ###
################################################################################
################################################################################


# Created Nov 16, 2022
# Last Updated Nov. 17, 2022
# Authors: Katie Goodwin


# code for bootstrap run copied from range_position_analysis.R for cluster


###########################################################################################
#### PACKAGES AND DATA ------------------------------------------------
###########################################################################################
#install.packages("glmmTMB", lib="~/R_Libs/")

#install.packages("glmmTMB", dependencies = T, repos="http://cran.r-project.org", lib="~/R_Libs/")
#install.packages("tidyverse", repos="http://cran.r-project.org", lib="~/R_Libs/")

getRversion()

library(glmmTMB, lib="~/R/x86_64-pc-linux-gnu-library/4.2/") # cluster run only
#library(glmmTMB, lib="~/R_Libs/")
library(tidyverse)

load('census8_forcluster.RData') # cluster run only

###########################################################################################
##### MODEL  -----------------------------------------------------
###########################################################################################

m1 <- glmmTMB(new_all ~  poly(temprange_pos, 2) + poly(pptrange_pos, 2) + canopy 
              + (1|site/rep) + (1|family), 
              ziformula = ~ poly(temprange_pos, 2) + poly(pptrange_pos, 2) + canopy 
              + (1|site/rep) + (1|family),
              data = census8, 
              family = nbinom2(link = "log"),
              na.action=na.exclude) 

summary(m1) 

##### BOOTSTRAPPING - HAS OPTIMAL THERMAL RECRUITMENT POSITION SHIFTED FROM ZERO? ---------------
# i.e., do 95% confidence intervals of peak of curve include zero?

nreps <- 20001 # number of reps - note that lots of these reps will not converge

# empty data frames: 
optima <- data.frame(opt_pptrange = rep(NA, nreps), opt_recruits = rep(NA, nreps)) # for optima
pred <- expand.grid(pptrange_pos = seq(-1.6,2.3, by = .01), temprange_pos = 0, canopy = "open", 
                    site = "NA", rep = "NA", family = "NA") # for predicted values. will take pop means for random effects. Since no interaction terms value of canopy or ppt shouldnt matter
warnings <- data.frame(warning = rep(NA, nreps)) # for bootstrap warnings


# the looooop
for (i in 1:nreps) {
  save(pred, file = 'pptrange_cond_bootstraps01DEC22.RData') # intermittent save if case crashes
  save(optima, file = 'pptrange_cond_optimaCI01DEC22.RData') #intermittent save if case crashes
  save(warnings, file = 'pptrange_cond_bootstrapwarnings01DEC22.RData') # intermittent save if crashes
  print(paste0("Running model ", i)) # keep track of progress
  newmodel <- quietly(safely((refit)))(m1, simulate(m1)[[1]]) # bootstrap run that saves warnings (quietly) AND errors (safely) - warnings are convergence issues and errors are fatal (e.g., entire simulated dataset are zeroes)
  if (is.null(newmodel$result$error) == TRUE) { # if there were no fatal errors...
    if(length(newmodel$warnings) == 0) { # and there were no warnings...
      print(paste0("model ", i, " converged!! yay!"))
      pred[ , i + 6] <- predict(newmodel$result$result, newdata = pred, type = "conditional", se = F,  allow.new.levels = T) # save fitted values from bootstrap run as new column
      optima[i, ] <- pred[which.max(pred[[i + 6]]), c(1, i + 6)] # save max fitted value from bootstrap run (i.e., peak of curve - optimal recruitment) as new row
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
