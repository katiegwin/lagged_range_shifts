################################################################################
################################################################################
###                                                                          ###
###                 RECRUITMENT WITH RANGE POSITION ANALYSIS                 ###
###                         PLOTS FOR PUBLICATION                            ###
###                                                                          ###
################################################################################
################################################################################


# Created Apr 5, 2023
# Last Updated Jan 24, 2024 
# Authors: Katie Goodwin

# Final plots for publication 

# analysis done in other scripts. Model outputs saved to use here. 

###########################################################################################
#### PACKAGES AND DATA ------------------------------------------------
###########################################################################################
library(ggplot2)
library(ggthemes)
library(tidyverse)
library("cowplot")

# range position of study sites
rangepos <- read.csv(file = "data/siterangepositions.csv")

# census data
load('data/recruit_census.RData') 

### outputs for model testing whether optimal recruitment has shifted from thermal centre. See recruitment_optima_analysis.R for how all these came to be
load('outputs/bootstraps/temp-cond/temprange_cond_optimaCI08DEC22.RData') # thermal optima from bootstraps
topt <- optima
load("outputs/models/optimarecruit_biascorrpredict.Rdata") # predicted trend line using predict function for main line and  bias corrected 95% confidence intervals of bootstraps

# fix classes
topt$opt_temprange <- as.numeric(as.character(topt$opt_temprange))
topt$opt_recruits <- as.numeric(as.character(topt$opt_recruits))
str(topt)

### outputs for model testing for successful recruitment beyond cool edge (see recruitment_edge_analysis.R for how all these came to be)
load("outputs/models/coolhypodata.RData") # data frame of species that recruited with sites beyond range with new within/beyond range factor of sites with thermal range position <0 (cool)
load("outputs/models/coldedge_biascorrpredict.Rdata") # predicted trend line using predict function for main line and bias corrected95% confidence intervals of bootstraps

### outputs for model testing for recruitment failure beyond new warm edge (see recruitment_edge_analysis.R for how all these came to be)
load("outputs/models/warmhypodata.Rdata") # data frame of species that recruited with sites beyond new warm limit with within/beyond factor and sites with thermal range position >0 (warm)
load("outputs/models/hotedge_biascorrpredict.Rdata") # Bias corrected CIs

### canopy cover plot (see canopy_analysis.R)
load("outputs/models/canopy_optima.Rdata") # thermal optima from bootstraps for mean value of open and closed canopy treatments

## zi part of optimal recruitment model
load("outputs/models/optimalzi_biascorrpredict.Rdata")

## canopy bias corrected bootstrap
load("outputs/models/canopy_biascorrpredict.Rdata")

## ppt prediction plot
load("outputs/models/ppt_biascorrpredict.Rdata")
load('outputs/bootstraps/ppt/pptrange_cond_bootstraps01DEC22.RData')
load('outputs/bootstraps/ppt/pptrange_cond_optimaCI01DEC22.RData')
popt <- optima

popt$opt_pptrange <- as.numeric(as.character(popt$opt_pptrange))
popt$opt_recruits <- as.numeric(as.character(popt$opt_recruits))
str(popt)


colcold <- c("#90caf9", "#1565c0")
colhot <- c("#e53935", "#ef9a9a")


###########################################################################################
#### RANGE POSITION OF SITES ------------------------------------------------
##########################################################################################
#### plot of each sites location within each thermal species range 
pdf(file = "outputs/final/thermal_ranges.pdf",   # The directory you want to save the file in
    width = 17, # The width of the plot in inches
    height = 10)
ggplot(data = rangepos) +
  geom_hline(yintercept=1, col = "red", size =2) +
  geom_hline(yintercept =0,  size = 2) +
  geom_hline(yintercept=-1, col="#0288D1", size = 2) +
  geom_point(mapping = aes (x = species, y = temprange_pos, color="darkolivegreen4"), size = 8, alpha=0.3) +
  scale_color_identity(guide=guide_none()) +
  scale_y_reverse() +
  ylab("Thermal range position") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size = 28)) +
  theme(axis.text = element_text(size = 22))  +
  theme(axis.title = element_text(size = 36))   +
  theme(plot.margin = margin(1,1,1,1, "cm"))
dev.off()

#### plot of each sites location within each precipitation species range 
pdf(file = "outputs/final/pptranges.pdf",   # The directory you want to save the file in
    width = 17, # The width of the plot in inches
    height = 10)
ggplot(data = rangepos) +
  geom_hline(yintercept=1, col = "#0288D1", size =2) +
  geom_hline(yintercept =0,  size = 2) +
  geom_hline(yintercept=-1, col="red", size = 2) +
  geom_point(mapping = aes (x = species, y = pptrange_pos, color="darkolivegreen4"), size = 8, alpha=0.3) +
  scale_color_identity(guide=guide_none()) +
  ylab("Precipitation range position") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size = 28)) +
  theme(axis.text = element_text(size = 22))  +
  theme(axis.title = element_text(size = 30))   +
  theme(plot.margin = margin(1,1,1,1, "cm"))
dev.off()

###########################################################################################
#### THERMAL OPTIMA PLOTS ------------------------------------------------
###########################################################################################

##### data frame for optima shift (species that recruited in at least 8 plots)
census8 <- census %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ABILAS" | species == "ANEOCC" |
           species == "ERIPER" | species == "PICENG" |
           species == "RUBURS" | species == "SORSIT" | species == "TELGRA" |
           species == "TOLMEN" | species == "VACDEL" | species == "VACPAR" |
           species == "MAHNER" | species == "ERILAN" | species == "MAHAQU") # remove species with indiv in less than 8 plots



# version for publication
pdf(file = "outputs/final/thermalrecruit.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 3)
ggplot(tpredictbias) +
  geom_vline(xintercept=1, col = "red", size =0.5) +
  geom_vline(xintercept = -0.26, linetype = 2, size = 0.5) +
  geom_vline(xintercept=-1, col="#0288D1", size = 0.5) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_point(data = census8, aes(x = temprange_pos, y = log(new_all + 1) ), col ="black", alpha = 0.4) +
  geom_line(aes(x = temprange_pos, y = log(pred + 1)), col = "black") +
  geom_ribbon(aes(x = temprange_pos, ymin = log(lower.ci + 1), ymax = log(upper.ci +1)), fill = "black", alpha=0.3) +
  scale_x_reverse() +
  ylab("ln (Seedling Count + 1)") + xlab("Thermal Range Position") + 
  theme_classic()
dev.off()


###########################################################################################
#### BOOTSTRAP OPTIMA MULTI PANEL PLOT ------------------------------------------------
###########################################################################################

##### optima recruitment from bootstrap runs 
# plot for botany symposium (transparent to be presented on dark background)
breaks <- seq(from = -2.6, to = 1.6, by = 0.1) # ensure there is a break at zero

# version for publication - temperature
# subset to 5000 rows 
topt1 <- na.omit(topt)
topt2 <- topt1[-c(5001:5057),] 

histt <- ggplot(topt2) +
  geom_histogram(aes(x = opt_temprange), binwidth = .1, breaks = breaks, fill = "gray", col = "black") + scale_x_reverse() +
  geom_vline(xintercept=1, col = "red", size = 0.5) +
  geom_vline(xintercept = -0.26, linetype = 2, size = 0.5) +
  geom_vline(xintercept=-1, col="#0288D1", size = 0.5) +
  geom_vline(xintercept = -0, size = 0.5) +
  xlab("Thermal Range Position") +
  ylab("Bootstrap Optimum Frequency")  + 
  theme_classic()

# version for publication - precipitation
# subset to 5000 bootstraps
popt1 <- na.omit(popt)
popt2 <- popt1[-c(5001:9948),] 

histp <- ggplot(popt2) +
  geom_histogram(aes(x = opt_pptrange), bins = 40, fill = "gray", col = "black") +
  geom_vline(aes(xintercept = 0), linetype = 2) + xlab("Precipitation Range Position") +
  geom_vline(xintercept= -1, col = "red", size = 0.5) +
  geom_vline(xintercept = 0.12, linetype = 2, size = 0.5) +
  geom_vline(xintercept= 1, col="#0288D1", size = 0.5) +
  geom_vline(xintercept = -0, size = 0.5) +
  theme_classic() +
  theme(axis.title.y=element_blank()) 

#version for publication - canopy cover
# subset to 5000 bootstraps (for each treatment so 10000)
toptcanopy1 <- na.omit(toptcanopy)
toptcanopy2 <- toptcanopy1[-c(10001:29288),] 
col <- c("#33691E", "#FFA000")

histc <- ggplot(toptcanopy2, aes(x = opt_temprange, fill = canopy)) +  # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.6, bins = 50, col = "black") +
  scale_fill_manual(values = col)  +
  theme_classic() + 
  ylab("Bootstrap Optima Frequency") +
  scale_x_reverse() +
  geom_vline(xintercept=1, col = "red", size = 0.5) +
  geom_vline(xintercept=-1, col="#0288D1", size = 0.5) +
  geom_vline(xintercept = -0, size = 0.5) +
  xlab("Thermal Range Position") +
  theme(legend.position = "none",
        axis.title.y=element_blank()) 
histc

##### panel plot with all the different bootstrap optima histogrames
pdf(file = "outputs/final/bootstraphists.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 3)
ggdraw() +
  draw_plot(histt, x = 0, y = 0, width = .33, height = 1) +
  draw_plot(histp, x = .33, y = 0, width = .33, height = 1) +
  draw_plot(histc, x = .66, y = 0, width = .33, height = 1) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0.005, 0.305, 0.635), y = c(1, 1, 1))
dev.off()


###########################################################################################
#### RECRUITMENT BEYOND EDGES PLOT ------------------------------------------------
###########################################################################################
# plot for manuscript
cpredictbias$temprange_fac1 <- factor(cpredictbias$temprange_fac,     # Reorder factor levels
                         c("within", "beyond"))
leading$temprange_fac1 <- factor(leading$temprange_fac,     # Reorder factor levels
                                      c("within", "beyond"))
cold <- ggplot(cpredictbias) +
  geom_violin(data = leading, aes(x = temprange_fac1, y = log(new_all+1), fill = temprange_fac1), alpha = 0.6, col = "gray25") +
  scale_fill_manual(values = colcold)  +
  geom_point(aes(x = temprange_fac1, y = log(pred +1)), size = 1, col = "black") +
  geom_linerange(aes(x = temprange_fac, ymin = log(lower.ci + 1), ymax = log(upper.ci + 1)), linewidth = .5, col = "black") +
  ylab("ln(Seedling Count + 1)") + xlab("Thermal range position") + 
  scale_x_discrete(labels = c('Cool Half Range', 'Beyond Cold Edge')) +
  #annotate("text", x = 2.3, y=5, label = "P + __", col = "black")+ 
  theme_classic() +
  theme(legend.position = "none")
cold

# plot for manuscript
hot <- ggplot(hotpredictbias) +
  geom_violin(data = trailing, aes(x = temprange_fac, y = log(new_all+1), fill = temprange_fac), alpha=1, col = "gray25") +
  scale_fill_manual(values = colhot)  +
  geom_point(aes(x = temprange_fac, y = log(pred +1)), size = 1, col = "black") +
  geom_linerange(aes(x = temprange_fac, ymin = log(lower.ci + 1), ymax = log(upper.ci + 1)), size = .5, col = "black") +
  ylab("ln(Seedling Count + 1)") + xlab("Thermal range position") + 
  scale_x_discrete(labels = c('Presumed Unsuitably Warm', 'Warm Half Range')) +
 # annotate("text", x = 2.3, y=4.3, label = "P = 0= __", col = "black")+ 
  theme_classic() +
  theme(legend.position = "none")
hot


##### panel plot with both edges
pdf(file = "outputs/final/edgespred.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
ggdraw() +
  draw_plot(hot, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(cold, x = .5, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1))
dev.off()



###########################################################################################
#### CANOPY PLOT ------------------------------------------------
###########################################################################################
censussuc <- census %>% 
  group_by(species) %>% 
  filter(species != "CARSPE" & species != "MAIDIL" & species !="MAIRAC" & species != "PINCON") 

census5 <- censussuc %>% #already subseted to species that recruited at least once
  group_by(species) %>% 
  filter(species != "CARSTI" & species != "SAMCER" & species !="RUBSPE" & species != "PINPON" & species != "SAMRAC" &
           species != "CARSPE" & species != "MAIDIL" & species !="MAIRAC" & species != "PINCON") 

census5a <- subset(census5, new_all != "NA")
census5 <- census5a

# where are all the optimas?
canpeaks <- canpredictbias %>%
  group_by(canopycont) %>%
  summarize(max = temprange_pos[which.max(pred)])
canpeaks

col <- c("#33691E", "#FFA000", "#CE93d8", "#4fc3f7")
toptcanopy$canopy <- plyr::revalue(toptcanopy$canopy, c("mean closed" = "closed",
                            "mean open" = "open"))
canpredictbias$canopycont <- as.factor(canpredictbias$canopycont)


pdf(file = "outputs/final/canopypred.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 3)
ggplot(data = canpredictbias) +
  geom_vline(xintercept=1, col = "red", size =0.5) +
  geom_vline(xintercept =0, size = 0.5) +
  geom_vline(xintercept=-1, col="#0288D1", size = 0.5) +
  geom_vline(xintercept=-0.23, linetype = 2, col="#33691E", size = 0.5) +
  geom_vline(xintercept=-0.38, linetype = 2, col="#FFA000", size = 0.5) +
  geom_vline(xintercept=-0.47, linetype = 2, col="#CE93d8", size = 0.5) +
  geom_vline(xintercept=-0.62, linetype = 2, col="#4fc3f7", size = 0.5) +
  geom_point(data = census5, aes(x = temprange_pos, y = log(new_all + 1),), alpha = 0.5) +
  geom_ribbon(aes(x = temprange_pos, ymin = log(lower.ci + 1), ymax = log(upper.ci + 1), fill = canopycont),  alpha = 0.2) +
  geom_line(mapping = aes (x = temprange_pos, y = log(pred + 1), col = canopycont)) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  scale_x_reverse() +
  theme_classic() +
  ylab("ln (Seedling Count + 1)") + xlab("Thermal Range Position") +
 theme(legend.position = "none")
dev.off()




###########################################################################################
#### ZERO INFLATION PART OF OPTIMAL MODEL ------------------------------------------------
###########################################################################################

# plot for manuscript
pdf(file = "outputs/final/zi_pred.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 4)
ggplot(tpredictzibias) +  
  geom_line(aes(x = temprange_pos, y = (1-pred)), col = "black") +
  geom_ribbon(aes(x = temprange_pos, ymin = (1-lower.ci), ymax = (1-upper.ci)), fill = "black", alpha=0.3) +
  ylab("Probability of Recruitment") + xlab("Thermal Range Position") +
  scale_x_reverse() +
  theme_classic() 
dev.off()


###########################################################################################
####  PPTOPTIMA SHIFT PLOTS ------------------------------------------------
###########################################################################################

##### data frame for optima shift (species that recruited in at least 8 plots)
census8 <- census %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ABILAS" | species == "ANEOCC" |
           species == "ERIPER" | species == "PICENG" |
           species == "RUBURS" | species == "SORSIT" | species == "TELGRA" |
           species == "TOLMEN" | species == "VACDEL" | species == "VACPAR" |
           species == "MAHNER" | species == "ERILAN" | species == "MAHAQU") # remove species with indiv in less than 8 plots

# version for publication
pdf(file = "outputs/final/pptpred.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
   height = 3)
ggplot(ppredictbias) +
  geom_vline(xintercept=-1, col = "red", size =0.5) +
  geom_vline(xintercept =0,  size = 0.5) +
  geom_vline(xintercept =0.12, linetype = 2, size = 0.5) +
  geom_vline(xintercept=1, col="#0288D1", size = 0.5) +
  geom_point(data = census8, aes(x = pptrange_pos, y = log(new_all + 1) ), col ="black", alpha = 0.4) +
  geom_line(aes(x = pptrange_pos, y = log(pred + 1)), col = "black") +
  geom_ribbon(aes(x = pptrange_pos, ymin = log(lower.ci + 1), ymax = log(upper.ci +1)), fill = "black", alpha=0.3) +
  ylab("ln (Seedling Count + 1)") + xlab("Precipitaiton Range Position") +
  theme_classic()
dev.off()


################################################################################
#### PRECIPITATION EDGE PLOTS ----
################################################################################
# data and model output
load("outputs/models/wetedge_biascorrpredict.Rdata") # predictions for plot
load("outputs/models/wethypodata.RData") # data for wet model

load("outputs/models/dryedge_biascorrpredict.Rdata") # predictions for plot
load("outputs/models/dryhypodata.RData") # data for dry model

# wet plot
leadingwet$pptrange_fac1 <- factor(leadingwet$pptrange_fac,     # Reorder factor levels
                                   c("within", "beyond"))
wetpredictbias$pptrange_fac1 <- factor(wetpredictbias$pptrange_fac,     # Reorder factor levels
                                   c("within", "beyond")) 


wetplot <- ggplot(wetpredictbias) +
  geom_violin(data = leadingwet,                      # adding the raw data 
              aes(x = pptrange_fac1, y = log(new_all + 1), fill = pptrange_fac1), alpha = 0.6,  scale = "width", col = 'gray25') +   
  scale_fill_manual(values = colcold)  +
  geom_point(aes(x = pptrange_fac1, y = log(pred + 1))) +
  geom_linerange(aes(x = pptrange_fac1, ymin = log(lower.ci + 1), ymax = log(upper.ci + 1)), size = 0.5, col = "black") + 
  ylab("ln(Seedling Count + 1") + xlab("Precipitation Range Position")+ theme_classic() +
  scale_x_discrete(labels = c('Wet Half Range', 'Beyond Wet Edge')) +
  theme(legend.position = "none")
wetplot

# dry plot

dryplot <- ggplot(drypredictbias) +
  geom_violin(data = trailingdry,                      # adding the raw data 
              aes(x = pptrange_fac, y = log(new_all + 1), fill = pptrange_fac), alpha = 1, col = "gray25", scale = "width") +
  scale_fill_manual(values = colhot)  +
  geom_point(aes(x = pptrange_fac, y = log(pred + 1))) +
  geom_linerange(aes(x = pptrange_fac, ymin = log(lower.ci + 1), ymax = log(upper.ci + 1)), size = 0.5, col = "black") + 
  ylab("ln(Seedling Count + 1") + xlab("Precipitation Range Position")+ theme_classic() +
  scale_x_discrete(labels = c('Presumed Unsuitably Dry', 'Dry Half Range')) +
  theme(legend.position = "none")
dryplot


# combining the two ppt plots
pdf(file = "outputs/final/pptedgespred.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
ggdraw() +
  draw_plot(dryplot, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(wetplot, x = .5, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1))
dev.off()
