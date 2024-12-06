# lagged_range_shifts
Code and data for Goodwin et al. 2024. Lagged climate-driven range shifts at species’ leading, but not trailing, range edges revealed by multispecies seed addition experiment in Ecography DOI: 10.1111/ecog.07331

# Compiled data
CLseed_sitesclim_Normal_1981_2010Y.csv = Climate NA 1981-2010 normals for each study site

metadata_variables.csv = description of all variables in other data frames

microclimate.RData = measurements of microclimate for each replicate plot (3 per site).

occur_biome_clim_Normal_1981_2010Y.RData = climate NA 1981-2010 normals for each GBIF adult occurrence record 

recruit_census.csv = yearly seedling census data at the replicate level for seed addition experiment (3 reps per site)

recruit_census.RData = yearly seedling census data at the replicate level for seed addition experiment (3 reps per site)

recruit_surv1yr.RData = data on whether or not recruits that germinated survived to their second year

recruit_surv2yrs.RData = data on whether or not recruits that germinated survived to their third year

recruit_survto2022.RData = data on whether or not recruits survived to 2022 (the end of the experiment)

siterangepositions.csv = description of where each seed addition site falls within each focal species climatic range

species_desc.csv = description of focal species

year2size.Rdata = data on mean height of each species in each plot in the first year of growth


# Scripts
1_climatic_range_position.R = quantify species climatic ranges from occurrence records

2_recruitment_optima_analysis.R = modelling where optimal recruitment is across climatic ranges

3_recruitment_edge_analysis.R = modelling recruitment at and beyond cold, war, dry, and wet edges

4_canopy_opt_analysis.R = modelling whether optimal recruitment changes with canopy cover

5_surv_height_analysis.R = modelling survival and growth following initial recruitment

6_speciesspecific.R = species-specific recruitment across the range

7_pretty_plots.R = putting all the figures together

bootstraps/ = folders of bootstrapping code to generate 95% confidence intervals for models to be run on a separate cluster server, since they take >24 hours.


# Outputs

/bootstraps = bootstrap outputs used to calculate 95% CIs includes predictions and, where relevant, optima for each bootstraps

/final = where final figures are saved

/model =  model predictions and 95% confidence intervals from analysis to be inputted into 7_pretty_plots.R to generate figures