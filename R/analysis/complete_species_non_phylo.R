## This is an R script to
## run the analysis on only species
## which have 'complete' data
## non phylogenetically
## these analyses will be compared with the imputed data analysis as well


# packages
library(dplyr)
library(ggplot2)
library(car)
library(ggcorrplot)
library(GGally)
library(tidyr)
library(broom)
library(tibble)
library(arm)
library(snow)
library(MuMIn)

# read in response variables
response <- readRDS("Data/response_variables.RDS")

# read in predictor variables
predictors <- readRDS("Data/predictor_variables.RDS")

# join data for 'analysis' dat
analysis <- response %>%
  left_join(., predictors)

# look at number of species
length(unique(analysis$COMMON_NAME))
length(unique(analysis$TipLabel))

# so there are four more 'species' than there are phylogenetic tips
# this isn't actually important for this part of the analysis
# but because we want to phylogenetically constrain the analysis later on
# will clean this up and select the species
# to include
more_than_1 <- analysis %>%
  group_by(TipLabel) %>%
  summarize(number_tips=length(unique(COMMON_NAME))) %>%
  dplyr::filter(number_tips > 1)

# show the species which have the same TipLabel  
analysis %>%
  dplyr::filter(TipLabel %in% more_than_1$TipLabel) %>%
  distinct(COMMON_NAME) %>%
  .$COMMON_NAME

# all of these are examples of species which have recently been split
# before Walter Jetz's tree etc.
# which explains why they have the same TipLabel but are treated as different
# by Clements/eBird
# in this instance, I'll select the nominate species for each of these pairs
# so I'll filter out those not of interest here
analysis <- analysis %>%
  dplyr::filter(!COMMON_NAME %in% c("Pacific Wren", "Ridgway's Rail", "Woodhouse's Scrub-Jay", "Bell's Sparrow"))

# now test the tiplabel and common name lengths
length(unique(analysis$COMMON_NAME))
length(unique(analysis$TipLabel))

# now they match! So we have a potential sample size of 486 species
# I might do some imputation for the missing species
# but in this R script lets just run the models for species which have data for each category
# that I am interested in
# so will first filter to the columns of interest and then
# filter out all species without complete data
analysis <- analysis %>%
  ungroup() %>%
  dplyr::select(1:10, brain_residual, diet_breadth, adult_body_mass_g, functional_diet,
                mean_flock_size, total_range_km2, habitat_generalism_scaled, clutch_size, migration_status) %>%
  dplyr::filter(complete.cases(.))

# now test the tiplabel and common name lengths
length(unique(analysis$COMMON_NAME))

# left with 245 species which have complete data that can be analyzed
# let's start by looking at correlation and distribution of response variable
hist(analysis$mean_urbanness)

# heavily skewed so will try log transform
hist(log(analysis$mean_urbanness))

# plot ggpairs of variables
analysis %>%
  dplyr::select(11:18) %>%
  ggpairs()

# doesn't look like too much colinearity problems
# some of the variables look pretty skewed, unsurprisingly
# so will need to log-transform these and they should look pretty normal
analysis <- analysis %>%
  mutate(response=log(mean_urbanness)) %>%
  mutate(log_body_size=log(adult_body_mass_g)) %>%
  mutate(log_flock_size=log(mean_flock_size)) %>%
  mutate(log_range_size=log(total_range_km2)) %>%
  mutate(weights=1/(sd_urbanness+0.00001))
           
# make a correlation plot figure
analysis %>%
  dplyr::select(log_body_size, log_flock_size, log_range_size, brain_residual,
                clutch_size, habitat_generalism_scaled, diet_breadth) %>%
  rename(`Body size (log)`=log_body_size) %>%
  rename(`Mean flock size (log)`=log_flock_size) %>%
  rename(`Range size (log km2)`=log_range_size) %>%
  rename(`Brain residual`=brain_residual) %>%
  rename(`Clutch size`=clutch_size) %>%
  rename(`Habitat generalism`=habitat_generalism_scaled) %>%
  rename(`Diet breadth`=diet_breadth) %>%
  cor(., use="pairwise.complete.obs") %>%
  ggcorrplot(lab=TRUE, 
             outline.col = "white",
             ggtheme = ggplot2::theme_classic,
             colors = c("#6D9EC1", "white", "#E46726"))+
  theme(axis.text=element_text(color="black"))

ggsave("Figures/correlation_among_predictors.png", width=6, height=7, units="in")

# alright - it generally looks pretty good
# now need to think about putting the model together
# but I want to run a separate model for every month
# so will need to functionalize it
# linear model function
modelling_function <- function(month) {
  
  # filter data for a given month
  dat <- analysis %>%
    dplyr::filter(MONTH==month) %>%
    mutate(z.log_body_size=rescale(log_body_size)) %>%
    mutate(z.log_flock_size=rescale(log_flock_size)) %>%
    mutate(z.log_range_size=rescale(log_range_size)) %>%
    mutate(z.brain_residual=rescale(brain_residual)) %>%
    mutate(z.clutch_size=rescale(clutch_size)) %>%
    mutate(z.habitat_generalism_scaled=rescale(habitat_generalism_scaled)) %>%
    mutate(z.diet_breadth=rescale(diet_breadth))
  
  # fit a linear global model
  # standardize the model
  lm.mod <- lm(response ~ z.log_body_size + z.log_flock_size + z.log_range_size + 
                 z.brain_residual + z.clutch_size + z.habitat_generalism_scaled + 
                 z.diet_breadth + functional_diet,
               data=dat, na.action="na.fail", weights=weights)
  
  # get a dataframe of variance inflation factors
  # from the arm package
  # which shows how (if) the collinearity among predictors influences the model results
  vif_df <- as.data.frame(vif(lm.mod)) %>%
    rownames_to_column(var="term") %>%
    rename(adjusted_GVIF=4) %>%
    mutate(MONTH=month) %>%
    mutate(adjusted_GVIF=round(adjusted_GVIF, digits=2))
  
  # now create a 'summary' dataframe
  summary_df <- tidy(lm.mod) %>%
    mutate(lwr_95_confint=confint(lm.mod)[,1]) %>%
    mutate(upr_95_confint=confint(lm.mod)[,2]) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="global_model")

  ### Now prepare for model averaging of the global model
  # this part just sets up a parallelizing which isn't necessary
  # but is increasingly necessary as more variables are included
  clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
  clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
  
  clusterExport(clust, "dat", envir = environment())
  
  # now uses the dredge function from MuMIn
  # but uses the 'pdredge' version which is for paralellizing
  model.set <- pdredge(lm.mod, m.lim=c(0, 8), cluster=clust, extra="R^2")
  
  # selects all models with deltaAic < 4
  top.models <- get.models(model.set, subset=delta<4) 
  
  # how many top models
  length(top.models)
  
  # Ranks these models based on AICc
  my.models <- model.sel(top.models, rank="AICc") 
  
  # actually does the model averaging part of the analysis
  averaged_models <- model.avg(my.models)
  
  # now get a summary of the averaged model results
  # using the 'full' set of models
  # not conditional
  # more details on that are here: https://onlinelibrary.wiley.com/doi/10.1111/j.1420-9101.2010.02210.x
  model_results <- as.data.frame(summary(averaged_models)$coefmat.full) %>%
    cbind(confint(averaged_models, full=TRUE)) %>%
    rownames_to_column(var="term") %>%
    left_join(., data.frame(importance=averaged_models$importance) %>%
                rownames_to_column(var="term"), by="term") %>%
    rename(estimate=Estimate) %>%
    rename(std.error=`Std. Error`) %>%
    rename(adjusted_std.error=`Adjusted SE`) %>%
    rename(z_value=`z value`) %>%
    rename(p.value=`Pr(>|z|)`) %>%
    rename(lwr_95_confint=`2.5 %`) %>%
    rename(upr_95_confint=`97.5 %`) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="model_averaging")
  
  modelling_results_list <- list(vif_df, summary_df, model_results)
  
}

# now apply this function for every month
# so we get 12 sets of 'results'

results_list <- lapply(as.character(unique(analysis$MONTH)), function(x){modelling_function(x)})

# now just need to get the three different
# dataframes in each list of lists
# out into a fashion to work with! 
# this is a lazy way to deal with it
# reflecting my lack of coding ability in regards to lists of lists
vif_summary <- bind_rows(results_list[[1]][[1]],
                         results_list[[2]][[1]],
                         results_list[[3]][[1]],
                         results_list[[4]][[1]],
                         results_list[[5]][[1]],
                         results_list[[6]][[1]],
                         results_list[[7]][[1]],
                         results_list[[8]][[1]],
                         results_list[[9]][[1]],
                         results_list[[10]][[1]],
                         results_list[[11]][[1]],
                         results_list[[12]][[1]])

non_phylo_global_model_results <- bind_rows(results_list[[1]][[2]],
                                            results_list[[2]][[2]],
                                            results_list[[3]][[2]],
                                            results_list[[4]][[2]],
                                            results_list[[5]][[2]],
                                            results_list[[6]][[2]],
                                            results_list[[7]][[2]],
                                            results_list[[8]][[2]],
                                            results_list[[9]][[2]],
                                            results_list[[10]][[2]],
                                            results_list[[11]][[2]],
                                            results_list[[12]][[2]])

non_phylo_model_averaging_results <- bind_rows(results_list[[1]][[3]],
                                               results_list[[2]][[3]],
                                               results_list[[3]][[3]],
                                               results_list[[4]][[3]],
                                               results_list[[5]][[3]],
                                               results_list[[6]][[3]],
                                               results_list[[7]][[3]],
                                               results_list[[8]][[3]],
                                               results_list[[9]][[3]],
                                               results_list[[10]][[3]],
                                               results_list[[11]][[3]],
                                               results_list[[12]][[3]])

rm(list=setdiff(ls(), c("vif_summary", "non_phylo_global_model_results", "non_phylo_model_averaging_results")))

save.image("Results/complete_species_non_phylo.RData")


## Now want to repeat the modelling function above, but for migrants and residents
migrant_dat <- analysis %>%
  dplyr::filter(migration_status=="Migrant")

modelling_function <- function(month) {
  
  # filter data for a given month
  dat <- migrant_dat %>%
    dplyr::filter(MONTH==month) %>%
    mutate(z.log_body_size=rescale(log_body_size)) %>%
    mutate(z.log_flock_size=rescale(log_flock_size)) %>%
    mutate(z.log_range_size=rescale(log_range_size)) %>%
    mutate(z.brain_residual=rescale(brain_residual)) %>%
    mutate(z.clutch_size=rescale(clutch_size)) %>%
    mutate(z.habitat_generalism_scaled=rescale(habitat_generalism_scaled)) %>%
    mutate(z.diet_breadth=rescale(diet_breadth))
  
  # fit a linear global model
  # standardize the model
  lm.mod <- lm(response ~ z.log_body_size + z.log_flock_size + z.log_range_size + 
                 z.brain_residual + z.clutch_size + z.habitat_generalism_scaled + 
                 z.diet_breadth + functional_diet,
               data=dat, na.action="na.fail", weights=weights)
  
  # get a dataframe of variance inflation factors
  # from the arm package
  # which shows how (if) the collinearity among predictors influences the model results
  vif_df <- as.data.frame(vif(lm.mod)) %>%
    rownames_to_column(var="term") %>%
    rename(adjusted_GVIF=4) %>%
    mutate(MONTH=month) %>%
    mutate(adjusted_GVIF=round(adjusted_GVIF, digits=2))
  
  # now create a 'summary' dataframe
  summary_df <- tidy(lm.mod) %>%
    mutate(lwr_95_confint=confint(lm.mod)[,1]) %>%
    mutate(upr_95_confint=confint(lm.mod)[,2]) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="global_model")
  
  ### Now prepare for model averaging of the global model
  # this part just sets up a parallelizing which isn't necessary
  # but is increasingly necessary as more variables are included
  clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
  clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
  
  clusterExport(clust, "dat", envir = environment())
  
  # now uses the dredge function from MuMIn
  # but uses the 'pdredge' version which is for paralellizing
  model.set <- pdredge(lm.mod, m.lim=c(0, 8), cluster=clust, extra="R^2")
  
  # selects all models with deltaAic < 4
  top.models <- get.models(model.set, subset=delta<4) 
  
  # how many top models
  length(top.models)
  
  # Ranks these models based on AICc
  my.models <- model.sel(top.models, rank="AICc") 
  
  # actually does the model averaging part of the analysis
  averaged_models <- model.avg(my.models)
  
  # now get a summary of the averaged model results
  # using the 'full' set of models
  # not conditional
  # more details on that are here: https://onlinelibrary.wiley.com/doi/10.1111/j.1420-9101.2010.02210.x
  model_results <- as.data.frame(summary(averaged_models)$coefmat.full) %>%
    cbind(confint(averaged_models, full=TRUE)) %>%
    rownames_to_column(var="term") %>%
    left_join(., data.frame(importance=averaged_models$importance) %>%
                rownames_to_column(var="term"), by="term") %>%
    rename(estimate=Estimate) %>%
    rename(std.error=`Std. Error`) %>%
    rename(adjusted_std.error=`Adjusted SE`) %>%
    rename(z_value=`z value`) %>%
    rename(p.value=`Pr(>|z|)`) %>%
    rename(lwr_95_confint=`2.5 %`) %>%
    rename(upr_95_confint=`97.5 %`) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="model_averaging")
  
  modelling_results_list <- list(vif_df, summary_df, model_results)
  
}

# now apply this function for every month
# so we get 12 sets of 'results'

results_list.migrants <- lapply(as.character(unique(analysis$MONTH)), function(x){modelling_function(x)})

# now just need to get the three different
# dataframes in each list of lists
# out into a fashion to work with! 
# this is a lazy way to deal with it
# reflecting my lack of coding ability in regards to lists of lists
vif_summary.migrants <- bind_rows(results_list.migrants[[1]][[1]],
                                  results_list.migrants[[2]][[1]],
                                  results_list.migrants[[3]][[1]],
                                  results_list.migrants[[4]][[1]],
                                  results_list.migrants[[5]][[1]],
                                  results_list.migrants[[6]][[1]],
                                  results_list.migrants[[7]][[1]],
                                  results_list.migrants[[8]][[1]],
                                  results_list.migrants[[9]][[1]],
                                  results_list.migrants[[10]][[1]],
                                  results_list.migrants[[11]][[1]],
                                  results_list.migrants[[12]][[1]])

non_phylo_global_model_results.migrants <- bind_rows(results_list.migrants[[1]][[2]],
                                                     results_list.migrants[[2]][[2]],
                                                     results_list.migrants[[3]][[2]],
                                                     results_list.migrants[[4]][[2]],
                                                     results_list.migrants[[5]][[2]],
                                                     results_list.migrants[[6]][[2]],
                                                     results_list.migrants[[7]][[2]],
                                                     results_list.migrants[[8]][[2]],
                                                     results_list.migrants[[9]][[2]],
                                                     results_list.migrants[[10]][[2]],
                                                     results_list.migrants[[11]][[2]],
                                                     results_list.migrants[[12]][[2]])

non_phylo_model_averaging_results.migrants <- bind_rows(results_list.migrants[[1]][[3]],
                                                        results_list.migrants[[2]][[3]],
                                                        results_list.migrants[[3]][[3]],
                                                        results_list.migrants[[4]][[3]],
                                                        results_list.migrants[[5]][[3]],
                                                        results_list.migrants[[6]][[3]],
                                                        results_list.migrants[[7]][[3]],
                                                        results_list.migrants[[8]][[3]],
                                                        results_list.migrants[[9]][[3]],
                                                        results_list.migrants[[10]][[3]],
                                                        results_list.migrants[[11]][[3]],
                                                        results_list.migrants[[12]][[3]])

save(non_phylo_global_model_results.migrants, 
     non_phylo_model_averaging_results.migrants,
     vif_summary.migrants, file="Results/complete_species_non_phylo_migrants_only.RData")

## Now want to repeat the modelling function above, but for migrants and residents
resident_dat <- analysis %>%
  dplyr::filter(migration_status=="Resident")

modelling_function <- function(month) {
  
  # filter data for a given month
  dat <- resident_dat %>%
    dplyr::filter(MONTH==month) %>%
    mutate(z.log_body_size=rescale(log_body_size)) %>%
    mutate(z.log_flock_size=rescale(log_flock_size)) %>%
    mutate(z.log_range_size=rescale(log_range_size)) %>%
    mutate(z.brain_residual=rescale(brain_residual)) %>%
    mutate(z.clutch_size=rescale(clutch_size)) %>%
    mutate(z.habitat_generalism_scaled=rescale(habitat_generalism_scaled)) %>%
    mutate(z.diet_breadth=rescale(diet_breadth))
  
  # fit a linear global model
  # standardize the model
  lm.mod <- lm(response ~ z.log_body_size + z.log_flock_size + z.log_range_size + 
                 z.brain_residual + z.clutch_size + z.habitat_generalism_scaled + 
                 z.diet_breadth + functional_diet,
               data=dat, na.action="na.fail", weights=weights)
  
  # get a dataframe of variance inflation factors
  # from the arm package
  # which shows how (if) the collinearity among predictors influences the model results
  vif_df <- as.data.frame(vif(lm.mod)) %>%
    rownames_to_column(var="term") %>%
    rename(adjusted_GVIF=4) %>%
    mutate(MONTH=month) %>%
    mutate(adjusted_GVIF=round(adjusted_GVIF, digits=2))
  
  # now create a 'summary' dataframe
  summary_df <- tidy(lm.mod) %>%
    mutate(lwr_95_confint=confint(lm.mod)[,1]) %>%
    mutate(upr_95_confint=confint(lm.mod)[,2]) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="global_model")
  
  ### Now prepare for model averaging of the global model
  # this part just sets up a parallelizing which isn't necessary
  # but is increasingly necessary as more variables are included
  clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
  clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
  
  clusterExport(clust, "dat", envir = environment())
  
  # now uses the dredge function from MuMIn
  # but uses the 'pdredge' version which is for paralellizing
  model.set <- pdredge(lm.mod, m.lim=c(0, 8), cluster=clust, extra="R^2")
  
  # selects all models with deltaAic < 4
  top.models <- get.models(model.set, subset=delta<4) 
  
  # how many top models
  length(top.models)
  
  # Ranks these models based on AICc
  my.models <- model.sel(top.models, rank="AICc") 
  
  # actually does the model averaging part of the analysis
  averaged_models <- model.avg(my.models)
  
  # now get a summary of the averaged model results
  # using the 'full' set of models
  # not conditional
  # more details on that are here: https://onlinelibrary.wiley.com/doi/10.1111/j.1420-9101.2010.02210.x
  model_results <- as.data.frame(summary(averaged_models)$coefmat.full) %>%
    cbind(confint(averaged_models, full=TRUE)) %>%
    rownames_to_column(var="term") %>%
    left_join(., data.frame(importance=averaged_models$importance) %>%
                rownames_to_column(var="term"), by="term") %>%
    rename(estimate=Estimate) %>%
    rename(std.error=`Std. Error`) %>%
    rename(adjusted_std.error=`Adjusted SE`) %>%
    rename(z_value=`z value`) %>%
    rename(p.value=`Pr(>|z|)`) %>%
    rename(lwr_95_confint=`2.5 %`) %>%
    rename(upr_95_confint=`97.5 %`) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="model_averaging")
  
  modelling_results_list <- list(vif_df, summary_df, model_results)
  
}

# now apply this function for every month
# so we get 12 sets of 'results'

results_list.residents <- lapply(as.character(unique(analysis$MONTH)), function(x){modelling_function(x)})

# now just need to get the three different
# dataframes in each list of lists
# out into a fashion to work with! 
# this is a lazy way to deal with it
# reflecting my lack of coding ability in regards to lists of lists
vif_summary.residents <- bind_rows(results_list.residents[[1]][[1]],
                                  results_list.residents[[2]][[1]],
                                  results_list.residents[[3]][[1]],
                                  results_list.residents[[4]][[1]],
                                  results_list.residents[[5]][[1]],
                                  results_list.residents[[6]][[1]],
                                  results_list.residents[[7]][[1]],
                                  results_list.residents[[8]][[1]],
                                  results_list.residents[[9]][[1]],
                                  results_list.residents[[10]][[1]],
                                  results_list.residents[[11]][[1]],
                                  results_list.residents[[12]][[1]])

non_phylo_global_model_results.residents <- bind_rows(results_list.residents[[1]][[2]],
                                                     results_list.residents[[2]][[2]],
                                                     results_list.residents[[3]][[2]],
                                                     results_list.residents[[4]][[2]],
                                                     results_list.residents[[5]][[2]],
                                                     results_list.residents[[6]][[2]],
                                                     results_list.residents[[7]][[2]],
                                                     results_list.residents[[8]][[2]],
                                                     results_list.residents[[9]][[2]],
                                                     results_list.residents[[10]][[2]],
                                                     results_list.residents[[11]][[2]],
                                                     results_list.residents[[12]][[2]])

non_phylo_model_averaging_results.residents <- bind_rows(results_list.residents[[1]][[3]],
                                                        results_list.residents[[2]][[3]],
                                                        results_list.residents[[3]][[3]],
                                                        results_list.residents[[4]][[3]],
                                                        results_list.residents[[5]][[3]],
                                                        results_list.residents[[6]][[3]],
                                                        results_list.residents[[7]][[3]],
                                                        results_list.residents[[8]][[3]],
                                                        results_list.residents[[9]][[3]],
                                                        results_list.residents[[10]][[3]],
                                                        results_list.residents[[11]][[3]],
                                                        results_list.residents[[12]][[3]])

save(non_phylo_global_model_results.residents, 
     non_phylo_model_averaging_results.residents,
     vif_summary.residents, file="Results/complete_species_non_phylo_residents_only.RData")
