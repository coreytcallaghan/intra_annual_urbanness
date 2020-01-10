## This is an R script to
## run the analysis on only species
## which have 'complete' data
## using a phylogenetic approach
## these analyses will be compared with the imputed data analysis as well
## the first part of this script is the same as the 
## script "complete_species_non_phylo.R"


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
library(ape)
library(phylolm)
library(phylosignal)
library(phylobase)

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
# I will do some imputation for the missing species
# but in this R script lets just run the models for species which have data for each category
# that I am interested in
# so will first filter to the columns of interest and then
# filter out all species without complete data
analysis <- analysis %>%
  ungroup() %>%
  dplyr::select(1:11, brain_residual, diet_breadth, adult_body_mass_g, functional_diet,
                mean_flock_size, total_range_km2, habitat_generalism_scaled, clutch_size, migration_status) %>%
  dplyr::filter(complete.cases(.))

# now test the tiplabel and common name lengths
length(unique(analysis$COMMON_NAME))

# left with 245 species which have complete data that can be analyzed
# let's start by looking at correlation and distribution of response variable
hist(analysis$mean_urbanness)

# doesn't look like too much colinearity problems
# some of the variables look pretty skewed, unsurprisingly
# so will need to log-transform these and they should look pretty normal
analysis <- analysis %>%
  mutate(response=log(mean_urbanness)) %>%
  mutate(log_body_size=log(adult_body_mass_g)) %>%
  mutate(log_flock_size=log(mean_flock_size)) %>%
  mutate(log_range_size=log(total_range_km2)) %>%
  mutate(weights=1/(sd_urbanness+0.00001))


##########################################
##########################################
#### Now start phylo stuff ###############

# function to read one tree in
read_one_tree<-function(path, x=1){
  
  one_bird_tree <- ape::read.tree(file = "Data/phylo/phy.tre")[[x]]

  return(one_bird_tree)
}

bird_tree <- read_one_tree()


# function to read all trees in
read_all_trees<-function(path){
  
  ape::read.tree(file = "Data/phylo/phy.tre")
  
}

#all_tress <- read_all_trees()

# a function to subset the tree to the tips of the 245 species
# described above
subset_tree <- function(bird_tree, dataset) {
  
  non_usa_sp <- bird_tree$tip.label[!bird_tree$tip.label %in% dataset$TipLabel]
  
  usa_bird_tree <- drop.tip(bird_tree, non_usa_sp)
  
  return(usa_bird_tree)
}

usa_tree <- subset_tree(bird_tree, analysis)

# a function to run many phylo models
run_many_phylo_models <- function(data, all_trees, n=1000){
  
  bird_trees_ss <- all_trees[1:n]
  
  ss_trees <- lapply(bird_trees_ss, subset_tree, data)
  
  list_o <- lapply(ss_trees, run_one_phylo_model, data=analysis)
  
  return(list_o)
}

# a function to extract a term
# in this case, brain residual
extract_brain <- function(mod, term_to_extract="brain_residual"){
  
  nn <-names(mod$coefficients)
  
  return(mod$coefficients[nn==term_to_extract])
  
}

# function to run a phylo model
# on non-scaled data
run_one_phylo_model <- function(usa_tree, analysis){
  
  row.names(analysis) <- analysis$TipLabel
  
  phy_mod <- phylolm(response ~ log_body_size + log_flock_size + log_range_size + brain_residual +
                       clutch_size + habitat_generalism_scaled + diet_breadth,
                     data=analysis, phy=usa_tree, na.action="na.fail", weights=weights)
  
  return(phy_mod)
}

# run a phylo model, but on standardized data
# using arm::rescale
standard_phylo_model <- function(usa_tree, analysis) {
  
  row.names(analysis) <- analysis$TipLabel
  
  phy_mod_rescaled <- phylolm(response ~ rescale(log_body_size) + rescale(log_flock_size) + 
                                rescale(log_range_size) + rescale(brain_residual) +
                                rescale(clutch_size) + rescale(habitat_generalism_scaled) + 
                                rescale(diet_breadth) + functional_diet,
                              data=analysis, phy=usa_tree, na.action="na.fail", weights=weights)
  
  return(phy_mod_rescaled)
  
}

## phylosignal analysis
phylosignal_analysis <- function(analysis, usa_tree){

  distinct_dat <- analysis %>%
    dplyr::select(COMMON_NAME, TipLabel, log_body_size, log_flock_size, 
                  log_range_size, brain_residual,
                  clutch_size, habitat_generalism_scaled, diet_breadth) %>%
    distinct()
  
  dat <- distinct_dat %>%
    dplyr::select(log_body_size, log_flock_size, log_range_size, brain_residual,
                    clutch_size, habitat_generalism_scaled, diet_breadth)
  
  row.names(dat) <- distinct_dat$TipLabel #name rows so that it matches the tree
  
  #dd$rand<-rnorm(dim(dd)[2]) #random numbers to test package
  
  p4d <- phylo4d(usa_tree, dat) #create phylobase object
  
  ps <- phyloSignal(p4d,reps = 9999) #run calculation, p values a bit unstable at 999 reps
  
  stats <- ps$stat %>%
    rownames_to_column(var="term") %>%
    mutate(value="Statistic")
  
  p_values <- ps$pvalue %>%
    rownames_to_column(var="term") %>%
    mutate(value="P value")
  
  phylo_sig_results <- bind_rows(stats, p_values)
  
}

phylosignal_results <- phylosignal_analysis(analysis, usa_tree)

saveRDS(phylosignal_results, "Results/complete_phylosignal_analysis.RDS")


###################################################################
###################################################################
############## Now start some of the modelling to look at phylogenetic model responses
## arm package is really tricky and won't be recognized to load
## into the cluster export. A similar issue to what we had in the
## remake workflow. Here, I pulled the arm::rescale function out and put
## it into the environment and then load it onto the clusterexport function
rescale <- function (x, binary.inputs = "center") 
{
  if (!is.numeric(x)) {
    x <- as.numeric(factor(x))
    x.obs <- x[!is.na(x)]
  }
  x.obs <- x[!is.na(x)]
  if (length(unique(x.obs)) == 2) {
    if (binary.inputs == "0/1") {
      x <- (x - min(x.obs))/(max(x.obs) - min(x.obs))
      return(x)
    }
    else if (binary.inputs == "-0.5,0.5") {
      return(x - 0.5)
    }
    else if (binary.inputs == "center") {
      return(x - mean(x.obs))
    }
    else if (binary.inputs == "full") {
      return((x - mean(x.obs))/(2 * sd(x.obs)))
    }
  }
  else {
    return((x - mean(x.obs))/(2 * sd(x.obs)))
  }
}

# alright - it generally looks pretty good on correlation and other things
# especially based on the other script
# I have now added the necessary function above
# now I will put these into one function
# which will allow me to
# run a separate model for every month
# so will need to functionalize it
# linear model function
modelling_function <- function(month){
  
  # filter data for a given month
  dat <- analysis %>%
    dplyr::filter(MONTH==month)
  
  # run a global phylogenetic model
  # on this filtered data
  phy_mod_rescaled <- standard_phylo_model(usa_tree, dat)
  
  # summarize the results of the model
  # I don't think there is compatibility with tidy
  # so will do this manually
  glob_mod_results <- as.data.frame(summary(phy_mod_rescaled)$coefficients) %>%
    cbind(confint(phy_mod_rescaled)) %>%
    rownames_to_column(var="term") %>%
    rename(estimate=Estimate) %>%
    rename(std.error=StdErr) %>%
    rename(lwr_95_confint=`2.5 %`) %>%
    rename(upr_95_confint=`97.5 %`) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="phylo_global_mod")
  
  # now do a model averaging approach for the phylogenetic model
  row.names(dat) <- dat$TipLabel
  
  phy_mod_rescaled <- phylolm(response ~ rescale(log_body_size) + rescale(log_flock_size) + 
                                rescale(log_range_size) + rescale(brain_residual) +
                                rescale(clutch_size) + rescale(habitat_generalism_scaled) + 
                                rescale(diet_breadth) + functional_diet,
                              data=dat, phy=usa_tree, na.action="na.fail", weights=weights)
  
  
  
  clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
  clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
  
  clusterExport(clust, c("dat", "phylolm", "rescale", "usa_tree"), envir = .GlobalEnv)
  
  # now uses the dredge function from MuMIn
  # but uses the 'pdredge' version which is for paralellizing
  model.set <- pdredge(phy_mod_rescaled, m.lim=c(0, 8), cluster=clust)
  
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
    rename(z_value=`z value`) %>%
    rename(p.value=`Pr(>|z|)`) %>%
    rename(lwr_95_confint=`2.5 %`) %>%
    rename(upr_95_confint=`97.5 %`) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="model_averaging")
  
  modelling_results_list <- list(glob_mod_results, model_results)
  
}

# now apply this function for every month
# so we get 12 sets of 'results'

results_list <- lapply(unique(analysis$MONTH), function(x){modelling_function(x)})

# now just need to get the three different
# dataframes in each list of lists
# out into a fashion to work with! 
# this is a lazy way to deal with it
# reflecting my lack of coding ability in regards to lists of lists
phylo_global_model_results <- bind_rows(results_list[[1]][[1]],
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

phylo_model_averaging_results <- bind_rows(results_list[[1]][[2]],
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

rm(list=setdiff(ls(), c("phylo_global_model_results", "phylo_model_averaging_results")))

save.image("Results/complete_species_phylo.RData")



###############################################################################
###############################################################################
###############################################################################
###### Repeat the analysis but x2 by migrants and residents
###############################################################################
###############################################################################
###############################################################################
# migrants
migrant_dat <- analysis %>%
  dplyr::filter(migration_status=="Migrant")

modelling_function <- function(month){
  
  tree <- subset_tree(bird_tree, migrant_dat)
  
  # filter data for a given month
  dat <- migrant_dat %>%
    dplyr::filter(MONTH==month)
  
  # run a global phylogenetic model
  # on this filtered data
  phy_mod_rescaled <- standard_phylo_model(tree, dat)
  
  # summarize the results of the model
  # I don't think there is compatibility with tidy
  # so will do this manually
  glob_mod_results <- as.data.frame(summary(phy_mod_rescaled)$coefficients) %>%
    cbind(confint(phy_mod_rescaled)) %>%
    rownames_to_column(var="term") %>%
    rename(estimate=Estimate) %>%
    rename(std.error=StdErr) %>%
    rename(lwr_95_confint=`2.5 %`) %>%
    rename(upr_95_confint=`97.5 %`) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="phylo_global_mod")
  
  # now do a model averaging approach for the phylogenetic model
  row.names(dat) <- dat$TipLabel
  
  phy_mod_rescaled <- phylolm(response ~ rescale(log_body_size) + rescale(log_flock_size) + 
                                rescale(log_range_size) + rescale(brain_residual) +
                                rescale(clutch_size) + rescale(habitat_generalism_scaled) + 
                                rescale(diet_breadth) + functional_diet,
                              data=dat, phy=tree, na.action="na.fail", weights=weights)
  
  
  
  clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
  clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
  
  clusterExport(clust, c("dat", "phylolm", "rescale", "tree"), envir = .GlobalEnv)
  
  # now uses the dredge function from MuMIn
  # but uses the 'pdredge' version which is for paralellizing
  model.set <- pdredge(phy_mod_rescaled, m.lim=c(0, 8), cluster=clust)
  
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
    rename(z_value=`z value`) %>%
    rename(p.value=`Pr(>|z|)`) %>%
    rename(lwr_95_confint=`2.5 %`) %>%
    rename(upr_95_confint=`97.5 %`) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="model_averaging")
  
  modelling_results_list <- list(glob_mod_results, model_results)
  
}

# now apply this function for every month
# so we get 12 sets of 'results'

results_list.migrants <- lapply(unique(analysis$MONTH), function(x){modelling_function(x)})

# now just need to get the three different
# dataframes in each list of lists
# out into a fashion to work with! 
# this is a lazy way to deal with it
# reflecting my lack of coding ability in regards to lists of lists
phylo_global_model_results.migrants <- bind_rows(results_list.migrants[[1]][[1]],
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

phylo_model_averaging_results.migrants <- bind_rows(results_list.migrants[[1]][[2]],
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

save(phylo_global_model_results.migrants, 
     phylo_model_averaging_results.migrants, 
     file="Results/complete_species_phylo_migrants_only.RData")


# residents
resident_dat <- analysis %>%
  dplyr::filter(migration_status=="Resident")

modelling_function <- function(month){
  
  tree <- subset_tree(bird_tree, resident_dat)
  
  # filter data for a given month
  dat <- resident_dat %>%
    dplyr::filter(MONTH==month)
  
  # run a global phylogenetic model
  # on this filtered data
  phy_mod_rescaled <- standard_phylo_model(tree, dat)
  
  # summarize the results of the model
  # I don't think there is compatibility with tidy
  # so will do this manually
  glob_mod_results <- as.data.frame(summary(phy_mod_rescaled)$coefficients) %>%
    cbind(confint(phy_mod_rescaled)) %>%
    rownames_to_column(var="term") %>%
    rename(estimate=Estimate) %>%
    rename(std.error=StdErr) %>%
    rename(lwr_95_confint=`2.5 %`) %>%
    rename(upr_95_confint=`97.5 %`) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="phylo_global_mod")
  
  # now do a model averaging approach for the phylogenetic model
  row.names(dat) <- dat$TipLabel
  
  phy_mod_rescaled <- phylolm(response ~ rescale(log_body_size) + rescale(log_flock_size) + 
                                rescale(log_range_size) + rescale(brain_residual) +
                                rescale(clutch_size) + rescale(habitat_generalism_scaled) + 
                                rescale(diet_breadth) + functional_diet,
                              data=dat, phy=tree, na.action="na.fail", weights=weights)
  
  
  
  clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
  clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
  
  clusterExport(clust, c("dat", "phylolm", "rescale", "tree"), envir = .GlobalEnv)
  
  # now uses the dredge function from MuMIn
  # but uses the 'pdredge' version which is for paralellizing
  model.set <- pdredge(phy_mod_rescaled, m.lim=c(0, 8), cluster=clust)
  
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
    rename(z_value=`z value`) %>%
    rename(p.value=`Pr(>|z|)`) %>%
    rename(lwr_95_confint=`2.5 %`) %>%
    rename(upr_95_confint=`97.5 %`) %>%
    mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
    mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
    mutate(MONTH=month) %>%
    mutate(model_type="model_averaging")
  
  modelling_results_list <- list(glob_mod_results, model_results)
  
}

# now apply this function for every month
# so we get 12 sets of 'results'

results_list.residents <- lapply(unique(analysis$MONTH), function(x){modelling_function(x)})

# now just need to get the three different
# dataframes in each list of lists
# out into a fashion to work with! 
# this is a lazy way to deal with it
# reflecting my lack of coding ability in regards to lists of lists
phylo_global_model_results.residents <- bind_rows(results_list.residents[[1]][[1]],
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

phylo_model_averaging_results.residents <- bind_rows(results_list.residents[[1]][[2]],
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

save(phylo_global_model_results.residents, 
     phylo_model_averaging_results.residents, 
     file="Results/complete_species_phylo_residents_only.RData")




