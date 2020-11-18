## This is an R script to
## run the analysis on only species
## which have 'complete' data
## using a phylogenetic approach
## these analyses will be compared with the imputed data analysis as well
## the first part of this script is the same as the 
## script "complete_species_non_phylo.R"


# packages
library(dplyr)
library(readr)
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
library(phytools)
library(patchwork)

source("R/global_functions.R")

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

# I redid analyses that will influence the number of species
# because I now only use data from 2014, so some species won't meet the 250
# threshold for all 12 months.
# get a list of these species
# and remove them from possible analyses
species_without_enough_dat <- analysis %>%
  dplyr::filter(number_obs<=250) %>%
  dplyr::select(COMMON_NAME) %>%
  distinct()

analysis <- analysis %>%
  dplyr::filter(!COMMON_NAME %in% species_without_enough_dat$COMMON_NAME)

# now test the tiplabel and common name lengths
length(unique(analysis$COMMON_NAME))
length(unique(analysis$TipLabel))

monthly_histograms <- ggplot(analysis, aes(x=mean_urbanness))+
  geom_histogram(color="black", fill="#4DAF4A", bins=50)+
  scale_x_log10()+
  facet_wrap(~MONTH)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_fill_brewer(palette="Set1")+
  xlab("Species-specific urban tolerance")+
  ylab("Number of species")+
  ggtitle("a)")

monthly_histograms

ggsave("Figures/monthly_urbanness_distributions.png", width=8.5, height=7.8, units="in")

analysis %>%
  group_by(MONTH) %>%
  summarize(mean=mean(mean_urbanness),
            sd=sd(mean_urbanness),
            N=n()) %>%
  mutate(se=sd/sqrt(N))

monthly_means <- analysis %>%
  group_by(MONTH) %>%
  summarize(mean=mean(mean_urbanness),
            sd=sd(mean_urbanness),
            N=n()) %>%
  mutate(se=sd/sqrt(N)) %>%
  ggplot(., aes(x=MONTH, y=mean, group=1))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.4)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_color_brewer(palette="Set1")+
  ylab("Mean of urban tolerances")+
  xlab("")+
  ggtitle("b)")

monthly_means

ggsave("Figures/monthly_urbanness_mean_line.png", width=7, height=5.8, units="in")


monthly_means_split <- analysis %>%
  group_by(MONTH, migration_status) %>%
  summarize(mean=mean(mean_urbanness),
            sd=sd(mean_urbanness),
            N=n()) %>%
  mutate(se=sd/sqrt(N)) %>%
  ggplot(., aes(x=MONTH, y=mean, group=migration_status, color=migration_status))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.4)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_color_brewer(palette="Set1")+
  ylab("Mean of urban tolerances")+
  xlab("")+
  ggtitle("c)")+
  labs(color="Migratory status")+
  theme(legend.position="bottom")

monthly_means_split

ggsave("Figures/monthly_urbanness_mean_line_split.png", width=7, height=5.8, units="in")



monthly_histograms + monthly_means + monthly_means_split + plot_layout(ncol=1)

ggsave("Figures/monthly_all_species_urbanness_summary.png", height=10, width=7.5, units="in")


# pick some species with different breeding periods
response %>%
  dplyr::filter(COMMON_NAME %in% c("Barred Owl", "American Goldfinch",
                                   "Northern Mockingbird", "Prairie Warbler",
                                   "Mallard", "Blackburnian Warbler",
                                   "Horned Lark", "Pine Grosbeak")) %>%
  ggplot(., aes(x=MONTH, y=urban_score, group=1))+
  geom_line()+
  geom_point()+
  facet_wrap(~COMMON_NAME, scales="free", ncol=2)+
  theme_bw()



# how many obs for results
sum(analysis$number_obs)

# doesn't look like too much colinearity problems
# some of the variables look pretty skewed, unsurprisingly
# so will need to log-transform these and they should look pretty normal
analysis <- analysis %>%
  mutate(response=log10(mean_urbanness)) %>%
  mutate(log_body_size=log(adult_body_mass_g)) %>%
  mutate(log_flock_size=log(mean_flock_size)) %>%
  mutate(log_range_size=log(total_range_km2)) %>%
  mutate(weights=1/(sd_urbanness)) %>%
  mutate(weights=ifelse(weights>50, 50, weights))


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

all_trees <- read_all_trees()

# a function to subset the tree to the tips of the 245 species
# described above
subset_tree <- function(bird_tree, dataset) {
  
  non_usa_sp <- bird_tree$tip.label[!bird_tree$tip.label %in% dataset$TipLabel]
  
  usa_bird_tree <- drop.tip(bird_tree, non_usa_sp)
  
  return(usa_bird_tree)
}

usa_tree <- subset_tree(bird_tree, analysis)

# need to get a consensus tree to run the phylogenetic analyses on
# first subset all trees to the 245 species
non_usa_sp <- bird_tree$tip.label[!bird_tree$tip.label %in% analysis$TipLabel]

subset_trees <- lapply(all_trees, drop.tip, tip=non_usa_sp)

con_tree <- consensus.edges(subset_trees,consensus.tree=consensus(subset_trees, p=0.5, check.labels=TRUE))



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
                     data=analysis, phy=con_tree, na.action="na.fail", weights=weights)
  
  return(phy_mod)
}

# run a phylo model, but on standardized data
# using arm::rescale
standard_phylo_model <- function(con_tree, analysis) {
  
  row.names(analysis) <- analysis$TipLabel
  
  phy_mod_rescaled <- phylolm(response ~ rescale(log_body_size) + rescale(log_flock_size) + 
                                rescale(log_range_size) + rescale(brain_residual) +
                                rescale(clutch_size) + rescale(habitat_generalism_scaled) + 
                                rescale(diet_breadth),
                              data=analysis, phy=con_tree, na.action="na.fail", weights=weights)
  
  return(phy_mod_rescaled)
  
}

## phylosignal analysis
## do this for every month
## because the signal for urbanness could change throughout the months
phylosig_month_function <- function(month) {
  
phylosignal_analysis <- function(analysis, con_tree){

  distinct_dat <- analysis %>%
    dplyr::select(COMMON_NAME, TipLabel, log_body_size, log_flock_size, 
                  log_range_size, brain_residual, response,
                  clutch_size, habitat_generalism_scaled, diet_breadth, MONTH) %>%
    distinct()
  
  dat <- distinct_dat %>%
    dplyr::select(log_body_size, log_flock_size, log_range_size, brain_residual, TipLabel,
                    clutch_size, habitat_generalism_scaled, diet_breadth, response, MONTH) %>%
    dplyr::filter(MONTH==month) %>%
    dplyr::select(-MONTH) %>%
    column_to_rownames(var="TipLabel")
  
  #dd$rand<-rnorm(dim(dd)[2]) #random numbers to test package
  
  p4d <- phylo4d(con_tree, dat) #create phylobase object
  
  ps <- phyloSignal(p4d,reps = 9999) #run calculation, p values a bit unstable at 999 reps
  
  stats <- ps$stat %>%
    rownames_to_column(var="term") %>%
    mutate(value="Statistic")
  
  p_values <- ps$pvalue %>%
    rownames_to_column(var="term") %>%
    mutate(value="P value")
  
  phylo_sig_results <- bind_rows(stats, p_values)
  
}

phylosignal_results <- phylosignal_analysis(analysis, con_tree) %>%
  mutate(MONTH=month)

return(phylosignal_results)
#saveRDS(phylosignal_results, paste0("Results/", month, "_complete_phylosignal_analysis.RDS"))

}


month_phylosignal_results <- bind_rows(lapply(unique(response$MONTH), phylosig_month_function))


response_phylosig <- month_phylosignal_results %>%
  dplyr::filter(term=="response") %>%
  round_df(., digits=4)

write_csv(response_phylosig, "Results/phylosignal_monthly_results.csv")

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
    dplyr::filter(MONTH==month) %>%
    column_to_rownames(var="TipLabel")
  
  # run a global phylogenetic model
  # on this filtered data
  phy_mod_rescaled <- phylolm(response ~ rescale(log_body_size) + rescale(log_flock_size) + 
                                rescale(log_range_size) + rescale(brain_residual) +
                                rescale(clutch_size) + rescale(habitat_generalism_scaled) + 
                                rescale(diet_breadth) + migration_status,
                              data=dat, phy=con_tree, na.action="na.fail", weights=weights)
  
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
  
  # now uses the dredge function from MuMIn
  # but uses the 'pdredge' version which is for paralellizing
  model.set <- dredge(phy_mod_rescaled, m.lim=c(0, 8))
  
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

save(phylo_global_model_results, 
     phylo_model_averaging_results, 
     file="Results/complete_species_phylo.RData")




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
  
  tree <- subset_tree(con_tree, migrant_dat)
  
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
                                rescale(diet_breadth),
                              data=dat, phy=tree, na.action="na.fail", weights=weights)
  
  
  
  # clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
  # clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
  # 
  # clusterExport(clust, c("dat", "phylolm", "rescale", "tree"), envir = .GlobalEnv)
  
  # now uses the dredge function from MuMIn
  # but uses the 'pdredge' version which is for paralellizing
  model.set <- dredge(phy_mod_rescaled, m.lim=c(0, 8))
  
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
  
  tree <- subset_tree(con_tree, resident_dat)
  
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
                                rescale(diet_breadth),
                              data=dat, phy=tree, na.action="na.fail", weights=weights)
  
  
  
  # clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
  # clust <- try(makeCluster(getOption("cl.cores", 10), type = clusterType))
  # 
  # clusterExport(clust, c("dat", "phylolm", "rescale", "tree"), envir = .GlobalEnv)
  
  # now uses the dredge function from MuMIn
  # but uses the 'pdredge' version which is for paralellizing
  model.set <- dredge(phy_mod_rescaled, m.lim=c(0, 8))
  
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




