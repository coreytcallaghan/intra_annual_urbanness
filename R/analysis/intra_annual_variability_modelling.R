## This is an R script to investigate
## what traits predict the standard deviation of the urbanness
## throughout the year
## will start with the dataset from the other
## 'analysis' scripts
## copy and pasting the initial part

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
library(forcats)
library(ggtree)
library(phylolm)
library(phylosignal)
library(phylobase)
library(ape)
library(phylosignal)

# read in response variables
response <- readRDS("Data/response_variables.RDS") %>%
  group_by(COMMON_NAME) %>%
  summarize(intra_annual_variability=sd(mean_urbanness),
            sd_of_intra_annual_variability=sd(sd_urbanness)) 

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
  dplyr::select(1:7, brain_residual, diet_breadth, adult_body_mass_g, functional_diet,
                mean_flock_size, total_range_km2, habitat_generalism_scaled, clutch_size,
                migration_status, migratory_status) %>%
  dplyr::filter(complete.cases(.))

# now test the tiplabel and common name lengths
length(unique(analysis$COMMON_NAME))

# left with 245 species which have complete data that can be analyzed
# let's start by looking at correlation and distribution of response variable
hist(analysis$intra_annual_variability)

# heavily skewed so will try log transform
hist(log(analysis$intra_annual_variability))

# plot ggpairs of variables
analysis %>%
  dplyr::select(8:17) %>%
  ggpairs()

# doesn't look like too much colinearity problems
# some of the variables look pretty skewed, unsurprisingly
# so will need to log-transform these and they should look pretty normal
analysis <- analysis %>%
  mutate(response=log(intra_annual_variability)) %>%
  mutate(log_body_size=log(adult_body_mass_g)) %>%
  mutate(log_flock_size=log(mean_flock_size)) %>%
  mutate(log_range_size=log(total_range_km2)) %>%
  mutate(weights=1/sd_of_intra_annual_variability)


# get dat for model
# and rescale variables
dat <- analysis %>%
  mutate(z.log_body_size=rescale(log_body_size)) %>%
  mutate(z.log_flock_size=rescale(log_flock_size)) %>%
  mutate(z.log_range_size=rescale(log_range_size)) %>%
  mutate(z.brain_residual=rescale(brain_residual)) %>%
  mutate(z.clutch_size=rescale(clutch_size)) %>%
  mutate(z.habitat_generalism_scaled=rescale(habitat_generalism_scaled)) %>%
  mutate(z.diet_breadth=rescale(diet_breadth))


ggplot(dat, aes(x=migration_status, y=intra_annual_variability))+
  geom_violin(position=position_dodge()) +
  geom_boxplot(width=0.1, color="black", position = position_dodge(width =0.9))+
  coord_flip()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Migration status")+
  ylab("Intra-annual urbanness variability")

# fit a linear global model
# standardize the model
simple_mig_mod <- lm(response ~ migration_status, data=dat, na.action="na.fail", weights=weights)
summary(simple_mig_mod)
AIC(simple_mig_mod)

complex_mig_mod <- lm(response ~ migratory_status, data=dat, na.action="na.fail", weights=weights)
summary(complex_mig_mod)
AIC(complex_mig_mod)


# since migration status is the simpler version here
# I will stick with this for now, and the change in AIC (above)
# is pretty marginal!
# I also might split the other models (by month)
# by migration status
# so makes sense to do this in a simpler version
lm.mod <- lm(response ~ migration_status + z.log_body_size + z.log_flock_size + z.log_range_size + 
               z.brain_residual + z.clutch_size + z.habitat_generalism_scaled + 
               z.diet_breadth  + functional_diet,
             data=dat, na.action="na.fail", weights=weights)

summary(lm.mod)

# reasonably high R2 value so that is promising
# make a quick plot of this
summary_df <- tidy(lm.mod) %>%
  mutate(lwr_95_confint=confint(lm.mod)[,1]) %>%
  mutate(upr_95_confint=confint(lm.mod)[,2]) %>%
  mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
  mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
  mutate(model_type="global_model") %>%
  arrange(estimate)

ggplot(summary_df, aes(x=fct_inorder(term), y=estimate))+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), 
                width=0.8, position=position_dodge(width=0.6))+
  geom_point(position=position_dodge(width=0.6), aes(color=significance))+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_hline(yintercept=0, color="red")+
  scale_color_brewer(palette="Set1")+
  xlab("")+
  ylab("Standardized parameter estimate")

# the fact range size is really positive is kind of unsurprising because
# it is possible that species with large ranges might
# just have a higher possibility of interacting with urban areas by chance
# so I'll try re-running the model without that term
# to see what happens
lm.mod.2 <- lm(response ~ migration_status + z.log_body_size + z.log_flock_size +
               z.brain_residual + z.clutch_size + z.habitat_generalism_scaled + 
               z.diet_breadth  + functional_diet,
             data=dat, na.action="na.fail", weights=weights)

summary(lm.mod.2)

# reasonably high R2 value so that is promising
# make a quick plot of this
summary_df <- tidy(lm.mod.2) %>%
  mutate(lwr_95_confint=confint(lm.mod.2)[,1]) %>%
  mutate(upr_95_confint=confint(lm.mod.2)[,2]) %>%
  mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
  mutate(trend=ifelse(.$estimate >0, "positive", "negative")) %>%
  mutate(model_type="global_model") %>%
  arrange(estimate)

ggplot(summary_df, aes(x=fct_inorder(term), y=estimate))+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), 
                width=0.8, position=position_dodge(width=0.6))+
  geom_point(position=position_dodge(width=0.6), aes(color=significance))+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_hline(yintercept=0, color="red")+
  scale_color_brewer(palette="Set1")+
  xlab("")+
  ylab("Standardized parameter estimate")

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
  mutate(model_type="model_averaging") %>%
  arrange(estimate)

ggplot(model_results, aes(x=fct_inorder(term), y=estimate))+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), 
                width=0.8, position=position_dodge(width=0.6))+
  geom_point(position=position_dodge(width=0.6), aes(color=significance))+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_hline(yintercept=0, color="red")+
  scale_color_brewer(palette="Set1")+
  xlab("")+
  ylab("Standardized parameter estimate")

# but we need to consider the potential of phylogenetic signal in intra-annual variability
# so first let's test the phylogenetic signal of this variable
# function to read one tree in
read_one_tree<-function(path, x=50){
  
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
subset_tree <- function(bird_tree) {
  
  non_usa_sp <- bird_tree$tip.label[!bird_tree$tip.label %in% analysis$TipLabel]
  
  usa_bird_tree <- drop.tip(bird_tree, non_usa_sp)
  
  return(usa_bird_tree)
}

usa_tree <- subset_tree(bird_tree)

# now do a phylogenetic signal analysis
phylo_dat <- dat %>%
  dplyr::select(COMMON_NAME, TipLabel, response) %>%
  distinct()

phylo_dat.2 <- phylo_dat %>%
  dplyr::select(response)

row.names(phylo_dat.2) <- phylo_dat$TipLabel #name rows so that it matches the tree

p4d <- phylo4d(usa_tree, phylo_dat.2) #create phylobase object

ps <- phyloSignal(p4d,reps = 9999) #run calculation, p values a bit unstable at 999 reps

stats <- ps$stat %>%
  rownames_to_column(var="term") %>%
  mutate(value="Statistic")

p_values <- ps$pvalue %>%
  rownames_to_column(var="term") %>%
  mutate(value="P value")

phylo_sig_results <- bind_rows(stats, p_values)

# plot the intra annual variability on a tree
ggtree(p4d, layout='circular', aes(color=response), 
       ladderize = FALSE, size=1)+
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue'))+
  geom_tiplab(aes(angle=angle), size=2.5)+
  theme(legend.position = c(.05, .85))

# looks like strong phylo signal in intra-annual variability, which probably isn't terrible surprising
# so let's now repeat the above models, but make them phylogenetic models

# now do a phylogenetic signal analysis
phylo_dat_3 <- dat %>%
  dplyr::select(TipLabel, response, weights, migration_status, z.log_body_size, z.log_flock_size,
                z.log_range_size, z.brain_residual, z.clutch_size, z.habitat_generalism_scaled,
                z.diet_breadth, functional_diet) %>%
  distinct() %>%
  column_to_rownames(var="TipLabel")

phy_mod_rescaled <- phylolm(response ~ z.log_body_size + z.log_flock_size + 
                              z.log_range_size + z.brain_residual +
                              z.clutch_size + z.habitat_generalism_scaled + 
                              z.diet_breadth + functional_diet + migration_status,
                            data=phylo_dat_3, phy=usa_tree, na.action="na.fail", weights=weights)

summary(phy_mod_rescaled)
