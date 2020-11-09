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
library(phytools)
library(RColorBrewer)
library(readr)

source("R/global_functions.R")

# read in response variables
response <- readRDS("Data/response_variables.RDS") 

# I redid analyses that will influence the number of species
# because I now only use data from 2014, so some species won't meet the 250
# threshold for all 12 months.
# get a list of these species
# and remove them from possible analyses
species_without_enough_dat <- response %>%
  dplyr::filter(number_obs<=250) %>%
  dplyr::select(COMMON_NAME) %>%
  distinct()

response <- response %>%
  dplyr::filter(!COMMON_NAME %in% species_without_enough_dat$COMMON_NAME) %>%
  group_by(COMMON_NAME) %>%
  summarize(intra_annual_variability=sd(mean_urbanness),
            range_of_variability=(max(mean_urbanness) - min(mean_urbanness)),
            sd_of_intra_annual_variability=mean(sd_urbanness)) 

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

# # First let's just look at migration as a function of the response
# # to summarize
# migration_analysis <- analysis %>%
#   dplyr::select(1:8, migration_status, migratory_status) %>%
#   dplyr::filter(complete.cases(migration_status))
# 
# # The interesting thing here is whether intra-annual variability is explained by migration status
# # which is a reasonable hypothesis as you would expect that migratory species 
# # are more likely to have higher variability
# # first plot this
# ggplot(migration_analysis, aes(x=factor(migration_status, levels=c("Resident", "Migrant")),
#                      y=intra_annual_variability, fill=migration_status))+
#   geom_violin(position=position_dodge()) +
#   geom_boxplot(width=0.1, color="black", position = position_dodge(width =0.9))+
#   coord_flip()+
#   scale_fill_brewer(palette="Set1")+
#   scale_y_log10()+
#   theme_bw()+
#   theme(axis.text=element_text(color="black"))+
#   xlab("")+
#   ylab("Intra-annual urbanness variability")+
#   guides(fill=FALSE)
# 
# ggsave("Figures/intra_annual_variability_vs_migration_status.png", width=6.5, height=4.5, units="in")
# 
# 
# # left with 245 species which have complete data that can be analyzed
# # let's start by looking at correlation and distribution of response variable
# hist(migration_analysis$intra_annual_variability)
# 
# # heavily skewed so will try log transform
# hist(log(migration_analysis$intra_annual_variability))
# 
# # make plot of this for supplementary material
# brewer.pal(name="Set1", n=1)
# 
# ggplot(migration_analysis, aes(x=intra_annual_variability))+
#   geom_histogram(color="black", fill="#E41A1C", bins=20)+
#   scale_x_log10()+
#   theme_bw()+
#   theme(axis.text=element_text(color="black"))+
#   xlab("Intra-annual variability of urban-tolerance")+
#   ylab("Number of species")
# 
# ggsave("Figures/intra_annual_variability_histogram.png", width=6.5, height=4.5, units='in')
# 
# # summary statistics for paper
# summary(migration_analysis$intra_annual_variability)
# mean(migration_analysis$intra_annual_variability)
# sd(migration_analysis$intra_annual_variability)
# 
# migration_analysis %>%
#   group_by(migration_status) %>%
#   summarize(mean=mean(intra_annual_variability),
#             sd=sd(intra_annual_variability))










# now they match! So we have a potential sample size of 486 species
# I might do some imputation for the missing species
# but in this R script lets just run the models for species which have data for each category
# that I am interested in
# so will first filter to the columns of interest and then
# filter out all species without complete data
analysis <- analysis %>%
  ungroup() %>%
  dplyr::select(1:8, brain_residual, diet_breadth, adult_body_mass_g, functional_diet,
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

# make plot of this for supplementary material
brewer.pal(name="Set1", n=1)

ggplot(analysis, aes(x=intra_annual_variability))+
  geom_histogram(color="black", fill="#E41A1C", bins=20)+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Intra-annual variability of urban tolerance")+
  ylab("Number of species")

ggsave("Figures/intra_annual_variability_histogram.png", width=6.5, height=4.5, units='in')

# summary statistics for paper
summary(analysis$intra_annual_variability)
mean(analysis$intra_annual_variability)
sd(analysis$intra_annual_variability)

analysis %>%
  group_by(migration_status) %>%
  summarize(mean=mean(intra_annual_variability),
            sd=sd(intra_annual_variability))

# plot the relationship between sd and range
ggplot(analysis, aes(x=intra_annual_variability, y=range_of_variability))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

# plot ggpairs of variables
analysis %>%
  dplyr::select(9:17) %>%
  ggpairs()

# doesn't look like too much colinearity problems
# some of the variables look pretty skewed, unsurprisingly
# so will need to log-transform these and they should look pretty normal
analysis <- analysis %>%
  mutate(response=log10(intra_annual_variability)) %>%
  mutate(log_body_size=log(adult_body_mass_g)) %>%
  mutate(log_flock_size=log(mean_flock_size)) %>%
  mutate(log_range_size=log(total_range_km2)) %>%
  mutate(weights=1/sd_of_intra_annual_variability)

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


##########################################
# prepare table S1
table_s1 <- analysis %>%
  dplyr::select(1, 5:8, 2, 9:17) %>%
  left_join(., readRDS("Data/response_variables.RDS") %>%
              dplyr::filter(COMMON_NAME %in% analysis$COMMON_NAME) %>%
              #dplyr::filter(!COMMON_NAME %in% c("Monk Parakeet", "Red-crowned Parrot")) %>%
              dplyr::select(COMMON_NAME, MONTH, mean_urbanness) %>%
              group_by(COMMON_NAME) %>%
              pivot_wider(names_from=MONTH, values_from=mean_urbanness)) %>%
  dplyr::select(1:5, 17:27, 6:16) %>%
  rename(SCIENTIFIC_NAME=ebird_SCIENTIFIC_NAME)

write_csv(table_s1, "Results/table_s1.csv")



# but we need to consider the potential of phylogenetic signal in intra-annual variability
# so first let's test the phylogenetic signal of this variable
# function to read one tree in
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

# now do a phylogenetic signal analysis
phylo_dat <- analysis %>%
  dplyr::select(COMMON_NAME, TipLabel, response) %>%
  distinct()

phylo_dat.2 <- phylo_dat %>%
  dplyr::select(response)

row.names(phylo_dat.2) <- phylo_dat$TipLabel #name rows so that it matches the tree

p4d <- phylo4d(con_tree, phylo_dat.2) #create phylobase object

ps <- phyloSignal(p4d,reps = 9999) #run calculation, p values a bit unstable at 999 reps

stats <- ps$stat %>%
  rownames_to_column(var="term") %>%
  mutate(value="Statistic")

p_values <- ps$pvalue %>%
  rownames_to_column(var="term") %>%
  mutate(value="P value")

phylo_sig_results <- bind_rows(stats, p_values) %>%
  round_df(., digits=4) %>%
  knitr::kable()

phylo_sig_results

# plot the intra annual variability on a tree
plasma_pal <- c(viridis::plasma(n = 12, direction=-1))

ggtree(p4d, layout='circular', aes(color=response), 
       ladderize = FALSE, size=1)+
  #scale_color_gradient(low = "yellow", high = "red", na.value = NA)+
  scale_color_gradientn(colors=plasma_pal, 
                        name = "Intra-annual variability:  ", breaks=c(-4.5,-2, 0), 
                        labels = c(0.000, 0.01, 1))+
  #scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue'))+
  geom_tiplab(aes(angle=angle), size=2.5)+
  theme(legend.position = "bottom")
  #guides(color=guide_legend(label.hjust = 0.01))

ggsave("Figures/phylo_tree_of_intra_annual_variability.png", width=10, height=10, units="in")

# looks like strong phylo signal in intra-annual variability, which probably isn't terrible surprising

# The interesting thing here is whether intra-annual variability is explained by migration status
# which is a reasonable hypothesis as you would expect that migratory species 
# are more likely to have higher variability
# first plot this
ggplot(analysis, aes(x=factor(migration_status, levels=c("Resident", "Migrant")),
                               y=intra_annual_variability, fill=migration_status))+
  geom_violin(position=position_dodge()) +
  geom_boxplot(width=0.1, color="black", position = position_dodge(width =0.9))+
  coord_flip()+
  scale_fill_brewer(palette="Set1")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("")+
  ylab("Intra-annual urban tolerance variability")+
  guides(fill=FALSE)

ggsave("Figures/intra_annual_variability_vs_migration_status.png", width=6.5, height=4.5, units="in")


dat <- analysis %>%
  mutate(z.log_body_size=rescale(log_body_size)) %>%
  mutate(z.log_flock_size=rescale(log_flock_size)) %>%
  mutate(z.log_range_size=rescale(log_range_size)) %>%
  mutate(z.brain_residual=rescale(brain_residual)) %>%
  mutate(z.clutch_size=rescale(clutch_size)) %>%
  mutate(z.habitat_generalism_scaled=rescale(habitat_generalism_scaled)) %>%
  mutate(z.diet_breadth=rescale(diet_breadth))


phylo_dat_3 <- dat %>%
  dplyr::select(TipLabel, response, weights, migration_status, z.log_body_size, z.log_flock_size,
                z.log_range_size, z.brain_residual, z.clutch_size, z.habitat_generalism_scaled,
                z.diet_breadth, functional_diet) %>%
  distinct() %>%
  column_to_rownames(var="TipLabel")


# Now run a phylo model
phy_mod_migration <- phylolm(response ~ migration_status,
                             data=phylo_dat_3, phy=con_tree, na.action="na.fail", weights=weights)

summary(phy_mod_migration)

# now run a non-phylogenetic model
simple_mig_mod <- lm(response ~ migration_status, data=dat, na.action="na.fail", weights=weights)
summary(simple_mig_mod)














############################## OLD STUFF!!! #############################################
######################## NO LONGER IN THE PAPER

phylo_dat_3 <- dat %>%
  dplyr::select(TipLabel, response, weights, migration_status, z.log_body_size, z.log_flock_size,
                z.log_range_size, z.brain_residual, z.clutch_size, z.habitat_generalism_scaled,
                z.diet_breadth, functional_diet) %>%
  distinct() %>%
  column_to_rownames(var="TipLabel")

phy_mod_rescaled <- phylolm(response ~ z.log_body_size + z.log_flock_size + z.brain_residual +
                              z.clutch_size + z.habitat_generalism_scaled + 
                              z.diet_breadth + migration_status,
                            data=phylo_dat_3, phy=con_tree, na.action="na.fail", weights=weights)

summary(phy_mod_rescaled)

phy_mod_results <- as.data.frame(summary(phy_mod_rescaled)$coefficients) %>%
  cbind(confint(phy_mod_rescaled)) %>%
  rownames_to_column(var="term") %>%
  rename(estimate=Estimate) %>%
  rename(std.error=StdErr) %>%
  rename(lwr_95_confint=`2.5 %`) %>%
  rename(upr_95_confint=`97.5 %`) %>%
  mutate(significance=ifelse(p.value <=0.05, "Significant", "Non-significant")) %>%
  mutate(trend=ifelse(.$estimate >0, "positive", "negative"))


phy_mod_results %>%
  dplyr::filter(term != "(Intercept)") %>%
  mutate(term=case_when(
    term == "z.log_body_size" ~ "Body size (log)",
    term == "z.log_flock_size" ~ "Mean flock size (log)",
    term == "z.log_range_size" ~ "Range size (log km2)",
    term == "z.brain_residual" ~ "Brain residual",
    term == "z.clutch_size" ~ "Clutch size",
    term == "z.habitat_generalism_scaled" ~ "Habitat generalism",
    term == "z.diet_breadth" ~ "Diet breadth",
    term == "migration_statusResident" ~ "Resident")) %>%
  ggplot(., aes(x=fct_inorder(term), y=estimate))+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), 
                width=0.8, position=position_dodge(width=0.6))+
  geom_point(position=position_dodge(width=0.6), aes(color=significance))+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_hline(yintercept=0, color="red")+
  scale_color_brewer(palette="Set1")+
  xlab("")+
  ylab("Parameter estimate")+
  labs(color="Significance")

ggsave("Figures/phylogenetic_global_model_of_intra_annual_variability.png", width=7.5, height=6.8, units="in")



# Let's try a k-means clustering approach to see if we can cluster species
# based on their seasonal responses
# i.e., those high in breeding season low otherwise, those high during migraton, etc.
all_response <- readRDS("Data/response_variables.RDS") %>%
  dplyr::filter(COMMON_NAME %in% analysis$COMMON_NAME) %>%
  #dplyr::filter(!COMMON_NAME %in% c("Monk Parakeet", "Red-crowned Parrot")) %>%
  dplyr::select(COMMON_NAME, MONTH, mean_urbanness) %>%
  group_by(COMMON_NAME) %>%
  mutate(mean_urbanness=scales::rescale(mean_urbanness)) %>%
  pivot_wider(names_from=MONTH, values_from=mean_urbanness) %>%
  column_to_rownames(var="COMMON_NAME")

library(gridExtra)
library(cluster)
library(factoextra)

distance <- get_dist(all_response)
fviz_dist(distance, gradient=list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
fviz_nbclust(all_response, kmeans, method="wss")
fviz_nbclust(all_response, kmeans, method="silhouette")

k2 <- kmeans(all_response, centers=2, nstart=30)

fviz_cluster(k2, data=all_response)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

# now get out the clusters and plot their seasonal trends to see if they make sense...
clusters <- data.frame(cluster=k2$cluster) %>%
  rownames_to_column(var="COMMON_NAME")

readRDS("Data/response_variables.RDS") %>%
  dplyr::filter(COMMON_NAME %in% analysis$COMMON_NAME) %>%
  #dplyr::filter(!COMMON_NAME %in% c("Monk Parakeet", "Red-crowned Parrot")) %>%
  dplyr::select(COMMON_NAME, MONTH, mean_urbanness) %>%
  group_by(COMMON_NAME) %>%
  mutate(mean_urbanness=scales::rescale(mean_urbanness)) %>%
  left_join(clusters) %>%
  mutate(cluster_name="Cluster") %>%
  unite(cluster, cluster_name, cluster, sep=" ") %>%
  mutate(cluster=as.character(cluster)) %>%
  ggplot()+
  geom_point(aes(x=MONTH, y=mean_urbanness, group=COMMON_NAME), size=0.6, color="gray90")+
  geom_line(aes(x=MONTH, y=mean_urbanness, group=COMMON_NAME), size=0.6, color="gray90")+
  #scale_y_continuous(labels = scales::percent)+
  geom_smooth(aes(x=MONTH, y=mean_urbanness, group=cluster, color=cluster), 
              method="gam", formula = y ~ s(x, bs = "cc"), size=1.2)+
  facet_wrap(~cluster, scales="free_y")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Month")+
  ylab("Scaled urbanness")+
  theme(legend.position="none")+
  scale_color_brewer(palette="Set1")+
  scale_x_discrete(breaks=c("Jan", "Feb", "Mar", "Apr", 
                            "May", "Jun", "Jul", "Aug",
                            "Sep", "Oct", "Nov", "Dec"), labels=c("Jan", "", "", "", "", "Jun", 
                                              "", "", "", "", "", "Dec"))

ggsave("Figures/clustering_results.png", width=7.5, height=4.8, units="in")



clusters %>%
  left_join(analysis) %>%
  group_by(cluster, migration_status) %>%
  summarize(N=n())


##########################################
# prepare table S1
table_s1 <- analysis %>%
  left_join(clusters) %>%
  dplyr::select(1, 5:8, 2, 9:17, 24) %>%
  left_join(., readRDS("Data/response_variables.RDS") %>%
              dplyr::filter(COMMON_NAME %in% analysis$COMMON_NAME) %>%
              #dplyr::filter(!COMMON_NAME %in% c("Monk Parakeet", "Red-crowned Parrot")) %>%
              dplyr::select(COMMON_NAME, MONTH, mean_urbanness) %>%
              group_by(COMMON_NAME) %>%
              pivot_wider(names_from=MONTH, values_from=mean_urbanness)) %>%
  dplyr::select(1:5, 17:28, 6:16) %>%
  rename(SCIENTIFIC_NAME=ebird_SCIENTIFIC_NAME)

write_csv(table_s1, "Results/table_s1.csv")

# fit a linear global model
# standardize the model
simple_mig_mod <- lm(response ~ migration_status, data=dat, na.action="na.fail", weights=weights)
summary(simple_mig_mod)
AIC(simple_mig_mod)







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


