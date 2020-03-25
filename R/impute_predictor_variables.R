### This is an R script to impute missing data
### because we have a lot of different traits
### compiled from different datasets
### we don't have a complete dataset
### after reading about a lot of different imputation procedures
### and talking with Shinichi Nakagawa, it sounds like mice is probably the best bet
### Blomberg also says that phylogenetic imputation can often cause weird things to happen
### because there is likely a strong phylogenetic signal in body size
### then the 'phylogeny' is --- in a way --- already accounted for in the imputation
### and including the phylogentic tree could be somewhat redundant
### would be good to show that body size shows a strong phylogentic signal

# packages
library(mice)
library(dplyr)
library(missCompare)


### read in predictor variables
predictors <- readRDS("Data/predictor_variables.RDS")

### for now I'll just select a few example columns
### as eventually I'll finalize the total list of traits to use in models
### I can also go through and manually fix the categorical missing data
### so will only focus on imputing continuous variables
data_to_impute <- predictors %>%
  dplyr::select(brain_size_mm3, brain_residual, diet_breadth, adult_body_mass_g,
                mean_flock_size, total_range_km2, habitat_generalism_scaled, clutch_size) %>%
  ungroup()


# set seed
set.seed(10)




