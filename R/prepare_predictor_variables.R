## This is an R script to prepare predictor variables
## and create one clean dataframe with common name
## and predictor variables

# packages
library(readr)
library(dplyr)
library(tidyr)

# read in clements clean data
clements_clean <- read_csv("Data/clements_clean.csv") %>%
  rename(COMMON_NAME=ebird_COMMON_NAME)

# first get a list of species which are considered for the response variables
# and join this with the clements taxonomy to bring in the scientific name
# and the "TipLabel"
species <- readRDS("Data/response_variables.RDS") %>%
  dplyr::select(COMMON_NAME) %>%
  distinct() %>%
  left_join(., clements_clean) %>%
  group_by(COMMON_NAME) %>%
  slice(1)

# but some species don't have a tip label
# for various reasons
# will manually fix these here
taxonomic_fix <- species %>%
  dplyr::filter(!complete.cases(TipLabel)) %>%
  dplyr::select(-TipLabel) %>%
  ungroup() %>%
  mutate(TipLabel=c("Gallinula_chloropus", "Porphyrio_porphyrio",
                    "Buteo_nitidus", "Troglodytes_troglodytes", 
                    "Rallus_longirostris", "Amphispiza_belli", 
                    "Charadrius_alexandrinus", "Gallinago_gallinago", "Aphelocoma_californica"))

species <- species %>%
  dplyr::filter(complete.cases(TipLabel)) %>%
  bind_rows(taxonomic_fix)

# read in body size data
# and prepare it to join up with
# the species data
body_size <- read_csv("Data/body_size_data/cleaned_body_size_data.csv")

body_size_matched.1 <- species %>%
  dplyr::select(COMMON_NAME) %>%
  left_join(., body_size) %>%
  dplyr::select(COMMON_NAME, adult_body_mass_g) %>%
  dplyr::filter(complete.cases(adult_body_mass_g))

body_size_matched.2 <- species %>%
  dplyr::filter(!COMMON_NAME %in% body_size_matched.1$COMMON_NAME) %>%
  dplyr::select(TipLabel) %>%
  left_join(., body_size, by="TipLabel") %>%
  dplyr::select(COMMON_NAME.x, adult_body_mass_g) %>%
  dplyr::filter(complete.cases(adult_body_mass_g)) %>%
  rename(COMMON_NAME=COMMON_NAME.x)

body_size_final <- bind_rows(body_size_matched.1,
                             body_size_matched.2) %>%
  distinct()


# Now repeat the same general process
# but for brain size data
# and other associated data in the same brain size
# data file
brain_size <- read_csv("Data/brain_size_data/brain_size_and_other_data.csv") %>%
  dplyr::select(-body_mass_g)

brain_size_matched.1 <- species %>%
  ungroup() %>%
  dplyr::select(ebird_SCIENTIFIC_NAME) %>%
  mutate(SCIENTIFIC_NAME=gsub(" ", "_", ebird_SCIENTIFIC_NAME)) %>%
  dplyr::select(SCIENTIFIC_NAME) %>% 
  left_join(., brain_size) %>%
  dplyr::filter(complete.cases(brain_residual)) %>%
  dplyr::select(-Order, -Family) %>%
  rename(ebird_SCIENTIFIC_NAME=SCIENTIFIC_NAME)

brain_size_matched.2 <- species %>%
  ungroup() %>%
  mutate(SCIENTIFIC_NAME=gsub(" ", "_", ebird_SCIENTIFIC_NAME)) %>%
  dplyr::filter(!SCIENTIFIC_NAME %in% brain_size_matched.1$ebird_SCIENTIFIC_NAME) %>%
  dplyr::select(TipLabel, SCIENTIFIC_NAME) %>%
  rename(ebird_SCIENTIFIC_NAME=SCIENTIFIC_NAME) %>%
  rename(SCIENTIFIC_NAME=TipLabel) %>%
  left_join(., brain_size, by="SCIENTIFIC_NAME") %>%
  dplyr::filter(complete.cases(brain_residual)) %>%
  dplyr::select(-Order, -Family, -SCIENTIFIC_NAME)
  
brain_size_final <- bind_rows(brain_size_matched.1,
                              brain_size_matched.2) %>%
  mutate(ebird_SCIENTIFIC_NAME=gsub("_", " ", ebird_SCIENTIFIC_NAME)) %>%
  distinct()

# Now repeat the same general process
# but for functional group categorical data
# which will require a bit of extra work as well
functional_data <- read_csv("Data/functional_data/BirdFuncDat.csv") %>%
  dplyr::select(8,9,20,26:30) %>%
  rename(SCIENTIFIC_NAME=Scientific) %>%
  mutate(SCIENTIFIC_NAME=gsub(" ", "_", .$SCIENTIFIC_NAME)) %>%
  rename(functional_diet=`Diet-5Cat`) %>%
  rename(ground=`ForStrat-ground`) %>%
  rename(understory=`ForStrat-understory`) %>%
  rename(mid_high=`ForStrat-midhigh`) %>%
  rename(canopy=`ForStrat-canopy`) %>%
  rename(aerial=`ForStrat-aerial`)
  
functional_matched.1 <- species %>%
  ungroup() %>%
  dplyr::select(ebird_SCIENTIFIC_NAME) %>%
  mutate(SCIENTIFIC_NAME=gsub(" ", "_", ebird_SCIENTIFIC_NAME)) %>%
  dplyr::select(SCIENTIFIC_NAME) %>% 
  left_join(., functional_data) %>%
  dplyr::filter(complete.cases(functional_diet)) %>%
  dplyr::select(-English) %>%
  rename(ebird_SCIENTIFIC_NAME=SCIENTIFIC_NAME)

functional_matched.2 <- species %>%
  ungroup() %>%
  mutate(SCIENTIFIC_NAME=gsub(" ", "_", ebird_SCIENTIFIC_NAME)) %>%
  dplyr::filter(!SCIENTIFIC_NAME %in% functional_matched.1$ebird_SCIENTIFIC_NAME) %>%
  dplyr::select(TipLabel, SCIENTIFIC_NAME) %>%
  rename(ebird_SCIENTIFIC_NAME=SCIENTIFIC_NAME) %>%
  rename(SCIENTIFIC_NAME=TipLabel) %>%
  left_join(., functional_data, by="SCIENTIFIC_NAME") %>%
  dplyr::filter(complete.cases(functional_diet)) %>%
  dplyr::select(-English, -SCIENTIFIC_NAME)

functional_final <- bind_rows(functional_matched.1,
                              functional_matched.2) %>%
  mutate(ebird_SCIENTIFIC_NAME=gsub("_", " ", ebird_SCIENTIFIC_NAME))

# now read in a different migration status delimitation
# with also some habitat characteristics
# first summarize and clean up the data
migration_habitat <- read_csv("Data/migration_status/migration_status_cleaned.csv") %>%
  group_by(COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarize(habitat_generalism=length(unique(IUCN_HABITATS))) %>%
  left_join(., migration_habitat <- read_csv("Data/migration_status/migration_status_cleaned.csv") %>%
              group_by(COMMON_NAME, SCIENTIFIC_NAME) %>%
              slice(1) %>%
              dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, MIGRATORY_STATUS)) %>%
  ungroup()

# now repeat the process above
# to match with the species list
mig_habitat_matched.1 <- species %>%
  dplyr::select(COMMON_NAME) %>%
  left_join(., migration_habitat) %>%
  dplyr::filter(complete.cases(MIGRATORY_STATUS))

mig_habitat_matched.2 <- species %>%
  ungroup() %>%
  dplyr::filter(!COMMON_NAME %in% mig_habitat_matched.1$COMMON_NAME) %>%
  dplyr::select(ebird_SCIENTIFIC_NAME, COMMON_NAME) %>%
  rename(SCIENTIFIC_NAME=ebird_SCIENTIFIC_NAME) %>%
  left_join(., migration_habitat, by="SCIENTIFIC_NAME") %>%
  dplyr::select(2, 1, 4, 5) %>%
  dplyr::filter(complete.cases(MIGRATORY_STATUS)) %>%
  rename(COMMON_NAME=COMMON_NAME.x)

mig_habitat_final <- bind_rows(mig_habitat_matched.1,
                               mig_habitat_matched.2) %>%
  dplyr::select(-SCIENTIFIC_NAME) %>%
  ungroup() %>%
  mutate(habitat_generalism_scaled=round(scales::rescale(habitat_generalism), digits=2)) %>%
  rename(migratory_status=MIGRATORY_STATUS)


# now get clutch size and mating system
# first read in data and clean it up to what is necessary
# and rename some column names
clutch_mating <- read_delim("Data/fecundity_and_life_history/avian_ssd_jan07.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(Species_name, English_name, Clutch_size, Mating_System) %>%
  rename(SCIENTIFIC_NAME=Species_name) %>%
  rename(COMMON_NAME=English_name) %>%
  mutate(Clutch_size=gsub("-999", NA, .$Clutch_size)) %>%
  mutate(Mating_System=gsub("-999", NA, .$Mating_System)) %>%
  rename(clutch_size=Clutch_size) %>%
  rename(mating_system=Mating_System)

# now join with the species list as above for previous datasets
clutch_mating_matched.1 <- species %>%
  dplyr::select(COMMON_NAME) %>%
  left_join(., clutch_mating) %>%
  dplyr::filter(complete.cases(clutch_size))

clutch_mating_matched.2 <- species %>%
  ungroup() %>%
  dplyr::filter(!COMMON_NAME %in% clutch_mating_matched.1$COMMON_NAME) %>%
  dplyr::select(ebird_SCIENTIFIC_NAME, COMMON_NAME) %>%
  rename(SCIENTIFIC_NAME=ebird_SCIENTIFIC_NAME) %>%
  left_join(., clutch_mating, by="SCIENTIFIC_NAME") %>%
  dplyr::select(2, 1, 4, 5) %>%
  dplyr::filter(complete.cases(clutch_size)) %>%
  rename(COMMON_NAME=COMMON_NAME.x)

clutch_mating_final <- bind_rows(clutch_mating_matched.1,
                                 clutch_mating_matched.2) %>%
  dplyr::select(-SCIENTIFIC_NAME)

# get flock size
# which was summarized monthly for each species
# so this is really the only one that can actually change by month
# and is repeated
flock_size <- readRDS("Data/flock_size/flock_size_per_month_for_each_species.RDS")

# get necessary data from flock size dataframe
# will take a mean for each species
# and can always write this out separately for a different model
# for each month
# which can be done in the analysis script when running models
# but this is robably related to the mean anyway
flock_size_final <- species %>%
  left_join(., flock_size) %>%
  group_by(COMMON_NAME) %>%
  summarize(mean_flock_size=mean(median_abund))

# now play around with range sizes
# which come from Federico
# which come from BirdLife international
# first read in a sort of key for later joining with range sizes and eBird
birdlife_key <- read_csv("Data/IUCN_data/HBW-BirdLife_Checklist_Version_3.csv") %>%
  dplyr::select(`Common name`, `Scientific name`) %>%
  distinct() %>%
  rename(COMMON_NAME=`Common name`) %>%
  rename(SCIENTIFIC_NAME=`Scientific name`) %>%
  dplyr::filter(complete.cases(.))

# now read in range size data and join with 'birdlife' 'key'
ranges <- read_delim("Data/range_size/bird.range.season041219.csv", 
                     ";", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(SCIENTIFIC_NAME=Species) %>%
  left_join(., birdlife_key)

# now follow the similar joining process as above
range_size_matched.1 <- species %>%
  dplyr::select(COMMON_NAME) %>%
  left_join(., ranges) %>%
  dplyr::filter(complete.cases(range.size.km2))

range_size_matched.2 <- species %>%
  ungroup() %>%
  dplyr::filter(!COMMON_NAME %in% range_size_matched.1$COMMON_NAME) %>%
  dplyr::select(ebird_SCIENTIFIC_NAME, COMMON_NAME) %>%
  rename(SCIENTIFIC_NAME=ebird_SCIENTIFIC_NAME) %>%
  left_join(., ranges, by="SCIENTIFIC_NAME") %>%
  dplyr::filter(complete.cases(range.size.km2)) %>%
  dplyr::select(-COMMON_NAME.y) %>%
  rename(COMMON_NAME=COMMON_NAME.x) %>%
  dplyr::select(2, 1, 3:7)

range_size_final <- bind_rows(range_size_matched.1,
                              range_size_matched.2) %>%
  dplyr::select(-SCIENTIFIC_NAME) %>%
  rename(resident_range_km2=resident.km2) %>%
  rename(breeding_range_km2=breeding.km2) %>%
  rename(non_breeding_range_km2=`non-breeding.km2`) %>%
  rename(passage_range_km2=passage.km2) %>%
  rename(total_range_km2=range.size.km2)


# Now put all of the above together
# into a dataframe that is for the predictor variables
predictor_variables <- species %>%
  left_join(., brain_size_final) %>%
  left_join(., body_size_final) %>%
  left_join(., functional_final) %>%
  left_join(., flock_size_final) %>%
  left_join(., range_size_final) %>%
  left_join(., mig_habitat_final) %>%
  left_join(., clutch_mating_final)

# check the class of the variables
skimr::skim(predictor_variables)

# clutch size is being treated as a character for some reason
# so I'll rewrite that to numeric
predictor_variables <- predictor_variables %>%
  mutate(clutch_size=as.numeric(as.character(clutch_size)))

# write out predictor variables
saveRDS(predictor_variables, "Data/predictor_variables.RDS")










