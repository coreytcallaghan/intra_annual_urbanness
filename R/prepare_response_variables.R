## This scrip reads in the summaries
## exported from the 'make_ggridges_and_get_urbanness.R' script
## and then saves out a couple RDSs of 'response variables'
## I also do some manual filtering of the species included
## ensuring tha they largely occur in the United States

# packages
library(dplyr)
library(purrr)
library(readr)

# read in ABA codes
aba_codes <- read_csv("Data/ABA_codes/ABA_Checklist-8.0.5.csv") %>%
  dplyr::select(2, 3, 5) %>%
  rename(COMMON_NAME=X2) %>%
  rename(SCIENTIFIC_NAME=X3) %>%
  rename(ABA_Code=X5) %>%
  dplyr::filter(complete.cases(ABA_Code))

# now join species with ABA code
potential_species <- read_csv("Data/list_of_potential_species.csv") %>%
  rename(ebird_COMMON_NAME=COMMON_NAME) %>%
  left_join(., read_csv("Data/clements_clean.csv"), by="ebird_COMMON_NAME") %>%
  rename(COMMON_NAME=ebird_COMMON_NAME) %>%
  dplyr::select(COMMON_NAME, ebird_SCIENTIFIC_NAME) %>%
  distinct() %>%
  left_join(aba_codes, by="COMMON_NAME")

# Now I want to remove any species that have a code 3, 4, or 5 which means they are 
# relatively uncommon occurrences in the ABA area
# some species don't join because they aren't on the ABA checklist and should also be removed
# also, I want to keep Gray-headed Swamphen as this species is on the checklist
# but as a different species name
species_to_keep <- potential_species %>%
  dplyr::filter(ABA_Code %in% c(1, 2)) %>%
  bind_rows(., potential_species %>%
              dplyr::filter(COMMON_NAME=="Gray-headed Swamphen")) %>%
  .$COMMON_NAME


# read in monthly dat
setwd("Data/species_monthly_summaries")
monthly_dat <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS) %>%
  dplyr::filter(COMMON_NAME %in% species_to_keep)

length(unique(monthly_dat$COMMON_NAME))

setwd("..")
setwd("..")
saveRDS(monthly_dat, "Data/response_variables.RDS")
  
  
  
  
  
# Old stuff when I did it before filtering by ABA code
  dplyr::filter(COMMON_NAME != "Amazon Kingfisher") %>%
  dplyr::filter(COMMON_NAME != "American Flamingo") %>%
  dplyr::filter(COMMON_NAME != "Aplomado Falcon") %>%
  dplyr::filter(COMMON_NAME != "Bananaquit") %>%
  dplyr::filter(COMMON_NAME != "Black-faced Grassquit") %>%
  dplyr::filter(COMMON_NAME != "Black-necked Swan") %>%
  dplyr::filter(COMMON_NAME != "Black-vented Oriole") %>%
  dplyr::filter(COMMON_NAME != "Blue-and-yellow Macaw") %>%
  dplyr::filter(COMMON_NAME != "Blue-crowned Parakeet") %>%
  dplyr::filter(COMMON_NAME != "Blue-throated Hummingbird") %>%
  dplyr::filter(COMMON_NAME != "Brown Jay") %>%
  dplyr::filter(COMMON_NAME != "Chestnut-fronted Macaw") %>%
  dplyr::filter(COMMON_NAME != "Collared Plover") %>%
  dplyr::filter(COMMON_NAME != "Coscoroba Swan") %>%
  dplyr::filter(COMMON_NAME != "Crimson-fronted Parakeet") %>%
  dplyr::filter(COMMON_NAME != "Green Parakeet") %>%
  dplyr::filter(COMMON_NAME != "Groove-billed Ani") %>%
  dplyr::filter(COMMON_NAME != "Mexican Jay") %>%
  dplyr::filter(COMMON_NAME != "Mitred Parakeet") %>%
  dplyr::filter(COMMON_NAME != "Monk Parakeet") %>%
  dplyr::filter(COMMON_NAME != "Morelet's Seedeater") %>%
  dplyr::filter(COMMON_NAME != "Northern Jacana") %>%
  dplyr::filter(COMMON_NAME != "Orange-winged Parrot") %>%
  dplyr::filter(COMMON_NAME != "Red-crowned Parrot") %>%
  dplyr::filter(COMMON_NAME != "Red-lored Parrot") %>%
  dplyr::filter(COMMON_NAME != "Red-masked Parakeet") %>%
  dplyr::filter(COMMON_NAME != "Rufous-backed Robin") %>%
  dplyr::filter(COMMON_NAME != "Rufous-capped Warbler") %>%
  dplyr::filter(COMMON_NAME != "Streak-backed Oriole") %>%
  dplyr::filter(COMMON_NAME != "Tropical Mockingbird") %>%
  dplyr::filter(COMMON_NAME != "Tropical Parula") %>%
  dplyr::filter(COMMON_NAME != "Turquoise-fronted Parrot") %>%
  dplyr::filter(COMMON_NAME != "White-cheeked Pintail") %>%
  dplyr::filter(COMMON_NAME != "White-eared Hummingbird") %>%
  dplyr::filter(COMMON_NAME != "White-eyed Parakeet") %>%
  dplyr::filter(COMMON_NAME != "White-faced Whistling-Duck") %>%
  dplyr::filter(COMMON_NAME != "White-fronted Parrot") %>%
  dplyr::filter(COMMON_NAME != "White-throated Thrush") %>%
  dplyr::filter(COMMON_NAME != "White-winged Parakeet") %>%
  dplyr::filter(COMMON_NAME != "Yellow-chevroned Parakeet") %>%
  dplyr::filter(COMMON_NAME != "Yellow-crowned Parrot") %>%
  dplyr::filter(COMMON_NAME != "Zenaida Dove")


  

# if I want to look at daily summaries
# but won't look at that quite yet
# will think about this more
setwd("Data/species_daily_summaries")
daily_dat <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)