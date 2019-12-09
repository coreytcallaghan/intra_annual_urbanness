## This scrip reads in the summaries
## exported from the 'make_ggridges_and_get_urbanness.R' script
## and then saves out a couple RDSs of 'response variables'
## I also do some manual filtering of the species included
## ensuring tha they largely occur in the United States

# packages
library(dplyr)
library(purrr)

setwd("Data/species_monthly_summaries")
monthly_dat <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS) %>%
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

setwd("..")
setwd("..")
saveRDS(monthly_dat, "Data/response_variables.RDS")
  


setwd("Data/species_daily_summaries")
daily_dat <- list.files(pattern = ".RDS") %>%
  map_dfr(readRDS)