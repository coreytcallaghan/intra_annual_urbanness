# This script interacts with the big query database
# to get all the data for a given species
# and saves it as an RDS
# then I will process the data from those RDSs in a separate analytical step


# packages
library(readr)
library(bigrquery)
library(dbplyr)
library(dplyr)
library(lubridate)

# create connectio with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa')

lights <- tbl(con, 'virrs_5000_qa')

get_data_for_a_species <- function(species_name) {
  
  # get data from bigquery
  dat <- ebird %>%
    dplyr::filter(COMMON_NAME == species_name) %>%
    dplyr::filter(STATE_CODE != "US-AK") %>%
    dplyr::filter(STATE_CODE != "US-HI") %>%
    dplyr::filter(OBSERVATION_DATE > "2010-01-01") %>%
    dplyr::filter(LONGITUDE < -27) %>%
    dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, OBSERVATION_COUNT,
                  LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE) %>%
    left_join(., lights, by="LOCALITY_ID") %>%
    dplyr::select(-"_geo") %>%
    dplyr::select(-system_index) %>%
    collect(n=Inf)
  
  # add month and day to each dataframe
  dat <- dat %>%
    mutate(DAY=yday(OBSERVATION_DATE)) %>%
    mutate(MONTH=month(OBSERVATION_DATE))
  
  file_species <- gsub(" ", "_", species_name)
  
  saveRDS(dat, file = paste0("Data/species_RDS/", file_species, "_ebird_May19.RDS"))
  
}

# get a list of species which was done in a separate R script
# titled "get_list_of_potential_us_species.R"
species_list <- read_csv("Data/list_of_potential_species.csv") %>%
  .$COMMON_NAME

# Now run the function to separate eBird data for each species into an RDS
lapply(species_list, function(x) {get_data_for_a_species(x)})


