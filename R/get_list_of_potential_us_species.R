## This is an R script to get a list of species
## which have at least > 250 observations in each month

# packages
library(readr)
library(bigrquery)
library(dbplyr)
library(dplyr)
library(lubridate)
library(tidyr)

# create connectio with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa')

# read in clements data
clements <- read_csv("Data/Clements-Checklist-v2018-August-2018.csv") %>%
  rename(COMMON_NAME = `English name`) %>%
  rename(SCIENTIFIC_NAME = `scientific name`) %>%
  dplyr::filter(category=="species") %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, order, family) %>%
  distinct()

# interact with bigquery db
# to get a potential initial list of species
# that may be possible to be included
# another level of filtering will occur down below
df <- ebird %>%
  dplyr::filter(Country == "United States") %>%
  dplyr::filter(STATE_CODE != "US-AK") %>%
  dplyr::filter(STATE_CODE != "US-HI") %>%
  dplyr::filter(OBSERVATION_DATE > "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, CATEGORY) %>%
  group_by(COMMON_NAME, CATEGORY) %>%
  summarise(N=n()) %>%
  collect(n=Inf)

ropi <- df %>%
  dplyr::filter(COMMON_NAME == "Rock Pigeon") %>%
  dplyr::filter(CATEGORY == "domestic")

df2 <- df %>%
  filter(CATEGORY %in% c("species", "issf", "form")) %>%
  bind_rows(ropi) %>%
  group_by(COMMON_NAME) %>%
  summarise(N=sum(N)) %>%
  dplyr::filter(N>100) %>%
  left_join(., clements) %>%
  dplyr::filter(family != "Stercorariidae (Skuas and Jaegers)") %>%
  dplyr::filter(family != "Alcidae (Auks, Murres, and Puffins)") %>%
  dplyr::filter(family != "Diomedeidae (Albatrosses)") %>%
  dplyr::filter(family != "Oceanitidae (Southern Storm-Petrels)") %>%
  dplyr::filter(family != "Hydrobatidae (Northern Storm-Petrels)") %>%
  dplyr::filter(family != "Procellariidae (Shearwaters and Petrels)") %>%
  dplyr::filter(family != "Fregatidae (Frigatebirds)") %>%
  dplyr::filter(family != "Sulidae (Boobies and Gannets)")

# now get this list of species
# and see how many observations each has
# by month
df3 <- ebird %>%
  dplyr::filter(STATE_CODE != "US-AK") %>%
  dplyr::filter(STATE_CODE != "US-HI") %>%
  dplyr::filter(OBSERVATION_DATE > "2010-01-01") %>%
  dplyr::filter(LONGITUDE < -27) %>%
  dplyr::filter(COMMON_NAME %in% local(unique(df2$COMMON_NAME))) %>%
  mutate(MONTH=month(OBSERVATION_DATE)) %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, MONTH) %>%
  group_by(COMMON_NAME, MONTH) %>%
  summarize(N=n()) %>%
  collect(n=Inf)

# create a dataframe
# to join with to fill in 'empty months'
species_months <- data.frame(COMMON_NAME=rep(unique(df2$COMMON_NAME), 12)) %>%
  arrange(COMMON_NAME) %>%
  mutate(MONTH=rep(c(1:12), length(unique(df2$COMMON_NAME))))

df4 <- df3 %>%
  right_join(., species_months) %>%
  replace_na(list(MONTH=0)) %>%
  mutate(greater_than_250=ifelse(N>=250, 1, 0)) %>%
  group_by(COMMON_NAME) %>%
  summarize(months_greater_than_250=sum(greater_than_250)) %>%
  dplyr::filter(months_greater_than_250 == 12) %>%
  dplyr::select(COMMON_NAME)

## save the file as a csv
write_csv(df4, "Data/list_of_potential_species.csv")

