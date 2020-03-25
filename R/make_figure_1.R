# This is an R script to make Figure 1 for the paper
# which will be four example species with their ggridges plotted
# for each month


# packages
library(dplyr)
library(ggplot2)
library(ggridges)
library(scales)
library(tidyr)
library(lubridate)


whwo <- readRDS(paste0("Data/species_RDS/White-headed_Woodpecker_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

rubl <- readRDS(paste0("Data/species_RDS/Rusty_Blackbird_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

acfl <- readRDS(paste0("Data/species_RDS/Acadian_Flycatcher_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

alhu <- readRDS(paste0("Data/species_RDS/Allen's_Hummingbird_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

caja <- readRDS(paste0("Data/species_RDS/Canada_Jay_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

haha <- readRDS(paste0("Data/species_RDS/Harris's_Hawk_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

weta <- readRDS(paste0("Data/species_RDS/Western_Tanager_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

ospr <- readRDS(paste0("Data/species_RDS/Osprey_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

oven <- readRDS(paste0("Data/species_RDS/Ovenbird_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

ambi <- readRDS(paste0("Data/species_RDS/American_Bittern_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

hosp <- readRDS(paste0("Data/species_RDS/House_Sparrow_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

all_dat <- bind_rows(haha, caja, ambi, oven, hosp, weta)

processed_dat <- readRDS("Data/response_variables.RDS") %>%
  dplyr::filter(COMMON_NAME %in% c("Harris's Hawk", "Canada Jay",
                                   "American Bittern", "Ovenbird",
                                   "House Sparrow", "Western Tanager")) %>%
  group_by(COMMON_NAME) %>%
  summarize(N=sum(number_obs),
            SD=round(sd(mean_urbanness), digits=2)) %>%
  unite(title, COMMON_NAME, N, sep="; N=") %>%
  unite(title, title, SD, sep="; SD=") %>%
  mutate(COMMON_NAME=c("American Bittern", "Canada Jay",
                       "Harris's Hawk", "House Sparrow", 
                       "Ovenbird", "Western Tanager"))

plot_dat <- all_dat %>%
  left_join(., processed_dat, by="COMMON_NAME")


# make a figure for visualization
# showing intra-annual differences for
# every species
# here is a vignette for ggridges: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
ggplot(plot_dat, aes(x=avg_rad, y=MONTH, fill=title, height=..density..))+
  geom_density_ridges(stat="density", color="black")+
  scale_x_log10(labels=comma)+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("")+
  xlab(bquote('Average radiance ('* 'nW' ~cm^-2~sr^-1*')'))+
  scale_y_discrete(labels=c("Dec", "Nov", "Oct", "Sep", "Aug", 
                            "Jul", "Jun", "May", "Apr", "Mar", "Feb", "Jan"), 
                   limits=c("Dec", "Nov", "Oct", "Sep", "Aug", 
                            "Jul", "Jun", "May", "Apr", "Mar", "Feb", "Jan"))+
  facet_wrap(~title, scales="free")+
  guides(fill=FALSE)


ggsave("Figures/example_ggridges_unedited.png", width=8.5, height=7.3, units="in")


