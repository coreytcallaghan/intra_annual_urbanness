# This is an R script to make Figure 1 for the paper
# which will be four example species with their ggridges plotted
# for each month


# packages
library(dplyr)
library(ggplot2)
library(ggridges)
library(scales)
library(tidyr)


whwo <- readRDS(paste0("Data/species_RDS/White-headed_Woodpecker_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

rubl <- readRDS(paste0("Data/species_RDS/Rusty_Blackbird_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

acfl <- readRDS(paste0("Data/species_RDS/Acadian_Flycatcher_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

alhu <- readRDS(paste0("Data/species_RDS/Allen's_Hummingbird_ebird_May19.RDS")) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE))

all_dat <- bind_rows(whwo, rubl, acfl, alhu)

processed_dat <- readRDS("Data/response_variables.RDS") %>%
  dplyr::filter(COMMON_NAME %in% c("White-headed Woodpecker", "Rusty Blackbird",
                                   "Acadian Flycatcher", "Allen's Hummingbird")) %>%
  group_by(COMMON_NAME) %>%
  summarize(N=sum(number_obs),
            SD=round(sd(mean_urbanness), digits=1)) %>%
  unite(title, COMMON_NAME, N, sep="; N=") %>%
  unite(title, title, SD, sep="; SD=") %>%
  mutate(COMMON_NAME=c("Acadian Flycatcher", "Allen's Hummingbird",
                        "Rusty Blackbird", "White-headed Woodpecker"))

plot_dat <- all_dat %>%
  left_join(., processed_dat, by="COMMON_NAME")


# make a figure for visualization
# showing intra-annual differences for
# every species
# here is a vignette for ggridges: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
ggplot(plot_dat, aes(x=avg_rad, y=MONTH, height=..density..))+
  geom_density_ridges(stat="density", color="black", fill="steelblue4")+
  scale_x_log10(labels=comma)+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  ylab("")+
  xlab(bquote('Average radiance ('* 'nW' ~cm^-2~sr^-1*')'))+
  scale_y_discrete(labels=c("Dec", "Nov", "Oct", "Sep", "Aug", 
                            "Jul", "Jun", "May", "Apr", "Mar", "Feb", "Jan"), 
                   limits=c("Dec", "Nov", "Oct", "Sep", "Aug", 
                            "Jul", "Jun", "May", "Apr", "Mar", "Feb", "Jan"))+
  facet_wrap(~title)


ggsave("Figures/example_ggridges_unedited.png", width=7, height=6.3, units="in")


