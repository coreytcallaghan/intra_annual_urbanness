# This is an R script to read in
# a dataset for each species
# and make a ggridges plot by month
# and export that ggridges plot to a folder
# and at the same time make a small dataframe of
# species-urban scores by month

# here is a vignette for ggridges: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html

# packages
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggridges)


# function to apply for a given species

make_ggridges_and_extract_urban_scores <- function(species_name) {}


file_names <- list.files("city_bird_urbanness_paper/Data/species_RDS/")

for (i in file_names) {

df <- readRDS(paste0("city_bird_urbanness_paper/Data/species_RDS/", i)) %>%
  mutate(month=month(OBSERVATION_DATE, label=TRUE, abbr=TRUE))


title_name <- unique(df$COMMON_NAME)


ggplot(df, aes(x=avg_rad_mean_5k, y=month, height=..density..))+
  geom_density_ridges(stat="density", fill="skyblue2")+
  scale_x_log10()+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  ggtitle(paste0(title_name))+
  ylab("")+
  xlab(bquote('Average radiance ('* 'nW' ~cm^-2~sr^-1*')'))+
  scale_y_discrete(labels=c("Dec", "Nov", "Oct", "Sep", "Aug", 
                            "Jul", "Jun", "May", "Apr", "Mar", "Feb", "Jan"), 
                   limits=c("Dec", "Nov", "Oct", "Sep", "Aug", 
                            "Jul", "Jun", "May", "Apr", "Mar", "Feb", "Jan"))

ggsave(filename = paste0("city_bird_urbanness_paper/Data/species_ggridges/", title_name, ".png"),
       width=4.6, height=3.8, units="in")

urbanness <- df %>%
  group_by(COMMON_NAME, month) %>%
  summarise(urban_score=median(avg_rad_mean_5k, na.rm=TRUE))

saveRDS(urbanness, file = paste0("city_bird_urbanness_paper/Data/species_monthly_summaries/", title_name, ".RDS"))

}


