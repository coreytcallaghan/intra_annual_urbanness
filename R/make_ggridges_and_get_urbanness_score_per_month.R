# This is an R script to read in
# a dataset for each species
# and make a ggridges plot by month
# and export that ggridges plot to a folder
# and at the same time make a small dataframe of
# species-urban scores by month
# and a species-urban scores by day
# the RDS files were too large to push all of them
# so I pushed 30 'example RDSs'
# which would allow someone to reproduce this code for the sake of completeness
# to know what is happening
# if you wanted to reproduce this, below where it says "species_RDS", you would replace with "species_RDS_examples"
# which would loop through the 30 example species
# but, as noted this code currently is on a list of files that are not pushed to the repository

# packages
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggridges)
library(mgcv)
library(scales)


# get list of file names
file_names <- list.files("Data/species_RDS/")

# loop through the filenames
# and do some summaries
for (i in file_names) {

df <- readRDS(paste0("Data/species_RDS/", i)) %>%
  mutate(MONTH=month(OBSERVATION_DATE, abbr=TRUE, label=TRUE)) %>%
  dplyr::filter(OBSERVATION_DATE > "2014-01-01")


title_name <- unique(df$COMMON_NAME)

# make a figure for visualization
# showing potential intra-annual differences for
# every species
# here is a vignette for ggridges: https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
ggplot(df, aes(x=avg_rad, y=MONTH, height=..density..))+
  geom_density_ridges(stat="density", fill="skyblue2")+
  scale_x_log10(labels=comma)+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  ggtitle(paste0(title_name))+
  ylab("")+
  xlab(bquote('Average radiance ('* 'nW' ~cm^-2~sr^-1*')'))+
  scale_y_discrete(labels=c("Dec", "Nov", "Oct", "Sep", "Aug", 
                            "Jul", "Jun", "May", "Apr", "Mar", "Feb", "Jan"), 
                   limits=c("Dec", "Nov", "Oct", "Sep", "Aug", 
                            "Jul", "Jun", "May", "Apr", "Mar", "Feb", "Jan"))

# save figure out
ggsave(filename = paste0("Figures/species_ggridges/", title_name, ".png"),
       width=4.6, height=3.8, units="in")

# write a function to resample the urbaness for each month
# using 100 random samples
# which helps to account for a lot of variation
resample_urbanness_month_function <- function(draw){
  
  dat <- df %>%
    group_by(MONTH) %>%
    sample_n(100) %>%
    group_by(COMMON_NAME, MONTH) %>%
    summarise(urban_score=median(avg_rad, na.rm=TRUE)) %>%
    mutate(month_numeric=1:12) %>%
    mutate(resample=draw)
  
}

# repeat the above function 1000 times
# to converge on a 'mean' urbanness score for each month
# will also be able to see which month have the highest variance around it
list_of_resamples <- lapply(c(1:1000), function(x){resample_urbanness_month_function(x)})
resample_month_dfs <- bind_rows(list_of_resamples)

# summarize the resampled dataframes into a single
# mean urbanness score for each month
resampled_urbanness_month <- resample_month_dfs %>%
  group_by(MONTH) %>%
  summarize(mean_urbanness=mean(urban_score),
            sd_urbanness=sd(urban_score))

# now get the static urbanness month score (without resampling)
# and join with the resampled urbanness month score above
urbanness_month <- df %>%
  group_by(COMMON_NAME, MONTH) %>%
  summarise(urban_score=median(avg_rad, na.rm=TRUE),
            number_obs=n()) %>%
  mutate(month_numeric=1:12) %>%
  left_join(., resampled_urbanness_month, by="MONTH")

# save the as an RDS
saveRDS(urbanness_month, file = paste0("Data/species_monthly_summaries/", title_name, ".RDS"))

}



############################################################################
# repeat the above process, but for 'day of year' as opposed to month
# write a function to resample the urbaness for each day
# using 10 random samples per day
# as it is possible that a lot of days of the year could have very few samples
# for some species
# will also add 'replace=TRUE' in case a species has <10 observations on any given day
# this should be fine and would likely only influence the species with lower observations anyway, which would
# likely be the ones with the least variance anyhow
resample_urbanness_day_function <- function(draw){
  
  dat <- df %>%
    group_by(DAY) %>%
    sample_n(10, replace=TRUE) %>%
    group_by(COMMON_NAME, DAY) %>%
    summarise(urban_score=median(avg_rad, na.rm=TRUE)) %>%
    mutate(resample=draw)
  
}

# repeat the above function 1000 times
# to converge on a 'mean' urbanness score for each day
list_of_resamples <- lapply(c(1:1000), function(x){resample_urbanness_day_function(x)})
resample_day_dfs <- bind_rows(list_of_resamples)

# summarize the resampled dataframes into a single
# mean urbanness score for each month
resampled_urbanness_day <- resample_day_dfs %>%
  group_by(DAY) %>%
  summarize(mean_urbanness=mean(urban_score),
            sd_urbanness=sd(urban_score))

urbanness_day <- df %>%
  group_by(COMMON_NAME, DAY) %>%
  summarise(urban_score=median(avg_rad, na.rm=TRUE),
            number_obs=n()) %>%
  left_join(., resampled_urbanness_day, by="DAY")

# save out this RDS
saveRDS(urbanness_day, file = paste0("Data/species_daily_summaries/", title_name, ".RDS"))

}


