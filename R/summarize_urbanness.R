## a script to explore and summarize urbanness of the species included in the analysis
## Still a preliminary script and will probably add to the 'analysis' script of some sort
## once everything is figured out

# packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(GGally)

# read in dat
dat <- readRDS("Data/response_variables.RDS")

# number of species possibly included
length(unique(dat$COMMON_NAME))

# plot the relationship between
# the total urbanness
# and the resampled urbanness
ggplot(dat, aes(x=urban_score, y=mean_urbanness))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")+
  ylab("Resampled urbanness")+
  xlab("Total urbanness")

# plot the relationship between each month's urbanness
# for each month, using ggpairs
# to see if they are correlated
jan <- dat %>%
  dplyr::filter(MONTH=="Jan") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Jan=mean_urbanness)
feb <- dat %>%
  dplyr::filter(MONTH=="Feb") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Feb=mean_urbanness)
mar <- dat %>%
  dplyr::filter(MONTH=="Mar") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Mar=mean_urbanness)
apr <- dat %>%
  dplyr::filter(MONTH=="Apr") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Apr=mean_urbanness)
may <- dat %>%
  dplyr::filter(MONTH=="May") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(May=mean_urbanness)
jun <- dat %>%
  dplyr::filter(MONTH=="Jun") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Jun=mean_urbanness)
jul <- dat %>%
  dplyr::filter(MONTH=="Jul") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Jul=mean_urbanness)
aug <- dat %>%
  dplyr::filter(MONTH=="Aug") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Aug=mean_urbanness)
sep <- dat %>%
  dplyr::filter(MONTH=="Sep") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Sep=mean_urbanness)
oct <- dat %>%
  dplyr::filter(MONTH=="Oct") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Oct=mean_urbanness)
nov <- dat %>%
  dplyr::filter(MONTH=="Nov") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Nov=mean_urbanness)
dec <- dat %>%
  dplyr::filter(MONTH=="Dec") %>%
  dplyr::select(COMMON_NAME, mean_urbanness) %>%
  rename(Dec=mean_urbanness)

jan %>%
  left_join(., feb) %>%
  left_join(., mar) %>%
  left_join(., apr) %>%
  left_join(., may) %>%
  left_join(., jun) %>%
  left_join(., jul) %>%
  left_join(., aug) %>%
  left_join(., sep) %>%
  left_join(., oct) %>%
  left_join(., nov) %>%
  left_join(., dec) %>%
  ggpairs(2:13)

# summarize intra-annual variability
summary <- dat %>%
  group_by(COMMON_NAME) %>%
  summarize(intra_annual_variance=sd(mean_urbanness),
            total_urbanness=mean(mean_urbanness))

ggplot(summary, aes(x=total_urbanness))+
  geom_histogram(bins=45, color="black", fill="orange")+
  scale_x_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Count")+
  xlab("Mean urban score (log-scale)")


ggplot(summary, aes(x=intra_annual_variance))+
  geom_histogram(bins=45, color="black", fill="orange")+
  #scale_x_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ylab("Count")+
  xlab("Intra-annual urbanness variability")


# plot total versus sd urbanness
ggplot(summary, aes(x=total_urbanness, y=intra_annual_variance))+
  geom_point()+
  theme_bw()+
  scale_x_log10(labels=comma)+
  theme(axis.text=element_text(color="black"))+
  #geom_smooth(method="lm")+
  ylab("Intra-annual urbanness variability")+
  xlab("Mean urban score (log-scale)")







