##

# packages
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(tidyr)
library(ggrepel)


# read in non phylogenetic modelling results
load("Results/complete_species_non_phylo.RData")


# read in phylogenetic modelling results
load("Results/complete_species_phylo.RData")


global_model_results <- bind_rows(non_phylo_global_model_results, 
                                  phylo_global_model_results)


model_averaging_results <- bind_rows(non_phylo_global_model_results,
                                     phylo_global_model_results %>%
                                       mutate(model_type="phylo_model_averaging"))


# make a plot showing the correlation between
# global model estimates and model averaged parameter estimates
# first for non-phylogenetic analysis
non_phylo_global_model_results %>%
  dplyr::select(term, estimate, lwr_95_confint,
                upr_95_confint, model_type, MONTH) %>%
  bind_rows(., non_phylo_model_averaging_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH)) %>%
  dplyr::filter(term != "(Intercept)") %>%
  pivot_wider(names_from=model_type, values_from=c(estimate, lwr_95_confint, upr_95_confint)) %>%
  ggplot(., aes(x=estimate_global_model, y=estimate_model_averaging))+
  geom_point(color="limegreen")+
  geom_smooth(method="lm", color="orchid3")+
  facet_wrap(~MONTH, scales="free")+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  xlab("Global model estimate")+
  ylab("Model averaing estimate")+
  ggtitle("Non-phylogenetic analyses")


# now for phylogenetic analyses
phylo_global_model_results %>%
  dplyr::select(term, estimate, lwr_95_confint,
                upr_95_confint, model_type, MONTH) %>%
  bind_rows(., phylo_model_averaging_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH)) %>%
  dplyr::filter(term != "(Intercept)") %>%
  pivot_wider(names_from=model_type, values_from=c(estimate, lwr_95_confint, upr_95_confint)) %>%
  ggplot(., aes(x=estimate_phylo_global_mod, y=estimate_model_averaging))+
  geom_point(color="limegreen")+
  geom_smooth(method="lm", color="orchid3")+
  facet_wrap(~MONTH, scales="free")+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  xlab("Global model estimate")+
  ylab("Model averaing estimate")+
  ggtitle("Phylogenetic analyses")


# make a plot to look at difference between parameter estimates for phylo
# and non-phylogenetic models
non_phylo_global_model_results %>%
  dplyr::select(term, estimate, lwr_95_confint,
                upr_95_confint, model_type, MONTH) %>%
  bind_rows(., phylo_global_model_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH)) %>%
  dplyr::filter(term != "(Intercept)") %>%
  mutate(term=gsub("rescale\\(", "z.", .$term)) %>%
  mutate(term=gsub("\\)", "", .$term)) %>%
  pivot_wider(names_from=model_type, values_from=c(estimate, lwr_95_confint, upr_95_confint)) %>%
  ggplot(., aes(x=estimate_global_model, y=estimate_phylo_global_mod))+
  geom_point(color="limegreen")+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  geom_abline(slope=1,intercept=0, color="orchid3")+
  xlab("Non-phylogenetic model")+
  ylab("Phylogenetic model")+
  geom_text_repel(aes(label = term), 
                  box.padding = unit(0.45, "lines"))+
  facet_wrap(~MONTH, scales="free")


# let's plot the parameter estimates through time
# for each of the methods
phylo_global_model_results %>%
  dplyr::filter(term != "(Intercept)") %>%
  mutate(term=gsub("rescale\\(", "z.", .$term)) %>%
  mutate(term=gsub("\\)", "", .$term)) %>%
  arrange(term, MONTH) %>%
  mutate(month_num=rep(1:12, 7)) %>%
  ggplot(., aes(x=month_num, y=estimate))+
  geom_point(color="limegreen")+
  facet_wrap(~term, scales="free")+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  ylab("Parameter estimate")+
  xlab("")+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=11))+
  scale_x_continuous(breaks=c(1:12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))






