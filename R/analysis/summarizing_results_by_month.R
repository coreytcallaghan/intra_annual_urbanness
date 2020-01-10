##

# packages
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(tidyr)
library(ggrepel)
library(wesanderson)


# read in non phylogenetic modelling results
load("Results/complete_species_non_phylo.RData")
load("Results/complete_species_non_phylo_migrants_only.RData")
load("Results/complete_species_non_phylo_residents_only.RData")


# read in phylogenetic modelling results
load("Results/complete_species_phylo.RData")
load("Results/complete_species_phylo_migrants_only.RData")
load("Results/complete_species_phylo_residents_only.RData")


global_model_results <- bind_rows(non_phylo_global_model_results, 
                                  phylo_global_model_results) %>%
  mutate(Species="All") %>%
  bind_rows(non_phylo_global_model_results.migrants %>% mutate(Species="Migrants"),
            phylo_global_model_results.migrants %>% mutate(Species="Migrants")) %>%
  bind_rows(non_phylo_global_model_results.residents %>% mutate(Species="Residents"),
            phylo_global_model_results.residents %>% mutate(Species="Residents"))


model_averaging_results <- bind_rows(non_phylo_model_averaging_results,
                                     phylo_model_averaging_results) %>%
  mutate(Species="All") %>%
  bind_rows(non_phylo_model_averaging_results.migrants %>% mutate(Species="Migrants"),
            phylo_model_averaging_results.migrants %>% mutate(Species="Migrants")) %>%
  bind_rows(non_phylo_model_averaging_results.residents %>% mutate(Species="Residents"),
            phylo_model_averaging_results.residents %>% mutate(Species="Residents"))


# make a plot showing the correlation between
# global model estimates and model averaged parameter estimates
# first for non-phylogenetic analysis
non_phylo_global_model_results %>%
  dplyr::select(term, estimate, lwr_95_confint,
                upr_95_confint, model_type, MONTH) %>%
  bind_rows(., non_phylo_model_averaging_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH)) %>%
  mutate(MONTH=factor(.$MONTH, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))) %>%
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
  mutate(MONTH=factor(.$MONTH, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))) %>%
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


# So it looks like that in general, the model averaging and 
# non-model averaging approaches showed good correlation!


# make a plot to look at difference between parameter estimates for phylo
# and non-phylogenetic models
non_phylo_global_model_results %>%
  dplyr::select(term, estimate, lwr_95_confint,
                upr_95_confint, model_type, MONTH) %>%
  bind_rows(., phylo_global_model_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH)) %>%
  mutate(MONTH=factor(.$MONTH, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))) %>%
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


# let's plot the parameter estimates through the year
# for each of the methods
# first the phylogenetic global model results
phylo_global_model_results %>%
  dplyr::filter(term != "(Intercept)") %>%
  mutate(term=gsub("rescale\\(", "z.", .$term)) %>%
  mutate(term=gsub("\\)", "", .$term)) %>%
  arrange(term, MONTH) %>%
  mutate(month_num=rep(1:12, 11)) %>%
  ggplot(., aes(x=month_num, y=estimate))+
  #geom_line(color="orchid3")+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), width=0.4)+
  geom_point(color="limegreen")+
  facet_wrap(~term, scales="free")+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  ylab("Parameter estimate")+
  xlab("")+
  #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=11), color="orchid3")+
  scale_x_continuous(breaks=c(1:12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  ggtitle("Phylogenetic global model results")


# Non-phylogenetic results
global_model_results %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::filter(model_type=="global_model") %>%
  dplyr::filter(Species=="All") %>%
  arrange(term, MONTH) %>%
  mutate(month_num=case_when(
    MONTH=="Jan" ~ 1,
    MONTH=="Feb" ~ 2,
    MONTH=="Mar" ~ 3,
    MONTH=="Apr" ~ 4,
    MONTH=="May" ~ 5,
    MONTH=="Jun" ~ 6,
    MONTH=="Jul" ~ 7,
    MONTH=="Aug" ~ 8,
    MONTH=="Sep" ~ 9,
    MONTH=="Oct" ~ 10,
    MONTH=="Nov" ~ 11,
    MONTH=="Dec" ~ 12)) %>%
  ggplot(., aes(x=month_num, y=estimate))+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), width=0.4)+
  geom_point(color="limegreen")+
  facet_wrap(~term, scales="free")+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  ylab("Parameter estimate")+
  xlab("")+
  #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=11), color="orchid3")+
  scale_x_continuous(breaks=c(1:12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  ggtitle("Non-phylogenetic global model results")


# phylogenetic model averaging results
phylo_model_averaging_results %>%
  dplyr::filter(term != "(Intercept)") %>%
  mutate(term=gsub("rescale\\(", "z.", .$term)) %>%
  mutate(term=gsub("\\)", "", .$term)) %>%
  arrange(term, MONTH) %>%
  mutate(month_num=case_when(
    MONTH=="Jan" ~ 1,
    MONTH=="Feb" ~ 2,
    MONTH=="Mar" ~ 3,
    MONTH=="Apr" ~ 4,
    MONTH=="May" ~ 5,
    MONTH=="Jun" ~ 6,
    MONTH=="Jul" ~ 7,
    MONTH=="Aug" ~ 8,
    MONTH=="Sep" ~ 9,
    MONTH=="Oct" ~ 10,
    MONTH=="Nov" ~ 11,
    MONTH=="Dec" ~ 12)) %>%
  ggplot(., aes(x=month_num, y=estimate))+
  #geom_line(color="orchid3")+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), width=0.4)+
  geom_point(color="limegreen")+
  facet_wrap(~term, scales="free_y")+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  ylab("Parameter estimate")+
  xlab("")+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=11), color="orchid3")+
  scale_x_continuous(breaks=c(1:12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  ggtitle("Phylogenetic model averaging results")

# Non-phylogenetic model averaging results
non_phylo_model_averaging_results %>%
  dplyr::filter(term != "(Intercept)") %>%
  arrange(term, MONTH) %>%
  mutate(month_num=case_when(
    MONTH=="Jan" ~ 1,
    MONTH=="Feb" ~ 2,
    MONTH=="Mar" ~ 3,
    MONTH=="Apr" ~ 4,
    MONTH=="May" ~ 5,
    MONTH=="Jun" ~ 6,
    MONTH=="Jul" ~ 7,
    MONTH=="Aug" ~ 8,
    MONTH=="Sep" ~ 9,
    MONTH=="Oct" ~ 10,
    MONTH=="Nov" ~ 11,
    MONTH=="Dec" ~ 12)) %>%
  ggplot(., aes(x=month_num, y=estimate))+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), width=0.4)+
  geom_point(color="limegreen")+
  facet_wrap(~term, scales="free")+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  ylab("Parameter estimate")+
  xlab("")+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=11), color="orchid3")+
  scale_x_continuous(breaks=c(1:12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  ggtitle("Non-phylogenetic model averaging results")


# plot the parameter estimates throughout the year
non_phylo_global_model_results %>%
  mutate(model_type="Non-phylo global") %>%
  dplyr::select(term, estimate, lwr_95_confint,
                upr_95_confint, model_type, MONTH) %>%
  bind_rows(., non_phylo_model_averaging_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH) %>%
              mutate(model_type="Non-phylo averaged")) %>%
  bind_rows(., phylo_global_model_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH) %>%
              mutate(model_type="Phylo global")) %>%
  bind_rows(., phylo_model_averaging_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH) %>%
              mutate(model_type="Phylo averaged")) %>%
  dplyr::filter(term != "(Intercept)") %>%
  mutate(term=gsub("rescale\\(", "z.", .$term)) %>%
  mutate(term=gsub("\\)", "", .$term)) %>%
  mutate(term=gsub("z.habitat_generalism_scaled", "z.habitat_generalism", .$term)) %>%
  mutate(MONTH=factor(.$MONTH, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))) %>%
  ggplot(., aes(x=term, y=estimate, color=model_type))+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), width=0.8, position=position_dodge(width=0.6))+
  geom_point(position=position_dodge(width=0.6))+
  coord_flip()+
  facet_wrap(~MONTH)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_hline(yintercept=0, color="red")+
  scale_color_brewer(palette="Set1")+
  xlab("")+
  ylab("Standardized parameter estimate")+
  guides(colour = guide_legend(title="           Model"))

# Same thing as above, but remove the categorical functional group
non_phylo_global_model_results %>%
  mutate(model_type="Non-phylo global") %>%
  dplyr::select(term, estimate, lwr_95_confint,
                upr_95_confint, model_type, MONTH) %>%
  bind_rows(., non_phylo_model_averaging_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH) %>%
              mutate(model_type="Non-phylo averaged")) %>%
  bind_rows(., phylo_global_model_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH) %>%
              mutate(model_type="Phylo global")) %>%
  bind_rows(., phylo_model_averaging_results %>%
              dplyr::select(term, estimate, lwr_95_confint,
                            upr_95_confint, model_type, MONTH) %>%
              mutate(model_type="Phylo averaged")) %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::filter(!term %in% c("functional_dietVertFishScav", "functional_dietPlantSeed",
                             "functional_dietOmnivore", "functional_dietInvertebrate")) %>%
  mutate(term=gsub("rescale\\(", "z.", .$term)) %>%
  mutate(term=gsub("\\)", "", .$term)) %>%
  mutate(term=gsub("z.habitat_generalism_scaled", "z.habitat_generalism", .$term)) %>%
  mutate(MONTH=factor(.$MONTH, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))) %>%
  ggplot(., aes(x=term, y=estimate, color=model_type))+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), width=0.8, position=position_dodge(width=0.6))+
  geom_point(position=position_dodge(width=0.6))+
  coord_flip()+
  facet_wrap(~MONTH)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_hline(yintercept=0, color="red")+
  scale_color_brewer(palette="Set1")+
  xlab("")+
  ylab("Standardized parameter estimate")+
  guides(colour = guide_legend(title="           Model"))

# now look at the difference between migrants and residents
# Non-phylogenetic model averaging results
global_model_results %>%
  dplyr::filter(model_type=="global_model") %>%
  dplyr::filter(Species != "All") %>%
  dplyr::filter(term != "(Intercept)") %>%
  arrange(term, MONTH) %>%
  mutate(month_num=case_when(
    MONTH=="Jan" ~ 1,
    MONTH=="Feb" ~ 2,
    MONTH=="Mar" ~ 3,
    MONTH=="Apr" ~ 4,
    MONTH=="May" ~ 5,
    MONTH=="Jun" ~ 6,
    MONTH=="Jul" ~ 7,
    MONTH=="Aug" ~ 8,
    MONTH=="Sep" ~ 9,
    MONTH=="Oct" ~ 10,
    MONTH=="Nov" ~ 11,
    MONTH=="Dec" ~ 12)) %>%
  ggplot(., aes(x=month_num, y=estimate, color=Species))+
  geom_errorbar(aes(ymin=lwr_95_confint, ymax=upr_95_confint), width=0.4)+
  geom_point()+
  facet_wrap(~term, scales="free")+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  ylab("Parameter estimate")+
  xlab("")+
  #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=11), color="orchid3")+
  scale_x_continuous(breaks=c(1:12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  ggtitle("Non-phylogenetic global model results")

