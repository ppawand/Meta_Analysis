#################
# library
##################
library(ggplot2)
library(metafor)
library(dplyr)
library(tidyr)
library(readxl)
library(party)
library(ggpubr)



# Summarize a data by a provide categorical variable and calculate bootstrapped 
# weighted mean and confidence interval for each category

source("~/Desktop/Meta-analysis/bootstrap mean and CI.R")
source("~/Desktop/Meta-analysis/var_importance.R")

##################
# data
##################
mbiomass <- read_excel("~/Desktop/MetaData.xlsx", sheet = 1)

mbiomass <- mbiomass %>% select(studyID, map, mat, env.type, ph, pH.level, irrigation.status,
                                ecosystem, Warming.technique, magnitude, magnitude.of.warming, 
                                duration, warming.duration, deg.duration,
                                LRR, var, weight, weight.adjusted, weighted.lnR, 
                                variable)

mbiomass[c("studyID", "env.type", "pH.level", "ecosystem", "Warming.technique",
           "irrigation.status", "magnitude.of.warming", "warming.duration")] <-
  lapply(mbiomass[c("studyID", "env.type", "pH.level", "ecosystem", "Warming.technique",
                    "irrigation.status","magnitude.of.warming","warming.duration")], factor)
sapply(mbiomass, class)

# Rename and reorder the factor levels
mbiomass$env.type <- factor(mbiomass$env.type,
                            levels = c("semi-arid", "dry-mesic"),
                            labels = c("Semi-arid", "Dry-mesic"))

mbiomass$pH.level <- factor(mbiomass$pH.level,
                            levels = c("Acidic", "Neutral", "Alkaline"),
                            labels = c("Acidic", "Neutral", "Alkaline"))

mbiomass$magnitude.of.warming <- factor(mbiomass$magnitude.of.warming,
                                    levels = c("< 2", "> 2"),
                                    labels = c("Low", "High"))

mbiomass$warming.duration <- factor(mbiomass$warming.duration,
                                    levels = c("<= 2", "2-5", ">=5"),
                                    labels = c("Short", "Medium", "Long"))

mbiomass$irrigation.status <- factor(mbiomass$irrigation.status,
                                     levels = c("yes", "no"),
                                     labels = c("Irrigated", "Non-irrigated"))

var_importance(mbiomass)

# overall effects of warming

overall_effect <- as.data.frame(bootstrap_mean_ci(data = mbiomass,
                  cat_var = "variable",
                  num_var = "LRR"))


overall_effect <- overall_effect %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                            lci_per = ((exp(lower_ci) - 1) * 100),
                                            uci_per = ((exp(upper_ci) - 1) * 100))

(overall_biomass_plot <- 
  ggplot(overall_effect, aes(factor(category, levels = c("MBC", "Gram Positive Bacteria",
                                                        "Gram Negative Bacteria", "Bacteria",
                                                        "Fungi", "Total Microbial Biomass")),
                            mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Microbial Groups", y = "Percentage Change (%)") +
  scale_x_discrete(breaks = c("MBC", "Gram Positive Bacteria",
                              "Gram Negative Bacteria", "Bacteria",
                              "Fungi", "Total Microbial Biomass"),
                   labels = c("MBC", "G+ Bacterial Lipid", "G- Bacterial Lipid",
                              "Bacterial Lipid", "Fungal Lipid",
                              "Total Lipid")) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


################
# MBC
################

mbc <- mbiomass %>% filter(variable == "MBC")


mbc_ecosystem <- bootstrap_mean_ci(data = mbc,
                  cat_var = "ecosystem",
                  num_var = "LRR")
mbc_ecosystem <- mbc_ecosystem %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                            lci_per = ((exp(lower_ci) - 1) * 100),
                                            uci_per = ((exp(upper_ci) - 1) * 100))

(mbc_ecosystem_plot <- 
  ggplot(mbc_ecosystem, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


mbc_envtype <- bootstrap_mean_ci(data = mbc,
                                   cat_var = "env.type",
                                   num_var = "LRR")
mbc_envtype <- mbc_envtype %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(mbc_envtype_plot <- 
    ggplot(mbc_envtype, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Environment Types", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())



mbc_technique <- bootstrap_mean_ci(data = mbc,
                  cat_var = "Warming.technique",
                  num_var = "LRR")

mbc_technique <- mbc_technique %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(mbc_technique_plot <- 
  ggplot(mbc_technique, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Warming Methods", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

mbc_ph <- bootstrap_mean_ci(data = mbc,
                  cat_var = "pH.level",
                  num_var = "LRR")
mbc_ph <- mbc_ph %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(mbc_ph_plot <- 
  ggplot(mbc_ph, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "pH Levels", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


mbc_magnitude <- bootstrap_mean_ci(data = mbc,
                  cat_var = "magnitude.of.warming",
                  num_var = "LRR")

mbc_magnitude <- mbc_magnitude %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(mbc_magnitude_plot <-
  ggplot(mbc_magnitude, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Magnitude", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


mbc_duration <- bootstrap_mean_ci(data = mbc,
                  cat_var = "warming.duration",
                  num_var = "LRR")

mbc_duration <- mbc_duration %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(mbc_duration_plot <-
  ggplot(mbc_duration, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())



mbc_irrigation <- bootstrap_mean_ci(data = mbc,
                                   cat_var = "irrigation.status",
                                   num_var = "LRR")
mbc_irrigation <- mbc_irrigation %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(mbc_irrigation_plot <- 
    ggplot(mbc_irrigation, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Irrigation", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())




##################
# Bacterial Biomass
#################

bacteria <- mbiomass %>% filter(variable == "Bacteria")


bacteria_ecosystem <- bootstrap_mean_ci(data = bacteria,
                                   cat_var = "ecosystem",
                                   num_var = "LRR")
bacteria_ecosystem <- bacteria_ecosystem %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(bacteria_ecosystem_plot <- 
  ggplot(bacteria_ecosystem, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


bacteria_envtype <- bootstrap_mean_ci(data = bacteria,
                                 cat_var = "env.type",
                                 num_var = "LRR")
bacteria_envtype <- bacteria_envtype %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                      lci_per = ((exp(lower_ci) - 1) * 100),
                                      uci_per = ((exp(upper_ci) - 1) * 100))

(bacteria_envtype_plot <- 
    ggplot(bacteria_envtype, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Environment Types", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())


bacteria_technique <- bootstrap_mean_ci(data = bacteria,
                                   cat_var = "Warming.technique",
                                   num_var = "LRR")

bacteria_technique <- bacteria_technique %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(bacteria_technique_plot <-
  ggplot(bacteria_technique, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Warming Methods", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

bacteria_ph <- bootstrap_mean_ci(data = bacteria,
                            cat_var = "pH.level",
                            num_var = "LRR")
bacteria_ph <- bacteria_ph %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                            lci_per = ((exp(lower_ci) - 1) * 100),
                            uci_per = ((exp(upper_ci) - 1) * 100))

(bacteria_ph_plot <-
  ggplot(bacteria_ph, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "pH Levels", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


bacteria_magnitude <- bootstrap_mean_ci(data = bacteria,
                                   cat_var = "magnitude.of.warming",
                                   num_var = "LRR")

bacteria_magnitude <- bacteria_magnitude %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(bacteria_magnitude_plot <-
  ggplot(bacteria_magnitude, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Magnitude", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


bacteria_duration <- bootstrap_mean_ci(data = bacteria,
                                  cat_var = "warming.duration",
                                  num_var = "LRR")

bacteria_duration <- bacteria_duration %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                        lci_per = ((exp(lower_ci) - 1) * 100),
                                        uci_per = ((exp(upper_ci) - 1) * 100))

(bacteria_duration_plot <-
  ggplot(bacteria_duration, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

bacteria_irrigation <- bootstrap_mean_ci(data = bacteria,
                                    cat_var = "irrigation.status",
                                    num_var = "LRR")
bacteria_irrigation <- bacteria_irrigation %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                            lci_per = ((exp(lower_ci) - 1) * 100),
                                            uci_per = ((exp(upper_ci) - 1) * 100))

(bacteria_irrigation_plot <- 
    ggplot(bacteria_irrigation, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Irrigation", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())

###########################
# Fungal Biomass
###########################
fungi <- mbiomass %>% filter(variable == "Fungi")


fungi_ecosystem <- bootstrap_mean_ci(data = fungi,
                                   cat_var = "ecosystem",
                                   num_var = "LRR")
fungi_ecosystem <- fungi_ecosystem %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(fungi_ecosystem_plot <-
  ggplot(fungi_ecosystem, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

fungi_envtype <- bootstrap_mean_ci(data = fungi,
                                 cat_var = "env.type",
                                 num_var = "LRR")
fungi_envtype <- fungi_envtype %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                      lci_per = ((exp(lower_ci) - 1) * 100),
                                      uci_per = ((exp(upper_ci) - 1) * 100))

(fungi_envtype_plot <- 
    ggplot(fungi_envtype, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Environment Types", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())



fungi_technique <- bootstrap_mean_ci(data = fungi,
                                   cat_var = "Warming.technique",
                                   num_var = "LRR")

fungi_technique <- fungi_technique %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(fungi_technique_plot <-
  ggplot(fungi_technique, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Warming Methods", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

fungi_ph <- bootstrap_mean_ci(data = fungi,
                            cat_var = "pH.level",
                            num_var = "LRR")
fungi_ph <- fungi_ph %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                            lci_per = ((exp(lower_ci) - 1) * 100),
                            uci_per = ((exp(upper_ci) - 1) * 100))

(fungi_ph_plot <- 
  ggplot(fungi_ph, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "pH Levels", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


fungi_magnitude <- bootstrap_mean_ci(data = fungi,
                                   cat_var = "magnitude.of.warming",
                                   num_var = "LRR")

fungi_magnitude <- fungi_magnitude %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(fungi_magnitude_plot <- 
  ggplot(fungi_magnitude, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Magnitude", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


fungi_duration <- bootstrap_mean_ci(data = fungi,
                                  cat_var = "warming.duration",
                                  num_var = "LRR")

fungi_duration <- fungi_duration %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                        lci_per = ((exp(lower_ci) - 1) * 100),
                                        uci_per = ((exp(upper_ci) - 1) * 100))

(fungi_duration_plot <- 
  ggplot(fungi_duration, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

fungi_irrigation <- bootstrap_mean_ci(data = fungi,
                                    cat_var = "irrigation.status",
                                    num_var = "LRR")
fungi_irrigation <- fungi_irrigation %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                            lci_per = ((exp(lower_ci) - 1) * 100),
                                            uci_per = ((exp(upper_ci) - 1) * 100))

(fungi_irrigation_plot <- 
    ggplot(fungi_irrigation, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Irrigation", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())


###########################
# Gram positive bacteria
###########################
Gpos.bacteria <- mbiomass %>% filter(variable == "Gram Positive Bacteria")

gpos_ecosystem <- bootstrap_mean_ci(data = Gpos.bacteria,
                                   cat_var = "ecosystem",
                                   num_var = "LRR")
gpos_ecosystem <- gpos_ecosystem %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(Gpos_ecosystem_plot <-
  ggplot(gpos_ecosystem, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


Gpos_envtype <- bootstrap_mean_ci(data = Gpos.bacteria,
                                 cat_var = "env.type",
                                 num_var = "LRR")
Gpos_envtype <- Gpos_envtype %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                      lci_per = ((exp(lower_ci) - 1) * 100),
                                      uci_per = ((exp(upper_ci) - 1) * 100))

(Gpos_envtype_plot <- 
    ggplot(Gpos_envtype, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Environment Types", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())



gpos_technique <- bootstrap_mean_ci(data = Gpos.bacteria,
                                   cat_var = "Warming.technique",
                                   num_var = "LRR")

gpos_technique <- gpos_technique %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(Gpos_technique_plot <-
  ggplot(gpos_technique, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Warming Methods", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

gpos_ph <- bootstrap_mean_ci(data = Gpos.bacteria,
                            cat_var = "pH.level",
                            num_var = "LRR")
gpos_ph <- gpos_ph %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                            lci_per = ((exp(lower_ci) - 1) * 100),
                            uci_per = ((exp(upper_ci) - 1) * 100))

(Gpos_ph_plot <- 
  ggplot(gpos_ph, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "pH Levels", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


gpos_magnitude <- bootstrap_mean_ci(data = Gpos.bacteria,
                                   cat_var = "magnitude.of.warming",
                                   num_var = "LRR")

gpos_magnitude <- gpos_magnitude %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(Gpos_magnitude_plot <-
  ggplot(gpos_magnitude, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Magnitude", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


gpos_duration <- bootstrap_mean_ci(data = Gpos.bacteria,
                                  cat_var = "warming.duration",
                                  num_var = "LRR")

gpos_duration <- gpos_duration %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                        lci_per = ((exp(lower_ci) - 1) * 100),
                                        uci_per = ((exp(upper_ci) - 1) * 100))

(Gpos_duration_plot <- 
  ggplot(gpos_duration, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

Gpos_irrigation <- bootstrap_mean_ci(data = Gpos.bacteria,
                                    cat_var = "irrigation.status",
                                    num_var = "LRR")
Gpos_irrigation <- Gpos_irrigation %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                            lci_per = ((exp(lower_ci) - 1) * 100),
                                            uci_per = ((exp(upper_ci) - 1) * 100))

(Gpos_irrigation_plot <- 
    ggplot(Gpos_irrigation, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Irrigation", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())


###########################
# Gram negative bacteria
###########################
Gneg.bacteria <- mbiomass %>% filter(variable == "Gram Negative Bacteria")

gneg_ecosystem <- bootstrap_mean_ci(data = Gneg.bacteria,
                                   cat_var = "ecosystem",
                                   num_var = "LRR")
gneg_ecosystem <- gneg_ecosystem %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(Gneg_ecosystem_plot <- 
  ggplot(gneg_ecosystem, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


Gneg_envtype <- bootstrap_mean_ci(data = Gneg.bacteria,
                                 cat_var = "env.type",
                                 num_var = "LRR")
Gneg_envtype <- Gneg_envtype %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                      lci_per = ((exp(lower_ci) - 1) * 100),
                                      uci_per = ((exp(upper_ci) - 1) * 100))

(Gneg_envtype_plot <- 
    ggplot(Gneg_envtype, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Environment Types", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())



gneg_technique <- bootstrap_mean_ci(data = Gneg.bacteria,
                                   cat_var = "Warming.technique",
                                   num_var = "LRR")

gneg_technique <- gneg_technique %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(Gneg_technique_plot <- 
  ggplot(gneg_technique, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Warming Methods", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

gneg_ph <- bootstrap_mean_ci(data = Gneg.bacteria,
                            cat_var = "pH.level",
                            num_var = "LRR")
gneg_ph <- gneg_ph %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                            lci_per = ((exp(lower_ci) - 1) * 100),
                            uci_per = ((exp(upper_ci) - 1) * 100))

(Gneg_ph_plot <-
  ggplot(gneg_ph, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "pH Levels", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


gneg_magnitude <- bootstrap_mean_ci(data = Gneg.bacteria,
                                   cat_var = "magnitude.of.warming",
                                   num_var = "LRR")

gneg_magnitude <- gneg_magnitude %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(Gneg_magnitude_plot <-
  ggplot(gneg_magnitude, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Magnitude", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


gneg_duration <- bootstrap_mean_ci(data = Gneg.bacteria,
                                  cat_var = "warming.duration",
                                  num_var = "LRR")

gneg_duration <- gneg_duration %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                        lci_per = ((exp(lower_ci) - 1) * 100),
                                        uci_per = ((exp(upper_ci) - 1) * 100))

(Gneg_duration_plot <-
  ggplot(gneg_duration, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

Gneg_irrigation <- bootstrap_mean_ci(data = Gneg.bacteria,
                                    cat_var = "irrigation.status",
                                    num_var = "LRR")
Gneg_irrigation <- Gneg_irrigation %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                            lci_per = ((exp(lower_ci) - 1) * 100),
                                            uci_per = ((exp(upper_ci) - 1) * 100))

(Gneg_irrigation_plot <- 
    ggplot(Gneg_irrigation, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Irrigation", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())


#######################
# total biomass
#####################
total.mbiomass <- mbiomass %>% filter(variable == "Total Microbial Biomass")


tbiomass_ecosystem <- bootstrap_mean_ci(data = total.mbiomass,
                                   cat_var = "ecosystem",
                                   num_var = "LRR")
tbiomass_ecosystem <- tbiomass_ecosystem %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(tbiomass_ecosystem_plot <-
  ggplot(tbiomass_ecosystem, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


tbiomass_envtype <- bootstrap_mean_ci(data = total.mbiomass,
                                 cat_var = "env.type",
                                 num_var = "LRR")
tbiomass_envtype <- tbiomass_envtype %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                      lci_per = ((exp(lower_ci) - 1) * 100),
                                      uci_per = ((exp(upper_ci) - 1) * 100))

(tbiomass_envtype_plot <- 
    ggplot(tbiomass_envtype, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Environment Types", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())



tbiomass_technique <- bootstrap_mean_ci(data = total.mbiomass,
                                   cat_var = "Warming.technique",
                                   num_var = "LRR")

tbiomass_technique <- tbiomass_technique %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(tbiomass_technique_plot <-
  ggplot(tbiomass_technique, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Warming Methods", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

tbiomass_ph <- bootstrap_mean_ci(data = total.mbiomass,
                            cat_var = "pH.level",
                            num_var = "LRR")
tbiomass_ph <- tbiomass_ph %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                            lci_per = ((exp(lower_ci) - 1) * 100),
                            uci_per = ((exp(upper_ci) - 1) * 100))

(tbiomass_ph_plot <-
  ggplot(tbiomass_ph, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "pH Levels", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


tbiomass_magnitude <- bootstrap_mean_ci(data = total.mbiomass,
                                   cat_var = "magnitude.of.warming",
                                   num_var = "LRR")

tbiomass_magnitude <- tbiomass_magnitude %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                          lci_per = ((exp(lower_ci) - 1) * 100),
                                          uci_per = ((exp(upper_ci) - 1) * 100))

(tbiomass_magnitude_plot <-
  ggplot(tbiomass_magnitude, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Magnitude", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


tbiomass_duration <- bootstrap_mean_ci(data = total.mbiomass,
                                  cat_var = "warming.duration",
                                  num_var = "LRR")

tbiomass_duration <- tbiomass_duration %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                        lci_per = ((exp(lower_ci) - 1) * 100),
                                        uci_per = ((exp(upper_ci) - 1) * 100))

(tbiomass_duration_plot <-
  ggplot(tbiomass_duration, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

tbiomass_irrigation <- bootstrap_mean_ci(data = total.mbiomass,
                                    cat_var = "irrigation.status",
                                    num_var = "LRR")
tbiomass_irrigation <- tbiomass_irrigation %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                            lci_per = ((exp(lower_ci) - 1) * 100),
                                            uci_per = ((exp(upper_ci) - 1) * 100))

(tbiomass_irrigation_plot <- 
    ggplot(tbiomass_irrigation, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Irrigation", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())


# best.model.total.biomass <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
#                                    mods = ~ ecosystem + ph + ahm + Warming.technique +
#                                      magnitude + duration,
#                                    random = ~ 1 | studyID,
#                                    data = total.mbiomass,
#                                    method = "ML")
# 
# summary(best.model.total.biomass)
# 
# 
# results.total.biomass <- data.frame(coeff = best.model.total.biomass$b, 
#                            lower.CI = best.model.total.biomass$ci.lb, 
#                            upper.CI = best.model.total.biomass$ci.ub,
#                            p.value = best.model.total.biomass$pval)
# 
# funnel(best.model.total.biomass, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

##################################################################
# activity
##################################################################
mfunction <- read_excel("~/Desktop/MetaData.xlsx", sheet = 2)



mfunction <- mfunction %>% select(studyID, map, mat, env.type, ph, pH.level, irrigation.status,
                                ecosystem, Warming.technique, magnitude, magnitude.of.warming, 
                                duration, warming.duration, deg.duration,
                                LRR, var, weight, weight.adjusted, weighted.lnR, 
                                variable)

mfunction[c("studyID", "env.type", "pH.level", "ecosystem", "Warming.technique",
           "irrigation.status", "magnitude.of.warming", "warming.duration")] <-
  lapply(mfunction[c("studyID", "env.type", "pH.level", "ecosystem", "Warming.technique",
                    "irrigation.status","magnitude.of.warming","warming.duration")], factor)
sapply(mfunction, class)

# Rename and reorder the factor levels
mfunction$env.type <- factor(mfunction$env.type,
                            levels = c("semi-arid", "dry-mesic"),
                            labels = c("Semi-arid", "Dry-mesic"))

mfunction$pH.level <- factor(mfunction$pH.level,
                            levels = c("Acidic", "Neutral", "Alkaline"),
                            labels = c("Acidic", "Neutral", "Alkaline"))

mfunction$magnitude.of.warming <- factor(mfunction$magnitude.of.warming,
                                        levels = c("< 2", "> 2"),
                                        labels = c("Low", "High"))

mfunction$warming.duration <- factor(mfunction$warming.duration,
                                    levels = c("<=2", "2-5", ">=5"),
                                    labels = c("Short", "Medium", "Long"))

mfunction$irrigation.status <- factor(mfunction$irrigation.status,
                                     levels = c("yes", "no"),
                                     labels = c("Irrigated", "Non-irrigated"))

var_importance(mfunction)


##################
# Soil Respiration
##################

total.respiration <- mfunction %>% filter(variable == "Total respiration")


tresp_ecosystem <- bootstrap_mean_ci(data = total.respiration,
                                        cat_var = "ecosystem",
                                        num_var = "LRR")
tresp_ecosystem <- tresp_ecosystem %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                                    lci_per = ((exp(lower_ci) - 1) * 100),
                                                    uci_per = ((exp(upper_ci) - 1) * 100))

(tresp_ecosystem_plot <-
  ggplot(tresp_ecosystem, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())



tresp_envtype <- bootstrap_mean_ci(data = total.respiration,
                                 cat_var = "env.type",
                                 num_var = "LRR")
tresp_envtype <- tresp_envtype %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                      lci_per = ((exp(lower_ci) - 1) * 100),
                                      uci_per = ((exp(upper_ci) - 1) * 100))

(tresp_envtype_plot <- 
    ggplot(tresp_envtype, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Environment Types", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())



tresp_technique <- bootstrap_mean_ci(data = total.respiration,
                                        cat_var = "Warming.technique",
                                        num_var = "LRR")

tresp_technique <- tresp_technique %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                                    lci_per = ((exp(lower_ci) - 1) * 100),
                                                    uci_per = ((exp(upper_ci) - 1) * 100))

(tresp_technique_plot <-
  ggplot(tresp_technique, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Warming Methods", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

tresp_ph <- bootstrap_mean_ci(data = total.respiration,
                                 cat_var = "pH.level",
                                 num_var = "LRR")
tresp_ph <- tresp_ph %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                      lci_per = ((exp(lower_ci) - 1) * 100),
                                      uci_per = ((exp(upper_ci) - 1) * 100))

(tresp_ph_plot <-
  ggplot(tresp_ph, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "pH Levels", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


tresp_magnitude <- bootstrap_mean_ci(data = total.respiration,
                                        cat_var = "magnitude.of.warming",
                                        num_var = "LRR")

tresp_magnitude <- tresp_magnitude %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                                    lci_per = ((exp(lower_ci) - 1) * 100),
                                                    uci_per = ((exp(upper_ci) - 1) * 100))

(tresp_magnitude_plot <-
  ggplot(tresp_magnitude, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Magnitude", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


tresp_duration <- bootstrap_mean_ci(data = total.respiration,
                                       cat_var = "warming.duration",
                                       num_var = "LRR")

tresp_duration <- tresp_duration %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                                  lci_per = ((exp(lower_ci) - 1) * 100),
                                                  uci_per = ((exp(upper_ci) - 1) * 100))

(tresp_duration_plot <-
  ggplot(tresp_duration, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


tresp_irrigation <- bootstrap_mean_ci(data = total.respiration,
                                    cat_var = "irrigation.status",
                                    num_var = "LRR")
tresp_irrigation <- tresp_irrigation %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                            lci_per = ((exp(lower_ci) - 1) * 100),
                                            uci_per = ((exp(upper_ci) - 1) * 100))

(tresp_irrigation_plot <- 
    ggplot(tresp_irrigation, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Irrigation", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())



# 
# 
# best.model.total.respiration <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
#                                        mods = ~ ecosystem + ph + ahm + Warming.technique +
#                                          magnitude + duration,
#                                        random = ~ 1 | studyID,
#                                        data = total.respiration,
#                                        method = "ML")
# 
# summary(best.model.total.respiration)
# 
# results.total.respiration <- data.frame(coeff = best.model.total.respiration$b, 
#                                         lower.CI = best.model.total.respiration$ci.lb, 
#                                         upper.CI = best.model.total.respiration$ci.ub,
#                                         p.value = best.model.total.respiration$pval)
# funnel(best.model.total.respiration,xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

############################
# Soil Microbial Respiration
############################

microbial.respiration <- mfunction %>% filter(variable == "Microbial respiration")


mresp_ecosystem <- bootstrap_mean_ci(data = microbial.respiration,
                                        cat_var = "ecosystem",
                                        num_var = "LRR")
mresp_ecosystem <- mresp_ecosystem %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                                    lci_per = ((exp(lower_ci) - 1) * 100),
                                                    uci_per = ((exp(upper_ci) - 1) * 100))

(mresp_ecosystem_plot <- 
  ggplot(mresp_ecosystem, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

mresp_envtype <- bootstrap_mean_ci(data = microbial.respiration,
                                 cat_var = "env.type",
                                 num_var = "LRR")
mresp_envtype <- mresp_envtype %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                      lci_per = ((exp(lower_ci) - 1) * 100),
                                      uci_per = ((exp(upper_ci) - 1) * 100))

(mresp_envtype_plot <- 
    ggplot(mresp_envtype, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Environment Types", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())



mresp_technique <- bootstrap_mean_ci(data = microbial.respiration,
                                        cat_var = "Warming.technique",
                                        num_var = "LRR")

mresp_technique <- mresp_technique %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                                    lci_per = ((exp(lower_ci) - 1) * 100),
                                                    uci_per = ((exp(upper_ci) - 1) * 100))

(mresp_technique_plot <- 
  ggplot(mresp_technique, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 3)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Warming Methods", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

mresp_ph <- bootstrap_mean_ci(data = microbial.respiration,
                                 cat_var = "pH.level",
                                 num_var = "LRR")
mresp_ph <- mresp_ph %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                      lci_per = ((exp(lower_ci) - 1) * 100),
                                      uci_per = ((exp(upper_ci) - 1) * 100))

(mresp_ph_plot <-
  ggplot(mresp_ph, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 3)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "pH Levels", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())


mresp_magnitude <- bootstrap_mean_ci(data = microbial.respiration,
                                        cat_var = "magnitude.of.warming",
                                        num_var = "LRR")

mresp_magnitude <- mresp_magnitude %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                                    lci_per = ((exp(lower_ci) - 1) * 100),
                                                    uci_per = ((exp(upper_ci) - 1) * 100))

(mresp_magnitude_plot <- 
  ggplot(mresp_magnitude, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Magnitude", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())



mresp_duration <- bootstrap_mean_ci(data = microbial.respiration,
                                       cat_var = "warming.duration",
                                       num_var = "LRR")

mresp_duration <- mresp_duration %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                                  lci_per = ((exp(lower_ci) - 1) * 100),
                                                  uci_per = ((exp(upper_ci) - 1) * 100))

(mresp_duration_plot <-
  ggplot(mresp_duration, aes(category, mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Ecosystem types", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())

mresp_irrigation <- bootstrap_mean_ci(data = microbial.respiration,
                                    cat_var = "irrigation.status",
                                    num_var = "LRR")
mresp_irrigation <- mresp_irrigation %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                            lci_per = ((exp(lower_ci) - 1) * 100),
                                            uci_per = ((exp(upper_ci) - 1) * 100))

(mresp_irrigation_plot <- 
    ggplot(mresp_irrigation, aes(category, mean_per)) + 
    geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
    geom_text(aes(label= n, vjust = 2)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                  width = 0.2, linewidth = 0.7,
                  position = position_dodge(width = 0.35)) +
    labs(x = "Irrigation", y = "Percentage Change (%)") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text= element_text(size = 11),
          axis.title.y = element_blank())  + coord_flip())


# 
# 
# best.model.microbial.respiration <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
#                                            mods = ~ ecosystem + ph + ahm + Warming.technique +
#                                              magnitude + duration,
#                                            random = ~ 1 | studyID,
#                                            data = microbial.respiration,
#                                            method = "ML")
# 
# summary(best.model.microbial.respiration)
# 
# results.microbial.respiration <- data.frame(coeff = best.model.microbial.respiration$b, 
#                                             lower.CI = best.model.microbial.respiration$ci.lb, 
#                                             upper.CI = best.model.microbial.respiration$ci.ub,
#                                             p.value = best.model.microbial.respiration$pval)
# 
# funnel(best.model.microbial.respiration,xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)


# overall effects of warming on microbial activity

overall_effect_activity <- as.data.frame(bootstrap_mean_ci(data = mfunction,
                                                  cat_var = "variable",
                                                  num_var = "LRR"))



overall_effect_activity <- overall_effect_activity %>% mutate(mean_per = ((exp(mean) - 1) * 100),
                                            lci_per = ((exp(lower_ci) - 1) * 100),
                                            uci_per = ((exp(upper_ci) - 1) * 100))

(overall_activity_plot <-
  ggplot(overall_effect_activity, aes(factor(category, 
                                           levels = c ("Oxidases","N-hydrolysis",
                                                       "C-hydrolysis",
                                                      "Microbial respiration", 
                                                      "Total respiration")),
                                    mean_per)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + theme_classic() +
  geom_text(aes(label= n, vjust = 2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lci_per, ymax = uci_per),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Microbial activity", y = "Percentage Change (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank())  + coord_flip())






















# ###########################
# # Enzyme activity
# #######################
# 
# 
# 
# beta.glucosidase <- mfunction %>% filter(variable == "BG")
# 
# 
# best.model.beta.glucosidase <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
#                                       mods = ~ ecosystem + ph + ahm + Warming.technique +
#                                         magnitude + duration,
#                                       random = ~ 1 | studyID,
#                                       data = beta.glucosidase,
#                                       method = "ML")
# 
# summary(best.model.beta.glucosidase)
# 
# results.beta.glucosidase <- data.frame(coeff = best.model.beta.glucosidase$b, 
#                                        lower.CI = best.model.beta.glucosidase$ci.lb, 
#                                        upper.CI = best.model.beta.glucosidase$ci.ub,
#                                        p.value = best.model.beta.glucosidase$pval)
# funnel(best.model.beta.glucosidase,xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)
# ## Cellobiohydrolase
# 
# cellobiohydrolase <- mfunction %>% filter(variable == "CBH")
# 
# 
# best.model.cellobiohydrolase <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
#                                        mods = ~ ecosystem + ph + ahm + Warming.technique +
#                                          magnitude + duration,
#                                        random = ~ 1 | studyID,
#                                        data = cellobiohydrolase,
#                                        method = "ML")
# 
# summary(best.model.cellobiohydrolase)
# 
# results.cellobiohyrolase <- data.frame(coeff = best.model.cellobiohydrolase$b, 
#                                        lower.CI = best.model.cellobiohydrolase$ci.lb, 
#                                        upper.CI = best.model.cellobiohydrolase$ci.ub,
#                                        p.value = best.model.cellobiohydrolase$pval)
# 
# funnel(best.model.cellobiohydrolase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)
# 
# ## N- acetyl-beta- Glucosaminidase
# 
# nag <- mfunction %>% filter(variable == "NAG")
# 
# 
# best.model.nag <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
#                          mods = ~ ecosystem + ph + ahm + Warming.technique +
#                            magnitude + duration,
#                          random = ~ 1 | studyID,
#                          data = nag,
#                          method = "ML")
# summary(best.model.nag)
# 
# results.nag <- data.frame(coeff = best.model.nag$b, 
#                           lower.CI = best.model.nag$ci.lb, 
#                           upper.CI = best.model.nag$ci.ub,
#                           p.value = best.model.nag$pval)
# funnel(best.model.nag, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)
# 
# ## beta-Xylosidase
# 
# beta.xylosidase <- mfunction %>% filter(variable == "BX")
# 
# 
# best.model.beta.xylosidase <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
#                                      mods = ~ ecosystem + ph + ahm + Warming.technique +
#                                        magnitude + duration,
#                                      random = ~ 1 | studyID,
#                                      data = beta.xylosidase,
#                                      method = "ML")
# 
# summary(best.model.beta.xylosidase)
# 
# results.beta.xylosidase<- data.frame(coeff = best.model.beta.xylosidase$b, 
#                                      lower.CI = best.model.beta.xylosidase$ci.lb, 
#                                      upper.CI = best.model.beta.xylosidase$ci.ub,
#                                      p.value = best.model.beta.xylosidase$pval)
# funnel(best.model.beta.xylosidase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)
# 
# ## Leucine aminopeptidase
# 
# lap <- mfunction %>% filter(variable == "LAP")
# 
# 
# best.model.lap <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
#                          mods = ~ ecosystem + ph + ahm + Warming.technique +
#                            magnitude + duration,
#                          random = ~ 1 | studyID,
#                          data = lap,
#                          method = "ML")
# 
# summary(best.model.lap)
# 
# results.lap <- data.frame(coeff = best.model.lap$b, 
#                           lower.CI = best.model.lap$ci.lb, 
#                           upper.CI = best.model.lap$ci.ub,
#                           p.value = best.model.lap$pval)
# funnel(best.model.lap, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)
# 
# 
# ## Phenol Oxidase
# 
# phenol.oxidase <- mfunction %>% filter(variable == "PHO")
# 
# 
# best.model.phenol.oxidase <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
#                                     mods = ~ ecosystem + ph + ahm + Warming.technique +
#                                       magnitude + duration,
#                                     random = ~ 1 | studyID,
#                                     data = phenol.oxidase,
#                                     method = "ML")
# 
# summary(best.model.phenol.oxidase)
# 
# results.phenol.oxidase <- data.frame(coeff = best.model.phenol.oxidase$b, 
#                                      lower.CI = best.model.phenol.oxidase$ci.lb, 
#                                      upper.CI = best.model.phenol.oxidase$ci.ub,
#                                      p.value = best.model.phenol.oxidase$pval)
# 
# 
# funnel(best.model.phenol.oxidase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)
# 
# ## perooxidase
# 
# perooxidase <- mfunction %>% filter(variable == "PO")
# 
# best.model.perooxidase <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
#                                  mods = ~ ecosystem + ph + ahm + Warming.technique +
#                                    magnitude + duration,
#                                  random = ~ 1 | studyID,
#                                  data = perooxidase,
#                                  method = "ML")
# 
# summary(best.model.perooxidase)
# 
# results.perooxidase <- data.frame(coeff = best.model.perooxidase$b, 
#                                   lower.CI = best.model.perooxidase$ci.lb, 
#                                   upper.CI = best.model.perooxidase$ci.ub,
#                                   p.value = best.model.perooxidase$pval)
# funnel(best.model.perooxidase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)
# 



















