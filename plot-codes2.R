library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(metafor)



activity <- read_excel("~/Desktop/MetaData.xlsx", sheet = 2)

#####################
# overall effects
#####################
warming_effects_activity <- activity %>% group_by(variable) %>% 
  summarise(meaneffect = (sum(LRR * weight.adjusted) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

warming_effects_activity <- warming_effects_activity %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))

ggplot(warming_effects_activity, aes(factor(variable, levels = c("PO","PHO", "LAP", "NAG", "CBH",
                                                                 "BX", "BG", "Microbial respiration", 
                                                                 "Total respiration")),
                                     meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) + theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Microbial Activity", y = "Effect Size") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text= element_text(size = 11)) +
  scale_x_discrete(breaks = c("PO","PHO", "LAP", "NAG", "CBH",
                              "BX", "BG", "Microbial respiration", 
                              "Total respiration"),
                   labels = c("PO","PHO", "LAP", "NAG", "CBH",
                              "BX", "BG", "Microbial\nRespiration", 
                              "Total\nRespiration")) +
  ylim(-0.4, 0.4) + coord_flip()



#################################
# effects of magnitude on activity
##################################

magnitude_activity <- activity %>% group_by(magnitude.of.warming, variable) %>% 
  summarise(meaneffect = (sum(LRR * weight.adjusted) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

magnitude_activity <- magnitude_activity %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))

ggplot(magnitude_activity, aes(factor(variable, levels = c("PO","PHO", "LAP", "NAG", "CBH",
                                                       "BX", "BG", "Microbial respiration", 
                                                       "Total respiration")),
                               meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.3)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) +
  theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.3)) + coord_flip() +
  facet_wrap(~ magnitude.of.warming) +
  labs(x = "Microbial Activity", y = "Effect Size") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_x_discrete(breaks = c("PO","PHO", "LAP", "NAG", "CBH",
                              "BX", "BG", "Microbial respiration", 
                              "Total respiration"),
                   labels = c("PO","PHO", "LAP", "NAG", "CBH",
                              "BX", "BG", "Microbial\nRespiration", 
                              "Total\nRespiration"))
################################
# ecosystem effects on mactivity 
################################

ecosystem_activity <- activity %>% group_by(ecosystem, variable) %>% 
  summarise(meaneffect = (sum(LRR * weight.adjusted) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

ecosystem_activity <- ecosystem_activity %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))

ggplot(ecosystem_activity, aes(factor(variable, levels = c("PO","PHO", "LAP", "NAG", "CBH",
                                                           "BX", "BG", "Microbial respiration", 
                                                           "Total respiration")),
                               meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.3)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) +
  theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.3)) + coord_flip() +
  facet_wrap(~ ecosystem) +
  labs(x = "Microbial Activity", y = "Effect Size") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_x_discrete(breaks = c("PO","PHO", "LAP", "NAG", "CBH",
                              "BX", "BG", "Microbial respiration", 
                              "Total respiration"),
                   labels = c("PO","PHO", "LAP", "NAG", "CBH",
                              "BX", "BG", "Microbial\nRespiration", 
                              "Total\nRespiration"))
####################
# PH effects
##################
pH_activity <- activity %>% group_by(pH.level, variable) %>% 
  summarise(meaneffect = (sum(LRR * weight.adjusted) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

pH_activity <- pH_activity %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))

ggplot(pH_activity, aes(factor(variable, levels = c("PO","PHO", "LAP", "NAG", "CBH",
                                                           "BX", "BG", "Microbial respiration", 
                                                           "Total respiration")),
                               meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.3)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) +
  theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.3)) + coord_flip() +
  facet_wrap(~ pH.level) +
  labs(x = "Microbial Activity", y = "Effect Size") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_x_discrete(breaks = c("PO","PHO", "LAP", "NAG", "CBH",
                              "BX", "BG", "Microbial respiration", 
                              "Total respiration"),
                   labels = c("PO","PHO", "LAP", "NAG", "CBH",
                              "BX", "BG", "Microbial\nRespiration", 
                              "Total\nRespiration"))
################
# warming methods
##################

methods_activity <- activity %>% group_by(Warming.technique, variable) %>% 
  summarise(meaneffect = (sum(LRR * weight.adjusted) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

methods_activity <- methods_activity %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))
methods_activity <-  methods_activity %>% filter(Warming.technique != "controlled using computers")
methods_activity <-  methods_activity %>% filter(Warming.technique != "Greenhouse")
methods_activity <-  methods_activity %>% filter(Warming.technique != "Laboratory Incubation")


ggplot(methods_activity, aes(factor(variable, levels = c("PO","PHO", "LAP", "NAG", "CBH",
                                                           "BX", "BG", "Microbial respiration", 
                                                           "Total respiration")),
                               meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.3)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) +
  theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.3)) + coord_flip() +
  facet_wrap(~ Warming.technique) +
  labs(x = "Microbial Activity", y = "Effect Size") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_x_discrete(breaks = c("PO","PHO", "LAP", "NAG", "CBH",
                              "BX", "BG", "Microbial respiration", 
                              "Total respiration"),
                   labels = c("PO","PHO", "LAP", "NAG", "CBH",
                              "BX", "BG", "Microbial\nRespiration", 
                              "Total\nRespiration"))


