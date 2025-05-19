library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
class(data$weight.adjusted)
 
#####################################
# Bootstrap mean and CI function
#####################################
 bootstrap_mean_ci <- function (data, cat_var, num_var, nIter = 5000, confprob = 0.05) {
 # get a unique values of categorical variable
   cat_levels <- unique(data[[cat_var]])
 
   # creating empty data frame to store results
   results <- data.frame(category = character(),
                         mean = numeric(),
                         lower_ci = numeric(),
                         upper_ci = numeric())
  
  for (cat in cat_levels) {
    # subset data for a current category
    
    cat_data <- data[data[[cat_var]] == cat, num_var]
    cat_data <- pull(cat_data, num_var)
    # perform bootstrapping
    boot_sample <- replicate(nIter, mean(sample(cat_data, replace = TRUE)))
    ci <- quantile(boot_sample, c(confprob/2, 1-confprob/2))
    
    results <- rbind(results, data.frame(category = cat,
                                         mean = mean(cat_data),
                                         lower_ci = ci[1],
                                         upper_ci = ci[2]))
                                         
  }
  
  return(results)
 }




data <- read_excel("~/Desktop/MetaData.xlsx", sheet = 1)

cat_data <- data[data[["ecosystem"]] == "Grassland", "weighted.lnR"]
replicate(1000, mean(sample(cat_data, replace = TRUE)))

bootstrap_mean_ci(data = mbc, cat_var = "ecosystem", num_var = "LRR")


##################
# overall effects
#################
warming_effects <- bacteria %>% group_by(Warming.technique) %>% 
  summarise(meaneffect = (sum(weighted.lnR) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

warming_effects <- warming_effects %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))

ggplot(warming_effects, aes(factor(variable, levels = c("MBC","FB ratio", "Gram Positive Bacteria",
                                                      "Gram Negative Bacteria", "Bacteria",
                                                      "Fungi", "Total Microbial Biomass")),
                          meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) + theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Microbial Biomass", y = "Effect Size") +
  scale_x_discrete(breaks = c("MBC","FB ratio", "Gram Positive Bacteria",
                              "Gram Negative Bacteria", "Bacteria",
                              "Fungi", "Total Microbial Biomass"),
                   labels = c("MBC", "FB ratio", "G+ Bacteria", "G- Bacteria", "Bacteria", "Fungi",
                              "Total Biomass")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank()) +
  ylim(-0.2, 0.3) + coord_flip()

##############################
# effects of warming magnitude
##############################

warming_level <- data %>% group_by(magnitude.of.warming, variable) %>% 
  summarise(meaneffect = (sum(LRR * weight.adjusted) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

warming_level <- warming_level %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))

ggplot(warming_level, aes(factor(variable, levels = c("MBC","FB ratio", "Gram Positive Bacteria",
                                                      "Gram Negative Bacteria", "Bacteria",
                                                      "Fungi", "Total Microbial Biomass")),
                                 meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) + theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.35)) +
  labs(x = "Microbial Biomass", y = "Effect Size") +
  scale_x_discrete(breaks = c("MBC","FB ratio", "Gram Positive Bacteria",
                                "Gram Negative Bacteria", "Bacteria",
                                "Fungi", "Total Microbial Biomass"),
                   labels = c("MBC", "FB ratio", "G+ Bacteria", "G- Bacteria", "Bacteria", "Fungi",
                              "Total Biomass")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  facet_wrap(~ magnitude.of.warming) + coord_flip()

#########################################
# warming effects in different ecosystem 
#########################################
ecosystem <- data %>% group_by(ecosystem, variable) %>% 
  summarise(meaneffect = (sum(LRR * weight.adjusted) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

ecosystem <- ecosystem %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))


ggplot(ecosystem, aes(factor(variable, levels = c("MBC","FB ratio", "Gram Positive Bacteria",
                                                      "Gram Negative Bacteria", "Bacteria",
                                                      "Fungi", "Total Microbial Biomass")),
                          meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.5)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) +
  theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.5)) +
  labs(x = "Microbial Biomass", y = "Effect Size") +
  scale_x_discrete(breaks = c("MBC","FB ratio", "Gram Positive Bacteria",
                              "Gram Negative Bacteria", "Bacteria",
                              "Fungi", "Total Microbial Biomass"),
                   labels = c("MBC", "FB ratio", "G+ Bacteria", "G- Bacteria", "Bacteria", "Fungi",
                              "Total Biomass")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text = element_text(size = 11),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  facet_wrap(~ ecosystem) + coord_flip()

##########################
# pH effects on microbes
#########################

pH <- data %>% group_by(pH.level, variable) %>% 
  summarise(meaneffect = (sum(LRR * weight.adjusted) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

pH<- pH %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))

ggplot(pH, aes(factor(variable, levels = c("MBC","FB ratio", "Gram Positive Bacteria",
                                                      "Gram Negative Bacteria", "Bacteria",
                                                      "Fungi", "Total Microbial Biomass")),
                          meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) + theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.35)) + coord_flip() +
  labs(x = "Microbial Biomass", y = "Effect Size") +
  scale_x_discrete(breaks = c("MBC","FB ratio", "Gram Positive Bacteria",
                              "Gram Negative Bacteria", "Bacteria",
                              "Fungi", "Total Microbial Biomass"),
                   labels = c("MBC", "FB ratio", "G+ Bacteria", "G- Bacteria", "Bacteria", "Fungi",
                              "Total Biomass")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  facet_wrap(~ pH.level)

#############################
# effects of warming duration
#############################

duration <- data %>% group_by(warming.duration, variable) %>% 
  summarise(meaneffect = (sum(LRR * weight.adjusted) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

duration <- duration %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))

ggplot(duration, aes(factor(variable, levels = c("MBC","FB ratio", "Gram Positive Bacteria",
                                           "Gram Negative Bacteria", "Bacteria",
                                           "Fungi", "Total Microbial Biomass")),
               meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) + theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.35)) + coord_flip() +
  labs(x = "Microbial Biomass", y = "Effect Size") +
  scale_x_discrete(breaks = c("MBC","FB ratio", "Gram Positive Bacteria",
                              "Gram Negative Bacteria", "Bacteria",
                              "Fungi", "Total Microbial Biomass"),
                   labels = c("MBC", "FB ratio", "G+ Bacteria", "G- Bacteria", "Bacteria", "Fungi",
                              "Total Biomass")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  facet_wrap(~ warming.duration)
#################################
# effects of warming technique
################################
techniques <- data %>% group_by(Warming.technique, variable) %>% 
  summarise(meaneffect = (sum(LRR * weight.adjusted) / sum(weight.adjusted)),
            vari = 1/sum(weight.adjusted),
            n = n(),
            lower.ci = meaneffect - 1.96 * sqrt(vari),
            upper.ci = meaneffect + 1.96 * sqrt (vari))

techniques <- techniques %>% mutate(percent_change = ((exp(meaneffect) - 1) * 100))

techniques <-  techniques %>% filter(Warming.technique != "controlled using computers")
techniques <-  techniques %>% filter(Warming.technique != "Greenhouse")
techniques <-  techniques %>% filter(Warming.technique != "Laboratory Incubation")

ggplot(techniques, aes(factor(variable, levels = c("MBC","FB ratio", "Gram Positive Bacteria",
                                                 "Gram Negative Bacteria", "Bacteria",
                                                 "Fungi", "Total Microbial Biomass")),
                     meaneffect)) + 
  geom_point(size = 3.0, position = position_dodge(width = 0.35)) + 
  geom_text(aes(label = n), size = 3, vjust = -1, hjust = 0.5) + theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.35)) + coord_flip() +
  labs(x = "Microbial Biomass", y = "Effect Size") +
  scale_x_discrete(breaks = c("MBC","FB ratio", "Gram Positive Bacteria",
                              "Gram Negative Bacteria", "Bacteria",
                              "Fungi", "Total Microbial Biomass"),
                   labels = c("MBC", "FB ratio", "G+ Bacteria", "G- Bacteria", "Bacteria", "Fungi",
                              "Total Biomass")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text= element_text(size = 11),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  facet_wrap(~ Warming.technique)
