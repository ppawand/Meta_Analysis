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


################
# function
################
# rank the predictor based on their relative predictive power

var_importance <- function (data) {
  cforest <- cforest(LRR ~ ahm + pH.level + ecosystem +
                       magnitude.of.warming + warming.duration + Warming.technique,
                     data = data,
                     controls = (cforest_control(ntree = 1000, replace = T)))
  
  varimp <- varimp(cforest)
  plot(varimp)
  ranked.importance <- as.data.frame(sort(varimp, decreasing = T))
  names(ranked.importance) <- "index"
  ranked.importance$predictor <- row.names(ranked.importance)
  ranked.importance$predictor <- factor(ranked.importance$predictor,
                                        levels = ranked.importance$predictor[order(ranked.importance$index)])
  ranked_plot <- ggplot(ranked.importance, aes(x = predictor, y = index)) +
    geom_bar(stat = "identity", fill = "black", col = "black") +
    coord_flip() +
    theme_bw()
  return(ranked_plot)
}


##################
# data
##################
mbiomass <- read_excel("~/Desktop/MetaData.xlsx", sheet = 1)

mbiomass <- mbiomass %>% select(studyID, map, mat, ahm, ph, pH.level, ecosystem,
                                Warming.technique, magnitude, magnitude.of.warming, 
                                duration, warming.duration, deg.duration,
                                LRR, var, weight, weight.adjusted, weighted.lnR, 
                                variable)

mbiomass[c("studyID", "pH.level", "ecosystem", "Warming.technique", "magnitude.of.warming",
           "warming.duration")] <-
  lapply(mbiomass[c("studyID", "pH.level", "ecosystem", "Warming.technique", "magnitude.of.warming",
                    "warming.duration")], factor)
sapply(mbiomass, class)


# mbc

mbc <- mbiomass %>% filter(variable == "MBC")
var_importance(mbc)


mod.mbc <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                mods = ~ ahm + pH.level + warming.duration + ecosystem + Warming.technique +
                  magnitude.of.warming,
                random = ~ 1 | studyID,
                data = mbc)

summary(mod.mbc)

funnel(mod.mbc, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)


# bacteria

bacteria <- mbiomass %>% filter(variable == "Bacteria")
var_importance(bacteria)


mod.bacteria <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                  mods = ~ warming.duration + Warming.technique + ahm + pH.level,
                  random = ~ 1 | studyID,
                  data = bacteria)

summary(mod.bacteria)

funnel(mod.bacteria, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)


# fungi

fungi  <- mbiomass %>% filter(variable == "Fungi")
var_importance(fungi)


mod.fungi <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                       mods = ~ pH.level + magnitude.of.warming + warming.duration,
                       random = ~ 1 | studyID,
                       data = fungi)

summary(mod.fungi)

funnel(mod.fungi, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)


# Gram Positive bacteria

Gpos.bacteria <- mbiomass %>% filter(variable == "Gram Positive Bacteria")

var_importance(Gpos.bacteria)
# all predictors have negative variable importance



# Gram Negative bacteria

Gneg.bacteria <- mbiomass %>% filter(variable == "Gram Negative Bacteria")

var_importance(Gneg.bacteria)

mod.Gneg <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                    mods = ~ ahm + magnitude.of.warming + ecosystem +
                     Warming.technique,
                    random = ~ 1 | studyID,
                    data = Gneg.bacteria)

summary(mod.Gneg)

funnel(mod.Gneg, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

# Total biomass

total.mbiomass <- mbiomass %>% filter(variable == "Total Microbial Biomass")
var_importance(total.mbiomass)

mod.total.biomass <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                   mods = ~ pH.level + warming.duration + magnitude.of.warming +
                     ahm,
                   random = ~ 1 | studyID,
                   data = total.mbiomass)

summary(mod.total.biomass)

funnel(mod.total.biomass, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

##################################################################
# activity
##################################################################
mfunction <- read_excel("~/Desktop/MetaData.xlsx", sheet = 2)

mfunction <- mfunction %>% select(studyID, map, mat, ahm, ph, pH.level, ecosystem,
                                  Warming.technique, magnitude, magnitude.of.warming, 
                                  duration, warming.duration, deg.duration,
                                  LRR, var, weight, weight.adjusted, weighted.lnR, 
                                  variable)

mfunction[c("studyID", "pH.level", "ecosystem", "Warming.technique", "magnitude.of.warming",
            "warming.duration")] <-
  lapply(mfunction[c("studyID", "pH.level", "ecosystem", "Warming.technique", "magnitude.of.warming",
                     "warming.duration")], factor)
sapply(mfunction, class)


# respiration

total.respiration <- mfunction %>% filter(variable == "Total respiration")

var_importance(total.respiration)

mod.total.respiration <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                            mods = ~ Warming.technique + ahm + pH.level + magnitude.of.warming +
                              ecosystem + warming.duration,
                            random = ~ 1 | studyID,
                            data = total.respiration)

summary(mod.total.respiration)

funnel(mod.total.respiration, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)


# microbial respiration
microbial.respiration <- mfunction %>% filter(variable == "Microbial respiration")

var_importance(microbial.respiration)

mod.mrespiration <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                                mods = ~ magnitude.of.warming + ahm,
                                random = ~ 1 | studyID,
                                data = microbial.respiration)

summary(mod.mrespiration)

funnel(mod.mrespiration, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

# oxidases

oxidases <- mfunction %>% filter(variable == "Oxidases")

var_importance(oxidases)

# not enough sample size

# C-hydrolysis

c_hydrolysis <- mfunction %>% filter(variable == "C-hydrolysis")
var_importance(c_hydrolysis)

mod.c_hydrolysis <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                           mods = ~ pH.level + ecosystem + warming.duration +
                             Warming.technique,
                           random = ~ 1 | studyID,
                           data = c_hydrolysis)

summary(mod.c_hydrolysis)

funnel(mod.c_hydrolysis, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

# N-hydrolysis
n_hydrolysis <- mfunction %>% filter(variable == "N-hydrolysis")
var_importance(n_hydrolysis)

mod.n_hydrolysis <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                           mods = ~ ahm + Warming.technique + magnitude.of.warming +
                            pH.level + warming.duration,
                           random = ~ 1 | studyID,
                           data = n_hydrolysis)

summary(mod.n_hydrolysis)

funnel(mod.n_hydrolysis, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)



# without converting numerical predictors to categories

var_importance <- function (data) {
  cforest <- cforest(LRR ~ ahm + ph + ecosystem +
                       magnitude + duration + Warming.technique,
                     data = data,
                     controls = (cforest_control(ntree = 1000, replace = T)))
  
  varimp <- varimp(cforest)
  plot(varimp)
  ranked.importance <- as.data.frame(sort(varimp, decreasing = T))
  names(ranked.importance) <- "index"
  ranked.importance$predictor <- row.names(ranked.importance)
  ranked.importance$predictor <- factor(ranked.importance$predictor,
                                        levels = ranked.importance$predictor[order(ranked.importance$index)])
  ranked_plot <- ggplot(ranked.importance, aes(x = predictor, y = index)) +
    geom_bar(stat = "identity", fill = "black", col = "black") +
    coord_flip() +
    theme_bw()
  return(ranked_plot)
}

# mbc

mbc <- mbiomass %>% filter(variable == "MBC")
var_importance(mbc)


mod.mbc <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                  mods = ~ ph + duration + ahm + magnitude + ecosystem,
                  random = ~ 1 | studyID,
                  data = mbc)

summary(mod.mbc)

# ph had a significant effect
regplot(mod.mbc, mod = "ph")

funnel(mod.mbc, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)


# bacteria

bacteria <- mbiomass %>% filter(variable == "Bacteria")
var_importance(bacteria)


mod.bacteria <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                       mods = ~ duration + Warming.technique + magnitude + ahm,
                       random = ~ 1 | studyID,
                       data = bacteria)

summary(mod.bacteria)

funnel(mod.bacteria, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)


# fungi

fungi  <- mbiomass %>% filter(variable == "Fungi")
var_importance(fungi)


mod.fungi <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                    mods = ~ ph + duration + magnitude,
                    random = ~ 1 | studyID,
                    data = fungi)

summary(mod.fungi)

funnel(mod.fungi, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)


# Gram Positive bacteria

Gpos.bacteria <- mbiomass %>% filter(variable == "Gram Positive Bacteria")

var_importance(Gpos.bacteria)
# all predictors have negative variable importance



# Gram Negative bacteria

Gneg.bacteria <- mbiomass %>% filter(variable == "Gram Negative Bacteria")

var_importance(Gneg.bacteria)

mod.Gneg <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                   mods = ~ ahm + ecosystem + duration +
                     Warming.technique,
                   random = ~ 1 | studyID,
                   data = Gneg.bacteria)

summary(mod.Gneg)

funnel(mod.Gneg, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

# Total biomass

total.mbiomass <- mbiomass %>% filter(variable == "Total Microbial Biomass")
var_importance(total.mbiomass)

mod.total.biomass <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                            mods = ~ ph + magnitude + ahm,
                            random = ~ 1 | studyID,
                            data = total.mbiomass)

summary(mod.total.biomass)

funnel(mod.total.biomass, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

# respiration

total.respiration <- mfunction %>% filter(variable == "Total respiration")

var_importance(total.respiration)

mod.total.respiration <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                                mods = ~ magnitude + ahm + Warming.technique + ph +
                                  duration + ecosystem,
                                random = ~ 1 | studyID,
                                data = total.respiration)

summary(mod.total.respiration)
# magnitude and ph is significant
regplot(mod.total.respiration, mod = "ph")
regplot(mod.total.respiration, mod = "magnitude")


funnel(mod.total.respiration, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)


# microbial respiration
microbial.respiration <- mfunction %>% filter(variable == "Microbial respiration") %>% 
  filter(duration < 7)

var_importance(microbial.respiration)

mod.mrespiration <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                           mods = ~magnitude + ahm + ecosystem + ph + duration,
                           random = ~ 1 | studyID,
                           data = microbial.respiration)

summary(mod.mrespiration)

# magnitude and ahm is significant

regplot(mod.mrespiration, mod = "magnitude")
regplot(mod.mrespiration, mod = "ahm")


funnel(mod.mrespiration, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

# oxidases

oxidases <- mfunction %>% filter(variable == "Oxidases")

var_importance(oxidases)

# not enough sample size

# C-hydrolysis

c_hydrolysis <- mfunction %>% filter(variable == "C-hydrolysis")
var_importance(c_hydrolysis)

mod.c_hydrolysis <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                           mods = ~ ecosystem + ph + Warming.technique + ahm,
                           random = ~ 1 | studyID,
                           data = c_hydrolysis)

summary(mod.c_hydrolysis)

# ecosystem, pg, and ahm is significant

regplot(mod.c_hydrolysis, mod = "ph")
regplot(mod.c_hydrolysis, mod = "ahm")



funnel(mod.c_hydrolysis, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

# N-hydrolysis
n_hydrolysis <- mfunction %>% filter(variable == "N-hydrolysis")
var_importance(n_hydrolysis)

mod.n_hydrolysis <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                           mods = ~ magnitude + a hm + Warming.technique + ph,
                           random = ~ 1 | studyID,
                           data = n_hydrolysis)

summary(mod.n_hydrolysis)

funnel(mod.n_hydrolysis, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)



























