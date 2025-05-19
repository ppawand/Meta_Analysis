#################
# library
##################
library(ggplot2)
library(metafor)
library(dplyr)
library(tidyr)
library(readxl)
library(party)

################
# function
################
# rank the predictor based on their relative predictive power

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

# Summarize a data by a provide categorical variable and calculate bootstrapped 
# mean and confidence interval for each category

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



##################
# data
##################
mbiomass <- read_excel("~/Desktop/MetaData.xlsx", sheet = 1)
mbiomass <- mbiomass %>% select(studyID, trt, map, mat, ahm, ph, ecosystem,
                                Warming.technique, magnitude, duration, deg.duration,
                                LRR, var, weight.adjusted, weighted.lnR, 
                                variable)

mbiomass[c("studyID", "ecosystem", "Warming.technique")] <-
  lapply(mbiomass[c("studyID","ecosystem", "Warming.technique")], factor)
sapply(mbiomass, class)

#######################
# variable importance
########################
var_importance(mbiomass)
################
# MBC
################

mbc <- mbiomass %>% filter(variable == "MBC")

mod.mbc <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                mods = ~ ecosystem + ph + ahm + Warming.technique +
                  magnitude + duration,
                random = ~ 1 | studyID,
                data = mbc,
                method = "ML")

summary(mod.mbc)

funnel(mod.mbc, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)



##################
# Bacterial Biomass
#################
bacteria <- mbiomass %>% filter(variable == "Bacteria")

mod.bacteria <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                              mods = ~ ecosystem + ph + ahm + Warming.technique +
                                magnitude + duration,
                              random = ~ 1 | studyID,
                              data = bacteria,
                              method = "ML")
summary(mod.bacteria)

mean_bacteria <- bootstrap_mean_ci(data = bacteria, cat_var = "Warming.technique", num_var = "LRR")
mean_bacteria %>% mutate(percent_change = ((exp(mean) - 1) * 100))

funnel(best.model.bacteria, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

###########################
# Fungal Biomass
###########################
fungi <- mbiomass %>% filter(variable == "Fungi")


best.model.fungi <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                           mods = ~ ecosystem + ph + ahm + Warming.technique +
                             magnitude + duration,
                           random = ~ 1 | studyID,
                           data = fungi,
                           method = "ML")

summary(best.model.fungi)


results.fungi <- data.frame(coeff = best.model.fungi$b, 
                               lower.CI = best.model.fungi$ci.lb, 
                               upper.CI = best.model.fungi$ci.ub,
                               p.value = best.model.fungi$pval)
funnel(best.model.fungi, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

###########################
# Gram positive bacteria
###########################
Gpos.bacteria <- mbiomass %>% filter(variable == "Gram Positive Bacteria")


best.model.Gpos <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                          mods = ~ ecosystem + ph + ahm + Warming.technique +
                            magnitude + duration,
                          random = ~ 1 | studyID,
                          data = Gpos.bacteria,
                          method = "ML")

summary(best.model.Gpos)

results.Gpos <- data.frame(coeff = best.model.Gpos$b, 
                            lower.CI = best.model.Gpos$ci.lb, 
                            upper.CI = best.model.Gpos$ci.ub,
                            p.value = best.model.Gpos$pval)
funnel(best.model.Gpos, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

###########################
# Gram negative bacteria
###########################
Gneg.bacteria <- mbiomass %>% filter(variable == "Gram Negative Bacteria")


best.model.Gneg <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                          mods = ~ ecosystem + ph + ahm + Warming.technique +
                            magnitude + duration,
                          random = ~ 1 | studyID,
                          data = Gneg.bacteria,
                          method = "ML")

summary(best.model.Gneg)

results.Gneg <- data.frame(coeff = best.model.Gneg$b, 
                           lower.CI = best.model.Gneg$ci.lb, 
                           upper.CI = best.model.Gneg$ci.ub,
                           p.value = best.model.Gneg$pval)
funnel(best.model.Gneg, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

#######################
# total biomass
#####################
total.mbiomass <- mbiomass %>% filter(variable == "Total Microbial Biomass")


best.model.total.biomass <- rma.mv(yi = LRR, V = var, W = weight.adjusted,
                                   mods = ~ ecosystem + ph + ahm + Warming.technique +
                                     magnitude + duration,
                                   random = ~ 1 | studyID,
                                   data = total.mbiomass,
                                   method = "ML")

summary(best.model.total.biomass)


results.total.biomass <- data.frame(coeff = best.model.total.biomass$b, 
                           lower.CI = best.model.total.biomass$ci.lb, 
                           upper.CI = best.model.total.biomass$ci.ub,
                           p.value = best.model.total.biomass$pval)

funnel(best.model.total.biomass, xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

##################################################################
# activity
##################################################################
mfunction <- read_excel("~/Desktop/MetaData.xlsx", sheet = 2)

mfunction <- mfunction %>% select(studyID, trt, map, mat, ahm, ph, ecosystem,
                                Warming.technique, magnitude, duration, deg.duration,
                                LRR, var, weight.adjusted, 
                                variable)

mfunction[c("studyID", "ecosystem", "Warming.technique")] <-
  lapply(mfunction[c("studyID","ecosystem", "Warming.technique")], factor)


############################
# variable importance
###########################
var_importance(mfunction)

##################
# Soil Respiration
##################

total.respiration <- mfunction %>% filter(variable == "Total respiration")


best.model.total.respiration <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
                                       mods = ~ ecosystem + ph + ahm + Warming.technique +
                                         magnitude + duration,
                                       random = ~ 1 | studyID,
                                       data = total.respiration,
                                       method = "ML")

summary(best.model.total.respiration)

results.total.respiration <- data.frame(coeff = best.model.total.respiration$b, 
                                        lower.CI = best.model.total.respiration$ci.lb, 
                                        upper.CI = best.model.total.respiration$ci.ub,
                                        p.value = best.model.total.respiration$pval)
funnel(best.model.total.respiration,xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)
############################
# Soil Microbial Respiration
############################

microbial.respiration <- mfunction %>% filter(variable == "Microbial respiration")


best.model.microbial.respiration <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
                                           mods = ~ ecosystem + ph + ahm + Warming.technique +
                                             magnitude + duration,
                                           random = ~ 1 | studyID,
                                           data = microbial.respiration,
                                           method = "ML")

summary(best.model.microbial.respiration)

results.microbial.respiration <- data.frame(coeff = best.model.microbial.respiration$b, 
                                            lower.CI = best.model.microbial.respiration$ci.lb, 
                                            upper.CI = best.model.microbial.respiration$ci.ub,
                                            p.value = best.model.microbial.respiration$pval)

funnel(best.model.microbial.respiration,xlim = c(-3, 3), ylim = c(0, 1.5), back = "white", pch = 1)

###########################
# Enzyme activity
#######################
## beta Gulcosidase

beta.glucosidase <- mfunction %>% filter(variable == "BG")


best.model.beta.glucosidase <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
                                      mods = ~ ecosystem + ph + ahm + Warming.technique +
                                        magnitude + duration,
                                      random = ~ 1 | studyID,
                                      data = beta.glucosidase,
                                      method = "ML")

summary(best.model.beta.glucosidase)

results.beta.glucosidase <- data.frame(coeff = best.model.beta.glucosidase$b, 
                                       lower.CI = best.model.beta.glucosidase$ci.lb, 
                                       upper.CI = best.model.beta.glucosidase$ci.ub,
                                       p.value = best.model.beta.glucosidase$pval)
funnel(best.model.beta.glucosidase,xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)
## Cellobiohydrolase

cellobiohydrolase <- mfunction %>% filter(variable == "CBH")


best.model.cellobiohydrolase <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
                                       mods = ~ ecosystem + ph + ahm + Warming.technique +
                                         magnitude + duration,
                                       random = ~ 1 | studyID,
                                       data = cellobiohydrolase,
                                       method = "ML")

summary(best.model.cellobiohydrolase)

results.cellobiohyrolase <- data.frame(coeff = best.model.cellobiohydrolase$b, 
                                       lower.CI = best.model.cellobiohydrolase$ci.lb, 
                                       upper.CI = best.model.cellobiohydrolase$ci.ub,
                                       p.value = best.model.cellobiohydrolase$pval)

funnel(best.model.cellobiohydrolase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)

## N- acetyl-beta- Glucosaminidase

nag <- mfunction %>% filter(variable == "NAG")


best.model.nag <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
                         mods = ~ ecosystem + ph + ahm + Warming.technique +
                           magnitude + duration,
                         random = ~ 1 | studyID,
                         data = nag,
                         method = "ML")
summary(best.model.nag)

results.nag <- data.frame(coeff = best.model.nag$b, 
                          lower.CI = best.model.nag$ci.lb, 
                          upper.CI = best.model.nag$ci.ub,
                          p.value = best.model.nag$pval)
funnel(best.model.nag, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)

## beta-Xylosidase

beta.xylosidase <- mfunction %>% filter(variable == "BX")


best.model.beta.xylosidase <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
                                     mods = ~ ecosystem + ph + ahm + Warming.technique +
                                       magnitude + duration,
                                     random = ~ 1 | studyID,
                                     data = beta.xylosidase,
                                     method = "ML")

summary(best.model.beta.xylosidase)

results.beta.xylosidase<- data.frame(coeff = best.model.beta.xylosidase$b, 
                                     lower.CI = best.model.beta.xylosidase$ci.lb, 
                                     upper.CI = best.model.beta.xylosidase$ci.ub,
                                     p.value = best.model.beta.xylosidase$pval)
funnel(best.model.beta.xylosidase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)

## Leucine aminopeptidase

lap <- mfunction %>% filter(variable == "LAP")


best.model.lap <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
                         mods = ~ ecosystem + ph + ahm + Warming.technique +
                           magnitude + duration,
                         random = ~ 1 | studyID,
                         data = lap,
                         method = "ML")

summary(best.model.lap)

results.lap <- data.frame(coeff = best.model.lap$b, 
                          lower.CI = best.model.lap$ci.lb, 
                          upper.CI = best.model.lap$ci.ub,
                          p.value = best.model.lap$pval)
funnel(best.model.lap, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)


## Phenol Oxidase

phenol.oxidase <- mfunction %>% filter(variable == "PHO")


best.model.phenol.oxidase <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
                                    mods = ~ ecosystem + ph + ahm + Warming.technique +
                                      magnitude + duration,
                                    random = ~ 1 | studyID,
                                    data = phenol.oxidase,
                                    method = "ML")

summary(best.model.phenol.oxidase)

results.phenol.oxidase <- data.frame(coeff = best.model.phenol.oxidase$b, 
                                     lower.CI = best.model.phenol.oxidase$ci.lb, 
                                     upper.CI = best.model.phenol.oxidase$ci.ub,
                                     p.value = best.model.phenol.oxidase$pval)


funnel(best.model.phenol.oxidase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)

## perooxidase

perooxidase <- mfunction %>% filter(variable == "PO")

best.model.perooxidase <- rma.mv(yi = LRR, V = var, W = weight.adjusted, 
                                 mods = ~ ecosystem + ph + ahm + Warming.technique +
                                   magnitude + duration,
                                 random = ~ 1 | studyID,
                                 data = perooxidase,
                                 method = "ML")

summary(best.model.perooxidase)

results.perooxidase <- data.frame(coeff = best.model.perooxidase$b, 
                                  lower.CI = best.model.perooxidase$ci.lb, 
                                  upper.CI = best.model.perooxidase$ci.ub,
                                  p.value = best.model.perooxidase$pval)
funnel(best.model.perooxidase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)




















