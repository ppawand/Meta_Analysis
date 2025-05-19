################
# Library
#################
library(metafor)
library(readxl)
library(dplyr)
library(tidyr)

################
# function
################
choose.best.model <- function(x) {
  
  mod1 <- rma.mv(LRR, LRR.var,
                 mods = ~ magnitude.of.warming + warming.duration,
                 random = ~ 1 | site.id,
                 data = x,
                 method = "ML", 
                 control = list(verbose = TRUE, optimizer="optim",
                                optmethod="Nelder-Mead"))
  mod2 <- rma.mv(LRR, LRR.var,
                 mods = ~ magnitude.of.warming + warming.duration + ecosystem,
                 random = ~ 1 | site.id,
                 data = x,
                 method = "ML", 
                 control = list(verbose = TRUE, optimizer="optim",
                                optmethod="Nelder-Mead"))
  mod3 <- rma.mv(LRR, LRR.var,
                 mods = ~ magnitude.of.warming + warming.duration + ecosystem + 
                   irrigation.status,
                 random = ~ 1 | site.id,
                 data = x,
                 method = "ML", 
                 control = list(verbose = TRUE, optimizer="optim",
                                optmethod="Nelder-Mead"))
  mod4 <- rma.mv(LRR, LRR.var,
                 mods = ~ magnitude.of.warming + warming.duration + ecosystem + 
                   irrigation.status + fertilization,
                 random = ~ 1 | site.id,
                 data = x,
                 method = "ML", 
                 control = list(verbose = TRUE, optimizer="optim",
                                optmethod="Nelder-Mead"))
  mod5 <- rma.mv(LRR, LRR.var,
                 mods = ~ magnitude.of.warming + warming.duration + ecosystem + 
                   irrigation.status + fertilization + pH,
                 random = ~ 1 | site.id,
                 data = x,
                 method = "ML", 
                 control = list(verbose = TRUE, optimizer="optim",
                                optmethod="Nelder-Mead"))
  mod6 <- rma.mv(LRR, LRR.var,
                 mods = ~ magnitude.of.warming + warming.duration + ecosystem + 
                   irrigation.status + fertilization + pH + map,
                 random = ~ 1 | site.id,
                 data = x,
                 method = "ML", 
                 control = list(verbose = TRUE, optimizer="optim",
                                optmethod="Nelder-Mead"))
  mod7 <- rma.mv(LRR, LRR.var,
                 mods = ~ magnitude.of.warming + warming.duration + ecosystem + 
                   irrigation.status + fertilization + pH + map + mat,
                 random = ~ 1 | site.id,
                 data = x,
                 method = "ML", 
                 control = list(verbose = TRUE, optimizer="optim",
                                optmethod="Nelder-Mead"))
  mod8 <- rma.mv(LRR, LRR.var,
                 mods = ~ magnitude.of.warming + ecosystem,
                 random = ~ 1 | site.id,
                 data = x,
                 method = "ML", 
                 control = list(verbose = TRUE, optimizer="optim",
                                optmethod="Nelder-Mead"))
  mod9 <- rma.mv(LRR, LRR.var,
                 mods = ~ magnitude.of.warming + warming.duration + ecosystem + pH,
                 random = ~ 1 | site.id,
                 data = x,
                 method = "ML", 
                 control = list(verbose = TRUE, optimizer="optim",
                                optmethod="Nelder-Mead"))
  aic <- AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9)
  best.model <- aic[which.min(aic$AIC), ]
  return(list(aic, best.model))
}

############
# Data
##############
mfunction <- read_excel("~/Desktop/Meta-analysis/MetaData.xlsx", sheet = 4)

mfunction <- mfunction %>%
  select(site.id, mat, map, pH, fertilization, irrigation.status, ecosystem, 
         magnitude.of.warming, warming.duration, control.mean, control.sd,
         control.n, trt.mean, trt.sd, trt.n, variable)
mfunction <- escalc(measure = "ROM",
                    m1i = trt.mean,
                    m2i = control.mean,
                    sd1i = trt.sd,
                    sd2i = control.sd,
                    n1i = trt.n,
                    n2i = control.n,
                    data = mfunction,
                    var.names = c("LRR", "LRR.var"))


na.pos <- is.na(mfunction$LRR.var)
mfunction$LRR.var[na.pos] <- 
  with(mfunction, 1 / ((trt.n[na.pos] * control.n[na.pos]) /  (trt.n[na.pos] + control.n[na.pos])))
mfunction <- mfunction %>% select(-control.mean, -control.sd, -control.n, -trt.mean, -trt.sd, -trt.n)

mfunction[c("fertilization", "irrigation.status", "ecosystem")] <-
  lapply(mfunction[c("fertilization", "irrigation.status", "ecosystem")], factor)
sapply(mfunction, class)

##################
# Soil Respiration
##################

total.respiration <- mfunction %>% filter(variable == "Total respiration")
  

choose.best.model(total.respiration)

best.model.total.respiration <- rma.mv(LRR, LRR.var,
                                   mods = ~ magnitude.of.warming + warming.duration +
                                     ecosystem + irrigation.status + fertilization +
                                     pH + map,
                                   random = ~ 1 | site.id,
                                   data = total.respiration,
                                   method = "ML", 
                                   control = list(verbose = TRUE, optimizer="optim",
                                                  optmethod="Nelder-Mead"))

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

choose.best.model(microbial.respiration)

best.model.microbial.respiration <- rma.mv(LRR, LRR.var,
                                       mods = ~ magnitude.of.warming + warming.duration +
                                         ecosystem + irrigation.status,
                                       random = ~ 1 | site.id,
                                       data = microbial.respiration,
                                       method = "ML", 
                                       control = list(verbose = TRUE, optimizer="optim",
                                                      optmethod="Nelder-Mead"))

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

choose.best.model(beta.glucosidase)

best.model.beta.glucosidase <- rma.mv(LRR, LRR.var,
                                       mods = ~ magnitude.of.warming + warming.duration +
                                         ecosystem + irrigation.status + fertilization,
                                       random = ~ 1 | site.id,
                                       data = beta.glucosidase,
                                       method = "ML", 
                                       control = list(verbose = TRUE, optimizer="optim",
                                                      optmethod="Nelder-Mead"))

summary(best.model.beta.glucosidase)

results.beta.glucosidase <- data.frame(coeff = best.model.beta.glucosidase$b, 
                                        lower.CI = best.model.beta.glucosidase$ci.lb, 
                                        upper.CI = best.model.beta.glucosidase$ci.ub,
                                        p.value = best.model.beta.glucosidase$pval)
funnel(best.model.beta.glucosidase,xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)
## Cellobiohydrolase

cellobiohydrolase <- mfunction %>% filter(variable == "CBH")

choose.best.model(cellobiohydrolase)

best.model.cellobiohydrolase <- rma.mv(LRR, LRR.var,
                                      mods = ~ magnitude.of.warming + warming.duration +
                                        ecosystem + irrigation.status + fertilization +
                                        pH,
                                      random = ~ 1 | site.id,
                                      data = cellobiohydrolase,
                                      method = "ML", 
                                      control = list(verbose = TRUE, optimizer="optim",
                                                     optmethod="Nelder-Mead"))

summary(best.model.cellobiohydrolase)

results.cellobiohyrolase <- data.frame(coeff = best.model.cellobiohydrolase$b, 
                                       lower.CI = best.model.cellobiohydrolase$ci.lb, 
                                       upper.CI = best.model.cellobiohydrolase$ci.ub,
                                       p.value = best.model.cellobiohydrolase$pval)

funnel(best.model.cellobiohydrolase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)
## N- acetyl-beta- Glucosaminidase

nag <- mfunction %>% filter(variable == "NAG")

choose.best.model(nag)

best.model.nag <- rma.mv(LRR, LRR.var,
                         mods = ~ magnitude.of.warming + warming.duration +
                         ecosystem + irrigation.status + fertilization +
                                         pH + map + mat,
                         random = ~ 1 | site.id,
                         data = nag,
                         method = "ML", 
                         control = list(verbose = TRUE, optimizer="optim",
                                                      optmethod="Nelder-Mead"))

summary(best.model.nag)

results.nag <- data.frame(coeff = best.model.nag$b, 
                          lower.CI = best.model.nag$ci.lb, 
                          upper.CI = best.model.nag$ci.ub,
                          p.value = best.model.nag$pval)
funnel(best.model.nag, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)

## beta-Xylosidase

beta.xylosidase <- mfunction %>% filter(variable == "BX")

choose.best.model(beta.xylosidase)

best.model.beta.xylosidase <- rma.mv(LRR, LRR.var,
                         mods = ~ magnitude.of.warming + warming.duration +
                           ecosystem,
                         random = ~ 1 | site.id,
                         data = beta.xylosidase,
                         method = "ML", 
                         control = list(verbose = TRUE, optimizer="optim",
                                        optmethod="Nelder-Mead"))

summary(best.model.beta.xylosidase)

results.beta.xylosidase<- data.frame(coeff = best.model.beta.xylosidase$b, 
                          lower.CI = best.model.beta.xylosidase$ci.lb, 
                          upper.CI = best.model.beta.xylosidase$ci.ub,
                          p.value = best.model.beta.xylosidase$pval)
funnel(best.model.beta.xylosidase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)

## Leucine aminopeptidase

lap <- mfunction %>% filter(variable == "LAP")

choose.best.model(lap)

best.model.lap <- rma.mv(LRR, LRR.var,
                         mods = ~ magnitude.of.warming + warming.duration +
                           ecosystem,
                         random = ~ 1 | site.id,
                         data = lap,
                         method = "ML", 
                         control = list(verbose = TRUE, optimizer="optim",
                                        optmethod="Nelder-Mead"))

summary(best.model.lap)

results.lap <- data.frame(coeff = best.model.lap$b, 
                          lower.CI = best.model.lap$ci.lb, 
                          upper.CI = best.model.lap$ci.ub,
                          p.value = best.model.lap$pval)
funnel(best.model.lap, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)


## Phenol Oxidase

phenol.oxidase <- mfunction %>% filter(variable == "PHO")

choose.best.model(phenol.oxidase)

best.model.phenol.oxidase <- rma.mv(LRR, LRR.var,
                         mods = ~ magnitude.of.warming + ecosystem,
                         random = ~ 1 | site.id,
                         data = phenol.oxidase,
                         method = "ML", 
                         control = list(verbose = TRUE, optimizer="optim",
                                        optmethod="Nelder-Mead"))

summary(best.model.phenol.oxidase)

results.phenol.oxidase <- data.frame(coeff = best.model.phenol.oxidase$b, 
                          lower.CI = best.model.phenol.oxidase$ci.lb, 
                          upper.CI = best.model.phenol.oxidase$ci.ub,
                          p.value = best.model.phenol.oxidase$pval)


funnel(best.model.phenol.oxidase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)

## perooxidase

perooxidase <- mfunction %>% filter(variable == "PO")
choose.best.model(perooxidase)

best.model.perooxidase <- rma.mv(LRR, LRR.var,
                         mods = ~ magnitude.of.warming + warming.duration,
                         random = ~ 1 | site.id,
                         data = perooxidase,
                         method = "ML", 
                         control = list(verbose = TRUE, optimizer="optim",
                                        optmethod="Nelder-Mead"))

summary(best.model.perooxidase)

results.perooxidase <- data.frame(coeff = best.model.perooxidase$b, 
                          lower.CI = best.model.perooxidase$ci.lb, 
                          upper.CI = best.model.perooxidase$ci.ub,
                          p.value = best.model.perooxidase$pval)
funnel(best.model.perooxidase, xlim = c(-4, 4), ylim = c(0, 1.8), back = "white", pch = 1)

###################
result.resp <- rbind(results.total.respiration, results.microbial.respiration,
                     results.beta.glucosidase, results.cellobiohyrolase,
                     results.nag, results.beta.xylosidase, results.lap,
                     results.phenol.oxidase, results.perooxidase)
write.csv(result.resp, "mfunction.results.csv")
