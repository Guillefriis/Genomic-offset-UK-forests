setwd('')

library(vegan)
library(tidyverse)

# Load data
##---------------------------------------------------------------------------------------------------------------------
## Load Clim and Freqs
load("Outliers0.RData")
load("ClimData1.RData")

Env_now_cov <- Env_now_all[row.names(freq.df), , drop=FALSE]
stopifnot(identical(rownames(freq.df), rownames(Env_now_cov)))

Env_now_cov <- Env_now_cov %>%
  as.matrix() %>%
  scale(center = TRUE, scale = TRUE)

Env_now_cov <- Env_now_cov[ , c("Bio5", "Bio6", "Bio8", "Bio14", "Bio15")]

# Redundancy analysis - Variable filtering with Ecological Model 4 (Lite)
##---------------------------------------------------------------------------------------------------------------------
# VIF
RDAclim <- rda(freq.df ~., data = Env_now_cov, scale = TRUE)

vifs <- vif.cca(RDAclim)
capture.output(vifs, file = "VIFS_LocalENM.txt")

rda_sum <- summary(RDAclim)
capture.output(rda_sum, file = "RDA_summary.txt")

anova.clim <- anova.cca(RDAclim, permutations = 999, parallel = 4)
capture.output(anova.clim, file = "Anova_LocalENM.txt")

gc()
save.image("Outliers1.RData")
