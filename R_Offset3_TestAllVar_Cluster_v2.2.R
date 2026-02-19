setwd('')

library(vegan)
library(tidyverse)

## Redundancy analysis
bio_var <- commandArgs(trailingOnly = TRUE)[1]

load("Outliers0.RData")
load("ClimData1.RData")

Env_now_cov <- Env_now_all[row.names(freq.df), bio_var, drop=FALSE]
stopifnot(identical(rownames(freq.df), rownames(Env_now_cov)))

Env_now_cov <- Env_now_cov %>%
  as.matrix() %>%
  scale(center = TRUE, scale = TRUE)

## Variable screening
rda_single <- rda(freq.df ~ ., data = Env_now_cov, scale = TRUE)
res <- anova.cca(rda_single, by = "term", parallel = 2, permutations = 999)

writeLines(capture.output(res), paste0("anova_", bio_var, "_999.txt"))


