setwd('')

library(tidyverse)
library(vegan)
library(robust)
library(qvalue)
library(data.table)

# Outliers - Capblancq method
##---------------------------------------------------------------------------------------------------------------------
## Load Freqs
load("Outliers1.RData")

## Function
rdadapt <- function(rda, K) {
  
  loadings <- rda$CCA$v[, 1:as.numeric(K), drop=FALSE]
  resscale <- apply(loadings, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda <- stats::median(resmaha) / stats::qchisq(0.5, df = K)
  reschi2test <- stats::pchisq(resmaha / lambda, K, lower.tail = FALSE)
  qval <- qvalue(reschi2test)
  data.frame(p.values = reschi2test, q.values = qval$qvalues)
  
}

K <- 2
rdadapt_env <- rdadapt(RDAclim, K)

## Q-value with FDR
thres_q <- 0.01

outliers <- data.frame(
  
  SNP   = colnames(freq.df)[which(rdadapt_env$q.values < thres_q)],
  q.value = rdadapt_env$q.values[which(rdadapt_env$q.values < thres_q)]
  
)

# Output outliers with CHROM and POS
outliers.df <- merge(outliers, chrom_pos.df, by = 'SNP')
colnames(outliers.df) <- c('SNP', 'Q-value', 'Chromosome', 'Position')
write.table(outliers.df, file = 'RDAclim_outliersQ001.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

## Clean up
rm(freq.df, K, rdadapt_env, rdadapt, thres_q, outliers)

save.image("Outliers2.RData")
