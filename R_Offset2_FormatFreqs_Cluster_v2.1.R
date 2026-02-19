setwd('')

library(vroom)
library(tibble)

## Freqs
freq.df <- vroom("SNPfreq_All.012", delim = "\t",
                 col_names = FALSE, progress = TRUE, num_threads = 4)

freq.df <- freq.df %>%
  rename(Pop = X1) %>%
  column_to_rownames("Pop") %>%
  as.data.frame()

colnames(freq.df) <- paste0("SNP", seq_len(ncol(freq.df)))

# POS table
chrom_pos.df <- vroom("SNPfreq_All.012.pos", delim = "\t",
                      col_names = c("CHROM", "POS"),
                      progress = TRUE, num_threads = 2)

chrom_pos.df <- chrom_pos.df %>%
  mutate(SNP = paste0("SNP", row_number())) %>%
  select(SNP, CHROM, POS) %>%
  as.data.frame()

save.image("Outliers0.RData")
