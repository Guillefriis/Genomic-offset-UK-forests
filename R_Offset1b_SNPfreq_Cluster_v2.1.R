setwd('')

library(vroom)

args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1] # chromosome name passed from the bash script (e.g., Bpe_Chr1)

cat("Chromosome being processed:", chrom, "\n")

# Load positions and subset by chromosome (global order matches .012 columns)
pos.vec <- read.table('birch_dataset.012.pos', header = FALSE)
chr_idx <- which(pos.vec$V1 == chrom)

# Load only the SNP columns for this chromosome from the .012 matrix
path_to_012 <- 'birch_dataset.012'
snps.dataset <- vroom(
  path_to_012,
  delim = '\t',
  col_names = FALSE,
  col_select = chr_idx + 1,
  col_types = cols(.default = col_integer())
)

# Individuals (row order for the .012)
inds <- read.table('birch_dataset.012.indv', stringsAsFactors = FALSE)$V1
row.names(snps.dataset) <- inds

# Sample → population
sample_pops.vec <- read.table('samplepops_birch.txt', sep = '\t', header = FALSE)
colnames(sample_pops.vec) <- c('Sample', 'Pop')
sample_pops.vec <- sample_pops.vec[match(inds, sample_pops.vec$Sample), ]
pop.vec <- as.factor(sample_pops.vec$Pop)

# Compute frequencies (ignore NAs) -- REPLACE aggregate() with pop-wise colMeans loop
cat("Computing SNP frequencies with pop-wise colMeans (NA-robust)...\n")

cn <- colnames(snps.dataset)

G <- as.matrix(snps.dataset)
rm(snps.dataset)
G[G == -1L] <- NA_integer_

pop_levels <- levels(pop.vec)
freqs.df <- matrix(NA_real_, nrow = length(pop_levels), ncol = ncol(G))
row.names(freqs.df) <- pop_levels
colnames(freqs.df) <- cn
rm(cn)

for (k in seq_along(pop_levels)) {
  rr <- pop.vec == pop_levels[k]
  freqs.df[k, ] <- colMeans(G[rr, , drop = FALSE], na.rm = TRUE) / 2
}

freqs.df <- as.data.frame(freqs.df)

# REMOVE sites invariant across populations BEFORE imputation
var_across_pops <- apply(freqs.df, 2, var, na.rm = TRUE)
keep_cols <- which(var_across_pops > 0)

dropped <- ncol(freqs.df) - length(keep_cols)
cat("Dropping", dropped, "invariant frequency sites across populations before imputation.\n")

if (dropped > 0) {
  freqs.df <- freqs.df[, keep_cols, drop = FALSE]
  pos_chr <- pos.vec[chr_idx, ][keep_cols, , drop = FALSE]

  outpos <- paste0('SNPfreq_variants_', chrom, '.012.pos')
  write.table(pos_chr, file = outpos, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Impute NAs per SNP with median across populations (must remove all NAs)
cat("Number of NAs before imputation:", sum(is.na(freqs.df)), "\n")

for (i in 1:ncol(freqs.df)){
  freqs.df[which(is.na(freqs.df[,i])), i] <- median(freqs.df[-which(is.na(freqs.df[,i])), i], na.rm = TRUE)
  if (i %% 10000 == 0) cat("Imputed SNP", i, "of", ncol(freqs.df), "\n")
}

cat("Number of NAs after imputation:", sum(is.na(freqs.df)), "\n")

# Write frequency table per chromosome
outname <- paste0("SNPfreq_", chrom, ".txt")
write.table(freqs.df, file = outname,
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# Cleanup
if (exists("pos_chr")) rm(pos_chr)
if (exists("outpos")) rm(outpos)
rm(G, k, pop_levels, rr, pos.vec, chr_idx, dropped, i, inds, keep_cols, outname, path_to_012, pop.vec, sample_pops.vec, var_across_pops)

# Save workspace
save.image(paste0('SNPfreq_', chrom, '.RData'))
cat('Chromosome', chrom, 'processing completed.\n')
