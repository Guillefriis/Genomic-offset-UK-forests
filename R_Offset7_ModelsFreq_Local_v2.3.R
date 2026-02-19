setwd('')

library(vroom)
library(dplyr)
library(readxl)
library(vegan)
library(terra)
library(tmap)
library(sf)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)

# SNP data
##---------------------------------------------------------------------------------------------------------------------
load("Outliers0.RData")
outliers_chrompos.df  <- vroom("BirchRDA_candsQ001.012.pos", col_names = c("CHR", "POS"))
outliers_chrompos.df  <- as.data.frame(outliers_chrompos.df)

str(chrom_pos.df)
str(outliers_chrompos.df)

idx_keep <- chrom_pos.df %>%
  mutate(idx = row_number()) %>%
  semi_join(outliers_chrompos.df, by = c("CHROM" = "CHR", "POS" = "POS")) %>%
  pull(idx)

rm(chrom_pos.df)
gc()

snp.df <- freq.df %>%
  select(all_of(idx_keep))

rm(freq.df, idx_keep)
gc()

colnames(snp.df) <- paste0("SNP", seq(1:ncol(snp.df)))
samps <- rownames(snp.df)

save.image("Outliers3.RData") ## Restart R, to actually free space
load("Outliers3.RData")

# Rebuilding terra objects
##---------------------------------------------------------------------------------------------------------------------
path_to_excel <- 'Birch_CompleteDataset_AnalysisFilt.xlsx'
pop_all.df <- read_excel(path_to_excel, sheet = "SamplesHF_filter_INDV20dip_drel")

pop_all.df <- pop_all.df %>%
  select(
    Sample = Sample,
    Pop = Population,
    Type = Type,
    Latitude = Lat,
    Longitude = Lon
  )

table(pop_all.df$Pop)
nlevels(as.factor(pop_all.df$Pop))

pop.df <- pop_all.df %>%
  add_count(Pop) %>%
  filter(n >= 9) %>%
  select(-n)

table(pop.df$Pop)
nlevels(as.factor(pop.df$Pop))

rm(path_to_excel)

# Load CHELSA current BIO climate layers (1981–2100)
##------------------------------------------------------------------------------------------------------------------
clim_now <- "E_CHELSA_Data/1981-2010/bio"
clim_fut <- "E_CHELSA_Data/2041-2070/MPI-ESM1-2-HR/ssp585/bio"
clim_far <- "E_CHELSA_Data/2071-2100/MPI-ESM1-2-HR/ssp585/bio"

uk_extent <- ext(-11, 2.5, 49.5, 61)
birch_coords <- vect(pop.df[, c(5, 4)], geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
land_poly <- ne_countries(scale = "medium", returnclass = "sf")

remove.NAs.stack <- function(rast.stack) {
  mask_layer <- sum(rast.stack, na.rm = TRUE)
  mask_layer[!is.na(mask_layer)] <- 1
  rast.stack <- mask(rast.stack, mask_layer)
  return(rast.stack)
}

## Current climate
bio_now <- list.files(clim_now, pattern = "bio.*1981-2010.*\\.tif$", full.names = TRUE)
bio_nstack <- rast(bio_now)
names(bio_nstack) <- c("Bio1", "Bio10", "Bio11", "Bio12", "Bio13",
                       "Bio14", "Bio15", "Bio16", "Bio17", "Bio18",
                       "Bio19", "Bio2", "Bio3", "Bio4", "Bio5",
                       "Bio6", "Bio7", "Bio8", "Bio9")

bio_ncrop <- crop(bio_nstack, uk_extent)
bio_ncrop <- remove.NAs.stack(bio_ncrop)
bio_n1 <- bio_ncrop[[1]]

land_rast <- rasterize(vect(land_poly), bio_n1, field = 1)
bion1_land <- mask(bio_n1, land_rast)

env_nvals <- extract(bio_ncrop, birch_coords)
birch_env_now <- cbind(Name = pop.df$Sample, env_nvals)

## Future climate
bio_fut <- list.files(clim_fut, pattern = "bio.*2041-2070.*\\.tif$", full.names = TRUE)
bio_fstack <- rast(bio_fut)
names(bio_fstack) <- c("Bio1", "Bio10", "Bio11", "Bio12", "Bio13",
                       "Bio14", "Bio15", "Bio16", "Bio17", "Bio18",
                       "Bio19", "Bio2", "Bio3", "Bio4", "Bio5",
                       "Bio6", "Bio7", "Bio8", "Bio9")


bio_fcrop <- crop(bio_fstack, uk_extent)
bio_fcrop <- remove.NAs.stack(bio_fcrop)
bio_f1 <- bio_fcrop[[1]]

land_rast <- rasterize(vect(land_poly), bio_f1, field = 1)
biof1_land <- mask(bio_f1, land_rast)

env_fvals <- extract(bio_fcrop, birch_coords)
birch_env_fut <- cbind(Name = pop.df$Sample, env_fvals)

## Far climate
bio_far <- list.files(clim_far, pattern = "bio.*2071-2100.*\\.tif$", full.names = TRUE)
bio_rstack <- rast(bio_far)
names(bio_rstack) <- c("Bio1", "Bio10", "Bio11", "Bio12", "Bio13",
                       "Bio14", "Bio15", "Bio16", "Bio17", "Bio18",
                       "Bio19", "Bio2", "Bio3", "Bio4", "Bio5",
                       "Bio6", "Bio7", "Bio8", "Bio9")


bio_rcrop <- crop(bio_rstack, uk_extent)
bio_rcrop <- remove.NAs.stack(bio_rcrop)
bio_r1 <- bio_rcrop[[1]]

land_rast <- rasterize(vect(land_poly), bio_r1, field = 1)
bior1_land <- mask(bio_r1, land_rast)

env_rvals <- extract(bio_rcrop, birch_coords)
birch_env_far <- cbind(Name = pop.df$Sample, env_rvals)

## Plot variables of interest
##---------------------------------------------------------------------------------------------------------------------
library(tmap)

tmap_mode("plot")
cols <- c("Bio5", "Bio6", "Bio8", "Bio14", "Bio15")

## UK polygon (same source as your example)
uk_sf <- ne_countries(country = "united kingdom",
                      scale = "medium",
                      returnclass = "sf")

uk_vect <- vect(uk_sf)

## mask function: crop + mask to UK borders
mask_to_uk <- function(r) {
  r <- crop(r, uk_vect)
  mask(r, uk_vect)
}

## apply to stacks
bio_n_uk <- mask_to_uk(bio_ncrop)
bio_f_uk <- mask_to_uk(bio_fcrop)
bio_r_uk <- mask_to_uk(bio_rcrop)

## plotting function: one figure per period, five maps
plot_period <- function(bio_stack, title_prefix) {
  
  maps <- lapply(cols, function(v) {
    tm_shape(bio_stack[[v]]) +
      tm_raster(
        title = v,
        palette = "viridis",
        style = "cont"
      ) +
      tm_layout(
        main.title = paste(title_prefix, v),
        legend.outside = TRUE
      )
  })
  
  tmap_arrange(maps, ncol = 3)
}

## produce figures
tmap_save(
  plot_period(bio_n_uk, "Current (1981–2010)"),
  filename = "Climate_Current_1981_2010.pdf",
  width = 10,
  height = 7,
  units = "in"
)

tmap_save(
  plot_period(bio_f_uk, "Future (2041–2070)"),
  filename = "Climate_Future_2041_2070.pdf",
  width = 10,
  height = 7,
  units = "in"
)

tmap_save(
  plot_period(bio_r_uk, "Far future (2071–2100)"),
  filename = "Climate_FarFuture_2071_2100.pdf",
  width = 10,
  height = 7,
  units = "in"
)

detach("package:tmap", unload = TRUE)

## Tidy up
rm(bio_f1, bio_fstack, bio_n1, bio_nstack, bio_r1, bio_rstack, biof1_land, bion1_land, bior1_land, env_fvals, env_nvals,
   env_rvals, land_rast, bio_far, bio_fut, bio_now, clim_far, clim_fut, clim_now, remove.NAs.stack, uk_extent, birch_coords,
   uk_sf, uk_vect, mask_to_uk, bio_n_uk, bio_f_uk, bio_r_uk, plot_period)
gc()

# Building Env datasets
##---------------------------------------------------------------------------------------------------------------------
# Current
Env_now_ind <- merge(pop.df, birch_env_now, by.x = 'Sample', by.y = 'Name')
Env_now_all <- Env_now_ind
row.names(Env_now_all) <- Env_now_all$Sample
Env_now_all <- Env_now_all[, -c(1, 3:6)]

Env_now_all <- Env_now_all %>%
  group_by(Pop) %>%
  summarise_all(mean, na.rm = TRUE)

Env_now_all <- as.data.frame(Env_now_all)
row.names(Env_now_all) <- Env_now_all$Pop
Env_now_all <- Env_now_all[, -1]

# Future
Env_fut_ind <- merge(pop.df, birch_env_fut, by.x = 'Sample', by.y = 'Name')
Env_fut_all <- Env_fut_ind
row.names(Env_fut_all) <- Env_fut_all$Sample
Env_fut_all <- Env_fut_all[, -c(1, 3:6)]

Env_fut_all <- Env_fut_all %>%
  group_by(Pop) %>%
  summarise_all(mean, na.rm = TRUE)

Env_fut_all <- as.data.frame(Env_fut_all)
row.names(Env_fut_all) <- Env_fut_all$Pop
Env_fut_all <- Env_fut_all[, -1]

# Far
Env_far_ind <- merge(pop.df, birch_env_far, by.x = 'Sample', by.y = 'Name')
Env_far_all <- Env_far_ind
row.names(Env_far_all) <- Env_far_all$Sample
Env_far_all <- Env_far_all[, -c(1, 3:6)]

Env_far_all <- Env_far_all %>%
  group_by(Pop) %>%
  summarise_all(mean, na.rm = TRUE)

Env_far_all <- as.data.frame(Env_far_all)
row.names(Env_far_all) <- Env_far_all$Pop
Env_far_all <- Env_far_all[, -1]

# Variables model
##---------------------------------------------------------------------------------------------------------------------
# scale NOW on a matrix, grab attrs, then convert to df
m_now <- Env_now_all[samps, cols, drop = FALSE]
m_now <- m_now %>%
  as.matrix() %>%
  scale(center = TRUE, scale = TRUE)

center_env <- attr(m_now, "scaled:center")
scale_env  <- attr(m_now, "scaled:scale")

Env_now <- m_now %>% as.data.frame()

# FUT and FAR use same cols and params, then convert to df
Env_fut <- Env_fut_all[samps, cols, drop = FALSE]
Env_fut <- Env_fut %>%
  select(all_of(cols)) %>%
  as.matrix() %>%
  scale(center = center_env[cols], scale = scale_env[cols]) %>%
  as.data.frame()

Env_far <- Env_far_all[samps, cols, drop = FALSE]
Env_far <- Env_far %>%
  select(all_of(cols)) %>%
  as.matrix() %>%
  scale(center = center_env[cols], scale = scale_env[cols]) %>%
  as.data.frame()

rm(m_now)
save.image('ClimData2.RData')

# Adaptive landscape
##---------------------------------------------------------------------------------------------------------------------
Env_now_cov <- Env_now
Env_fut_cov <- Env_fut
Env_far_cov <- Env_far

# pop info aligned to SNP order
pop_aligned <- pop.df %>%
  filter(Sample %in% samps) %>%
  slice(match(samps, Sample))

identical(row.names(snp.df), row.names(Env_now_cov))
identical(row.names(snp.df), row.names(Env_fut_cov))
identical(row.names(snp.df), row.names(Env_far_cov))
RDAout <- rda(snp.df ~., data = Env_now_cov)

RsquareAdj(RDAout)

# Plotting SNPs
TAB_loci <- as.data.frame(scores(RDAout, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDAout, choices=c(1:2), display="bp"))
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*3, y=RDA2*3), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))


rm(TAB_loci, TAB_var)

## Capblancq adaptive index, weighting:
adaptive_index2 <- function(RDA, K, env_pres, range = NULL, scale_env, center_env) {
  
  # Environmental variables used in RDA
  env_vars <- rownames(RDA$CCA$biplot)
  
  # Extract raster values + coords
  var_env_proj_pres <- as.data.frame(env_pres[[env_vars]], xy = TRUE)
  
  # Standardize environmental variables
  var_env_proj_RDA <- as.data.frame(scale(
    var_env_proj_pres[, env_vars],
    center = center_env[env_vars],
    scale = scale_env[env_vars]
  ))
  
  # Empty list for rasters
  Proj_pres <- list()
  
  # Generate rasters for each RDA axis
  for (i in 1:K) {
    values <- apply(var_env_proj_RDA[, env_vars], 1, function(x) sum(x * RDA$CCA$biplot[, i]))
    coords <- var_env_proj_pres[, c("x", "y")]
    raster_df <- cbind(coords, values)
    ras_pres <- terra::rast(raster_df, type = "xyz", crs = terra::crs(env_pres))
    names(ras_pres) <- paste0("RDA_pres_", i)
    Proj_pres[[paste0("RDA", i)]] <- ras_pres
  }
  
  # Optionally mask
  if (!is.null(range)) {
    Proj_pres <- lapply(Proj_pres, function(x) terra::mask(x, range))
  }
  
  # Compute combined weighted raster
  weights <- RDA$CCA$eig / sum(RDA$CCA$eig)
  weights <- weights[1:K]
  
  combined_vals <- Reduce(`+`, lapply(1:K, function(i) Proj_pres[[i]] * weights[i]))
  names(combined_vals) <- "RDA_weighted"
  
  # Add weighted raster to list
  Proj_pres[["RDA_weighted"]] <- combined_vals
  
  return(Proj_pres)
}


## Run the updated adaptive_index function
uk_shape <- land_poly[land_poly$name %in% "United Kingdom", ]
uk_range <- vect(uk_shape)

res_RDA_proj_current <- adaptive_index2(RDA = RDAout, K = 2, env_pres = bio_ncrop, range = uk_range,
                                        scale_env = scale_env, center_env = center_env)

# Convert each raster to data.frame, label it, and combine
plot_data_list <- lapply(names(res_RDA_proj_current), function(name) {
  df <- as.data.frame(res_RDA_proj_current[[name]], xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "value"
  df$value <- (df$value - min(df$value)) / (max(df$value) - min(df$value))  # Scale to [0–1]
  df$variable <- name
  return(df)
})

TAB_RDA <- do.call(rbind, plot_data_list)
TAB_RDA$variable <- factor(TAB_RDA$variable, levels = names(res_RDA_proj_current))
summary(TAB_RDA)

# Plot
ggplot(data = TAB_RDA) +
  geom_sf(data = land_poly, fill = "grey90", color = NA) +
  geom_raster(aes(x = x, y = y, fill = value), alpha = 0.8) +
  scale_fill_viridis_c(direction = -1, option = "A", guide = "none") +
  geom_sf(data = land_poly, fill = NA, color = "black", size = 0.2) +
  coord_sf(xlim = c(-11, 3), ylim = c(49.5, 61), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  facet_wrap(~ variable, ncol = 3) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size = 11)
  )

rm(RDA_proj, res_RDA_proj_current, plot_data_list, TAB_RDA)

# Genomic offset
##---------------------------------------------------------------------------------------------------------------------
# Function to predict genomic offset from a RDA model
genomic_offset <- function(RDA, K, env_pres, env_fut, range = NULL, scale_env, center_env){
  
  # Mask with the range if supplied
  if(!is.null(range)){
    env_pres <- mask(env_pres, range)
    env_fut <- mask(env_fut, range)
  }
  
  # Formatting and scaling environmental rasters for projection
  env_vars <- row.names(RDA$CCA$biplot)
  
  df_pres <- as.data.frame(env_pres[[env_vars]], xy = T, na.rm = TRUE)
  df_fut <- as.data.frame(env_fut[[env_vars]], xy = T, na.rm = TRUE)
  
  var_env_proj_pres <- as.data.frame(scale(df_pres[, env_vars], center = center_env[env_vars], scale = scale_env[env_vars]))
  var_env_proj_fut  <- as.data.frame(scale(df_fut[,  env_vars], center = center_env[env_vars], scale = scale_env[env_vars]))
  
  # Projection for each RDA axis
  Proj_pres <- list()
  Proj_fut <- list()
  Proj_offset <- list()
    
  for(i in 1:K){
      
      # Current climates
      ras_pres <- env_pres[[1]]
      ras_pres[!is.na(ras_pres)] <- as.vector(apply(var_env_proj_pres[,env_vars], 1, function(x) sum( x * RDA$CCA$biplot[,i])))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
      
      # Future climates
      ras_fut <- env_fut[[1]]
      ras_fut[!is.na(ras_fut)] <- as.vector(apply(var_env_proj_fut[,env_vars], 1, function(x) sum( x * RDA$CCA$biplot[,i])))
      Proj_fut[[i]] <- ras_fut
      names(ras_fut) <- paste0("RDA_fut_", as.character(i))
      names(Proj_fut)[i] <- paste0("RDA", as.character(i))
      
      # Single axis genetic offset 
      Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
      names(Proj_offset)[i] <- paste0("RDA", as.character(i))
      
    }
  
  # Weights based on axis eigenvalues
  weights <- RDA$CCA$eig/sum(RDA$CCA$eig)
  
  # Weighing the current and future adaptive indices based on the eigen values of the associated axes
  Proj_offset_pres <- do.call(cbind, lapply(1:K, function(x) as.data.frame(Proj_pres[[x]], xy = FALSE)))
  Proj_offset_pres <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_pres[,x]*weights[x])))
  Proj_offset_fut  <- do.call(cbind, lapply(1:K, function(x) as.data.frame(Proj_fut[[x]], xy = FALSE)))
  Proj_offset_fut  <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut[,x]*weights[x])))
  
  # Predict a global genetic offset, incorporating the K first axes weighted by their eigenvalues
  ras <- Proj_offset[[1]]
  ras[!is.na(ras)] <- unlist(lapply(1:nrow(Proj_offset_pres), function(x)
    dist(rbind(Proj_offset_pres[x,], Proj_offset_fut[x,]), method = "euclidean")))
  names(ras) <- "Global_offset"
  Proj_offset_global <- ras
  
  # Return projections for current and future climates for each RDA axis, prediction of genetic offset for each RDA axis and a global genetic offset
  
  return(list(
    Proj_pres = Proj_pres,
    Proj_fut = Proj_fut,
    Proj_offset = Proj_offset,
    Proj_offset_global = Proj_offset_global,
    weights = weights[1:K]
  ))
  
}

## Projections
res_RDA_proj2041 <- genomic_offset(RDAout, K = 2, env_pres = bio_ncrop, env_fut = bio_fcrop,
                                   range = uk_range, scale_env = scale_env, center_env = center_env)


res_RDA_proj2071 <- genomic_offset(RDAout, K = 2, env_pres = bio_ncrop, env_fut = bio_rcrop,
                                   range = uk_range, scale_env = scale_env, center_env = center_env)

# Convert offset rasters to data frames with coordinates
df_2041 <- as.data.frame(res_RDA_proj2041$Proj_offset_global, xy = TRUE)
df_2041$Date <- "2041-2070"
df_2071 <- as.data.frame(res_RDA_proj2071$Proj_offset_global, xy = TRUE)
df_2071$Date <- "2071-2100"

# Combine into one table
RDA_proj_offset <- rbind(df_2041, df_2071)

summary(RDA_proj_offset)
aggregate(Global_offset ~ Date, RDA_proj_offset, summary)

# Define color palette
colors <- c(
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[6:5])(2),
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[4:3])(2),
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[2:1])(3)
)

# Plot
ggplot(data = RDA_proj_offset) +
  geom_sf(data = land_poly[land_poly$name == "United Kingdom", ], fill = gray(.9), size = 0) +
  geom_raster(aes(x = x, y = y, fill = Global_offset)) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral")[6:1],
                       name = "Genomic offset") +
  geom_sf(data = land_poly[land_poly$name == "United Kingdom", ], fill = NA, size = 0.1) +
  coord_sf(xlim = c(-11, 2.5), ylim = c(49.5, 61), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~Date) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 11))

rm(df_2041, df_2071, RDA_df, RDA_proj_offset, res_RDA_proj2041, res_RDA_proj2071, colors)


# Climate-resilient assisted gene flow
##---------------------------------------------------------------------------------------------------------------------
adaptive_translocation <- function(RDA, K, env_pres, env_fut, target_coords, range = NULL, scale_env, center_env) {
  library(terra)
  
  # Extract variables used in RDA
  env_vars <- rownames(RDA$CCA$biplot)
  
  # Mask layers if range is provided
  if (!is.null(range)) {
    env_pres <- terra::mask(env_pres, range)
    env_fut  <- terra::mask(env_fut, range)
  }
  
  # Extract values and coordinates
  df_pres <- as.data.frame(env_pres[[env_vars]], xy = TRUE, na.rm = TRUE)
  df_fut  <- as.data.frame(env_fut[[env_vars]],  xy = TRUE, na.rm = TRUE)
  
  # Standardize
  var_pres <- scale(df_pres[, env_vars], center = center_env[env_vars], scale = scale_env[env_vars])
  var_fut  <- scale(df_fut[, env_vars],  center = center_env[env_vars], scale = scale_env[env_vars])
  
  # RDA space projection using loadings
  pres_proj <- sapply(1:K, function(i) apply(var_pres, 1, function(x) sum(x * RDA$CCA$biplot[, i])))
  fut_proj  <- sapply(1:K, function(i) apply(var_fut,  1, function(x) sum(x * RDA$CCA$biplot[, i])))
  
  # Identify the index of the target pixel in future projection
  fut_coords <- df_fut[, c("x", "y")]
  target_idx <- which.min((fut_coords$x - target_coords[1])^2 + (fut_coords$y - target_coords[2])^2)
  
  # Get the RDA space projection at the target site
  target_point <- fut_proj[target_idx, ]
  
  # Compute Euclidean distance from each present pixel to the target future RDA projection
  offset_vals <- apply(pres_proj, 1, function(row) {
    dist(rbind(row, target_point), method = "euclidean")
  })
  
  # Rebuild raster
  coords <- df_pres[, c("x", "y")]
  raster_df <- cbind(coords, offset_vals)
  ras_offset <- terra::rast(raster_df, type = "xyz", crs = terra::crs(env_pres))
  names(ras_offset) <- "Offset_to_Target"
  
  return(ras_offset)
}

#target_coords <- c(-3.64, 56.86)
target_coords <- c(-3.77, 50.62)

ras_offset_2041 <- adaptive_translocation(RDA = RDAout, K = 2, env_pres = bio_ncrop, env_fut = bio_fcrop,
                                   target_coords = target_coords, range = uk_range, scale_env = scale_env,
                                   center_env = center_env)

ras_offset_2071 <- adaptive_translocation(RDA = RDAout, K = 2, env_pres = bio_ncrop, env_fut = bio_rcrop,
                                          target_coords = target_coords, range = uk_range, scale_env = scale_env,
                                          center_env = center_env)



# Convert offset rasters to data frames with coordinates
df_2041 <- as.data.frame(ras_offset_2041, xy = TRUE)
df_2041$Date <- "2041–2070"

df_2071 <- as.data.frame(ras_offset_2071, xy = TRUE)
df_2071$Date <- "2071–2100"

# Combine into one data frame
df_offset <- rbind(df_2041, df_2071)

summary(df_offset)
aggregate(Offset_to_Target ~ Date, df_offset, summary)

# Plot using ggplot2
colors <- c(
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[6:5])(4),
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[4:3])(4),
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[2:1])(4)
)


ggplot(data = df_offset) + 
  geom_sf(data = land_poly[land_poly$name == "United Kingdom", ], fill = gray(.9), size = 0) +
  geom_raster(aes(x = x, y = y, fill = Offset_to_Target), alpha = 1) + 
  scale_fill_viridis_c(
    direction = -1,
    option = "A",
    name = paste0("Genomic offset\nto target site\n(", target_coords[2], ", ", target_coords[1], ")")
  ) +
  geom_point(data = data.frame(x = target_coords[1], y = target_coords[2]), 
             aes(x = x, y = y), 
             shape = 8, size = 3, color = "black", stroke = 1.2) +
  geom_sf(data = land_poly[land_poly$name == "United Kingdom", ], fill = NA, size = 0.1) +
  coord_sf(xlim = c(-11, 2.5), ylim = c(49.5, 61), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~Date) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size = 11)
  )

rm(df_2041, df_2071, RDA_df, ras_offset_2041, ras_offset_2071, colors, target_coords)

# Climate-resilient afforestation coverage
##---------------------------------------------------------------------------------------------------------------------
adaptive_coverage <- function(RDA, K, env_pres, env_fut, target_coords_list, range = NULL, scale_env, center_env) {
  
  # Extract variables used in RDA
  env_vars <- rownames(RDA$CCA$biplot)
  
  # Mask layers if range is provided
  if (!is.null(range)) {
    env_pres <- terra::mask(env_pres, range)
    env_fut  <- terra::mask(env_fut, range)
  }
  
  # Extract and standardize values and coordinates
  df_pres <- as.data.frame(env_pres[[env_vars]], xy = TRUE, na.rm = TRUE)
  df_fut  <- as.data.frame(env_fut[[env_vars]],  xy = TRUE, na.rm = TRUE)
  
  var_pres <- scale(df_pres[, env_vars], center = center_env[env_vars], scale = scale_env[env_vars])
  var_fut  <- scale(df_fut[, env_vars],  center = center_env[env_vars], scale = scale_env[env_vars])
  
  # RDA space projection using loadings
  pres_proj <- sapply(1:K, function(i) apply(var_pres, 1, function(x) sum(x * RDA$CCA$biplot[, i])))
  fut_proj  <- sapply(1:K, function(i) apply(var_fut,  1, function(x) sum(x * RDA$CCA$biplot[, i])))
  
  # Initialize raster stack
  offset_stack <- list()
  coords <- df_pres[, c("x", "y")]
  
  for (j in seq_along(target_coords_list)) {
    target_coords <- target_coords_list[[j]]
    
    # Identify the index of the target pixel in future projection
    fut_coords <- df_fut[, c("x", "y")]
    target_idx <- which.min((fut_coords$x - target_coords[1])^2 + (fut_coords$y - target_coords[2])^2)
    
    # Get the RDA projection of the target site
    target_point <- fut_proj[target_idx, ]
    
    # Compute Euclidean distance from each present pixel to the target projection
    offset_vals <- apply(pres_proj, 1, function(row) {
      dist(rbind(row, target_point), method = "euclidean")
    })
    
    # Create raster for each target point
    raster_df <- cbind(coords, offset_vals)
    ras_offset <- terra::rast(raster_df, type = "xyz", crs = terra::crs(env_pres))
    names(ras_offset) <- paste0("Offset_to_Target_", j)
    
    offset_stack[[j]] <- ras_offset
  }
  
  # Normalize each raster to 0–1 and find per-cell minimum
  norm_stack <- lapply(offset_stack, function(ras) {
    r_min <- terra::global(ras, fun = "min", na.rm = TRUE)[[1]]
    r_max <- terra::global(ras, fun = "max", na.rm = TRUE)[[1]]
    norm_ras <- (ras - r_min) / (r_max - r_min)
    return(norm_ras)
  })
  
  combined_offset <- do.call(terra::app, c(norm_stack, list(fun = min, na.rm = TRUE)))
  names(combined_offset) <- "Min_Normalised_Offset"
  
  return(combined_offset)
}


coords_plustrees <- read.table('registered_seedsources.txt', header = T, sep = '\t', dec = '.')
coords_plustrees <- coords_plustrees[, c('Lon', 'Lat')]

ras_coverage_2041 <- adaptive_coverage(RDA = RDAout, K = 2, env_pres = bio_ncrop, env_fut = bio_fcrop,
                                          target_coords_list = coords_plustrees, range = uk_range, scale_env = scale_env,
                                          center_env = center_env)

ras_coverage_2071 <- adaptive_coverage(RDA = RDAout, K = 2, env_pres = bio_ncrop, env_fut = bio_rcrop,
                                          target_coords_list = coords_plustrees, range = uk_range, scale_env = scale_env,
                                          center_env = center_env)

# Dataframes
df_cov_2041 <- as.data.frame(ras_coverage_2041, xy = TRUE)
df_cov_2041$Date <- "2041–2070"

df_cov_2071 <- as.data.frame(ras_coverage_2071, xy = TRUE)
df_cov_2071$Date <- "2071–2100"

df_cov_all <- rbind(df_cov_2041, df_cov_2071)

summary(df_cov_all)
aggregate(Min_Normalised_Offset ~ Date, df_cov_all, summary)

# Plot
ggplot(data = df_cov_all) +
  geom_sf(data = land_poly[land_poly$name == "United Kingdom", ], fill = gray(.9), size = 0) +
  geom_raster(aes(x = x, y = y, fill = Min_Normalised_Offset), alpha = 1) +
  scale_fill_viridis_c(
    direction = -1,
    option = "A",
    name = "Minimum normalised\ngenomic offset"
  ) +
  geom_point(data = coords_plustrees, aes(x = Lon, y = Lat),
             color = "black", fill = "red", size = 2, shape = 23, stroke = 0.5) +
  geom_sf(data = land_poly[land_poly$name == "United Kingdom", ], fill = NA, size = 0.1) +
  coord_sf(xlim = c(-11, 2.5), ylim = c(49.5, 61), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~Date) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size = 11)
  )

rm(coords_plustrees, df_cov_2041, df_cov_2071, df_cov_all, ras_coverage_2041, ras_coverage_2071)


# Maladaptation in saplings using observed vs predicted genotype counts (FIXED .012 parsing)
##---------------------------------------------------------------------------------------------------------------------
# Keep planted and colonisers
pop_placol.df <- pop_all.df %>%
  filter(Type %in% c("Pla", "Col"))

pop_placol.df <- pop_placol.df %>%
  group_by(Pop) %>%
  summarise(
    Type = first(Type),
    .groups = "drop"
  )  

# Predicted genotypes (genotype counts from RDA)
Env_now_sap <- Env_now_all[rownames(Env_now_all) %in% pop_placol.df$Pop, , drop = FALSE]

Env_now_sap <- Env_now_sap %>%
  select(all_of(cols)) %>%
  as.matrix() %>%
  scale(center = center_env[cols], scale = scale_env[cols]) %>%
  as.data.frame()

freqs_pred.df <- as.data.frame(predict(RDAout, newdata = Env_now_sap, type = "response"))

### SAP FREQS ###

path_to_012 <- 'BirchSAP_candsQ001.012'
snps.dataset <- vroom(
  path_to_012,
  delim = '\t',
  col_names = FALSE,
  col_types = cols(.default = col_integer())
)

snps.dataset <- as.data.frame(snps.dataset[,-1])

# Individuals (row order for the .012)
inds <- read.table('BirchSAP_candsQ001.012.indv', stringsAsFactors = FALSE)$V1
row.names(snps.dataset) <- inds

# Sample → population
sample_pops.vec <- read.table('samplepops20_birch_Sap.txt', sep = '\t', header = FALSE)
colnames(sample_pops.vec) <- c('Sample', 'Pop')
sample_pops.vec <- sample_pops.vec[match(inds, sample_pops.vec$Sample), ]
pop.vec <- as.factor(sample_pops.vec$Pop)

# Compute frequencies
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

freqs_sap.df <- as.data.frame(freqs.df)
colnames(freqs_sap.df) <- paste0("SNP", seq_len(ncol(freqs_sap.df)))
rm(G, k, pop_levels, rr, inds, path_to_012, pop.vec, sample_pops.vec, freqs.df)

## Sanity checks
stopifnot(
  rownames(freqs_sap.df) == rownames(freqs_pred.df),
  colnames(freqs_sap.df) == colnames(freqs_pred.df)
)

## Euclidean maladaptation per population
maladapt_euclid <- sapply(seq_len(nrow(freqs_pred.df)), function(i) {
  dist(rbind(freqs_sap.df[i, ], freqs_pred.df[i, ]))
})

results <- data.frame(
  Pop = rownames(freqs_pred.df),
  Type = tapply(pop_placol.df$Type, pop_placol.df$Pop, unique)[rownames(freqs_pred.df)],
  Maladapt = as.numeric(maladapt_euclid),
  row.names = NULL
)

results

ggplot(results, aes(x = Type, y = Maladapt, fill = Type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 2) +
  theme_bw() +
  labs(y = "Euclidean distance (obs vs pred)",
       x = "Population type")

# Cute plot
results$Type <- factor(results$Type, levels = c("Col", "Pla"))

results_sorted <- results %>%
  arrange(Type, Maladapt) %>%
  mutate(Pop = factor(Pop, levels = rev(Pop)))

ggplot(results_sorted, aes(x = Maladapt, y = Pop, fill = Type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(Maladapt, 3)), 
            hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = c("Col" = "#56B4E9", "Pla" = "#E69F00")) +
  theme_bw() +
  xlab("Maladaptation distance") +
  ylab("Population") +
  coord_cartesian(xlim = c(0, max(results$Maladapt) * 1.15)) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "top"
  )

## Significance permutation test: Observed statistic: difference in means (Pla – Col)
obs_stat <- with(results, mean(Maladapt[Type == "Pla"]) - mean(Maladapt[Type == "Col"]))
set.seed(1)

perm_stats <- replicate(10000, {
  perm_type <- sample(results$Type)
  mean(results$Maladapt[perm_type == "Pla"]) -
    mean(results$Maladapt[perm_type == "Col"])
})

# two-sided p-value
p_perm <- mean(abs(perm_stats) >= abs(obs_stat))
p_perm

p_to_stars <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

sig_label <- p_to_stars(p_perm)
sig_label

y_max <- max(results$Maladapt)

ggplot(results_sorted, aes(x = Maladapt, y = Pop, fill = Type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(Maladapt, 3)),
            hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = c("Col" = "#56B4E9", "Pla" = "#E69F00")) +
  theme_bw() +
  xlab("Maladaptation distance") +
  ylab("Population") +
  coord_cartesian(xlim = c(0, y_max * 1.25)) +
  annotate("text",
           x = y_max * 1.15,
           y = 3.5,
           label = paste0("Pla vs Col: ", sig_label, "\n(p = ",
                          signif(p_perm, 3), ")"),
           size = 4)

## CI whiskers with bootstrap
set.seed(1)

boot_maladapt <- function(obs_freqs, pred_freqs, B = 1000) {
  n_snp <- length(obs_freqs)
  replicate(B, {
    idx <- sample.int(n_snp, replace = TRUE)
    as.numeric(dist(rbind(obs_freqs[idx], pred_freqs[idx])))
  })
}

boot_df <- lapply(seq_len(nrow(freqs_sap.df)), function(i) {
  
  obs  <- as.numeric(freqs_sap.df[i, ])
  pred <- as.numeric(freqs_pred.df[i, ])
  
  b <- boot_maladapt(obs, pred)
  
  data.frame(
    Pop = rownames(freqs_sap.df)[i],
    Maladapt = as.numeric(dist(rbind(obs, pred))),
    LCI = quantile(b, 0.025),
    UCI = quantile(b, 0.975)
  )
}) %>%
  bind_rows() %>%
  left_join(results[, c("Pop", "Type")], by = "Pop") %>%
  arrange(Type, Maladapt) %>%
  mutate(Pop = factor(Pop, levels = rev(Pop)))

ggplot(boot_df, aes(x = Maladapt, y = Pop, fill = Type)) +
  geom_col(width = 0.7) +
  geom_errorbarh(
    aes(xmin = LCI, xmax = UCI),
    height = 0.25,
    linewidth = 0.35
  ) +
  geom_text(
    aes(x = UCI, label = round(Maladapt, 3)),
    hjust = -0.25,
    size = 3.5
  ) +
  scale_fill_manual(values = c("Col" = "#56B4E9", "Pla" = "#E69F00")) +
  theme_bw() +
  xlab("Maladaptation distance") +
  ylab("Population") +
  coord_cartesian(xlim = c(0, max(boot_df$UCI) * 1.3)) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "top"
  )
