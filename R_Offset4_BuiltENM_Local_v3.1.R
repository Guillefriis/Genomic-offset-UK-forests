setwd('')

library(robust)
library(tidyverse)

# Variable testing
##---------------------------------------------------------------------------------------------------------------------
## Load Clim and Freqs
load("ClimData1.RData")
pops_freq.df <- read.table("pops_freqdf.txt", header = F)

Env_now_set <- Env_now_all[pops_freq.df$V1, , drop=FALSE]
identical(pops_freq.df$V1, rownames(Env_now_set))

Env_now_set <- Env_now_set %>%
  as.matrix() %>%
  scale(center = TRUE, scale = TRUE)

# Redundancy analysis - Variable filtering with Ecological Model 4 (Lite)
##---------------------------------------------------------------------------------------------------------------------
## Correlation
Env_now_cov <- as.data.frame(Env_now_set[ , c("Bio1", "Bio5", "Bio6", "Bio14")])
psych::pairs.panels(Env_now_cov,
                    method = "pearson",
                    hist.col = "#00AFBB",
                    density = TRUE,
                    ellipses = TRUE,
                    lm = TRUE,
                    ci = TRUE)

cor_matrix <- cor(Env_now_cov, use = "pairwise.complete.obs", method = 'pearson')
hc <- caret::findCorrelation(cor_matrix, cutoff = 0.8)
colnames(Env_now_cov[, setdiff(seq_len(ncol(Env_now_cov)), hc), drop = FALSE])

# Plot selected variables
##------------------------------------------------------------------------------------------------------------------
library(terra)
library(tmap)
library(sf)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)

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

## Plot
library(tmap)

tmap_mode("plot")
cols <- c("Bio1", "Bio5", "Bio6", "Bio14")

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
