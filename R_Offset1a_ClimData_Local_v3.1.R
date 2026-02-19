setwd('')

library(readxl)
library(tidyverse)
library(terra)
library(rnaturalearth)
library(psych)
library(tmap)

# Metadata
##---------------------------------------------------------------------------------------------------------------------
path_to_excel <- 'Birch_CompleteDataset_AnalysisFilts.xlsx'
pop_all.df <- read_excel(path_to_excel, sheet = "DatasetAnalysis_RDA")

pop_all.df <- pop_all.df %>%
  select(
    Sample = Sample,
    Pop = Population,
    Type = Type,
    Latitude = Lat,
    Longitude = Lon,
    Platform = Platform
)

table(pop_all.df$Pop)
nlevels(as.factor(pop_all.df$Pop))

pop.df <- pop_all.df %>%
  add_count(Pop) %>%
  filter(n >= 9) %>%
  select(-n)

table(pop.df$Pop)
nlevels(as.factor(pop.df$Pop))

# Check platform
platform_check <- pop_all.df %>%
  count(Pop, Platform) %>%
  group_by(Pop) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(Pop, desc(freq))

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

tmap_mode("plot")

tm_shape(bion1_land) +
  tm_raster(style = "cont", palette = "-RdYlBu", title = "BIO1 (°C)", legend.reverse = TRUE) +
  tm_shape(land_poly) +
  tm_borders(lwd = 0.5, col = "black") +
  tm_shape(birch_coords) +
  tm_symbols(shape = 21, col = "darkblue", size = 0.4, border.col = "black", border.lwd = 0.5) +
  tm_layout(main.title = "Annual Mean Temperature (1981–2010)",
            legend.outside = TRUE,
            frame = FALSE)

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

tmap_mode("plot")

tm_shape(biof1_land) +
  tm_raster(style = "cont", palette = "-RdYlBu", title = "BIO1 (°C)", legend.reverse = TRUE) +
  tm_shape(land_poly) +
  tm_borders(lwd = 0.5, col = "black") +
  tm_shape(birch_coords) +
  tm_symbols(shape = 21, col = "darkblue", size = 0.4, border.col = "black", border.lwd = 0.5) +
  tm_layout(main.title = "Annual Mean Temperature (2041–2070)",
            legend.outside = TRUE,
            frame = FALSE)

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

tmap_mode("plot")

tm_shape(bior1_land) +
  tm_raster(style = "cont", palette = "-RdYlBu", title = "BIO1 (°C)", legend.reverse = TRUE) +
  tm_shape(land_poly) +
  tm_borders(lwd = 0.5, col = "black") +
  tm_shape(birch_coords) +
  tm_symbols(shape = 21, col = "darkblue", size = 0.4, border.col = "black", border.lwd = 0.5) +
  tm_layout(main.title = "Annual Mean Temperature (2071–2100)",
            legend.outside = TRUE,
            frame = FALSE)

rm(bio_f1, bio_fstack, bio_n1, bio_nstack, bio_r1, bio_rstack, biof1_land, bion1_land, bior1_land, env_fvals, env_nvals,
   env_rvals, land_rast, bio_far, bio_fut, bio_now, clim_far, clim_fut, clim_now, remove.NAs.stack, uk_extent, birch_coords)

# Validating data
##---------------------------------------------------------------------------------------------------------------------
## Current
Env_now_ind <- merge(pop.df, birch_env_now, by.x = 'Sample', by.y = 'Name')
Env_now_all <- Env_now_ind
row.names(Env_now_all) <- Env_now_all$Sample
Env_now_all <- Env_now_all[, -c(1, 3:5, 7)]

Env_now_all <- Env_now_all %>%
  group_by(Pop) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    PropILL = mean(Platform == "ILLUMINA", na.rm = TRUE)
  )

Env_now_all <- as.data.frame(Env_now_all)
row.names(Env_now_all) <- Env_now_all$Pop
Env_now_all <- Env_now_all[, -1]

# Future
Env_fut_ind <- merge(pop.df, birch_env_fut, by.x = 'Sample', by.y = 'Name')
Env_fut_all <- Env_fut_ind
row.names(Env_fut_all) <- Env_fut_all$Sample
Env_fut_all <- Env_fut_all[, -c(1, 3:5, 7)]

Env_fut_all <- Env_fut_all %>%
  group_by(Pop) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    PropILL = mean(Platform == "ILLUMINA", na.rm = TRUE)
  )

Env_fut_all <- as.data.frame(Env_fut_all)
row.names(Env_fut_all) <- Env_fut_all$Pop
Env_fut_all <- Env_fut_all[, -1]

# Far
Env_far_ind <- merge(pop.df, birch_env_far, by.x = 'Sample', by.y = 'Name')
Env_far_all <- Env_far_ind
row.names(Env_far_all) <- Env_far_all$Sample
Env_far_all <- Env_far_all[, -c(1, 3:5, 7)]

Env_far_all <- Env_far_all %>%
  group_by(Pop) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    PropILL = mean(Platform == "ILLUMINA", na.rm = TRUE)
  )

Env_far_all <- as.data.frame(Env_far_all)
row.names(Env_far_all) <- Env_far_all$Pop
Env_far_all <- Env_far_all[, -1]

rm(bio_fcrop, bio_rcrop, bio_ncrop) # These need to be rebuild downstreams

save.image('ClimData1b.RData')

