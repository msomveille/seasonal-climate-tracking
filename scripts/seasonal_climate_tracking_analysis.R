library(tidyverse)
library(rworldmap)
library(sf)
library(raster)
library(terra)
library(fields)
library(ggpubr)
library(ebirdst)
library(tidyterra)
library(MASS)
library(pracma)
library(ggnewscale)
library(gridExtra)

## Load data for sampled individuals
## This is the metadata obtained from the previously conducted genoscape analyses. It has 5 columns: species, season, population, longitude, and latitude. Each row is an individual sample.
data_for_analysis <- read_csv("resources/data_for_analysis_final.csv") %>%
  filter(!(population %in% c("KER", "SCA", "WM", "SouthDakota"))) #remove populations with less than 3 samples during at least one season
species.names <- unique(data_for_analysis$species)
data_for_analysis$longitude[which(data_for_analysis$longitude == 104.36700)] <- -104.36700 # correct a mistake

## Map of the Americas
newmap <- getMap(resolution = "low")
newmap <- spTransform(newmap, '+proj=longlat +datum=WGS84')
newmap@data$world <- rep(1,length(newmap@data$SOVEREIGNT))
newmap <- vect(newmap)
newmap <- terra::buffer(newmap, 0)
newmap <- terra::aggregate(newmap)
newmap <- terra::crop(newmap, ext(-180,-30,-60,90))

## Load ecoregion polygons (obtained from Dinerstein et al. 2017 Biosciences — freely available)
ecoregions <- terra::vect("resources/Ecoregions2017/Ecoregions2017.shp")
crs(ecoregions) <- "+proj=longlat +datum=WGS84"
ecoregions <- terra::crop(ecoregions, newmap)
ecoregions <- ecoregions[-96,] # remove an ecoregion whose polygon is causing problems with spatial manipulations

## Download ebird seasonal relative abundance surfaces for all species (freely available)
species.names.ebird <- c("wilfly", "yelwar", "wlswar", "comyel", "amered", "paibun", "herthr")
sp_path <- "resources/ebird_st/2022"
for(i in 1:length(species.names.ebird)){
  ebirdst::ebirdst_download_status(species.names.ebird[i], path="resources/ebird_st")
}

## Calculate the species relative abundance in each ecoregion and for each species
ecoregions_abund_presence <- list()
for(i in 1:length(species.names)){
  data_for_analysis2 <- data_for_analysis %>% filter(species == species.names[i])
  pops <- unique(data_for_analysis2$population)
  data_for_analysis2_W <- data_for_analysis2 %>% filter(season == "Wintering")
  data_for_analysis2_B <- data_for_analysis2 %>% filter(season == "Breeding")

  # Load eBird seasonal abundance surfaces
  abunds <- terra::rast(paste0(sp_path, "/", species.names.ebird[i], "/seasonal/", species.names.ebird[i], "_abundance_seasonal_mean_9km_2022.tif"))
  abunds <- terra::project(abunds, '+proj=longlat +datum=WGS84', method = "ngb")
  abunds[is.na(abunds)] <- 0
  abunds_B <- abunds[[1]]
  abunds_W <- abunds[[2]]
  
  # Calculate relative abundance in each ecoregion
  ecoregions_abund <- terra::extract(abunds_W, ecoregions, sum) %>%
    left_join(terra::extract(abunds_B, ecoregions, sum))

  # For each population, get occupied ecoregions 
  points_W_ecoregions <- points_B_ecoregions <- list()
  pops_name <- vector()
  for(j in 1:length(pops)){
    # Requires at least 3 distinct locations per season
    if(length(which(data_for_analysis2_W$population == pops[j])) > 0 & length(which(data_for_analysis2_B$population == pops[j])) > 0){
      points_W <- SpatialPoints(data_for_analysis2_W[which(data_for_analysis2_W$population == pops[j]), 4:5])
      points_B <- SpatialPoints(data_for_analysis2_B[which(data_for_analysis2_B$population == pops[j]), 4:5])
      crs(points_W) <- crs(points_B) <- "+proj=longlat +datum=WGS84"
      points_W_ecoregions[[j]] <- apply(relate(ecoregions, terra::vect(points_W), "contains"), 1, sum)
      points_B_ecoregions[[j]] <- apply(relate(ecoregions, terra::vect(points_B), "contains"), 1, sum)
      pops_name <- c(pops_name, pops[j])
    }
  }
  ecoregions_abund_presence[[i]] <- cbind(ecoregions_abund, do.call(cbind, points_W_ecoregions), do.call(cbind, points_B_ecoregions))
  colnames(ecoregions_abund_presence[[i]]) <- c("ID", "wintering_abund", "breeding_abund", paste0(pops_name, "_wintering_points"), paste0(pops_name, "_breeding_points"))
}

# write abundance distributions as csv files
seasonalAbundances_B <- data.frame(WIFL = ecoregions_abund_presence[[1]]$breeding_abund,
                                   YEWA = ecoregions_abund_presence[[2]]$breeding_abund,
                                   WIWA = ecoregions_abund_presence[[3]]$breeding_abund,
                                   COYE = ecoregions_abund_presence[[4]]$breeding_abund,
                                   AMRE = ecoregions_abund_presence[[5]]$breeding_abund,
                                   PABU = ecoregions_abund_presence[[6]]$breeding_abund,
                                   HETH = ecoregions_abund_presence[[7]]$breeding_abund)
seasonalAbundances_W <- data.frame(WIFL = ecoregions_abund_presence[[1]]$wintering_abund,
                                   YEWA = ecoregions_abund_presence[[2]]$wintering_abund,
                                   WIWA = ecoregions_abund_presence[[3]]$wintering_abund,
                                   COYE = ecoregions_abund_presence[[4]]$wintering_abund,
                                   AMRE = ecoregions_abund_presence[[5]]$wintering_abund,
                                   PABU = ecoregions_abund_presence[[6]]$wintering_abund,
                                   HETH = ecoregions_abund_presence[[7]]$wintering_abund)
seasonalAbundances_B <- apply(seasonalAbundances_B, 2, function(x) ifelse(is.na(x)==T, 0, x))
seasonalAbundances_W <- apply(seasonalAbundances_W, 2, function(x) ifelse(is.na(x)==T, 0, x))
seasonalAbundances_B <- apply(seasonalAbundances_B, 2, function(x) x/sum(x))
seasonalAbundances_W <- apply(seasonalAbundances_W, 2, function(x) x/sum(x))
write.csv(seasonalAbundances_B, "results/seasonalAbundances_B.csv", row.names = F)
write.csv(seasonalAbundances_W, "results/seasonalAbundances_W.csv", row.names = F)


## Fig S8: map of ecoregions with samples
# Ecoregions occupied by the species used in the study
occupied_ecoregions <- which(apply(do.call(cbind, lapply(ecoregions_abund_presence, function(x) x$wintering_abund + x$breeding_abund)), 1, sum) > 0)
pdf("results/Figure_S8.pdf", width = 12, height = 11)
ggplot() +
  geom_sf(data=newmap, col="grey60", fill="grey60") +
  geom_spatvector(data=ecoregions[occupied_ecoregions,], aes(fill=as.character(ECO_NAME)), col="grey20", lwd=0.05) +
  theme_void() + theme(legend.position="none") +
  ylim(c(-15, 90))
dev.off()


## Extract climate for each ecoregion. Climate was dowloaded from Chelsa (freely available)
years <- 2000:2018
months <- c("01","02","03","04","05","06","07","08","09","10","11","12")
months_days <- c("001", "032", "060", "091", "121", "152", "182", "213", "244", "274", "305", "335")
Temp_files <- list.files("resources/chelsa/chelsa_V2/GLOBAL/monthly/tas")
Temp_files <- Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[4])) %in% years]
Prec_files <- list.files("resources/chelsa/chelsa_V2/GLOBAL/monthly/pr")
Prec_files <- Prec_files[unlist(lapply(strsplit(Prec_files, "_"), function(x) x[4])) %in% years]
summer <- months[6:7]
summer_days <- months_days[6:7]
winter <- months[c(1:2, 12)]
winter_days <- months_days[c(1:2, 12)]

temp_BR <- matrix(Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[3])) %in% summer], ncol=length(summer))
temp_BR <- terra::rast(apply(temp_BR, 1, function(x) terra::rast(paste0("resources/chelsa/chelsa_V2/GLOBAL/monthly/tas/", x)))) # load temperature rasters
temp_BR <- (mean(temp_BR) / 10) - 273.15
prec_BR <- matrix(Prec_files[unlist(lapply(strsplit(Prec_files, "_"), function(x) x[3])) %in% summer], ncol=length(summer))
prec_BR <- terra::rast(apply(prec_BR, 1, function(x) terra::rast(paste0("resources/chelsa/chelsa_V2/GLOBAL/monthly/pr/", x)))) # load precipitation rasters
prec_BR <- mean(prec_BR)
temp_NB <- matrix(Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[3])) %in% winter], ncol=length(winter))
temp_NB <- terra::rast(apply(temp_NB, 1, function(x) terra::rast(paste0("resources/chelsa/chelsa_V2/GLOBAL/monthly/tas/", x)))) # load temperature rasters
temp_NB <- (mean(temp_NB) / 10) - 273.15
prec_NB <- matrix(Prec_files[unlist(lapply(strsplit(Prec_files, "_"), function(x) x[3])) %in% winter], ncol=length(winter))
prec_NB <- terra::rast(apply(prec_NB, 1, function(x) terra::rast(paste0("resources/chelsa/chelsa_V2/GLOBAL/monthly/pr/", x)))) # load precipitation rasters
prec_NB <- mean(prec_NB)
crs(temp_BR) <- crs(prec_BR) <- crs(temp_NB) <- crs(prec_NB) <- "+proj=longlat +datum=WGS84"

# crop around the Americas
temp_BR <- terra::crop(temp_BR, newmap, mask=T)
temp_NB <- terra::crop(temp_NB, newmap, mask=T)
prec_BR <- terra::crop(prec_BR, newmap, mask=T)
prec_NB <- terra::crop(prec_NB, newmap, mask=T)

# zscores
temp_mean <- mean(c(as.vector(temp_BR$mean), as.vector(temp_NB$mean)), na.rm=T)
temp_sd <- sd(c(as.vector(temp_BR$mean), as.vector(temp_NB$mean)), na.rm=T)
prec_mean <- mean(c(as.vector(prec_BR$mean), as.vector(prec_NB$mean)), na.rm=T)
prec_sd <- sd(c(as.vector(prec_BR$mean), as.vector(prec_NB$mean)), na.rm=T)
temp_zscore_BR <- (temp_BR - temp_mean) / temp_sd
prec_zscore_BR <- (prec_BR - prec_mean) / prec_sd
temp_zscore_NB <- (temp_NB - temp_mean) / temp_sd
prec_zscore_NB <- (prec_NB - prec_mean) / prec_sd

# extract climate for all ecoregions
ecoregions_climate_NB <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions, fun="mean", na.rm=T)
ecoregions_climate_BR <- terra::extract(c(temp_zscore_BR, prec_zscore_BR), ecoregions, fun="mean", na.rm=T)
colnames(ecoregions_climate_NB) <- colnames(ecoregions_climate_BR) <- c("ID", "temperature", "precipitation")

# write distance matrix as a csv file
distanceMat <- as.matrix(st_distance(st_as_sf(ecoregions)) / 1000)
write.csv(distanceMat, "results/distanceMatrix.csv", row.names = F, col.names = F)

# Run ORSIM using orsim.py

# calculate and write climate distance matrices as csv files
climate_distanceMat <- rdist(as.matrix(ecoregions_climate_BR[,2:3]), as.matrix(ecoregions_climate_NB[,2:3]))
thermal_distanceMat <- abs(outer(ecoregions_climate_BR$temperature, ecoregions_climate_NB$temperature, "-"))
precipitation_distanceMat <- abs(outer(ecoregions_climate_BR$precipitation, ecoregions_climate_NB$precipitation, "-"))

##  Calculate population variables used in the analysis (easonal climate overlap and migration distance)
niche_size <- niche_size_orsim <- niche_size_climate <- niche_size_thermal <- niche_size_precipitation <- list()
niche_overlap <- niche_overlap_orsim <- niche_overlap_climate <- niche_overlap_thermal <- niche_overlap_precipitation <- migration_distance_climate <- migration_distance_thermal <- migration_distance_precipitation <- migration_distance_orsim <- migration_distance <- list()
niche_overlap_temperature <- niche_overlap_precipitation <- niche_overlap_temperature_orsim <- niche_overlap_precipitation_orsim <- niche_overlap_temperature_climate <- niche_overlap_precipitation_climate <- niche_overlap_temperature_thermal <- niche_overlap_precipitation_thermal <- niche_overlap_temperature_precipitation <- niche_overlap_precipitation_precipitation <- population_names <- list()
niche_overlap_null <- niche_overlap_temperature_null <- niche_overlap_precipitation_null <- pop_centroid_breeding <- pop_centroid_wintering <- pop_centroid_breeding_orsim <- pop_centroid_wintering_orsim <- pop_centroid_wintering_climate <- pop_centroid_wintering_thermal <- pop_centroid_wintering_precipitation <- list()
for(k in 1:length(species.names)){
  # Population names
  pops_seas <- do.call(rbind, strsplit(colnames(ecoregions_abund_presence[[k]])[4:ncol(ecoregions_abund_presence[[k]])], "_"))
  pops_names <- unique(pops_seas[,1])
  population_names[[k]] <- pops_names
  
  # Which ecoregions are occupied by the species
  breeding_ecoregions <- which(ecoregions_abund_presence[[k]]$breeding_abund > 0)
  wintering_ecoregions <- which(ecoregions_abund_presence[[k]]$wintering_abund > 0)
  
  # Percentage of individuals of each population in each ecoregion
  ecoregions_abund_presence_W <- ecoregions_abund_presence[[k]][,4:(3+length(pops_names))]
  ss <- apply(ecoregions_abund_presence_W, 1, sum)
  ss <- ifelse(ss==0, 1, ss)
  ecoregions_abund_presence_W <- apply(ecoregions_abund_presence_W, 2, function(x) x/ss)
  ecoregions_abund_presence_B <- ecoregions_abund_presence[[k]][,(4+length(pops_names)):ncol(ecoregions_abund_presence[[k]])]
  ss <- apply(ecoregions_abund_presence_B, 1, sum)
  ss <- sum(ecoregions_abund_presence_B)
  ss <- ifelse(ss==0, 1, ss)
  ecoregions_abund_presence_B <- apply(ecoregions_abund_presence_B, 2, function(x) x/ss)
  
  # load ORSIM results for the species
  ORSIM_results <- read.csv(paste0("results/ORSIM_results_", species.names[k], ".csv"), header=F)
  
  # climate distances for the focal species
  climate_distanceMat_2 <- climate_distanceMat[breeding_ecoregions, wintering_ecoregions]
  thermal_distanceMat_2 <- thermal_distanceMat[breeding_ecoregions, wintering_ecoregions]
  precipitation_distanceMat_2 <- precipitation_distanceMat[breeding_ecoregions, wintering_ecoregions]
  
  # geographical centroids of wintering ecoregions
  centroid_wintering_ecoregions_spp <- geom(terra::centroids(ecoregions[wintering_ecoregions,]))[,3:4]
  
  # Calculate climate niches (empirical, simulated from ORSIM, simulated from seasonal climate tracking)
  niche_climate_wintering <- niche_temp_wintering <- niche_prec_wintering <- niche_climate_breeding <- niche_temp_breeding <- niche_prec_breeding <- list()
  niche_climate_wintering_orsim <- niche_temp_wintering_orsim <- niche_prec_wintering_orsim <- niche_climate_breeding_orsim <- niche_temp_breeding_orsim <- niche_prec_breeding_orsim <- list()
  niche_climate_wintering_climate <- niche_temp_wintering_climate <- niche_prec_wintering_climate <- niche_climate_breeding_climate <- niche_temp_breeding_climate <- niche_prec_breeding_climate <- list()
  niche_climate_wintering_thermal <- niche_temp_wintering_thermal <- niche_prec_wintering_thermal <- niche_climate_breeding_thermal <- niche_temp_breeding_thermal <- niche_prec_breeding_thermal <- list()
  niche_climate_wintering_precipitation <- niche_temp_wintering_precipitation <- niche_prec_wintering_precipitation <- niche_climate_breeding_precipitation <- niche_temp_breeding_precipitation <- niche_prec_breeding_precipitation <- pop_resampling <- list()
  niche_climate_wintering_null <- niche_temp_wintering_null <- niche_prec_wintering_null <- list()
  centroid_wintering_ecoregions_orsim <- centroid_breeding_ecoregions_orsim <- centroid_wintering_ecoregions <- centroid_breeding_ecoregions <- centroid_wintering_ecoregions_climate <- centroid_wintering_ecoregions_thermal <- centroid_wintering_ecoregions_precipitation <- vector()
  for(j in 4:(3+length(unique(pops_seas[,1])))){
    # which ecoregions are seasonally occupied by a given population
    breeding_pop_ecoregions <- which(ecoregions_abund_presence[[k]][breeding_ecoregions,][,j+length(pops_names)] > 0)
    wintering_pop_ecoregions <- which(ecoregions_abund_presence[[k]][wintering_ecoregions,][,j] > 0)
    
    # which wintering ecoregions are predicted to be occupied by ORSIM
    wintering_pop_ecoregions_orsim <- which(apply(ORSIM_results[breeding_pop_ecoregions,], 2, sum) > 0)
    
    # which wintering ecoregions are predicted to be occupied by seasonal climate tracking
    if(length(breeding_pop_ecoregions) > 1){
      wintering_pop_ecoregions_climate <- apply(climate_distanceMat_2[breeding_pop_ecoregions,], 1, which.min)
      wintering_pop_ecoregions_thermal <- apply(thermal_distanceMat_2[breeding_pop_ecoregions,], 1, which.min)
      wintering_pop_ecoregions_precipitation <- apply(precipitation_distanceMat_2[breeding_pop_ecoregions,], 1, which.min)
    }else{
      wintering_pop_ecoregions_climate <- which.min(climate_distanceMat_2[breeding_pop_ecoregions,])
      wintering_pop_ecoregions_thermal <- which.min(thermal_distanceMat_2[breeding_pop_ecoregions,])
      wintering_pop_ecoregions_precipitation <- which.min(precipitation_distanceMat_2[breeding_pop_ecoregions,])
    }
    
    # Ecoregions weights, based on population relative abundance and the fraction of individuals of the population in each occupied ecoregion 
    ecoregions_weights_wintering <- ecoregions_abund_presence[[k]]$wintering_abund[wintering_ecoregions][wintering_pop_ecoregions] * ecoregions_abund_presence_W[,j-3][wintering_ecoregions][wintering_pop_ecoregions]
    ecoregions_weights_wintering <- ecoregions_weights_wintering / sum(ecoregions_weights_wintering)
    ecoregions_weights_breeding <- ecoregions_abund_presence[[k]]$breeding_abund[breeding_ecoregions][breeding_pop_ecoregions] * ecoregions_abund_presence_B[,j-3][breeding_ecoregions][breeding_pop_ecoregions]
    ecoregions_weights_breeding <- ecoregions_weights_breeding / sum(ecoregions_weights_breeding)
    if(length(wintering_pop_ecoregions_orsim) > 1){
      ecoregions_weights_wintering_orsim <- apply(ORSIM_results[breeding_pop_ecoregions,][,wintering_pop_ecoregions_orsim], 2, sum)
      ecoregions_weights_breeding_orsim <- apply(ORSIM_results[breeding_pop_ecoregions,][,wintering_pop_ecoregions_orsim], 1, sum)
    }else{
      ecoregions_weights_wintering_orsim <- ORSIM_results[breeding_pop_ecoregions,][,wintering_pop_ecoregions_orsim]
      ecoregions_weights_breeding_orsim <- ORSIM_results[breeding_pop_ecoregions,][,wintering_pop_ecoregions_orsim]
    }
    ecoregions_weights_wintering_orsim <- ecoregions_weights_wintering_orsim / sum(ecoregions_weights_wintering_orsim)
    ecoregions_weights_breeding_orsim <- ecoregions_weights_breeding_orsim / sum(ecoregions_weights_breeding_orsim)
    
    # centroids of seasonally occupied ecoregions
    if(length(ecoregions[wintering_ecoregions,][wintering_pop_ecoregions,]) > 1){
      centroid_wintering_ecoregions <- rbind(centroid_wintering_ecoregions, 
                                             apply(matrix(geom(terra::centroids(ecoregions[wintering_ecoregions,][wintering_pop_ecoregions,]))[,3:4], ncol=2), 2, function(x) weighted.mean(x, w=ecoregions_weights_wintering))
      )
    }else{
      centroid_wintering_ecoregions <- rbind(centroid_wintering_ecoregions, 
                                             matrix(geom(terra::centroids(ecoregions[wintering_ecoregions,][wintering_pop_ecoregions,]))[,3:4], ncol=2)
      )
    }
    if(length(ecoregions[breeding_ecoregions,][breeding_pop_ecoregions,]) > 1){
      centroid_breeding_ecoregions <- rbind(centroid_breeding_ecoregions, 
                                            apply(matrix(geom(terra::centroids(ecoregions[breeding_ecoregions,][breeding_pop_ecoregions,]))[,3:4], ncol=2), 2, function(x) weighted.mean(x, w=ecoregions_weights_breeding))
      )
    }else{
      centroid_breeding_ecoregions <- rbind(centroid_breeding_ecoregions, 
                                            matrix(geom(terra::centroids(ecoregions[breeding_ecoregions,][breeding_pop_ecoregions,]))[,3:4], ncol=2)
      )
    }
    
    # centroids of seasonally occupied ecoregions from ORSIM
    if(length(ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_orsim,]) > 1){
      centroid_wintering_ecoregions_orsim <- rbind(centroid_wintering_ecoregions_orsim, 
                                                   apply(matrix(geom(terra::centroids(ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_orsim,]))[,3:4], ncol=2), 2, function(x) weighted.mean(x, w=ecoregions_weights_wintering_orsim))
      )
    }else{
      centroid_wintering_ecoregions_orsim <- rbind(centroid_wintering_ecoregions_orsim, 
                                                   matrix(geom(terra::centroids(ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_orsim,]))[,3:4], ncol=2)
      )
    }
    
    if(length(ecoregions[breeding_ecoregions,][breeding_pop_ecoregions,]) > 1){
      centroid_breeding_ecoregions_orsim <- rbind(centroid_breeding_ecoregions_orsim, 
                                                  apply(matrix(geom(terra::centroids(ecoregions[breeding_ecoregions,][breeding_pop_ecoregions,]))[,3:4], ncol=2), 2, function(x) weighted.mean(x, w=ecoregions_weights_breeding_orsim))
      )
    }else{
      centroid_breeding_ecoregions_orsim <- rbind(centroid_breeding_ecoregions_orsim, 
                                                  matrix(geom(terra::centroids(ecoregions[breeding_ecoregions,][breeding_pop_ecoregions,]))[,3:4], ncol=2)
      )
    }
    
    # centroids of seasonally occupied ecoregions from the model based on seasonal climate tracking
    centroid_wintering_ecoregions_climate <- rbind(centroid_wintering_ecoregions_climate, 
                                                   apply(matrix(geom(terra::centroids(ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_climate,]))[,3:4], ncol=2), 2, mean))
    centroid_wintering_ecoregions_thermal <- rbind(centroid_wintering_ecoregions_thermal, 
                                                   apply(matrix(geom(terra::centroids(ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_thermal,]))[,3:4], ncol=2), 2, mean))
    centroid_wintering_ecoregions_precipitation <- rbind(centroid_wintering_ecoregions_precipitation, 
                                                         apply(matrix(geom(terra::centroids(ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_precipitation,]))[,3:4], ncol=2), 2, mean))
    
    # null model resampling wintering ecoregions around orsim prediction
    # mean centroid of the wintering ecoregions predicted by ORSIM
    centroid_wintering_ecoregions_pop_orsim <- geom(terra::centroids(ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_orsim,]))[,3:4]
    if(length(wintering_pop_ecoregions_orsim) > 1){
      centroid_wintering_ecoregions_pop_orsim_mean <- apply(centroid_wintering_ecoregions_pop_orsim, 2, function(x) weighted.mean(x, w=ecoregions_weights_wintering_orsim))
    }else{
      centroid_wintering_ecoregions_pop_orsim_mean <- centroid_wintering_ecoregions_pop_orsim
    }
    # distance between mean centroid of the wintering ecoregions predicted by ORSIM and empirical wintering ecoregions
    ctrs_W_orsim_dists <- rdist.earth(matrix(centroid_wintering_ecoregions_pop_orsim_mean, ncol=2), centroid_wintering_ecoregions_spp, miles=F)[1,]
    if(length(wintering_pop_ecoregions_orsim) > 1){
      max_dist <- max(rdist.earth(centroid_wintering_ecoregions_pop_orsim, miles=F))
    }else{
      max_dist <- 1500
    }
    ecoreg_subset <- which(ctrs_W_orsim_dists < (max_dist/2))
    # Function to randomly sample ecoregions 
    sampling_fct <- function(x){
      wintering_ecoregions[sample(ecoreg_subset, length(wintering_pop_ecoregions_orsim), replace = F)]
    }
    pop_resampling <- t(sapply(rep(j,1000), function(x) sampling_fct(max_dist)))
    if(length(wintering_pop_ecoregions_orsim) == 1){ pop_resampling <- t(pop_resampling) }
    
    # Extract pixels of climate data from climate rasters in each ecoregion
    ecoregions_climate_wintering_orsim <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_orsim,])
    ecoregions_climate_wintering <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_ecoregions,][wintering_pop_ecoregions,])
    ecoregions_climate_breeding <- terra::extract(c(temp_zscore_BR, prec_zscore_BR), ecoregions[breeding_ecoregions,][breeding_pop_ecoregions,])
    colnames(ecoregions_climate_wintering_orsim) <- colnames(ecoregions_climate_wintering) <- colnames(ecoregions_climate_breeding) <- c("ID", "temp", "prec")
    
    ecoregions_climate_climate <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_climate,])
    ecoregions_climate_thermal <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_thermal,])
    ecoregions_climate_precipitation <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_ecoregions,][wintering_pop_ecoregions_precipitation,])
    colnames(ecoregions_climate_climate) <- colnames(ecoregions_climate_thermal) <- colnames(ecoregions_climate_precipitation) <- c("ID", "temp", "prec")
    
    # Resample climate pixels across the species geographical distribution based on ecoregion weights
    ecoregions_climate_resample_wintering <- ecoregions_climate_resample_wintering_orsim <- ecoregions_climate_resample_breeding <- ecoregions_climate_resample_breeding_orsim <- vector()
    to_sample_wintering <- round(ecoregions_weights_wintering * 10000)
    for(h in 1:length(to_sample_wintering)){
      ecoregions_climate_resample_wintering <- rbind(ecoregions_climate_resample_wintering, ecoregions_climate_wintering[which(ecoregions_climate_wintering$ID == h),][sample(1:length(which(ecoregions_climate_wintering$ID == h)), to_sample_wintering[h], replace=T),])
    }
    ecoregions_climate_resample_wintering <- ecoregions_climate_resample_wintering %>% mutate(population = pops_names[j-3])
    to_sample_breeding <- round(ecoregions_weights_breeding * 10000)
    for(h in 1:length(to_sample_breeding)){
      ecoregions_climate_resample_breeding <- rbind(ecoregions_climate_resample_breeding, ecoregions_climate_breeding[which(ecoregions_climate_breeding$ID == h),][sample(1:length(which(ecoregions_climate_breeding$ID == h)), to_sample_breeding[h], replace=T),])
    }
    ecoregions_climate_resample_breeding <- ecoregions_climate_resample_breeding %>% mutate(population = pops_names[j-3])
    to_sample_wintering_orsim <- round(ecoregions_weights_wintering_orsim * 10000)
    for(h in 1:length(to_sample_wintering_orsim)){
      ecoregions_climate_resample_wintering_orsim <- rbind(ecoregions_climate_resample_wintering_orsim, ecoregions_climate_wintering_orsim[which(ecoregions_climate_wintering_orsim$ID == h),][sample(1:length(which(ecoregions_climate_wintering_orsim$ID == h)), to_sample_wintering_orsim[h], replace=T),])
    }
    ecoregions_climate_resample_wintering_orsim <- ecoregions_climate_resample_wintering_orsim %>% mutate(population = pops_names[j-3])
    to_sample_breeding_orsim <- round(ecoregions_weights_breeding_orsim * 10000)
    for(h in 1:length(to_sample_breeding_orsim)){
      ecoregions_climate_resample_breeding_orsim <- rbind(ecoregions_climate_resample_breeding_orsim, ecoregions_climate_breeding[which(ecoregions_climate_breeding$ID == h),][sample(1:length(which(ecoregions_climate_breeding$ID == h)), to_sample_breeding_orsim[h], replace=T),])
    }
    ecoregions_climate_resample_breeding_orsim <- ecoregions_climate_resample_breeding_orsim %>% mutate(population = pops_names[j-3])
    
    # Calculate seasonal climate niches
    niche_climate_wintering[[j-3]] <- nicheDensityRaster(ecoregions_climate_resample_wintering %>% dplyr::select(temp, prec))
    niche_temp_wintering[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_resample_wintering %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_prec_wintering[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_resample_wintering %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_climate_breeding[[j-3]] <- nicheDensityRaster(ecoregions_climate_resample_breeding %>% dplyr::select(temp, prec))
    niche_temp_breeding[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_resample_breeding %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_prec_breeding[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_resample_breeding %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_climate_wintering_orsim[[j-3]] <- nicheDensityRaster(ecoregions_climate_resample_wintering_orsim %>% dplyr::select(temp, prec))
    niche_temp_wintering_orsim[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_resample_wintering_orsim %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_prec_wintering_orsim[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_resample_wintering_orsim %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_climate_breeding_orsim[[j-3]] <- nicheDensityRaster(ecoregions_climate_resample_breeding_orsim %>% dplyr::select(temp, prec))
    niche_temp_breeding_orsim[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_resample_breeding_orsim %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_prec_breeding_orsim[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_resample_breeding_orsim %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_climate_wintering_climate[[j-3]] <- nicheDensityRaster(ecoregions_climate_climate %>% dplyr::select(temp, prec))
    niche_temp_wintering_climate[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_climate %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_prec_wintering_climate[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_climate %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_climate_wintering_thermal[[j-3]] <- nicheDensityRaster(ecoregions_climate_thermal %>% dplyr::select(temp, prec))
    niche_temp_wintering_thermal[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_thermal %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_prec_wintering_thermal[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_thermal %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_climate_wintering_precipitation[[j-3]] <- nicheDensityRaster(ecoregions_climate_precipitation %>% dplyr::select(temp, prec))
    niche_temp_wintering_precipitation[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_precipitation %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_prec_wintering_precipitation[[j-3]] <- stats::density(unlist(as.vector(ecoregions_climate_precipitation %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_climate_wintering_n <- niche_temp_wintering_n <- niche_prec_wintering_n <- list()
    for(i in 1:1000){
      ecoregions_climate_null <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[pop_resampling[i,],])
      colnames(ecoregions_climate_null) <- c("ID", "temp", "prec")
      ecoregions_weights_null <- ecoregions_abund_presence[[k]]$wintering_abund[pop_resampling[i,]]
      to_sample_null <- round((ecoregions_weights_null/sum(ecoregions_weights_null)) * 10000)
      ecoregions_climate_resample_wintering_null <- vector()
      for(h in 1:length(to_sample_null)){
        ecoregions_climate_resample_wintering_null <- rbind(ecoregions_climate_resample_wintering_null, ecoregions_climate_null[which(ecoregions_climate_null$ID == h),][sample(1:length(which(ecoregions_climate_null$ID == h)), to_sample_null[h], replace=T),])
      }
      ecoregions_climate_resample_wintering_null <- ecoregions_climate_resample_wintering_null %>% mutate(population = pops_names[j-3])
      
      niche_climate_wintering_n[[i]] <- nicheDensityRaster(ecoregions_climate_resample_wintering_null %>% dplyr::select(temp, prec))
      niche_temp_wintering_n[[i]] <- stats::density(unlist(as.vector(ecoregions_climate_resample_wintering_null %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
      niche_prec_wintering_n[[i]] <- stats::density(unlist(as.vector(ecoregions_climate_resample_wintering_null %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    }
    niche_climate_wintering_null[[j-3]] <- niche_climate_wintering_n
    niche_temp_wintering_null[[j-3]] <- niche_temp_wintering_n
    niche_prec_wintering_null[[j-3]] <- niche_prec_wintering_n
  }
  
  # population centroids
  pop_centroid_breeding[[k]] <- centroid_breeding_ecoregions
  pop_centroid_wintering[[k]] <- centroid_wintering_ecoregions
  pop_centroid_breeding_orsim[[k]] <- centroid_breeding_ecoregions_orsim
  pop_centroid_wintering_orsim[[k]] <- centroid_wintering_ecoregions_orsim
  pop_centroid_wintering_climate[[k]] <- centroid_wintering_ecoregions_climate
  pop_centroid_wintering_thermal[[k]] <- centroid_wintering_ecoregions_thermal
  pop_centroid_wintering_precipitation[[k]] <- centroid_wintering_ecoregions_precipitation
  
  # migration distance
  migration_distance[[k]] <- diag(rdist.earth(centroid_wintering_ecoregions, centroid_breeding_ecoregions, miles=F))
  migration_distance_orsim[[k]] <- diag(rdist.earth(centroid_wintering_ecoregions_orsim, centroid_breeding_ecoregions_orsim, miles=F))
  migration_distance_climate[[k]] <- diag(rdist.earth(centroid_wintering_ecoregions_climate, centroid_breeding_ecoregions, miles=F))
  migration_distance_thermal[[k]] <- diag(rdist.earth(centroid_wintering_ecoregions_thermal, centroid_breeding_ecoregions, miles=F))
  migration_distance_precipitation[[k]] <- diag(rdist.earth(centroid_wintering_ecoregions_precipitation, centroid_breeding_ecoregions, miles=F))
  
  # Seeasonal climate overlap
  density_wintering <- rasterToPoints(raster::stack(niche_climate_wintering))
  density_breeding <- rasterToPoints(raster::stack(niche_climate_breeding))
  density_wintering_orsim <- rasterToPoints(raster::stack(niche_climate_wintering_orsim))
  density_breeding_orsim <- rasterToPoints(raster::stack(niche_climate_breeding_orsim))
  density_wintering_climate <- rasterToPoints(raster::stack(niche_climate_wintering_climate))
  density_wintering_thermal <- rasterToPoints(raster::stack(niche_climate_wintering_thermal))
  density_wintering_precipitation <- rasterToPoints(raster::stack(niche_climate_wintering_precipitation))
  niche_overlap[[k]] <- 1 - (0.5 * apply(abs(density_breeding[,-c(1,2)] - density_wintering[,-c(1,2)]), 2, sum))
  niche_overlap_orsim[[k]] <- 1 - (0.5 * apply(abs(density_breeding_orsim[,-c(1,2)] - density_wintering_orsim[,-c(1,2)]), 2, sum))
  niche_overlap_climate[[k]] <- 1 - (0.5 * apply(abs(density_breeding[,-c(1,2)] - density_wintering_climate[,-c(1,2)]), 2, sum))
  niche_overlap_thermal[[k]] <- 1 - (0.5 * apply(abs(density_breeding[,-c(1,2)] - density_wintering_thermal[,-c(1,2)]), 2, sum))
  niche_overlap_precipitation[[k]] <- 1 - (0.5 * apply(abs(density_breeding[,-c(1,2)] - density_wintering_precipitation[,-c(1,2)]), 2, sum))
  
  niche_overlap_null_pops <- list()
  for(i in 1:length(unique(pops_names))){
    density_wintering_null <- rasterToPoints(raster::stack(niche_climate_wintering_null[[i]]))
    niche_overlap_null_pops[[i]] <- 1 - (0.5 * apply(abs(density_breeding_orsim[,-c(1,2)][,i] - density_wintering_null[,-c(1,2)]), 2, sum))
  }
  niche_overlap_null[[k]] <- do.call(cbind, niche_overlap_null_pops)
  
  # Seasonal 1-D niche overlap
  niche_overlap_temp <- niche_overlap_prec <- vector()
  for(i in 1:length(unique(pops_names))){
    # temperature
    X <- niche_temp_breeding[[i]]$x
    Y1 <- niche_temp_breeding[[i]]$y
    Y2 <- niche_temp_wintering[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche_overlap_temp[i]  <- trapz ( X, Overlap ) / Total
    # precipitation
    X <- niche_prec_breeding[[i]]$x
    Y1 <- niche_prec_breeding[[i]]$y
    Y2 <- niche_prec_wintering[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche_overlap_prec[i] <- trapz ( X, Overlap ) / Total
  }
  niche_overlap_temperature[[k]] <- niche_overlap_temp
  niche_overlap_precipitation[[k]] <- niche_overlap_prec
  
  # Seasonal 1-D niche overlap for Orsim
  niche_overlap_temp_orsim <- niche_overlap_prec_orsim <- vector()
  for(i in 1:length(unique(pops_names))){
    # temperature
    X <- niche_temp_breeding_orsim[[i]]$x
    Y1 <- niche_temp_breeding_orsim[[i]]$y
    Y2 <- niche_temp_wintering_orsim[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche_overlap_temp_orsim[i]  <- trapz ( X, Overlap ) / Total
    # precipitation
    X <- niche_prec_breeding_orsim[[i]]$x
    Y1 <- niche_prec_breeding_orsim[[i]]$y
    Y2 <- niche_prec_wintering_orsim[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche_overlap_prec_orsim[i] <- trapz ( X, Overlap ) / Total
  }
  niche_overlap_temperature_orsim[[k]] <- niche_overlap_temp_orsim
  niche_overlap_precipitation_orsim[[k]] <- niche_overlap_prec_orsim
  
  # Seasonal 1-D niche overlap for climate tracking model
  niche_overlap_temp_climate <- niche_overlap_prec_climate <- vector()
  for(i in 1:length(unique(pops_names))){
    # temperature
    X <- niche_temp_breeding[[i]]$x
    Y1 <- niche_temp_breeding[[i]]$y
    Y2 <- niche_temp_wintering_climate[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche_overlap_temp_climate[i]  <- trapz ( X, Overlap ) / Total
    # precipitation
    X <- niche_prec_breeding[[i]]$x
    Y1 <- niche_prec_breeding[[i]]$y
    Y2 <- niche_prec_wintering_climate[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche_overlap_prec_climate[i] <- trapz ( X, Overlap ) / Total
  }
  niche_overlap_temperature_climate[[k]] <- niche_overlap_temp_climate
  niche_overlap_precipitation_climate[[k]] <- niche_overlap_prec_climate
  
  # Seasonal 1-D niche overlap for thermal tracking model
  niche_overlap_temp_thermal <- niche_overlap_prec_thermal <- vector()
  for(i in 1:length(unique(pops_names))){
    # temperature
    X <- niche_temp_breeding[[i]]$x
    Y1 <- niche_temp_breeding[[i]]$y
    Y2 <- niche_temp_wintering_thermal[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche_overlap_temp_thermal[i]  <- trapz ( X, Overlap ) / Total
    # precipitation
    X <- niche_prec_breeding[[i]]$x
    Y1 <- niche_prec_breeding[[i]]$y
    Y2 <- niche_prec_wintering_thermal[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche_overlap_prec_thermal[i] <- trapz ( X, Overlap ) / Total
  }
  niche_overlap_temperature_thermal[[k]] <- niche_overlap_temp_thermal
  niche_overlap_precipitation_thermal[[k]] <- niche_overlap_prec_thermal
  
  # Seasonal 1-D niche overlap for precipitation tracking model
  niche_overlap_temp_precipitation <- niche_overlap_prec_precipitation <- vector()
  for(i in 1:length(unique(pops_names))){
    # temperature
    X <- niche_temp_breeding[[i]]$x
    Y1 <- niche_temp_breeding[[i]]$y
    Y2 <- niche_temp_wintering_precipitation[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche_overlap_temp_precipitation[i]  <- trapz ( X, Overlap ) / Total
    # precipitation
    X <- niche_prec_breeding[[i]]$x
    Y1 <- niche_prec_breeding[[i]]$y
    Y2 <- niche_prec_wintering_precipitation[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche_overlap_prec_precipitation[i] <- trapz ( X, Overlap ) / Total
  }
  niche_overlap_temperature_precipitation[[k]] <- niche_overlap_temp_precipitation
  niche_overlap_precipitation_precipitation[[k]] <- niche_overlap_prec_precipitation
  
  # Seasonal 1-D niche overlap for null model
  niche_overlap_temp_null <- niche_overlap_prec_null <- list()
  for(i in 1:length(unique(pops_names))){
    niche_overlap_temp_null_pops <- niche_overlap_prec_null_pops <- vector()
    for(repet in 1:length(niche_temp_wintering_null[[i]])){
      # temperature
      X <- niche_temp_breeding_orsim[[i]]$x
      Y1 <- niche_temp_breeding_orsim[[i]]$y
      Y2 <- niche_temp_wintering_null[[i]][[repet]]$y
      Overlap <- pmin ( Y1, Y2 )
      Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
      niche_overlap_temp_null_pops[repet]  <- trapz ( X, Overlap ) / Total
      # precipitation
      X <- niche_prec_breeding_orsim[[i]]$x
      Y1 <- niche_prec_breeding_orsim[[i]]$y
      Y2 <- niche_prec_wintering_null[[i]][[repet]]$y
      Overlap <- pmin ( Y1, Y2 )
      Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
      niche_overlap_prec_null_pops[repet] <- trapz ( X, Overlap ) / Total
    }
    niche_overlap_temp_null[[i]] <- niche_overlap_temp_null_pops
    niche_overlap_prec_null[[i]] <- niche_overlap_prec_null_pops
  }
  niche_overlap_temperature_null[[k]] <- do.call(cbind, niche_overlap_temp_null)
  niche_overlap_precipitation_null[[k]] <- do.call(cbind, niche_overlap_prec_null)
  print(k)
}

pop_centroid_breeding <- do.call(rbind, pop_centroid_breeding)
pop_centroid_wintering <- do.call(rbind, pop_centroid_wintering)
pop_centroid_breeding_orsim <- do.call(rbind, pop_centroid_breeding_orsim)
pop_centroid_wintering_orsim <- do.call(rbind, pop_centroid_wintering_orsim)
pop_centroid_wintering_climate <- do.call(rbind, pop_centroid_wintering_climate)
pop_centroid_wintering_thermal <- do.call(rbind, pop_centroid_wintering_thermal)
pop_centroid_wintering_precipitation <- do.call(rbind, pop_centroid_wintering_precipitation)

# Put all the explanatory variables in the same data frame
niche_data <- data.frame(
  species = rep(species.names, unlist(lapply(population_names, length))),
  population = unlist(population_names),
  niche_overlap = unlist(niche_overlap),
  niche_overlap_temp = unlist(niche_overlap_temperature),
  niche_overlap_prec = unlist(niche_overlap_precipitation),
  niche_overlap_orsim = unlist(niche_overlap_orsim),
  niche_overlap_temp_orsim = unlist(niche_overlap_temperature_orsim),
  niche_overlap_prec_orsim = unlist(niche_overlap_precipitation_orsim),
  niche_overlap_climate = unlist(niche_overlap_climate),
  niche_overlap_temp_climate = unlist(niche_overlap_temperature_climate),
  niche_overlap_prec_climate = unlist(niche_overlap_precipitation_climate),
  niche_overlap_thermal = unlist(niche_overlap_thermal),
  niche_overlap_temp_thermal = unlist(niche_overlap_temperature_thermal),
  niche_overlap_prec_thermal = unlist(niche_overlap_precipitation_thermal),
  niche_overlap_precipitation = unlist(niche_overlap_precipitation),
  niche_overlap_temp_precipitation = unlist(niche_overlap_temperature_precipitation),
  niche_overlap_prec_precipitation = unlist(niche_overlap_precipitation_precipitation),
  migration_distance = unlist(migration_distance),
  migration_distance_orsim = unlist(migration_distance_orsim),
  migration_distance_climate = unlist(migration_distance_climate),
  migration_distance_thermal = unlist(migration_distance_thermal),
  migration_distance_precipitation = unlist(migration_distance_precipitation),
  long_breeding = pop_centroid_breeding[,1],
  lat_breeding = pop_centroid_breeding[,2],
  long_wintering = pop_centroid_wintering[,1],
  lat_wintering = pop_centroid_wintering[,2],
  long_breeding_orsim = pop_centroid_breeding_orsim[,1],
  lat_breeding_orsim = pop_centroid_breeding_orsim[,2],
  long_wintering_orsim = pop_centroid_wintering_orsim[,1],
  lat_wintering_orsim = pop_centroid_wintering_orsim[,2],
  long_wintering_climate = pop_centroid_wintering_climate[,1],
  lat_wintering_climate = pop_centroid_wintering_climate[,2],
  long_wintering_thermal = pop_centroid_wintering_thermal[,1],
  lat_wintering_thermal = pop_centroid_wintering_thermal[,2],
  long_wintering_precipitation = pop_centroid_wintering_precipitation[,1],
  lat_wintering_precipitation = pop_centroid_wintering_precipitation[,2]
)
niche_data <- niche_data %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "PABU_Louisiana", "SWTH_PNW"))


# rename populations
niche_data$population <- c("W", "SR", "CE", "SW",
                           "NW", "NWC", "SR", "C",
                           "NW", "R", "W", "SW1", "SW2", "NE",
                           "RC", "SW", "CE", "E",
                           "NWR", "SE", "NCE",
                           "SC",
                           "NWCE", "SR", "R", "NW", "SW")

# population centroids
pop_centroids_df <- niche_data[,c("species", "long_breeding", "lat_breeding", "long_wintering", "lat_wintering")]
pop_centroids_sf <- rbind(pop_centroids_df %>% dplyr::select(c("species", "long_breeding", "lat_breeding")) %>% mutate(season = "breeding") %>% rename(long = long_breeding, lat = lat_breeding),
                          pop_centroids_df %>% dplyr::select(c("species", "long_wintering", "lat_wintering")) %>% mutate(season = "wintering") %>% rename(long = long_wintering, lat = lat_wintering)) %>% 
  st_as_sf(coords = c("long", "lat"), crs='+proj=longlat +datum=WGS84')
pop_centroids_orsim_df <- niche_data[,c("species", "long_breeding_orsim", "lat_breeding_orsim", "long_wintering_orsim", "lat_wintering_orsim")]
pop_centroids_orsim_sf <- rbind(pop_centroids_orsim_df %>% dplyr::select(c("species", "long_breeding_orsim", "lat_breeding_orsim")) %>% mutate(season = "breeding") %>% rename(long = long_breeding_orsim, lat = lat_breeding_orsim),
                                pop_centroids_orsim_df %>% dplyr::select(c("species", "long_wintering_orsim", "lat_wintering_orsim")) %>% mutate(season = "wintering") %>% rename(long = long_wintering_orsim, lat = lat_wintering_orsim)) %>% 
  st_as_sf(coords = c("long", "lat"), crs='+proj=longlat +datum=WGS84')
pop_centroids_climate_df <- niche_data[,c("species", "long_breeding", "lat_breeding", "long_wintering_climate", "lat_wintering_climate")]
pop_centroids_climate_sf <- rbind(pop_centroids_climate_df %>% dplyr::select(c("species", "long_breeding", "lat_breeding")) %>% mutate(season = "breeding") %>% rename(long = long_breeding, lat = lat_breeding),
                                  pop_centroids_climate_df %>% dplyr::select(c("species", "long_wintering_climate", "lat_wintering_climate")) %>% mutate(season = "wintering") %>% rename(long = long_wintering_climate, lat = lat_wintering_climate)) %>% 
  st_as_sf(coords = c("long", "lat"), crs='+proj=longlat +datum=WGS84')
pop_centroids_thermal_df <- niche_data[,c("species", "long_breeding", "lat_breeding", "long_wintering_thermal", "lat_wintering_thermal")]
pop_centroids_thermal_sf <- rbind(pop_centroids_thermal_df %>% dplyr::select(c("species", "long_breeding", "lat_breeding")) %>% mutate(season = "breeding") %>% rename(long = long_breeding, lat = lat_breeding),
                                  pop_centroids_thermal_df %>% dplyr::select(c("species", "long_wintering_thermal", "lat_wintering_thermal")) %>% mutate(season = "wintering") %>% rename(long = long_wintering_thermal, lat = lat_wintering_thermal)) %>% 
  st_as_sf(coords = c("long", "lat"), crs='+proj=longlat +datum=WGS84')
pop_centroids_precipitation_df <- niche_data[,c("species", "long_breeding", "lat_breeding", "long_wintering_precipitation", "lat_wintering_precipitation")]
pop_centroids_precipitation_sf <- rbind(pop_centroids_precipitation_df %>% dplyr::select(c("species", "long_breeding", "lat_breeding")) %>% mutate(season = "breeding") %>% rename(long = long_breeding, lat = lat_breeding),
                                        pop_centroids_precipitation_df %>% dplyr::select(c("species", "long_wintering_precipitation", "lat_wintering_precipitation")) %>% mutate(season = "wintering") %>% rename(long = long_wintering_precipitation, lat = lat_wintering_precipitation)) %>% 
  st_as_sf(coords = c("long", "lat"), crs='+proj=longlat +datum=WGS84')


## FIGURES ##

## Figure 1: plot patterns of winter and summer climate

g_temp_BR <- ggplot() +
  geom_spatraster(data=temp_zscore_BR) + ylim(c(35, 72)) + xlim(c(-150, -50)) +
  theme_void() + scale_fill_hypso_tint_c(palette = "colombia_hypso", limits = c(-0.05, 1.35), direction=-1) +
  theme(legend.position="none")

g_temp_NB <-ggplot() +
  geom_spatraster(data=temp_zscore_NB) + ylim(c(5, 33.5)) + xlim(c(-135, -60)) + 
  theme_void() + scale_fill_hypso_tint_c(palette = "colombia_hypso", limits = c(-0.05, 1.3), direction=-1) +
  theme(legend.position="none")

pdf("results/Fig_1_temp_BR.pdf", width = 4, height = 4)
g_temp_BR
dev.off()
pdf("results/Fig_1_temp_NB.pdf", width = 4, height = 4)
g_temp_NB
dev.off()

                               
# Figure S10: examine relationships between climate and species relative abundance

ecoregions_temp_winter_mean <- ecoregions_prec_winter_mean <- ecoregions_temp_summer_mean <- ecoregions_prec_summer_mean <- vector()
for(i in 1:length(ecoregions)){
  ecoregions_temp_winter_mean[i] <- mean(ecoregions_climate_NB$temperature[ecoregions_climate_NB$ID == i])
  ecoregions_prec_winter_mean[i] <- mean(ecoregions_climate_NB$precipitation[ecoregions_climate_NB$ID == i])
  ecoregions_temp_summer_mean[i] <- mean(ecoregions_climate_BR$temperature[ecoregions_climate_BR$ID == i])
  ecoregions_prec_summer_mean[i] <- mean(ecoregions_climate_BR$precipitation[ecoregions_climate_BR$ID == i])
  print(i)
}

# Plot for all 7 species (k = 1,...,7) 
k=7
ecoregion_occupied_W <- which(ecoregions_abund_presence[[k]]$wintering_abund > 0)
ecoregion_occupied_B <- which(ecoregions_abund_presence[[k]]$breeding_abund > 0)
data_for_plot_W <- data.frame(wintering_abund = ecoregions_abund_presence[[k]]$wintering_abund[ecoregion_occupied_W],
                              temp_winter = ecoregions_temp_winter_mean[ecoregion_occupied_W],
                              prec_winter = ecoregions_prec_winter_mean[ecoregion_occupied_W])
data_for_plot_B <- data.frame(breeding_abund = ecoregions_abund_presence[[k]]$breeding_abund[ecoregion_occupied_B],
                              temp_summer = ecoregions_temp_summer_mean[ecoregion_occupied_B],
                              prec_summer = ecoregions_prec_summer_mean[ecoregion_occupied_B])
g_T_W <- ggplot() +
  geom_point(data = data_for_plot_W, aes(x=temp_winter, y=wintering_abund)) +
  geom_smooth(data = data_for_plot_W, aes(x=temp_winter, y=wintering_abund), se=T, span=1) +
  xlab("") + ylab("")
g_P_W <- ggplot() +
  geom_point(data = data_for_plot_W, aes(x=prec_winter, y=wintering_abund)) +
  geom_smooth(data = data_for_plot_W, aes(x=prec_winter, y=wintering_abund), se=T, span=1) +
  xlab("") + ylab("")
g_T_B <- ggplot() +
  geom_point(data = data_for_plot_B, aes(x=temp_summer, y=breeding_abund)) +
  geom_smooth(data = data_for_plot_B, aes(x=temp_summer, y=breeding_abund), se=T, span=1) +
  xlab("") + ylab("")
g_P_B <- ggplot() +
  geom_point(data = data_for_plot_B, aes(x=prec_summer, y=breeding_abund)) +
  geom_smooth(data = data_for_plot_B, aes(x=prec_summer, y=breeding_abund), se=T, span=1) +
  xlab("") + ylab("")
pdf(paste0("results/Fig_SM_climate_abund_", k, ".pdf"), width = 12, height = 3)
ggarrange(g_T_B, g_P_B, g_T_W, g_P_W, nrow=1)
dev.off()


## Figure S11 and S12: Check for climate bias in sampling

ecoregions_W_all_climate <- ecoregions_B_all_climate <- ecoregions_W_samples_climate <- ecoregions_B_samples_climate <- list()
for(k in 1:length(species.names)){
  breeding_ecoregions <- which(ecoregions_abund_presence[[k]]$breeding_abund > 0)
  wintering_ecoregions <- which(ecoregions_abund_presence[[k]]$wintering_abund > 0)
  pops_seas <- do.call(rbind, strsplit(colnames(ecoregions_abund_presence[[k]])[4:ncol(ecoregions_abund_presence[[k]])], "_"))
  pops_names <- unique(pops_seas[,1])
  wintering_ecoregions_samples <- which(apply(ecoregions_abund_presence[[k]][,4:(3+length(pops_names))], 1, sum) > 0)
  breeding_ecoregions_samples <- which(apply(ecoregions_abund_presence[[k]][,(4+length(pops_names)):ncol(ecoregions_abund_presence[[k]])], 1, sum) > 0)
  ecoregions_W_all_climate[[k]] <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_ecoregions,], fun=mean)
  ecoregions_B_all_climate[[k]] <- terra::extract(c(temp_zscore_BR, prec_zscore_BR), ecoregions[breeding_ecoregions,], fun=mean)
  ecoregions_W_samples_climate[[k]] <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_ecoregions_samples,], fun=mean)
  ecoregions_B_samples_climate[[k]] <- terra::extract(c(temp_zscore_BR, prec_zscore_BR), ecoregions[breeding_ecoregions_samples,], fun=mean)
  colnames(ecoregions_W_all_climate[[k]]) <- colnames(ecoregions_B_all_climate[[k]]) <- colnames(ecoregions_W_samples_climate[[k]]) <- colnames(ecoregions_B_samples_climate[[k]]) <- c("ID", "temp", "prec")
}

# Fig S11: Plot breeding temperature and precipitation density of all ecoregions occupied versus sampled ecoregions
k=1
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]])))) 
g_temp_wifl <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(a) WIFL temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_wifl <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(b) WIFL precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=2
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_yewa <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(c) YEWA temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_yewa <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(d) YEWA precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=3
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_wiwa <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(e) WIWA temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_wiwa <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(f) WIWA precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=4
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_coye <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(g) COYE temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_coye <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(h) COYE precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=5
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_amre <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(i) AMRE temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_amre <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(j) AMRE precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=6
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_pabu <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(k) PABU temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_pabu <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(l) PABU precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=7
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_heth <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(m) HETH temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_heth <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(n) HETH precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())

pdf("results/Fig_S11_new2.pdf", width = 8, height = 15)
ggarrange(g_temp_wifl, g_prec_wifl,
          g_temp_yewa, g_prec_yewa,
          g_temp_wiwa, g_prec_wiwa,
          g_temp_coye, g_prec_coye,
          g_temp_amre, g_prec_amre, 
          g_temp_pabu, g_prec_pabu, 
          g_temp_heth, g_prec_heth, nrow=7, ncol=2, common.legend = TRUE, legend="bottom")
dev.off()

# Fig S12: Plot wintering temperature and precipitation density of all ecoregions occupied versus sampled ecoregions
k=1
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]])))) 
g_temp_wifl <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(a) WIFL temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_wifl <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(b) WIFL precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=2
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_yewa <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(c) YEWA temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_yewa <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(d) YEWA precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=3
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_wiwa <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(e) WIWA temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_wiwa <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(f) WIWA precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=4
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_coye <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(g) COYE temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_coye <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(h) COYE precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=5
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_amre <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(i) AMRE temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_amre <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(j) AMRE precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=6
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_pabu <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(k) PABU temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_pabu <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(l) PABU precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=7
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_heth <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-0.8, 1.5)) + ggtitle("(m) HETH temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_heth <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) + theme_classic() +
  xlim(c(-2, 4)) + ggtitle("(n) HETH precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())

pdf("results/Fig_S12_new2.pdf", width = 8, height = 15)
ggarrange(g_temp_wifl, g_prec_wifl,
          g_temp_yewa, g_prec_yewa,
          g_temp_wiwa, g_prec_wiwa,
          g_temp_coye, g_prec_coye,
          g_temp_amre, g_prec_amre,
          g_temp_pabu, g_prec_pabu, 
          g_temp_heth, g_prec_heth, nrow=7, ncol=2, common.legend = TRUE, legend="bottom")
dev.off()


## Figures 2 and S9: patterns of simulated and empirical migratory connectivity and seasonal climate tracking
# Migratory connectivity
g_conn <- ggplot() + geom_sf(data = newmap) + ylim(c(0,68)) + xlim(c(-158,-51)) +
  geom_sf(data = pop_centroids_sf, aes(col=species, shape=season)) + 
  geom_segment(data = pop_centroids_df, aes(x=long_breeding, y=lat_breeding, xend=long_wintering, yend=lat_wintering, col=species)) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  ggtitle("(a) Empirical migratory connectivity") + theme_void()
g_conn_orsim <- ggplot() + geom_sf(data = newmap) + ylim(c(0,68)) + xlim(c(-158,-51)) +
  geom_sf(data = pop_centroids_orsim_sf, aes(col=species, shape=season)) + 
  geom_segment(data = pop_centroids_orsim_df, aes(x=as.numeric(long_breeding_orsim), y=as.numeric(lat_breeding_orsim), xend=as.numeric(long_wintering_orsim), yend=as.numeric(lat_wintering_orsim), col=species)) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  ggtitle("(f) Simulated connectivity (ORSIM)") + theme_void()
g_conn_climate <- ggplot() + geom_sf(data = newmap) + ylim(c(0,68)) + xlim(c(-158,-51)) +
  geom_sf(data = pop_centroids_climate_sf, aes(col=species, shape=season)) + 
  geom_segment(data = pop_centroids_climate_df, aes(x=as.numeric(long_breeding), y=as.numeric(lat_breeding), xend=as.numeric(long_wintering_climate), yend=as.numeric(lat_wintering_climate), col=species)) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  ggtitle("(k) Simulated connectivity (climate)") + theme_void()
g_conn_thermal <- ggplot() + geom_sf(data = newmap) + ylim(c(0,68)) + xlim(c(-158,-51)) +
  geom_sf(data = pop_centroids_thermal_sf, aes(col=species, shape=season)) + 
  geom_segment(data = pop_centroids_thermal_df, aes(x=as.numeric(long_breeding), y=as.numeric(lat_breeding), xend=as.numeric(long_wintering_thermal), yend=as.numeric(lat_wintering_thermal), col=species)) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  ggtitle("(k) Simulated connectivity (thermal)") + theme_void()
g_conn_precipitation <- ggplot() + geom_sf(data = newmap) + ylim(c(0,68)) + xlim(c(-158,-51)) +
  geom_sf(data = pop_centroids_precipitation_sf, aes(col=species, shape=season)) + 
  geom_segment(data = pop_centroids_precipitation_df, aes(x=as.numeric(long_breeding), y=as.numeric(lat_breeding), xend=as.numeric(long_wintering_precipitation), yend=as.numeric(lat_wintering_precipitation), col=species)) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  ggtitle("(p) Simulated connectivity (precipitation)") + theme_void()
# Overlap versus distance plots
g_overlap_dist <- ggplot(niche_data, aes(x = migration_distance, y = niche_overlap, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.95)) + xlim(c(0,7000)) + ggtitle("(c)") + theme_classic() +
  xlab("Empirical migration distance (km)") + ylab("Empirical 2D climate overlap")
g_overlap_temp_dist <- ggplot(niche_data, aes(x = migration_distance, y = niche_overlap_temp, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.6)) + xlim(c(0,7000)) + ggtitle("(d)") + theme_classic() +
  xlab("Empirical migration distance (km)") + ylab("Empirical thermal overlap")
g_overlap_prec_dist <- ggplot(niche_data, aes(x = migration_distance, y = niche_overlap_prec, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.6)) + xlim(c(0,7000)) + ggtitle("(e)") + theme_classic() +
  xlab("Empirical migration distance (km)") + ylab("Empirical precipitation overlap")
g_overlap_dist_orsim <- ggplot(niche_data, aes(x = migration_distance_orsim, y = niche_overlap_orsim, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.95)) + xlim(c(0,7000)) + ggtitle("(h)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated 2D climate overlap (ORSIM)")
g_overlap_temp_dist_orsim <- ggplot(niche_data, aes(x = migration_distance_orsim, y = niche_overlap_temp_orsim, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.6)) + xlim(c(0,7000)) + ggtitle("(i)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated thermal overlap (ORSIM)")
g_overlap_prec_dist_orsim <- ggplot(niche_data, aes(x = migration_distance_orsim, y = niche_overlap_prec_orsim, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.6)) + xlim(c(0,7000)) + ggtitle("(j)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated precipitation overlap (ORSIM)")
g_overlap_dist_climate <- ggplot(niche_data, aes(x = migration_distance_climate, y = niche_overlap_climate, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.95)) + xlim(c(0,7000)) + ggtitle("(m)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated 2D climate overlap (climate)")
g_overlap_temp_dist_climate <- ggplot(niche_data, aes(x = migration_distance_climate, y = niche_overlap_temp_climate, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.6)) + xlim(c(0,7000)) + ggtitle("(n)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated thermal overlap (climate)")
g_overlap_prec_dist_climate <- ggplot(niche_data, aes(x = migration_distance_climate, y = niche_overlap_prec_climate, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.6)) + xlim(c(0,7000)) + ggtitle("(o)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated precipitation overlap (climate)")
g_overlap_dist_thermal <- ggplot(niche_data, aes(x = migration_distance_thermal, y = niche_overlap_thermal, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.95)) + xlim(c(0,7000)) + ggtitle("(m)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated 2D climate overlap (thermal)")
g_overlap_temp_dist_thermal <- ggplot(niche_data, aes(x = migration_distance_thermal, y = niche_overlap_temp_thermal, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.6)) + xlim(c(0,7000)) + ggtitle("(n)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated thermal overlap (thermal)")
g_overlap_prec_dist_thermal <- ggplot(niche_data, aes(x = migration_distance_thermal, y = niche_overlap_prec_thermal, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.6)) + xlim(c(0,7000)) + ggtitle("(o)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated precipitation overlap (thermal)")
g_overlap_dist_precipitation <- ggplot(niche_data, aes(x = migration_distance_precipitation, y = niche_overlap_precipitation, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.95)) + xlim(c(0,7000)) + ggtitle("(r)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated 2D climate overlap (precipitation)")
g_overlap_temp_dist_precipitation <- ggplot(niche_data, aes(x = migration_distance_precipitation, y = niche_overlap_temp_precipitation, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.6)) + xlim(c(0,7000)) + ggtitle("(s)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated thermal overlap (precipitation)")
g_overlap_prec_dist_precipitation <- ggplot(niche_data, aes(x = migration_distance_precipitation, y = niche_overlap_prec_precipitation, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1.5) +
  ylim(c(0,0.6)) + xlim(c(0,7000)) + ggtitle("(t)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated precipitation overlap (precipitation)")
g_overlap_density <- ggplot(niche_data, aes(x = niche_overlap)) +
  geom_density(adjust = 2, fill="grey") + theme_classic() +
  xlim(c(0, 1)) + ggtitle("(b)") + theme(legend.title=element_blank()) +
  xlab("Empirical 2D climate overlap")
g_overlap_density_orsim <- ggplot(niche_data, aes(x = niche_overlap_orsim)) +
  geom_density(adjust = 2, fill="grey") + theme_classic() +
  xlim(c(0, 1)) + ggtitle("(g)") + theme(legend.title=element_blank()) +
  xlab("Simulated 2D climate overlap (ORSIM)")
g_overlap_density_climate <- ggplot(niche_data, aes(x = niche_overlap_climate)) +
  geom_density(adjust = 2, fill="grey") + theme_classic() +
  xlim(c(0, 1)) + ggtitle("(l)") + theme(legend.title=element_blank()) +
  xlab("Simulated 2D climate overlap (climate)")

# Figure 2
pdf("results/Fig_2_new3.pdf", width = 16, height = 10)
ggarrange(g_conn, g_overlap_density, g_overlap_dist, g_overlap_temp_dist, g_overlap_prec_dist, g_conn_orsim, g_overlap_density_orsim, g_overlap_dist_orsim, g_overlap_temp_dist_orsim, g_overlap_prec_dist_orsim, g_conn_climate, g_overlap_density_climate, g_overlap_dist_climate, g_overlap_temp_dist_climate, g_overlap_prec_dist_climate, ncol=5, nrow=3, common.legend = TRUE, legend="bottom")
dev.off()

# Figure S9
pdf("results/Fig_S9_new3.pdf", width = 13, height = 13)
ggarrange(g_conn, g_overlap_dist, g_overlap_temp_dist, g_overlap_prec_dist, g_conn_orsim, g_overlap_dist_orsim, g_overlap_temp_dist_orsim, g_overlap_prec_dist_orsim, g_conn_thermal, g_overlap_dist_thermal, g_overlap_temp_dist_thermal, g_overlap_prec_dist_thermal, g_conn_precipitation, g_overlap_dist_precipitation, g_overlap_temp_dist_precipitation, g_overlap_prec_dist_precipitation, ncol=4, nrow=4, common.legend = TRUE, legend="bottom")
dev.off()


# Correlation between precipitation overal and migration distance
cor_prec_migr <- cor.test(niche_data$migration_distance, niche_data$niche_overlap_prec) 

# Statistical comparison between empirical seasonal climate overlap and the values simulated by ORSIM and the climate tracking model 
ks.test(niche_data$niche_overlap, niche_data$niche_overlap_orsim)
ks.test(niche_data$niche_overlap, niche_data$niche_overlap_climate)

## Testing correlations between simulated and empirical variables
# Correlations with ORSIM
cor_niche <- cor.test(niche_data$niche_overlap, niche_data$niche_overlap_orsim) # seasonal climate overlap
cor_temp <- cor.test(niche_data$niche_overlap_temp, niche_data$niche_overlap_temp_orsim) # seasonal thermal overlap
cor_prec <- cor.test(niche_data$niche_overlap_prec, niche_data$niche_overlap_prec_orsim) # seasonal precipitation overlap
cor_migr <- cor.test(niche_data$migration_distance, niche_data$migration_distance_orsim) # migration distance
# Correlations with climate tracking model
cor_niche_climate <- cor.test(niche_data$niche_overlap, niche_data$niche_overlap_climate) # seasonal climate overlap
cor_temp_climate <- cor.test(niche_data$niche_overlap_temp, niche_data$niche_overlap_temp_climate) # seasonal thermal overlap
cor_prec_climate <- cor.test(niche_data$niche_overlap_prec, niche_data$niche_overlap_prec_climate) # seasonal precipitation overlap
cor_migr_climate <- cor.test(niche_data$migration_distance, niche_data$migration_distance_climate) # migration distance
# Correlations with thermal tracking model
cor_niche_thermal <- cor.test(niche_data$niche_overlap, niche_data$niche_overlap_thermal) # seasonal climate overlap
cor_temp_thermal <- cor.test(niche_data$niche_overlap_temp, niche_data$niche_overlap_temp_thermal) # seasonal thermal overlap
cor_prec_thermal <- cor.test(niche_data$niche_overlap_prec, niche_data$niche_overlap_prec_thermal) # seasonal precipitation overlap
cor_migr_thermal <- cor.test(niche_data$migration_distance, niche_data$migration_distance_thermal) # migration distance
# Correlations with precipitation tracking model
cor_niche_precipitation <- cor.test(niche_data$niche_overlap, niche_data$niche_overlap_precipitation) # seasonal climate overlap
cor_temp_precipitation <- cor.test(niche_data$niche_overlap_temp, niche_data$niche_overlap_temp_precipitation) # seasonal thermal overlap
cor_prec_precipitation <- cor.test(niche_data$niche_overlap_prec, niche_data$niche_overlap_prec_precipitation) # seasonal precipitation overlap
cor_migr_precipitation <- cor.test(niche_data$migration_distance, niche_data$migration_distance_precipitation) # migration distance

##  Are residuals of the relationship with ORSIM skewed towards observed > simulated? ##
orsim_niche_tracking <- niche_data$niche_overlap - niche_data$niche_overlap_orsim
orsim_temp_tracking <- niche_data$niche_overlap_temp - niche_data$niche_overlap_temp_orsim
orsim_prec_tracking <- niche_data$niche_overlap_prec - niche_data$niche_overlap_prec_orsim
hist(orsim_niche_tracking, breaks=10)
hist(orsim_temp_tracking, breaks=10)
hist(orsim_prec_tracking, breaks=10)
ks.test(orsim_niche_tracking,"pnorm", mean=0, sd=sd(orsim_niche_tracking), alternative="less")
ks.test(orsim_temp_tracking,"pnorm", mean=0, sd=sd(orsim_temp_tracking), alternative="less")
ks.test(orsim_prec_tracking,"pnorm", mean=0, sd=sd(orsim_prec_tracking), alternative="less")


# Correlations with ORSIM without HETH SW
niche_data_woHETHSW <- niche_data[-which(niche_data$species == "HETH" & niche_data$population == "SW"),]
cor_niche <- cor.test(niche_data_woHETHSW$niche_overlap, niche_data_woHETHSW$niche_overlap_orsim) 
cor_temp <- cor.test(niche_data_woHETHSW$niche_overlap_temp, niche_data_woHETHSW$niche_overlap_temp_orsim) 
cor_prec <- cor.test(niche_data_woHETHSW$niche_overlap_prec, niche_data_woHETHSW$niche_overlap_prec_orsim)


### Processing the results of the null model randomizing the seasonal grounds after ORSIM

# Comparing null and observed, treating each population as independent
pvals <- pvals_temp <- pvals_prec <- vector()
for(k in 1:length(species.names)){
  niche_data_spp <- niche_data %>% filter(species == species.names[k])
  for(i in 1:nrow(niche_data_spp)){
    pvals <- c(pvals, length(which(niche_overlap_null[[k]][,i] > niche_data_spp$niche_overlap[i])) / nrow(niche_overlap_null[[k]]))
    pvals_temp <- c(pvals_temp, length(which(niche_overlap_temperature_null[[k]][,i] > niche_data_spp$niche_overlap_temp[i])) / nrow(niche_overlap_temperature_null[[k]]))
    pvals_prec <- c(pvals_prec, length(which(niche_overlap_precipitation_null[[k]][,i] > niche_data_spp$niche_overlap_prec[i])) / nrow(niche_overlap_precipitation_null[[k]]))
  }
}
niche_data2 <- cbind(niche_data, pvals, pvals_temp, pvals_prec)

# Are distributions of scaled ranks skewed towards low values? ##
ks.test(niche_data2$pvals, "punif",0,1, alternative="greater") # no
ks.test(niche_data2$pvals_temp, "punif",0,1, alternative="greater") # no
ks.test(niche_data2$pvals_prec, "punif",0,1, alternative="greater") # yes

# Figure 3: Relationship between empirical climate niche overlap and null expectations
g_orsim_1 <- ggplot(niche_data2, aes(x = migration_distance, y = migration_distance_orsim)) +
  geom_point(col = "black", size=2.5) + theme_classic() +
  geom_abline(aes(slope=1, intercept=0)) + ggtitle("(a)") + xlim(c(0,6500)) + ylim(c(0,6500)) +
  xlab("Empirical migration distance (km)") + ylab("Simulated migration distance (km)")+
  scale_size(name="1-rank") + scale_colour_discrete(name="")
g_orsim_2 <- ggplot(niche_data2, aes(x = niche_overlap, y = niche_overlap_orsim)) +
  geom_abline(aes(slope=1, intercept=0)) + ggtitle("(b)") + xlim(c(0,0.9)) + ylim(c(0,0.9)) +
  geom_point(aes(col = species, size = 1-pvals), alpha=0.75) + theme_classic() +
  xlab("Empirical climate overlap") + ylab("Simulated climate overlap")+
  scale_size(name="1-rank") + scale_colour_discrete(name="")
g_orsim_3 <- ggplot(niche_data2, aes(x = niche_overlap_temp, y = niche_overlap_temp_orsim)) +
  geom_abline(aes(slope=1, intercept=0)) + ggtitle("(c)") + xlim(c(0.05,0.5)) + ylim(c(0.05,0.5)) +
  geom_point(aes(col = species, size = 1-pvals_temp), alpha=0.75) + theme_classic() +
  xlab("Empirical thermal overlap") + ylab("Simulated thermal overlap") +
  scale_size(name="1-rank") + scale_colour_discrete(name="")
g_orsim_4 <- ggplot(niche_data2, aes(x = niche_overlap_prec, y = niche_overlap_prec_orsim)) +
  geom_abline(aes(slope=1, intercept=0)) + ggtitle("(d)") + xlim(c(0.05,0.5)) + ylim(c(0.05,0.5)) +
  geom_point(aes(col = species, size = 1-pvals_prec), alpha=0.75) + theme_classic() +
  xlab("Empirical precipitation overlap") + ylab("Simulated precipitation overlap") +
  scale_size(name="1-rank") + scale_colour_discrete(name="")

pdf("results/Fig_3_new3.pdf", width = 8, height = 8)
ggarrange(g_orsim_1, g_orsim_2, g_orsim_3, g_orsim_4, nrow=2, ncol=2, common.legend = TRUE, legend="bottom")
dev.off()


##  Figures S1 — S7: distribution of samples  ##

##  Willow Flycatcher  ##
data_for_analysis_spp <- data_for_analysis %>% filter(species == "WIFL")
species_climate_spp <- species_climate %>% filter(species == species.names[1])
ecoregions_weights_spp <- ecoregions_weights %>% filter(species == species.names[1])
pop.names <- unique(species_climate_spp$population)
# Pop 1: Pacific North West — change to West
# Pop 2: Interior North West — change to "South Rockies"
# Pop 3: East — change to "Central-East"
# Pop 4: South West — change to "South West"
k=1 # Run for all populations (k= 1 to 4)
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_spatvector(data = newmap, col="grey90", fill="grey90") + ylim(c(-5,60)) + xlim(c(-135,-65)) +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"])) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"])) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf(paste0("results/Fig_S1_", k,".pdf"), width = 10, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()

##  Yellow Warbler  ##
data_for_analysis_spp <- data_for_analysis %>% filter(species == "YEWA")
species_climate_spp <- species_climate %>% filter(species == species.names[2])
ecoregions_weights_spp <- ecoregions_weights %>% filter(species == species.names[2])
pop.names <- unique(species_climate_spp$population)
# Pop 1: Alaska — change to NW
# Pop 2: SouthWest — mistake! Change to NWC
# Pop 3: PNW — Change to "South Rockies"
# Pop 4: East — Keep "East"
# Pop 5: MidWest — Change to "Central"
k=1 # Run for all populations (k= 1 to 5)
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_spatvector(data = newmap, col="grey90", fill="grey90") + ylim(c(5,70)) + xlim(c(-165,-57)) +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"])) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"])) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf(paste0("results/Fig_S2_", k,".pdf"), width = 10, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()

##  Wilson's Warbler  ##
data_for_analysis_spp <- data_for_analysis %>% filter(species == "WIWA")
species_climate_spp <- species_climate %>% filter(species == species.names[3])
ecoregions_weights_spp <- ecoregions_weights %>% filter(species == species.names[3])
pop.names <- unique(species_climate_spp$population)
# Pop 1: AK2Alberta — change to North West
# Pop 2: RockyMtn - change to Rockies
# Pop 3: PNW — Change to West
# Pop 4: CoastalCA — Change to South West 1
# Pop 5: Sierra — Change to South West 2
# Pop 6: Eastern — Change to North East
k=1 # Run for all populations (k= 1 to 6)
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_spatvector(data = newmap, col="grey90", fill="grey90") + ylim(c(6,67)) + xlim(c(-159,-60)) +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"])) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"])) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf(paste0("results/Fig_S3_", k,".pdf"), width = 10, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()

##  Common Yellowthroat  ##
data_for_analysis_spp <- data_for_analysis %>% filter(species == "COYE")
species_climate_spp <- species_climate %>% filter(species == species.names[4])
ecoregions_weights_spp <- ecoregions_weights %>% filter(species == species.names[4])
pop.names <- unique(species_climate_spp$population)
# Pop 1: West — change to Rockies—Central
# Pop 2: CA - change to California
# Pop 3: Southwest — Change to South West
# Pop 4: Midwest — Change to Central—East
# Pop 5: NewEngland — Change to East
k=1 # Run for all populations (k= 1 to 5)
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_spatvector(data = newmap, col="grey90", fill="grey90") + ylim(c(8,60)) + xlim(c(-135,-58)) +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"])) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"])) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf(paste0("results/Fig_S4_", k,".pdf"), width = 10, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()

##  American Redstart  ##
data_for_analysis_spp <- data_for_analysis %>% filter(species == "AMRE")
species_climate_spp <- species_climate %>% filter(species == species.names[5])
ecoregions_weights_spp <- ecoregions_weights %>% filter(species == species.names[5])
pop.names <- unique(species_climate_spp$population)
# Pop 1: Northwest — change to North West-Rockies
# Pop 2: South - change to South East
# Pop 3: Northeast — Change to North Central-East
k=1 # Run for all populations (k= 1 to 3)
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_spatvector(data = newmap, col="grey90", fill="grey90") + ylim(c(0,62)) + xlim(c(-135,-58)) +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"])) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"])) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf(paste0("results/Fig_S5_", k,".pdf"), width = 10, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()

##  Painted Bunting  ##
data_for_analysis_spp <- data_for_analysis %>% filter(species == "PABU")
species_climate_spp <- species_climate %>% filter(species == species.names[6])
ecoregions_weights_spp <- ecoregions_weights %>% filter(species == species.names[6])
pop.names <- unique(species_climate_spp$population)
# Pop 1: Central — change to South Central
# Pop 2: Louisiana - change to Louisiana
k=1 # Run for all populations (k= 1 to 2)
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_spatvector(data = newmap, col="grey90", fill="grey90") + ylim(c(5,55)) + xlim(c(-130,-65)) +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"])) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"])) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf(paste0("results/Fig_S6_", k,".pdf"), width = 10, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()

##  Hermit Thrush  ##
data_for_analysis_spp <- data_for_analysis %>% filter(species == "HETH")
species_climate_spp <- species_climate %>% filter(species == species.names[7])
ecoregions_weights_spp <- ecoregions_weights %>% filter(species == species.names[7])
pop.names <- unique(species_climate_spp$population)
# Pop 1: EasternTaiga — change to North West-Central-East
# Pop 2: InteriorWest - change to South Rockies
# Pop 3: PacificCentral — Change to Rockies
# Pop 4: PacificNorth - change to North West
# Pop 5: PacificSouth — Change to South West 2
k=1 # Run for all populations (k= 1 to 5)
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_spatvector(data = newmap, col="grey90", fill="grey90") + ylim(c(10,67)) + xlim(c(-160,-53)) +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="wintering"])) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_spatvector(data = ecoregions[ecoregions_weights_spp$ecoregion[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"],], aes(fill=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"], col=ecoregions_weights_spp$weight[ecoregions_weights_spp$population == pop.names[k] & ecoregions_weights_spp$season=="breeding"])) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf(paste0("results/Fig_S7_", k,".pdf"), width = 11, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()
