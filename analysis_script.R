library(tidyverse)
library(rworldmap)
library(rgeos)
library(sf)
library(terra)
library(fields)
library(ggpubr)

## Get data for sampled individuals
data_for_analysis <- read.csv("results/output/data_for_analysis.csv")
species.names <- unique(data_for_analysis$species)

## Map of the Americas
newmap <- getMap(resolution = "low")
newmap <- spTransform(newmap, '+proj=longlat +datum=WGS84')
newmap@data$world <- rep(1,length(newmap@data$SOVEREIGNT))
newmap <- gBuffer(newmap, byid=TRUE, width=0)
newmap <- gUnaryUnion(newmap, id=newmap@data$world)
newmap <- st_as_sf(newmap)
newmap <- st_crop(newmap, xmin = -180, xmax = -30, ymin = -60, ymax = 90)

## Get ecoregion polygons
ecoregions <- terra::vect("resources/Ecoregions/wwf_terr_ecos.shp")
crs(ecoregions) <- "+proj=longlat +datum=WGS84"
ecoregions <- terra::crop(ecoregions, newmap)

## Spatial locations of samples and UDs
species.names.ebird <- c("wilfly", "yelwar", "wlswar", "comyel", "amered", "paibun")
sp_path <- "seas_abund_rasters"

ecoregions_abund_presence <- list()
for(i in 1:length(species.names)){
  data_for_analysis2 <- data_for_analysis %>% filter(species == species.names[i])
  pops <- unique(data_for_analysis2$population)
  data_for_analysis2_W <- data_for_analysis2 %>% filter(season == "Wintering")
  data_for_analysis2_B <- data_for_analysis2 %>% filter(season == "Breeding")
  
  # Load eBird seasonal abundance surfaces
  
  if(i==6){
    abunds_B <- terra::rast(paste0(sp_path, "/2021/", species.names.ebird[i], "/seasonal/", species.names.ebird[i], "_abundance_seasonal_breeding_mean_2021.tif"))
    abunds_B <- terra::project(abunds_B, '+proj=longlat +datum=WGS84', method = "ngb")
    abunds_B[is.na(abunds_B)] <- 0
    abunds_W <- terra::rast(paste0(sp_path, "/2021/", species.names.ebird[i], "/seasonal/", species.names.ebird[i], "_abundance_seasonal_nonbreeding_mean_2021.tif"))
    abunds_W <- terra::project(abunds_W, '+proj=longlat +datum=WGS84', method = "ngb")
    abunds_W[is.na(abunds_W)] <- 0
  }else{
    abunds <- terra::rast(paste0(sp_path, "/2021/", species.names.ebird[i], "/seasonal/", species.names.ebird[i], "_abundance_seasonal_mean_mr_2021.tif"))
    abunds <- terra::project(abunds, '+proj=longlat +datum=WGS84', method = "ngb")
    abunds[is.na(abunds)] <- 0
    abunds_B <- abunds[[1]]
    abunds_W <- abunds[[2]]
  }
  
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

## Extract climate

years <- 2000:2018
months <- c("01","02","03","04","05","06","07","08","09","10","11","12")
months_days <- c("001", "032", "060", "091", "121", "152", "182", "213", "244", "274", "305", "335")
Temp_files <- list.files("chelsa/tas")
Temp_files <- Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[4])) %in% years]
Prec_files <- list.files("chelsa/pr")
Prec_files <- Prec_files[unlist(lapply(strsplit(Prec_files, "_"), function(x) x[4])) %in% years]
summer <- months[6:7]
summer_days <- months_days[6:7]
winter <- months[c(1:2, 12)]
winter_days <- months_days[c(1:2, 12)]

temp_BR <- matrix(Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[3])) %in% summer], ncol=length(summer))
temp_BR <- terra::rast(apply(temp_BR, 1, function(x) terra::rast(paste0("chelsa/tas/", x)))) # load temperature rasters
prec_BR <- matrix(Prec_files[unlist(lapply(strsplit(Prec_files, "_"), function(x) x[3])) %in% summer], ncol=length(summer))
prec_BR <- terra::rast(apply(prec_BR, 1, function(x) terra::rast(paste0("chelsa/pr/", x)))) # load rasters
temp_NB <- matrix(Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[3])) %in% winter], ncol=length(winter))
temp_NB <- terra::rast(apply(temp_NB, 1, function(x) terra::rast(paste0("chelsa/tas/", x)))) # load temperature rasters
prec_NB <- matrix(Prec_files[unlist(lapply(strsplit(Prec_files, "_"), function(x) x[3])) %in% winter], ncol=length(winter))
prec_NB <- terra::rast(apply(prec_NB, 1, function(x) terra::rast(paste0("chelsa/pr/", x)))) # load rasters
temp_BR_mean <- (mean(temp_BR) / 10) - 273.15
prec_BR_mean <- mean(prec_BR)
temp_NB_mean <- (mean(temp_NB) / 10) - 273.15
prec_NB_mean <- mean(prec_NB)
crs(temp_BR_mean) <- crs(prec_BR_mean) <- crs(temp_NB_mean) <- crs(prec_NB_mean) <- "+proj=longlat +datum=WGS84"

# crop around the Americas
temp_BR_mean_2 <- terra::extract(temp_BR_mean, terra::vect(newmap))
temp_NB_mean_2 <- terra::extract(temp_NB_mean, terra::vect(newmap))
prec_BR_mean_2 <- terra::extract(prec_BR_mean, terra::vect(newmap))
prec_NB_mean_2 <- terra::extract(prec_NB_mean, terra::vect(newmap))

# zscores
temp_mean <- mean(c(temp_BR_mean_2$mean, temp_NB_mean_2$mean), na.rm=T)
temp_sd <- sd(c(temp_BR_mean_2$mean, temp_NB_mean_2$mean), na.rm=T)
prec_mean <- mean(c(prec_BR_mean_2$mean, prec_NB_mean_2$mean), na.rm=T)
prec_sd <- sd(c(prec_BR_mean_2$mean, prec_NB_mean_2$mean), na.rm=T)

temp_zscore_BR <- (temp_BR_mean - temp_mean) / temp_sd
prec_zscore_BR <- (prec_BR_mean - prec_mean) / prec_sd
temp_zscore_NB <- (temp_NB_mean - temp_mean) / temp_sd
prec_zscore_NB <- (prec_NB_mean - prec_mean) / prec_sd

# Check for climate bias in sampling 
ecoregions_W_all_climate <- ecoregions_B_all_climate <- ecoregions_W_samples_climate <- ecoregions_B_samples_climate <- list()
for(k in 1:length(species.names)){
  pops_seas <- do.call(rbind, strsplit(colnames(ecoregions_abund_presence[[k]])[4:ncol(ecoregions_abund_presence[[k]])], "_"))
  pops_names <- unique(pops_seas[,1])
  ecoregions_W_all <- which(ecoregions_abund_presence[[k]]$wintering_abund > 0)
  ecoregions_B_all <- which(ecoregions_abund_presence[[k]]$breeding_abund > 0)
  ecoregions_W_samples <- which(apply(ecoregions_abund_presence[[k]][,4:(3+length(pops_names))], 1, sum) > 0)
  ecoregions_B_samples <- which(apply(ecoregions_abund_presence[[k]][,(4+length(pops_names)):ncol(ecoregions_abund_presence[[k]])], 1, sum) > 0)
  ecoregions_W_all_climate[[k]] <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[ecoregions_W_all,], fun=mean)
  ecoregions_B_all_climate[[k]] <- terra::extract(c(temp_zscore_BR, prec_zscore_BR), ecoregions[ecoregions_B_all,], fun=mean)
  ecoregions_W_samples_climate[[k]] <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[ecoregions_W_samples,], fun=mean)
  ecoregions_B_samples_climate[[k]] <- terra::extract(c(temp_zscore_BR, prec_zscore_BR), ecoregions[ecoregions_B_samples,], fun=mean)
  colnames(ecoregions_W_all_climate[[k]]) <- colnames(ecoregions_B_all_climate[[k]]) <- colnames(ecoregions_W_samples_climate[[k]]) <- colnames(ecoregions_B_samples_climate[[k]]) <- c("ID", "temp", "prec")
}


## Fig S8: Plot breeding temperature and precipitation density of all ecoregions occupied versus sampled ecoregions
k=1
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]])))) 
g_temp_wifl <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(a) WIFL temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_wifl <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(b) WIFL precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=2
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_yewa <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(c) YEWA temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_yewa <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(d) YEWA precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=3
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_wiwa <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(e) WIWA temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_wiwa <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(f) WIWA precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=4
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_coye <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(g) COYE temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_coye <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(h) COYE precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=5
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_amre <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(i) AMRE temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_amre <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(j) AMRE precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=6
ecoregions_clim <- rbind(ecoregions_B_all_climate[[k]], ecoregions_B_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_B_all_climate[[k]])), rep("samples", nrow(ecoregions_B_samples_climate[[k]]))))
g_temp_pabu <- ggplot(ecoregions_clim, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(k) PABU temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_pabu <- ggplot(ecoregions_clim, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(l) PABU precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())

pdf("results/figures/Fig_S8.pdf", width = 8, height = 13)
ggarrange(g_temp_wifl, g_prec_wifl,
          g_temp_yewa, g_prec_yewa,
          g_temp_wiwa, g_prec_wiwa,
          g_temp_coye, g_prec_coye,
          g_temp_amre, g_prec_amre, 
          g_temp_pabu, g_prec_pabu, nrow=6, ncol=2, common.legend = TRUE, legend="bottom")
dev.off()


## Fig S9: Plot wintering temperature and precipitation density of all ecoregions occupied versus sampled ecoregions
k=1
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]])))) 
g_temp_wifl <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(a) WIFL temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_wifl <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(b) WIFL precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=2
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_yewa <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(c) YEWA temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_yewa <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(d) YEWA precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=3
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_wiwa <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(e) WIWA temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_wiwa <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(f) WIWA precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=4
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_coye <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(g) COYE temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_coye <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(h) COYE precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=5
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_amre <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(i) AMRE temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_amre <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(j) AMRE precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())
k=6
ecoregions_W_climate <- rbind(ecoregions_W_all_climate[[k]], ecoregions_W_samples_climate[[k]]) %>%
  mutate(cat = c(rep("all", nrow(ecoregions_W_all_climate[[k]])), rep("samples", nrow(ecoregions_W_samples_climate[[k]]))))
g_temp_pabu <- ggplot(ecoregions_W_climate, aes(x = temp)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-0.8, 1.5)) + ggtitle("(k) PABU temperature") + theme(legend.title=element_blank(), axis.title.x=element_blank())
g_prec_pabu <- ggplot(ecoregions_W_climate, aes(x = prec)) +
  geom_density(aes(fill = cat), adjust = 1, alpha=0.4) +
  xlim(c(-2, 4)) + ggtitle("(l) PABU precipitation") + theme(legend.title=element_blank(), axis.title.x=element_blank())

pdf("results/figures/Fig_S9.pdf", width = 8, height = 13)
ggarrange(g_temp_wifl, g_prec_wifl,
          g_temp_yewa, g_prec_yewa,
          g_temp_wiwa, g_prec_wiwa,
          g_temp_coye, g_prec_coye,
          g_temp_amre, g_prec_amre,
          g_temp_pabu, g_prec_pabu, nrow=6, ncol=2, common.legend = TRUE, legend="bottom")
dev.off()



# Extract climate for each population based on seasonally occupied ecoregions

species_climate <- ecoregions_wintering <- ecoregions_breeding <- list()
for(k in 1:length(species.names)){
  pops_seas <- do.call(rbind, strsplit(colnames(ecoregions_abund_presence[[k]])[4:ncol(ecoregions_abund_presence[[k]])], "_"))
  pops_names <- unique(pops_seas[,1])
  ecoregions_abund_presence_W <- ecoregions_abund_presence[[k]][,4:(3+length(pops_names))]
  ss <- apply(ecoregions_abund_presence_W, 1, sum)
  ss <- ifelse(ss==0, 1, ss)
  ecoregions_abund_presence_W <- apply(ecoregions_abund_presence_W, 2, function(x) x/ss)
  ecoregions_abund_presence_B <- ecoregions_abund_presence[[k]][,(4+length(pops_names)):ncol(ecoregions_abund_presence[[k]])]
  ss <- apply(ecoregions_abund_presence_B, 1, sum)
  ss <- ifelse(ss==0, 1, ss)
  ecoregions_abund_presence_B <- apply(ecoregions_abund_presence_B, 2, function(x) x/ss)
  
  ecoregions_climate_resampled <- ecoregions_W_w <- ecoregions_B_w <- list()
  for(j in 1:length(pops_names)){
    ecoregions_W <- which(ecoregions_abund_presence_W[,j] > 0)
    ecoregions_W_weight <- ecoregions_abund_presence[[k]]$wintering_abund[ecoregions_W] * ecoregions_abund_presence_W[,j][ecoregions_W]
    ecoregions_W_w[[j]] <- cbind(ecoregions_W, ecoregions_W_weight)
    ecoregions_B <- which(ecoregions_abund_presence_B[,j] > 0)
    ecoregions_B_weight <- ecoregions_abund_presence[[k]]$breeding_abund[ecoregions_B] * ecoregions_abund_presence_B[,j][ecoregions_B]
    ecoregions_B_w[[j]] <- cbind(ecoregions_B, ecoregions_B_weight)
    
    # Extract climate for wintering 
    ecoregions_climate <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[ecoregions_W,])
    colnames(ecoregions_climate) <- c("ID", "temp", "prec")
    to_sample <- round((ecoregions_W_weight/sum(ecoregions_W_weight)) * 10000)
    ecoregions_W_climate_resample <- vector()
    for(i in 1:length(to_sample)){
      ecoregions_W_climate_resample <- rbind(ecoregions_W_climate_resample, ecoregions_climate[which(ecoregions_climate$ID == i),][sample(1:length(which(ecoregions_climate$ID == i)), to_sample[i], replace=T),])
    }
    ecoregions_W_climate_resample <- ecoregions_W_climate_resample %>% mutate(population = pops_seas[j,1], season = "wintering")
    
    # Extract climate for breeding 
    ecoregions_climate <- terra::extract(c(temp_zscore_BR, prec_zscore_BR), ecoregions[ecoregions_B,])
    colnames(ecoregions_climate) <- c("ID", "temp", "prec")
    to_sample <- round((ecoregions_B_weight/sum(ecoregions_B_weight)) * 10000)
    ecoregions_B_climate_resample <- vector()
    for(i in 1:length(to_sample)){
      ecoregions_B_climate_resample <- rbind(ecoregions_B_climate_resample, ecoregions_climate[which(ecoregions_climate$ID == i),][sample(1:length(which(ecoregions_climate$ID == i)), to_sample[i], replace=T),])
    }
    ecoregions_B_climate_resample <- ecoregions_B_climate_resample %>% mutate(population = pops_seas[j,1], season = "breeding")
    
    ecoregions_climate_resampled[[j]] <- rbind(ecoregions_W_climate_resample, ecoregions_B_climate_resample)
    
  }
  species_climate[[k]] <- do.call(rbind, ecoregions_climate_resampled) %>%
    mutate(species = species.names[k])
  ecoregions_wintering[[k]] <- ecoregions_W_w
  ecoregions_breeding[[k]] <- ecoregions_B_w
}
species_climate <- do.call(rbind, species_climate)

## Function to estimate seasonal 2-D climatic niche
nicheDensityRaster <- function(seasonalNiche){
  niche.kernel <- kde2d(seasonalNiche[,1], seasonalNiche[,2], n=50, h=1, lims=c(-2,3, -2,3))
  niche.kernel$z = niche.kernel$z/max(niche.kernel$z)
  niche.raster <- raster(niche.kernel)
  threshold=0; i=0
  while(threshold <= 0.95 * sum(niche.kernel$z)){
    i=i+1
    threshold = threshold + sort(as.vector(niche.raster), decreasing=T)[i]
  }
  niche.raster[which(as.vector(niche.raster) < sort(as.vector(niche.raster), decreasing=T)[i])] = 0
  niche.raster = niche.raster / sum(as.vector(niche.raster))
  return(niche.raster)
}

## Calculate seasonal climate niches
# Declare lists to store the niches
breeding.niche <- wintering.niche <- breeding.niche.temp <- wintering.niche.temp <- breeding.niche.prec <- wintering.niche.prec <- vector("list", length=length(species.names))
names(breeding.niche) <- names(wintering.niche) <- names(breeding.niche.temp) <- names(wintering.niche.temp) <- names(breeding.niche.prec) <- names(wintering.niche.prec) <- species.names
pop.names <- list()
# Estimate niches for every population
for(k in 1:length(species.names)){
  species_climate_spp <- species_climate %>% filter(species == species.names[k])
  pop.names[[k]] <- unique(species_climate_spp$population)
  niche_B <- niche_W <- niche_temp_B <- niche_temp_W <- niche_prec_B <- niche_prec_W <- list() #vector("list", length=length(climate_pop_breeding[[k]]))
  for(i in 1:length(pop.names[[k]])){
    niche_B[[i]] <- nicheDensityRaster(species_climate_spp %>% filter(population == pop.names[[k]][[i]] & season == "breeding") %>% dplyr::select(temp, prec))
    niche_W[[i]] <- nicheDensityRaster(species_climate_spp %>% filter(population == pop.names[[k]][[i]] & season == "wintering") %>% dplyr::select(temp, prec))
    niche_temp_B[[i]] <- stats::density(unlist(as.vector(species_climate_spp %>% filter(population == pop.names[[k]][[i]] & season == "breeding") %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_temp_W[[i]] <- stats::density(unlist(as.vector(species_climate_spp %>% filter(population == pop.names[[k]][[i]] & season == "wintering") %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_prec_B[[i]] <- stats::density(unlist(as.vector(species_climate_spp %>% filter(population == pop.names[[k]][[i]] & season == "breeding") %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_prec_W[[i]] <- stats::density(unlist(as.vector(species_climate_spp %>% filter(population == pop.names[[k]][[i]] & season == "wintering") %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
  }
  breeding.niche[[k]] <- niche_B
  wintering.niche[[k]] <- niche_W
  breeding.niche.temp[[k]] <- niche_temp_B
  wintering.niche.temp[[k]] <- niche_temp_W
  breeding.niche.prec[[k]] <- niche_prec_B
  wintering.niche.prec[[k]] <- niche_prec_W
}

## Calculate 2-D niche density
density.breeding <- density.winter <- vector("list", length=length(species.names))
names(density.breeding) <- names(density.winter) <- species.names
for(i in 1:length(species.names)){
  density.breeding[[i]] <- rasterToPoints(raster::stack(breeding.niche[[i]]))
  density.winter[[i]] <- rasterToPoints(raster::stack(wintering.niche[[i]]))
}

## Calculate 2-D niche size
breeding.niche.size <- wintering.niche.size <- vector("list", length=length(species.names))
names(breeding.niche.size) <- names(wintering.niche.size) <- species.names
for(i in 1:length(species.names)){
  breeding.niche.size[[i]] <- apply(density.breeding[[i]][,-c(1,2)], 2, function(x){length(which(x>0))})
  wintering.niche.size[[i]] <- apply(density.winter[[i]][,-c(1,2)], 2, function(x){length(which(x>0))})
}

# Seasonal 2-D niche overlap
niche.overlap <- vector("list", length=length(species.names))
names(niche.overlap) <- species.names
for(i in 1:length(species.names)){
  niche.overlap[[i]] <- 1 - (0.5 * apply(abs(density.breeding[[i]][,-c(1,2)] - density.winter[[i]][,-c(1,2)]), 2, sum))
  names(niche.overlap[[i]]) <- unique((species_climate %>% filter(species == species.names[i]))$population)
}

# Seasonal 1-D niche overlap
niche.overlap.temp <- niche.overlap.prec <- vector("list", length=length(species.names))
names(niche.overlap.temp) <- names(niche.overlap.prec) <- species.names
for(k in 1:length(breeding.niche.temp)){
  niche.overlap.T <- niche.overlap.P <- vector("list", length=length(breeding.niche.temp[[k]]))
  names(niche.overlap.T) <- names(niche.overlap.P) <- pop.names[[k]]
  for(i in 1:length(breeding.niche.temp[[k]])){
    # temperature
    X <- breeding.niche.temp[[k]][[i]]$x
    Y1 <- breeding.niche.temp[[k]][[i]]$y
    Y2 <- wintering.niche.temp[[k]][[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche.overlap.T[[i]]  <- trapz ( X, Overlap ) / Total
    # precipitation
    X <- breeding.niche.prec[[k]][[i]]$x
    Y1 <- breeding.niche.prec[[k]][[i]]$y
    Y2 <- wintering.niche.prec[[k]][[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche.overlap.P[[i]]  <- trapz ( X, Overlap ) / Total
  }
  niche.overlap.temp[[k]] <- niche.overlap.T
  niche.overlap.prec[[k]] <- niche.overlap.P
}

# Population centroids
pop_centroids <- vector()
for(i in 1:length(species.names)){
  pops_seas <- do.call(rbind, strsplit(colnames(ecoregions_abund_presence[[i]])[4:ncol(ecoregions_abund_presence[[i]])], "_"))
  centroid_occupied_ecoregions_W <- centroid_occupied_ecoregions_B <- vector()
  for(j in 1:length(unique(pops_seas[,1]))){
    centroid_occupied_ecoregions_W <- rbind(centroid_occupied_ecoregions_W, apply(matrix(geom(terra::centroids(ecoregions[ecoregions_wintering[[i]][[j]][,1],]))[,3:4], ncol=2), 2, function(x) weighted.mean(x, w=ecoregions_wintering[[i]][[j]][,2])))
    centroid_occupied_ecoregions_B <- rbind(centroid_occupied_ecoregions_B, apply(matrix(geom(terra::centroids(ecoregions[ecoregions_breeding[[i]][[j]][,1],]))[,3:4], ncol=2), 2, function(x) weighted.mean(x, w=ecoregions_breeding[[i]][[j]][,2])))
  }
  pop_centroids <- rbind(pop_centroids, as.data.frame(cbind(species.names[i], pops_seas[,1:2], rbind(centroid_occupied_ecoregions_W, centroid_occupied_ecoregions_B))))
}
colnames(pop_centroids) <- c("species", "population", "season", "longitude", "latitude") 
pop_centroids <- pop_centroids %>% unite("species_pop", species:population, remove=F)

# migration distances
migration_distances <- vector()
for(i in 1:length(unique(pop_centroids$species_pop))){
  ctrs <- pop_centroids %>% filter(species_pop == unique(species_pop)[i]) %>% dplyr::select(longitude, latitude)
  migration_distances[i] <- rdist.earth(matrix(as.numeric(as.matrix(ctrs)), ncol=2, byrow=F), miles=F)[1,2]
}
migration_distances <- data.frame(species_pop = unique(pop_centroids$species_pop),
                                  migration_distance = migration_distances) %>%
  separate(species_pop, c("species", "population"), remove=F)


# Put all the explanatory variables in the same data frame
niche_data <- data.frame(
  species = migration_distances$species,
  population = migration_distances$population,
  niche_overlap = unlist(niche.overlap),
  niche_overlap_temp = unlist(niche.overlap.temp),
  niche_overlap_prec = unlist(niche.overlap.prec),
  breeding_niche_size = unlist(breeding.niche.size),
  wintering_niche_size = unlist(wintering.niche.size),
  migration_distance = migration_distances$migration_distance
)
niche_data <- niche_data %>% unite("species_pop", species:population, remove=F)

#load("results/output/niche_data_new5.RData")

# rename populations
niche_data$population <- c("W", "SR", "CE", "SW",
                           "NW", "NWC", "SR", "E", "C",
                           "NW", "R", "W", "SW1", "SW2", "NE",
                           "RC", "CA", "SW", "CE", "E",
                           "NWR", "SC", "NCE", "C",
                           "SC", "Lsn")

# remove population with poor data 
'%!in%' <- function(x,y)!('%in%'(x,y))
niche_data_2 <- niche_data %>% 
  unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_E", "COYE_CA", "AMRE_C", "PABU_Lsn"))


## Simulating migratory connectivity using ORSIM ##

# write abundance distributions as csv files
seasonalAbundances_B <- data.frame(WIFL = ecoregions_abund_presence[[1]]$breeding_abund,
                                   YEWA = ecoregions_abund_presence[[2]]$breeding_abund,
                                   WIWA = ecoregions_abund_presence[[3]]$breeding_abund,
                                   COYE = ecoregions_abund_presence[[4]]$breeding_abund,
                                   AMRE = ecoregions_abund_presence[[5]]$breeding_abund,
                                   PABU = ecoregions_abund_presence[[6]]$breeding_abund)
seasonalAbundances_W <- data.frame(WIFL = ecoregions_abund_presence[[1]]$wintering_abund,
                                   YEWA = ecoregions_abund_presence[[2]]$wintering_abund,
                                   WIWA = ecoregions_abund_presence[[3]]$wintering_abund,
                                   COYE = ecoregions_abund_presence[[4]]$wintering_abund,
                                   AMRE = ecoregions_abund_presence[[5]]$wintering_abund,
                                   PABU = ecoregions_abund_presence[[6]]$wintering_abund)
seasonalAbundances_B <- apply(seasonalAbundances_B, 2, function(x) ifelse(is.na(x)==T, 0, x))
seasonalAbundances_W <- apply(seasonalAbundances_W, 2, function(x) ifelse(is.na(x)==T, 0, x))
seasonalAbundances_B <- apply(seasonalAbundances_B, 2, function(x) x/sum(x))
seasonalAbundances_W <- apply(seasonalAbundances_W, 2, function(x) x/sum(x))
write.csv(seasonalAbundances_B, "results/output/seasonalAbundances_B.csv", row.names = F)
write.csv(seasonalAbundances_W, "results/output/seasonalAbundances_W.csv", row.names = F)

# write distance matrix as a csv file
distanceMat <- as.matrix(st_distance(st_as_sf(ecoregions)) / 1000)
write.csv(distanceMat, "results/output/distanceMatrix.csv", row.names = F, col.names = F)

# Run orsim.py

# Load orsim results
ORSIM_results <- list()
ORSIM_results[[1]] <- read.csv("results/output/ORSIM_results_WIFL.csv", header=F)
ORSIM_results[[2]] <- read.csv("results/output/ORSIM_results_YEWA.csv", header=F)
ORSIM_results[[3]] <- read.csv("results/output/ORSIM_results_WIWA.csv", header=F)
ORSIM_results[[4]] <- read.csv("results/output/ORSIM_results_COYE.csv", header=F)
ORSIM_results[[5]] <- read.csv("results/output/ORSIM_results_AMRE.csv", header=F)
ORSIM_results[[6]] <- read.csv("results/output/ORSIM_results_PABU.csv", header=F)


# Map ORSIM results and quantify how good they are at capturing the observed patterns
pop_centroids_orsim <- vector()
for(k in 1:length(species.names)){
  wintering_ecoregions <- which(ecoregions_abund_presence[[k]]$wintering_abund > 0)
  breeding_ecoregions <- which(ecoregions_abund_presence[[k]]$breeding_abund > 0)
  pops_seas <- do.call(rbind, strsplit(colnames(ecoregions_abund_presence[[k]])[4:ncol(ecoregions_abund_presence[[k]])], "_"))
  centroid_wintering_ecoregions <- vector()
  for(j in 4:(3 + length(unique(pops_seas[,1])))){
    breeding_pop_ecoregions <- which(ecoregions_abund_presence[[k]][breeding_ecoregions,][,j+length(unique(pops_seas[,1]))] > 0)
    wintering_pop_ecoregions_orsim <- which(apply(ORSIM_results[[k]][breeding_pop_ecoregions,], 2, sum) > 0)
    
    if(length(ecoregions[wintering_ecoregions][wintering_pop_ecoregions_orsim]) > 1){
      centroid_wintering_ecoregions <- rbind(centroid_wintering_ecoregions, 
                                             apply(matrix(geom(terra::centroids(ecoregions[wintering_ecoregions][wintering_pop_ecoregions_orsim]))[,3:4], ncol=2), 2, function(x) weighted.mean(x, w=apply(ORSIM_results[[k]][breeding_pop_ecoregions,][,wintering_pop_ecoregions_orsim], 2, sum)))
      )
    }else{
      centroid_wintering_ecoregions <- rbind(centroid_wintering_ecoregions, 
                                             matrix(geom(terra::centroids(ecoregions[wintering_ecoregions][wintering_pop_ecoregions_orsim]))[,3:4], ncol=2)
      )
    }
  }
  pop_centroids_orsim <- rbind(pop_centroids_orsim, 
                               as.data.frame(cbind(species.names[k], pops_seas[1:length(unique(pops_seas[,1])),1:2], centroid_wintering_ecoregions)) %>% rename(species=V1, population=V2, season=V3, longitude=V4, latitude=V5),
                               pop_centroids %>% filter(species == species.names[k] & season=="breeding") %>% dplyr::select(species:latitude))
}
colnames(pop_centroids_orsim) <- c("species", "population", "season", "longitude", "latitude")
pop_centroids_orsim <- pop_centroids_orsim %>% unite("species_pop", species:population, remove=F)

pop_centroids_sf <- pop_centroids %>% st_as_sf(coords = c("longitude", "latitude"), crs='+proj=longlat +datum=WGS84')
pop_centroids_df <- pop_centroids %>% filter(season=="wintering") %>% left_join(pop_centroids %>% filter(season=="breeding"), by=c("species_pop", "species", "population"))
pop_centroids_B <- pop_centroids %>% filter(season=="breeding")

pop_centroids_orsim_sf <- pop_centroids_orsim %>% st_as_sf(coords = c("longitude", "latitude"), crs='+proj=longlat +datum=WGS84')
pop_centroids_orsim_df <- pop_centroids_orsim %>% filter(season=="wintering") %>% left_join(pop_centroids_orsim %>% filter(season=="breeding"), by=c("species_pop", "species", "population"))
pop_centroids_orsim_B <- pop_centroids_orsim %>% filter(season=="breeding")

pop_centroids_sf <- pop_centroids_sf %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))
pop_centroids_df <- pop_centroids_df %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))
pop_centroids_orsim_sf <- pop_centroids_orsim_sf %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))
pop_centroids_orsim_df <- pop_centroids_orsim_df %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))

##  Extract climate for ORSIM simulations, and calculate niches and migration distances
niche_overlap_null_2 <- niche_overlap_temp_null_2 <- niche_overlap_prec_null_2 <- wintering_niche_size_null_2 <- list()
for(k in 1:length(species.names)){
  wintering_ecoregions <- which(ecoregions_abund_presence[[k]]$wintering_abund > 0)
  breeding_ecoregions <- which(ecoregions_abund_presence[[k]]$breeding_abund > 0)
  pops_seas <- do.call(rbind, strsplit(colnames(ecoregions_abund_presence[[k]])[4:ncol(ecoregions_abund_presence[[k]])], "_"))
  niche_W_null <- niche_temp_W_null <- niche_prec_W_null <- list()
  for(j in 4:(3+length(unique(pops_seas[,1])))){
    breeding_pop_ecoregions <- which(ecoregions_abund_presence[[k]][breeding_ecoregions,][,j+length(unique(pops_seas[,1]))] > 0)
    wintering_pop_ecoregions_orsim <- which(apply(ORSIM_results[[k]][breeding_pop_ecoregions,], 2, sum) > 0)
    ecoregions_climate <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_pop_ecoregions_orsim,])
    colnames(ecoregions_climate) <- c("ID", "temp", "prec")
    if(length(breeding_pop_ecoregions) == 1 & length(wintering_pop_ecoregions_orsim)==1){
      ecoregions_weights <- ORSIM_results[[k]][breeding_pop_ecoregions,][,wintering_pop_ecoregions_orsim]
    }else{
      ecoregions_weights <- apply(ORSIM_results[[k]][breeding_pop_ecoregions,][,wintering_pop_ecoregions_orsim], 2, sum)
    }
    to_sample <- round((ecoregions_weights/sum(ecoregions_weights)) * 10000)
    ecoregions_W_climate_resample <- vector()
    for(h in 1:length(to_sample)){
      ecoregions_W_climate_resample <- rbind(ecoregions_W_climate_resample, ecoregions_climate[which(ecoregions_climate$ID == h),][sample(1:length(which(ecoregions_climate$ID == h)), to_sample[h], replace=T),])
    }
    ecoregions_W_climate_resample <- ecoregions_W_climate_resample %>% mutate(population = pops_seas[j-3,1])
    niche_W_null[[j-3]] <- nicheDensityRaster(ecoregions_W_climate_resample %>% dplyr::select(temp, prec))
    niche_temp_W_null[[j-3]] <- stats::density(unlist(as.vector(ecoregions_W_climate_resample %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    niche_prec_W_null[[j-3]] <- stats::density(unlist(as.vector(ecoregions_W_climate_resample %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
  }
  density_W_null <- rasterToPoints(raster::stack(niche_W_null))
  wintering_niche_size_null_2[[k]] <- apply(density_W_null[,-c(1,2)], 2, function(x){length(which(x>0))})
  niche_overlap_null_2[[k]] <- 1 - (0.5 * apply(abs(density.breeding[[k]][,-c(1,2)] - density_W_null[,-c(1,2)]), 2, sum))
  # Seasonal 1-D niche overlap
  niche.overlap.T <- niche.overlap.P <- vector()
  for(i in 1:length(unique(pops_seas[,1]))){
    # temperature
    X <- breeding.niche.temp[[k]][[i]]$x
    Y1 <- breeding.niche.temp[[k]][[i]]$y
    Y2 <- niche_temp_W_null[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche.overlap.T[i]  <- trapz ( X, Overlap ) / Total
    # precipitation
    X <- breeding.niche.prec[[k]][[i]]$x
    Y1 <- breeding.niche.prec[[k]][[i]]$y
    Y2 <- niche_prec_W_null[[i]]$y
    Overlap <- pmin ( Y1, Y2 )
    Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
    niche.overlap.P[i] <- trapz ( X, Overlap ) / Total
  }
  niche_overlap_temp_null_2[[k]] <- niche.overlap.T
  niche_overlap_prec_null_2[[k]] <- niche.overlap.P
}
#load("results/output/null_model_2_results_new3.RData")

# Migration distances predicted by ORSIM
migration_distance_orsim <- diag(rdist.earth(apply(as.matrix(pop_centroids_orsim_df[,5:6]), 2, as.numeric), 
                                             apply(as.matrix(pop_centroids_orsim_df[,8:9]), 2, as.numeric), 
                                             miles=F))
niche_data_null_orsim <- cbind(niche_data, unlist(niche_overlap_null_2), unlist(niche_overlap_temp_null_2), unlist(niche_overlap_prec_null_2))
colnames(niche_data_null_orsim) <- c(colnames(niche_data), "niche_overlap_orsim", "niche_overlap_temp_orsim", "niche_overlap_prec_orsim")
niche_data_null_orsim_2 <- niche_data_null_orsim %>% 
  unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_E", "COYE_CA", "AMRE_C", "PABU_Lsn")) %>%
  mutate(migration_distance_orsim = migration_distance_orsim)



###  Climate tracking simulation  ###

# extract climate for all ecoregions
ecoregions_climate_NB <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions, fun="mean", na.rm=T)
ecoregions_climate_BR <- terra::extract(c(temp_zscore_BR, prec_zscore_BR), ecoregions, fun="mean", na.rm=T)
colnames(ecoregions_climate_NB) <- colnames(ecoregions_climate_BR) <- c("ID", "temperature", "precipitation")

# calculate and write climate distance matrices as csv files
climate_distanceMat <- rdist(as.matrix(ecoregions_climate_BR[,2:3]), as.matrix(ecoregions_climate_NB[,2:3]))
thermal_distanceMat <- abs(outer(ecoregions_climate_BR$temperature, ecoregions_climate_NB$temperature, "-"))
precipitation_distanceMat <- abs(outer(ecoregions_climate_BR$precipitation, ecoregions_climate_NB$precipitation, "-"))
write.csv(climate_distanceMat, "results/output/climate_distanceMatrix.csv", row.names = F, col.names = F)
write.csv(thermal_distanceMat, "results/output/thermal_distanceMatrix.csv", row.names = F, col.names = F)
write.csv(precipitation_distanceMat, "results/output/precipitation_distanceMatrix.csv", row.names = F, col.names = F)

# Optimal climate tracking
overlap_fct <- function(X, Y1, Y2){
  Overlap <- pmin ( Y1, Y2 )
  Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
  return(trapz ( X, Overlap ) / Total)
}
pop_centroids_climate <- pop_centroids_thermal <- pop_centroids_precipitation <- vector()
niche_overlap_climate <- niche_overlap_temp_climate <- niche_overlap_prec_climate <- wintering_niche_size_climate <- list()
niche_overlap_thermal <- niche_overlap_temp_thermal <- niche_overlap_prec_thermal <- wintering_niche_size_thermal <- list()
niche_overlap_precipitation <- niche_overlap_temp_precipitation <- niche_overlap_prec_precipitation <- wintering_niche_size_precipitation <- list()
for(k in 1:length(species.names)){
  wintering_ecoregions <- which(ecoregions_abund_presence[[k]]$wintering_abund > 0)
  breeding_ecoregions <- which(ecoregions_abund_presence[[k]]$breeding_abund > 0)
  climate_distanceMat_2 <- climate_distanceMat[breeding_ecoregions, wintering_ecoregions]
  thermal_distanceMat_2 <- thermal_distanceMat[breeding_ecoregions, wintering_ecoregions]
  precipitation_distanceMat_2 <- precipitation_distanceMat[breeding_ecoregions, wintering_ecoregions]
  
  pops_seas <- do.call(rbind, strsplit(colnames(ecoregions_abund_presence[[k]])[4:ncol(ecoregions_abund_presence[[k]])], "_"))
  centroid_wintering_ecoregions_climate <- centroid_wintering_ecoregions_thermal <- centroid_wintering_ecoregions_precipitation <- vector()
  clim_W_null_climate <- temp_W_null_climate <- prec_W_null_climate <- list()
  clim_W_null_thermal <- temp_W_null_thermal <- prec_W_null_thermal <- list()
  clim_W_null_precipitation <- temp_W_null_precipitation <- prec_W_null_precipitation <- list()
  for(j in 4:(3 + length(unique(pops_seas[,1])))){
    breeding_pop_ecoregions <- which(ecoregions_abund_presence[[k]][breeding_ecoregions,][,j+length(unique(pops_seas[,1]))] > 0)
    if(length(breeding_pop_ecoregions) > 1){
      wintering_pop_ecoregions_climate <- apply(climate_distanceMat_2[breeding_pop_ecoregions,], 1, which.min)
      wintering_pop_ecoregions_thermal <- apply(thermal_distanceMat_2[breeding_pop_ecoregions,], 1, which.min)
      wintering_pop_ecoregions_precipitation <- apply(precipitation_distanceMat_2[breeding_pop_ecoregions,], 1, which.min)
    }else{
      wintering_pop_ecoregions_climate <- which.min(climate_distanceMat_2[breeding_pop_ecoregions,])
      wintering_pop_ecoregions_thermal <- which.min(thermal_distanceMat_2[breeding_pop_ecoregions,])
      wintering_pop_ecoregions_precipitation <- which.min(precipitation_distanceMat_2[breeding_pop_ecoregions,])
    }
    centroid_wintering_ecoregions_climate <- rbind(centroid_wintering_ecoregions_climate, 
                                                   apply(matrix(geom(terra::centroids(ecoregions[wintering_ecoregions][wintering_pop_ecoregions_climate]))[,3:4], ncol=2), 2, mean))
    centroid_wintering_ecoregions_thermal <- rbind(centroid_wintering_ecoregions_thermal, 
                                                   apply(matrix(geom(terra::centroids(ecoregions[wintering_ecoregions][wintering_pop_ecoregions_thermal]))[,3:4], ncol=2), 2, mean))
    centroid_wintering_ecoregions_precipitation <- rbind(centroid_wintering_ecoregions_precipitation, 
                                                         apply(matrix(geom(terra::centroids(ecoregions[wintering_ecoregions][wintering_pop_ecoregions_precipitation]))[,3:4], ncol=2), 2, mean))
    
    ecoregions_clim_climate <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_pop_ecoregions_climate,])
    ecoregions_clim_thermal <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_pop_ecoregions_thermal,])
    ecoregions_clim_precipitation <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_pop_ecoregions_precipitation,])
    colnames(ecoregions_clim_climate) <- colnames(ecoregions_clim_thermal) <- colnames(ecoregions_clim_precipitation) <- c("ID", "temp", "prec")
    clim_W_null_climate[[j-3]] <- nicheDensityRaster(ecoregions_clim_climate %>% dplyr::select(temp, prec))
    temp_W_null_climate[[j-3]] <- stats::density(unlist(as.vector(ecoregions_clim_climate %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    prec_W_null_climate[[j-3]] <- stats::density(unlist(as.vector(ecoregions_clim_climate %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    clim_W_null_thermal[[j-3]] <- nicheDensityRaster(ecoregions_clim_thermal %>% dplyr::select(temp, prec))
    temp_W_null_thermal[[j-3]] <- stats::density(unlist(as.vector(ecoregions_clim_thermal %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    prec_W_null_thermal[[j-3]] <- stats::density(unlist(as.vector(ecoregions_clim_thermal %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    clim_W_null_precipitation[[j-3]] <- nicheDensityRaster(ecoregions_clim_precipitation %>% dplyr::select(temp, prec))
    temp_W_null_precipitation[[j-3]] <- stats::density(unlist(as.vector(ecoregions_clim_precipitation %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
    prec_W_null_precipitation[[j-3]] <- stats::density(unlist(as.vector(ecoregions_clim_precipitation %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
  }
  
  pop_centroids_climate <- rbind(pop_centroids_climate, 
                                 as.data.frame(cbind(species.names[k], pops_seas[1:length(unique(pops_seas[,1])),1:2], centroid_wintering_ecoregions_climate)) %>% rename(species=V1, population=V2, season=V3, longitude=V4, latitude=V5),
                                 pop_centroids %>% filter(species == species.names[k] & season=="breeding") %>% dplyr::select(species:latitude))
  pop_centroids_thermal <- rbind(pop_centroids_thermal, 
                                 as.data.frame(cbind(species.names[k], pops_seas[1:length(unique(pops_seas[,1])),1:2], centroid_wintering_ecoregions_thermal)) %>% rename(species=V1, population=V2, season=V3, longitude=V4, latitude=V5),
                                 pop_centroids %>% filter(species == species.names[k] & season=="breeding") %>% dplyr::select(species:latitude))
  pop_centroids_precipitation <- rbind(pop_centroids_precipitation, 
                                       as.data.frame(cbind(species.names[k], pops_seas[1:length(unique(pops_seas[,1])),1:2], centroid_wintering_ecoregions_precipitation)) %>% rename(species=V1, population=V2, season=V3, longitude=V4, latitude=V5),
                                       pop_centroids %>% filter(species == species.names[k] & season=="breeding") %>% dplyr::select(species:latitude))
  
  density_W_null_climate <- rasterToPoints(raster::stack(clim_W_null_climate))
  density_W_null_thermal <- rasterToPoints(raster::stack(clim_W_null_thermal))
  density_W_null_precipitation <- rasterToPoints(raster::stack(clim_W_null_precipitation))
  wintering_niche_size_climate[[k]] <- apply(density_W_null_climate[,-c(1,2)], 2, function(x){length(which(x>0))})
  wintering_niche_size_thermal[[k]] <- apply(density_W_null_thermal[,-c(1,2)], 2, function(x){length(which(x>0))})
  wintering_niche_size_precipitation[[k]] <- apply(density_W_null_precipitation[,-c(1,2)], 2, function(x){length(which(x>0))})
  niche_overlap_climate[[k]] <- 1 - (0.5 * apply(abs(density.breeding[[k]][,-c(1,2)] - density_W_null_climate[,-c(1,2)]), 2, sum))
  niche_overlap_thermal[[k]] <- 1 - (0.5 * apply(abs(density.breeding[[k]][,-c(1,2)] - density_W_null_thermal[,-c(1,2)]), 2, sum))
  niche_overlap_precipitation[[k]] <- 1 - (0.5 * apply(abs(density.breeding[[k]][,-c(1,2)] - density_W_null_precipitation[,-c(1,2)]), 2, sum))
  # Seasonal 1-D niche overlap
  niche_overlap_T_climate <- niche_overlap_P_climate <- niche_overlap_T_thermal <- niche_overlap_P_thermal <- niche_overlap_T_precipitation <- niche_overlap_P_precipitation <- vector()
  for(i in 1:length(unique(pops_seas[,1]))){
    niche_overlap_T_climate[i] <- overlap_fct(breeding.niche.temp[[k]][[i]]$x, breeding.niche.temp[[k]][[i]]$y, temp_W_null_climate[[i]]$y)
    niche_overlap_P_climate[i] <- overlap_fct(breeding.niche.prec[[k]][[i]]$x, breeding.niche.prec[[k]][[i]]$y, prec_W_null_climate[[i]]$y)
    niche_overlap_T_thermal[i] <- overlap_fct(breeding.niche.temp[[k]][[i]]$x, breeding.niche.temp[[k]][[i]]$y, temp_W_null_thermal[[i]]$y)
    niche_overlap_P_thermal[i] <- overlap_fct(breeding.niche.prec[[k]][[i]]$x, breeding.niche.prec[[k]][[i]]$y, prec_W_null_thermal[[i]]$y)
    niche_overlap_T_precipitation[i] <- overlap_fct(breeding.niche.temp[[k]][[i]]$x, breeding.niche.temp[[k]][[i]]$y, temp_W_null_precipitation[[i]]$y)
    niche_overlap_P_precipitation[i] <- overlap_fct(breeding.niche.prec[[k]][[i]]$x, breeding.niche.prec[[k]][[i]]$y, prec_W_null_precipitation[[i]]$y)
  }
  niche_overlap_temp_climate[[k]] <- niche_overlap_T_climate
  niche_overlap_prec_climate[[k]] <- niche_overlap_P_climate
  niche_overlap_temp_thermal[[k]] <- niche_overlap_T_thermal
  niche_overlap_prec_thermal[[k]] <- niche_overlap_P_thermal
  niche_overlap_temp_precipitation[[k]] <- niche_overlap_T_precipitation
  niche_overlap_prec_precipitation[[k]] <- niche_overlap_P_precipitation
}
#load("results/output/null_model_2_results_climate_thermal_precipitation.RData")

pop_centroids_climate <- pop_centroids_climate %>% unite("species_pop", species:population, remove=F)
pop_centroids_thermal <- pop_centroids_thermal %>% unite("species_pop", species:population, remove=F)
pop_centroids_precipitation <- pop_centroids_precipitation %>% unite("species_pop", species:population, remove=F)

##  Plot connectivity map of population migration simulated by ORSIM climate ##
pop_centroids_sf <- pop_centroids %>% st_as_sf(coords = c("longitude", "latitude"), crs='+proj=longlat +datum=WGS84')
pop_centroids_df <- pop_centroids %>% filter(season=="wintering") %>% left_join(pop_centroids %>% filter(season=="breeding"), by=c("species_pop", "species", "population"))
pop_centroids_B <- pop_centroids %>% filter(season=="breeding")
pop_centroids_climate_sf <- pop_centroids_climate %>% st_as_sf(coords = c("longitude", "latitude"), crs='+proj=longlat +datum=WGS84')
pop_centroids_climate_df <- pop_centroids_climate %>% filter(season=="wintering") %>% left_join(pop_centroids_climate %>% filter(season=="breeding"), by=c("species_pop", "species", "population"))
pop_centroids_climate_B <- pop_centroids_climate %>% filter(season=="breeding")
pop_centroids_thermal_sf <- pop_centroids_thermal %>% st_as_sf(coords = c("longitude", "latitude"), crs='+proj=longlat +datum=WGS84')
pop_centroids_thermal_df <- pop_centroids_thermal %>% filter(season=="wintering") %>% left_join(pop_centroids_thermal %>% filter(season=="breeding"), by=c("species_pop", "species", "population"))
pop_centroids_thermal_B <- pop_centroids_thermal %>% filter(season=="breeding")
pop_centroids_precipitation_sf <- pop_centroids_precipitation %>% st_as_sf(coords = c("longitude", "latitude"), crs='+proj=longlat +datum=WGS84')
pop_centroids_precipitation_df <- pop_centroids_precipitation %>% filter(season=="wintering") %>% left_join(pop_centroids_precipitation %>% filter(season=="breeding"), by=c("species_pop", "species", "population"))
pop_centroids_precipitation_B <- pop_centroids_precipitation %>% filter(season=="breeding")

pop_centroids_sf <- pop_centroids_sf %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))
pop_centroids_df <- pop_centroids_df %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))
pop_centroids_climate_sf <- pop_centroids_climate_sf %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))
pop_centroids_climate_df <- pop_centroids_climate_df %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))
pop_centroids_thermal_sf <- pop_centroids_thermal_sf %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))
pop_centroids_thermal_df <- pop_centroids_thermal_df %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))
pop_centroids_precipitation_sf <- pop_centroids_precipitation_sf %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))
pop_centroids_precipitation_df <- pop_centroids_precipitation_df %>% unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_East", "COYE_CA", "AMRE_SouthDakota", "PABU_Louisiana"))



# Migration distance predicted by climate tracking models
migration_distance_climate <- diag(rdist.earth(apply(as.matrix(pop_centroids_climate_df[,5:6]), 2, as.numeric), 
                                               apply(as.matrix(pop_centroids_climate_df[,8:9]), 2, as.numeric), 
                                               miles=F))
migration_distance_thermal <- diag(rdist.earth(apply(as.matrix(pop_centroids_thermal_df[,5:6]), 2, as.numeric), 
                                               apply(as.matrix(pop_centroids_thermal_df[,8:9]), 2, as.numeric), 
                                               miles=F))
migration_distance_precipitation <- diag(rdist.earth(apply(as.matrix(pop_centroids_precipitation_df[,5:6]), 2, as.numeric), 
                                                     apply(as.matrix(pop_centroids_precipitation_df[,8:9]), 2, as.numeric), 
                                                     miles=F))

niche_data_climate <- cbind(niche_data, 
                            unlist(niche_overlap_climate), unlist(niche_overlap_temp_climate), unlist(niche_overlap_prec_climate),
                            unlist(niche_overlap_thermal), unlist(niche_overlap_temp_thermal), unlist(niche_overlap_prec_thermal),
                            unlist(niche_overlap_precipitation), unlist(niche_overlap_temp_precipitation), unlist(niche_overlap_prec_precipitation))
colnames(niche_data_climate) <- c(colnames(niche_data),
                                  "niche_overlap_climate", "niche_overlap_temp_climate", "niche_overlap_prec_climate",
                                  "niche_overlap_thermal", "niche_overlap_temp_thermal", "niche_overlap_prec_thermal",
                                  "niche_overlap_precipitation", "niche_overlap_temp_precipitation", "niche_overlap_prec_precipitation")

niche_data_climate <- niche_data_climate %>% 
  unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_E", "COYE_CA", "AMRE_C", "PABU_Lsn")) %>%
  mutate(migration_distance_climate = migration_distance_climate,
         migration_distance_thermal = migration_distance_thermal,
         migration_distance_precipitation = migration_distance_precipitation)

niche_data_all_climate <- cbind(niche_data_all, niche_data_climate[,10:21])


##  Fig. 1: niche tracking patterns: observed and simulated  ##

g_conn <- ggplot() + geom_sf(data = newmap) + ylim(c(0,70)) + xlim(c(-158,-53)) +
  geom_sf(data = pop_centroids_sf, aes(col=species, shape=season)) + 
  geom_segment(data = pop_centroids_df, aes(x=as.numeric(longitude.x), y=as.numeric(latitude.x), xend=as.numeric(longitude.y), yend=as.numeric(latitude.y), col=species)) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  ggtitle("(a) Empirical migratory connectivity") + theme_void()

g_conn_orsim <- ggplot() + geom_sf(data = newmap) + ylim(c(0,70)) + xlim(c(-158,-53)) +
  geom_sf(data = pop_centroids_orsim_sf, aes(col=species, shape=season)) + 
  geom_segment(data = pop_centroids_orsim_df, aes(x=as.numeric(longitude.x), y=as.numeric(latitude.x), xend=as.numeric(longitude.y), yend=as.numeric(latitude.y), col=species)) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  ggtitle("(e) Simulated connectivity (ORSIM)") + theme_void()

g_conn_climate <- ggplot() + geom_sf(data = newmap) + ylim(c(0,70)) + xlim(c(-158,-53)) +
  geom_sf(data = pop_centroids_climate_sf, aes(col=species, shape=season)) + 
  geom_segment(data = pop_centroids_climate_df, aes(x=as.numeric(longitude.x), y=as.numeric(latitude.x), xend=as.numeric(longitude.y), yend=as.numeric(latitude.y), col=species)) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  ggtitle("(i) Simulated connectivity (climate)") + theme_void()

g_conn_thermal <- ggplot() + geom_sf(data = newmap) + ylim(c(0,70)) + xlim(c(-158,-53)) +
  geom_sf(data = pop_centroids_thermal_sf, aes(col=species, shape=season)) + 
  geom_segment(data = pop_centroids_thermal_df, aes(x=as.numeric(longitude.x), y=as.numeric(latitude.x), xend=as.numeric(longitude.y), yend=as.numeric(latitude.y), col=species)) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  ggtitle("(i) Simulated connectivity (thermal)") + theme_void()

g_conn_precipitation <- ggplot() + geom_sf(data = newmap) + ylim(c(0,70)) + xlim(c(-158,-53)) +
  geom_sf(data = pop_centroids_precipitation_sf, aes(col=species, shape=season)) + 
  geom_segment(data = pop_centroids_precipitation_df, aes(x=as.numeric(longitude.x), y=as.numeric(latitude.x), xend=as.numeric(longitude.y), yend=as.numeric(latitude.y), col=species)) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  ggtitle("(m) Simulated connectivity (precipitation)") + theme_void()

g_overlap_dist <- ggplot(niche_data_all_climate, aes(x = migration_distance, y = niche_overlap, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.95)) + ggtitle("(b)") + theme_classic() +
  xlab("Empirical migration distance (km)") + ylab("Empirical 2D climate overlap")

g_overlap_temp_dist <- ggplot(niche_data_all_climate, aes(x = migration_distance, y = niche_overlap_temp, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.6)) + ggtitle("(c)") + theme_classic() +
  xlab("Empirical migration distance (km)") + ylab("Empirical thermal overlap")

g_overlap_prec_dist <- ggplot(niche_data_all_climate, aes(x = migration_distance, y = niche_overlap_prec, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.6)) + ggtitle("(d)") + theme_classic() +
  xlab("Empirical migration distance (km)") + ylab("Empirical precipitation overlap")

g_overlap_dist_orsim <- ggplot(niche_data_all_climate, aes(x = migration_distance_orsim, y = niche_overlap_orsim, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.95)) + ggtitle("(f)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated 2D climate overlap (ORSIM)")

g_overlap_temp_dist_orsim <- ggplot(niche_data_all_climate, aes(x = migration_distance_orsim, y = niche_overlap_temp_orsim, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.6)) + ggtitle("(g)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated thermal overlap (ORSIM)")

g_overlap_prec_dist_orsim <- ggplot(niche_data_all_climate, aes(x = migration_distance_orsim, y = niche_overlap_prec_orsim, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.6)) + ggtitle("(h)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated precipitation overlap (ORSIM)")

g_overlap_dist_climate <- ggplot(niche_data_all_climate, aes(x = migration_distance_climate, y = niche_overlap_climate, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.95)) + ggtitle("(j)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated 2D climate overlap (climate)")

g_overlap_temp_dist_climate <- ggplot(niche_data_all_climate, aes(x = migration_distance_climate, y = niche_overlap_temp_climate, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.6)) + ggtitle("(k)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated thermal overlap (climate)")

g_overlap_prec_dist_climate <- ggplot(niche_data_all_climate, aes(x = migration_distance_climate, y = niche_overlap_prec_climate, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.6)) + ggtitle("(l)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated precipitation overlap (climate)")

g_overlap_dist_thermal <- ggplot(niche_data_all_climate, aes(x = migration_distance_thermal, y = niche_overlap_thermal, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.95)) + ggtitle("(j)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated 2D climate overlap (thermal)")

g_overlap_temp_dist_thermal <- ggplot(niche_data_all_climate, aes(x = migration_distance_thermal, y = niche_overlap_temp_thermal, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.6)) + ggtitle("(k)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated thermal overlap (thermal)")

g_overlap_prec_dist_thermal <- ggplot(niche_data_all_climate, aes(x = migration_distance_thermal, y = niche_overlap_prec_thermal, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.6)) + ggtitle("(l)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated precipitation overlap (thermal)")

g_overlap_dist_precipitation <- ggplot(niche_data_all_climate, aes(x = migration_distance_precipitation, y = niche_overlap_precipitation, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.95)) + ggtitle("(n)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated 2D climate overlap (precipitation)")

g_overlap_temp_dist_precipitation <- ggplot(niche_data_all_climate, aes(x = migration_distance_precipitation, y = niche_overlap_temp_precipitation, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.6)) + ggtitle("(o)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated thermal overlap (precipitation)")

g_overlap_prec_dist_precipitation <- ggplot(niche_data_all_climate, aes(x = migration_distance_precipitation, y = niche_overlap_prec_precipitation, label = population)) +
  geom_point(aes(col = species), size=1.5) + 
  geom_text(hjust=0.5, vjust=-0.5, aes(col = species), size=3) + 
  geom_smooth(se = FALSE, col="black", span = 1) +
  ylim(c(0,0.6)) + ggtitle("(p)") + theme_classic() +
  xlab("Simulated migration distance (km)") + ylab("Simulated precipitation overlap (precipitation)")

# Figure 1
pdf("results/figures/Fig_1_new_with_orsim_and_climate.pdf", width = 13, height = 10)
ggarrange(g_conn, g_overlap_dist, g_overlap_temp_dist, g_overlap_prec_dist, g_conn_orsim, g_overlap_dist_orsim, g_overlap_temp_dist_orsim, g_overlap_prec_dist_orsim, g_conn_climate, g_overlap_dist_climate, g_overlap_temp_dist_climate, g_overlap_prec_dist_climate, ncol=4, nrow=3, common.legend = TRUE, legend="bottom")
dev.off()

# Figure S7
pdf("results/figures/Fig_S7_new_with_orsim_and_thermal_precipitation.pdf", width = 13, height = 13)
ggarrange(g_conn, g_overlap_dist, g_overlap_temp_dist, g_overlap_prec_dist, g_conn_orsim, g_overlap_dist_orsim, g_overlap_temp_dist_orsim, g_overlap_prec_dist_orsim, g_conn_thermal, g_overlap_dist_thermal, g_overlap_temp_dist_thermal, g_overlap_prec_dist_thermal, g_conn_precipitation, g_overlap_dist_precipitation, g_overlap_temp_dist_precipitation, g_overlap_prec_dist_precipitation, ncol=4, nrow=4, common.legend = TRUE, legend="bottom")
dev.off()


# Does ORSIM predict niche overlap?
cor_niche <- cor.test(niche_data_all_climate$niche_overlap, niche_data_all_climate$niche_overlap_orsim) 
cor_temp <- cor.test(niche_data_all_climate$niche_overlap_temp, niche_data_all_climate$niche_overlap_temp_orsim) 
cor_prec <- cor.test(niche_data_all_climate$niche_overlap_prec, niche_data_all_climate$niche_overlap_prec_orsim) 

cor_niche2 <- cor.test(niche_data_all_climate$niche_overlap, niche_data_all_climate$niche_overlap_climate) 
cor_temp2 <- cor.test(niche_data_all_climate$niche_overlap_temp, niche_data_all_climate$niche_overlap_temp_climate)
cor_prec2 <- cor.test(niche_data_all_climate$niche_overlap_prec, niche_data_all_climate$niche_overlap_prec_climate)

cor_niche3 <- cor.test(niche_data_all_climate$niche_overlap, niche_data_all_climate$niche_overlap_thermal) 
cor_temp3 <- cor.test(niche_data_all_climate$niche_overlap_temp, niche_data_all_climate$niche_overlap_temp_thermal) 
cor_prec3 <- cor.test(niche_data_all_climate$niche_overlap_prec, niche_data_all_climate$niche_overlap_prec_thermal) 

cor_niche4 <- cor.test(niche_data_all_climate$niche_overlap, niche_data_all_climate$niche_overlap_precipitation)
cor_temp4 <- cor.test(niche_data_all_climate$niche_overlap_temp, niche_data_all_climate$niche_overlap_temp_precipitation)
cor_prec4 <- cor.test(niche_data_all_climate$niche_overlap_prec, niche_data_all_climate$niche_overlap_prec_precipitation)


##  Are residuals skewed towards observed > simulated? ##
orsim_niche_tracking <- niche_data_null_orsim_2$niche_overlap - niche_data_null_orsim_2$niche_overlap_orsim
orsim_temp_tracking <- niche_data_null_orsim_2$niche_overlap_temp - niche_data_null_orsim_2$niche_overlap_temp_orsim
orsim_prec_tracking <- niche_data_null_orsim_2$niche_overlap_prec - niche_data_null_orsim_2$niche_overlap_prec_orsim
hist(orsim_niche_tracking, breaks=10)
hist(orsim_temp_tracking, breaks=10)
hist(orsim_prec_tracking, breaks=10)
ks.test(orsim_niche_tracking,"pnorm", mean=0, sd=sd(orsim_niche_tracking), alternative="less")
ks.test(orsim_temp_tracking,"pnorm", mean=0, sd=sd(orsim_temp_tracking), alternative="less")
ks.test(orsim_prec_tracking,"pnorm", mean=0, sd=sd(orsim_prec_tracking), alternative="less")



###  Null model randomizing the seasonal grounds after ORSIM

niche_overlap_null_4 <- migration_distances_null_4 <- niche_overlap_temp_null_4 <- niche_overlap_prec_null_4 <- wintering_niche_size_null_4 <- list()
for(k in 1:length(species.names)){
  
  pops_seas <- do.call(rbind, strsplit(colnames(ecoregions_abund_presence[[k]])[4:ncol(ecoregions_abund_presence[[k]])], "_"))
  pops_names <- unique(pops_seas[,1])
  ecoregions_abund_presence_W <- ecoregions_abund_presence[[k]][,4:(3+length(pops_names))]
  ss <- apply(ecoregions_abund_presence_W, 1, sum)
  ss <- ifelse(ss==0, 1, ss)
  ecoregions_abund_presence_W <- apply(ecoregions_abund_presence_W, 2, function(x) x/ss)
  ecoregions_abund_presence_B <- ecoregions_abund_presence[[k]][,(4+length(pops_names)):ncol(ecoregions_abund_presence[[k]])]
  ss <- apply(ecoregions_abund_presence_B, 1, sum)
  ss <- ifelse(ss==0, 1, ss)
  ecoregions_abund_presence_B <- apply(ecoregions_abund_presence_B, 2, function(x) x/ss)
  
  # Randomizing wintering grounds by selecting ecoregions within the wintering range of the species and around the ORSIM prediction
  species_winter_ecoregions <- which(ecoregions_abund_presence[[k]]$wintering_abund > 0)
  ctrs_W_spp <- matrix(geom(terra::centroids(ecoregions[species_winter_ecoregions,]))[,3:4], ncol=2)
  species_winter_ecoregions_dists <- rdist.earth(apply(as.matrix(ctrs_W_spp), 2, as.numeric), miles=F)
  ctrs_W_orsim <- pop_centroids_orsim %>% filter(season == "wintering" & species == species.names[k]) %>% dplyr::select(longitude, latitude)
  ctrs_W_orsim_dists <- rdist.earth(apply(as.matrix(ctrs_W_orsim), 2, as.numeric), apply(as.matrix(ctrs_W_spp), 2, as.numeric), miles=F)
  
  sampling_fct <- function(x, d){
    ecoreg_subset <- which(ctrs_W_orsim_dists[x,] < (d/2))
    species_winter_ecoregions[sample(ecoreg_subset, length(which(ecoregions_abund_presence_W[,x] > 0)), replace = F)]
  }
  
  pop_resampling <- list()
  for(j in 1:ncol(ecoregions_abund_presence_W)){
    ctrs_W_pop <- matrix(geom(terra::centroids(ecoregions[which(ecoregions_abund_presence_W[,j] > 0),]))[,3:4], ncol=2)
    if(nrow(ctrs_W_pop) > 1){
      max_dist <- max(rdist.earth(apply(as.matrix(ctrs_W_pop), 2, as.numeric), miles=F))
      pop_resampling[[j]] <- t(sapply(rep(j,1000), function(x) sampling_fct(x,max_dist)))
    }else if(nrow(ctrs_W_pop) == 1){
      pop_resampling[[j]] <- matrix(sample(1:length(species_winter_ecoregions), 1000), ncol=1)
    }
  }
  
  # Extract climate, and calculate niches and migration distances
  niche.overlap.null <- niche.overlap.temp.null <- niche.overlap.prec.null <- wintering.niche.size.null <- migra_dists_null <- list()
  for(j in 1:1000){
    niche_W_null <- niche_temp_W_null <- niche_prec_W_null <- list()
    for(i in 1:length(pop_resampling)){
      
      ecoregions_climate <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[pop_resampling[[i]][j,],])
      colnames(ecoregions_climate) <- c("ID", "temp", "prec")
      ecoregions_weights <- ecoregions_abund_presence[[k]]$wintering_abund[pop_resampling[[i]][j,]] + 0.1
      to_sample <- round((ecoregions_weights/sum(ecoregions_weights)) * 10000)
      ecoregions_W_climate_resample <- vector()
      for(h in 1:length(to_sample)){
        ecoregions_W_climate_resample <- rbind(ecoregions_W_climate_resample, ecoregions_climate[which(ecoregions_climate$ID == h),][sample(1:length(which(ecoregions_climate$ID == h)), to_sample[h], replace=T),])
      }
      ecoregions_W_climate_resample <- ecoregions_W_climate_resample %>% mutate(population = unique(pops_seas[,1])[i])
      niche_W_null[[i]] <- nicheDensityRaster(ecoregions_W_climate_resample %>% dplyr::select(temp, prec))
      niche_temp_W_null[[i]] <- stats::density(unlist(as.vector(ecoregions_W_climate_resample %>% dplyr::select(temp))), bw=0.25, kernel="gaussian", from=-2, to=3)
      niche_prec_W_null[[i]] <- stats::density(unlist(as.vector(ecoregions_W_climate_resample %>% dplyr::select(prec))), bw=0.25, kernel="gaussian", from=-2, to=3)
    }
    density_W_null <- rasterToPoints(raster::stack(niche_W_null))
    wintering.niche.size.null[[j]] <- apply(density_W_null[,-c(1,2)], 2, function(x){length(which(x>0))})
    niche.overlap.null[[j]] <- 1 - (0.5 * apply(abs(density.breeding[[k]][,-c(1,2)] - density_W_null[,-c(1,2)]), 2, sum))
    # Seasonal 1-D niche overlap
    niche.overlap.T <- niche.overlap.P <- vector()
    for(i in 1:length(pop_resampling)){
      # temperature
      X <- breeding.niche.temp[[k]][[i]]$x
      Y1 <- breeding.niche.temp[[k]][[i]]$y
      Y2 <- niche_temp_W_null[[i]]$y
      Overlap <- pmin ( Y1, Y2 )
      Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
      niche.overlap.T[i]  <- trapz ( X, Overlap ) / Total
      # precipitation
      X <- breeding.niche.prec[[k]][[i]]$x
      Y1 <- breeding.niche.prec[[k]][[i]]$y
      Y2 <- niche_prec_W_null[[i]]$y
      Overlap <- pmin ( Y1, Y2 )
      Total <- trapz ( X, Y1 ) + trapz ( X, Y2 )
      niche.overlap.P[i] <- trapz ( X, Overlap ) / Total
    }
    niche.overlap.temp.null[[j]] <- niche.overlap.T
    niche.overlap.prec.null[[j]] <- niche.overlap.P
    # migration distances
    ctrs_W <- vector()
    for(i in 1:length(pop_resampling)){
      ctrs_W <- rbind(ctrs_W, apply(matrix(geom(terra::centroids(ecoregions[pop_resampling[[i]][j,],]))[,3:4], ncol=2), 2, mean))
    }
    ctrs_B <- pop_centroids %>% filter(species == species.names[k]) %>% filter(season == "breeding") %>% dplyr::select(longitude, latitude)
    migra_dists_null[[j]] <- diag(rdist.earth(matrix(as.numeric(as.matrix(ctrs_B)), ncol=2, byrow=F), ctrs_W, miles=F))
  }
  migration_distances_null_4[[k]] <- do.call(rbind, migra_dists_null)
  niche_overlap_null_4[[k]] <- do.call(rbind, niche.overlap.null)
  niche_overlap_temp_null_4[[k]] <- do.call(rbind, niche.overlap.temp.null)
  niche_overlap_prec_null_4[[k]] <- do.call(rbind, niche.overlap.prec.null)
  wintering_niche_size_null_4[[k]] <- do.call(rbind, wintering.niche.size.null)
  print(k)
}
#load("results/output/null_model_4_results_4.RData")

# Comparing null and observed, treating each population as independent
pvals <- pvals_temp <- pvals_prec <- vector()
for(k in 1:length(species.names)){
  niche_data_spp <- niche_data %>% filter(species == species.names[k])
  for(i in 1:nrow(niche_data_spp)){
    pvals <- c(pvals, length(which(niche_overlap_null_4[[k]][,i] > niche_data_spp$niche_overlap[i])) / nrow(niche_overlap_null_4[[k]]))
    pvals_temp <- c(pvals_temp, length(which(niche_overlap_temp_null_4[[k]][,i] > niche_data_spp$niche_overlap_temp[i])) / nrow(niche_overlap_temp_null_4[[k]]))
    pvals_prec <- c(pvals_prec, length(which(niche_overlap_prec_null_4[[k]][,i] > niche_data_spp$niche_overlap_prec[i])) / nrow(niche_overlap_prec_null_4[[k]]))
  }
}
niche_data3 <- cbind(niche_data, pvals, pvals_temp, pvals_prec) %>% 
  unite("species_pop", species:population, remove=F) %>%
  filter(species_pop %!in% c("YEWA_E", "COYE_CA", "AMRE_C", "PABU_Lsn"))

##  Are distributions of scaled ranks skewed towards low values? ##
ks.test(niche_data3$pvals, "punif",0,1, alternative="greater") # yes
ks.test(niche_data3$pvals_temp, "punif",0,1, alternative="greater") # no
ks.test(niche_data3$pvals_prec, "punif",0,1, alternative="greater") # yes

# Without YEWA
niche_data_3b <- niche_data_3 %>% filter(species != "YEWA")
ks.test(niche_data_3b$pvals_prec, "punif",0,1, alternative="greater") # yes


##  Figure 2: Relationship between empirical climate niche overlap and null expectations   ##

niche_data_all <- niche_data_null_orsim_2 %>%
  left_join(niche_data3)

g_orsim_1 <- ggplot(niche_data_all, aes(x = migration_distance, y = migration_distance_orsim)) +
  geom_point(col = "black", size=2.5) + theme_classic() +
  geom_abline(aes(slope=1, intercept=0)) + ggtitle("(a)") + xlim(c(0,6500)) + ylim(c(0,6500)) +
  xlab("Empirical migration distance (km)") + ylab("Simulated migration distance (km)")+
  scale_size(name="1-rank") + scale_colour_discrete(name="")
g_orsim_2 <- ggplot(niche_data_all, aes(x = niche_overlap, y = niche_overlap_orsim)) +
  geom_point(aes(col = species, size = 1-pvals)) + theme_classic() +
  geom_abline(aes(slope=1, intercept=0)) + ggtitle("(b)") + xlim(c(0,0.9)) + ylim(c(0,0.9)) +
  xlab("Empirical climate overlap") + ylab("Simulated climate overlap")+
  scale_size(name="1-rank") + scale_colour_discrete(name="")
g_orsim_3 <- ggplot(niche_data_all, aes(x = niche_overlap_temp, y = niche_overlap_temp_orsim)) +
  geom_point(aes(col = species, size = 1-pvals_temp)) + theme_classic() +
  geom_abline(aes(slope=1, intercept=0)) + ggtitle("(c)") + xlim(c(0.05,0.5)) + ylim(c(0.05,0.5)) +
  xlab("Empirical thermal overlap") + ylab("Simulated thermal overlap") +
  scale_size(name="1-rank") + scale_colour_discrete(name="")
g_orsim_4 <- ggplot(niche_data_all, aes(x = niche_overlap_prec, y = niche_overlap_prec_orsim)) +
  geom_point(aes(col = species, size = 1-pvals_prec)) + theme_classic() +
  geom_abline(aes(slope=1, intercept=0)) + ggtitle("(d)") + xlim(c(0.05,0.5)) + ylim(c(0.05,0.5)) +
  xlab("Empirical precipitation overlap") + ylab("Simulated precipitation overlap") +
  scale_size(name="1-rank") + scale_colour_discrete(name="")

pdf("results/figures/Fig_2_full_new_6.pdf", width = 8, height = 8)
ggarrange(g_orsim_1, g_orsim_2, g_orsim_3, g_orsim_4, nrow=2, ncol=2, common.legend = TRUE, legend="bottom")
dev.off()




##  Fig. S1 — S6: distribution of samples  ##

##  WIFL  ##

data_for_analysis_spp <- data_for_analysis %>% filter(species == "WIFL")
species_climate_spp <- species_climate %>% filter(species == species.names[1])
pop.names <- unique(species_climate_spp$population)
# Pop 1: Pacific North West — change to West
# Pop 2: Interior North West — change to "South Rockies"
# Pop 3: East — change to "Central-East"
# Pop 4: South West — change to "South West"
k=1
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_sf(data = newmap, col="grey90", fill="grey90") + ylim(c(-5,60)) + xlim(c(-135,-65)) +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_wintering[[1]][[k]][,1],]), aes(fill=abs(-(ecoregions_wintering[[1]][[k]][,2]/sum(ecoregions_wintering[[1]][[k]][,2]))), col=abs(-(ecoregions_wintering[[1]][[k]][,2]/sum(ecoregions_wintering[[1]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_breeding[[1]][[k]][,1],]), aes(fill=abs(-(ecoregions_breeding[[1]][[k]][,2]/sum(ecoregions_breeding[[1]][[k]][,2]))), col=abs(-(ecoregions_breeding[[1]][[k]][,2]/sum(ecoregions_breeding[[1]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf("results/figures/Fig_S1_1.pdf", width = 10, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()


##  YEWA  ##

data_for_analysis_spp <- data_for_analysis %>% filter(species == "YEWA")
species_climate_spp <- species_climate %>% filter(species == species.names[2])
pop.names <- unique(species_climate_spp$population)

# Pop 1: Alaska — change to NW
# Pop 2: SouthWest — mistake! Change to NWC
# Pop 3: PNW — Change to "South Rockies"
# Pop 4: East — Keep "East"
# Pop 5: MidWest — Change to "Central"
k=1
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_sf(data = newmap, col="grey90", fill="grey90") + ylim(c(5,68)) + xlim(c(-160,-57)) +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_wintering[[2]][[k]][,1],]), aes(fill=abs(-(ecoregions_wintering[[2]][[k]][,2]/sum(ecoregions_wintering[[2]][[k]][,2]))), col=abs(-(ecoregions_wintering[[2]][[k]][,2]/sum(ecoregions_wintering[[2]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_breeding[[2]][[k]][,1],]), aes(fill=abs(-(ecoregions_breeding[[2]][[k]][,2]/sum(ecoregions_breeding[[2]][[k]][,2]))), col=abs(-(ecoregions_breeding[[2]][[k]][,2]/sum(ecoregions_breeding[[2]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf("results/figures/Fig_S2_1.pdf", width = 10, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()


##  WIWA  ##

data_for_analysis_spp <- data_for_analysis %>% filter(species == "WIWA")
species_climate_spp <- species_climate %>% filter(species == species.names[3])
pop.names <- unique(species_climate_spp$population)

# Pop 1: AK2Alberta — change to North West
# Pop 2: RockyMtn - change to Rockies
# Pop 3: PNW — Change to West
# Pop 4: CoastalCA — Change to South West 1
# Pop 5: Sierra — Change to South West 2
# Pop 6: Eastern — Change to North East
k=6
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_sf(data = newmap, col="grey90", fill="grey90") + ylim(c(6,67)) + xlim(c(-159,-60)) +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_wintering[[3]][[k]][,1],]), aes(fill=abs(-(ecoregions_wintering[[3]][[k]][,2]/sum(ecoregions_wintering[[3]][[k]][,2]))), col=abs(-(ecoregions_wintering[[3]][[k]][,2]/sum(ecoregions_wintering[[3]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_breeding[[3]][[k]][,1],]), aes(fill=abs(-(ecoregions_breeding[[3]][[k]][,2]/sum(ecoregions_breeding[[3]][[k]][,2]))), col=abs(-(ecoregions_breeding[[3]][[k]][,2]/sum(ecoregions_breeding[[3]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf("results/figures/Fig_S3_6.pdf", width = 11, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()


##  COYE  ##

data_for_analysis_spp <- data_for_analysis %>% filter(species == "COYE")
species_climate_spp <- species_climate %>% filter(species == species.names[4])
pop.names <- unique(species_climate_spp$population)

# Pop 1: West — change to Rockies—Central
# Pop 2: CA - change to California
# Pop 3: Southwest — Change to South West
# Pop 4: Midwest — Change to Central—East
# Pop 5: NewEngland — Change to East
k=1
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_sf(data = newmap, col="grey90", fill="grey90") + ylim(c(8,60)) + xlim(c(-135,-58)) +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_wintering[[4]][[k]][,1],]), aes(fill=abs(-(ecoregions_wintering[[4]][[k]][,2]/sum(ecoregions_wintering[[4]][[k]][,2]))), col=abs(-(ecoregions_wintering[[4]][[k]][,2]/sum(ecoregions_wintering[[4]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_breeding[[4]][[k]][,1],]), aes(fill=abs(-(ecoregions_breeding[[4]][[k]][,2]/sum(ecoregions_breeding[[4]][[k]][,2]))), col=abs(-(ecoregions_breeding[[4]][[k]][,2]/sum(ecoregions_breeding[[4]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf("results/figures/Fig_S4_1.pdf", width = 11, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()


##  AMRE  ##

data_for_analysis_spp <- data_for_analysis %>% filter(species == "AMRE")
species_climate_spp <- species_climate %>% filter(species == species.names[5])
pop.names <- unique(species_climate_spp$population)

# Pop 1: Northwest — change to North West-Rockies
# Pop 2: South - change to South Central
# Pop 3: Northeast — Change to North Central-East
# Pop 4: SouthDakota — Change to Central
k=4
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_sf(data = newmap, col="grey90", fill="grey90") + ylim(c(5,62)) + xlim(c(-135,-58)) +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_wintering[[5]][[k]][,1],]), aes(fill=abs(-(ecoregions_wintering[[5]][[k]][,2]/sum(ecoregions_wintering[[5]][[k]][,2]))), col=abs(-(ecoregions_wintering[[5]][[k]][,2]/sum(ecoregions_wintering[[5]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_breeding[[5]][[k]][,1],]), aes(fill=abs(-(ecoregions_breeding[[5]][[k]][,2]/sum(ecoregions_breeding[[5]][[k]][,2]))), col=abs(-(ecoregions_breeding[[5]][[k]][,2]/sum(ecoregions_breeding[[5]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf("results/figures/Fig_S5_4.pdf", width = 11, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()


##  PABU  ##

data_for_analysis_spp <- data_for_analysis %>% filter(species == "PABU")
species_climate_spp <- species_climate %>% filter(species == species.names[6])
pop.names <- unique(species_climate_spp$population)

# Pop 1: Central — change to South Central
# Pop 2: Louisiana - change to Louisiana
k=2
theme_set(theme_classic())
g_map_pop <- ggplot() + 
  geom_sf(data = newmap, col="grey90", fill="grey90") + ylim(c(5,55)) + xlim(c(-130,-65)) +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_wintering[[6]][[k]][,1],]), aes(fill=abs(-(ecoregions_wintering[[6]][[k]][,2]/sum(ecoregions_wintering[[6]][[k]][,2]))), col=abs(-(ecoregions_wintering[[6]][[k]][,2]/sum(ecoregions_wintering[[6]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("azure1","skyblue3")) + scale_fill_gradientn(colours = c("azure1","skyblue3")) +
  ggnewscale::new_scale_colour() + ggnewscale::new_scale_fill() +
  geom_sf(data = st_as_sf(ecoregions[ecoregions_breeding[[6]][[k]][,1],]), aes(fill=abs(-(ecoregions_breeding[[6]][[k]][,2]/sum(ecoregions_breeding[[6]][[k]][,2]))), col=abs(-(ecoregions_breeding[[6]][[k]][,2]/sum(ecoregions_breeding[[6]][[k]][,2]))))) +
  scale_color_gradientn(colours = c("cornsilk","coral2")) + scale_fill_gradientn(colours = c("cornsilk","coral2")) +
  geom_point(data = data_for_analysis_spp %>% filter(population == pop.names[k]), aes(x=longitude, y=latitude), col="black", size=0.3) + theme_void() + 
  theme(legend.title=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
g_temp_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = temp, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Temperature") + ylab("Density") + theme(legend.position="none")
g_prec_pop <- ggplot(data = species_climate_spp %>% filter(population == pop.names[k])) + 
  geom_density(aes(x = prec, fill = season), bw=0.25, kernel="gaussian", alpha=0.4) + xlim(c(-2, 4)) +
  xlab("Precipitation") + ylab("Density") + theme(legend.position="none")

pdf("results/figures/Fig_S6_2.pdf", width = 11, height = 3)
grid.arrange(g_map_pop, g_temp_pop, g_prec_pop, ncol=3)
dev.off()



