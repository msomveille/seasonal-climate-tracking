###  Preparing data for analysis  ###

library(tidyverse)
library(dggridR)
library(rworldmap)
library(rgeos)
library(terra)
library(sf)
library(adehabitatHR)
library(ebirdst)
library(rnaturalearth)
library(raster)
library(sfsmisc)

###
###   Location and population assignment of sampled individuals   ###
###

# Species names
species.names <- c("WIFL", "YEWA", "WIWA", "COYE", "AMRE", "PABU")

##  Breeding season
# Get breeding assignments 
breeding.assignments.files <- list(
  WIFL = read.table("resources/WIFL/Appendix_1_Indiv_Assignments.csv", header=T, sep=";"),
  YEWA = readRDS("resources/YEWA/YWAR_cleaned.rds"),
  WIWA = read.table("resources/WIWA/NEW/WIWA.Breeding.location_popassignment_structure.txt", header=T),
  COYE = read.table("resources/COYE/COYE_breeding.assignment.csv", header=T, sep=";"),
  AMRE = read.csv("resources/AMRE/amre_breeding_assignments.csv"),
  PABU = read.csv("resources/PABU/NEW/PABU_breeding_assignments.csv")
)
breeding.assignments.files$WIFL <- breeding.assignments.files$WIFL[which(breeding.assignments.files$WIFL$Stage == "Breeding"),]
breeding.assignments.files$YEWA <- breeding.assignments.files$YEWA[which(breeding.assignments.files$YEWA$Stage == "Breeding"),]
# Which columns indicate assignment probabilities for populations
breeding.assignments.pops <- list(
  WIFL = 31:37,
  YEWA = NA,
  WIWA = 3:8,
  COYE = 7:11,
  AMRE = 2:6,
  PABU = 2:5
)
# Put all the data in a data frame
data_for_analysis <- data.frame()
for(k in 1:length(species.names)){
  # Assignment of breeding individuals to breeding populations
  breeding.assignments <- breeding.assignments.files[[k]]
  if(k==2){
    popBR <- as.character(breeding.assignments$Pop)
  }else{
    popBR <- apply(breeding.assignments[,breeding.assignments.pops[[k]]], 1, function(x) colnames(breeding.assignments[,breeding.assignments.pops[[k]]])[which(x>=0.7)])
    popBR[which(lapply(popBR, length) == 0)] <- NA
    popBR <- unlist(popBR)
    breeding.assignments <- breeding.assignments[which(is.na(popBR)==F),]
    popBR <- popBR[which(is.na(popBR)==F)]
  }
  # Location of breeding individuals
  if(k==3){
    ptsBR <- as.matrix(cbind(breeding.assignments$Longitude, breeding.assignments$Latitude))
  }else{
    ptsBR <- as.matrix(cbind(breeding.assignments$Long, breeding.assignments$Lat))
  }
  # Store data in a data frame (with each row being a sampled individual)
  output <- data.frame(
    species = rep(species.names[k], nrow(ptsBR)),
    season = rep("Breeding", nrow(ptsBR)),
    population = popBR,
    longitude = ptsBR[,1],
    latitude = ptsBR[,2]
  )
  data_for_analysis <- rbind(data_for_analysis, output)
}

##  Wintering season
# Get wintering assignments 
wintering.assignments.files <- list(
  WIFL = readRDS("resources/WIFL/winter-migrant-assignments.rds"),
  YEWA = read_csv("resources/YEWA/NEW/yewa.winter.df.csv"),
  WIWA = read.table("resources/WIWA/NEW/WIWA.wint_birds.pop_assignment.txt", header=T, sep="\t"),
  COYE = read.table("resources/COYE/NEW/Coye.Wint_only.inclPlate29.inclTaylor.rm_miss10.rep_indiv_estimates.meta.Appendix2.redo_noCAs.txt", header=T, sep="\t"),
  AMRE = read.csv("resources/AMRE/wintering.assignment.consistent.135.csv"),
  PABU = read.table("resources/PABU/NEW/PABU.WintMigrVag.rm_miss20.rep_indiv_est.4map.lat_long.txt", header=T, sep="\t")
)
wintering.assignments.files$WIFL <- wintering.assignments.files$WIFL[which(wintering.assignments.files$WIFL$Stage == "Wintering"),]
wintering.assignments.files$COYE <- wintering.assignments.files$COYE[which(wintering.assignments.files$COYE$Notes == "Wintering"),]
wintering.assignments.files$PABU <- wintering.assignments.files$PABU[which(wintering.assignments.files$PABU$Stage == "WINTER"),]
# Which columns indicate assignment probabilities for populations
wintering.assignments.pops <- list(
  WIFL = NA,
  YEWA = NA,
  WIWA = 5:10,
  COYE = 5:9,
  AMRE = NA,
  PABU = NA
)
# Add all the data in the data frame with breeding individuals
for(k in 1:length(species.names)){
  # Assignment of breeding individuals to breeding populations
  wintering.assignments <- wintering.assignments.files[[k]]
  if(k %in% c(1,6)){
    popNB <- wintering.assignments$repunit
    popNB[which(wintering.assignments$PofZ < 0.7)] <- NA
    wintering.assignments <- wintering.assignments[which(is.na(popNB)==F),]
    popNB <- popNB[which(is.na(popNB)==F)]
  }
  if(k==2){
    popNB <- as.character(wintering.assignments$Pop)
  }
  if(k==5){
    popNB <- as.character(wintering.assignments$AssignedPop_Full)
  }
  if(k %in% 3:4){
    popNB <- apply(wintering.assignments[,wintering.assignments.pops[[k]]], 1, function(x) colnames(wintering.assignments[,wintering.assignments.pops[[k]]])[which(x>=0.7)])
    popNB[which(lapply(popNB, length) == 0)] <- NA
    popNB <- unlist(popNB)
    wintering.assignments <- wintering.assignments[which(is.na(popNB)==F),]
    popNB <- popNB[which(is.na(popNB)==F)]
  }
  # Location of breeding individuals
  if(k==1){
    ptsNB <- as.matrix(cbind(wintering.assignments$lon, wintering.assignments$lat))
  }
  if(k==3){
    ptsNB <- as.matrix(cbind(wintering.assignments$x, wintering.assignments$y))
  }
  if(k==5){
    ptsNB <- as.matrix(cbind(wintering.assignments$Lon, wintering.assignments$Lat))
  }
  if(k %in% c(2,4,6)){
    ptsNB <- as.matrix(cbind(wintering.assignments$Long, wintering.assignments$Lat))
  }
  # Store data in a data frame (with each row being a sampled individual)
  output <- data.frame(
    species = rep(species.names[k], nrow(ptsNB)),
    season = rep("Wintering", nrow(ptsNB)),
    population = popNB,
    longitude = ptsNB[,1],
    latitude = ptsNB[,2]
  )
  data_for_analysis <- rbind(data_for_analysis, output)
}

data_for_analysis$population[which(data_for_analysis$population=="PacNorthwest")] <- "PNW"
data_for_analysis$population[which(data_for_analysis$population=="1")] <- "SouthWest"
data_for_analysis$population[which(data_for_analysis$population=="2")] <- "MidWest"
data_for_analysis$population[which(data_for_analysis$population=="3")] <- "PNW"
data_for_analysis$population[which(data_for_analysis$population=="4")] <- "East"
data_for_analysis$population[which(data_for_analysis$population=="5")] <- "Alaska"
data_for_analysis$population[which(data_for_analysis$population=="Cluster1")] <- "Louisiana"
data_for_analysis$population[which(data_for_analysis$population=="Cluster2")] <- "South"
data_for_analysis$population[which(data_for_analysis$population=="Cluster3")] <- "East"
data_for_analysis$population[which(data_for_analysis$population=="Cluster4")] <- "Central"

write.csv(data_for_analysis, "results/output/data_for_analysis.csv", row.names = F)
