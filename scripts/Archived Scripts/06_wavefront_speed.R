###############################################################################
### This function is used to approximate the leading edge of an infection wave 
### and key properties, such as length and estimated speed of expansion. 
### It will identify hexes above a defined threshold (infThreshold) at time 
### (t) for a set of infections per capita estimates distributed
### across a hexgrid (hexObs)
###############################################################################

rm(list = ls())
gc()

library(tidyverse)
library(sf)
library(units)

hexes <- sf::st_read("Data/data-sources/hexes.geojson") |> 
  filter(as.integer(hexid) < 7662)

hexes_centroids <- as.data.frame(st_coordinates(st_cast(st_centroid(hexes), "MULTIPOINT"))) |> 
  rename(hexid = L1) |> 
  mutate(hexid = as.character(hexid)) |> 
  left_join(hexes) |> 
  mutate(geometry = st_centroid(geometry)) |> 
  st_as_sf()

## US States border lines
us_states <- tigris::states(cb = T, resolution = "500k") |> 
  st_transform(crs = 26915) |> 
  filter(!NAME %in% c("Alaska",
                      "Hawaii",
                      "Puerto Rico", 
                      "Guam", 
                      "Commonwealth of the Northern Mariana Islands",
                      "American Samoa",
                      "United States Virgin Islands"))

## hexgrids
dataset <- "preomicron"

## Pre-Omicron
CAR_df_preomicron <- vroom::vroom(paste0("Data/data-sources/tsa_",
                                         dataset, 
                                         ".csv"))

## peak date
alpha_peak <- as.Date("2020-11-19")
delta_peak <- as.Date("2021-09-04")

## Threshold for filtering given the distribution, any value of trend that it is above the 3rd Quarter of the trend distribution
threshold_mean <- 165

## Pre-Omicron

## Alpha
CAR_lag_alpha <- CAR_df_preomicron |> 
  mutate(hexid = as.character(hexid)) |> 
  left_join(hexes) |> 
  dplyr::select(hexid, date, mean, sd, geometry) |>
  st_as_sf() |> 
  filter(date %in% seq.Date(to = alpha_peak,
                            from = (alpha_peak-63),
                            length.out = 63))|>
  st_transform(crs = 26915)

## delta
CAR_lag_delta <- CAR_df_preomicron |> 
  mutate(hexid = as.character(hexid)) |> 
  left_join(hexes) |> 
  dplyr::select(hexid, date, mean, sd, geometry) |>
  st_as_sf() |> 
  filter(date %in% seq.Date(to = delta_peak,
                            from = (delta_peak-63),
                            length.out = 63))|>
  st_transform(crs = 26915)

### testing parameters
infThreshold <- 165 ##daily infections per capita
t <- as.Date("2021-08-18")
hexObs <- CAR_lag_delta
hexgrid <- sf::st_read("Data/data-sources/hexes.geojson")
###

idWavefront <- function(hexObs,
                        hexgrid,
                        infThreshold,
                        t){
  ###############################################################################
  ### Load packages
  ###############################################################################  
  library(dplyr)
  ###############################################################################
  ### Filter the dataset
  ###############################################################################
  filtHexObs <- hexObs %>% dplyr::filter(date == t) %>%
    mutate(wavefront = ifelse(mean > infThreshold, 
                              TRUE, FALSE),
           ### Find the centroid of each hex  
           centroid = st_centroid(geometry))
  
  ### Check that the data set is nonempty: 
  if (sum(filtHexObs$wavefront) == 0) {
    print(paste("No hexes identified above threshold of", infThreshold, 
                "infections per day. Maximum infections estimated on", t,
                "is", 
                round(max((hexObs %>% dplyr::filter(date == t))[,"infectionsPC"])*1e5)))
    break()
  } else{
    print(paste("Identified", sum(filtHexObs$wavefront), 
                "hexes with daily infections per capita greater than", 
                infThreshold,".", 
                "This represents", round(sum(filtHexObs$wavefront)/length(unique(filtHexObs$hexid)),2)*100, 
                "% of all hexes."))
  }
  
  ### Check with a plot
  # ggplot() + 
  #   geom_sf(data=filtHexObs, mapping=aes(fill=wavefront)) +
  #   geom_sf(data=filtHexObs %>% filter(wavefront==TRUE),
  #           mapping=aes(geometry= centroid), size =.01)
  
  wave <- filtHexObs %>% filter(wavefront == TRUE) %>%
    mutate(waveNumb = ifelse(as.numeric(hexid) > 2400, 1, 1)) %>% 
    filter(waveNumb == 1)
  
  notWave <- filtHexObs %>% filter(wavefront == FALSE) %>%
    mutate(waveNumb = ifelse(as.numeric(hexid) > 2400, 1, 1)) %>% 
    filter(waveNumb == 1)
  
  boundary <- st_intersection(wave, notWave) %>% st_simplify()
  
  # plotBound <- ggplot() + 
  #   geom_sf(data=wave, mapping=aes(geometry= geometry), fill="salmon") +
  #   geom_sf(data=notWave, mapping=aes(geometry= geometry), fill="orange") +
  #   geom_sf(data=boundary, mapping=aes(geometry= geometry), color="dodgerblue") + 
  #   ggtitle(paste(t))
  
  
  ### Calculate the length of the boundary 
  boundLength <- sum(st_length(boundary))
  
  ### Create a list of the relevant outputs 
  waveFrontList <- list("wave" = wave, 
                        "not wave" = notWave, 
                        "boundary" = boundary, 
                        "boundary length" = boundLength,
                        "boundary plot" = plotBound)
  
  return(waveFrontList)
  
}

wavesurface_list <- CAR_lag_delta |> 
  group_split(date)

dates_wave <- unique(CAR_lag_delta$date)

waveboundary_list <- lapply(wavesurface_list, idWavefront, hexes, 165, unique(CAR_lag_delta$date))
