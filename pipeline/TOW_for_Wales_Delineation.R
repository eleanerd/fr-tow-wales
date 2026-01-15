##################
# Freddie Hunter
# 03 Dec 2025
# Wales TOW
#################

# Review by Eleanor Downer on 15 Jan 2026

# Assumptions
# This code assumes that the input is a VOM with a minimum pixel cluster size of 5pixel.
# The code also assumes that the VOM has been filtered with a minimum of 130cm 
# The code also assumes that the VOM has NA values not zero values

rm(list = ls())

#############################################################

# Input files
# There are many parameters imbedded in the code. 
working_dir <- "//forestresearch.gov.uk/shares/IFOS/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/"
setwd(working_dir)

chm_raster_path <- '04_Spatial Analysis/4_Processing/CHM_FILT/SS79_CHM_FILT.tif'
nfi_layer_path  <- "1_Reference_Data/6_Wales_NFI_2024/NFI_Wales_2024_WoodlandOnly.gpkg"

#############################################################
# Libraries
library(terra)
library(dplyr)
library(sf)
library(raster)
library(glue)
library(deldir)
library(lidR)
library(lwgeom)
library(nngeo)
library(smoothr)
library(purrr)
library(foreach)
library(doParallel)

#################################################
# Processing

# Load CHM raster
chm_raster <- terra::rast(chm_raster_path)

# Split into tiles
chm_tiles <- terra::makeTiles(chm_raster, c(2000, 2000))  
tile_file_list <- list.files(
  "04_Spatial Analysis/4_Processing/CHM_FILT/",
  pattern = "^tile_[0-9]{1,3}\\.tif$",
  full.names = TRUE
)

###############################
### Run first part in parallel
###############################

num_cores <- parallel::detectCores()

# Create a cluster
cluster <- makeCluster(num_cores - 15)
registerDoParallel(cluster)

# Run method per tile
tile_results <- foreach(
  tile_index = 1:length(tile_file_list),
  .packages = c(
    "sf",
    "dplyr",
    "units",
    "nngeo",
    "terra",
    "raster",
    "glue",
    "deldir",
    "lidR",
    "lwgeom",
    "smoothr",
    "purrr"
  ),
  .export = c("tile_file_list", "nfi_layer_path", "working_dir")
) %dopar% {
  
  tile_path <- tile_file_list[tile_index]
  
  if (file.exists(glue("{working_dir}/filt_max_class_{tile_path}.tif"))) {
    print('processed')
    next()
  }
  
  tile_chm <- rast(tile_path)
  
  if (is.nan(minmax(tile_chm)[2])) {
    print('no data exists')
    return(NULL)
  } else { 
    print('data exists')
  }
  
  tile_name <- basename(tile_path) %>% stringr::str_remove('.tif')
  
  # Step 1: 5x5 local maximum filter
  local_max_raster <- terra::focal(
    tile_chm,
    w = matrix(1, 5, 5),
    fun = max,
    na.rm = TRUE
  )
  
  # Restore NA mask
  local_max_raster[is.na(tile_chm)] <- NA
  
  # Step 2: Initial height-based classification
  reclass_matrix <- matrix(
    c(0.49, 10, 3,
      10, 15, 10,
      15, 20, 15,
      20, 30, 20,
      30, Inf, 30),
    ncol = 3,
    byrow = TRUE
  )
  
  local_max_raster <- as(local_max_raster, "SpatRaster")
  height_class_raster <- classify(local_max_raster, reclass_matrix)
  
  # First sieve (5 pixels)
  height_class_raster[is.na(height_class_raster)] <- 1
  sieved_class_small <- sieve(
    height_class_raster,
    threshold = 5,
    directions = 4
  )
  
  # Second sieve (25 pixels)
  sieved_class_small[sieved_class_small == 1] <- NA
  filtered_height_class <- sieve(
    sieved_class_small,
    threshold = 25,
    directions = 4
  )
  
  # Tree detection model
  window_heights <- c(3, 5, 10, 20)
  window_sizes   <- c(2, 2.6, 3.2, 5)
  
  size_model <- lm(log(window_sizes) ~ log(window_heights))
  model_a <- exp(coef(size_model)[1])
  model_b <- coef(size_model)[2]
  
  adaptive_window <- function(x) model_a * x^model_b
  
  # CHM smoothing
  smoothing_window <- matrix(1, nrow = 3, ncol = 3)
  
  tile_chm[is.na(tile_chm)] <- 0
  smoothed_chm <- terra::focal(
    x = tile_chm,
    w = smoothing_window,
    fun = mean,
    na.rm = TRUE
  )
  
  smoothed_chm[is.na(filtered_height_class)] <- NA
  smoothed_chm[smoothed_chm < 0.5] <- NA
  
  treetops <- locate_trees(
    smoothed_chm,
    lmf(ws = adaptive_window, shape = 'circular', hmin = 3)
  )
  
  # Transform treetops CRS
  nfi_layers <- st_layers(nfi_layer_path)
  nfi_meta   <- st_read(
    nfi_layer_path,
    query = glue("SELECT * from {nfi_layers$name} LIMIT 1 OFFSET 1")
  )
  
  treetops <- st_transform(treetops, crs = st_crs(nfi_meta))
  st_write(
    treetops,
    glue('04_Spatial Analysis/4_Processing/TOW_POLYGONS/ttops_{tile_name}.gpkg'),
    append = FALSE
  )
  
  # Crown delineation
  crown_raster <- dalponte2016(
    chm = smoothed_chm,
    treetops = treetops
  )()
  
  crown_raster[smoothed_chm >= 3] <- 2
  crown_raster[!is.na(crown_raster)] <- 2
  crown_raster[is.na(crown_raster)] <- 1
  
  crown_raster_filled <- sieve(
    as(crown_raster, 'SpatRaster'),
    threshold = 5,
    directions = 4
  )
  
  filtered_height_class[crown_raster_filled == 1] <- NA
  height_class_bins <- filtered_height_class
  
  writeRaster(
    height_class_bins,
    glue("04_Spatial Analysis/4_Processing/TOW_POLYGONS/filt_max_class_bins_{tile_name}.tif"),
    overwrite = TRUE
  )
  
  filtered_height_class[!is.na(filtered_height_class)] <- 1
  
  # Memory cleanup
  objects_to_remove <- ls()
  objects_to_remove <- objects_to_remove[!objects_to_remove %in% c(
    "filtered_height_class",
    "tile_name",
    "working_dir", 
    "nfi_layer_path", 
    "tile_file_list"
  )]  
  rm(list = objects_to_remove)
  
  # Polygonize canopy
  canopy_polygons <- st_as_sf(
    stars::st_as_stars(filtered_height_class),
    point = FALSE,
    merge = TRUE,
    connect8 = TRUE
  ) %>%
    st_union() %>%
    st_cast("POLYGON") %>%
    st_transform(st_crs(nfi_meta))
  
  # ---- NFI OHC SECTION ----
  
  canopy_bbox <- st_bbox(canopy_polygons, crs = st_crs(canopy_polygons)) %>%
    st_as_sfc() %>%
    st_transform(st_crs(nfi_meta))
  
  bbox_wkt <- st_as_text(canopy_bbox)
  
  nfi_crop <- st_read(
    nfi_layer_path,
    wkt_filter = bbox_wkt
  ) %>%
    st_collection_extract("POLYGON") %>%
    st_union() %>%
    st_cast('POLYGON') %>%
    st_as_sf()
  
  if (nrow(nfi_crop) != 0) {
    nfi_crop$ROWID <- seq_len(nrow(nfi_crop))
  } else {
    nfi_crop$ROWID <- NULL
  }
  
  canopy_minus_nfi <- st_difference(canopy_polygons, nfi_crop)
  nfi_buffer_10m  <- nfi_crop %>% st_buffer(10)
  
  nfi_ohc <- st_intersection(canopy_minus_nfi, nfi_buffer_10m) %>%
    st_collection_extract("POLYGON") %>%
    st_union() %>%
    st_cast('POLYGON')
  
  nfi_ohc_final <- st_difference(nfi_ohc, nfi_crop) %>%
    st_union() %>%
    st_cast('POLYGON') %>%
    st_as_sf()
  
  if (nrow(nfi_ohc_final) != 0) {
    nfi_ohc_final$id <- seq_len(nrow(nfi_ohc_final))
  } else {
    nfi_ohc_final$id <- NULL
  }
  
  overlaps <- st_overlaps(nfi_ohc_final, nfi_crop) %>% as.data.frame()
  touches  <- st_touches(nfi_ohc_final, nfi_crop) %>% as.data.frame()
  
  keep_ids <- rbind(overlaps, touches) %>% dplyr::select(row.id)
  
  if (nrow(keep_ids) != 0) {
    nfi_ohc_final <- nfi_ohc_final %>%
      filter(id %in% unlist(keep_ids))
  }
  
  nfi_ohc_final$area <- st_area(nfi_ohc_final) %>% units::drop_units()
  nfi_ohc_final <- nfi_ohc_final[nfi_ohc_final$area >= 5, ]
  
  if (nrow(keep_ids) != 0) {
    nfi_ohc_final$Woodland_Type <- 'NFI OHC'
  }
  
  # Remove OHC from TOW
  nfi_ohc_final <- st_transform(nfi_ohc_final, st_crs(canopy_polygons))
  canopy_polygons <- st_make_valid(canopy_polygons)
  nfi_ohc_final   <- st_make_valid(nfi_ohc_final)
  
  intersect_idx <- st_intersects(canopy_polygons, nfi_ohc_final)
  
  tow_geoms <- map2(
    st_geometry(canopy_polygons),
    intersect_idx,
    \(geom, ids) {
      geom <- st_sfc(geom, crs = st_crs(canopy_polygons))
      if (length(ids) == 0) return(geom)
      st_difference(geom, st_union(nfi_ohc_final[ids, ]))
    }
  )
  
  canopy_tow <- do.call(c, tow_geoms) %>%
    st_union() %>%
    st_cast('POLYGON') %>%
    st_as_sf()
  
  st_write(
    nfi_ohc_final,
    glue('04_Spatial Analysis/4_Processing/TOW_POLYGONS/NFI_OHC_F_{tile_name}.gpkg'),
    append = FALSE
  )
  
  st_write(
    canopy_tow,
    glue('04_Spatial Analysis/4_Processing/TOW_POLYGONS/canopy_area_tow_{tile_name}.gpkg'),
    append = FALSE
  )
}

stopCluster(cluster)

# End part 1
