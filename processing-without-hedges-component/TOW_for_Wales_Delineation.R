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
wd <- "//forestresearch.gov.uk/shares/IFOS/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/"
setwd(wd)

chm_raster_path <- "04_Spatial Analysis/4_Processing/CHM_FILT/SS79_CHM_FILT.tif"
nfi_layer_path  <- "04_Spatial Analysis/1_Reference_Data/6_Wales_NFI_2024/NFI_Wales_2024_WoodlandOnly.gpkg"

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

################################################
# Functions


#################################################
# Processing

# Load CHM raster
chm_raster <- terra::rast(chm_raster_path)
tile_name <- "SS79"
  
# Step 1: 5x5 local maximum filter
local_max_raster <- terra::focal(
  chm_raster,
  w = matrix(1, 5, 5),
  fun = max,
  na.rm = TRUE
)
  
# Restore NA mask
local_max_raster[is.na(chm_raster)] <- NA
  
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
  
chm_raster[is.na(chm_raster)] <- 0
smoothed_chm <- terra::focal(
  x = chm_raster,
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
treetops <- st_transform(treetops, crs = st_crs(27700))

# Load into memory
chm_r <- raster(smoothed_chm)
crs(chm_r) <- crs(27700)

# Crown delineation
crown_raster <- dalponte2016(chm = smoothed_chm, treetops = treetops)()

crown_raster[smoothed_chm >= 3] <- 2
crown_raster[!is.na(crown_raster)] <- 2
crown_raster[is.na(crown_raster)] <- 1
  
crown_raster_filled <- sieve(
  as(crown_raster, 'SpatRaster'),
  threshold = 5,
  directions = 4
)
  
filtered_height_class[crown_raster_filled == 1] <- NA
filtered_height_class_2 <- filtered_height_class
filtered_height_class[!is.na(filtered_height_class)] <- 1

# Polygonize canopy
canopy_polygons <- st_as_sf(
  stars::st_as_stars(filtered_height_class),
  point = FALSE,
  merge = TRUE,
  connect8 = TRUE
) %>%
  st_union() %>%
  st_cast("POLYGON") %>%
  st_as_sf() %>%
  st_transform(st_crs(27700))

# Can remove:
#rm(height_class_bins)
#rm(crown_raster)
#rm(crown_raster_filled)
#gc()

# End part 1

########################################################################
# Start Part 2

# ---- Load NFI Woodland  ----

# Convert bounding box to polygon, then get bbox as WKT 
dat_bbox <- st_as_sfc(st_bbox(canopy_polygons, crs = st_crs(27700)))
bbox_wkt <- st_as_text(dat_bbox)

nfi <- st_read(nfi_layer_path,
               wkt_filter = bbox_wkt) %>%
  filter(units::drop_units(st_area(.)) > 0) %>%
  st_union() %>%
  st_cast("MULTIPOLYGON") %>%
  st_as_sf()

# --- Separate TOW polygons by intersection with NFI ---
polygons_in_nfi <- st_filter(canopy_polygons, nfi, .predicate = st_intersects)
polygons_out_nfi <- st_filter(canopy_polygons, nfi, .predicate = st_disjoint)

# =========================================
# Classify NFI OHC
# =========================================

# Clip polygons_in_nfi by NFI (remove NFI from TOW)
tow_nfi_diff <- polygons_in_nfi %>%
  st_difference(nfi) %>%
  st_collection_extract("POLYGON") %>%
  st_cast("POLYGON") %>%
  st_as_sf %>%
  mutate(
    area = round(units::drop_units(st_area(.)), 2) # Calculate area
  ) %>%
  filter(area > 10) # Remove polygons under 10m2

# Buffer NFI by 10 m
nfi_buf10 <- st_buffer(nfi, 10)
nfi_buf_only <- st_difference(nfi_buf10, nfi) %>%
  st_collection_extract("POLYGON") %>%
  st_cast("POLYGON") %>%
  st_as_sf()

# NFI OHC (polygons intersecting buffer only)
NFI_OHC <- tow_nfi_diff %>%
  st_intersection(nfi_buf_only) %>%
  st_collection_extract("POLYGON") %>%
  st_cast("POLYGON") %>%
  st_make_valid() %>%
  st_as_sf() %>%
  mutate(
    area = round(units::drop_units(st_area(.)), 2), # Calculate area
    Woodland_Type = "NFI OHC"
  ) %>%
  filter(area > 5)

# Keep polygons that intersect or touch NFI Woodland
NFI_OHC <- NFI_OHC %>%
  st_filter(nfi, .predicate = function(a, b) {
    intersects_matrix <- st_intersects(a, b, sparse = FALSE)
    touches_matrix    <- st_touches(a, b, sparse = FALSE)
    
    # Keep polygons that either intersect OR touch
    apply(intersects_matrix | touches_matrix, 1, any)
  })


nfi_buf_only <- nfi_buf_only %>% st_union()

# Non-intersect polygons_out_nfi classification by are
connecting_polygons <- tow_nfi_diff %>%
  st_difference(nfi_buf_only) %>%
  st_collection_extract("POLYGON") %>%
  st_cast("POLYGON") %>%
  st_make_valid() %>%
  mutate(
    area = round(units::drop_units(st_area(.)), 2),
    Woodland_Type = case_when(
      area < 1000 ~ "Group of Trees",
      area >= 1000 ~ "Small Woodland"
    )
  ) %>%
  filter(area > 5)

NFI_OHC$trct <- lengths(st_intersects(NFI_OHC, treetops))
connecting_polygons$trct <- lengths(st_intersects(connecting_polygons, treetops))

# ---- Calculate canopy attributes ----

polygons_out_nfi$area <- st_area(polygons_out_nfi) %>% units::drop_units() %>% round(2)
polygons_out_nfi <- polygons_out_nfi[polygons_out_nfi$area >= 5, ] # 
polygons_out_nfi$mbc <- st_area(lwgeom::st_minimum_bounding_circle(polygons_out_nfi$x)) %>% units::drop_units()
polygons_out_nfi$trct <- lengths(st_intersects(polygons_out_nfi, treetops))
polygons_out_nfi$trct[is.na(polygons_out_nfi$trct)] <- 0
polygons_out_nfi$mbc_p <- polygons_out_nfi$area / polygons_out_nfi$mbc * 100

# ---- Assign Woodland Types ----

polygons_out_nfi <-polygons_out_nfi %>%
  mutate(Woodland_Type = case_when(
      area <= 350 & trct == 0 ~ 'Lone Tree',
      area <= 350 & mbc_p >= 55 ~ 'Lone Tree',
      trct == 1 ~ 'Lone Tree',
      area <= 350 & trct == 2 & mbc_p >= 50 ~ 'Lone Tree',
      area <= 1000 & trct > 1 ~ 'Group of Trees',
      area > 1000 ~ 'Small Woodland'
    ))

# ---- Merge datasets ----

# Ensure all datasets have the same columns
common_cols <- Reduce(
  intersect,
  list(names(NFI_OHC),
       names(connecting_polygons),
       names(polygons_out_nfi)))


polygons_joined <- dplyr::bind_rows(
  NFI_OHC %>% dplyr::select(all_of(common_cols)),
  connecting_polygons %>% dplyr::select(all_of(common_cols)),
  polygons_out_nfi %>% dplyr::select(all_of(common_cols))
)

# Add TOW ID
num_digits <- nchar(length(polygons_joined$Woodland_Type))
polygons_joined$TOW_ID <- paste0(
  formatC(1:length(polygons_joined$Woodland_Type), width = num_digits, flag = "0"),
  glue("_SS79")
)

# ---- Merge height-class rasters ----

height_bin_polygons <- st_as_sf(
  stars::st_as_stars(filtered_height_class_2),
  point = FALSE,
  connect8 = TRUE,
  merge = TRUE
) %>%
  st_cast("POLYGON") %>%
  st_transform(st_crs(27700))

# ---- Burn height bins into canopy polygons ----

canopy_with_bins <- st_intersection(
  polygons_joined,
  height_bin_polygons
) %>%
  st_collection_extract("POLYGON") %>%
  st_cast('POLYGON')

canopy_with_bins$Bin_Area <- st_area(canopy_with_bins) %>% units::drop_units()

# ---- Add CHM height statistics ----

height_stats <- terra::extract(
  chm_raster,
  canopy_with_bins,
  fun = function(x) {
    c(
      min  = round(min(x, na.rm = TRUE), 2),
      mean = round(mean(x, na.rm = TRUE), 2),
      max  = round(max(x, na.rm = TRUE), 2),
      sd   = round(sd(x, na.rm = TRUE), 2)
    )
  }
)

canopy_with_bins$min_ht  <- height_stats[, 2]
canopy_with_bins$mean_ht <- height_stats[, 3]
canopy_with_bins$max_ht  <- height_stats[, 4]
canopy_with_bins$sd_ht   <- height_stats[, 5]

# Assign sub_id attribute

canopy_with_bins <- canopy_with_bins %>%
  group_by(TOW_ID) %>%
  mutate(sub_id = row_number()) %>%
  ungroup()

# Re-order cols
canopy_with_bins <- canopy_with_bins[, c("TOW_ID", "Woodland_Type", "area", "spat_76f434aa2110_30452_l9EFoqzFKoIF0vD.tif", 
                                         "Bin_Area", "min_ht", "mean_ht", "max_ht", "sd_ht", "sub_id", "x")]


colnames(canopy_with_bins)[colnames(canopy_with_bins) == "spat_76f434aa2110_30452_l9EFoqzFKoIF0vD.tif"] <- "ht_bin"
colnames(canopy_with_bins)[colnames(canopy_with_bins) == "area"] <- "Area_m"
colnames(canopy_with_bins)[colnames(canopy_with_bins) == "min_ht"] <- "Min_Ht"
colnames(canopy_with_bins)[colnames(canopy_with_bins) == "mean_ht"] <- "Mean_Ht"
colnames(canopy_with_bins)[colnames(canopy_with_bins) == "max_ht"] <- "Max_Ht"
colnames(canopy_with_bins)[colnames(canopy_with_bins) == "sd_ht"] <- "Stdv_Ht"
colnames(canopy_with_bins)[colnames(canopy_with_bins) == "sub_id"] <- "Sub_ID"
colnames(canopy_with_bins)[colnames(canopy_with_bins) == "ht_bin"] <- "Ht_Bin"

st_write(
  canopy_with_bins,
  glue('04_Spatial Analysis/4_Processing/TOW_POLYGONS/TOW_Wales_SS79.gpkg'),
  append = FALSE
)

# End of script
