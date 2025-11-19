####################
# Freddie Hunter
# 18 Sept 2025
# Wales TOW
###################

# Reviewed 16 Oct 2025 by Eleanor Downer

#################################################
# Settings and Libraries
#################################################

rm(list = ls())
gc()

library(terra)
library(dplyr)
library(sf)
library(glue)
library(deldir)
library(dplyr)
library(centerline)
library(smoothr)
library(raster)
library(foreach)
library(doParallel)
library(lwgeom)
library(purrr)
library(stringr)

#################################################
# Functions
#################################################

# --- Function to compute adjusted straightness
get_adjusted_straightness <- function(geom, penalty = 2) {
  coords <- st_coordinates(geom)
  if (nrow(coords) < 2) return(NA_real_)
  # Start and end coordinates
  start <- coords[1, c("X", "Y")]
  end <- coords[nrow(coords), c("X", "Y")]
  # Euclidean distance
  euclid_dist <- sqrt(sum((end - start)^2))
  # Total line length
  total_length <- as.numeric(st_length(geom))
  if (total_length == 0) return(NA_real_)
  # Base straightness
  straightness <- euclid_dist / total_length
  # Penalized straightness
  adjusted <- straightness^penalty
  return(c(adjusted))
}

get_euclidian_dist <- function(geom) {
  coords <- st_coordinates(geom)
  if (nrow(coords) < 2) return(NA_real_)
  # Start and end coordinates
  start <- coords[1, c("X", "Y")]
  end <- coords[nrow(coords), c("X", "Y")]
  # Euclidean distance
  euclid_dist <- sqrt(sum((end - start)^2))
  # Total line length
  total_length <- as.numeric(st_length(geom))
  if (total_length == 0) return(NA_real_)
  return(euclid_dist)
}

# Function to calculate maximum area
max_rectangle_area <- function(length) {
  width <- 5  # Maximum allowed width
  area <- length * width
  return(area)
}


#################################################
# Processing
#################################################

wd <- "//forestresearch.gov.uk/shares/IFOS/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data"
setwd(wd)

files_list <- list.files("//forestresearch.gov.uk/shares/IFOS/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/4_Processing/VOM_Processing/VOM_Pt1/",
                        pattern = "_VOM_with_NFI_NDVI_BUA.tif$")
tiles_list <- str_remove(files_list, "_VOM_with_NFI_NDVI_BUA.tif")

# Register parallel backend
cl <- makeCluster(detectCores() - 22) # Leave one core free
registerDoParallel(cl)

for (tile_of_interest in tiles_list) {
  
  print(tile_of_interest)
  
  if (file.exists(glue("C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/hedge_processing/hedges/{tile_of_interest}_hedges_v4.gpkg"))) {
    next
  }
  
  if (file.exists(glue("C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/hedge_processing/hedges/{tile_of_interest}_hedges.gpkg"))) {
    next
  }
    
  # Load CHM raster
  chm_ndvi_filt <- rast(glue("//forestresearch.gov.uk/shares/IFOS/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/4_Processing/VOM_Processing/VOM_Pt1/{tile_of_interest}_VOM_with_NFI_NDVI_BUA.tif") )
  chm_ndvi_filt[chm_ndvi_filt < 0.5] <- NA # Set values below 0.5m to NA
  
  ncol_chm <- ncol(chm_ndvi_filt)
  nrow_chm <- nrow(chm_ndvi_filt)
  
  # Section parameters
  section_size <- 1000
  overlap <- 10
    
  # Lists
  hedges_tiles <- list()
  
  # Iterate over sections (~100 per 10km tile)
  for (row_start in seq(1, nrow_chm, by = section_size - overlap)) {
    for (col_start in seq(1, ncol_chm, by = section_size - overlap)) {
      print(glue("row {row_start}, col {col_start}"))
      
      # Define section boundaries
      row_end <- min(row_start + section_size - 1, nrow_chm)
      col_end <- min(col_start + section_size - 1, ncol_chm)
      # Get cell numbers for corners
      cell_top_left <- cellFromRowCol(chm_ndvi_filt, row_start, col_start)
      cell_bottom_right <- cellFromRowCol(chm_ndvi_filt, row_end, col_end)
      # Get coordinates of those cells
      coords_top_left <- xyFromCell(chm_ndvi_filt, cell_top_left)
      coords_bottom_right <- xyFromCell(chm_ndvi_filt, cell_bottom_right)
      # Create extent from coordinates
      ext_section <- ext(
        coords_top_left[1], coords_bottom_right[1],
        coords_bottom_right[2], coords_top_left[2]
      )
      # Crop section
      chm <- crop(chm_ndvi_filt, ext_section)
      #chm_full_crop <- crop(chm_full, ext_section)
  
      # Skip if section is entirely NA
      if (all(is.na(values(chm)))) {
        print("Section is all NA, skipping...")
        next
      }
  
      # Step 1: 5x5 m local maximum filter
      print("Step 1: 5x5 m local maximum filter")
      local_max <- terra::focal(chm, w = matrix(1, 5, 5), fun = max, na.rm = TRUE)
      # set all NA values in original data back to NA to stop expansion by window
      local_max[is.na(chm)] <- NA
      local_max <- as(local_max, "SpatRaster")
  
      # Step 2: Initial height-based classification
      print("Step 2: Initial height-based classification")
      # Reclassification matrix: start, end, new class
      # Using integers for class labels - faster
      rcl <- matrix(c(
        0.49, 6, 1,
        6, Inf, 2
      ), ncol=3, byrow=TRUE)
      class_raster <- classify(local_max, rcl)
  
      #rm(local_max) # free up memory
  
      # Apply modal filter
      fw <- matrix(1, 3, 3)
      class_raster_mod <- terra::focal(
        x = class_raster,
        w = fw,
        fun = modal,
        na.rm = T
      )
      class_raster_mod[is.na(class_raster)] <- 0
      #rm(class_raster) # free up memory
      #plot(class_raster_mod)
  
      # seive based on islands under 5 pixels
      #raster_fp <- glue("0_VOM/Hedges/{tile_of_interest}_class_raster_mod_{row_start}_{col_start}.tif")
      raster_fp <- glue("C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/hedge_processing/{tile_of_interest}_class_raster_mod_{row_start}_{col_start}.tif")
      terra::writeRaster(class_raster_mod, raster_fp, overwrite = TRUE)
  
      # GDAL Sieve
      py_path <- "C:\\Program Files\\QGIS 3.44.1\\apps\\Python312\\python.exe"
      gdal_sieve <- "C:\\Program Files\\QGIS 3.44.1\\apps\\Python312\\Scripts\\gdal_sieve.exe"
      #sieved_raster_fp <- glue("0_VOM/Hedges/{tile_of_interest}_sieve_5m_class_raster_mod_{row_start}_{col_start}.tif")
      sieved_raster_fp <- glue("C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/hedge_processing/{tile_of_interest}_sieve_5m_class_raster_mod_{row_start}_{col_start}.tif")
      system(glue('"{gdal_sieve}" -st 5 -8 -of GTiff "{raster_fp}" "{sieved_raster_fp}"'))
      
      #rm(class_raster_mod) # free up memory
      
      filt_max_class <- terra::rast(sieved_raster_fp)
      filt_max_class[filt_max_class == 0] <- NA
      filt_max_class[filt_max_class != 0] <- 1
      
      canopy_area <- stars:::st_as_sf.stars(
        stars::st_as_stars(filt_max_class), 
        point = FALSE, 
        merge = TRUE, 
        connect8 = TRUE) %>% 
        st_buffer(0.01) %>% 
        st_union() %>% 
        st_cast('POLYGON') %>% 
        st_as_sf() 
      
      # Compute polygon areas for all canopy polygons
      canopy_area$area_m <- as.numeric(st_area(canopy_area))
      canopy_area <- canopy_area %>% filter(area_m >= 20)
      
      # remove pixels over 6m
      filt6 <- terra::rast(sieved_raster_fp)
      msk6 <- filt6 == 2
      msk6[isFALSE(msk6)] <- NA 
      mask_poly <- as.polygons(msk6) %>% st_as_sf() %>% st_buffer(0.2) # this is used later
      
      file.remove(raster_fp)
      file.remove(sieved_raster_fp)
      
      print("Finding hedges in polygons...")
  
      # Find hedges{
      hedges_list <- foreach(
        pol = seq_along(canopy_area$x),
        .packages = c("sf",
                      "smoothr",
                      "centerline",
                      "dplyr",
                      "lwgeom",
                      "purrr"),
        .errorhandling = "pass"
      ) %dopar% {
  
        poly <- canopy_area$x[pol]
  
        # Buffer inward by a small amount
        shrunk <- st_buffer(poly, -2.49, endCapStyle = "FLAT")
        # Buffer back outward
        regrown <- st_buffer(shrunk, 2.6, endCapStyle = "FLAT")
  
        # Subtract to get narrow border regions
        narrow_parts1 <- st_difference(poly, regrown) %>%
          st_union() %>%
          st_cast("POLYGON") %>%
          st_as_sf()
  
        # Skip if no narrow parts
        if (isTRUE(dim(narrow_parts1)[1] == 0)) {
          return(NULL)
        }
  
        narrow_parts1$id <- c(1:nrow(narrow_parts1))
        narrow_parts1$area_m <- as.numeric(st_area(narrow_parts1))
        narrow_parts1 <- narrow_parts1 %>% filter(area_m >= 20)
  
        # Skip if no narrow parts
        if (isTRUE(dim(narrow_parts1)[1] == 0)) {
          return(NULL)
        }
        ## What if mask_poly is empty? Aka there's no polygons over 6m in section
  
        # Remove pixels over 6m
        poly_cleaned <- st_difference(narrow_parts1, mask_poly) %>%
          st_union() %>% 
          st_cast("POLYGON")
  
        # Skip if no narrow parts
        if (isTRUE(dim(poly_cleaned)[1] == 0)) {
          return(NULL)
        }
  
        # Rerun previous part
        shrunk2 <- st_buffer(poly_cleaned, -2.49, endCapStyle = "FLAT") # Buffer inward by a small amount
        regrown2 <- st_buffer(shrunk2, 2.6, endCapStyle = "FLAT") # Buffer back outward
  
        # Subtract to get narrow border regions
        narrow_parts2 <- st_difference(poly_cleaned, regrown2) %>%
          st_union() %>%
          st_cast("POLYGON") %>%
          st_as_sf()
  
        # Skip if no narrow parts
        if (isTRUE(dim(narrow_parts2)[1] == 0)) {
          return(NULL)
        }
  
        narrow_parts2$id <- c(1:nrow(narrow_parts2))
        narrow_parts2$area_m <- narrow_parts2 %>%
          st_area() %>%
          units::drop_units()
        narrow_parts2 <- narrow_parts2 %>% filter(area_m >= 20)
  
        if (isTRUE(dim(narrow_parts2)[1] == 0)) {
          return(NULL)
        }
        #plot(narrow_parts$x, add = T, col = 'blue')
  
        # Smooth the polygon
        narrow_parts <- smoothr::smooth(
          narrow_parts2,
          method = "ksmooth",
          smoothness = 15
        ) %>%
          st_make_valid()
  
        narrow_parts <- st_simplify(
          narrow_parts,
          dTolerance = 0.05,
          preserveTopology = TRUE
        )
  
        # Create skeleton of polygon
        skel_dense <- centerline::cnt_skeleton(narrow_parts, keep = 1.5)
  
        # Was getting error when skel_dense was made of multiple geometries
        if (nrow(skel_dense) == 0) {
          return(NULL)
        }
  
        # Get centerline or return NULL if fails
        skel_dense <- st_make_valid(skel_dense)
        plot(skel_dense)
        skel_cent <- centerline::cnt_path_guess(
          input = narrow_parts,
          skeleton = skel_dense
        )
  
        if (nrow(skel_cent) == 0) {
          return(NULL)
        }
        
        # Filter out short lines
        skel_cent$cnt_len_adj <- as.numeric(st_length(skel_cent)) * 0.92
        skel_cent <- skel_cent %>% filter(cnt_len_adj >= 20)
        if (dim(skel_cent)[1] == 0) {
          return(NULL)
        }
  
        # Ggt straightness, auclid and length
        skel_cent <- skel_cent %>%
          mutate(
            straightness = purrr::map_dbl(geometry, get_adjusted_straightness),
            euclid = purrr::map_dbl(geometry, get_euclidian_dist)
          )
  
        # Filter out wiggly lines
        filtered_skel_cent <- skel_cent %>%
          filter(straightness > 0.35)
        if (dim(filtered_skel_cent)[1] == 0) {
          return(NULL)
        }
  
        other_dat <- filtered_skel_cent %>% st_drop_geometry()
        
        # Select hedges that are straight
        hedges_sel <- filter(narrow_parts, id %in% filtered_skel_cent$id)
        
        # Minimum bounding circle
        mbc_a <- lwgeom::st_minimum_bounding_circle(hedges_sel) %>%
          st_area() %>%
          units::drop_units()
        mbc_p <- hedges_sel$area_m/mbc_a * 100
        hedges_sel$mbcp <- mbc_p
        # Remove 'compact' polygons; those with mbcp > 30%
        hedges_sel <- hedges_sel[hedges_sel$mbcp < 30, ]
        
        hedges_sel <- merge(hedges_sel, other_dat, by = "id")
        
        # Get max rectangle area
        hedges_sel$max_rec_area  <-max_rectangle_area(hedges_sel$cnt_len_adj) 
  
        # Computes polygon area as a percentage of the estimated max rectangle area
        # Remove polygons where area exceeds the rectangle (likely 'fat' polygons)
        hedges_sel$max_area_perc <- hedges_sel$area_m.x / hedges_sel$max_rec_area * 100
        hedges_sel <- hedges_sel[hedges_sel$max_area_perc < 100, ]
        
        # Calculate length-area-ratio
        hedges_sel$lar = st_length(st_boundary(hedges_sel)) / st_area(hedges_sel) %>%units::drop_units()
  
        #hedges <- rbind(hedges, hedges_sel)
        return(hedges_sel)
      }
      
      # Combine only valid sf results
      hedges_list <- Filter(function(x) inherits(x, "sf") && nrow(x) > 0, hedges_list)
      
      if (length(hedges_list) > 0) {
        hedges <- do.call(rbind, hedges_list)
      } else {
        hedges <- NULL
      }
    
      # Save outputs to lists
      if (!is.null(hedges) && nrow(hedges) > 0) {
        
        hedges_tiles[[length(hedges_tiles) + 1]] <- hedges
        
      } else {
        print(glue("No hedges, skipping section"))
      }
      
      rm(hedges)
      
    }
  }
  
  # Merge all hedges
  merged_hedges <- do.call(rbind, hedges_tiles)
  
  if (is.null(merged_hedges)) {
    next
  }
  # Update attribute values to have 2 decimal places
  merged_hedges <- merged_hedges %>% mutate(across(where(is.numeric), ~ round(., 2)))
  st_write(merged_hedges, glue("C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/hedge_processing/hedges/{tile_of_interest}_hedges_v4.gpkg"), delete_dsn = TRUE)
}

# Clean up cluster after parallel work
stopCluster(cl)
registerDoSEQ()