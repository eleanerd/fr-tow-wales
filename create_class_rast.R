####################################
# Freddie Hunter
# 18 Sept 2025
# Wales TOW - VOM Delineation
####################################

rm(list = ls())
gc()

#############################################################
#Libraries
library(terra)
library(dplyr)
library(sf)
library(raster)
library(glue)
library(deldir)
library(dplyr)

#############################################################
# Input files
# There are many parameters imbedded in the code. 
#wd <- "//forestresearch.gov.uk/shares/IFOS/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/3_Test_Square"
#raz_files <- 'SS79_CHM_1.3m_masked_NFI_OSSM_Cables.tif'
wd <- 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/3_Test_Square'
raz_files <- glue('{wd}/VOM/SS79_VOM_final.tif')

#################################################
#Processing

# Load CHM raster
chm_full <- rast(raz_files)
chm_full[chm_full < 2.5] <- NA

ncol_chm <- ncol(chm_full)
nrow_chm <- nrow(chm_full)

section_size <- 2000
overlap <- 10

for (row_start in seq(1, nrow_chm, by = section_size - overlap)) {
  for (col_start in seq(1, ncol_chm, by = section_size - overlap)) {
    
    print(glue("row {row_start}, col {col_start}"))
    
    row_end <- min(row_start + section_size - 1, nrow_chm)
    col_end <- min(col_start + section_size - 1, ncol_chm)
    
    # Get cell numbers for corners
    cell_top_left <- cellFromRowCol(chm_full, row_start, col_start)
    cell_bottom_right <- cellFromRowCol(chm_full, row_end, col_end)
    
    # Get coordinates of those cells
    coords_top_left <- xyFromCell(chm_full, cell_top_left)
    coords_bottom_right <- xyFromCell(chm_full, cell_bottom_right)
    
    # Create extent from coordinates
    ext_section <- ext(
      coords_top_left[1], coords_bottom_right[1],
      coords_bottom_right[2], coords_top_left[2]
    )
    
    # Crop section
    chm <- crop(chm_full, ext_section)
    chm[chm < 1.3] <- NA
    
    
    # Step 1: 5x5 m local maximum filter
    print('Step 1: 5x5 m local maximum filter')
    local_max <- terra::focal(chm, w = matrix(1, 5, 5), fun = max, na.rm = T)
    local_max[is.na(chm)] <- NA # set all NA values in original data back to NA to stop expansion by window
    
    # Step 2: Initial height-based classification
    print('Step 2: Initial height-based classification')
    # Reclassification matrix: start, end, new class
    rcl <- matrix(c(
      1.3, 3,      1.3,
      3, 10,       3,
      10, 15,     10,
      15, 20,     15,
      20, 30,    20,
      30, Inf,   30
    ), ncol=3, byrow=TRUE)
    
    local_max <- as(local_max, "SpatRaster")
    class_raster <- classify(local_max, rcl)
    
    rm(local_max)
    gc()
    
    # Step 3: Patch size calculation for class 0.75–1.3 (needs special handling)
    print('Step 3: Patch size calculation for class 0.75–1.3')
    # Extract to process patches but exclude 1.3-3 bin.
    target_patch <- as.numeric(class_raster > 1.3)
    target_patch[target_patch == 0] <- NA # assumes NA gets value of 1
    
    # Use terra's patches function to get connected regions
    patches <- patches(target_patch, directions = 4)
    rm(target_patch)
    gc()
    patchFreq <- terra::freq(patches)
    
    # modal filter
    fw <- matrix(nrow = 3, ncol = 3)
    fw[is.na(fw)] <- 1
    class_raster_mod <- terra::focal(
      x = class_raster,
      w = fw,
      fun = modal,
      na.rm = T
    )
    class_raster_mod[is.na(class_raster)] <- 0
    
    # seive based on islands under 
    writeRaster(class_raster_mod, glue('{wd}/Delineation/SS79_class_raster_mod_{row_start}_{col_start}.tif'), overwrite = TRUE)
    system(glue('"C:\\Program Files\\QGIS 3.44.1\\apps\\Python312\\python" "C:\\Program Files\\QGIS 3.44.1\\apps\\Python312\\Scripts\\gdal_sieve.py" -st 5 -4 -of GTiff "{wd}/Delineation/SS79_class_raster_mod_{row_start}_{col_start}.tif" "{wd}/Delineation/SS79_sieve_5m_class_raster_mod_{row_start}_{col_start}.tif"'))
    
  }
}

#########################################
# Merge rasters with GDAL in R
#########################################

gdal_merge_path <- "C:\\Program Files\\QGIS 3.44.1\\apps\\Python312\\Scripts\\gdal_merge.py"
python_path <- "C:\\Program Files\\QGIS 3.44.1\\apps\\Python312\\python.exe"

output_file <- "Y:\\Forest Inventory\\0700_NonCore_Funded\\0726_TOW_Wales\\04_Spatial Analysis\\3_Test_Square\\Delineation\\SS79_class_raster.tif"

input_pattern <- "Y:\\Forest Inventory\\0700_NonCore_Funded\\0726_TOW_Wales\\04_Spatial Analysis\\3_Test_Square\\Delineation\\SS79_*.tif"
input_files <- Sys.glob(input_pattern)

args <- c(
  gdal_merge_path,
  "-o", shQuote(output_file),
  "-of", "GTiff",
  "-a_nodata", "-9999",
  shQuote(input_files)
)
system2(command = python_path, args = args, wait = TRUE)


# 4 Decide what to do with the remaining tree canopy data.
#     Probably polygonise the remaining clusters of pixels at the class level
#     Give them an ID, a height class, some height attributes (mean, min etc) and a woodland type
#     Woodland types should probably be classifed with similar but improved... criteria as TOW england. 
