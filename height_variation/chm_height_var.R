#################################
# Create height variation raster
#################################

library(terra)
library(glue)
library(plyr)
library(sf)

tile_of_interest <- commandArgs(trailingOnly = TRUE)
wd <- "Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/3_Test_Square"

# read in files
chm_dir <- glue("{wd}/VOM/{tile_of_interest}_VOM_with_NFI.tif")
chm_fns <- list.files(chm_dir, pattern = "*_VOM_with_NFI_NDVI_masked.tif$", full.names=TRUE)

window_size <- 3
window <- matrix(1, nrow = window_size, ncol = window_size)

for (chm_fn in chm_fns) {

  chm <- terra::rast(chm_fn)
  r_var <- focal(chm, w = window, fun = sd, na.policy = "omit", na.rm = TRUE)

  filepath <- glue("{chm_dir}/height_var/SS79_hv3x3.tif")
  writeRaster(r_var, filepath, overwrite=TRUE)

}
