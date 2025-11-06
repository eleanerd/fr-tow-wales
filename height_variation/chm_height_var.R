#################################
# Create height variation raster
#################################

library(terra)
library(glue)
library(plyr)
library(sf)

wd <- "Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/0_VOM/with_nfi"

# read in files
chm_fns <- list.files(wd, pattern = "*with_NFI.tif$", full.names=TRUE)

window_size <- 3
window <- matrix(1, nrow = window_size, ncol = window_size)

for (chm_fn in chm_fns) {
  
  tile <- strsplit(basename(chm_fn), "_")[[1]][[1]]
  filepath <- glue("{wd}/height_var/{tile}_height_var.tif")
  
  if (file.exists(filepath) {
    next
  }
  
  chm <- terra::rast(chm_fn)
  r_var <- focal(chm, w = window, fun = sd, na.policy = "omit", na.rm = TRUE)

  writeRaster(r_var, filepath, overwrite=TRUE)

}
 