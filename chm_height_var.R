library(terra)
library(glue)
library(plyr)
library(sf)

# read in files
chm_dir <- 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/3_Test_Square'
chm_fns <- list.files(chm_dir, pattern = '*_CHM_2.5_masked_NFI.tif$', full.names=TRUE)

print(chm_fns)

window_size <- 9

window <- matrix(1, nrow = window_size, ncol = window_size)

for (chm_fn in chm_fns) {
  
  chm <- terra::rast(chm_fn)
  
  r_var <- focal(chm, w = window, fun = sd, na.policy = "omit", na.rm = TRUE)
  #r_mean <- focal(chm, w = window_mean, fun = mean, na.policy = "omit", na.rm = TRUE)
  
  filepath <- glue('{chm_dir}/CHM_variation/SS79_var_{window_size}x{window_size}.tif')
  writeRaster(r_var, filepath, overwrite=TRUE)
  
}
