
library(sf)
library(dplyr)

rm(list = ls())

# File paths
gpkg_fn <- "C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/Results/TOW_Wales.gpkg"
out_fn <- "C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/Results/TOW_Wales/TOW_Wales.gpkg"

# Columns to keep
#columns_to_keep <- c('area_m2','woodland_type','tow_id','focal_mean','bin_area',
#                     'min_ht','mean_ht','max_ht','sd_ht','sub_id','lidar_year','lidar_date', 'Tilename')

# Read only those columns (sf passes to GDAL)
gpkg <- st_read(gpkg_fn)

gpkg <- gpkg %>% st_cast("POLYGON")

# Rename
# names(gpkg) <- c("Lidar_Date", "Lidar_Year", "Sub_ID", "Stdv_Height", "Height_Class", "Mean_Height", "TOW_Area", "Woodland_Type",
#                  "Max_Height", "Min_Height", "Bin_Area", "Tilename", "TOW_ID", "geom" )


#st_geometry(gpkg) <- "geometry"

# Reorder if needed (geometry last is conventional)
gpkg <- gpkg[, c("TOW_ID","TOW_Area","Woodland_Type","Height_Class","Bin_Area",
                 "Min_Height","Mean_Height","Max_Height","Stdv_Height","Sub_ID",
                 "Lidar_Year","Lidar_Date", "Tilename", "geom")]

gpkg <- gpkg %>% filter(Bin_Area != 0)

colnames(gpkg)[colnames(gpkg) == "Tilename"] <- "British_Gr"


st_write(gpkg, out_fn, delete_dsn = TRUE)

summary <- gpkg %>%
  group_by(Woodland_Type) %>%
  summarise(total_Bin_Area = sum(Bin_Area, na.rm = TRUE))