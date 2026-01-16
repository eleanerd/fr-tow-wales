#################################
# Wales VOM Part 1.2
# Eleanor R. Downer
# 19 November 2025
#################################

###############################
# Import libraries
###############################

import os
import sys
import rasterio
from rasterio import features
import numpy as np
import geopandas as gpd

# Suppress warnings from opening 'wales_osm_clean.gpkg' with multiple layers
import warnings
warnings.filterwarnings("ignore", message="More than one layer found")

###############################
# Settings / Inputs 
###############################

tile_of_interest = sys.argv[1]
print('Processing tile:', tile_of_interest)

wd = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis'

output_path = f'C:\\Users\\eleanor.downer\\OneDrive - Forest Research\\Documents\\TOW_Wales\\hedge_processing\\vom_pt1\\{tile_of_interest}_VOM_with_NFI_NDVI_BUA.tif'
if os.path.exists(output_path):
    print(f'Output for tile {tile_of_interest} already exists. Exiting.')
    sys.exit()

###############################
# Get extent of tile 
###############################

fp = f'{wd}/1_Reference_Data/4_Wales_10km_Squares/Wales_10km_Squares_Planet_overlap.gpkg'
wales_10km_footprint = gpd.read_file(fp)
tile_footprint = wales_10km_footprint[wales_10km_footprint['tileref'] == tile_of_interest]
tile_footprint = tile_footprint.set_crs('EPSG:27700')

if tile_footprint.empty:
    print(f'Tile {tile_of_interest} not found in Wales 10km squares footprint. Exiting.')
    sys.exit()
    
#########################################
# Load VOM Pt1 CHM with NFI & NDVI masked
#########################################

print('Loading VOM Pt1 rasters')

chm_fp = f'{wd}/4_Processing/VOM_Processing/VOM_Pt1/{tile_of_interest}_VOM_with_NFI_NDVI.tif'

with rasterio.open(chm_fp) as src:

    chm_data = src.read(1).astype('float32')

    out_transform = src.transform
    out_meta = src.meta.copy()
    chm_crs = src.crs
    
    out_shape = (int(chm_data.shape[0]), int(chm_data.shape[1]))

###########################################
# Mask out built up areas
###########################################

print('Masking out built up areas')

built_up_areas_fp = f'{wd}\\1_Reference_Data\\13_OS_Built_Up_Areas\\OS_Open_Built_Up_Areas_GeoPackage\\os_open_built_up_areas.gpkg'
built_up_areas = gpd.read_file(built_up_areas_fp, layer='os_open_built_up_areas', bbox=tile_footprint)
built_up_areas = built_up_areas.set_crs(chm_crs)
built_up_areas["geometry"] = built_up_areas.buffer(-10)

if built_up_areas.empty:
    print('No built up areas found in tile, skipping cable masking.')
else:
    if built_up_areas.crs != chm_crs:
        built_up_areas = built_up_areas.to_crs(chm_crs)

    # Rasterise cables
    built_up_areas_mask = features.rasterize(
        [(geom, 1) for geom in built_up_areas.geometry],
        out_shape=out_shape,
        transform=out_transform,
        fill=0,
        dtype="uint8"
    )

    built_up_areas_mask = built_up_areas_mask.astype(bool)

    chm_data[built_up_areas_mask] = np.nan

####################################
# Save to file - without NFI
####################################

print('Saving masked raster to file')

out_meta.update({
    'width' : out_shape[1],
    'height' : out_shape[0],
    'crs' : chm_crs,
    'transform' : out_transform,
    "dtype": "float32",
    "compress": "LZW",
    "tiled": True
})

chm_data = np.where(np.isnan(chm_data), -9999.0, chm_data)

output_path = f'C:\\Users\\eleanor.downer\\OneDrive - Forest Research\\Documents\\TOW_Wales\\hedge_processing\\vom_pt1\\{tile_of_interest}_VOM_with_NFI_NDVI_BUA.tif'
with rasterio.open(output_path, "w", **out_meta) as dst:
    dst.write(chm_data, 1)

print(f'Processing complete for tile {tile_of_interest}.')