#################################
# Wales VOM Part 2
# Eleanor R. Downer
# 14 October 2025
#################################

###############################
# Import libraries
###############################

import os
import sys
import pandas as pd
import rasterio
from rasterio import features
from rasterio.mask import mask
import numpy as np
import geopandas as gpd
from scipy import ndimage
from skimage import morphology

###############################
# Settings / Inputs 
###############################

tile_of_interest = sys.argv[1] 
wd = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data'

print('Processing tile:', tile_of_interest)

output_path = f'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/4_Processing/VOM_Processing/VOM_Pt2/{tile_of_interest}_VOM.tif'

if os.path.exists(output_path):
    print(f'Output for tile {tile_of_interest} already exists. Exiting.')
    sys.exit()

###############################
# Get extent of tile 
###############################

fp = f'{wd}/4_Wales_10km_Squares/Wales_10km_Squares.shp'
wales_10km_footprint = gpd.read_file(fp)

tile_footprint = wales_10km_footprint[wales_10km_footprint['tileref'] == tile_of_interest]
tile_footprint = tile_footprint.set_crs('EPSG:27700')

#############################################
# Load in CHM
#############################################

chm_path = f'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/4_Processing/VOM_Processing/VOM_Pt1/{tile_of_interest}_VOM_with_NFI.tif'

with rasterio.open(chm_path) as src:
    chm_data = src.read(1).astype('float32')
    chm_transform = src.transform
    chm_meta = src.meta
    chm_crs = src.crs

#############################################
# Mask out hedges
#############################################

print('Masking out hedges')

hedge_gpkg_fp = f'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/4_Processing/VOM_Processing/Hedges/{tile_of_interest}_hedges_dissolved.gpkg'
hedge_gdf = gpd.read_file(hedge_gpkg_fp, bbox=tile_footprint)

# Reproject hedges to match CHM
if hedge_gdf.crs != chm_crs:
    hedge_gdf = hedge_gdf.to_crs(chm_crs)

# Rasterise hedges
hedge_mask = features.rasterize(
    [(geom, 1) for geom in hedge_gdf.geometry],
    out_shape=chm_data.shape, #out_shape=out_shape,
    transform=chm_transform, #transform=out_transform,
    fill=0,
    dtype="uint8"
)

# Mask hedges out of VOM
hedge_mask = hedge_mask.astype(bool)
chm_data[hedge_mask] = np.nan

################################################
# Mask out large, flat objects
# (flat buildings/structures and gorse/shrub)
################################################

print('Masking out flat top objects')

chm_var_path = f"{wd}/0_VOM/with_nfi/height_var/{tile_of_interest}_height_var.tif"

with rasterio.open(chm_var_path) as var_src:
    chm_var_data = var_src.read(1).astype(np.float16)
    chm_var_data = np.where(chm_var_data > 0.3, np.nan, chm_var_data) # Set pixels over 0.3 to NaN

    var_meta = var_src.meta.copy()
    chm_crs = var_src.crs

mask = ~np.isnan(chm_var_data) # Create binary mask of valid pixels
structure = ndimage.generate_binary_structure(2, 1) # Label connected components
labeled, num = ndimage.label(mask, structure=structure)
sizes = np.bincount(labeled.ravel()) # Count component sizes

# Build keep mask
keep_mask = np.zeros_like(sizes, dtype=bool)
keep_mask[sizes >= 100] = True
keep_mask[0] = False

chm_var_data = np.where(keep_mask[labeled], chm_var_data, np.nan)
mask = (~np.isnan(chm_var_data))

# Dilate mask by 3x3 window
selem = morphology.disk(1)  # Create a disk-shaped element for erosion
dilated_mask = morphology.dilation(mask, selem)

# Mask out low variation areas from VOM
chm_data[dilated_mask] = np.nan

####################################
# NDVI Thresholding (0.2)
####################################

print('NDVI thresholding')

ndvi_fp = f'{wd}/8_NDVI_and_NDWI/ndvi/{tile_of_interest}_ndvi_composite_resampled.tif'
with rasterio.open(ndvi_fp) as src:
    ndvi_array = src.read(1).astype('float32')

ndvi_array[ndvi_array == -9999] = np.nan
ndvi_mask = ndvi_array < 0.2

chm_data[ndvi_mask] = np.nan

####################################
# NDVI Thresholding (2023 Update)
####################################

thresh = 0.5
mask_data = np.ones(ndvi_array.shape)
mask_data[ndvi_array <= thresh] = 0 # Set pixels with mean yearly NDVI greater than threshold to 1
mask_data = mask_data.astype(np.int16)  # Should print [0, 1]

shrink_size = 3  # Number of pixels to shrink by
selem = morphology.disk(shrink_size)  # Create a disk-shaped element for erosion
shrunken_mask = morphology.dilation(mask_data, selem)

#####################################
# Apply to pixel clusters over 175m2
#####################################

# Create small pixel clump mask

mask = ~np.isnan(chm_data) # Create binary mask of valid pixels
structure = ndimage.generate_binary_structure(2, 1) # Label connected components
labeled, num = ndimage.label(mask, structure=structure)
sizes = np.bincount(labeled.ravel()) # Count component sizes

keep_mask = np.zeros_like(sizes, dtype=bool) # Build keep mask
keep_mask[sizes <= 175] = True # Pixel clumps smaller than 175m2
keep_mask[0] = False 

# Apply 'shrunken_mask' to chm_data
chm_data_ndvi_masked = np.where(shrunken_mask, chm_data, np.nan)
# Add back data from small pixel clump
chm_data_ndvi_masked = np.where(keep_mask[labeled], chm_data, chm_data_ndvi_masked)

####################################
# Remove pixels lower than 2.5m
####################################

#chm_data[chm_data < 2.5] = np.nan
chm_data_ndvi_masked[chm_data_ndvi_masked < 2.5] = np.nan

####################################
# Mask out NFI Woodland
####################################

print('Masking out NFI Woodland from CHM')

NFI_woodland_fp = f'{wd}/6_Wales_NFI_2024/NFI_Wales_2024_WoodlandOnly.gpkg'
nfi_gdf = gpd.read_file(NFI_woodland_fp, bbox=tile_footprint)

# Reproject NFI to match CHM 
if nfi_gdf.crs != chm_crs:
    nfi_gdf = nfi_gdf.to_crs(chm_crs)

# Rasterise NFI
nfi_mask = features.rasterize(
    [(geom, 1) for geom in nfi_gdf.geometry],
    out_shape=chm_data.shape, #out_shape=out_shape,
    transform=chm_transform, #transform=out_transform,
    fill=0,
    dtype="uint8"
)

# Mask NFI out of VOM
nfi_mask = nfi_mask.astype(bool)
#chm_data[nfi_mask] = np.nan
chm_data_ndvi_masked[nfi_mask] = np.nan

####################################
# Remove Power Cables (OpenStreetMap)
# Missed in Pt 1
####################################

power_lines_fp = f'{wd}/11_OpenStreetMap/osm_power_lines.gpkg'
power_lines_gdf = gpd.read_file(power_lines_fp, bbox=tile_footprint)
minor_lines_gdf = power_lines_gdf[power_lines_gdf['power'].isin(['minor_line'])]

# Reproject to match CHM
if minor_lines_gdf.crs != chm_crs:
    minor_lines_gdf = minor_lines_gdf.to_crs(chm_crs)

if minor_lines_gdf.empty:
    print('No minor power lines to mask')
else:
    # Rasterise power lines
    power_lines_mask = features.rasterize(
        [(geom, 1) for geom in minor_lines_gdf.geometry],
        out_shape=chm_data.shape, #out_shape=out_shape,
        transform=chm_transform, #transform=out_transform,
        fill=0,
        dtype="uint8"
    )

    # Mask power lines out of VOM
    power_lines_mask = power_lines_mask.astype(bool)
    chm_data_ndvi_masked[power_lines_mask] = np.nan

####################################
# Remove clusters smaller than 5m2
####################################

# Create binary mask of valid pixels
mask = ~np.isnan(chm_data_ndvi_masked)

structure = ndimage.generate_binary_structure(2, 1) # Label connected components
labeled, num = ndimage.label(mask, structure=structure)
sizes = np.bincount(labeled.ravel()) # Count component sizes

# Build keep mask
keep_mask = np.zeros_like(sizes, dtype=bool)
keep_mask[sizes >= 5] = True
keep_mask[0] = False  # background always False

# Apply filter
chm_data_ndvi_masked = np.where(
    keep_mask[labeled], 
    chm_data_ndvi_masked, 
    np.nan
)
####################################
# Save to file
####################################

print('Saving to file')

chm_meta.update({
    'width' : (chm_data.shape)[1], # out_shape[1],
    'height' : (chm_data.shape)[0], # out_shape[0],
    'crs' : chm_crs,
    'transform' : chm_transform # out_transform,
})

chm_data_ndvi_masked = np.where(
    np.isnan(chm_data_ndvi_masked),
    -9999.0,
    chm_data_ndvi_masked
)

with rasterio.open(output_path, "w", **chm_meta) as dst:
    dst.write(chm_data_ndvi_masked, 1)