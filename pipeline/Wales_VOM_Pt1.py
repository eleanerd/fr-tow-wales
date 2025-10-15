#################################
# Wales VOM Part 1
# Eleanor R. Downer
# 14 October 2025
#################################

###############################
# Import libraries
###############################

#import os
import sys
import pandas as pd
import rasterio
from rasterio import features
from rasterio.mask import mask
#from rasterio.crs import CRS
from rasterio.windows import from_bounds
#from shapely.geometry import shape
import numpy as np
import geopandas as gpd
from scipy import ndimage
from skimage import morphology

###############################
# Settings / Inputs 
###############################

tile_of_interest = sys.argv[1]
wd = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/3_Test_Square'

###############################
# Get extent of tile 
###############################

fp = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/4_Wales_10km_Squares/Wales_10km_Squares.shp'
wales_10km_footprint = gpd.read_file(fp)

tile_footprint = wales_10km_footprint[wales_10km_footprint['tileref'] == tile_of_interest]
tile_footprint = tile_footprint.set_crs('EPSG:27700')

#####################################
# Remove pixels below 1.3m from CHM
#####################################

print('Removing pixels below 1.3m from CHM')

chm_fp = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/1_Wales_LiDAR/Wales_CHM_2020_22_32bit.tif'

with rasterio.open(chm_fp) as src:

    window = from_bounds(*tile_footprint.total_bounds, transform=src.transform)
    chm_data = src.read(1, window=window).astype('float32')

    out_transform = src.window_transform(window)
    out_meta = src.meta.copy()
    chm_crs = src.crs
    
    chm_data = np.where(chm_data < 1.3, np.nan, chm_data)
    
out_shape = (int(window.height), int(window.width))

#################################
# Create dictionary OSMM terms
#################################

osm_desc_terms = pd.read_csv('Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/5_Wales_OSMM_2023/descterm_categories.csv')
osm_desc_terms.columns = osm_desc_terms.columns.str.strip()
osm_desc_terms = osm_desc_terms.map(lambda x: x.strip() if isinstance(x, str) else x)

osm_terms = {
    "tree_terms": osm_desc_terms.loc[osm_desc_terms["category"] == "tree_term", "descterm"].tolist(),
    "mask_terms": osm_desc_terms.loc[osm_desc_terms["category"] == "mask_term", "descterm"].tolist(),
    "structure_terms": osm_desc_terms.loc[osm_desc_terms["category"] == "structure_term", "descterm"].tolist(),
    "water_terms": osm_desc_terms.loc[osm_desc_terms["category"] == "water_body", "descterm"].tolist()
}

for k, v in osm_terms.items():
    print(f"{k}: {v}")

##############################
# Categorise OSMM
##############################

print('Categorising OSMM features')

OSMM_path = f'{wd}/OSMM_Square_{tile_of_interest}/OSMM_Square_{tile_of_interest}.shp'
OSMM = gpd.read_file(OSMM_path)

if OSMM.crs != chm_crs:
    OSMM = OSMM.to_crs(chm_crs)

# Update features with blank 'descterm' field
OSMM["descterm"].replace('', np.nan, inplace=True)

for value in OSMM["descgroup"].unique():
    OSMM.loc[
        (OSMM["descterm"].isna()) & (OSMM["descgroup"] == value),
        "descterm"
    ] = value

# Create Surf_Obj field
OSMM["Surf_Obj"] = None

# Assign Surf_Obj values based on descterm categories
OSMM.loc[OSMM["descterm"].isin(osm_terms['tree_terms']), "Surf_Obj"] = "No"
OSMM.loc[OSMM["descterm"].isin(osm_terms['structure_terms']), "Surf_Obj"] = "Yes"
OSMM.loc[OSMM["descterm"].isin(osm_terms['mask_terms']), "Surf_Obj"] = "Mask"
OSMM.loc[OSMM["descterm"].isin(osm_terms['water_terms']), "Surf_Obj"] = "Water"

######################################
# Mask water bodies
######################################

print('Masking water bodies')

water_gdf = OSMM[OSMM["Surf_Obj"] == "Water"].copy()
water_gdf = water_gdf.to_crs(chm_crs)

# Rasterise mask_terms
water_mask = features.rasterize(
    [(geom, 1) for geom in water_gdf.geometry],
    out_shape=out_shape,
    transform=out_transform,
    fill=0,
    dtype="uint8"
)

# Apply water mask to CHM
chm_data = np.where(water_mask == 1, np.nan, chm_data)

###########################################
# Mask out non-tree OSMM features from CHM
###########################################

print('Masking out non-tree OSMM features from CHM')

# Buffer and reproject surface_terms
surface_gdf = OSMM[OSMM["Surf_Obj"] == "Yes"].copy()
surface_gdf = surface_gdf.to_crs(chm_crs)
surface_gdf["geometry"] = surface_gdf.buffer(2)

# Get mask_terms and reproject
mask_gdf = OSMM[OSMM["Surf_Obj"] == "Mask"].copy()
mask_gdf = mask_gdf.to_crs(chm_crs)

# Rasterise buffered surface_terms
surface_mask = features.rasterize(
        [(geom, 1) for geom in surface_gdf.geometry],
        out_shape=out_shape,
        transform=out_transform,
        fill=0,
        dtype="uint8"
)

# Rasterise mask_terms
mask_mask = features.rasterize(
        [(geom, 1) for geom in mask_gdf.geometry],
        out_shape=out_shape,
        transform=out_transform,
        fill=0,
        dtype="uint8"
)

# Apply both masks
combined_mask = (surface_mask == 1) | (mask_mask == 1)
chm_data[combined_mask] = np.nan

###########################################
# Mask out buildings from OpenStreetMap
###########################################

print('Masking out buildings from OpenStreetMap')

osm_fp = f'C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/wales_osm_clean.gpkg'
osm_polygons = gpd.read_file(osm_fp, layer='multipolygons', bbox=tile_footprint)
osm_buildings = osm_polygons[osm_polygons['building'].notna()].copy()
osm_buildings = osm_buildings.to_crs(chm_crs)
osm_buildings["geometry"] = osm_buildings.buffer(2)

# Rasterise buildings
osm_buildings_mask = features.rasterize(
    [(geom, 1) for geom in osm_buildings.geometry],
    out_shape=out_shape,
    transform=out_transform,
    fill=0,
    dtype="uint8"
)

osm_buildings_mask = osm_buildings_mask.astype(bool)
chm_data[osm_buildings_mask] = np.nan

###########################################
# Mask power cables
###########################################

print('Masking out power cables from CHM')

cables_fp = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/10_Power_Lines/OHL_buffered.gpkg'
cables = gpd.read_file(cables_fp)
cables = cables.set_crs(chm_crs) 

# Rasterise cables
cables_mask = features.rasterize(
    [(geom, 1) for geom in cables.geometry],
    out_shape=out_shape,
    transform=out_transform,
    fill=0,
    dtype="uint8"
)

cables_mask = cables_mask.astype(bool)

over_10m_mask = chm_data > 10
cables_mask = cables_mask & over_10m_mask

chm_data[cables_mask] = np.nan

###########################################
# Mask power structures from OpenStreetMap
###########################################

print('Masking out power structures from CHM')

power_structures_fp = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/11_OpenStreetMap/osm_power_structures_v2.gpkg'
power_structures = gpd.read_file(power_structures_fp, bbox=tile_footprint)

power_structures = power_structures[power_structures['power_source'] != 'solar'].copy()

if power_structures.crs != chm_crs:
    power_structures = power_structures.to_crs(chm_crs)

power_structures_mask = features.rasterize(
    [(geom, 1) for geom in power_structures.geometry],
    out_shape=out_shape,
    transform=out_transform,
    fill=0,
    dtype="uint8"
)

power_structures_mask = power_structures_mask.astype(bool)
chm_data[power_structures_mask] = np.nan

####################################
# Save to file - with NFI still
####################################

out_meta.update({
    'width' : out_shape[1],
    'height' : out_shape[0],
    'crs' : chm_crs,
    'transform' : out_transform,
})

chm_data = np.where(np.isnan(chm_data), -9999.0, chm_data)

output_path = f'{wd}/VOM/{tile_of_interest}_VOM_with_NFI.tif'
with rasterio.open(output_path, "w", **out_meta) as dst:
    dst.write(chm_data, 1)

####################################
# NDVI Thresholding
####################################

chm_data_ndvi = chm_data.copy()

ndvi_fp = f'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/3_Test_Square/NDVI/{tile_of_interest}_ndvi_composite_resampled.tif'
with rasterio.open(ndvi_fp) as src:
    ndvi_array = src.read(1).astype('float32')

ndvi_array[ndvi_array == -9999] = np.nan
ndvi_mask = ndvi_array < 0.2

chm_data_ndvi[ndvi_mask] = np.nan

thresh = 0.35
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

mask = ~np.isnan(chm_data_ndvi) # Create binary mask of valid pixels
structure = ndimage.generate_binary_structure(2, 1) # Label connected components
labeled, num = ndimage.label(mask, structure=structure)
sizes = np.bincount(labeled.ravel()) # Count component sizes

keep_mask = np.zeros_like(sizes, dtype=bool) # Build keep mask
keep_mask[sizes <= 250] = True # Pixel clumps smaller than 250m2
keep_mask[0] = False 

# Apply 'shrunken_mask' to chm_data
chm_data_ndvi_masked = np.where(shrunken_mask, chm_data_ndvi, np.nan)
# Add back data from small pixel clump
chm_data_ndvi_masked = np.where(keep_mask[labeled], chm_data_ndvi, chm_data_ndvi_masked)

####################################
# Save to file - without NFI
####################################

out_meta.update({
    'width' : out_shape[1],
    'height' : out_shape[0],
    'crs' : chm_crs,
    'transform' : out_transform,
})

chm_data_ndvi_masked = np.where(np.isnan(chm_data_ndvi_masked), -9999.0, chm_data_ndvi_masked)

output_path = f'{wd}/VOM/{tile_of_interest}_VOM_with_NFI_NDVI_masked.tif'
with rasterio.open(output_path, "w", **out_meta) as dst:
    dst.write(chm_data_ndvi_masked, 1)