#################################
# Wales VOM Part 1
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
from rasterio.windows import from_bounds
from rasterio.warp import reproject, Resampling
import numpy as np
import geopandas as gpd
from scipy import ndimage
from skimage import morphology

# Suppress warnings from opening 'wales_osm_clean.gpkg' with multiple layers
import warnings
warnings.filterwarnings("ignore", message="More than one layer found")

# To Do
# Don't rasterize if gdf is empty - apply to all masks 

###############################
# Settings / Inputs 
###############################

tile_of_interest = sys.argv[1] 
wd = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis'

print('Processing tile:', tile_of_interest)

output_path = f'{wd}/1_Reference_Data/0_VOM/with_nfi/{tile_of_interest}_VOM_with_NFI_NDVI.tif'

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
    
#####################################
# Remove pixels below 1.3m from CHM
#####################################

print('Removing pixels below 0.5m from CHM')

chm_fp = f'C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/data/{tile_of_interest}_CHM.tif'

with rasterio.open(chm_fp) as src:

    #window = from_bounds(*tile_footprint.total_bounds, transform=src.transform)
    chm_data = src.read(1).astype('float32')

    out_transform = src.transform
    out_meta = src.meta.copy()
    chm_crs = src.crs
    
    out_shape = (int(chm_data.shape[0]), int(chm_data.shape[1]))

    chm_data = np.where(chm_data < 0.5, np.nan, chm_data)

#####################################
# Remove pixels not on land
#####################################

wales_boundary_fp = f'{wd}/1_Reference_Data/2_Wales_Boundary/Wales_Boundary.shp'
wales_boundary = gpd.read_file(wales_boundary_fp)

if wales_boundary.crs != chm_crs:
    wales_boundary = wales_boundary.to_crs(chm_crs)

# Mask CHM to Wales boundary
print('Masking CHM to Wales boundary')

land_mask = features.rasterize(
    [(geom, 1) for geom in wales_boundary.geometry],
    out_shape=out_shape,
    transform=out_transform,
    fill=0,
    dtype="uint8"
)

# Apply land mask to CHM
chm_data = np.where(land_mask == 0, np.nan, chm_data)

if np.all(np.isnan(chm_data)):
    print('No valid CHM data remains, skipping tile.')
    sys.exit()

#################################
# Skip tile if no data remains
#################################

if np.all(np.isnan(chm_data)):
    print('No valid CHM data remains, skipping tile.')
    sys.exit()

#################################
# Create dictionary OSMM terms
#################################

osm_desc_terms = pd.read_csv(f'{wd}/1_Reference_Data/5_Wales_OSMM_2023/descterm_categories.csv')
osm_desc_terms.columns = osm_desc_terms.columns.str.strip()
osm_desc_terms = osm_desc_terms.map(lambda x: x.strip() if isinstance(x, str) else x)

osm_terms = {
    "tree_terms": osm_desc_terms.loc[osm_desc_terms["category"] == "tree_term", "descterm"].tolist(),
    "mask_terms": osm_desc_terms.loc[osm_desc_terms["category"] == "mask_term", "descterm"].tolist(),
    "structure_terms": osm_desc_terms.loc[osm_desc_terms["category"] == "structure_term", "descterm"].tolist(),
    "water_terms": osm_desc_terms.loc[osm_desc_terms["category"] == "water_body", "descterm"].tolist()
}

##############################
# Categorise OSMM
##############################

print('Categorising OSMM features')

OSMM_path = f'{wd}/1_Reference_Data/5_Wales_OSMM_2023/GPKG/OS_Map_{tile_of_interest}.gpkg'
OSMM = gpd.read_file(OSMM_path)

if OSMM.crs != chm_crs:
    OSMM = OSMM.to_crs(chm_crs)

# Update features with blank 'descterm' field
#OSMM["descterm"].replace('', np.nan, inplace=True)
OSMM["descterm"] = OSMM["descterm"].replace('', np.nan)

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
if water_gdf.empty:
    water_mask = np.zeros(out_shape, dtype="uint8")
else:
    water_mask = features.rasterize(
        [(geom, 1) for geom in water_gdf.geometry],
        out_shape=out_shape,
        transform=out_transform,
        fill=0,
        dtype="uint8"
    )

    # Apply water mask to CHM
    chm_data = np.where(water_mask == 1, np.nan, chm_data)

    if np.all(np.isnan(chm_data)):
        print('No valid CHM data remains, skipping tile.')
        sys.exit()

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

if surface_gdf.empty:
    print('No surface_terms found in tile, skipping surface_term masking.')
    surface_mask = np.zeros(out_shape, dtype="uint8")
else:
    # Rasterise buffered surface_terms
    surface_mask = features.rasterize(
        [(geom, 1) for geom in surface_gdf.geometry],
        out_shape=out_shape,
        transform=out_transform,
        fill=0,
        dtype="uint8"
    )

if mask_gdf.empty:
    print('No mask_terms found in tile, skipping mask_term masking.')
    mask_mask = np.zeros(out_shape, dtype="uint8")
else:
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

if np.all(np.isnan(chm_data)):
    print('No valid CHM data remains, skipping tile.')
    sys.exit()

###########################################
# Mask out buildings from OpenStreetMap
###########################################

print('Masking out buildings from OpenStreetMap')

osm_fp = f'{wd}/1_Reference_Data/11_OpenStreetMap/wales_osm_clean.gpkg'
osm_polygons = gpd.read_file(osm_fp,
                             layer='multipolygons',
                             bbox=tile_footprint)
osm_buildings = osm_polygons[osm_polygons['building'].notna()].copy()

if osm_buildings.empty:
    print('No buildings found in tile, skipping building masking.')
else:

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

    if np.all(np.isnan(chm_data)):
        print('No valid CHM data remains, skipping tile.')
        sys.exit()

###########################################
# Mask power cables
###########################################

print('Masking out power cables from CHM')

cables_fp = f'{wd}/1_Reference_Data/10_Power_Lines/OHL_buffered.gpkg'
cables = gpd.read_file(cables_fp)
cables = cables.set_crs(chm_crs) 

if cables.empty:
    print('No power cables found in tile, skipping cable masking.')
else:
    if cables.crs != chm_crs:
        cables = cables.to_crs(chm_crs)

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

    if np.all(np.isnan(chm_data)):
        print('No valid CHM data remains, skipping tile.')
        sys.exit()

###########################################
# Mask power structures from OpenStreetMap
###########################################

print('Masking out power structures from CHM')

power_structures_fp = f'{wd}/1_Reference_Data/11_OpenStreetMap/osm_power_structures_v2.gpkg'
power_structures = gpd.read_file(power_structures_fp, bbox=tile_footprint)

power_structures = power_structures[power_structures['power_source'] != 'solar'].copy()

if power_structures.empty:
    print('No power structures found in tile, skipping power structure masking.')
else:
        
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

    if np.all(np.isnan(chm_data)):
        print('No valid CHM data remains, skipping tile.')
        sys.exit()

###########################################
# Mask out coastline
###########################################

coastline_fp = f'{wd}/1_Reference_Data/11_OpenStreetMap/osm_coastline_buffered_50m.gpkg'
coastline = gpd.read_file(coastline_fp, bbox=tile_footprint)

# Check to see if gpkg is empty
if coastline.empty:
    print('No coastline features found in tile, skipping coastline masking.')
else:
    if coastline.crs != chm_crs:
        coastline = coastline.to_crs(chm_crs)

    coastline_mask = features.rasterize(
        [(geom, 1) for geom in coastline.geometry],
        out_shape=out_shape,
        transform=out_transform,
        fill=0,
        dtype="uint8"
    )

    coastline_mask = coastline_mask.astype(bool)
    chm_data[coastline_mask] = np.nan

    if np.all(np.isnan(chm_data)):
        print('No valid CHM data remains, skipping tile.')
        sys.exit()


# ####################################
# # Save to file - with NFI still
# ####################################

out_meta.update({
    'width' : out_shape[1],
    'height' : out_shape[0],
    'crs' : chm_crs,
    'transform' : out_transform,
})

chm_data = np.where(np.isnan(chm_data), -9999.0, chm_data)

output_path = f'{wd}/1_Reference_Data/0_VOM/with_nfi/{tile_of_interest}_VOM_with_NFI.tif'
with rasterio.open(output_path, "w", **out_meta) as dst:
    dst.write(chm_data, 1)

####################################
# NDVI Thresholding
####################################

chm_data_ndvi = chm_data.copy()

print('Resampling NDVI')

ndvi_fp = f'{wd}/1_Reference_Data/8_NDVI_and_NDWI/ndvi/{tile_of_interest}_ndvi_composite.tif'
ndvi_resampled_fp = f'{wd}/1_Reference_Data/8_NDVI_and_NDWI/ndvi/{tile_of_interest}_ndvi_composite_resampled.tif'

# Resample NDVI raster to match CHM raster (1m2)
with rasterio.open(ndvi_fp) as src:
    dst_meta = out_meta.copy()
    dst_meta.update({
        "dtype": src.meta["dtype"],  # keep original dtype
        "count": 1
    })

    with rasterio.open(ndvi_resampled_fp, 'w', **dst_meta) as dst:
        reproject(
            source=rasterio.band(src,1),
            destination=rasterio.band(dst,1),
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=out_transform,
            dst_crs=out_meta['crs'],
            dst_width=out_meta['width'],
            dst_height=out_meta['height'],
            resampling=Resampling.bilinear
        )

print('Applying NDVI threshold')

with rasterio.open(ndvi_resampled_fp) as src:
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

print('Saving VOM CHM with NDVI mask to file')

out_meta.update({
    'width' : out_shape[1],
    'height' : out_shape[0],
    'crs' : chm_crs,
    'transform' : out_transform,
})

chm_data_ndvi_masked = np.where(np.isnan(chm_data_ndvi_masked), -9999.0, chm_data_ndvi_masked)

output_path = f'{wd}/1_Reference_Data/0_VOM/with_nfi/{tile_of_interest}_VOM_with_NFI_NDVI.tif'
with rasterio.open(output_path, "w", **out_meta) as dst:
    dst.write(chm_data_ndvi_masked, 1)

print(f'Processing complete for tile {tile_of_interest}.')