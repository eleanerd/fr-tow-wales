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

###############################
# Settings / Inputs 
###############################

#tile_of_interest = sys.argv[1] 
tile_of_interest = 'SH57'  # For testing
wd = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis'

output_path = f'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/4_Processing/CHM_FILT/{tile_of_interest}_CHM_FILT.tif'
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

print('Removing pixels below 1.3m from CHM')

chm_fp = f'C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/data/{tile_of_interest}_CHM.tif'

with rasterio.open(chm_fp) as src:

    #window = from_bounds(*tile_footprint.total_bounds, transform=src.transform)
    chm_data = src.read(1).astype('float32')

    out_transform = src.transform
    out_meta = src.meta.copy()
    chm_crs = src.crs
    
    out_shape = (int(chm_data.shape[0]), int(chm_data.shape[1]))

    chm_data = np.where(chm_data < 1.3, np.nan, chm_data)

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
    
    water_gdf["geometry"] = water_gdf.buffer(-6)
    water_gdf["geometry"] = water_gdf.buffer(4)
    water_gdf = water_gdf[~water_gdf.is_empty]

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

####################################
# Remove Power Cables (OpenStreetMap)
# Missed in Pt 1
####################################

power_lines_fp = f'{wd}/1_Reference_Data/11_OpenStreetMap/osm_power_lines.gpkg'
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
        out_shape=out_shape, #out_shape=out_shape,
        transform=out_transform, #transform=out_transform,
        fill=0,
        dtype="uint8"
    )

    # Mask power lines out of VOM
    power_lines_mask = power_lines_mask.astype(bool)
    chm_data[power_lines_mask] = np.nan


###########################################
# Mask out coastline
###########################################

print('Masking out coastline')

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

################################################
# Mask out large, flat objects
# (flat buildings/structures and gorse/shrub)
################################################

print('Masking out flat top objects')

chm_var_path = f"{wd}/1_Reference_Data/0_VOM/with_nfi/height_var/{tile_of_interest}_height_var.tif"

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

ndvi_fp = f'{wd}/1_Reference_Data/8_NDVI_and_NDWI/ndvi/{tile_of_interest}_ndvi_composite_resampled.tif'
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

out_meta.update({
    'width' : (chm_data.shape)[1],
    'height' : (chm_data.shape)[0],
    'crs' : chm_crs,
    'transform' : out_transform,
    "nodata": -9999.0,
    "dtype": "float32",
    "compress": "LZW",
    "tiled": True
})

chm_data_ndvi_masked = np.where(
    np.isnan(chm_data_ndvi_masked),
    -9999.0,
    chm_data_ndvi_masked
)

with rasterio.open(output_path, "w", **out_meta) as dst:
    dst.write(chm_data_ndvi_masked, 1)