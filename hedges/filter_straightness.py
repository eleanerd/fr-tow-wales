#########################################
# Filter hedge polygons by straightness
#########################################

import geopandas as gpd
import glob
import os

# --- Inputs ---

print("Loading hedge files...")
hedge_dir = 'C:\\Users\\eleanor.downer\\OneDrive - Forest Research\\Documents\\TOW_Wales\\hedge_processing\\hedges'
hedge_files = glob.glob(f'{hedge_dir}\\*_hedge*.gpkg')

print("Loading built-up areas...")
built_up_areas_gpkg = 'Y:\\Forest Inventory\\0700_NonCore_Funded\\0726_TOW_Wales\\04_Spatial Analysis\\1_Reference_Data\\13_OS_Built_Up_Areas\\OS_Open_Built_Up_Areas_GeoPackage\\os_open_built_up_areas_wales.gpkg'
built_up_areas = gpd.read_file(built_up_areas_gpkg, layer='os_open_built_up_areas')
built_up_areas["geometry"] = built_up_areas.buffer(-10)
built_up_areas_union = built_up_areas.union_all()

# --- Process each tile ---

for hedge_file in hedge_files:

    tile = hedge_file.split('\\')[-1].split('_')[0]
    print('Processing tile', tile)

    output_path = f'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/4_Processing/VOM_Processing/Hedges/{tile}_hedges.gpkg'
    output_layer_name = f'{tile}_hedges'
    if os.path.exists(output_path):
        print(f'{tile} already processed. Skipping.')
        continue

    hedges = gpd.read_file(hedge_file)

    # --- Filter by straightness ---

    print("Filtering hedges by straightness...")

    hedges = hedges[~((hedges['straightness'] < 0.6) & (hedges['mbcp'] > 10))]
    hedges['geometry'] = hedges['geometry'].buffer(0)

    # --- Remove hedges in built-up areas --
    # Only run for hegdes_v4 files (which contain built-up area hedges)
    if 'hedges_v4' in hedge_file:
        print("Removing hedges in built-up areas...")
        hedges = hedges.overlay(built_up_areas, how="difference")
        hedges = hedges[~hedges.geometry.is_empty]

    # --- Dissolve touching/overlapping polygons ---

    print("Dissolving touching/overlapping polygons...")

    dissolved = gpd.GeoDataFrame(
        geometry=[hedges.geometry.union_all()],
        crs=hedges.crs
    ).explode(index_parts=False).reset_index(drop=True)

    dissolved['area_m2'] = dissolved.geometry.area.round(2)

    # --- Remove small polygons ---
    
    dissolved = dissolved[dissolved['area_m2'] >= 20]

    # --- Save updated ---

    print("Saving updated hedge file...")
    
    dissolved.to_file(
        output_path,
        driver="GPKG",
        layer=output_layer_name,
        overwrite=True
    )

    print(f"Saved updated hedge file: {output_path}")