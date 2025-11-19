#########################################
# Filter hedge polygons by straightness
#########################################

import geopandas as gpd
import glob
import os

# --- Inputs ---

hedge_dir = 'C:\\Users\\eleanor.downer\\OneDrive - Forest Research\\Documents\\TOW_Wales\\hedge_processing\\hedges'
hedge_files = glob.glob(f'{hedge_dir}\\*_hedge*.gpkg')

built_up_areas_gpkg = gpd.read_file('Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/13_OS_Built_Up_Areas/OS_Open_Built_Up_Areas_GeoPackage/os_open_built_up_areas.gpkg')
built_up_areas_layer = 'os_open_built_up_areas'
built_up_areas = gpd.read_file(built_up_areas_gpkg, layer=built_up_areas_layer)
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

    hedges = hedges[~((hedges['straightness'] < 0.6) & (hedges['mbcp'] > 10))]
    hedges['geometry'] = hedges['geometry'].buffer(0)

    # --- Remove hedges in built-up areas ---

    hedges["geometry"] = hedges.geometry.apply(lambda g: g.difference(built_up_areas_union))
    hedges = hedges[~hedges.geometry.is_empty]

    # --- Dissolve touching/overlapping polygons ---

    dissolved = gpd.GeoDataFrame(
        geometry=[hedges.geometry.union_all()],
        crs=hedges.crs
    ).explode(index_parts=False).reset_index(drop=True)

    dissolved['area_m2'] = dissolved.geometry.area.round(2)

    # --- Save updated ---
    
    dissolved.to_file(
        output_path,
        driver="GPKG",
        layer=output_layer_name,
        overwrite=True
    )

    print(f"Saved updated hedge file: {output_path}")