###############################
# Eleanor R Downer
# 17 October 2025
# Clip Wales CHM to 10km tiles - using gdal_translate
###############################

import os
import subprocess
import geopandas as gpd

# Paths
CHM_PATH="Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/1_Wales_LiDAR/Wales_CHM_2020_22_32bit.tif"
SHP_PATH="Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/4_Wales_10km_Squares/Wales_10km_Squares.shp"
#OUTPUT_DIR="Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/1_Wales_LiDAR/CHM/"
OUTPUT_DIR="C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/data/"

print("Loading shapefile...")
gdf = gpd.read_file(SHP_PATH)

for idx, row in gdf.iterrows():

    tileref = row["tileref"]
    minx, miny, maxx, maxy = row.geometry.bounds
    output_path = f"{OUTPUT_DIR}{tileref}_CHM.tif"

    if os.path.exists(output_path):
        print(f"File {output_path} already exists. Skipping...")
        continue

    cmd = [
        r"C:\Program Files\QGIS 3.44.1\bin\gdal_translate.exe",
        "-of", "GTiff",
        "-projwin", str(minx), str(maxy), str(maxx), str(miny),
        "-a_nodata", "-9999",
        "-co", "COMPRESS=DEFLATE",
        "-co", "TILED=YES",
        "-co", "BIGTIFF=YES",
        CHM_PATH,
        output_path
    ]

    print(f"Processing {tileref}...")
    subprocess.run(cmd, check=True)
