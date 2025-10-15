###############################
# Eleanor R Downer
# 15 October 2025
# Clip Wales CHM to 10km tiles
###############################

import os
import subprocess
import geopandas as gpd
import glob

CHM_PATH="Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/1_Wales_LiDAR/Wales_CHM_2020_22_32bit.tif"
SHP_PATH="Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/4_Wales_10km_Squares/Wales_10km_Squares.shp"
OUTPUT_DIR="C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/data/clipped_chms/"

# Load the shapefile with 10km grid squares
gdf = gpd.read_file(SHP_PATH)

# Loop through each grid square
for idx, row in gdf.iterrows():
    tileref = row['tileref']  # Adjust based on actual column name
    print(f"Processing square: {tileref}")

    # Save the geometry to a temporary shapefile
    temp_shp = f"{OUTPUT_DIR}{tileref}.shp"
    gdf_temp = gpd.GeoDataFrame([row], crs=gdf.crs)
    gdf_temp.to_file(temp_shp)

    output_path = f"{OUTPUT_DIR}{tileref}_CHM.tif"

    minx, miny, maxx, maxy = row['geometry'].bounds

    cmd = [
        r"C:\Program Files\QGIS 3.44.1\bin\gdalwarp.exe",
        "-overwrite",
        "-t_srs", "EPSG:27700",
        "-te", str(minx), str(miny), str(maxx), str(maxy),
        "-cutline", temp_shp,
        "-crop_to_cutline",
        "-dstnodata", "-9999.0",
        "-ot", "Float32",
        "-multi",
        "-co", "COMPRESS=DEFLATE",
        "-co", "TILED=YES",
        CHM_PATH,
        output_path
    ]

    subprocess.run(cmd, check=True)

    # Remove all files for a shapefile (shp, shx, dbf, prj, cpg, etc.)
    shapefile_base = temp_shp.rsplit('.', 1)[0]  # remove extension
    for f in glob.glob(f"{shapefile_base}.*"):
        try:
            os.remove(f)
        except Exception as e:
            print(f"Error removing {f}: {e}")

    print(f"Saved clipped CHM for square {tileref} to {output_path}")