
############################
# Merge rasters with GDAL
############################

import subprocess
import glob

# Path to gdal_merge.py
gdal_merge_path = r"C:\Program Files\QGIS 3.44.1\apps\Python312\Scripts\gdal_merge.py"

# Output filename
#output_file = r"Y:\Forest Inventory\0700_NonCore_Funded\0726_TOW_Wales\04_Spatial Analysis\3_Test_Square\VOM\SS79_VOM_NFI_hedges_extracted.tif"
output_file = r"Y:\Forest Inventory\0700_NonCore_Funded\0726_TOW_Wales\04_Spatial Analysis\3_Test_Square\Delineation\SS79_class_raster.tif"

# Input file pattern
#input_pattern = r"Y:\Forest Inventory\0700_NonCore_Funded\0726_TOW_Wales\04_Spatial Analysis\3_Test_Square\Hedges\CHMs\section_*.tif"
#input_files = glob.glob(input_pattern) # Create list of input files

input_pattern = r"Y:\Forest Inventory\0700_NonCore_Funded\0726_TOW_Wales\04_Spatial Analysis\3_Test_Square\Delineation\SS79_*.tif"
input_files = glob.glob(input_pattern) # Create list of input files

cmd = [
    r"C:\Program Files\QGIS 3.44.1\apps\Python312\python.exe",
    gdal_merge_path,
    "-o", output_file,
    "-of", "GTiff",
    "-a_nodata", "-9999" # no data value
] + input_files

subprocess.run(cmd, check=True)