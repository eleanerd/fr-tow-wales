###################################
# Convert cable lines into polygons
###################################

import geopandas as gpd

lines = gpd.read_file('Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/10_Power_Lines/OHL.shp')

if lines.crs != 'EPSG:27700':
    lines = lines.to_crs('EPSG:27700') 

# Buffer the lines
buffered = lines.buffer(12)

# Turn buffered geometries into polygons (buffer returns polygons already)
polygons = gpd.GeoDataFrame(geometry=buffered, crs=lines.crs)

# Optionally dissolve overlapping buffers into single polygons
polygons = polygons.dissolve()

polygons.to_file('Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/1_Reference_Data/10_Power_Lines/OHL_buffered.gpkg')
