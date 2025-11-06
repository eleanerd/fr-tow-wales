#########################################
# Filter hedge polygons by straightness
#########################################

import geopandas as gpd
import glob

hedge_dir = 'C:\\Users\\eleanor.downer\\OneDrive - Forest Research\\Documents\\TOW_Wales\\hedge_processing\\hedges'
hedge_files = glob.glob(f'{hedge_dir}\\*_hedges_v3.gpkg')

for hedge_file in hedge_files:

    tile = hedge_file.split('\\')[-1].split('_')[0]
    print('Processing tile', tile)

    hedges = gpd.read_file(hedge_file)

    hedges = hedges[~((hedges['straightness'] < 0.6) & (hedges['mbcp'] > 10))]

    # Save updated
    hedges.to_file(f'{hedge_dir}\\{tile}_hedges_filt.gpkg', driver="GPKG", overwrite=True)