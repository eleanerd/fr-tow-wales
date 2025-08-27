####################################
# NDVI and NDWI for TOW Wales
####################################

import os
import rasterio
import glob
import numpy as np

planet_dir = 'Y:\\Forest Inventory\\0700_NonCore_Funded\\0726_TOW_Wales\\04_Spatial Analysis\\1_Reference_Data\\7_Planet_Basemaps\\'
planet_tiles = glob.glob(planet_dir + '\\*\\L15-1002E-1368N.tif')

output_dir = 'Y:\\Forest Inventory\\0700_NonCore_Funded\\0726_TOW_Wales\\04_Spatial Analysis\\1_Reference_Data\\8_NDVI_and_NDWI'

for planet_tile in planet_tiles:

    tile = planet_tile.split('\\')[-1].split('.')[0]
    print(f'Processing {tile}...')

    month = planet_tile.split('\\')[-2].split('_2023')[0]
    print(month)

    if os.path.exists(f'{output_dir}//NDVI//{tile}_{month}_NDVI.tif') and os.path.exists(f'{output_dir}//NDWI//{tile}_{month}_NDWI.tif'):
        print(f'NDVI and NDWI for {tile} in {month} already exist. Skipping...')
        continue
    
    with rasterio.open(planet_tile) as src:

        meta = src.meta

        # Read the red, green, and NIR bands
        red = src.read(6).astype(float)
        green = src.read(3).astype(float)
        nir = src.read(8).astype(float)

        # Calculate NDVI
        ndvi = (nir - red) / (nir + red)
        
        # Calculate NDWI
        ndwi = (green - nir) / (green + nir)

    meta.update({
        'dtype': 'float32',
        'count': 1
    })

    # Save the NDVI and NDWI results
    with rasterio.open(f'{output_dir}//NDVI//{tile}_{month}_NDVI.tif', 'w', **meta) as dst:
        dst.write(ndvi.astype(np.float32), 1)

    with rasterio.open(f'{output_dir}//NDWI//{tile}_{month}_NDWI.tif', 'w', **meta) as dst:
        dst.write(ndwi.astype(np.float32), 1)

print('Calculate yearly and seasonal means of NDWI')

spring_months = ['March', 'April', 'May']
summer_months = ['June', 'July', 'August']
winter_months = ['December', 'January', 'Feb'] 
autumn_months = ['September', 'October', 'November']

ndwi_files = glob.glob(f'{output_dir}//NDWI//L15-1002E-1368N_*.tif')
ndwi_stack = np.array([rasterio.open(f).read(1) for f in ndwi_files])
ndwi_stack[ndwi_stack == -9999] = np.nan

ndwi_mean = np.nanmean(ndwi_stack, axis=0)
ndwi_spring = np.nanmean(np.array([rasterio.open(f).read(1) for f in ndwi_files if any(m in f for m in spring_months)]), axis=0)
ndwi_summer = np.nanmean(np.array([rasterio.open(f).read(1) for f in ndwi_files if any(m in f for m in summer_months)]), axis=0)
ndwi_winter = np.nanmean(np.array([rasterio.open(f).read(1) for f in ndwi_files if any(m in f for m in winter_months)]), axis=0)
ndwi_autumn = np.nanmean(np.array([rasterio.open(f).read(1) for f in ndwi_files if any(m in f for m in autumn_months)]), axis=0)

with rasterio.open(ndwi_files[0]) as src:
    meta = src.meta

meta.update({
    'dtype': 'float32',
    'count': 5
})

with rasterio.open(f'{output_dir}//NDWI//L15-1002E-1368N_2023_NDWI_yearly_and_seasonal_means.tif', 'w', **meta) as dst:
    dst.write(ndwi_mean.astype(np.float32), 1)
    dst.set_band_description(1, 'Yearly Mean NDWI')
    dst.write(ndwi_spring.astype(np.float32), 2)
    dst.set_band_description(2, 'Spring Mean NDWI')
    dst.write(ndwi_summer.astype(np.float32), 3)
    dst.set_band_description(3, 'Summer Mean NDWI')
    dst.write(ndwi_winter.astype(np.float32), 4)
    dst.set_band_description(4, 'Winter Mean NDWI')
    dst.write(ndwi_autumn.astype(np.float32), 5)
    dst.set_band_description(5, 'Autumn Mean NDWI')


ndwi_median = np.nanmedian(ndwi_stack, axis=0)
ndwi_spring_median = np.nanmedian(np.array([rasterio.open(f).read(1) for f in ndwi_files if any(m in f for m in spring_months)]), axis=0)

meta.update({
    'count': 2
})

with rasterio.open(f'{output_dir}//NDWI//L15-1002E-1368N_2023_NDWI_median.tif', 'w', **meta) as dst:
    dst.write(ndwi_median.astype(np.float32), 1)
    dst.set_band_description(1, 'Yearly Median NDWI')
    dst.write(ndwi_spring_median.astype(np.float32), 2)
    dst.set_band_description(2, 'Spring Median NDWI')

ndvi_files = glob.glob(f'{output_dir}//NDVI//L15-1002E-1368N_*.tif')
ndvi_stack = np.array([rasterio.open(f).read(1) for f in ndvi_files])
ndvi_stack[ndvi_stack == -9999] = np.nan

ndvi_median = np.nanmedian(ndvi_stack, axis=0)
ndvi_spring_median = np.nanmedian(np.array([rasterio.open(f).read(1) for f in ndvi_files if any(m in f for m in spring_months)]), axis=0)

with rasterio.open(f'{output_dir}//NDVI//L15-1002E-1368N_2023_NDVI_median.tif', 'w', **meta) as dst:
    dst.write(ndvi_median.astype(np.float32), 1)
    dst.set_band_description(1, 'Yearly Median NDVI')
    dst.write(ndvi_spring_median.astype(np.float32), 2)
    dst.set_band_description(2, 'Spring Median NDVI')