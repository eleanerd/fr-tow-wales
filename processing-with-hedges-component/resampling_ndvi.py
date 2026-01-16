###########################################
# Resample NDVI raster to match CHM raster
###########################################

import rasterio
from rasterio.warp import reproject, Resampling

ndvi_raster = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/3_Test_Square/NDVI/SS79_ndvi_composite.tif'
chm_raster = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/3_Test_Square/VOM/SS79_VOM_no_NFI.tif'
out_raster = 'Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/3_Test_Square/NDVI/SS79_ndvi_composite_resampled.tif'

with rasterio.open(chm_raster) as ref:
    ref_meta = ref.meta.copy()
    ref_transform = ref.transform
    ref_crs = ref.crs
    ref_shape = (ref.height, ref.width)

with rasterio.open(ndvi_raster) as src:
    dst_meta = ref_meta.copy()
    dst_meta.update({
        "dtype": src.meta["dtype"],  # keep original dtype
        "count": 1
    })

    with rasterio.open(out_raster, 'w', **dst_meta) as dst:
        reproject(
            source=rasterio.band(src,1),
            destination=rasterio.band(dst,1),
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=ref_transform,
            dst_crs=ref_crs,
            dst_width=ref_shape[1],
            dst_height=ref_shape[0],
            resampling=Resampling.bilinear
        )