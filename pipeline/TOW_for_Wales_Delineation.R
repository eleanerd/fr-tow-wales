##################
# Freddie Hunter
# 03 Dec 2025
# Wales TOW V3
#################

# Assumptions
# This code assumes that the input is a VOM with a minimum pixel cluster size of 5pixel.
# The code also assumes that the VOM has been filtered with a minimum of 50cm 
# The code also asumes that the VOM has NA values not zero values

rm(list = ls())
#############################################################
# Input files
# There are many parameters imbedded in the code. 
wd <- "//forestresearch.gov.uk/shares/IFOS/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/3_Test_Squares/test/"
wd <- "C:/Users/freddie.hunter/Downloads/New folder/"
raz_files <- 'SS79_CHM_1.3m_masked_NFI_OSSM_Cables.tif'
NFI_Layer = "C:/Users/freddie.hunter/Downloads/NFI_Wales_2024_WoodlandOnly.gpkg"
setwd(wd)

#############################################################
#Libraries
library(terra)
library(dplyr)
library(sf)
library(raster)
library(glue)
library(deldir)
library(dplyr)
library(lidR)
library(lwgeom)
library(nngeo)
library(smoothr)
library(purrr)
library(foreach)
library(doParallel)

#################################################
#Processing

# Load CHM raster
chm <- terra::rast(raz_files)

# split into tiles
tiles <- terra::makeTiles(chm, c(2000,2000))  
flist <- list.files(wd, pattern = "^tile_[0-9]{1,3}\\.tif$", full.names = T)

###############################
### Run first part in parallel
###############################

numCores <- parallel::detectCores()

# Create a cluster
cl <- makeCluster(numCores - 15)   # keep 1 core free
registerDoParallel(cl)

# run method per sw
new_data <- foreach(
  i = 1:length(flist),
  .packages = c(
    "sf",
    "dplyr",
    "units",
    "nngeo",
    "terra",
    "raster",
    "glue",
    "deldir",
    "lidR",
    "lwgeom",
    "smoothr",
    "purrr"
  ),
  .export = c("flist", "NFI_Layer", "wd")
) %dopar% {
  
  raz_l <- flist[i]
  
  if(file.exists(glue("{wd}/filt_max_class_{raz_l}.tif"))){
    print('processed')
    next()
  }
  
  chm <- rast(raz_l)
  
  if (is.nan(minmax(chm)[2])){
    print('no data exists')
    return(NULL)
  } else { 
    print('data exists')
  }
  
  raz_l <- basename(raz_l) %>% stringr::str_remove('.tif')
  # Step 1: 5x5 local maximum filter
  # Assume cell size = 1m → 6x6 window = 6 cells
  local_max <- terra::focal(chm,
                            w = matrix(1, 5, 5),
                            fun = max,
                            na.rm = T)
  
  # set all NA values in original data back to NA to stop expansion by window
  local_max[is.na(chm)] <- NA
  
  # Step 2: Initial height-based classification
  # Reclassification matrix: start, end, new class
  rcl <- matrix(
    c(0.49, 10, 3, 10, 15, 10, 15, 20, 15, 20, 30, 20, 30, Inf, 30),
    ncol = 3,
    byrow = TRUE
  )
  
  local_max <- as(local_max, "SpatRaster")
  class_raster <- classify(local_max, rcl)
  
  # second round of filter for 5 clus with na set to 1
  class_raster[is.na(class_raster)] <- 1
  r_sieved_1 <- sieve(class_raster, threshold = 5, # minimum number of connected cells
                      directions = 4)   # 4-connected (use 8 for diagonals))
  
  # second round of filter for 25clus size with 1 set na
  r_sieved_1[r_sieved_1 == 1] <- NA
  filt_max_class <- sieve(r_sieved_1, threshold = 25, # minimum number of connected cells
                          directions = 4)   # 4-connected (use 8 for diagonals))
  #writeRaster(filt_max_class, glue('{wd}/filt_max_class_{raz_files}'), overwrite=TRUE)
  
  # Tree detection
  # get tree detection function
  x <- c(3, 5, 10, 20)
  y <- c(2, 2.6, 3.2, 5)
  
  model <- lm(log(y) ~ log(x))
  a <- exp(coef(model)[1])
  b <- coef(model)[2]
  
  # Define function
  f <- function(x) a * x^b
  
  # Detect treetops using adaptive window size
  # clean up CHM
  fw <- matrix(nrow = 3, ncol = 3)
  fw[is.na(fw)] <- 1
  
  chm[is.na(chm)] <- 0
  chm_s <- terra::focal(
    x = chm,
    w = fw,
    fun = mean,
    na.rm = T
  )
  # set na values in filtt_max_class as na in chm_s
  chm_s[is.na(filt_max_class)] <- NA
  chm_s[chm_s < 0.5] <- NA
  
  ttops <- locate_trees(chm_s, lmf(ws = f, shape = 'circular', hmin = 3))
  
  # change ttops crs to NFI crs #####
  layer <- st_layers(NFI_Layer)
  meta <- st_read(NFI_Layer, query = glue("SELECT * from {layer$name} LIMIT 1 OFFSET 1"))
  
  # change tops crs
  ttops <- st_transform(ttops, crs = st_crs(meta))
  st_write(ttops, glue('{wd}/ttops_{raz_l}.gpkg'), append = F)
  
  # Get crown extent for each tree.
  # The crown extent height values can go as low as possible, this avoids cropping crowns at the 2.5m height, attributing total crown area to trees with height >= 3m.
  # only pixels below 3m which do not directly correspond to a 3m or above tree's canopy.
  crowns <- dalponte2016(chm = chm_s, treetops = ttops)()
  # but we do want to ensure all VOM pixels equal and above 3 are included in the output.
  crowns[chm_s >= 3] <- 2
  
  # Make crowns binary
  # avoid using 0 as gdal treats it as no data.
  crowns[!is.na(crowns)] <- 2
  crowns[is.na(crowns)] <- 1
  
  # write crowns to file
  #writeRaster(crowns, glue('{wd}/crowns_{raz_files}'), overwrite = T)
  
  # filter out any pixel gaps smaller than 5m.
  crowns_fllt <- sieve(as(crowns, 'SpatRaster'),
                       threshold = 5,
                       directions = 4)   # 4-connected (use 8 for diagonals))
  #writeRaster(crowns_fllt, glue("{wd}/crowns_clus_5_{raz_files}"), overwrite = T)
  
  # set areas that are not part of 3m tree crowns to zero in binary layer
  filt_max_class[crowns_fllt == 1] <- NA
  # Keep version with bins
  filt_max_class_bins <- filt_max_class
  
  writeRaster(filt_max_class_bins,
              glue("{wd}/filt_max_class_bins_{raz_l}.tif"),
              overwrite = T)
  
  # Now define woodland as binary
  filt_max_class[!is.na(filt_max_class)] <- 1
  
  # empty ram
  vars_all = ls()
  vars_all <- vars_all[!vars_all %in% c(
    "filt_max_class",
    "raz_l",
    "wd", 
    "NFI_Layer", 
    "flist"
  )]  
  rm(list = vars_all)
  
  #writeRaster(filt_max_class,
  #            glue("{wd}/filt_max_class_{raz_l}.tif"),
  #            overwrite = T)
  
  # set crs
  layer <- st_layers(NFI_Layer)
  meta <- st_read(NFI_Layer, query = glue("SELECT * from {layer$name} LIMIT 1 OFFSET 1"))
  
  # delineate canopy
  canopy_area <- st_as_sf(
    stars::st_as_stars(filt_max_class),
    point = FALSE,
    merge = TRUE,
    connect8 = TRUE
  ) %>% st_union() %>%
    st_cast("POLYGON") %>%
    st_transform(st_crs(meta))
  
  # get bbx for NFI crop
  dat_bbox <- st_bbox(canopy_area, crs = st_crs(canopy_area)) %>% st_as_sfc() %>% st_transform(st_crs(meta))
  
  # CHANGE TO LATEST NFI!!!
  bbx_cutout <-
    st_read(NFI_Layer) %>%
    st_intersection(dat_bbox) %>%
    st_as_sf() 
  bbx_cutout <- bbx_cutout %>% st_collection_extract("POLYGON") %>% st_union() %>% st_cast('POLYGON') %>% st_as_sf()
  
  if (isTRUE(dim(bbx_cutout)[0] != 0)){
    bbx_cutout$ROWID <- c(1:dim(bbx_cutout)[1])
  } else {
    bbx_cutout$ROWID <- NULL
  }
  
  # cut out nfi from TOW
  canopy_area_nfirm <- st_difference(canopy_area, bbx_cutout)
  
  # get nfi 10m buffer
  NFI_10BUF <- bbx_cutout %>% st_buffer(10) 
  
  # set NFI OHC
  NFI_OHC <- st_intersection(canopy_area_nfirm, NFI_10BUF) %>% st_collection_extract("POLYGON") %>% st_union() %>% st_cast('POLYGON')
  NFI_OHC_F <- st_difference(NFI_OHC, bbx_cutout)  %>% st_union() %>% st_cast('POLYGON') %>% st_as_sf()
  
  if (isTRUE(dim(NFI_OHC_F)[0] != 0)){
    NFI_OHC_F$id <- c(1:dim(NFI_OHC_F)[1])
  } else {
    NFI_OHC_F$id <- NULL
  }
  
  # Keep if touching or overlapping original NFI
  OL <- st_overlaps(NFI_OHC_F, bbx_cutout) %>% as.data.frame()
  TC <- st_touches(NFI_OHC_F, bbx_cutout) %>% as.data.frame()
  ids_keep <- rbind(OL, TC) %>% dplyr::select(row.id)
  if (isTRUE(dim(ids_keep)[0] != 0)){
    NFI_OHC_F <- NFI_OHC_F %>% filter(id %in% unlist(ids_keep))
  }
  NFI_OHC_F$area <- st_area(NFI_OHC_F) %>% units::drop_units()
  NFI_OHC_F <- NFI_OHC_F[NFI_OHC_F$area >= 5,]
  if (isTRUE(dim(ids_keep)[0] != 0)){
    NFI_OHC_F$Woodland_Type <- 'NFI OHC'
    
  }
  
  # Remove OHC from TOW
  NFI_OHC_F <- st_transform(NFI_OHC_F, st_crs(canopy_area))
  canopy_area <- st_make_valid(canopy_area)
  NFI_OHC_F   <- NFI_OHC_F %>% st_make_valid()
  idx <- st_intersects(canopy_area, NFI_OHC_F) 
  
  result_geom <- map2(
    st_geometry(canopy_area),
    idx,
    \(g, id) {
      g <- st_sfc(g, crs = st_crs(canopy_area))  
      if (length(id) == 0) {
        return(g)   # keep original geometry
      }
      d <- st_difference(g, st_union(NFI_OHC_F[id, ]))
      return(d)
    }
  )
  
  # 4. Bind results back with correct CRS
  canopy_area_tow <- do.call(c, result_geom) %>% st_union() %>% st_cast('POLYGON') %>% st_as_sf()
  
  # Write final NFI OHC area per tile
  st_write(NFI_OHC_F, glue('{wd}/NFI_OHC_F_{raz_l}.gpkg'), append = F)
  
  # Write final TOW area per tile
  st_write(canopy_area_tow, glue('{wd}/canopy_area_tow_{raz_l}.gpkg'), append = F)
  
}

stopCluster(cl)
# END part 1
########################################################################

# Start Part 2

# set crs
layer <- st_layers(NFI_Layer)
meta <- st_read(NFI_Layer, query = glue("SELECT * from {layer$name} LIMIT 1 OFFSET 1"))

# merge tow canopy area
cat <- list.files(wd, pattern = 'canopy_area_tow_tile_[0-9]{1,3}\\.gpkg$', full.names = T)
cat_list <- lapply(cat, st_read, quiet = TRUE)
cat_merged <- do.call(rbind, cat_list) %>% st_union() %>% st_cast('POLYGON') %>% st_as_sf()
cat_merged$id <- 1:dim(cat_merged)[1]
# Merge ttops
ttops <- list.files(wd, pattern = 'ttops_tile_[0-9]{1,3}\\.gpkg$', full.names = T)
ttops_list <- lapply(ttops, st_read, quiet = TRUE)
ttops_merged <- do.call(rbind, ttops_list) 

# Merge NFI_OHC_F
nfi_ohc <- list.files(wd, pattern = 'NFI_OHC_F_tile_[0-9]{1,3}\\.gpkg$', full.names = T)
nfi_ohc_list <- lapply(nfi_ohc, st_read, quiet = TRUE)
nfi_ohc_merged <- do.call(rbind, nfi_ohc_list) %>% st_union() %>% st_cast('POLYGON')  %>% st_as_sf()

# add trct to NFI OHC merged
nfi_ohc_merged$id <- c(1:dim(nfi_ohc_merged)[1])
NFI_OHC_F_trct <- st_intersection(nfi_ohc_merged,ttops_merged) %>% st_drop_geometry() %>% group_by(id) %>% summarise(trct = n())
nfi_ohc_merged <- base::merge(nfi_ohc_merged, NFI_OHC_F_trct, by = 'id', all = TRUE )
nfi_ohc_merged$trct[is.na(nfi_ohc_merged$trct)] <- 0

# get atts for WT
ttops_merged <- ttops_merged %>% st_transform(st_crs(meta))
cat_merged$area <- st_area(cat_merged) %>% units::drop_units()
cat_merged <- cat_merged[cat_merged$area >=5,]
cat_merged$mbc <-  st_area(lwgeom::st_minimum_bounding_circle(cat_merged$x)) %>% units::drop_units()
trct_dat <- st_intersection(cat_merged,ttops_merged) %>% st_drop_geometry() %>% group_by(id) %>% summarise(trct = n())
cat_merged <- base::merge(cat_merged, trct_dat, by = 'id', all = TRUE )
cat_merged$trct[is.na(cat_merged$trct)] <- 0
cat_merged$mbc_p <- cat_merged$area/cat_merged$mbc * 100

# set Woodland types
wt_canopy <-
  cat_merged %>% mutate(Woodland_Type = case_when((area >= 5 & area <= 350 & trct == 0 ~ 'Lone Tree'),
                                                  (area >= 5 & area <= 350 & mbc_p >= 55 ~ 'Lone Tree'),
                                                  (area >= 5 & trct == 1 ~ 'Lone Tree'),
                                                  (area >= 5 & area <= 350 & trct == 2 & mbc_p >= 50 ~ 'Lone Tree'),
                                                  (area >= 5 & area <= 1000 & trct > 1 ~ 'Group of Trees'),
                                                  (area > 1000 ~ 'Small Woodland')))

# Set any lone trees touching NFI OHC to Group of trees. 
NFI_OHC_F_LT <- wt_canopy %>% filter(Woodland_Type == 'Lone Tree') %>% st_touches(nfi_ohc_merged) %>% as.data.frame()
wt_canopy$Woodland_Type[wt_canopy$Woodland_Type == 'Lone Tree'][NFI_OHC_F_LT$row.id] <-  'Group of Trees'

# Add back in NFI_OHC_F
nfi_ohc_merged$mbc <- NA
nfi_ohc_merged$mbc_p <- NA
nfi_ohc_merged$Woodland_Type <- 'NFI OHC'
nfi_ohc_merged$area <- st_area(nfi_ohc_merged) %>% units::drop_units()
wt_canopy_f <- rbind(wt_canopy, nfi_ohc_merged)

# clear out the ram
#rm(list = c('ttops_merged', 'nfi_ohc_merged', 'cat_merged'))

# add back in height bins
# Merge rasters and polygonize
filt_max_class_bins_list  <- list.files(wd, pattern = 'filt_max_class_bins_tile_[0-9]{1,3}\\.tif$', full.names = T)
filt_max_class_bins_rlist <- lapply(filt_max_class_bins_list, rast)
filt_max_class_bins_m <- do.call(merge, filt_max_class_bins_rlist) 

# set areas that are not part of 3m tree crowns to zero. 
filt_max_class_bins_p <- st_as_sf(
  stars::st_as_stars(filt_max_class_bins_m),
  point = FALSE,
  connect8 = TRUE,
  merge = TRUE) %>%
  st_cast("POLYGON") %>%
  st_transform(st_crs(meta))

# now burn back in the bins to the main delineation
# empty ram
vars_all = ls()
vars_all <- vars_all[!vars_all %in% c(
  "wd",
  "raz_l",
  "wt_canopy_f", 
  "filt_max_class_bins_p", 
  "raz_files",
  "vars_all",
  "NFI_Layer"
)]  
rm(list = vars_all)

all_canopy_bins <- st_intersection(wt_canopy_f, filt_max_class_bins_p) %>% st_collection_extract("POLYGON") %>% st_cast('POLYGON')
all_canopy_bins$Bin_Area <- st_area(all_canopy_bins) %>% units::drop_units()

# add height raster stats
chm <- rast(glue("{wd}/{raz_files}"))
ht_dat <- terra::extract(
  chm,
  all_canopy_bins,
  fun = function(x) {
    c(
      min  = min(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      max  = max(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE)
    )
  }
)

all_canopy_bins$min_ht <- ht_dat[,2]
all_canopy_bins$mean_ht <- ht_dat[,3]
all_canopy_bins$max_ht <- ht_dat[,4]
all_canopy_bins$sd_ht <- ht_dat[,5]
names(all_canopy_bins) <- c('TOW_ID', 'Area_M', 'MBC', 'TRCT', 'MBC_P', 'Woodland_Type', 'HT_BIN', 'x', 'HT_BIN_AREA_M', 'HT_MIN', 'HT_MEAN', 'HT_MAX', 'HT_SD')
st_write(all_canopy_bins, glue('{wd}/TOW_Wales_Canopy_{raz_files}.gpkg'), append = F)

# end