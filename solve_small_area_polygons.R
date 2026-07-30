
library(sf)
library(dplyr)
library(purrr)
library(nngeo)   # for st_nn (nearest neighbor / touching polygons)
library(glue)
library(doParallel)
library(foreach)

merge_small_polygons <- function(tile_sf, min_area = 5) {
  
  repeat {
    
    # Find small polygons
    small_polys <- tile_sf %>% filter(Bin_Area < min_area)
    if (nrow(small_polys) == 0) break
    
    message(nrow(small_polys))
    
    merged_any <- FALSE
    
    for (i in seq_len(nrow(small_polys))) {
      
      small_poly <- small_polys[i, ]
      
      # First filter by TOW_ID and size
      candidates <- tile_sf %>%
        filter(TOW_ID == small_poly$TOW_ID & Bin_Area >= min_area)
      
      if (nrow(candidates) == 0) next
      
      # Then do spatial predicate only on candidates
      neighbors_idx <- st_touches(small_poly, candidates)[[1]]
      
      if (length(neighbors_idx) == 0) next
      
      neighbors <- candidates[neighbors_idx, ]
      
      if (nrow(neighbors) == 0) next
      
      # Choose neighbor by closest Height_Class, then smaller Bin_Area
      neighbors <- neighbors %>% mutate(Height_diff = abs(Height_Class - small_poly$Height_Class))
      target <- neighbors %>% arrange(Height_diff, Bin_Area) %>% slice(1)
      
      # Merge geometries
      merged_geom <- st_union(st_geometry(small_poly), st_geometry(target))
      
      # Update target polygon
      target_row <- which(st_geometry(tile_sf) == st_geometry(target))
      tile_sf$geom[target_row] <- merged_geom
      tile_sf$Bin_Area[target_row] <- as.numeric(st_area(merged_geom)) %>% as.integer()
      
      # Remove small polygon
      small_row <- which(st_geometry(tile_sf) == st_geometry(small_poly))
      tile_sf <- tile_sf[-small_row, ]
      
      merged_any <- TRUE
      break  # start over because row indices changed
    }
    
    if (!merged_any) break
  }
  
  return(tile_sf)
}

# ---- Paths ----
dir <- "Z:/0700_NonCore_Funded/0726_TOW_Wales/05_Results/simp_chunks/TOW_Wales_3.gpkg"
output_file <- "Z:/0700_NonCore_Funded/0726_TOW_Wales/05_Results/final_chunks/TOW_Wales_3.gpkg"

# ---- Read ----
gpkg <- st_read(dir, quiet = TRUE) |> 
  st_make_valid()

if (!"geometry" %in% colnames(gpkg)) {
  st_geometry(gpkg) <- gpkg$geom
}

tile_list <- split(gpkg, gpkg$Tilename)

# ---- Parallel setup ----
ncores <- max(1, parallel::detectCores() - 5)
cl <- makeCluster(ncores)
registerDoParallel(cl)

# ---- Process tiles in parallel ----
result_list <- foreach(tile = tile_list,
                       .packages = c("sf", "dplyr", "glue"),
                       .export = c("merge_small_polygons")) %dopar% {
                         
                         fn <- tile$Tilename[[1]]
                         output_file <- glue("Z:/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/4_Processing/final_polys/{fn}.gpkg")
                         
                         if (file.exists(output_file)) {
                           return(NULL)
                          }
                         
                         tile <- st_make_valid(tile)
                         
                         new_tile <- merge_small_polygons(tile)
                         
                         new_tile <- new_tile %>%
                           group_by(TOW_ID) %>%
                           mutate(Sub_ID = row_number()) %>%
                           ungroup() %>%
                           st_cast("POLYGON")
                         
                         final_gpkg <- st_make_valid(new_tile)
                         
                         
                         st_write(final_gpkg,
                                  output_file,
                                  delete_dsn = TRUE,
                                  quiet = TRUE)
                         
                         
                         return(1)
                           
                       }

stopCluster(cl)

# ---- Combine all tiles ----
final_gpkg <- do.call(rbind, result_list)

# Ensure valid geometry
final_gpkg <- st_make_valid(final_gpkg)

# ---- Write final GPKG ----
st_write(final_gpkg,
         output_file,
         delete_dsn = TRUE,
         quiet = TRUE)

cat("Final GPKG saved to:", output_file, "\n")
