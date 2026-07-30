#!/bin/bash

# List files in dir
dir="Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/4_Processing/VOM_Processing/Hedges/"
files=$(ls "$dir")

# Extract tile names
tiles=()
for file in $files; do
    if [[ $file == *_hedges_dissolved.gpkg ]]; then
        tile_name=${file%%_*}
        tiles+=("$tile_name")
        
        python Wales_VOM_Pt2.py "$tile_name"
    fi
done