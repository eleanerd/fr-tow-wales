#!/bin/bash

# List files in dir
dir="Y:/Forest Inventory/0700_NonCore_Funded/0726_TOW_Wales/04_Spatial Analysis/4_Processing/VOM_Processing/VOM_Pt1"
files=$(ls "$dir")

# Extract tile names
tiles=()
for file in $files; do
    if [[ $file == *_VOM_with_NFI.tif ]]; then
        tile_name=${file%%_*}
        tiles+=("$tile_name")
        
        python Wales_VOM_Pt1_2.py "$tile_name"
    fi
done