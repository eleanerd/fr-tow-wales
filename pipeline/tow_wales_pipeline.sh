#!/bin/bash

# List files in dir
dir="C:/Users/eleanor.downer/OneDrive - Forest Research/Documents/TOW_Wales/data/"
files=$(ls "$dir")

# Extract tile names
tiles=()
for file in $files; do
    if [[ $file == *_CHM.tif ]]; then
        tile_name=${file%%_*}
        tiles+=("$tile_name")
        
        python Wales_VOM_Pt2.py "$tile_name"
    fi
done