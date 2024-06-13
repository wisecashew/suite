#!/bin/bash

dir1="."
dir2="hulls"

# Iterate over files in dir1
for file in "$dir1"/*.hull.pkl; do
    # Extract the filename from the path
    filename=$(basename "$file")
    
    # Check if the file exists in dir2
    if [ ! -e "$dir2/$filename" ]; then
        echo "$filename does not exist in dir2"
    fi
done
