#!/bin/bash

set -e

# Loop through all .binodal files in the directory
for file in vs*.binodal; do
    echo "Processing $file"
    
    # Perform the operation to remove duplicates and update the file in-place
    awk '!seen[$0]++' "$file" > temp.txt && mv temp.txt "$file"
    
    echo "Finished processing $file"
done
