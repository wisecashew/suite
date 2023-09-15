#!/bin/bash

set -e

# Loop through all .binodal files in the directory
for file in trim*.binodal; do
    echo "Processing $file..."
    
    lcount=`cat $file | wc -l`

    echo "line count = $lcount"
    echo "Finished processing $file!"
done
