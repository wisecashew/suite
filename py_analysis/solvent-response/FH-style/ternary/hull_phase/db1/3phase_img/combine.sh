#!/bin/bash

# Output file for combined data
set -e 
rm combined.db || true 
touch combined.db 

output_file="combined.db"  # Change this if needed
# Extract header from the first line
header="vs | vc | vp | chi_sc | chi_ps | chi_pc | phi_s | phi_p | label0 | label1 | label2 | phi_s1 | phi_p1 | phi_s2 | phi_p2 | phi_s3 | phi_p3 | w1 | w2 | w3"
echo "$header" >> "$output_file"
# Loop through all files with .db extension
for file in chisc*.db; do
  # Append data from line 2 onwards
  tail -n +2 "$file" >> "$output_file"
done

echo "Combined data saved to: $output_file"
