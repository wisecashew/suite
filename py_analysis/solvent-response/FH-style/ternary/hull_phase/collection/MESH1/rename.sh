#!/bin/bash

for file in *.pkl; do
    new_name=$(echo "$file" | sed 's/\.\././g')
    mv "$file" "$new_name"
done
