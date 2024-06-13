#!/bin/bash

# Define a function to check the conditions for each database

# Iterate over all databases
for database in chisc*.db; do
    # Call the function for each database
    awk -F'|' -f db_parser.awk $database
done
