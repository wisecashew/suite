#!/bin/bash

set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?." ' EXIT
rm wu.out || true 
rm pi.out || true 

for file in *wu-idx*; do
	# echo "$file"
	idx=$(echo "$file" | sed -n 's/.*idx_\([0-9]\+\).*/\1/p')
	if `echo "$file" | grep -q "0\.1"`; then
		echo "$file with 0.1"
		python pickle_parser.py -i $file -o with_wu_$idx >> wu.out 2>&1 
	else
		echo "$file without 0.1"
		python pickle_parser.py -i $file -o without_wu_$idx >> wu.out 2>&1 
	fi
done


for file in *pi-idx*; do
	# echo "$file"
	idx=$(echo "$file" | sed -n 's/.*idx_\([0-9]\+\).*/\1/p')
	if `echo "$file" | grep -q "0\.1"`; then
		echo "$file with 0.1"
		python pickle_parser.py -i $file -o with_pi_$idx >> pi.out 2>&1 
	else
		echo "$file without 0.1"
		python pickle_parser.py -i $file -o without_pi_$idx >> pi.out 2>&1 
	fi
done
