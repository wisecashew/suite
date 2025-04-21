#!/bin/bash

set -e

./sct -D 32 -T 1.0 -N 208 -p 0.25 -M -1 -m -1 -S -1 -s 0 --alpha 6.64515447931444 --beta 67.06061106873511 -d thermo.dump
