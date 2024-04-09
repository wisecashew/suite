#!/bin/bash

set -e 

python database.py --hull collection/MESH_5001/chisc_1-chips_1-chipc_1-vs_2-vp_3-vc_3.hull.pkl 
python compare.py  --hull collection/MESH_5001/chisc_1-chips_1-chipc_1-vs_2-vp_3-vc_3.hull.pkl --tget triple.db --tdump trips_valid.db

