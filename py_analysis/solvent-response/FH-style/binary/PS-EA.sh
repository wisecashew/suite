#!/bin/bash

set -e
python nuanced_arr.py --label LOOP --skips 500 --draw-bin \
--xlim 0 0.4 --ylim 350 500 \
--Tbot 350 350 350 --Ttop 500 500 500 --T 200 500 --Tcrit 350 550 \
--vs 2.15333589 2.1617160717778257 2.17391301 --vm 169.97627999 350.9686598655294 600.96728852 \
--pv 1 1 1 --pwms 0.253633 0.252362 0.253864 --pwmm 0 0 0 --pwss 0 0 0 \
--emsa -241.84367196 -241.824049887406 -241.81661662 --emsn -65.12996552 -65.13687424103456 -65.11821106 \
--emma -260.26034119 -259.53263708252217 -260.21919075 --emmn -260.26034119 -259.53263708252217 -260.21919075 \
--essa 0 0 0 --essn 0 0 0 \
--csv PS-ethylacetate-100k.csv PS-ethylacetate_223k.csv PS-ethylacetate_600k.csv \
--img PS-EA.svg --dump PS-EA-100k.dump PS-EA-223k.dump PS-EA-600k.dump 

# this is the 100k good!
# python constricted-cmaes-with_labels.py --label LOOP --skip 1000 \
# --T 10 500 --T-loop 450 --vm 169.97627999 --vs 2.15333589 \
# --xfrmpv 24.17572319 --xfrmpw -1.07932993 \
# --emsa -241.84367196 --emsn -65.12996552 \
# --emm -260.26034119 --ess 0 \
# --img PS-EA-cons-100k --csv PS-ethylacetate-100k.csv --maxiter 10 # > PS-EA-cons-100k.out 


# this is the 223k good!
# python constricted-cmaes-with_labels.py --label LOOP --skip 1000 \
# --T 350 500 --T-loop 500 --vm 350.98078888 --vs 2.16532315 \
# --xfrmpv 24.17044419 --xfrmpw -1.08755902 \
# --emsa -241.83914792 --emsn -65.12509132 \
# --emm -259.55454803 --ess 0 \
# --img PS-EA-cons-223k --csv PS-ethylacetate_223k.csv --maxiter 10 --optimize # > PS-EA-cons-100k.out 

# this is the 600k good!
# python constricted-cmaes-with_labels.py --label LOOP --skip 1000 \
# --T 10 500 --T-loop 500 --vm 600.96728852 --vs 2.17391301 \
# --xfrmpv 24.1587254 --xfrmpw -1.07810775 \
# --emsa -241.81661662 --emsn -65.11821106 \
# --emm -260.21919075 --ess 0 \
# --img PS-EA-cons-600k --csv PS-ethylacetate_600k.csv --maxiter 10 --optimize # > PS-EA-cons-100k.out 
