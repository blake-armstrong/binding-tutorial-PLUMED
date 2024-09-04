#!/bin/bash
tmp=$(plumed kt --temp 300 | awk '{print $11" " $12}')
echo "kBT is ${tmp}"
kt=$(echo ${tmp} | awk '{print $1}')
cat BIAS/HILLS* > HILLS
plumed sum_hills --hills HILLS --mintozero --kt ${kt} --outfile fes_1D.dat
radius=0.1
python3 ../../scripts/mtd_analysis.py fes_1D.dat 14473.0 300 ${radius} | tee dG.dat
