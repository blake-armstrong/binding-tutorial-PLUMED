#!/bin/bash

error_exit() {
    echo "Error: $1" >&2
    exit 1
}

if [ $# -gt 0 ]; then
    if ! [[ $1 =~ ^[0-9]+$ ]]; then
        error_exit "Number of walkers should be a positive integer."
    fi
    
    if [ $1 -le 1 ]; then
        error_exit "Need more than one walker for multiple walkers."
    fi
    
fi

numwalkers=${1:-4}

names=("radius_0.0" "radius_0.1" "radius_0.2" "radius_0.3" "radius_0.4" "radius_0.5")
radiuses=(0.0 0.1 0.2 0.3 0.4 0.5)





cwd=$(pwd)
mkdir -p data
cd data
for i in ${!names[@]}; do
  name=${names[$i]}
  radius=${radiuses[$i]}
  mkdir -p ${name}
  cd ${name}
  sed "s/REPLACERADIUS/${radius}/" ${cwd}/simulation_files/plumed_template.inp > plumed.inp
  sed "s/REPLACERADIUS/${radius}/" ${cwd}/scripts/analysis.sh > analysis.sh
  cp ${cwd}/scripts/run.sh .
  mkdir -p BIAS
  if [ -n "$(ls -A BIAS/HILLS* 2>/dev/null)" ]; then
    cd BIAS
    for b in HILLS*; do
      mv ${b} bck.${b}
    done
    cd ..
  fi
  for j in $(seq 0 $((numwalkers-1))); do
    mkdir -p Walker_${j}
    cd Walker_${j}/
    rand=$((25 + RANDOM%(14)))
    sed "s/YY/${rand}/" ${cwd}/simulation_files/ref.pdb > coord.pdb
    # cp ${cwd}/simulation_files/${j}.pdb coord.pdb
    cp ../plumed.inp .
    cp ${cwd}/simulation_files/argon.xml .
    sed -i.bak "s/REPLACETOTALWALKERS/${numwalkers}/" plumed.inp
    sed -i.bak "s/REPLACEWALKERID/${j}/" plumed.inp && rm plumed.inp.bak
    touch ../BIAS/HILLS.${j}
    cd ..
  done
  cd ..
done
cd ..
