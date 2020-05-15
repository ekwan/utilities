#!/bin/bash

# runs one pair of reactant/product geometries for a GSM TS search with xtb

# - if molecule has non-standard charge or multiplicity, modify the ograd file
# - if breaking bonds (or not), adjust the TS_FINAL_TYPE parameter in inpfileq 
# - ensure atom numbers are the same in both xyz files
# - ideally xyz files should be minimized with xtb first
# - job will run locally, assumes the current session will not die (use "screen")

# usage:
#
# ./submit_one.sh jobname reactant.xyz product.xyz 

if [ "$#" -ne 3 ]; then
    echo Invalid number of command-line arguments.
    echo
    echo usage:
    echo ./submit_one.sh jobname reactant.xyz product.xyz
    exit 1
fi

# extract command line arguments
jobname=$1
reactant=$2
product=$3

# check if directory already exists
if [ -d $jobname ]; then
    echo Job directory $jobname already exists!  Aborting.
    exit 1
fi

# check that reactant and product are present
if [ ! -f $reactant ]; then
    echo Reactant file $reactant does not exist!  Aborting.
    exit 1
fi
if [ ! -f $product ]; then
    echo Product file $product does not exist!  Aborting.
    exit 1
fi

# setup job directory
cp -r template $jobname
mv $reactant $jobname
mv $product $jobname
cd $jobname
mkdir scratch
cat ${reactant} ${product} > scratch/initial0001.xyz
cp xtb_to_orca.py scratch

# load necessary modules
module load cmake/3.10.2 intel/mkl/1.1 git openmpi/4.0.1

# 4 is the number of cores that xtb will use
# stacksize might be important for huge inputs
export OMP_STACKSIZE=2G
export OMP_NUM_THREADS=4,1
export OMP_MAX_ACTIVE_LEVELS=1
export MKL_NUM_THREADS=4
export XTBHOME=/SFS/user/bm/kwaneu/sw/anaconda2/envs/xtb/share/xtb
export XTBPATH=/SFS/user/bm/kwaneu/sw/anaconda2/envs/xtb/share/xtb

# 8 is supposed to be how many cores GSM will use, but doesn't seem to have much effect
#./gsm.orca.exe 1 8
date >> gsm.log 
stdbuf -o 0 ./gsm.orca.exe 1 8 >> gsm.log 2>&1 &
echo Now running job \"${jobname}\".
