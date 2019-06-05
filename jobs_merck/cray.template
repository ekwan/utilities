#!/bin/bash

#PBS -j oe
#PBS -q @QUEUE
#PBS -l select=1:ncpus=@CPU:nppus=@CPU:vntype=cray_compute:proctype=@PROCTYPE
#PBS -N @NAME

cd $PBS_O_WORKDIR

export STMP='/mnt/lustre2/crayxc/scratch'
mkdir $STMP/$USER.$PBS_JOBID || exit_on_error $? "could not create directory"
module load @MODULE
export GAUSS_SCRDIR=$STMP/$USER.$PBS_JOBID
. $g16root/g16/bsd/g16.profile

date
echo running g16

aprun -n 1 -d @CPU g16 @NAME.gjf @NAME.out

date

rm -rf $STMP/$USER.$PBS_JOBID
mv @NAME.out ../output
date > finished.txt
