#!/bin/bash 

#SBATCH -N 1        # number of nodes

#SBATCH -n @CORES   # number of cores 

#SBATCH -p @QUEUE   # partition to submit to 

#SBATCH --mem=64000  # memory per node in MB (see also --mem-per-cpu)

#SBATCH -t @RUNTIME # expected runtime in minutes

#SBATCH -J @JOBNAME # name of this job

# This is a template script for submitting Gaussian jobs
# to SLURM.  Tags starting with @ will be replaced with
# a script.  Eugene Kwan, May 2014

# prevent core dumps on job failure
ulimit -c 0

# set scratch
mkdir $LOCAL_SCRATCH/@USERNAME_@JOBNAME/
export GAUSS_SCRDIR=$LOCAL_SCRATCH/@USERNAME_@JOBNAME

# write out when and where the job started

echo "*************************************" > log.txt
echo "Scratch is: " $LOCAL_SCRATCH/@USERNAME_@JOBNAME >> log.txt
echo "Running on host:" >> log.txt
hostname >> log.txt
echo "Job @JOBNAME started at..." >> log.txt
date >> log.txt

# run job
g09 @JOBNAME.gjf @JOBNAME.out

# remove scratch
rm -rf $LOCAL_SCRATCH/@USERNAME_@JOBNAME

# analyze the result
./analyze.sh @JOBNAME.out >> log.txt

# add it to the master log
cat log.txt >> ../output/output.txt

# move completed job to output directory
mv @JOBNAME.out "../output/"


