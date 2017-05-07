#!/bin/bash 

#SBATCH -N 1        # number of nodes

#SBATCH -n @CORES   # number of cores 

#SBATCH -p @QUEUE   # partition to submit to 

#SBATCH --mem=@MEM  # memory per node in MB (see also --mem-per-cpu)

#SBATCH -t @RUNTIME # expected runtime in minutes

#SBATCH -J g09_@JOBNAME # name of this job

# This is a template script for submitting Gaussian jobs
# to SLURM.  Tags starting with @ will be replaced with
# a script.  Eugene Kwan, May 2014

# swap to local node
mkdir /scratch/@USERNAME_@JOBNAME/
export GAUSS_SCRDIR=/scratch/@USERNAME_@JOBNAME

# prevent core dumps on job failure
ulimit -c 0

# write out when and where the job started

echo "*************************************" > log.txt
echo "Running on host:"
hostname >> log.txt
echo "Job @JOBNAME started at..." >> log.txt
date >> log.txt

# run job
g09 @JOBNAME.gjf @JOBNAME.out

# remove temporary directory
rm -rf /scratch/@USERNAME_@JOBNAME

# analyze the result
./analyze.sh @JOBNAME.out >> log.txt

# add it to the master log
cat log.txt >> ../output/output.txt

# move completed job to output directory
mv @JOBNAME.out "../output/"


