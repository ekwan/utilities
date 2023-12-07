#!/bin/bash

# Gaussian/SLURM Job Submission Script
# Eugene Kwan, May 2014

# Description:
# 
# - creates a directory for your job
# - moves the job input file from the current directory into this new job directory
# - submits the job to SLURM with the specified option
# - shared memory (RAM) on local node is used

# Syntax:
# You must have this file (eek.sh) and the template file (template.sh) in the same directory.
# Execute the script as follows:
#
# . ./eek.sh input_filename number_proc queue_name runtime_in_mins
#
# input_filename: assumes .gjf extension (e.g., write "test01" instead of "test01.gjf")
# number_proc: the number of processors you want (should match input file)
# queue_name: the queue you want to submit to (e.g., jacobsen)
#
# If the intended directory for this job already exists, the script will try to rename the job.

# ensure this gjf file exists

inputFilename=${1}.gjf
if [ -f $inputFilename ]; then

    jobdir=${1}

    # if directory already exists, ask for user input
    while [ -d $jobdir ]; do
        echo Directory $jobdir already exists!
        echo "(O)verwrite, (r)ename, or (s)kip?"
        read choice
        if [ $choice == "o" -o $choice == "O" ]; then
            rm -rf $jobdir
            echo Existing directory deleted.
            echo Note that any existing output files in the output directory will be overwritten.
        elif [ $choice == "r" -o $choice == "R" ]; then
            echo "Enter new directory name: "
            read jobdir
            if [ -z $choice ]; then
                echo Invalid selection!
                echo
                jobdir=${1}
            fi
        else
            echo Skipped.
            return
        fi
    done
    
    # create job directory
    mkdir $jobdir

    # gather information
    username=$(whoami)
    cores=$2
    queue=$3
    # add 2 GB to the memory requested in g09
    mem=`grep %mem $inputFilename | cut -f 2 -d"=" | sed s/GB/000/ | sed s/MB// | awk '{print $1+2000}'`
    if [ -z "$4" ]; then
        runtime=10000 # default runtime is one week
    else
        runtime=$4
    fi
    jobname=$jobdir

    # check all fields present
    if [ -z "$mem" ]; then
        echo Error in memory specification.
        exit 1
    fi

    # create template script
    sed s/@CORES/$cores/g template.sh | sed s/@QUEUE/$queue/g | sed s/@MEM/$mem/g | sed s/@RUNTIME/$runtime/g | sed s/@JOBNAME/$jobname/g | sed s/@USERNAME/$username/g > $jobname.sh

    # copy files to directory
    mv $jobname.sh $jobdir
    mv $inputFilename $jobdir/$jobdir.gjf
    cp analyze.sh $jobdir

    # submit job
    cd $jobdir
    sbatch $jobname.sh
    cd ..

    # report back
    echo "Job submitted with a maximum runtime of "$runtime" minutes."

else
    echo "File "$inputFilename" not found -- job not started."
fi
