#!/bin/bash

# Gaussian Job Submission Script
# Eugene Kwan, May 2014
#
# This script submits every .gjf file in directory to SLURM.
#
# Usage:
# ./submit.sh number_of_processors desired_queue expected_runtime(mins)
#
# Suggested number of cores:
# jacobsen 24, set nprocshared=12
# general 64, set nprocshared=64

# check number of command line arguments
if (( $# < 2 )); then
    echo Invalid number of command-line arguments.
    echo
    echo Usage: ./submit.sh number_of_processors desired_queue expected_runtime_in_mins
    exit 1
fi

echo "Beginning batch submission process..."
echo "Each job is being submitted to ${1} processors in the queue ${2}."
echo

for i in *.gjf; do
 
 echo "Submitting job file ${i%%.gjf}..."
 . ./eek.sh ${i%%.gjf} ${1} ${2} ${3}
 echo
done

echo "Job submission complete."
