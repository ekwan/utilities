#!/bin/bash

# Gaussian Job Submission Script (Merck)
# Eugene Kwan, May 2019
#
# This script submits every .gjf file in directory to PBS.
#
# Usage:
# ./submit.sh number_of_processors desired_queue memory_in_GB
#
# Suggested resources:
# 16 / 256 GB (cluster)
# 24 / 256 GB (crayxc)
#
# fast   10 min
# small   2 h
# medium 12 h
# large  48 h
# huge   infinite
# long   ?

# check number of command line arguments
if (( $# != 3 )); then
    echo Invalid number of command-line arguments.
    echo
    echo Usage: ./submit.sh number_of_processors desired_queue[@processor_type] memory_in_GB
    echo
    echo cray: 256 GB
    echo 16 crayxc@sandybridge
    echo 24 crayxc@haswell
    echo
    echo "cluster: 16 cores/256 GB"
    echo
    echo "fast     10 min"
    echo "small     2 h"
    echo "medium   12 h"
    echo "large    48 h"
    echo "huge    inf  "
    echo "long      ?  "
    echo
    exit 1
fi

# ensure there are input files present
count=`ls -1 *.gjf 2>/dev/null | wc -l`
if (( $count == 0 )); then
    echo Error: No .gjf files found!
    exit 1
fi

# submit all jobs
echo "Beginning batch submission process..."
echo "Each job is being submitted to ${1} processors in the queue ${2}."
echo

for i in *.gjf; do
 
    jobdir=${i%%.gjf}
    echo "Submitting $jobdir..."
 
    # ensure target directory does not exist
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
    
    # make target directory
    mkdir $jobdir

    # determine which template to use
    # and make replacements in template
    cray=`echo ${2} | grep -c cray`
    if (( $cray == 0 )); then
        # this is a cluster job
        sed s/@NAME/$jobdir/g cluster.template | sed s/@CPU/${1}/ | sed s/@QUEUE/${2}/ | sed s/@MEM/${3}/ > temp.sh
    else
        # this is a cray job
        cray_grep=`echo ${2} | grep -c @`
        if (( cray_grep != 1 )); then
            echo No cray architecture specified.  Assuming sandybridge.
            cray_queue=${2}
            cray_architecture="sandybridge"
            cray_module="gaussian/g16.a03-avx3"
        else
            cray_queue=`echo ${2} | awk '{n=split($1,fields,"@"); print fields[1]}'`
            cray_architecture=`echo ${2} | awk '{n=split($1,fields,"@"); print fields[2]}'`
            if [ "$cray_architecture" == "sandybridge" ]; then
                cray_module="gaussian/g16.a03-avx"
            elif [ "$cray_architecture" == "haswell" ]; then
                cray_module="gaussian/g16.a03-avx2"
            else
                echo "Unrecognized cray architecture:" ${cray_architecture}
                echo "Falling back to sandybridge."
                cray_module="gaussian/g16.a03-avx"
                cray_architecture="sandybridge"
            fi
        fi
        sed s/@NAME/$jobdir/g cray.template | sed s/@CPU/${1}/g | sed s/@QUEUE/${cray_queue}/ | sed s/@PROCTYPE/${cray_architecture}/ | sed sh@MODULEh${cray_module}h > temp.sh

    fi

    # submit job
    mv $i ${jobdir}
    mv temp.sh ${jobdir}/${jobdir}.sh
    cd ${jobdir}
    qsub $jobdir.sh
    cd ..
    echo
done

echo "Job submission complete."
