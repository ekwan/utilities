#!/bin/bash

# Takes modredundant jobs with freq checkpoints and starts ts jobs.
# Failed jobs or duplicates will be ignored.

##################

# where the output files are
output_file_location=/share/PI/nburns/jobs4/output

# where the checkpoints with the frequencies are
checkpoint_file_directory=/share/PI/nburns/jobs4

# all output files in output_file_location matching this will be submitted
output_file_mask=indocat_mannosyl_diax_BnOH_???.out

# all checkpoints will be called this
checkpoint_standard_filename=checkpoint.chk

# slurm settings
username=$(whoami)
cores=16
queue=normal
mem=60000
runtime=2880

###################

echo "Starting job submission process..."
for i in ${output_file_location}/${output_file_mask}; do
    echo
    echo Submitting job for $i.
    
    # make sure a directory is not already present
    rawdir=`echo $i | awk 'BEGIN {FS="/"} {gsub(".out", "", $0); print $NF }'`
    jobdir=${rawdir}_ts
    if [ -d $jobdir ]; then
        echo Directory already present.
        continue
    fi
    
    # make sure a checkpoint is available
    source_checkpoint_filename=${checkpoint_file_directory}/${rawdir}/${checkpoint_standard_filename}
    if [ ! -f ${source_checkpoint_filename} ]; then
        echo Checkpoint file ${source_checkpoint_filename} not found!
        continue
    fi

    # make sure a frequency job was completed
    freqs=`grep -c Frequencies $i`
    if [ $freqs -lt 10 ]; then
        echo Frequencies not found.
        continue
    fi
 
    # create job directory and copy checkpoint
    mkdir $jobdir
    destination_checkpoint_filename=${jobdir}/${checkpoint_standard_filename}
    cp ${source_checkpoint_filename} ${destination_checkpoint_filename}

    # copy in job file
    cp freq_to_ts.gjf ${jobdir}/${jobdir}.gjf
    cp analyze.sh ${jobdir}

    # create slurm script
    slurm_script=${jobdir}.sh
    sed s/@CORES/$cores/g template.sh | sed s/@QUEUE/$queue/g | sed s/@MEM/$mem/g | sed s/@RUNTIME/$runtime/g | sed s/@JOBNAME/$jobdir/g | sed s/@USERNAME/$username/g >    $slurm_script
    mv ${slurm_script} ${jobdir}

    # start job
    cd ${jobdir}
    sbatch ${slurm_script}
    cd ..

done

echo "Job submission complete."
