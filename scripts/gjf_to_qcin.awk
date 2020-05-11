FNR == 1 {
    # keep track of number of files and filenames
    n_files++
    filenames[n_files]=FILENAME

    # initialize
    last_blank = 0
    blanks = 0
    charge[n_files] = "X"
    multiplicity[n_files] = "Y"
    n_atoms[n_files] = 0
    this_n_atoms = 0

    # loop while true
    while (1==1) {
        if (NF==0) {
            # skip header
            if ( last_blank == 0 )
                blanks++
            last_blank = 1
        }
        else
            last_blank = 0
        if ( blanks == 2 ) {
            if ( charge[n_files] == "X" ) {
                # read charge and multiplicity
                getline
                charge[n_files] = $1
                multiplicity[n_files] = $2
            }
            else {
                # read geometry
                this_n_atoms++
                n_atoms[n_files] = this_n_atoms
                geometry[n_files,this_n_atoms] = $0
            }
        }
        else if ( blanks > 2 )
            # skip anything at the end
            nextfile
        getline
    }
}

END {
    # iterate through files
    for (i=1; i <= n_files; i++) {
        # construct output filename
        input_filename = filenames[i]
        output_filename = input_filename
        n=split(output_filename, fields, "/")
        output_filename = fields[n]
        gsub(".gjf","",output_filename)
        output_filename = output_filename ".qcin"
        print output_filename

		# determine if this is a ts job
        is_ts_job = match(output_filename,"-ts")

        # job step 1
        printf "$molecule\n" > output_filename
        printf "%s %s\n", charge[i], multiplicity[i] >> output_filename
        for (j=1; j <= n_atoms[i]; j++) {
            split(geometry[i,j],fields," ")
            printf "%2s %16.8f %16.8f %16.8f\n", fields[1], fields[2], fields[3], fields[4] >> output_filename
        }
        print "$end\n" >> output_filename
        print "$pcm\nTHEORY CPCM\nMETHOD SWIG\n$end\n" >> output_filename
        print "$solvent\nDIELECTRIC 78.39\n$end\n" >> output_filename
        print "$rem\nBASIS 6-31G*\nDFT_D D3_BJ\nMETHOD B3LYP\nMEM_TOTAL 3000" >> output_filename
        print "SOLVENT_METHOD PCM" >> output_filename
        if ( is_ts_job > 0 )
            print "JOBTYPE FREQ" >> output_filename
        else
            print "JOBTYPE OPT" >> output_filename
        print "$end\n" >> output_filename

        # job step 2
        if ( is_ts_job > 0 ) {
            print "@@@\n\n$molecule\nread\n$end\n" >> output_filename
            print "$pcm\nTHEORY CPCM\nMETHOD SWIG\n$end\n" >> output_filename
            print "$solvent\nDIELECTRIC 78.39\n$end\n" >> output_filename
            print "$rem\nBASIS 6-31G*\nDFT_D D3_BJ\nMETHOD B3LYP\nMEM_TOTAL 3000" >> output_filename
            print "SOLVENT_METHOD PCM" >> output_filename
            print "SCF_GUESS READ\nJOBTYPE TS\nGEOM_OPT_HESSIAN READ" >> output_filename
            print "$end\n" >> output_filename
        }

        # job step 3
        print "@@@\n\n$molecule\nread\n$end\n" >> output_filename
        print "$pcm\nTHEORY CPCM\nMETHOD SWIG\n$end\n" >> output_filename
        print "$solvent\nDIELECTRIC 78.39\n$end\n" >> output_filename
        print "$rem\nBASIS 6-31G*\nDFT_D D3_BJ\nMETHOD B3LYP\nMEM_TOTAL 3000" >> output_filename
        print "SOLVENT_METHOD PCM" >> output_filename
        print "SCF_GUESS READ\nJOBTYPE FREQ" >> output_filename
        print "$end\n" >> output_filename
    }
}
