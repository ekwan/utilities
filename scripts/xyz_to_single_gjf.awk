# converts multiple xyz files to one gjf file
#
# usage: awk -v output_filename="$var" -f xyz_to_single_gjf.awk *.xyz

BEGIN {
    if (length(output_filename) == 0) {
        print "must specify output_filename"
        print
        print "usage: awk -v output_filename=\"$var\" -f xyz_to_single_gjf.awk *.xyz"
        exit_code = 1
        exit
    }
}

FNR == 1 {
    n_files++
    n_atoms[n_files]=$1
}

FNR == 2 {
    titles[n_files]=$0
}

FNR > 2 && NF == 4 {
    geometry[n_files,FNR-2] = $0
}

END {
    if (exit_code != 0)
        exit exit_code
    printf "" > output_filename
    for (i=1; i <= n_files; i++) {
        if ( i > 1 )
            print "--Link1--" >> output_filename
        print "%nprocshared=4\n%mem=3GB\n#t b3lyp empiricaldispersion=gd3bj gen NMR pop=none int=finegrid\n" >> output_filename
        printf "%s\n\n", titles[i] >> output_filename
        print "0 1" >> output_filename
        for (j=1; j <= n_atoms[i]; j++)
            print geometry[i,j] >> output_filename
        print "\n@/n/jacobsen_lab/ekwan/solvent/basis/pcSseg-1.bas\n" >> output_filename
    }
    printf "\n\n" >> output_filename
}
