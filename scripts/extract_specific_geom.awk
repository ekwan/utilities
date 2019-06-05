# retrieves the geometry of a particular step from an optimization
#
# usage: awk -f extract_specific_geom.awk g16.out 25
# (extract geometry number 25)

BEGIN {
    # check number of command line arguments
    if ( ARGC != 3 ) {
        exit_code = 1
        print "Check number of command line arguments."
        print
        print "usage: awk -f extract_specific_geom.awk g16.out 25"
        print "(extract geometry number 25)"
        exit
    }

    input_filename = ARGV[1]
    geometry_number = ARGV[2]
    ARGV[2] = ""
}

# the "second" file is just the geometry number, so don't read that
NR != FNR {
    exit
}

$1 ~ /NAtoms=/ {
    number_of_atoms = $2
}

# keep track of iterations
/Standard orientation:/ {
    iterations++
}

# get geometry
/Standard orientation/,/Rotational constants/ {
    if ( NF == 6 && match($0,"[a-zA-Z]") == 0 ) {
            symbol[$1]=$2
            x[iterations,$1]=$4
            y[iterations,$1]=$5
            z[iterations,$1]=$6
        }
}

# get charge and multiplicity
$1 == "Charge" && $4 == "Multiplicity" {
    charge = $3
    multiplicity = $6
}

END {
    # quit if there was a problem
    if ( exit_code != 0 )
        exit
    
    # generate output filename
    n_fields = split(input_filename,fields,"/")
    abbreviation = fields[n_fields]
    n_fields = split(abbreviation,fields,"-")
    output_filename = ""
    for (i=1; i <= n_fields; i++) {
        if ( i == 6 )
            output_filename = output_filename "ts"
        else
            output_filename = output_filename fields[i]
        if ( i < n_fields )
            output_filename = output_filename "-"
    }
    gsub(".out",".gjf",output_filename)

    #abc/clay_new-glu_gluC6-AADL-Na_1Me2O-phenylate-245-b3lyp_d3bj-631gd-pcm_thf.out
    #clay_new-glu_gluC6-AADL-Na_1Me2O-phenylate-ts-b3lyp_d3bj-631gd-pcm_thf.out

    # print out requested geometry
    print "%nprocshared=16\n%mem=24GB\n#p b3lyp empiricaldispersion=gd3bj 6-31g* opt=(ts,calcfc,noeigentest,maxstep=3) scrf=(pcm,solvent=thf) freq=noraman\n\ntitle\n" > output_filename
    printf "%s %s\n", charge, multiplicity >> output_filename
    for (i=1; i <= number_of_atoms; i++)
        printf "%2d %15.8f %15.8f %15.8f\n", symbol[i], x[geometry_number,i], y[geometry_number,i], z[geometry_number,i] >> output_filename
    printf "\n\n" >> output_filename
    printf "Wrote geometry %d to %s.\n", geometry_number, output_filename
}
