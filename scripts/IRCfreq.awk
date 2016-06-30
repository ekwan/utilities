# extracts each optimized point from an IRC job and creates a single point frequency job
#
# usage:
# awk -f IRCfreq.awk IRC.out
#
# should only be used on one file

FNR != FNR { 
    print "warning: ignoring additional input files beyond the first"
    exit
}

$1 == "Charge" && $4 == "Multiplicity" {
    charge = $3
    multiplicity = $5
}

/Cartesian Coordinates/,/CHANGE IN THE REACTION COORDINATE/ {
    if ( match($1,"[0-9]") > 0 )
        {
            atomNumber = $1
            if ( atomNumber == 1 )
                points++
            symbol[atomNumber] = $2
            numberOfAtoms = $1
            x[points,atomNumber] = $3
            y[points,atomNumber] = $4
            z[points,atomNumber] = $5i
        }
}

END {
    for (i=1; i <= points; i++)
        {
            filename = sprintf("gjf/IRC_reverse_%03d.gjf", i)
            print "%nprocshared=4" > filename
            print "%mem=3GB" >> filename
            print "#p b3lyp/6-31g* empiricaldispersion=gd3bj scrf=(pcm,solvent=water) freq=noraman pop=none" >> filename
            print "\ntitle\n" >> filename
            printf "%d %d\n", charge, multiplicity >> filename
            for (atom=1; atom <= numberOfAtoms; atom++)
                printf "%3d %15.8f %15.8f %15.8f\n", symbol[atom], x[i,atom], y[i,atom], z[i,atom] >> filename
            printf "\n\n" >> filename
        }
    printf "%d gjf files written.\n", points
}
