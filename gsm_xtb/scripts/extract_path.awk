# usage: awk -f extract_path.awk stringfile.xyz0001

{
    if ( NF == 1 ) {
        getline
        n_geoms++
        n_atoms=0
        energy[n_geoms] = $1
    }
    else {
        n_atoms++
        geom[n_geoms,n_atoms]=$0
    }
}

END {
    for (i=1; i <= n_geoms; i++) {
        filename=sprintf("analysis/path_%03d.gjf",i)
        #print filename
        print "#\n" > filename
        print energy[i] >> filename
        print "\n0 1" > filename
        for (j=1; j <= n_atoms; j++) {
            print geom[i,j] >> filename
        }
        print "\n" >> filename
    }
}

