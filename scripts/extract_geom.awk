FNR == 1 {
    fileCount++
    filenames[fileCount]=FILENAME
}

/Standard orientation/,/Rotational constants/ {
    if ( NF == 6 && match($0,"[a-zA-Z]") == 0 )
        {
            symbol[fileCount,$1]=$2
            x[fileCount,$1]=$4
            y[fileCount,$1]=$5
            z[fileCount,$1]=$6
            atoms[fileCount]=$1
        }
}

END {
    for (i=1; i <= fileCount; i++)
        {
            filename = filenames[i]
            gsub(".out",".gjf",filename)
            print filename
            print "%nprocshared=16\n%mem=50GB\n#p m062x/6-31g* opt freq=noraman scrf=(pcm,solvent=toluene)\n\ntitle\n\n0 1" > filename
            for (j=1; j <= atoms[i]; j++)
                printf "%2d %15.8f %15.8f %15.8f\n", symbol[i,j], x[i,j], y[i,j], z[i,j]  >> filename
            printf "\n\n" >> filename
        }
}
