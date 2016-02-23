# shows how to epimerize two groups
# awk -f epimerize.awk *.out
# write out a bunch of gjf files to a gjf folder

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
            filename = "gjf/" filenames[i]
            gsub(".out","_proper.gjf",filename)
            print filename
            print "%nprocshared=4\n%mem=3GB\n#p m062x/cc-pVTZ opt scrf=(pcm,solvent=toluene)\n\ntitle\n\n0 1" > filename

            # get C6-H10 bond length
            CH_bond_length = distance(i,6,10)

            # get C6-Cl1 bond length
            CCl_bond_length = distance(i,6,1)

            # switch the positions of atom H10 and Cl1
            H_x = x[i,10]
            H_y = y[i,10]
            H_z = z[i,10]
            x[i,10] = x[i,1]
            y[i,10] = y[i,1]
            z[i,10] = z[i,1]
            x[i,1] = H_x
            y[i,1] = H_y
            z[i,1] = H_z

            # set C6-H10 bond length
            set_distance(i, 6, 10, CH_bond_length)

            # set C6-Cl1 bond length
            set_distance(i, 6, 1, CCl_bond_length)

            for (j=1; j <= atoms[i]; j++)
                printf "%2d %15.8f %15.8f %15.8f\n", symbol[i,j], x[i,j], y[i,j], z[i,j]  >> filename
            printf "\n\n" >> filename
        }
}

function set_distance(file_number, atom_number_1, atom_number_2, requested_distance,       scaling, d_x, d_y, d_z)
{
    current_distance = distance(file_number, atom_number_1, atom_number_2)
    scaling = (requested_distance - current_distance) / current_distance
    d_x = x[file_number,atom_number_2] + ((x[file_number,atom_number_2] - x[file_number,atom_number_1]) * scaling)
    d_y = y[file_number,atom_number_2] + ((y[file_number,atom_number_2] - y[file_number,atom_number_1]) * scaling)
    d_z = z[file_number,atom_number_2] + ((z[file_number,atom_number_2] - z[file_number,atom_number_1]) * scaling)   
    x[file_number,atom_number_2] = d_x
    y[file_number,atom_number_2] = d_y
    z[file_number,atom_number_2] = d_z
}

function distance(thisFilename,atomNumber1,atomNumber2) {
    x1=x[thisFilename,atomNumber1]
    y1=y[thisFilename,atomNumber1]
    z1=z[thisFilename,atomNumber1]
    x2=x[thisFilename,atomNumber2]
    y2=y[thisFilename,atomNumber2]
    z2=z[thisFilename,atomNumber2]
    temp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
    temp=sqrt(temp)
    return temp
}
