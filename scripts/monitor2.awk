# monitors multiple g16 jobs
# assumes the jobs are all of related structures
# assumes #p used in route card
#
# usage:
# awk -f monitor2.awk job*.out

FNR == 1 {
    number_of_files++
    filenames[number_of_files]=FILENAME
    iterations=0
}

# get geometries
/Standard orientation:/ { iterations++ }

/Standard orientation/,/Rotational constants/ {
    if (match($0,"[a-zA-Z]") == 0 && NF==6) {
        # get atom symbols or
        # ensure atom symbols stay the same
        if ( atomSymbol[$1] == 0 )
            atomSymbol[$1]=$2
        else if ( atomSymbol[$1] != $2 ) {
            printf "Expected atom %d to have symbol %s but found symbol %s.\n", $1, atomSymbol[$1], $2
            exit_code=1
            exit
        }

        # get geometries
        number_of_iterations[number_of_files]=iterations
        x[number_of_files,iterations,$1]=$4
        y[number_of_files,iterations,$1]=$5
        z[number_of_files,iterations,$1]=$6
    }
}

$1 ~ /NAtoms=/ {
	# get number of atoms or
	# ensure number of atoms is constant
	if ( number_of_atoms == 0 )
		number_of_atoms = $2
	else if ( iterations > 1 && number_of_atoms != $2 ) {
		printf "Expected %d atoms but found %d.\n", number_of_atoms, $2
		exit_code=1 
		exit
	}
}

# get electronic energies
/SCF Done/ {
    energy[number_of_files,iterations]=$5
    if ( lowest_energy == 0.0 )
        lowest_energy = $5
    else if ( $5 < lowest_energy && (lowest_energy-$5)*627.509469 < 1000.0 )
        lowest_energy = $5
}

/Cartesian Forces:/ {
    cartesian_forces[number_of_files,iterations]=$6
}

END {
    if ( exit_code == 1 ) {
        printf "Error: mismatch in number of atoms for %s.\n", FILENAME
        print "Check that all of these files are for the same structure."
        exit
    }

    # print out energies and gradients for each iteration of each file
    print "filename                                                                                              iteration   rel_energy   cart_force   bond1  bond2"
    
    for (i=1; i <= number_of_files; i++) {
        current_filename = filenames[i]
        if ( length(current_filename) > 100 )
            current_filename = substr(current_filename, length(current_filename)-99, 100)

        gsub(".out","",current_filename)
        for (j=1; j <= number_of_iterations[i]; j++) {
            # calculate relative energy in kcal/mol
            relative_energy = (energy[i,j] - lowest_energy)*627.509469

            # get cartesian force
            current_force = cartesian_forces[i,j]*1E5
            if ( current_force == 0.0 )
                current_force = -999.9

            # calculate geometric parameters
            bond1 = distance(i, j, 19, 34)
            bond2 = distance(i, j, 34, 43)
        
            # print out requested stuff
            current_iteration = j
            if ( energy[i,j] != 0.0 && abs(relative_energy) < 1000.0 ) 
                printf "%100s    %3d  %12.1f    %9.1f      %5.3f  %5.3f\n", current_filename, current_iteration, relative_energy, current_force, bond1, bond2
            else
                printf "%100s    %3d  error\n", current_filename, current_iteration
        }
    }
}

function min(value1,value2) {
    return (value1 < value2 ? value1 : value2)
}

function abs(value) {
  return (value<0?-value:value);
}

function distance(file_number,iteration_number,atomNumber1,atomNumber2,   x1,x2,y1,y2,z1,z2,temp) {
    x1=x[file_number,iteration_number,atomNumber1]
    y1=y[file_number,iteration_number,atomNumber1]
    z1=z[file_number,iteration_number,atomNumber1]
    x2=x[file_number,iteration_number,atomNumber2]
    y2=y[file_number,iteration_number,atomNumber2]
    z2=z[file_number,iteration_number,atomNumber2]
    temp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
    temp=sqrt(temp)
    return temp
}

function angle(file_number,iteration_number,atomNumber1,atomNumber2,atomNumber3,  value) {
    value=( ( -distance(file_number,iteration_number,atomNumber1,atomNumber3)^2 + distance(file_number,iteration_number,atomNumber1,atomNumber2)^2 + distance(file_number,iteration_number,atomNumber2,atomNumber3)^2   ) / ( 2*distance(file_number,iteration_number,atomNumber1,atomNumber2) * distance(file_number,iteration_number,atomNumber2,atomNumber3)   )   )
   return acos(value)
}

function asin(x) { return (180/3.141592)*atan2(x, sqrt(1-x*x)) }

function acos(x) { return (180/3.141592)*atan2(sqrt(1-x*x), x) }

function atan(x) { return (180/3.141592)*atan2(x,1) }

function dihedral(file_number,iteration_number,atomNumber1,atomNumber2,atomNumber3,atomNumber4,   
                     x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,B1x,B1y,B1z,B2x,B2y,B2z,B3x,B3y,B3z,
                     modB2,yAx,yAy,yAz,CP2,CP2y,CP2z,termY,CPx,CPy,CPz,termX,dihed4) {
    x1=x[file_number,iteration_number,atomNumber1]
    y1=y[file_number,iteration_number,atomNumber1]
    z1=z[file_number,iteration_number,atomNumber1]
    x2=x[file_number,iteration_number,atomNumber2]
    y2=y[file_number,iteration_number,atomNumber2]
    z2=z[file_number,iteration_number,atomNumber2]
    x3=x[file_number,iteration_number,atomNumber3]
    y3=y[file_number,iteration_number,atomNumber3]
    z3=z[file_number,iteration_number,atomNumber3]
    x4=x[file_number,iteration_number,atomNumber4]
    y4=y[file_number,iteration_number,atomNumber4]
    z4=z[file_number,iteration_number,atomNumber4]
   
    B1x=x2-x1
    B1y=y2-y1
    B1z=z2-z1
    
    B2x=x3-x2
    B2y=y3-y2
    B2z=z3-z2

    B3x=x4-x3
    B3y=y4-y3
    B3z=z4-z3
   
   modB2=sqrt((B2x^2)+(B2y^2)+(B2z^2))
# yAx is x-coord. etc of modulus of B2 times B1
   yAx=modB2*(B1x)
   yAy=modB2*(B1y)
   yAz=modB2*(B1z)
# CP2 is the crossproduct of B2 and B3
   CP2x=(B2y*B3z)-(B2z*B3y)
   CP2y=(B2z*B3x)-(B2x*B3z)
   CP2z=(B2x*B3y)-(B2y*B3x)
   termY=((yAx*CP2x)+(yAy*CP2y)+(yAz*CP2z))
# CP is the crossproduct of B1 and B2
   CPx=(B1y*B2z)-(B1z*B2y)
   CPy=(B1z*B2x)-(B1x*B2z)
   CPz=(B1x*B2y)-(B1y*B2x)
   termX=((CPx*CP2x)+(CPy*CP2y)+(CPz*CP2z))
   dihed4=(180/3.141592)*atan2(termY,termX)
   return dihed4
}
