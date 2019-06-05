# pulls out last geometry and generates new file with same
# name and route card
#
# only works for single route jobs and %chk lines not read yet
#
# usage: awk -f extract_geom_with_route.awk *.out

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

$1 == "Charge" && $4 == "Multiplicity" {
    charge[fileCount] = $3
    multiplicity[fileCount] = $6
}

/nprocshared=/ && NF == 1 { processors[fileCount]=$0 }
/mem=/ && NF == 1 { memory[fileCount]=$0 }
/#p/ && route_read[fileCount] == 0 {
    route_read[fileCount] = 1
    while ( match($0,"-----") == 0 ) {
        route_card[fileCount] = route_card[fileCount] $0
        getline
    }
}

END {
    for (i=1; i <= fileCount; i++)
        {
            filename=filenames[i]
            gsub(".out",".gjf",filename)
            print filename
            print trim(processors[i]) > filename
            print trim(memory[i]) >> filename
            print trim(route_card[i]) >> filename
            print "\ntitle\n" >> filename
            printf "%s %s\n", charge[i], multiplicity[i] >> filename
            for (j=1; j <= atoms[i]; j++)
                printf "%2d %15.8f %15.8f %15.8f\n", symbol[i,j], x[i,j], y[i,j], z[i,j] >> filename
            printf "\n\n" >> filename
        }
}

function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function trim(s) { return rtrim(ltrim(s)); }
