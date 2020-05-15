# converts gaussian input files (.gjf) to xyz files
#
# usage:
# mawk -f gjf_to_xyz.awk *.gjf

FNR == 1 {
	n_files++
	filenames[n_files]=FILENAME
    blanks = 0
    last_blank = 0
    while (1 == 1) {
        if ( blanks == 2 && NF == 4 ) {
            n_atoms[n_files]++
            geometry[n_files,n_atoms[n_files]]=$0
        }
        else if ( blanks == 3 )
            break
		is_blank = length(strip($0)) == 0
        if ( is_blank ) {
            if ( last_blank == 0 )
                blanks++
            last_blank = 1
        }
        else
            last_blank = 0
        getline
    }

    nextfile
}

END {
    for (i=1; i <= n_files; i++) {
        input_filename = filenames[i]
        n = split(input_filename, fields, "/")
        input_filename = fields[n]
        output_filename = input_filename
        gsub("gjf$","xyz",output_filename)
        print output_filename
        printf "%d\ntitle\n", n_atoms[i] > output_filename
        for (j=1; j <= n_atoms[i]; j++)
            print geometry[i,j] >> output_filename
    }
}

function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function strip(s) { return rtrim(ltrim(s)); }

