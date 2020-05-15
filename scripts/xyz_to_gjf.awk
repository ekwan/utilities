# converts xyz files to gjf files
#
# usage: awk -f xyz-to_gjf.awk *.xyz

FNR == 1 {
    n_files++
    filenames[n_files]=FILENAME
    i=0
}

FNR > 2 {
    i++
    n_lines[n_files]=i
    line[n_files,i] = $0
}

END {
    for (i=1; i <= n_files; i++) {
        filename = filenames[i]
        n=split(filename,fields,"/")
        filename = fields[n]
        gsub("xyz$","gjf",filename)
        print "#\n" > filename
        print "title\n" >> filename
        print "0 1" >> filename
        for (j=1; j <= n_lines[i]; j++)
            print line[i,j] >> filename
        print "\n\n" >> filename
    }
}
