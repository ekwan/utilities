# converts maestro xyz files to gjf files
#
# usage: awk -f xyz-to_gjf.awk *.xyz

FNR == 1 {
    file_count++
    i=0
}

FNR > 2 {
    i++
    n_lines[file_count]=i
    line[file_count,i] = $0
}

END {
    for (i=1; i <= file_count; i++) {
        filename = sprintf("gjf/indocat-conformer_%04d-b3lyp_d3bj-631gd-gas.gjf", i)
        print filename
        print "%nprocshared=32" > filename
        print "%mem=48GB" >> filename
        print "#p b3lyp 6-31g* empiricaldispersion=gd3bj opt pop=none" >> filename
        print "\ntitle\n" >> filename
        print "0 1" >> filename
        for (j=1; j <= n_lines[i]; j++)
            print line[i,j] >> filename
        print "\n\n" >> filename
    }
}
