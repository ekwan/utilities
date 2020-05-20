# awk -f fix_header.awk head.txt abc.gjf tail.txt

FNR == 1 {
    fileCount++
    spaces = 0
}

fileCount == 1 || fileCount == 3 {print}

fileCount == 2 {
    if ( length($0) == 0 )
        spaces++
    if ( spaces == 2 )
        {
            getline
            printf "\n"
            spaces += 2
        }
    if ( spaces == 4)
        print $0
}
