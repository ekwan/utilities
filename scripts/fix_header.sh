for i in *.gjf; do
    awk -f fix_header.awk head.txt $i tail.txt > temp.gjf
    mv temp.gjf $i
    echo $i
done
