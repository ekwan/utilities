for i in *.gjf; do
    echo $i
    awk -f fix.awk $i > temp.gjf
    mv temp.gjf $i
done
