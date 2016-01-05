#!/bin/bash

# changes nproc in every file to specified number

for i in *.gjf; do
 
 echo "Fixing job file $i..."
 sed 3s/.*/"%nprocshared=8"/  $i > temp.txt
 rm $i
 mv temp.txt $i

done

echo "Done."
