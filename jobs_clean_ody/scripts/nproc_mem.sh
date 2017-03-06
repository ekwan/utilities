#!/bin/bash

# changes nproc in every file to specified number

for i in *.gjf; do
 
 echo "Fixing job file $i..."
 awk '{
        if ( match($0,"%nproc") > 0 )
            $0 = "%nprocshared=8"
        else if ( match($0,"%mem") > 0 )
            $0 = "%mem=5GB" 
        print $0
      }
 ' $i > temp.txt
 rm $i
 mv temp.txt $i

done

echo "Done."
