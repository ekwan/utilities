#!/bin/bash

echo
echo "filename, title_card, electronic_energy, free_energy, num_imag_freq" > batch.txt

for i in *.out; do
 echo -ne "."
 . ./read.sh $i 
done

echo
awk 'BEGIN {FS=", "} {printf "%50-s %40-s %15s %15s %3s\n", $1, $2, $3, $4, $5}' batch.txt > batch2.txt
mv batch2.txt batch.txt
more batch.txt

echo
