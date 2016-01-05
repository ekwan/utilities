#!/bin/bash

echo
echo "Job finished at:"
date
echo

SUCCESS=-1

grep -q "Normal termination" ${1}
let SUCCESS=SUCCESS*$?

grep -q "Stationary point" ${1}
let SUCCESS=SUCCESS*$?

if [ $SUCCESS -eq 0 ]
then
 echo "Job terminated normally..."
 echo
 echo "Final energy:"
 grep -i "SCF Done" ${1} > temp.txt
 energy=`tail -n 1 temp.txt`
 rm temp.txt
 echo $energy
 echo

 grep -q "freq" ${1}
 
 if [ $? -eq 0 ]
 then
  grep -q "imaginary" ${1}
  if [ $? -eq 1 ]
   then
    echo "No imaginary frequencies found." 
   else
    grep imaginary ${1}
  fi
 grep -i "electronic and thermal free" ${1}
 else
  echo "No frequency analysis was performed."
 fi

else
 echo "Abort--Job did not terminate normally."
fi

echo "*************************************"
echo


