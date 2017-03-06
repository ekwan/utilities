#!/bin/bash

title=`awk 'BEGIN {count=0}; /l101.exe/,/Symbolic/ {count++; if (count==3) {print} }' ${1}`

grep -q "Normal termination" ${1}

if [ $? -eq 0 ]
then
 name1=`echo ${1} | cut -f1 -d"."`
 energy=`grep -i "SCF Done" ${1} | tail -n 1 | cut -f8 -d" "`
 free=`grep -i "electronic and thermal free" ${1} | cut -f2 -d"="`
 if [ -z $free ]
  then
   free="0"
  fi
 imaginary=`grep -i frequencies ${1} | grep ignored | awk '{print $1}'`
 if [ -z $imaginary  ]
  then
   imaginary="0"
  fi
 echo $name1, $title, $energy, $free, $imaginary >> batch.txt
else
 echo $name1, $title, error, error, error >> batch.txt
fi



