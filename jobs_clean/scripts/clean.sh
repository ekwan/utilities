#!/bin/bash

# This script cleans the checkpoint files from any directory that has the done___.txt file in it.  It doesn't remove anything else.

echo " "
echo "Beginning cleaning process..."
echo "------ "

for directory in *; do
 if [ -d $directory ]
  then
   echo "File $directory is a directory."
   if [ -f ${directory}/done_${directory}.txt ]
    then   
     echo "File done_${directory}.txt exists."
     if [ -f ${directory}/checkpoint.chk ]
      then
       rm -f ${directory}/checkpoint.chk
       rm -f ${directory}/core.*
       echo "Checkpoint file deleted.  Core files, if any, deleted."
      else
       echo "Checkpoint file not found.  No action taken."
     fi
    else
     echo "File done_${directory}.txt DOES NOT exist.  No action taken."
   fi
  else
   echo "File $directory is NOT a directory.  No action taken."
 fi
 echo "-------"
done

echo "Cleaning process complete."
echo " "


