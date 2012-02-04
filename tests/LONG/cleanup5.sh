#!/bin/bash

# change rain to your user name on blade or altix wherever you are going to use it
ps aux | sed -n '/maerzkek/ p' > tmp
# change the program name whatever name you generated with the make file. I am assuing you are running the program with the same name. 
cat tmp | sed -n '/NpTgaro/ p'  > tmp1

# Dont need to do anything here.

for PIDS in $(awk '{print $2}' < tmp1 )
do

 kill -9 $PIDS

done

exit 0
