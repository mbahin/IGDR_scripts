#!/bin/bash

# Mathieu Bahin, 30/06/14

# Script to rename files within directory and sub-directories.
# The inputs are the path from where the script has to find all the files, directories and sub-directories, the current name to change and the new name.
# There is no output.
# When the name to change is, for example 'LUPA1' and there may be files with names like 'LUPA10', the user has to be careful not to change the second kind of file when processing the first one. This is why there is the second condition in the find. The first is for the cases where a directory has the exact name the user want to change.

# Getting parameters back
#path='/home/genouest/umr6061/recomgen/dog/hitte/JurasHic/FASTQ.dir/BROAD.dir/BROAD.trim_paired.dir'
#old='BROAD_BLOOD'
#new='BROAD-ANN-2011_BLOOD_BEGL'
root=$1
old=$2
new=$3

# List files and directories and move them
for f in $(find $root -depth -name "*$old" -o -name "*$old[^0-9]*"); do
	new_path=$(echo $f | sed "s|\(.*\)/\(.*\)$old\(.*\)|\1/\2$new\3|g")
	mv $f $new_path
done