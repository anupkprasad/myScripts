#!/bin/bash

for file in *.pdb; do          # loop through the matching files
                               # no need to worry about white spaces or glob characters in file names
  echo "file=$file"
  sed "s/filename.pdb/$file/" hexdock_template.mac > temp.mac  # put pdb file name at the placeholder filename.pdb in hex_template code
  echo "input changed: $file"
  IFS='.'
  read -ra pdbprifix <<< "$file"
  sed "s/filename_matrix.hex/${pdbprifix}_matrix.hex/" temp.mac > hexinput.mac ###  replace with output *hex at placeholder
  echo "matrix_name changed: ${pdbprifix}_matrix.hex"
  hex <hexinput.mac >hex.log ## call hex with inputfile
done


#
#INFILE=./test.mac
##read inputfile
#while read -r LINE
#do
	#printf '%s\n' "$LINE"
#done < "$INFILE"
#
