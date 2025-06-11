#!/bin/bash
mkdir hexdock

for file in ./pdbdir/*.pdb; do          # loop through the matching files
                               # no need to worry about white spaces or glob characters in file names
  echo "file=$file"
  IFS='/'
  read -ra pdbfull <<< "$file"
  IFS='.'
  read -ra pdbprifix <<< "${pdbfull[2]}"

  sed "s/381.pdb/${pdbfull[2]}/" hexdock_template.mac > temp.mac
  echo "input changed: ${pdbfull[2]}"
  sed "s/381_matrix.hex/${pdbprifix}_matrix.hex/" temp.mac > hexinput.mac
  echo "matrix_name changed: ${pdbprifix}_matrix.hex"
  hex <hexinput.mac >hex.log
done


#
#INFILE=./test.mac
##read inputfile
#while read -r LINE
#do
	#printf '%s\n' "$LINE"
#done < "$INFILE"
#
