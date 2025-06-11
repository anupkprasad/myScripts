#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 13:53:15 2020

@author: anup
"""

import os
#filesRename = ['AccountStatement-2022-06-30.pdf' ]
path = "/home/anupkumar/Downloads/"

file_folder_list = os.listdir(path)
# Iterate
for file in os.listdir(path):
	# Checking if the file is present in the list
	#if file in filesRename:
  	oldName = os.path.join(path, file)
  	#n = os.path.splitext(file)[0]
  	b =   file[:34] + "1min_for10" + file[42:] ##'GUVs_peptide0.5uM' + file[12:]
  	newName = os.path.join(path, b)
  	# Rename the file uncomment below carefully
  	#os.rename(oldName, newName)



for file in file_folder_list:
    if file.endswith(".docx"):
        print(file)