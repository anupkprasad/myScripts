#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 13:53:15 2020
@author: anup
"""
import os
#filesRename = ['AccountStatement-2022-06-30.pdf' ]
path = "/mnt/mydrive/alpfaFold3/dynamics/DCLK/DCLK_KD_tail/dclk1_amsa_atemp/allPdbs/"

file_folder_list = os.listdir(path)
# Iterate
for file in os.listdir(path):
	# Checking if the file is present in the list
    oldName = os.path.join(path, file)
    #n = os.path.splitext(file)[0]
    new_name  =   "dclk1_" + file[13:]
    newName = os.path.join(path, new_name)
    # Rename the file uncomment below carefully
    os.rename(oldName, newName)
    print(f"Renamed '{oldName}' to '{newName}'")