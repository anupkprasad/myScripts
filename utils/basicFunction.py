#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-29 (Y/M/D) at 15:28
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

from tsmoothie.smoother import *

# This function calculates the running mean of a given array `y` using a convolution smoother.
def runningMean(y, N):
    smoother1 = ConvolutionSmoother(window_len=N, window_type='ones')
    smoother1.smooth(y)
    return smoother1.smooth_data[0]

def renameFiles(path, prefix="dclk1_"):
    """
    Renames files in the specified directory by adding a prefix to each file name.
    
    Args:
        path (str): The directory containing the files to rename.
        prefix (str): The prefix to add to each file name (default is "dclk1_").
    """
    import os
    
    file_folder_list = os.listdir(path)
    
    # Iterate through each file in the directory
    for i, file in enumerate(file_folder_list):
        if i < 10:  # Skip the first 10 files
            continue
        oldName = os.path.join(path, file)
        new_name = f"6G{i}_H" ##prefix + file[13:]  # Adjust the slicing as needed
        newName = os.path.join(path, new_name)
        # Rename the file
        os.rename(oldName, newName)
        print(f"Renamed '{oldName}' to '{newName}'")


if __name__ == "__main__":
    # Example usage
    # renameFiles("/path/to/files/", prefix="my_prefix_")
    pass
