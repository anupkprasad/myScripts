#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-20 (Y/M/D) at 15:44
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

def multi2SingleLineFasta(input_file: str, output_file: str):
    """
    Convert a multi-line FASTA file to a single-line format.
        input_file (str): Path to the input multi-line FASTA file.
        output_file (str): Path to save the single-line FASTA file.
    """
    data = []
    with open(input_file, 'r') as infile:
        for line in infile:
            if not line.startswith('>'):
                line = line.strip()
                line = line.replace(' ', '')  # Remove spaces
                if line:  # Skip empty lines
                    data.append(line)
            else:
                data.append('\n' + line.strip() + '\n')  # Write header lines
    # Write to output file in single-line format
    with open(output_file, 'w') as outfile:
        for line in data:
            if line.startswith('>'):
                outfile.write(line)
            else:
                outfile.write(line)
        # Ensure the output file ends with a newline; needed for fa2cma
        if data and not data[-1].endswith('\n'):
            outfile.write('\n')
    print(f"Converted {input_file} to single-line format in {output_file}.")



path = "/mnt/mydrive/alpfaFold3/dynamics/PRKACA(PKA)"
multi2SingleLineFasta(input_file=path + "/PRKACA_KD_tail_query.fasta", output_file= path + "/PRKACA_KD_tail_query.fasta")