#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-20 (Y/M/D) at 15:44
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""
import argparse

def multi2SingleLineFasta(input_file: str, output_file: str):
    """
    Convert a multi-line FASTA file to a single-line format.
        input_file (str): Path to the input multi-line FASTA file.
        output_file (str): Path to save the single-line FASTA file.
    """
    data = []
    with open(input_file, 'r') as infile:
        first_header = True
        for line in infile:
            if not line.startswith('>'):
                line = line.strip()
                line = line.replace(' ', '')  # Remove spaces
                if line:  # Skip empty lines
                    data.append(line)
            else:
                if first_header:
                    data.append(line.strip() + '\n')  # First header without leading newline
                    first_header = False
                else:
                    data.append('\n' + line.strip() + '\n')  # Subsequent headers with leading newline
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

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Convert multiline fasta to single line")
	parser.add_argument("input", help="input fasta file")
	parser.add_argument("output", help="output fasta file")
	args = parser.parse_args()
	multi2SingleLineFasta(args.input, args.output)
	print("done")
