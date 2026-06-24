#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-20 (Y/M/D) at 15:44
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""
import argparse
import re
from textwrap import fill


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



def a3mtoFasta(inputfile, outputfile, insert=True):
    """
    Converts an A3M MSA file to FASTA format by removing lowercase letters.
    Args:
        inputfile (str): Path to the A3M file.
        outputfile (str): Path to the output FASTA file.
    Returns:
        None
    """
    with open(inputfile, 'r') as infile, open(outputfile, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                outfile.write(line)
            else:
                if insert:
                    cleaned_line = ''.join(c for c in line if not c.islower())
                else:
                    cleaned_line = line
                outfile.write(cleaned_line)
    print(f"Converted {inputfile} to {outputfile} in FASTA format.")


def cmatoFasta(infile, outfile, width=60, header_token_first=False):
    """
    Convert a CMA alignment file to FASTA format in a single pass.
    
    Args:
        infile (str): Path to input CMA file
        outfile (str): Path to output FASTA file
        width (int): Line width for wrapped sequences (default: 60)
        header_token_first (bool): If True, use only the first token of the header
    
    Extracts sequences from CMA blocks and writes them to FASTA.
    """
    # Read the input CMA file
    with open(infile, 'r', encoding='utf-8', errors='ignore') as fh:
        text = fh.read()
    
    # Extract sequences from CMA text
    seqs = []
    for m in re.finditer(r'>([^\n\r]+)', text):
        header = m.group(1).strip()
        start_pos = m.end()
        
        # Find the next '{' that starts the sequence block
        open_brace = text.find('{', start_pos)
        if open_brace == -1:
            seqs.append((header, ""))
            continue
        
        # Prefer sequence end pattern '()}'; fallback to next '}'
        close_marker = text.find('()}', open_brace)
        if close_marker != -1:
            seq_block = text[open_brace+1:close_marker]
        else:
            close_brace = text.find('}', open_brace)
            if close_brace == -1:
                seq_block = text[open_brace+1:]
            else:
                seq_block = text[open_brace+1:close_brace]
        
        # Extract only letters and uppercase
        seq_letters = re.sub(r'[^A-Za-z]', '', seq_block).upper()
        seqs.append((header, seq_letters))
    
    # Write FASTA file
    n_empty = 0
    with open(outfile, 'w') as out:
        for header, seq in seqs:
            if header_token_first:
                fasta_header = header.split()[0]
            else:
                fasta_header = header.replace('>', '').strip()
            
            out.write(f'>{fasta_header}\n')
            if seq:
                out.write(fill(seq, width=width) + '\n')
            else:
                n_empty += 1
                out.write('\n')
    
    # Report results
    if n_empty:
        print(f"Warning: {n_empty} sequence(s) had no letters extracted.")
    print(f"Wrote {len(seqs)} entries to {outfile}. {n_empty} empty sequences.")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert multiline fasta to single line")
    parser.add_argument("input", help="input fasta file")
    parser.add_argument("output", help="output fasta file")
    parser.add_argument("--task", help="task to perform: multi2single or cma2fasta", default="multi2single")

    args = parser.parse_args()
    if args.task == "cma2fasta":
        cmatoFasta(args.input, args.output)
    else:
        multi2SingleLineFasta(args.input, args.output)
    print("done")



if __name__ == "__main__":
    # Example usage of the functions above
    # cmatoFasta("input.cma", "output.fasta")
    # multi2SingleLineFasta("input.fasta", "output.fasta")
    # a3mtoFasta("input.a3m", "output.fasta")
    pass