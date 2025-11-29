#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-23 (Y/M/D) at 17:47
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""
 
import json



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




def getAF3GeneratedMSA(inputJson, outputfile, MSAformat = "fasta", insert=True):
    """
    Converts a single-line MSA string to a FASTA file format.
    Args:
        file (str): Path to the MSA file.
    Returns:
        str: Path to the output FASTA file.
    """
    # Open and read the JSON file
    with open(inputJson, 'r') as file:
        jdata = json.load(file)
    seqjdata = jdata["sequences"][0]["protein"]
    unpairedMsa = seqjdata["unpairedMsa"].split("\n")
    # open file in write mode
    MSAOutput = []
    with open(outputfile, 'w') as fp:
        for line in unpairedMsa:
            # write each item on a new line
            if not line.startswith('>'):
                if MSAformat == "fasta" and insert:
                    line = line
                else:
                    line = ''.join(c for c in line if not c.islower())
                print(len(line))
            fp.write("%s\n" % line)
            MSAOutput.append(line)
        print('Done')
    return MSAOutput



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


path = "/mnt/mydrive/labTool/darc/darc_on_AF3_ensemble/dclk1_full_length"
inputJson = path + "/dclk1fulllength_amsa_atemp_data.json"
outputfile = path + "/dclk1_new.fasta"
getAF3GeneratedMSA(inputJson, outputfile, MSAformat="fasta", insert=True)

#a3mtoFasta(path + "/dcx2.a3m", path + "/dcx2.fasta")





"""
Convert .cma alignment file to FASTA.

Usage:
    python cma_to_fasta_robust.py input.cma output.fasta

This parser is tolerant of the CMA structure shown in your example:
  - headers start with '>'
  - sequence content appears inside a following {...} block, often between '{()' and '()}'.
  - the script strips any non-letter characters from the extracted sequence.
"""

import sys
import re
from textwrap import fill

def extract_sequences_from_cma(text):
    """
    Returns a list of (header, sequence) tuples extracted from the CMA text.
    The algorithm:
      - find each header line '>' and use that as the header
      - for each header, locate the next '{' after the header and try to find a sequence
        end marker '()}' (preferred). If not found, fall back to the nearest '}'.
      - extract the substring inside the braces and strip any non-letter characters.
    """
    results = []
    # find all header positions
    for m in re.finditer(r'>([^\n\r]+)', text):
        header = m.group(1).strip()
        start_pos = m.end()

        # find the next '{' that starts the sequence block
        open_brace = text.find('{', start_pos)
        if open_brace == -1:
            # nothing found for this header
            results.append((header, ""))
            continue

        # prefer sequence end pattern '()}' as seen in examples
        close_marker = text.find('()}', open_brace)
        if close_marker != -1:
            seq_block = text[open_brace+1:close_marker]   # content inside { ... before ()}
        else:
            # fallback: find the next '}' after open_brace
            close_brace = text.find('}', open_brace)
            if close_brace == -1:
                seq_block = text[open_brace+1:]
            else:
                seq_block = text[open_brace+1:close_brace]

        # Remove any bracketed annotations like "(...)" that may surround the sequence marker,
        # then remove any non-letter characters. This keeps only A-Z/a-z and '-' if desired.
        # Here we keep only letters A-Z.
        seq_letters = re.sub(r'[^A-Za-z]', '', seq_block)

        # uppercase
        seq_letters = seq_letters.upper()

        results.append((header, seq_letters))

    return results

def write_fasta(seq_tuples, outpath, width=60, header_token_first=False):
    """
    Write (header, seq) tuples to FASTA file.
    If header_token_first is True, the FASTA ID will be the first token of the header
    (useful if header lines are long).
    """
    with open(outpath, 'w') as out:
        for header, seq in seq_tuples:
            if header_token_first:
                fasta_header = header.split()[0]
            else:
                # sanitize header to avoid '>' characters (they shouldn't be there)
                fasta_header = header.replace('>', '').strip()
            out.write(f'>{fasta_header}\n')
            if seq:
                # wrap lines
                out.write(fill(seq, width=width) + '\n')
            else:
                out.write('>')  # write an empty line (will be caught below)
    return

def cmatoFasta(infile, outfile):
 
    with open(infile, 'r', encoding='utf-8', errors='ignore') as fh:
        text = fh.read()

    seqs = extract_sequences_from_cma(text)

    # warn about empty sequences and optionally filter them out
    n_empty = sum(1 for h,s in seqs if not s)
    if n_empty:
        print(f"Warning: {n_empty} sequence(s) had no letters extracted. They will be written as empty entries.")
        # You might want to remove entries with empty sequences:
        # seqs = [(h,s) for (h,s) in seqs if s]

    write_fasta(seqs, outfile, width=60, header_token_first=False)
    print(f"Wrote {len(seqs)} entries to {outfile}. {n_empty} empty sequences.")





cmatoFasta("/mnt/mydrive/others/rezwan/seq_rhl.cma",
             "/mnt/mydrive/others/rezwan/seq_rhl.fasta")


multi2SingleLineFasta("/mnt/mydrive/others/rezwan/seq_rhl.fasta",
                      "/mnt/mydrive/others/rezwan/seq_rhl_singleline.fasta")