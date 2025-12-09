#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-23 (Y/M/D) at 17:47
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""
 
import json
import argparse

def getAF3GeneratedMSA(inputJson, outputfile, MSAformat="fasta", nolowcase=False, nogap=False):
    """
    Converts a single-line MSA string to a FASTA file format.
    Args:
        inputJson (str): Path to the AF3 JSON file.
        outputfile (str): Path to output FASTA file.
        MSAformat (str): Output format (default: fasta).
        nolowcase (bool): Remove lowercase letters if True.
        nogap (bool): Remove gap characters (-) if True.
    Returns:
        list: List of output sequences.
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
                # Remove lowercase if flag is set
                if nolowcase:
                    line = ''.join(c for c in line if not c.islower())
                # Remove gaps if flag is set
                if nogap:
                    line = line.replace("-", "")
                print(len(line))
            fp.write("%s\n" % line)
            MSAOutput.append(line)
        print('Done')
    return MSAOutput



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process AF3 generated MSA JSON to FASTA")
    parser.add_argument("inputJson", help="Input JSON file containing MSA data")
    parser.add_argument("outputfile", help="Output FASTA file")
    parser.add_argument("--MSAformat", default="fasta", help="MSA format, default is fasta")
    parser.add_argument("--nolowcase", action='store_true', help="Remove lowercase letters (inserts) from sequence")
    parser.add_argument("--nogap", action='store_true', help="Remove gaps (-) from the sequence")
    
    args = parser.parse_args()

    getAF3GeneratedMSA(args.inputJson, args.outputfile, MSAformat=args.MSAformat, 
                       nolowcase=args.nolowcase, nogap=args.nogap)






# Example usages:

# # Default: keep lowercase and gaps
# python3 outputAF3.py input.json output.fasta

# # Remove lowercase only
# python3 outputAF3.py input.json output.fasta --nolowcase

# # Remove gaps only
# python3 outputAF3.py input.json output.fasta --nogap

# # Remove both lowercase and gaps
# python3 outputAF3.py input.json output.fasta --nolowcase --nogap
