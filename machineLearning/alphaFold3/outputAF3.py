#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-23 (Y/M/D) at 17:47
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""
 
import json

def getAF3GeneratedMSA(inputJson, outputfile, no_insert = True):
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
                if no_insert:
                    line = ''.join(c for c in line if not c.islower())
                print(len(line))
            fp.write("%s\n" % line)
            MSAOutput.append(line)
        print('Done')
    return MSAOutput


path = "/mnt/mydrive/alpfaFold3/dynamics/DCLK/DCLK_KD_tail/dclk3_dclkmsa_atemp"
inputJson = path + "/dclk3_dclkmsa_atemp_data.json"
outputfile = path + "/dclk3_dclkmsa_atemp_dataMSA.fasta"
getAF3GeneratedMSA(inputJson, outputfile, no_insert=True)