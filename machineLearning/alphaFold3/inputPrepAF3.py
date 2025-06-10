#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AF3_inputPrep.py
Description: Script for preparing input data for AlphaFold3.
""" 
import json
from typing import List

def generate_af3_input_json(name: str, sequence: str, PTMs: List[dict] ,model_seeds: List[int],
                            unpaired_msa: str, paired_msa: str,templates: List[str],
                            output_file: str):
    """
    Generates an AlphaFold3 input JSON file.
    Args:
        name (str): AF3 job Name.
        sequence (str): Amino acid sequence of the protein.
        model_seeds (List[int]): List of model seeds.

        All foloowing three fields are required together if want use any of them.
        unpaired_msa (str): Unpaired MSA string in single line format.
                           "" ; MSA-free mode.
        paired_msa (str): Paired MSA string. This is for complex proteins.
        templates (List[str]):  [mmcif, id]; custom template,
                                null or ommition of this field; AF3-searched template,
                                []; template free

        output_file (str): Path to save the output JSON file.
    """
    if PTMs is None:
        PTMs = []
    else:
        PTMs = [{"ptmType": ptm["type"], "ptmPosition": ptm["position"]} for ptm in PTMs]

    af3_input = {
        "dialect": "alphafold3",
        "version": 1,
        "name": name,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": sequence,
                    "modifications": PTMs,
                    "unpairedMsa": unpaired_msa,
                    "pairedMsa": paired_msa,
                    "templates": templates,
                }
            }
        ],
        "modelSeeds": model_seeds
    }

    with open(output_file, "w") as f:
        json.dump(af3_input, f, indent=2)
    print(f" AlphaFold 3 input JSON saved to: {output_file}")



def customMSA (file):
    """
    Reads a MSA file and processes it to create a single-line MSA string.
    Args:
        file (str): Path to the custom MSA file.
    Returns: (str, str)
        str: Processed sequence with X as a gap.
        str: Processed MSA string with sequences concatenated.
    """
    custom_msa = open(file, "r").readlines()
    oneline_msa = ""
    for i,l in enumerate(custom_msa):
        if len(l) > 1:
            if l[0] != ">":
                no_lower_seq = ''.join(c for c in l if not c.islower())
                if i == 1:
                    sequence = no_lower_seq.replace("-", "X")
                    count = 0
                    for c in sequence:
                        if c == "X":
                            count += 1
                    oneline_msa = oneline_msa + sequence
                    sequence = sequence[:-1] # remove last character(a newline character)
                    print("First and query sequence: ", sequence)
                    print("Query sequence length: ", len(sequence))
                    print("Number of gaps/deletion in query sequence: ", count)
                else:
                    oneline_msa = oneline_msa + no_lower_seq
                    print(len(no_lower_seq[:-1]))# remove last character(the newline character)

            else:
                oneline_msa = oneline_msa + l
    return sequence, oneline_msa



def seqModification(file, modRemove = ["lowercase"]):
    """
    Modifies the sequence by replacing certain characters with 'X'.
    Args:
        sequence (str): Input sequence.
    Returns:
        str: Modified sequence.
    """
    sequences = open(file, "r").readlines()
    modified_sequece = []
    for i, l in enumerate(sequences):
        if len(l) > 1:
            if l[0] != ">":
                # Remove unwanted characters
                for char in modRemove:
                    if char == "lowercase":
                        l = ''.join(c for c in l if not c.islower())
                    else:
                        l = l.replace(char, "")
                modified_sequece.append(l)
            else:
                modified_sequece.append(l)
    
    with open(file[:-6] + "Mod.fasta", "w") as f:
        for seq in modified_sequece:
            f.write(seq)
    
    print("Modified sequence saved to: ", file[:-6] + "Mod.fasta")
    return




def AF3InputMSAtoFatsa(inputfile, outputfile):
    """
    Converts a single-line MSA string to a FASTA file format.
    Args:
        file (str): Path to the MSA file.
    Returns:
        str: Path to the output FASTA file.
    """
    # Open and read the JSON file
    with open(inputfile, 'r') as file:
        jdata = json.load(file)
    seqjdata = jdata["sequences"][0]["protein"]
    unpairedMsa = seqjdata["unpairedMsa"].split("\n")
    # open file in write mode
    with open(outputfile, 'w') as fp:
        for line in unpairedMsa:
            # write each item on a new line
            fp.write("%s\n" % line)
        print('Done')
    return unpairedMsa

