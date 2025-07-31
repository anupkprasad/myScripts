#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-25 (Y/M/D) at 11:11
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""
"""
AF3_inputPrep.py
Description: Script for preparing inputs for AlphaFold3.
""" 
import json
from typing import List

def generateAF3InputJson(
    name: str,
    model_seeds: List[int],
    sequence: str,
    output_file: str,
    PTMs: List[dict]= [],
    unpaired_msa: str= None,
    paired_msa: str= None,
    templates: List[str]=None,
    chainID: str = "A",
    ligand: str = None
):
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
    if PTMs != []:
        PTMs = [{"ptmType": ptm["type"], "ptmPosition": ptm["position"]} for ptm in PTMs]

    if ligand:
        af3_input = {
            "dialect": "alphafold3",
            "version": 1,
            "name": name,
            "sequences": [
                {
                    "protein": {
                        "id": chainID,
                        "sequence": sequence,
                        "modifications": PTMs,
                        "unpairedMsa": unpaired_msa,
                        "pairedMsa": paired_msa,
                        "templates": templates,
                    }
                },
                {
                    "ligand": {
                        "id": ["L"],  # Example ligand ID
                        "ccdCodes": [ligand]  # Example CCD code for ligands
                    }
                }
            ],
            "modelSeeds": model_seeds
        }
    else:
        af3_input = {
            "dialect": "alphafold3",
            "version": 1,
            "name": name,
            "sequences": [
                {
                    "protein": {
                        "id": chainID,
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




def generateMultichainAF3InputJson(
    name: str,
    model_seeds: List[int],
    sequence: List[str],
    output_file: str,
    PTMs: List[dict]= [],
    unpaired_msa: str= None,
    paired_msa: str= None,
    templates: List[str]=None,
    chainID: List[str] = ["A", "B"]
):
    """
    Generates an AlphaFold3 input JSON file for multiple chains.
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
                    "id": chainID[0],
                    "sequence": sequence[0],
                    "modifications": PTMs,
                    "unpairedMsa": unpaired_msa,
                    "pairedMsa": paired_msa,
                    "templates": templates,
                }
            },
            {
                "protein": {
                    "id": chainID[1],
                    "sequence": sequence[1],
                    "modifications": PTMs,
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

