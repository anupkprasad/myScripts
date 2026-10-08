#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-25 (Y/M/D) at 11:11
AF3_inputPrep.py
Description: Script for preparing inputs for AlphaFold3.
more information: https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
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
        #Phosphoserine → "SEP", Phosphothreonine → "TPO", Phosphotyrosine → "PTR"

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
                },
                {
                    "ligand": {
                        "id": ["I"],
                        "ccdCodes": ["MG"]
                }
                },
                {
                    "ligand": {
                        "id": ["J"],
                        "ccdCodes": ["MG"]
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
    PTMs: List = None,
    unpaired_msa: str = None,
    paired_msa: str = None,
    templates: List[str] = None,
    chainID: List[str] = ["A", "B"]
):
    """
    Generates an AlphaFold3 input JSON file for multiple chains.

    PTMs: per-chain modifications. Either:
      - list of lists aligned with chainID, e.g.
        [[], [{"type": "SEP", "position": 11}]]  # phospho only on chain B
      - flat list applied to every chain (legacy), e.g.
        [{"type": "SEP", "position": 11}]
    Positions are 1-based (residue 1 = first amino acid). CCD codes:
      Phosphoserine SEP, Phosphothreonine TPO, Phosphotyrosine PTR.
    """
    def _format_ptms(ptm_list):
        if not ptm_list:
            return []
        return [{"ptmType": ptm["type"], "ptmPosition": ptm["position"]} for ptm in ptm_list]

    if PTMs is None:
        PTMs = [[] for _ in sequence]
    elif PTMs and isinstance(PTMs[0], dict):
        # Legacy: same PTMs on all chains
        PTMs = [_format_ptms(PTMs) for _ in sequence]
    else:
        PTMs = [_format_ptms(chain_ptms) for chain_ptms in PTMs]

    if len(PTMs) != len(sequence):
        raise ValueError(f"PTMs length ({len(PTMs)}) must match sequence length ({len(sequence)})")

    sequences_block = []
    for i, (cid, seq) in enumerate(zip(chainID, sequence)):
        protein = {
            "id": cid,
            "sequence": seq,
            "modifications": PTMs[i],
            "templates": templates,
        }
        if i == 0:
            protein["unpairedMsa"] = unpaired_msa
            protein["pairedMsa"] = paired_msa
        sequences_block.append({"protein": protein})

    af3_input = {
        "dialect": "alphafold3",
        "version": 1,
        "name": name,
        "sequences": sequences_block,
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



USAGE = """
Usage: inputPrepAF3.py
======================

This file is a function module. Import functions and call them from Python.

Run:
    python3 your_script.py

1) Single-chain AF3 JSON
------------------------
from inputPrepAF3 import generateAF3InputJson

generateAF3InputJson(
    name="job_single_chain",
    model_seeds=[0, 100, 200],
    sequence="MSSIGTGYDLSASTF",
    output_file="single_chain.json",
    PTMs=[{"type": "SEP", "position": 11}],  # 1-based index
    unpaired_msa=None,
    paired_msa=None,
    templates=None,
    chainID="A",
)

2) Single-chain with ligand
---------------------------
generateAF3InputJson(
    name="job_with_ligand",
    model_seeds=[0, 100, 200],
    sequence="MSSIGTGYDLSASTF",
    output_file="single_chain_ligand.json",
    PTMs=[],
    ligand="ATP",  # CCD code
)

3) Multi-chain AF3 JSON (different PTMs per chain)
---------------------------------------------------
from inputPrepAF3 import generateMultichainAF3InputJson

generateMultichainAF3InputJson(
    name="job_multichain",
    model_seeds=[0, 100, 200],
    sequence=[
        "AAAAA",
        "MSSIGTGYDLSASTF",
    ],
    output_file="multichain.json",
    chainID=["A", "B"],
    PTMs=[
        [],  # chain A
        [{"type": "SEP", "position": 11}],  # chain B
    ],
    unpaired_msa=None,
    paired_msa=None,
    templates=None,
)

4) Multi-chain legacy PTM mode (same PTM applied to all chains)
----------------------------------------------------------------
generateMultichainAF3InputJson(
    name="job_multichain_legacy",
    model_seeds=[0, 100, 200],
    sequence=["AAAAA", "SSSSS"],
    output_file="multichain_legacy.json",
    PTMs=[{"type": "SEP", "position": 2}],
)

5) Convert custom MSA file to one-line AF3 MSA string
------------------------------------------------------
from inputPrepAF3 import customMSA

sequence, oneline_msa = customMSA("my_alignment.a3m")

# Then pass these values to AF3 JSON generator:
generateAF3InputJson(
    name="job_with_custom_msa",
    model_seeds=[0],
    sequence=sequence,
    output_file="msa_input.json",
    PTMs=[],
    unpaired_msa=oneline_msa,
    paired_msa="",
    templates=[],
)

6) Clean FASTA/MSA sequence characters
--------------------------------------
from inputPrepAF3 import seqModification

seqModification("input.fasta", modRemove=["lowercase", "-"])

Notes:
- PTM position is 1-based (first residue = 1).
- PTM type uses CCD code, e.g. SEP, TPO, PTR.
- For each modified site, residue identity should match PTM chemistry.
- In multi-chain mode, PTMs should usually be provided per chain as list-of-lists.
"""

