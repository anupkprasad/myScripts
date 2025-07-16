#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-20 (Y/M/D) at 15:44
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""
import os
from dash import Dash, html
import dash_bio as dashbio

def viewAlignment(fasta_file: str, port: int = 8050):
    """
    Launch an interactive MSA alignment viewer using Dash Bio.
        fasta_file (str): Path to the input FASTA file.
        port (int): Port to run the Dash server on (default 8050).
        Use multiple separate port to compare multiple MSA files.
    """
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    with open(fasta_file) as f:
        data = f.read()

    app = Dash(__name__)
    app.layout = html.Div([
        html.H3(f"MSA Viewer - {os.path.basename(fasta_file)}"),
        dashbio.AlignmentChart(
            id='alignment-viewer',
            data=data
        )
    ])

    print(f"Launching viewer at http://127.0.0.1:{port}/")
    app.run(port=port, debug=True)


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



path = "/mnt/mydrive/alpfaFold3/dynamics/EGFR/KD_tail"
multi2SingleLineFasta(input_file=path + "/PKB3_KD_tail_query.fasta", output_file= path + "/EGFR_KD_tail_query.fasta")