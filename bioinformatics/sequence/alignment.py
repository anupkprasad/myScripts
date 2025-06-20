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

