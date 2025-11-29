#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-11-20 (Y/M/D) at 15:57
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

"""
clean_msa_with_report.py

Clean an MSA (FASTA) for coevolution analysis and produce QC plots & summary.

Usage:
    python clean_msa_with_report.py input_msa.fasta cleaned_msa.fasta --outdir results \
        [--remove_paralogs] [--min_cov 0.7] [--max_col_gap 0.3] [--theta 0.8]
        [--use_jackhmmer /usr/bin/jackhmmer --query query.fasta]

Notes:
 - Default parameters follow common DCA practice:
     min_cov = 0.70, max_col_gap = 0.30, theta = 0.80
 - If --use_jackhmmer is provided, the script will attempt to run jackhmmer to realign/expand
   around the provided --query file. This is optional; if you do not want jackhmmer, skip it.
"""

import argparse
import subprocess
import json
import os
import re
from collections import defaultdict
from math import log2
from datetime import datetime

from Bio import AlignIO, SeqIO
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

# -----------------------------
# Helper functions
# -----------------------------

def get_species_id(header):
    """
    Try to extract a species identifier from header.
    Heuristic; adapt if your headers have different format.
    """
    h = header.strip()
    # Uniprot-like: >sp|P12345|NAME OS=Homo sapiens ...
    m = re.search(r"OS=([A-Za-z_ ]+)", h)
    if m:
        return m.group(1).strip().replace(' ', '_')
    # Try last token uppercase (HUMAN, YEAST)
    m = re.search(r"([A-Z]{3,6})$", h)
    if m:
        return m.group(1)
    # Fall back to first token
    return h.split()[0]

def read_fasta_msa(path):
    """Read fasta MSA and return list of SeqRecord"""
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) == 0:
        raise ValueError(f"No sequences found in {path}")
    return records

def write_fasta_msa(records, path):
    SeqIO.write(records, path, "fasta")

def seq_to_strs(records):
    return [str(r.seq) for r in records]

def coverage_fraction(seq):
    L = len(seq)
    return 1.0 - (seq.count('-') / L)

def filter_by_coverage(records, min_cov=0.7):
    L = len(records[0].seq)
    kept = []
    for r in records:
        cov = coverage_fraction(str(r.seq))
        if cov >= min_cov:
            kept.append(r)
    return kept

def remove_paralogs(records):
    species_map = defaultdict(list)
    for r in records:
        sp = get_species_id(r.id)
        species_map[sp].append(r)
    reps = []
    for sp, recs in species_map.items():
        # choose sequence with best coverage (non-gaps)
        rep = max(recs, key=lambda x: len(str(x.seq)) - str(x.seq).count('-'))
        reps.append(rep)
    return reps

def trim_gappy_columns(records, max_col_gap=0.3):
    arr = np.array([list(str(r.seq)) for r in records])
    gap_frac = (arr == '-').mean(axis=0)
    keep_mask = gap_frac <= max_col_gap
    if keep_mask.sum() == 0:
        raise ValueError("No columns remain after trimming; relax max_col_gap.")
    trimmed = ["".join(arr[i, keep_mask]) for i in range(arr.shape[0])]
    for i, r in enumerate(records):
        r.seq = trimmed[i]
    return records, gap_frac  # return raw gap_frac for plotting original profile

def seq_identity(a, b):
    # Expect same length
    matches = sum(1 for x,y in zip(a,b) if x == y)
    return matches / len(a)

def compute_reweighting_weights(records, theta=0.8):
    seqs = seq_to_strs(records)
    N = len(seqs)
    weights = np.zeros(N, dtype=float)
    for i in range(N):
        cnt = 0
        si = seqs[i]
        for j in range(N):
            sj = seqs[j]
            if seq_identity(si, sj) >= theta:
                cnt += 1
        weights[i] = 1.0/cnt
    Neff = weights.sum()
    return weights, Neff

def shannon_entropy_column(column):
    # column: list/array of residues, include gaps
    # compute freq excluding gaps? Usually include gaps as a state; here we exclude gaps for entropy calc.
    counts = {}
    total = 0
    for a in column:
        if a == '-':
            continue
        counts[a] = counts.get(a, 0) + 1
        total += 1
    if total == 0:
        return 0.0
    ent = 0.0
    for v in counts.values():
        p = v/total
        ent -= p * log2(p)
    return ent

def column_entropies(records):
    arr = np.array([list(str(r.seq)) for r in records])
    ent = [shannon_entropy_column(arr[:,i]) for i in range(arr.shape[1])]
    return np.array(ent)

def make_seqid_matrix(records, out_png):
    seqs = seq_to_strs(records)
    N = len(seqs)
    mat = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(i, N):
            v = seq_identity(seqs[i], seqs[j])
            mat[i,j] = v
            mat[j,i] = v
    fig, ax = plt.subplots(figsize=(6,6))
    im = ax.imshow(mat, vmin=0, vmax=1)
    ax.set_title("Sequence identity matrix")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xlabel("sequence index")
    ax.set_ylabel("sequence index")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return mat

def plot_gap_profile(gap_frac, out_png):
    fig, ax = plt.subplots(figsize=(8,3))
    ax.plot(np.arange(1, len(gap_frac)+1), gap_frac)
    ax.set_xlabel("alignment column (1-based)")
    ax.set_ylabel("gap fraction")
    ax.set_title("Gap fraction per column (original alignment)")
    ax.set_xlim(1, len(gap_frac))
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

def plot_entropy(entropies, out_png):
    fig, ax = plt.subplots(figsize=(8,3))
    ax.plot(np.arange(1, len(entropies)+1), entropies)
    ax.set_xlabel("alignment column (1-based)")
    ax.set_ylabel("Shannon entropy (bits)")
    ax.set_title("Per-column Shannon entropy (no gaps)")
    ax.set_xlim(1, len(entropies))
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

def plot_seqid_histogram(mat, out_png):
    # use upper triangle excluding diagonal
    triu = mat[np.triu_indices_from(mat, k=1)]
    fig, ax = plt.subplots(figsize=(5,3))
    ax.hist(triu, bins=40)
    ax.set_xlabel("pairwise sequence identity")
    ax.set_ylabel("count")
    ax.set_title("Pairwise sequence identity distribution")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

def run_jackhmmer(query_fasta, database_fasta, output_sto, iterations, jackhmmer_path='/usr/bin/jackhmmer'):
    cmd = [
        jackhmmer_path,
        '-N', str(iterations),
        '-A', output_sto,
        query_fasta,
        database_fasta
    ]
    print("Running jackhmmer: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    # convert sto to fasta using Biopython AlignIO
    aln = AlignIO.read(output_sto, "stockholm")
    fasta_path = output_sto + ".fasta"
    AlignIO.write(aln, fasta_path, "fasta")
    return fasta_path

# -----------------------------
# Main pipeline
# -----------------------------

def pipeline(args):
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    start_time = datetime.utcnow().isoformat() + "Z"
    report = {
        "input_msa": args.input,
        "cleaned_msa": args.output,
        "parameters": {
            "min_cov": args.min_cov,
            "max_col_gap": args.max_col_gap,
            "theta": args.theta,
            "remove_paralogs": args.remove_paralogs,
            "used_jackhmmer": bool(args.use_jackhmmer)
        },
        "start_time": start_time
    }

    # Optionally run jackhmmer if requested (realign/expand around query)
    working_msa = args.input
    if args.use_jackhmmer:
        if not args.query:
            raise ValueError("When using --use_jackhmmer you must provide --query query.fasta")
        sto_out = os.path.join(outdir, "jackhmmer_out.sto")
        fasta_from_sto = run_jackhmmer(args.query, working_msa, sto_out, iterations=2, jackhmmer_path=args.use_jackhmmer)
        working_msa = fasta_from_sto
        report["jackhmmer_sto"] = sto_out
        report["working_msa_after_jackhmmer"] = working_msa

    # Read MSA
    records = read_fasta_msa(working_msa)
    report["n_initial_sequences"] = len(records)
    L0 = len(records[0].seq)
    report["initial_alignment_length"] = L0

    # 1) Filter by coverage
    records = filter_by_coverage(records, min_cov=args.min_cov)
    report["n_after_coverage"] = len(records)
    if len(records) == 0:
        raise ValueError("No sequences survive coverage filter; relax --min_cov")

    # 2) Optional paralog removal
    if args.remove_paralogs:
        records = remove_paralogs(records)
    report["n_after_paralog_filter"] = len(records)

    # 3) Trim highly gappy columns (we still keep original gap profile for plotting)
    # compute original gap profile for plotting
    arr0 = np.array([list(str(r.seq)) for r in records])
    orig_gap_frac = (arr0 == '-').mean(axis=0)

    records, raw_gap_frac = trim_gappy_columns(records, max_col_gap=args.max_col_gap)

    # 4) Compute weights & Neff
    weights, Neff = compute_reweighting_weights(records, theta=args.theta)
    L = len(records[0].seq)
    report["n_final_sequences"] = len(records)
    report["final_alignment_length"] = L
    report["Neff"] = float(Neff)
    report["Neff_over_L"] = float(Neff / L)

    # 5) Entropy profile
    ent = column_entropies(records)

    # 6) Sequence identity matrix and histogram
    seqid_mat = make_seqid_matrix(records, os.path.join(outdir, "seqid_matrix.png"))
    plot_seqid_histogram(seqid_mat, os.path.join(outdir, "seqid_hist.png"))

    # 7) Gap profile plot (original)
    plot_gap_profile(orig_gap_frac, os.path.join(outdir, "gap_profile_original.png"))

    # 8) Entropy plot
    plot_entropy(ent, os.path.join(outdir, "entropy_profile.png"))

    # 9) Save cleaned msa
    write_fasta_msa(records, args.output)

    # 10) Save report JSON
    end_time = datetime.utcnow().isoformat() + "Z"
    report["end_time"] = end_time
    report["plots"] = {
        "seqid_matrix": "seqid_matrix.png",
        "seqid_histogram": "seqid_hist.png",
        "gap_profile_original": "gap_profile_original.png",
        "entropy_profile": "entropy_profile.png"
    }
    report_path = os.path.join(outdir, "cleaning_report.json")
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)

    # 11) Write a short human-readable summary
    summary_txt = os.path.join(outdir, "summary.txt")
    with open(summary_txt, "w") as fh:
        fh.write(f"MSA cleaning report\n")
        fh.write(f"Input MSA: {args.input}\n")
        fh.write(f"Cleaned MSA: {args.output}\n")
        fh.write(f"Initial sequences: {report['n_initial_sequences']}\n")
        fh.write(f"After coverage filter: {report['n_after_coverage']}\n")
        fh.write(f"After paralog filter: {report['n_after_paralog_filter']}\n")
        fh.write(f"Final sequences: {report['n_final_sequences']}\n")
        fh.write(f"Initial alignment length: {report['initial_alignment_length']}\n")
        fh.write(f"Final alignment length: {report['final_alignment_length']}\n")
        fh.write(f"Neff: {report['Neff']:.2f}\n")
        fh.write(f"Neff/L: {report['Neff_over_L']:.3f}\n")
        fh.write(f"Parameters: min_cov={args.min_cov}, max_col_gap={args.max_col_gap}, theta={args.theta}\n")
        fh.write("\nPlots saved to directory: " + outdir + "\n")
    print("Done. Results in:", outdir)
    print("Summary written to:", summary_txt)
    print("JSON report written to:", report_path)

# -----------------------------
# CLI
# -----------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean MSA FASTA and produce QC plots & report")
    parser.add_argument("input", help="input MSA in FASTA format")
    parser.add_argument("output", help="output cleaned MSA FASTA")
    parser.add_argument("--outdir", default="msa_cleaning_results", help="output directory for plots and report")
    parser.add_argument("--remove_paralogs", action="store_true", help="remove paralogs (one representative per species)")
    parser.add_argument("--min_cov", type=float, default=0.70, help="minimum sequence coverage (fraction non-gaps)")
    parser.add_argument("--max_col_gap", type=float, default=0.30, help="maximum allowed gap fraction per column")
    parser.add_argument("--theta", type=float, default=0.80, help="sequence identity threshold for reweighting")
    parser.add_argument("--use_jackhmmer", type=str, default="", help="path to jackhmmer executable to run (optional)")
    parser.add_argument("--query", type=str, default="", help="query FASTA to use with jackhmmer (required if using jackhmmer)")
    args = parser.parse_args()
    # Normalize empty string to False
    if args.use_jackhmmer == "":
        args.use_jackhmmer = None
    pipeline(args)

