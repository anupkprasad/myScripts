#!/usr/bin/env python3

"""
Pipeline:
1. Parse domtblout file to extract PFAM hits.
2. Apply E-value cutoff to select confident family hits.
3. Map each sequence → best family.
4. Convert sequences to embeddings using k-mer frequency.
5. UMAP visualization.
"""

import pandas as pd
from Bio import SeqIO
import numpy as np
from sklearn.preprocessing import normalize
import umap
import matplotlib.pyplot as plt

# -----------------------------
# 1. Parse domtblout
# -----------------------------
def parse_domtblout(file, evalue_cutoff=1e-5):
    rows = []
    with open(file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            if len(fields) < 22:
                continue

            family = fields[0]
            accession = fields[1]
            query = fields[3]
            evalue = float(fields[6])

            if evalue > evalue_cutoff:
                continue

            rows.append({
                "query": query,
                "family": family,
                "accession": accession,
                "evalue": evalue,
                "score": float(fields[7])
            })

    df = pd.DataFrame(rows)
    return df


# -----------------------------
# 2. Choose best family hit per sequence
# -----------------------------
def select_best_hits(df):
    return df.sort_values("evalue").groupby("query").first().reset_index()


# -----------------------------
# 3. Convert sequences → k-mer embedding
# -----------------------------
def kmer_embedding(seqs, k=3):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    kmer_list = [a+b+c for a in amino_acids for b in amino_acids for c in amino_acids]
    kmer_index = {k:i for i,k in enumerate(kmer_list)}

    embeddings = np.zeros((len(seqs), len(kmer_list)))

    for i, seq in enumerate(seqs):
        s = str(seq).upper()
        for j in range(len(s)-k+1):
            kmer = s[j:j+k]
            if kmer in kmer_index:
                embeddings[i, kmer_index[kmer]] += 1

    # Normalize for comparison
    return normalize(embeddings)


# -----------------------------
# 4. Make UMAP plot
# -----------------------------
def make_umap_plot(embeddings, labels, out="umap_plot.png"):
    reducer = umap.UMAP(n_neighbors=20, min_dist=0.1, metric='euclidean')
    X = reducer.fit_transform(embeddings)

    plt.figure(figsize=(10,7))
    for f in sorted(set(labels)):
        mask = labels == f
        plt.scatter(X[mask,0], X[mask,1], s=20, label=f, alpha=0.7)

    plt.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title("UMAP clustering of kinase sequences by PFAM family")
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    print(f"Saved UMAP plot → {out}")


# -----------------------------
# MAIN
# -----------------------------
def main(domtblout, fasta):
    print("Parsing domtblout...")
    df = parse_domtblout(domtblout)

    print("Selecting best family hit per sequence...")
    best = select_best_hits(df)

    print("Loading FASTA...")
    seq_dict = {rec.id: rec.seq for rec in SeqIO.parse(fasta, "fasta")}

    print("Matching sequences to family annotations...")
    best["sequence"] = best["query"].map(seq_dict)

    best = best.dropna()

    print("Building k-mer embeddings...")
    embeddings = kmer_embedding(list(best["sequence"]))

    print("Generating UMAP visualization...")
    make_umap_plot(embeddings, best["family"].values)

    best.to_csv("family_assignments.tsv", sep="\t", index=False)
    print("Saved table → family_assignments.tsv")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--domtblout", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    args = parser.parse_args()

    main(args.domtblout, args.fasta)
