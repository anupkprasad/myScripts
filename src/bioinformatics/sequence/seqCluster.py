#!/usr/bin/env python3
"""
seqCluster — cluster / embed sequences (HMMER domtblout → k-mer → UMAP).

When to use
    You have a FASTA of sequences and a HMMER ``domtblout`` labeling families;
    make a 2-D UMAP colored by best PFAM/family hit.

Public API
    ParseDomtblout(file, evalueCutoff=1e-5) → DataFrame
    SelectBestHits(df) → DataFrame
    KmerEmbedding(seqs, k=3) → ndarray
    PlotUmap(embeddings, labels, outFile=...)
    ClusterByDomtblout(domtblout, fasta, outDir=...)

Requires (optional extras)
    pandas, scikit-learn, umap-learn

Example
    from bioinformatics.sequence import ClusterByDomtblout
    ClusterByDomtblout("hits.domtblout", "seqs.fasta", outDir="cluster_out")

CLI
    python -m bioinformatics.sequence.seqCluster -domtblout hits.domtblout -fasta seqs.fasta
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np

PathLike = Union[str, Path]


def _requirePandas():
    try:
        import pandas as pd
    except ImportError as exc:
        raise SystemExit("seqCluster needs pandas") from exc
    return pd


def ParseDomtblout(file: PathLike, evalueCutoff: float = 1e-5):
    """Parse HMMER domtblout → DataFrame (query, family, accession, evalue, score)."""
    pd = _requirePandas()
    path = Path(file).expanduser().resolve()
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            if len(fields) < 22:
                continue
            evalue = float(fields[6])
            if evalue > evalueCutoff:
                continue
            rows.append(
                {
                    "query": fields[3],
                    "family": fields[0],
                    "accession": fields[1],
                    "evalue": evalue,
                    "score": float(fields[7]),
                }
            )
    return pd.DataFrame(rows)


def SelectBestHits(df):
    """Best (lowest E-value) family hit per query."""
    if df is None or len(df) == 0:
        return df
    return df.sort_values("evalue").groupby("query", as_index=False).first()


def KmerEmbedding(seqs: Sequence, k: int = 3) -> np.ndarray:
    """Normalized k-mer frequency embedding over standard amino acids."""
    try:
        from sklearn.preprocessing import normalize
    except ImportError as exc:
        raise SystemExit("seqCluster needs scikit-learn") from exc

    aminoAcids = "ACDEFGHIKLMNPQRSTVWY"
    if k != 3:
        raise ValueError("only k=3 is implemented")
    kmerList = [a + b + c for a in aminoAcids for b in aminoAcids for c in aminoAcids]
    kmerIndex = {kmer: i for i, kmer in enumerate(kmerList)}
    embeddings = np.zeros((len(seqs), len(kmerList)), dtype=float)
    for i, seq in enumerate(seqs):
        s = str(seq).upper()
        for j in range(len(s) - k + 1):
            kmer = s[j : j + k]
            idx = kmerIndex.get(kmer)
            if idx is not None:
                embeddings[i, idx] += 1
    return normalize(embeddings)


def PlotUmap(
    embeddings: np.ndarray,
    labels: Sequence,
    outFile: PathLike = "umap_plot.png",
    title: str = "UMAP of sequences by family",
):
    try:
        import umap
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "PlotUmap needs umap-learn and matplotlib:\n"
            "  pip install umap-learn"
        ) from exc

    reducer = umap.UMAP(n_neighbors=20, min_dist=0.1, metric="euclidean")
    xy = reducer.fit_transform(embeddings)
    outPath = Path(outFile).expanduser().resolve()
    outPath.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 7))
    labels = np.asarray(labels)
    for family in sorted(set(labels)):
        mask = labels == family
        plt.scatter(xy[mask, 0], xy[mask, 1], s=20, label=family, alpha=0.7)
    plt.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outPath, dpi=300)
    plt.close()
    print(f"seqCluster: wrote {outPath}")
    return outPath


def ClusterByDomtblout(
    domtblout: PathLike,
    fasta: PathLike,
    outDir: PathLike = ".",
    evalueCutoff: float = 1e-5,
) -> Path:
    """End-to-end: parse hits → embed → UMAP → write assignments TSV."""
    from Bio import SeqIO

    outPath = Path(outDir).expanduser().resolve()
    outPath.mkdir(parents=True, exist_ok=True)

    print("seqCluster: parsing domtblout…")
    df = ParseDomtblout(domtblout, evalueCutoff=evalueCutoff)
    print("seqCluster: selecting best hits…")
    best = SelectBestHits(df)
    seqDict = {rec.id: rec.seq for rec in SeqIO.parse(str(Path(fasta).expanduser()), "fasta")}
    best["sequence"] = best["query"].map(seqDict)
    best = best.dropna(subset=["sequence"])
    if len(best) == 0:
        raise ValueError("no queries matched between domtblout and FASTA")

    print("seqCluster: k-mer embeddings…")
    embeddings = KmerEmbedding(list(best["sequence"]))
    print("seqCluster: UMAP…")
    PlotUmap(embeddings, best["family"].values, outFile=outPath / "umap_plot.png")
    tsvPath = outPath / "family_assignments.tsv"
    best.drop(columns=["sequence"]).to_csv(tsvPath, sep="\t", index=False)
    print(f"seqCluster: wrote {tsvPath}")
    return tsvPath


def main(argv=None):
    parser = argparse.ArgumentParser(description="Cluster sequences via domtblout + UMAP")
    parser.add_argument("-domtblout", required=True)
    parser.add_argument("-fasta", required=True)
    parser.add_argument("-outdir", default=".")
    parser.add_argument("-evalue", type=float, default=1e-5)
    args = parser.parse_args(argv)
    ClusterByDomtblout(
        args.domtblout,
        args.fasta,
        outDir=args.outdir,
        evalueCutoff=args.evalue,
    )


if __name__ == "__main__":
    main()
