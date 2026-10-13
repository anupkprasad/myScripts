#!/usr/bin/env python3
"""
msaStats — column / sequence statistics on an MSA (foundation for conservation plots).

When to use
    After LoadMsa / CleanMsa, compute per-column conservation, entropy, gaps,
    or pairwise identity. Plot helpers write simple PNG profiles you can extend
    later for richer MSA / conservation figures.

Public API
    GapFraction(records) → (nCol,)
    ColumnEntropy(records) → (nCol,) bits (gaps excluded)
    ColumnConservation(records) → (nCol,) frequency of consensus residue
    SeqIdentityMatrix(records) → (nSeq, nSeq)
    EffectiveSeqCount(records, theta=0.8) → (weights, Neff)
    PlotColumnProfile(values, outFile, ...)  generic 1-D column plot
    PlotConservation(records, outFile)       consensus frequency vs column
    PlotEntropy(records, outFile)
    PlotGapProfile(records, outFile)
    PlotSeqIdMatrix(records, outFile)

Example
    from bioinformatics.sequence import LoadMsa, PlotConservation, ColumnConservation

    records = LoadMsa("aln.fasta")
    scores = ColumnConservation(records)
    PlotConservation(records, outFile="conservation.png")
"""

from __future__ import annotations

import argparse
from collections import Counter
from math import log2
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np

from .msaIO import LoadMsa

PathLike = Union[str, Path]


def _asCharArray(records: Sequence) -> np.ndarray:
    seqs = [str(r.seq) for r in records]
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"unequal sequence lengths: {sorted(lengths)[:5]}")
    return np.asarray([list(s) for s in seqs])


def GapFraction(records: Sequence) -> np.ndarray:
    """Per-column gap fraction (0–1)."""
    arr = _asCharArray(records)
    return (arr == "-").mean(axis=0)


def ColumnEntropy(records: Sequence) -> np.ndarray:
    """Per-column Shannon entropy in bits (gaps excluded)."""
    arr = _asCharArray(records)
    nCol = arr.shape[1]
    ent = np.zeros(nCol, dtype=float)
    for i in range(nCol):
        counts: dict = {}
        total = 0
        for a in arr[:, i]:
            if a == "-":
                continue
            counts[a] = counts.get(a, 0) + 1
            total += 1
        if total == 0:
            continue
        e = 0.0
        for v in counts.values():
            p = v / total
            e -= p * log2(p)
        ent[i] = e
    return ent


def ColumnConservation(records: Sequence) -> np.ndarray:
    """
    Per-column conservation = frequency of the most common non-gap residue (0–1).

    Good default score for a conservation strip above an MSA plot.
    """
    arr = _asCharArray(records)
    nCol = arr.shape[1]
    scores = np.zeros(nCol, dtype=float)
    for i in range(nCol):
        col = [a for a in arr[:, i] if a != "-"]
        if not col:
            continue
        scores[i] = Counter(col).most_common(1)[0][1] / len(col)
    return scores


def SeqIdentity(a: str, b: str) -> float:
    """Fraction of identical columns (including gaps as characters)."""
    if len(a) != len(b):
        raise ValueError("sequences must be equal length")
    if not a:
        return 0.0
    return sum(x == y for x, y in zip(a, b)) / len(a)


def SeqIdentityMatrix(records: Sequence) -> np.ndarray:
    """Pairwise sequence-identity matrix (nSeq × nSeq)."""
    seqs = [str(r.seq) for r in records]
    n = len(seqs)
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        mat[i, i] = 1.0
        for j in range(i + 1, n):
            v = SeqIdentity(seqs[i], seqs[j])
            mat[i, j] = v
            mat[j, i] = v
    return mat


def EffectiveSeqCount(
    records: Sequence,
    theta: float = 0.8,
) -> Tuple[np.ndarray, float]:
    """Henikoff-style reweighting: weight_i = 1 / #{j : id(i,j) ≥ theta}; Neff = Σ w."""
    seqs = [str(r.seq) for r in records]
    n = len(seqs)
    weights = np.zeros(n, dtype=float)
    for i in range(n):
        cnt = sum(1 for j in range(n) if SeqIdentity(seqs[i], seqs[j]) >= theta)
        weights[i] = 1.0 / max(cnt, 1)
    return weights, float(weights.sum())


def PlotColumnProfile(
    values: Sequence[float],
    outFile: PathLike,
    ylabel: str = "score",
    title: str = "Column profile",
    xlabel: str = "alignment column (1-based)",
    ylim: Optional[Tuple[float, float]] = None,
    dpi: int = 150,
) -> Path:
    """Generic per-column line plot (used by conservation / entropy / gaps)."""
    import matplotlib.pyplot as plt

    values = np.asarray(values, dtype=float)
    outPath = Path(outFile).expanduser().resolve()
    outPath.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 3))
    x = np.arange(1, len(values) + 1)
    ax.plot(x, values, color="#1f4e79", lw=1.0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(1, max(len(values), 1))
    if ylim is not None:
        ax.set_ylim(*ylim)
    fig.tight_layout()
    fig.savefig(outPath, dpi=dpi)
    plt.close(fig)
    print(f"msaStats: wrote {outPath}")
    return outPath


def PlotConservation(records: Sequence, outFile: PathLike = "conservation.png", **kwargs) -> Path:
    scores = ColumnConservation(records)
    return PlotColumnProfile(
        scores,
        outFile,
        ylabel="conservation (consensus freq)",
        title=kwargs.pop("title", "Per-column conservation"),
        ylim=(0, 1),
        **kwargs,
    )


def PlotEntropy(records: Sequence, outFile: PathLike = "entropy.png", **kwargs) -> Path:
    ent = ColumnEntropy(records)
    return PlotColumnProfile(
        ent,
        outFile,
        ylabel="Shannon entropy (bits)",
        title=kwargs.pop("title", "Per-column Shannon entropy"),
        **kwargs,
    )


def PlotGapProfile(records: Sequence, outFile: PathLike = "gap_profile.png", **kwargs) -> Path:
    gaps = GapFraction(records)
    return PlotColumnProfile(
        gaps,
        outFile,
        ylabel="gap fraction",
        title=kwargs.pop("title", "Gap fraction per column"),
        ylim=(0, 1),
        **kwargs,
    )


def PlotSeqIdMatrix(records: Sequence, outFile: PathLike = "seqid_matrix.png", dpi: int = 150) -> np.ndarray:
    import matplotlib.pyplot as plt

    mat = SeqIdentityMatrix(records)
    outPath = Path(outFile).expanduser().resolve()
    outPath.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(mat, vmin=0, vmax=1, cmap="viridis")
    ax.set_title("Sequence identity matrix")
    ax.set_xlabel("sequence index")
    ax.set_ylabel("sequence index")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(outPath, dpi=dpi)
    plt.close(fig)
    print(f"msaStats: wrote {outPath}")
    return mat


def main(argv=None):
    parser = argparse.ArgumentParser(description="MSA column statistics / profiles")
    parser.add_argument("msaFile")
    parser.add_argument("-fmt", default="fasta", choices=["fasta", "a3m", "clustal"])
    parser.add_argument("-out", default="msa_stats")
    parser.add_argument(
        "-plot",
        default="conservation",
        choices=["conservation", "entropy", "gap", "seqid", "all"],
    )
    args = parser.parse_args(argv)
    records = LoadMsa(args.msaFile, fmt=args.fmt)
    outDir = Path(args.out)
    outDir.mkdir(parents=True, exist_ok=True)
    plots = [args.plot] if args.plot != "all" else ["conservation", "entropy", "gap", "seqid"]
    for name in plots:
        if name == "conservation":
            PlotConservation(records, outDir / "conservation.png")
        elif name == "entropy":
            PlotEntropy(records, outDir / "entropy.png")
        elif name == "gap":
            PlotGapProfile(records, outDir / "gap_profile.png")
        elif name == "seqid":
            PlotSeqIdMatrix(records, outDir / "seqid_matrix.png")


if __name__ == "__main__":
    main()
