#!/usr/bin/env python3
"""
phyloTree — distance-based phylogenetic tree from an MSA (NJ / UPGMA).

When to use
    Quick tree from cleaned MSA identity distances. Not a substitute for
    IQ-TREE / RAxML on large datasets, but good for QC and clade overview.

Public API
    DistanceFromMsa(records) → (names, distMatrix)
    BuildDistanceTree(records, method="nj"|"upgma") → Bio.Phylo tree
    PlotPhyloTree(tree|records, outFile=..., method="nj")
    WriteNewick(tree, outFile)

Example
    from bioinformatics.sequence import LoadMsa, PlotPhyloTree, BuildDistanceTree

    records = LoadMsa("clean.fasta")
    PlotPhyloTree(records[:80], outFile="tree.png", method="nj")
    tree = BuildDistanceTree(records[:80], method="nj")

CLI
    python -m bioinformatics.sequence.phyloTree clean.fasta -out tree.png -maxSeqs 80
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np

from .msaIO import LoadMsa
from .msaStats import SeqIdentityMatrix
from .msaView import SortBySimilarity

PathLike = Union[str, Path]


def _uniqueNames(records: Sequence, maxLen: int = 40) -> List[str]:
    """Phylo needs unique tip labels; shorten and disambiguate collisions."""
    names: List[str] = []
    seen = {}
    for i, rec in enumerate(records):
        raw = re.sub(r"[^\w.|/-]+", "_", rec.id).strip("_") or f"seq{i}"
        if len(raw) > maxLen:
            raw = raw[: maxLen - 1] + "~"
        if raw in seen:
            seen[raw] += 1
            raw = f"{raw}_{seen[raw]}"
        else:
            seen[raw] = 0
        names.append(raw)
    return names


def DistanceFromMsa(records: Sequence) -> Tuple[List[str], np.ndarray]:
    """
    Names + full distance matrix with d_ij = 1 − sequence identity.

    Diagonal is 0.
    """
    if len(records) < 2:
        raise ValueError("need at least 2 sequences for a tree")
    names = _uniqueNames(records)
    idMat = SeqIdentityMatrix(records)
    dist = np.clip(1.0 - idMat, 0.0, 1.0)
    np.fill_diagonal(dist, 0.0)
    return names, dist


def BuildDistanceTree(
    records: Sequence,
    method: str = "nj",
):
    """
    Build a Bio.Phylo tree from MSA identity distances.

    method: ``nj`` (neighbor-joining) or ``upgma``.
    """
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    method = method.lower()
    if method not in {"nj", "upgma"}:
        raise ValueError("method must be 'nj' or 'upgma'")

    names, dist = DistanceFromMsa(records)
    # Bio.Phylo wants lower-triangular rows including diagonal
    lower = []
    for i in range(len(names)):
        lower.append([float(dist[i, j]) for j in range(i + 1)])
    dm = DistanceMatrix(names, matrix=lower)
    constructor = DistanceTreeConstructor()
    if method == "nj":
        tree = constructor.nj(dm)
    else:
        tree = constructor.upgma(dm)
    return tree


def WriteNewick(tree, outFile: PathLike) -> Path:
    """Write Newick (.nwk / .tree)."""
    from Bio import Phylo

    outPath = Path(outFile).expanduser().resolve()
    outPath.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(outPath), "newick")
    print(f"phyloTree: wrote {outPath}")
    return outPath


def PlotPhyloTree(
    treeOrRecords,
    outFile: PathLike = "phylo_tree.png",
    method: str = "nj",
    maxSeqs: Optional[int] = 80,
    queryId: Optional[str] = None,
    sortByQuery: bool = False,
    dpi: int = 150,
    title: Optional[str] = None,
    writeNewick: bool = True,
):
    """
    Plot a rectangular phylogenetic tree.

    Pass either a Bio.Phylo tree or an MSA record list (tree is built then).
    For large MSAs, set maxSeqs (default 80).
    """
    from Bio import Phylo
    import matplotlib.pyplot as plt

    tree = treeOrRecords
    newickPath = None
    if not hasattr(tree, "clade") and not hasattr(tree, "root"):
        # treat as records
        records = list(treeOrRecords)
        if sortByQuery:
            records = SortBySimilarity(records, queryId=queryId)
        if maxSeqs is not None and len(records) > maxSeqs:
            print(f"phyloTree: using first {maxSeqs} of {len(records)} sequences")
            records = records[:maxSeqs]
        tree = BuildDistanceTree(records, method=method)
        if writeNewick:
            outPath = Path(outFile).expanduser().resolve()
            newickPath = outPath.with_suffix(".nwk")
            WriteNewick(tree, newickPath)

    nTips = tree.count_terminals()
    height = max(4.0, min(24.0, nTips * 0.22 + 1.5))
    width = max(8.0, min(16.0, 8.0 + nTips * 0.02))
    fig, ax = plt.subplots(figsize=(width, height))
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)
    ax.set_ylabel("")
    ax.set_xlabel("distance (1 − sequence identity)")
    ax.set_title(title or f"Distance tree ({method.upper()}, {nTips} tips)")
    # Phylo.draw leaves a lot of whitespace sometimes
    fig.tight_layout()
    outPath = Path(outFile).expanduser().resolve()
    outPath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outPath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"phyloTree: wrote {outPath}")
    return outPath


def main(argv=None):
    parser = argparse.ArgumentParser(description="MSA distance tree (NJ/UPGMA)")
    parser.add_argument("msaFile")
    parser.add_argument("-fmt", default="fasta", choices=["fasta", "a3m", "clustal"])
    parser.add_argument("-out", default="phylo_tree.png")
    parser.add_argument("-method", default="nj", choices=["nj", "upgma"])
    parser.add_argument("-maxSeqs", type=int, default=80)
    parser.add_argument("-queryId", default=None)
    parser.add_argument("-sortByQuery", action="store_true")
    args = parser.parse_args(argv)
    records = LoadMsa(args.msaFile, fmt=args.fmt)
    PlotPhyloTree(
        records,
        outFile=args.out,
        method=args.method,
        maxSeqs=args.maxSeqs,
        queryId=args.queryId,
        sortByQuery=args.sortByQuery,
    )


if __name__ == "__main__":
    main()
