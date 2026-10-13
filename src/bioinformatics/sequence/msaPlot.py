#!/usr/bin/env python3
"""
msaPlot — publication-style MSA figure (conservation strip + colored residues).

When to use
    After LoadMsa / CleanMsa, make a PNG of the alignment with a conservation
    track. Prefer a subset for large MSAs (maxSeqs).

Public API
    PlotMsa(records, outFile=..., maxSeqs=40, queryId=None, showConservation=True)

Example
    from bioinformatics.sequence import LoadMsa, PlotMsa

    records = LoadMsa("clean.fasta")
    PlotMsa(records, outFile="msa.png", maxSeqs=40, queryId="query")

CLI
    python -m bioinformatics.sequence.msaPlot clean.fasta -out msa.png -maxSeqs 40
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np

from .msaIO import LoadMsa
from .msaStats import ColumnConservation
from .msaView import SortBySimilarity, _resolveQuery

PathLike = Union[str, Path]

# ClustalX-like residue colors (hex)
_AA_COLORS = {
    "A": "#80a0f0",
    "I": "#80a0f0",
    "L": "#80a0f0",
    "M": "#80a0f0",
    "F": "#80a0f0",
    "W": "#80a0f0",
    "V": "#80a0f0",  # hydrophobic
    "K": "#f01505",
    "R": "#f01505",  # positive
    "D": "#c048c0",
    "E": "#c048c0",  # negative
    "S": "#15c015",
    "T": "#15c015",  # polar
    "N": "#15c015",
    "Q": "#15c015",
    "C": "#f08080",  # cysteine
    "G": "#f09048",  # glycine
    "P": "#c0c000",  # proline
    "H": "#15a4a4",
    "Y": "#15a4a4",  # aromatic polar
    "B": "#ffffff",
    "Z": "#ffffff",
    "X": "#ffffff",
    "-": "#f0f0f0",
    ".": "#f0f0f0",
}


def _hexToRgb(h: str) -> tuple:
    h = h.lstrip("#")
    return tuple(int(h[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


def _recordsToRgb(records: Sequence) -> np.ndarray:
    """(nSeq, nCol, 3) RGB array from residue letters."""
    arr = np.asarray([list(str(r.seq).upper()) for r in records])
    nSeq, nCol = arr.shape
    rgb = np.ones((nSeq, nCol, 3), dtype=float)
    default = _hexToRgb("#ffffff")
    for aa, hexCol in _AA_COLORS.items():
        mask = arr == aa
        if not mask.any():
            continue
        rgb[mask] = _hexToRgb(hexCol)
    # unknown letters → white
    known = set(_AA_COLORS)
    for i in range(nSeq):
        for j in range(nCol):
            if arr[i, j] not in known:
                rgb[i, j] = default
    return rgb


def _shortLabel(rec, maxLen: int = 28) -> str:
    label = rec.id
    if len(label) > maxLen:
        return label[: maxLen - 1] + "…"
    return label


def PlotMsa(
    records: Sequence,
    outFile: PathLike = "msa.png",
    maxSeqs: int = 40,
    queryId: Optional[str] = None,
    showConservation: bool = True,
    sortByQuery: bool = True,
    dpi: int = 150,
    title: Optional[str] = None,
):
    """
    Plot MSA as a colored residue matrix with optional conservation strip.

    - maxSeqs: keep query + nearest neighbors (if sortByQuery) or first N rows
    - showConservation: top panel = ColumnConservation
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    if not records:
        raise ValueError("empty alignment")

    work = list(records)
    if sortByQuery:
        work = SortBySimilarity(work, queryId=queryId)
    else:
        # still resolve query to put first if requested
        if queryId is not None:
            q = _resolveQuery(work, queryId)
            work = [q] + [r for r in work if r.id != q.id]

    if maxSeqs is not None and len(work) > maxSeqs:
        work = work[:maxSeqs]

    nSeq = len(work)
    nCol = len(str(work[0].seq))
    rgb = _recordsToRgb(work)
    cons = ColumnConservation(work) if showConservation else None

    # figure size scales with MSA dimensions
    width = max(8.0, min(22.0, nCol / 18.0))
    height = max(4.0, min(14.0, nSeq / 4.0 + (1.2 if showConservation else 0.0)))

    if showConservation:
        fig, axes = plt.subplots(
            2,
            1,
            figsize=(width, height),
            gridspec_kw={"height_ratios": [1, max(3, nSeq / 8)], "hspace": 0.08},
            sharex=True,
            constrained_layout=True,
        )
        axCons, axMsa = axes
        x = np.arange(1, nCol + 1)
        axCons.fill_between(x, cons, color="#1f4e79", alpha=0.85, step="mid")
        axCons.set_ylim(0, 1)
        axCons.set_ylabel("cons.")
        axCons.set_xlim(0.5, nCol + 0.5)
        axCons.set_yticks([0, 0.5, 1.0])
        axCons.tick_params(labelbottom=False)
        if title:
            axCons.set_title(title)
        else:
            axCons.set_title(f"MSA ({nSeq} sequences × {nCol} columns)")
    else:
        fig, axMsa = plt.subplots(figsize=(width, height), constrained_layout=True)
        if title:
            axMsa.set_title(title)

    axMsa.imshow(rgb, aspect="auto", interpolation="nearest", origin="upper")
    axMsa.set_xlabel("alignment column (1-based)")
    axMsa.set_ylabel("sequence")
    axMsa.set_yticks(np.arange(nSeq))
    axMsa.set_yticklabels([_shortLabel(r) for r in work], fontsize=7)
    # sparse x ticks
    xStep = max(1, nCol // 12)
    xTicks = list(range(0, nCol, xStep))
    if xTicks[-1] != nCol - 1:
        xTicks.append(nCol - 1)
    axMsa.set_xticks(xTicks)
    axMsa.set_xticklabels([str(t + 1) for t in xTicks], fontsize=8)

    # small legend of residue classes
    legendHandles = [
        Patch(facecolor=_AA_COLORS["A"], label="hydrophobic"),
        Patch(facecolor=_AA_COLORS["K"], label="positive"),
        Patch(facecolor=_AA_COLORS["D"], label="negative"),
        Patch(facecolor=_AA_COLORS["S"], label="polar"),
        Patch(facecolor=_AA_COLORS["G"], label="G"),
        Patch(facecolor=_AA_COLORS["P"], label="P"),
        Patch(facecolor=_AA_COLORS["C"], label="C"),
        Patch(facecolor=_AA_COLORS["-"], label="gap"),
    ]
    axMsa.legend(
        handles=legendHandles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        fontsize=7,
        frameon=False,
    )

    outPath = Path(outFile).expanduser().resolve()
    outPath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outPath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"msaPlot: wrote {outPath} ({nSeq} × {nCol})")
    return outPath


def main(argv=None):
    parser = argparse.ArgumentParser(description="Plot colored MSA + conservation")
    parser.add_argument("msaFile")
    parser.add_argument("-fmt", default="fasta", choices=["fasta", "a3m", "clustal"])
    parser.add_argument("-out", default="msa.png")
    parser.add_argument("-maxSeqs", type=int, default=40)
    parser.add_argument("-queryId", default=None)
    parser.add_argument("-noCons", action="store_true")
    parser.add_argument("-noSort", action="store_true")
    args = parser.parse_args(argv)
    records = LoadMsa(args.msaFile, fmt=args.fmt)
    PlotMsa(
        records,
        outFile=args.out,
        maxSeqs=args.maxSeqs,
        queryId=args.queryId,
        showConservation=not args.noCons,
        sortByQuery=not args.noSort,
    )


if __name__ == "__main__":
    main()
