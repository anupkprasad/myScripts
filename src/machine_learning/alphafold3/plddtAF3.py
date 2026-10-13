#!/usr/bin/env python3
"""
plddtAF3 — mean CA pLDDT across AF3 predicted models.

When to use
    You have an AF3 job folder (seed-*_sample-*/model.cif) and want per-residue
    confidence for a selection, averaged over the ensemble.

Important
    In MDAnalysis/VMD, multi-CIF trajectories keep B-factors from *frame 0*
    only. This module re-reads the B-factor (pLDDT) from *each* CIF.

Public API
    MeanCAPlddt(universe|directory, sel=...) → resids, chainIds, mean, plddt
    PlotCAPlddt(...)                         → PNG of mean ± std vs CA index

Related modules
    utils.parseFiles.LoadAf3  → load seed CIFs as one Universe
    extractAF3                → MSA/templates from *_data.json (not pLDDT)

Example (Python)
    from machine_learning.alphafold3 import PlotCAPlddt, MeanCAPlddt
    from utils.parseFiles import LoadAf3

    u = LoadAf3("/path/to/af3_job")
    PlotCAPlddt(
        u,
        sel="protein and name CA and chainID A",
        outFile="plddt.png",
    )
    # or without loading coords first:
    # MeanCAPlddt(directory="/path/to/af3_job", sel="protein and name CA")

Example (CLI)
    python -m machine_learning.alphafold3.plddtAF3 \\
        -dir ./af3_job -sel "protein and name CA and chainID A" -out plddt.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np

PathLike = Union[str, Path]

HELP = """
plddtAF3 — mean CA pLDDT across AF3 models

  python -m machine_learning.alphafold3.plddtAF3 -dir ./af3_job \\
      -sel "protein and name CA and chainID A" -out plddt.png

pLDDT is read from each CIF B-factor (not frozen MDA tempfactors).
"""


def _resolveCifPaths(
    universe=None,
    directory: Optional[PathLike] = None,
    pattern: str = "seed",
) -> List[Path]:
    from utils.parseFiles import Af3CifPaths

    if universe is not None and getattr(universe, "_af3CifPaths", None):
        return list(universe._af3CifPaths)
    if directory is None:
        raise ValueError("provide directory= or a Universe from LoadAf3(...)")
    return Af3CifPaths(directory, pattern=pattern)


def MeanCAPlddt(
    universe=None,
    directory: Optional[PathLike] = None,
    sel: str = "protein and name CA",
    pattern: str = "seed",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Mean pLDDT of selected CA atoms across AF3 models.

    Reads B-factors from **each** CIF (trajectory tempfactors stay at frame 0).

    Returns ``resids, chainIds, meanPlddt, plddt`` with
    ``plddt`` shape ``(nFrames, nCA)``.
    """
    from utils.parseFiles import _bfactorsFromPath, _universeFromStructureFile

    paths = _resolveCifPaths(universe=universe, directory=directory, pattern=pattern)
    if universe is None:
        universe = _universeFromStructureFile(paths[0])

    atoms = universe.select_atoms(sel)
    if atoms.n_atoms == 0:
        raise ValueError(f"empty selection: {sel!r}")
    if not np.all(atoms.names == "CA"):
        atoms = atoms.select_atoms("name CA")
        if atoms.n_atoms == 0:
            raise ValueError(f"selection has no CA atoms: {sel!r}")

    indices = atoms.indices
    nFrames = len(paths)
    plddt = np.empty((nFrames, atoms.n_atoms), dtype=np.float32)
    for i, path in enumerate(paths):
        bf = _bfactorsFromPath(path)
        if bf.shape[0] <= int(indices.max()):
            raise ValueError(f"atom count mismatch in {path}")
        plddt[i] = bf[indices]

    meanPlddt = plddt.mean(axis=0)
    resids = atoms.resids.copy()
    try:
        chainIds = np.asarray(atoms.chainIDs)
    except Exception:
        chainIds = np.asarray(["?"] * atoms.n_atoms)
    return resids, chainIds, meanPlddt, plddt


def PlotCAPlddt(
    universe=None,
    directory: Optional[PathLike] = None,
    sel: str = "protein and name CA",
    pattern: str = "seed",
    outFile: Optional[PathLike] = None,
    title: Optional[str] = None,
    show: bool = False,
):
    """Plot mean CA pLDDT (± std across models). Returns ``(resids, meanPlddt, fig)``."""
    import matplotlib.pyplot as plt

    resids, chainIds, meanPlddt, plddt = MeanCAPlddt(
        universe=universe,
        directory=directory,
        sel=sel,
        pattern=pattern,
    )
    stdPlddt = plddt.std(axis=0)
    x = np.arange(len(resids))

    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(x, meanPlddt, color="#1f4e79", lw=1.2, label="mean pLDDT")
    ax.fill_between(
        x,
        meanPlddt - stdPlddt,
        meanPlddt + stdPlddt,
        color="#1f4e79",
        alpha=0.2,
        label="±1σ",
    )
    ax.set_ylim(0, 100)
    ax.set_xlabel("CA index (selection order)")
    ax.set_ylabel("pLDDT")
    ax.set_title(title or f"Mean CA pLDDT ({plddt.shape[0]} models)\n{sel}")
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()

    if outFile is not None:
        outPath = Path(outFile).expanduser().resolve()
        outPath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outPath, dpi=150)
        print(f"plddtAF3: wrote {outPath}")
    if show:
        plt.show()
    elif outFile is not None:
        plt.close(fig)
    return resids, meanPlddt, fig


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("command", nargs="?", default="plddt")
    parser.add_argument("-dir", default=".")
    parser.add_argument("-out", default="af3_ca_plddt.png")
    parser.add_argument("-sel", default="protein and name CA")
    parser.add_argument("-pattern", default="seed")
    parser.add_argument("-h", "--help", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.command in {"help", "-h", "--help"} or args.help:
        print(HELP)
        return 0

    resids, meanPlddt, _ = PlotCAPlddt(
        directory=args.dir,
        sel=args.sel,
        pattern=args.pattern,
        outFile=args.out,
    )
    print(
        f"plddtAF3: {len(resids)} CA atoms, "
        f"mean pLDDT {float(meanPlddt.mean()):.2f} "
        f"(min {float(meanPlddt.min()):.2f}, max {float(meanPlddt.max()):.2f})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
