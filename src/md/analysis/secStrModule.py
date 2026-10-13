#!/usr/bin/env python3
"""
Secondary structure: compute (MDAnalysis DSSP) and/or plot (VMD letter files).

VMD writes per-frame letter strings (e.g. ss_traj.dat). This module can plot
those, or compute DSSP from scratch on a trajectory / existing Universe.

    # Plot data already calculated in VMD
    secstr-module from-vmd -file ss_traj.dat -resNum 17 -out ss.png

    # Compute with MDAnalysis DSSP then plot heatmap
    secstr-module from-traj -top top.pdb -traj traj.xtc \\
        -sel "protein" -stride 5 -out ss.png -matrix ss.dat

    # From Python with an existing Universe (package path after pip install -e .)
    from md.analysis.secStrModule import computeSecStr, plotSecStrHeatmap, SecStr
    frames, resids, codes = computeSecStr(universe=u, sel="protein")
    plotSecStrHeatmap(codes, frames, resids, outFile="ss.png")
    SecStr("ss_traj.dat", resNum=17).plotBasicSs("ss_panels.png")
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Union

import numpy as np

# DSSP letter order for heatmaps (MDAnalysis / VMD-compatible)
SS_ORDER = ("H", "G", "I", "E", "B", "T", "S", "-")
SS_LABELS = (
    r"$\alpha$-helix",
    r"$3_{10}$-helix",
    r"$\pi$-helix",
    r"$\beta$-sheet",
    r"$\beta$-bridge",
    "turn",
    "bend",
    "coil",
)
SS_COLORS = (
    "red",
    "limegreen",
    "magenta",
    "blue",
    "tan",
    "cyan",
    "orange",
    "white",
)
_SS_ALIASES = {
    "C": "-",
    " ": "-",
    "~": "-",
    "NA": "-",
    "X": "-",
}

PathLike = Union[str, Path]

HELP = """
SecStrModule — secondary structure (VMD file plot or MDAnalysis DSSP)

Usage:
  secstr-module <command> [options]
  python -m md.analysis.secStrModule help

Commands:
  from-vmd     Plot panels from a VMD letter file (e.g. ss_traj.dat)
  from-traj    Compute DSSP with MDAnalysis, write optional matrix, plot heatmap
  help         Show this message

from-vmd options:
  -file <path>     VMD secondary-structure letter file (default: ss_traj.dat)
  -pepIdx all|list peptide indices (default: all)
  -resNum <n>      residues per peptide copy (default: 17)
  -first <i>       first frame index in file (default: 0)
  -last <i>        last frame index exclusive; -1 = end (default: -1)
  -out <png>       output figure (default: ss_vmd.png)

from-traj options:
  -top <file>      topology (required unless you only use the Python API)
  -traj <file>     trajectory (optional for single structure)
  -sel <text>      MDAnalysis selection (default: protein)
  -first <i>       first trajectory frame (default: 0)
  -last <i|end>    last frame (default: end)
  -stride <n>      stride (default: 1)
  -simplified      helix / strand / coil only
  -matrix <file>   write residue×frame letter matrix
  -out <png>       heatmap PNG (default: ss_heatmap.png)
  -dpi <n>         figure DPI (default: 200)
  -title <text>    plot title

Python API:
  from md.analysis.secStrModule import SecStr, computeSecStr, plotSecStrHeatmap
  SecStr("ss_traj.dat", resNum=17).plotBasicSs("out.png")
  frames, resids, codes = computeSecStr(universe=u, sel="protein")
  plotSecStrHeatmap(codes, frames, resids, outFile="ss.png")
"""


# ---------------------------------------------------------------------------
# VMD letter-file analysis (existing workflow)
# ---------------------------------------------------------------------------

def runningMean(y, windowLen):
    """Smooth 1-D series; uses tsmoothie if available, else numpy convolution."""
    try:
        from tsmoothie.smoother import ConvolutionSmoother

        smoother = ConvolutionSmoother(window_len=windowLen, window_type="ones")
        smoother.smooth(y)
        return smoother.smooth_data[0]
    except ImportError:
        y = np.asarray(y, dtype=float)
        if windowLen < 1:
            return y
        kernel = np.ones(windowLen) / float(windowLen)
        return np.convolve(y, kernel, mode="same")


class SecStr:
    """
    Analyse / plot secondary structure from a VMD letter file.

    Each line is one frame: a string of DSSP-like letters (H, G, I, E, B, T, C, …),
    optionally concatenating several identical peptide copies.
    """

    def __init__(self, fileName="ss_traj.dat", pepIdx="all", resNum=17):
        self.fileName = fileName
        self.pepIdx = pepIdx
        self.resNum = resNum
        self.ssOrder = {
            r"$\alpha$-helix": "H",
            r"$3_{10}$-helix": "G",
            r"$\pi$-helix": "I",
            r"$\beta$-sheet": "E",
            r"$\beta$-bridge": "B",
            "turn": "T",
            "coil": "C",
        }
        with open(fileName) as handle:
            self.dataAll = handle.read().replace(" ", "").split("\n")[:-1]

    def dataPepIdx(self):
        """Return array of shape (nFrames, nResSelected) of letter codes."""
        pepIdx = self.pepIdx
        dataAll = self.dataAll
        resNum = self.resNum

        if pepIdx == "all":
            data = dataAll
        else:
            arr = []
            for i in range(len(dataAll)):
                data1 = ""
                for j in pepIdx:
                    c1 = resNum * j
                    if c1 >= len(dataAll[0]):
                        raise SystemExit("pepIdx index out of range")
                    c2 = resNum + resNum * j
                    data1 = data1 + dataAll[i][c1:c2]
                arr.append(data1)
            data = arr
        return np.array([list(row) for row in data])

    def ssPercent(self):
        """Per-frame secondary-structure percentages."""
        array = self.dataPepIdx()
        arraySp = np.zeros((len(array), 8))
        sso = list(self.ssOrder.values())
        allRes = len(array[0])
        for j in range(len(array)):
            arraySp[j, 0] = j + 1
            for i, st in enumerate(sso):
                arraySp[j, i + 1] = (np.count_nonzero(array[j] == st) / allRes) * 100
        return arraySp

    def ssFreq(self, timeI=0, timeF=-1):
        """Per-residue frequency of each SS type over frames."""
        array = self.dataPepIdx()
        resNum = self.resNum
        if not (timeI == 0 and timeF == -1):
            array = array[timeI:timeF]
        arraySsf = np.zeros((resNum, 8))
        allRes = len(array[0])
        divisorFactor = len(array) * (allRes / resNum)
        array = array.T
        sso = list(self.ssOrder.values())
        for j in range(resNum):
            arraySsf[j, 0] = j + 1
            for i, st in enumerate(sso):
                arraySsf[j, i + 1] = (
                    np.count_nonzero(array[j:allRes:resNum] == st) / divisorFactor
                ) * 100
        return arraySsf

    def ssTimeAvg(self, timeI=0, timeF=-1):
        """Average of ssPercent over a frame window."""
        arraySp = self.ssPercent()
        if not (timeI == 0 and timeF == -1):
            arraySp = arraySp[timeI:timeF]
        div = len(arraySp)
        arraySp = np.sum(arraySp, axis=0) / div
        arraySp[0] = timeF - timeI
        return arraySp

    def ssHeatmap(self, timeI=0, timeF=-1):
        """Integer matrix (nFrames, nRes) for pcolormesh (x=residue, y=frame)."""
        array = self.dataPepIdx()
        if not (timeI == 0 and timeF == -1):
            array = array[timeI:timeF]
        sso = list(self.ssOrder.values())
        # work on a copy so we do not corrupt letter data
        work = array.copy()
        for i, st in enumerate(sso):
            work[work == st] = i + 1
        return work.astype(int)

    def plotBasicSs(self, outFile, timeI=0, timeF=-1):
        """Three-panel figure: heatmap, % vs time, per-residue frequency."""
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator

        fig, axs = plt.subplots(3, 1, figsize=(3.5, 6), constrained_layout=True)
        ssf = self.ssFreq(timeI=timeI, timeF=timeF)
        ssp = self.ssPercent()
        ssHmp = self.ssHeatmap(timeI=timeI, timeF=timeF)

        axs[1].plot(ssp[:, 1], color="red", alpha=0.3)
        axs[1].plot(runningMean(ssp[:, 1], 10), label=r"$\alpha$-helix", color="red")
        axs[1].plot(runningMean(ssp[:, 2], 10), label=r"$3_{10}$-helix", color="green")
        axs[1].plot(ssp[:, 2], color="green", alpha=0.3)
        axs[1].plot(ssp[:, 4], label=r"$\beta$-sheet", color="blue")
        axs[2].plot(ssf[:, 1], label=r"$\alpha$-helix", color="red", marker="o")
        axs[2].plot(ssf[:, 2], label=r"$3_{10}$-helix", color="green", marker="o")
        axs[2].plot(ssf[:, 4], label=r"$\beta$-sheet", color="blue", marker="o")

        cmap = matplotlib.colors.ListedColormap(
            ["red", "limegreen", "magenta", "blue", "tan", "cyan", "white"], 7
        )
        mesh = axs[0].pcolormesh(ssHmp, cmap=cmap)
        evenly = (ssHmp.max() - ssHmp.min()) / 7.0
        ticks = np.arange(
            ssHmp.min() + evenly / 2, ssHmp.max() + evenly / 2, evenly
        )
        cbar = fig.colorbar(
            mesh, ax=axs[0], ticks=ticks, location="top", fraction=1, aspect=50
        )
        cbar.ax.set_xticklabels(list(self.ssOrder.keys()), fontsize=6)
        yLabel = ["Frame", "Secondary structure (%)", "Frequency (%)"]
        xLabel = ["Residue index", "Time (ns)", "Residue index"]
        for i in range(3):
            axs[i].set_ylabel(yLabel[i])
            axs[i].set_xlabel(xLabel[i])
            if i != 0:
                axs[i].legend(fontsize=6)
        axs[0].yaxis.set_major_locator(MaxNLocator(integer=True))
        axs[0].invert_yaxis()  # frame 0 at top
        fig.savefig(outFile, dpi=600)
        plt.close(fig)
        print(f"SecStrModule: wrote {outFile}")


# Backward-compatible names
Sec_Str = SecStr


def running_mean(y, N):
    return runningMean(y, N)


# ---------------------------------------------------------------------------
# MDAnalysis DSSP compute + heatmap
# ---------------------------------------------------------------------------

def _requireMda():
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise SystemExit(
            "MDAnalysis is required. Activate ComBio, e.g.\n"
            "  conda activate ComBio"
        ) from exc
    return mda


def loadUniverse(top, traj=None, **kwargs):
    mda = _requireMda()
    if traj is None:
        return mda.Universe(top, **kwargs)
    return mda.Universe(top, traj, **kwargs)


def frameIndices(u, first=0, last=None, stride=1) -> List[int]:
    n = len(u.trajectory)
    if last is None or last == "end":
        last = n - 1
    first, last, stride = int(first), int(last), int(stride)
    if first < 0 or last >= n or first > last or stride < 1:
        raise ValueError(
            f"bad frame range first={first} last={last} stride={stride} (n={n})"
        )
    return list(range(first, last + 1, stride))


def _residuesForSelection(u, sel):
    picked = u.select_atoms(sel)
    if picked.n_atoms == 0:
        raise ValueError(f"selection matched no atoms: {sel!r}")
    atoms = picked.residues.atoms.select_atoms("protein")
    if atoms.n_residues == 0:
        raise ValueError(
            f"selection {sel!r} has no protein residues; DSSP needs protein atoms"
        )
    return atoms


def _normalizeCodes(codes, simplified=False):
    flat = np.asarray(codes, dtype=object).reshape(-1)
    out = np.empty(flat.shape[0], dtype="<U1")
    for i, raw in enumerate(flat):
        c = str(raw).strip().upper()
        c = _SS_ALIASES.get(c, c)
        if c not in SS_ORDER:
            c = "-"
        out[i] = c
    out = out.reshape(np.asarray(codes).shape)
    if not simplified:
        return out
    simp = np.empty_like(out)
    helix = np.isin(out, ("H", "G", "I"))
    strand = np.isin(out, ("E", "B"))
    simp[:] = "-"
    simp[helix] = "H"
    simp[strand] = "E"
    return simp


def computeSecStr(
    top=None,
    traj=None,
    sel="protein",
    first=0,
    last=None,
    stride=1,
    simplified=False,
    guessHydrogens=True,
    universe=None,
    frames: Optional[Sequence[int]] = None,
):
    """
    Compute DSSP secondary structure with MDAnalysis.

    Pass an existing ``universe=u`` **or** ``top`` (+ optional ``traj``).

    Returns
    -------
    frames : list[int]
    resids : np.ndarray (nRes,)
    codes  : np.ndarray (nFrames, nRes) of DSSP letters
    """
    from MDAnalysis.analysis.dssp import DSSP

    if universe is None:
        if top is None:
            raise ValueError("computeSecStr: pass universe=... or top=...")
        u = loadUniverse(top, traj)
    else:
        u = universe

    if frames is None:
        frames = frameIndices(u, first=first, last=last, stride=stride)
    frames = list(frames)
    if not frames:
        raise ValueError("no frames to analyze")

    atoms = _residuesForSelection(u, sel)

    if len(frames) == 1:
        step = 1
        regular = True
    else:
        step = frames[1] - frames[0]
        regular = step > 0 and frames == list(range(frames[0], frames[-1] + 1, step))

    if regular:
        dssp = DSSP(atoms, guess_hydrogens=guessHydrogens).run(
            start=frames[0], stop=frames[-1] + 1, step=step
        )
        codes = np.asarray(dssp.results.dssp)
        resids = np.asarray(dssp.results.resids)
    else:
        dssp = DSSP(atoms, guess_hydrogens=guessHydrogens).run(
            start=min(frames), stop=max(frames) + 1, step=1
        )
        allCodes = np.asarray(dssp.results.dssp)
        resids = np.asarray(dssp.results.resids)
        offset = min(frames)
        codes = np.stack([allCodes[f - offset] for f in frames], axis=0)

    codes = _normalizeCodes(codes, simplified=simplified)
    return frames, resids, codes


def codesToInt(codes, simplified=False):
    order = ("H", "E", "-") if simplified else SS_ORDER
    table = {c: i for i, c in enumerate(order)}
    flat = np.asarray(codes).reshape(-1)
    out = np.array([table.get(str(c), table["-"]) for c in flat], dtype=int)
    return out.reshape(codes.shape), order


def writeSecStrMatrix(path, frames, resids, codes):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("# SecStrModule DSSP matrix (rows=frames, cols=residues)\n")
        handle.write(f"# frames: {' '.join(str(f) for f in frames)}\n")
        handle.write(f"# resids: {' '.join(str(r) for r in resids)}\n")
        handle.write("frame " + " ".join(str(r) for r in resids) + "\n")
        for fi, frame in enumerate(frames):
            row = " ".join(str(codes[fi, j]) for j in range(len(resids)))
            handle.write(f"{frame} {row}\n")
    print(f"SecStrModule: wrote {path}")


def readSecStrMatrix(path):
    path = Path(path)
    frames, rows = [], []
    resids = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0] == "frame":
                resids = [int(float(x)) for x in parts[1:]]
                continue
            frames.append(int(float(parts[0])))
            rows.append(parts[1:])
    if resids is None or not rows:
        raise ValueError(f"no secondary-structure matrix in {path}")
    return frames, np.asarray(resids, dtype=int), np.asarray(rows, dtype="<U1")


def plotSecStrHeatmap(
    codes,
    frames,
    resids,
    outFile="ss_heatmap.png",
    simplified=False,
    dpi=200,
    title=None,
):
    """Plot frame (y) × residue (x) secondary-structure heatmap."""
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import BoundaryNorm, ListedColormap
        from matplotlib.ticker import MaxNLocator
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required. Activate ComBio, e.g.\n"
            "  conda activate ComBio"
        ) from exc

    codes = np.asarray(codes)
    if codes.ndim != 2:
        raise ValueError(f"codes must be 2-D (nFrames, nRes), got {codes.shape}")

    intMat, order = codesToInt(codes, simplified=simplified)
    # pcolormesh C shape is (ny, nx) → y=frame, x=residue
    z = intMat

    if simplified:
        colors = ("red", "blue", "white")
        labels = (r"helix", r"strand", "coil")
    else:
        colors = SS_COLORS
        labels = SS_LABELS

    cmap = ListedColormap(list(colors), name="dssp")
    bounds = np.arange(-0.5, len(order) + 0.5, 1.0)
    norm = BoundaryNorm(bounds, cmap.N)

    nFrames, nRes = codes.shape
    width = max(6.0, min(16.0, nRes / 25.0))
    height = max(4.0, min(12.0, nFrames / 40.0))
    fig, ax = plt.subplots(figsize=(width, height))

    x = np.arange(nRes + 1)
    y = np.arange(nFrames + 1)
    pcm = ax.pcolormesh(x, y, z, cmap=cmap, norm=norm, shading="flat")

    cbar = fig.colorbar(pcm, ax=ax, ticks=np.arange(len(order)), fraction=0.03, pad=0.02)
    cbar.ax.set_yticklabels(list(labels), fontsize=8)
    cbar.set_label("Secondary structure")

    ax.set_xlabel("Residue")
    ax.set_ylabel("Frame")
    ax.set_title(title or f"Secondary structure ({nRes} residues, {nFrames} frames)")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.invert_yaxis()  # frame 0 at top

    xStep = max(1, nRes // 12)
    xTicks = list(range(0, nRes, xStep))
    if xTicks[-1] != nRes - 1:
        xTicks.append(nRes - 1)
    ax.set_xticks([t + 0.5 for t in xTicks])
    ax.set_xticklabels([str(resids[t]) for t in xTicks], rotation=45, ha="right")

    yStep = max(1, nFrames // 10)
    yTicks = list(range(0, nFrames, yStep))
    if yTicks[-1] != nFrames - 1:
        yTicks.append(nFrames - 1)
    ax.set_yticks([t + 0.5 for t in yTicks])
    ax.set_yticklabels([str(frames[t]) for t in yTicks])

    fig.tight_layout()
    outFile = Path(outFile)
    outFile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outFile, dpi=dpi)
    plt.close(fig)
    print(f"SecStrModule: wrote {outFile}")
    return outFile


# Aliases used by older secStrHeatmap imports
compute_secstr = computeSecStr
plot_secstr_heatmap = plotSecStrHeatmap
write_secstr_matrix = writeSecStrMatrix
read_secstr_matrix = readSecStrMatrix
load_universe = loadUniverse


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parsePepIdx(text):
    if text is None or text == "all":
        return "all"
    return [int(x) for x in text.replace(",", " ").split()]


def main(argv=None):
    parser = argparse.ArgumentParser(prog="secStrModule.py", add_help=False)
    parser.add_argument("command", nargs="?", default="help")
    parser.add_argument("-file", default="ss_traj.dat")
    parser.add_argument("-pepIdx", default="all")
    parser.add_argument("-resNum", type=int, default=17)
    parser.add_argument("-top", default=None)
    parser.add_argument("-traj", default=None)
    parser.add_argument("-sel", default="protein")
    parser.add_argument("-first", default=0)
    parser.add_argument("-last", default="end")
    parser.add_argument("-stride", type=int, default=1)
    parser.add_argument("-simplified", action="store_true")
    parser.add_argument("-matrix", default=None)
    parser.add_argument("-out", default=None)
    parser.add_argument("-dpi", type=int, default=200)
    parser.add_argument("-title", default=None)
    parser.add_argument("-h", "--help", action="store_true")
    args = parser.parse_args(argv)

    if args.command in {"help", "-h", "--help"} or args.help:
        print(HELP)
        return 0

    if args.command in {"from-vmd", "vmd", "plot-vmd"}:
        timeF = -1 if str(args.last) in {"end", "-1"} else int(args.last)
        timeI = int(args.first)
        out = args.out or "ss_vmd.png"
        sec = SecStr(
            fileName=args.file,
            pepIdx=_parsePepIdx(args.pepIdx),
            resNum=args.resNum,
        )
        sec.plotBasicSs(out, timeI=timeI, timeF=timeF)
        return 0

    if args.command in {"from-traj", "traj", "compute", "heatmap"}:
        if args.top is None:
            print("SecStrModule from-traj: -top is required", file=sys.stderr)
            return 1
        last = None if args.last in (None, "end") else int(args.last)
        frames, resids, codes = computeSecStr(
            top=args.top,
            traj=args.traj,
            sel=args.sel,
            first=args.first,
            last=last,
            stride=args.stride,
            simplified=args.simplified,
        )
        print(f"SecStrModule: {len(frames)} frames × {len(resids)} residues (MDAnalysis DSSP)")
        if args.matrix:
            writeSecStrMatrix(args.matrix, frames, resids, codes)
        plotSecStrHeatmap(
            codes,
            frames,
            resids,
            outFile=args.out or "ss_heatmap.png",
            simplified=args.simplified,
            dpi=args.dpi,
            title=args.title,
        )
        return 0

    print(f"SecStrModule: unknown command '{args.command}'", file=sys.stderr)
    print(HELP)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
