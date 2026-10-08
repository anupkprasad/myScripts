#!/usr/bin/env python3
"""
Secondary-structure heatmap for a trajectory selection.

Uses MDAnalysis.analysis.dssp (default). If mdtraj is installed, you can
request it with -backend mdtraj.

CLI:

    conda activate ComBio
    python secStrHeatmap.py -top top.pdb -traj traj.xtc \\
        -sel "protein" -out ss_heatmap.png

    python secStrHeatmap.py -top top.pdb -traj traj.xtc \\
        -sel "protein and resid 1-120" -stride 5 \\
        -matrix ss_matrix.dat -out ss_heatmap.png

From Python:

    from secStrHeatmap import compute_secstr, plot_secstr_heatmap
    frames, resids, codes = compute_secstr("top.pdb", "traj.xtc", sel="protein")
    plot_secstr_heatmap(codes, frames, resids, outfile="ss_heatmap.png")
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# DSSP letter -> display index (matches colormap order below)
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
# Map alternate coil / unknown letters onto '-'
_SS_ALIASES = {
    "C": "-",
    " ": "-",
    "~": "-",
    "NA": "-",
    "X": "-",
}


HELP = """
SecStrHeatmap — residue × frame secondary-structure heatmap

Usage:
  python secStrHeatmap.py [options]
  python secStrHeatmap.py help

Options:
  -top <file>       topology (PDB, GRO, TPR, ...) [required]
  -traj <file>      trajectory (XTC, DCD, TRR, ...); omit for single structure
  -sel <text>       MDAnalysis selection of residues to analyze
                    (default: protein). Backbone atoms of those residues are used.
  -first <i>        first frame (default: 0)
  -last <i|end>     last frame (default: end)
  -stride <n>       take every n-th frame (default: 1)
  -backend <name>   mda | mdtraj (default: mda)
  -simplified       collapse DSSP to helix / strand / coil (3 classes)
  -matrix <file>    write residue×frame letter matrix (.dat)
  -out <file>       PNG heatmap (default: ss_heatmap.png)
  -dpi <n>          figure DPI (default: 200)
  -title <text>     plot title

Examples:
  python secStrHeatmap.py -top top.pdb -traj traj.xtc -sel "protein" -out ss.png
  python secStrHeatmap.py -top top.pdb -traj traj.xtc \\
      -sel "protein and resid 40-200" -stride 10 -matrix ss.dat -out ss.png
"""


def _require_mda():
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise SystemExit(
            "MDAnalysis is required. Activate ComBio, e.g.\n"
            "  conda activate ComBio"
        ) from exc
    return mda


def _require_mdtraj():
    try:
        import mdtraj as md
    except ImportError as exc:
        raise SystemExit(
            "mdtraj is not installed. Use -backend mda, or:\n"
            "  conda activate ComBio && pip install mdtraj"
        ) from exc
    return md


def load_universe(top, traj=None, **kwargs):
    mda = _require_mda()
    if traj is None:
        return mda.Universe(top, **kwargs)
    return mda.Universe(top, traj, **kwargs)


def frame_indices(u, first=0, last=None, stride=1):
    n = len(u.trajectory)
    if last is None or last == "end":
        last = n - 1
    first, last, stride = int(first), int(last), int(stride)
    if first < 0 or last >= n or first > last or stride < 1:
        raise ValueError(
            f"bad frame range first={first} last={last} stride={stride} (n={n})"
        )
    return list(range(first, last + 1, stride))


def _residues_for_selection(u, sel):
    """Residues touched by sel; DSSP needs protein atoms of those residues."""
    picked = u.select_atoms(sel)
    if picked.n_atoms == 0:
        raise ValueError(f"selection matched no atoms: {sel!r}")
    residues = picked.residues
    atoms = residues.atoms.select_atoms("protein")
    if atoms.n_residues == 0:
        raise ValueError(
            f"selection {sel!r} has no protein residues; DSSP needs protein atoms"
        )
    return atoms


def _normalize_codes(codes, simplified=False):
    """codes: (n_frames, n_res) array of str → canonical DSSP letters."""
    flat = np.asarray(codes, dtype=object).reshape(-1)
    out = np.empty(flat.shape[0], dtype="<U1")
    for i, raw in enumerate(flat):
        c = str(raw).strip().upper()
        c = _SS_ALIASES.get(c, c)
        if c not in SS_ORDER:
            c = "-"
        out[i] = c
    out = out.reshape(codes.shape)
    if not simplified:
        return out
    simp = np.empty_like(out)
    helix = np.isin(out, ("H", "G", "I"))
    strand = np.isin(out, ("E", "B"))
    simp[:] = "-"
    simp[helix] = "H"
    simp[strand] = "E"
    return simp


def compute_secstr_mda(u, sel="protein", frames=None, simplified=False, guess_hydrogens=True):
    """
    Run MDAnalysis DSSP on residues in sel.

    Returns
    -------
    frames : list[int]
    resids : np.ndarray shape (n_res,)
    codes  : np.ndarray shape (n_frames, n_res) of DSSP letters
    """
    from MDAnalysis.analysis.dssp import DSSP

    atoms = _residues_for_selection(u, sel)
    if frames is None:
        frames = list(range(len(u.trajectory)))
    if not frames:
        raise ValueError("no frames to analyze")

    # Prefer a single run with start/stop/step when the frame list is arithmetic.
    if len(frames) == 1:
        step = 1
        regular = True
    else:
        step = frames[1] - frames[0]
        regular = step > 0 and frames == list(range(frames[0], frames[-1] + 1, step))

    if regular:
        dssp = DSSP(atoms, guess_hydrogens=guess_hydrogens).run(
            start=frames[0], stop=frames[-1] + 1, step=step
        )
        codes = np.asarray(dssp.results.dssp)
        resids = np.asarray(dssp.results.resids)
    else:
        dssp = DSSP(atoms, guess_hydrogens=guess_hydrogens).run(
            start=min(frames), stop=max(frames) + 1, step=1
        )
        all_codes = np.asarray(dssp.results.dssp)
        resids = np.asarray(dssp.results.resids)
        offset = min(frames)
        codes = np.stack([all_codes[f - offset] for f in frames], axis=0)

    codes = _normalize_codes(codes, simplified=simplified)
    return frames, resids, codes


def compute_secstr_mdtraj(top, traj, sel="protein", frames=None, simplified=False):
    """
    DSSP via mdtraj. Residue list is taken from an MDAnalysis selection
    (same resid numbers), then columns are sliced from mdtraj's DSSP array.
    """
    md = _require_mdtraj()
    mda = _require_mda()

    t = md.load(top) if traj is None else md.load(traj, top=top)
    u = mda.Universe(top) if traj is None else mda.Universe(top, traj)
    atoms = _residues_for_selection(u, sel)
    selected_resids = [int(r.resid) for r in atoms.residues]
    want = set(selected_resids)

    res_indices, resids = [], []
    for i, res in enumerate(t.topology.residues):
        rid = int(res.resSeq)
        if rid in want:
            res_indices.append(i)
            resids.append(rid)
    if not res_indices:
        raise ValueError("could not map selection residues onto mdtraj topology")

    if frames is None:
        frames = list(range(t.n_frames))
    t = t[frames]

    dssp = md.compute_dssp(t, simplified=bool(simplified))
    codes = _normalize_codes(np.asarray(dssp[:, res_indices]), simplified=simplified)
    return frames, np.asarray(resids, dtype=int), codes


def compute_secstr(
    top,
    traj=None,
    sel="protein",
    first=0,
    last=None,
    stride=1,
    backend="mda",
    simplified=False,
    guess_hydrogens=True,
    universe=None,
):
    """
    Compute DSSP secondary structure for residues in sel.

    Returns frames, resids, codes with codes.shape == (n_frames, n_res).
    """
    backend = backend.lower()
    if backend == "mdtraj":
        u = universe or load_universe(top, traj)
        frames = frame_indices(u, first=first, last=last, stride=stride)
        return compute_secstr_mdtraj(
            top, traj, sel=sel, frames=frames, simplified=simplified
        )

    if backend != "mda":
        raise ValueError("backend must be 'mda' or 'mdtraj'")

    u = universe or load_universe(top, traj)
    frames = frame_indices(u, first=first, last=last, stride=stride)
    return compute_secstr_mda(
        u,
        sel=sel,
        frames=frames,
        simplified=simplified,
        guess_hydrogens=guess_hydrogens,
    )


def codes_to_int(codes, simplified=False):
    """Map DSSP letters to integers for pcolormesh."""
    order = ("H", "E", "-") if simplified else SS_ORDER
    table = {c: i for i, c in enumerate(order)}
    flat = np.asarray(codes).reshape(-1)
    out = np.array([table.get(str(c), table["-"]) for c in flat], dtype=int)
    return out.reshape(codes.shape), order


def write_secstr_matrix(path, frames, resids, codes):
    """Write frames as rows, residues as columns (letter codes)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("# SecStrHeatmap DSSP matrix (rows=frames, cols=residues)\n")
        handle.write(f"# frames: {' '.join(str(f) for f in frames)}\n")
        handle.write(f"# resids: {' '.join(str(r) for r in resids)}\n")
        handle.write("frame " + " ".join(str(r) for r in resids) + "\n")
        for fi, frame in enumerate(frames):
            row = " ".join(str(codes[fi, j]) for j in range(len(resids)))
            handle.write(f"{frame} {row}\n")
    print(f"SecStrHeatmap: wrote {path}")


def read_secstr_matrix(path):
    """Read a matrix written by write_secstr_matrix."""
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
    codes = np.asarray(rows, dtype="<U1")
    return frames, np.asarray(resids, dtype=int), codes


def plot_secstr_heatmap(
    codes,
    frames,
    resids,
    outfile="ss_heatmap.png",
    simplified=False,
    dpi=200,
    title=None,
):
    """
    Plot residue (y) × frame (x) secondary-structure heatmap.

    codes : (n_frames, n_res) DSSP letters
    """
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
        raise ValueError(f"codes must be 2-D (n_frames, n_res), got {codes.shape}")

    int_mat, order = codes_to_int(codes, simplified=simplified)
    # pcolormesh expects (n_res, n_frames) with residue on y
    z = int_mat.T

    if simplified:
        colors = ("red", "blue", "white")
        labels = (r"helix", r"strand", "coil")
    else:
        colors = SS_COLORS
        labels = SS_LABELS

    cmap = ListedColormap(list(colors), name="dssp")
    bounds = np.arange(-0.5, len(order) + 0.5, 1.0)
    norm = BoundaryNorm(bounds, cmap.N)

    n_frames, n_res = codes.shape
    # Wide figure for long trajectories
    width = max(6.0, min(16.0, n_frames / 40.0))
    height = max(4.0, min(12.0, n_res / 25.0))
    fig, ax = plt.subplots(figsize=(width, height))

    # Coordinates so each cell is one frame × one residue
    x = np.arange(n_frames + 1)
    y = np.arange(n_res + 1)
    pcm = ax.pcolormesh(x, y, z, cmap=cmap, norm=norm, shading="flat")

    cbar = fig.colorbar(pcm, ax=ax, ticks=np.arange(len(order)), fraction=0.03, pad=0.02)
    cbar.ax.set_yticklabels(list(labels), fontsize=8)
    cbar.set_label("Secondary structure")

    ax.set_xlabel("Frame index (analyzed)")
    ax.set_ylabel("Residue")
    ax.set_title(title or f"Secondary structure ({n_res} residues, {n_frames} frames)")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Tick labels: analyzed-frame index → trajectory frame number; y → resid
    x_step = max(1, n_frames // 10)
    x_ticks = list(range(0, n_frames, x_step))
    if x_ticks[-1] != n_frames - 1:
        x_ticks.append(n_frames - 1)
    ax.set_xticks([t + 0.5 for t in x_ticks])
    ax.set_xticklabels([str(frames[t]) for t in x_ticks], rotation=45, ha="right")

    y_step = max(1, n_res // 12)
    y_ticks = list(range(0, n_res, y_step))
    if y_ticks[-1] != n_res - 1:
        y_ticks.append(n_res - 1)
    ax.set_yticks([t + 0.5 for t in y_ticks])
    ax.set_yticklabels([str(resids[t]) for t in y_ticks])

    fig.tight_layout()
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)
    print(f"SecStrHeatmap: wrote {outfile}")
    return outfile


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="secStrHeatmap.py",
        description="Secondary-structure heatmap from a trajectory",
        add_help=False,
    )
    parser.add_argument("command", nargs="?", default="run")
    parser.add_argument("-top", default=None)
    parser.add_argument("-traj", default=None)
    parser.add_argument("-sel", default="protein")
    parser.add_argument("-first", default=0)
    parser.add_argument("-last", default="end")
    parser.add_argument("-stride", type=int, default=1)
    parser.add_argument("-backend", default="mda", choices=("mda", "mdtraj"))
    parser.add_argument("-simplified", action="store_true")
    parser.add_argument("-matrix", default=None)
    parser.add_argument("-out", default="ss_heatmap.png")
    parser.add_argument("-dpi", type=int, default=200)
    parser.add_argument("-title", default=None)
    parser.add_argument("-h", "--help", action="store_true")
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    if args.command in {"help", "-h", "--help"} or args.help:
        print(HELP)
        return 0
    if args.top is None:
        print("SecStrHeatmap: -top is required", file=sys.stderr)
        print(HELP)
        return 1

    last = None if args.last in (None, "end") else int(args.last)
    frames, resids, codes = compute_secstr(
        args.top,
        traj=args.traj,
        sel=args.sel,
        first=args.first,
        last=last,
        stride=args.stride,
        backend=args.backend,
        simplified=args.simplified,
    )
    print(
        f"SecStrHeatmap: {len(frames)} frames × {len(resids)} residues "
        f"(backend={args.backend})"
    )
    if args.matrix:
        write_secstr_matrix(args.matrix, frames, resids, codes)
    plot_secstr_heatmap(
        codes,
        frames,
        resids,
        outfile=args.out,
        simplified=args.simplified,
        dpi=args.dpi,
        title=args.title,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
