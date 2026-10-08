#!/usr/bin/env python3
"""
Cluster a trajectory by pairwise RMSD of one atom selection.

Standalone (needs MDAnalysis + scipy):

    python trajClustering.py help
    python trajClustering.py cluster -top top.pdb -traj traj.xtc \\
        -sel "protein and name CA" -ncluster 5 -outdir traj_clusters

From another Python script:

    from trajClustering import load_universe, cluster_trajectory
    u = load_universe("top.pdb", "traj.xtc")
    labels, matrix = cluster_trajectory(u, sel="protein and name CA", ncluster=5)
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np


HELP = """
TrajCluster (Python) — RMSD clustering of one selection

Commands:
  rmsd      RMSD of each frame to a reference (default: first frame)
  matrix    Pairwise RMSD matrix after fitting the selection
  cluster   Hierarchical clustering from that matrix; write PDBs
  plot      Heatmap of a written RMSD matrix (.dat from VMD or Python)
  help      Show this message

Common options:
  -top <file>       topology (PDB, GRO, TPR, ...)
  -traj <file>      trajectory (XTC, DCD, TRR, ...)
  -sel <text>       MDAnalysis selection (default: protein and name CA)
  -first <i>        first frame (default: 0)
  -last <i|end>     last frame (default: end)
  -stride <n>       take every n-th frame (default: 1)
  -out <file>       output file for rmsd / matrix / plot

cluster also accepts one of:
  -ncluster <k>     cut the tree into k groups
  -cutoff <A>       cut where average-linkage RMSD exceeds this value
  -outdir <dir>     output directory (default: traj_clusters)
  -method <name>    linkage method: average|complete|single (default: average)

plot options:
  -matrix <file>    RMSD matrix .dat (required for plot)
  -out <file>       PNG path (default: rmsd_heatmap.png)
  -vmax <A>         color scale upper limit (default: auto)
  -cmap <name>      matplotlib colormap (default: viridis)

Examples:
  python trajClustering.py rmsd -top top.pdb -traj traj.xtc -sel "backbone" -out rmsd.dat
  python trajClustering.py matrix -top top.pdb -traj traj.xtc -sel "name CA" -out matrix.dat
  python trajClustering.py cluster -top top.pdb -traj traj.xtc -sel "name CA" -ncluster 4
  python trajClustering.py plot -matrix traj_clusters/rmsd_matrix.dat -out heatmap.png
"""


def _require_mda():
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis.rms import rmsd
    except ImportError as exc:
        raise SystemExit(
            "MDAnalysis is required. Activate an env that has it, e.g.\n"
            "  conda activate ComBio && pip install MDAnalysis scipy"
        ) from exc
    return mda, rmsd


def load_universe(top, traj=None, **kwargs):
    """Load a topology and optional trajectory."""
    mda, _ = _require_mda()
    if traj is None:
        return mda.Universe(top, **kwargs)
    return mda.Universe(top, traj, **kwargs)


def frame_indices(u, first=0, last=None, stride=1):
    n = len(u.trajectory)
    if last is None or last == "end":
        last = n - 1
    first = int(first)
    last = int(last)
    stride = int(stride)
    if first < 0 or last >= n or first > last or stride < 1:
        raise ValueError(f"bad frame range first={first} last={last} stride={stride} (n={n})")
    return list(range(first, last + 1, stride))


def cal_rmsd(u, sel="protein and name CA", ref_frame=0, frames=None):
    """RMSD of each frame to ref_frame after fitting the selection."""
    _, rmsd = _require_mda()
    atoms = u.select_atoms(sel)
    if atoms.n_atoms < 3:
        raise ValueError(f"selection has {atoms.n_atoms} atoms; need at least 3")
    if frames is None:
        frames = list(range(len(u.trajectory)))
    u.trajectory[ref_frame]
    ref = atoms.positions.copy()
    values = np.zeros(len(frames), dtype=float)
    for i, frame in enumerate(frames):
        u.trajectory[frame]
        values[i] = rmsd(atoms.positions, ref, superposition=True)
    return values


def cal_rmsd_matrix(u, sel="protein and name CA", frames=None):
    """
    Pairwise RMSD matrix for one selection.

    For each pair (i, j) the selection in frame i is fitted onto frame j,
    then the RMSD is computed. The matrix is symmetric.
    """
    _, rmsd = _require_mda()
    atoms = u.select_atoms(sel)
    if atoms.n_atoms < 3:
        raise ValueError(f"selection has {atoms.n_atoms} atoms; need at least 3")
    if frames is None:
        frames = list(range(len(u.trajectory)))
    n = len(frames)
    coords = []
    for frame in frames:
        u.trajectory[frame]
        coords.append(atoms.positions.copy())
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            value = rmsd(coords[i], coords[j], superposition=True)
            matrix[i, j] = value
            matrix[j, i] = value
        if i == 0 or (i + 1) % 25 == 0 or i == n - 1:
            print(f"TrajCluster: RMSD row {i + 1}/{n}", flush=True)
    return matrix


def hierarchical_cluster(matrix, ncluster=None, cutoff=None, method="average"):
    """Average/complete/single linkage on a square RMSD matrix."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    if (ncluster is None) == (cutoff is None):
        raise ValueError("specify exactly one of ncluster or cutoff")
    if method not in {"average", "complete", "single"}:
        raise ValueError("method must be average, complete, or single")
    condensed = squareform(matrix, checks=False)
    tree = linkage(condensed, method=method)
    if ncluster is not None:
        labels = fcluster(tree, t=int(ncluster), criterion="maxclust")
    else:
        labels = fcluster(tree, t=float(cutoff), criterion="distance")
    return labels, tree


def write_rmsd_table(path, frames, values):
    path = Path(path)
    with path.open("w") as handle:
        handle.write("# TrajCluster RMSD to reference\n")
        handle.write("# pandas: read_csv(path, sep=r'\\s+', comment='#')\n")
        handle.write("frame rmsd_A\n")
        for frame, value in zip(frames, values):
            handle.write(f"{frame} {value:.4f}\n")
    print(f"TrajCluster: wrote {path}")


def write_matrix(path, frames, matrix):
    path = Path(path)
    n = len(frames)
    with path.open("w") as handle:
        handle.write("# TrajCluster pairwise RMSD matrix\n")
        handle.write(f"# frame_order: {' '.join(str(f) for f in frames)}\n")
        handle.write("# pandas: read_csv(path, sep=r'\\s+', comment='#', index_col=0)\n")
        handle.write("frame " + " ".join(str(f) for f in frames) + "\n")
        for i in range(n):
            row = " ".join(f"{matrix[i, j]:.4f}" for j in range(n))
            handle.write(f"{frames[i]} {row}\n")
    print(f"TrajCluster: wrote {path} ({n} x {n})")


def read_matrix(path):
    """Read a TrajCluster RMSD matrix .dat (VMD or Python). Returns frames, matrix."""
    path = Path(path)
    frames = []
    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0] == "frame":
                continue
            frames.append(int(float(parts[0])))
            rows.append([float(x) for x in parts[1:]])
    if not rows:
        raise ValueError(f"no matrix data in {path}")
    matrix = np.asarray(rows, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"expected square matrix in {path}, got shape {matrix.shape}")
    if len(frames) != matrix.shape[0]:
        raise ValueError(f"frame labels ({len(frames)}) != matrix size ({matrix.shape[0]})")
    return frames, matrix


def plot_rmsd_heatmap(path, outfile="rmsd_heatmap.png", vmax=None, cmap="viridis"):
    """Plot a heatmap from an RMSD matrix .dat file."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required for plot. Activate ComBio, e.g.\n"
            "  conda activate ComBio"
        ) from exc

    frames, matrix = read_matrix(path)
    n = len(frames)
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(
        matrix,
        origin="lower",
        aspect="equal",
        interpolation="nearest",
        cmap=cmap,
        vmin=0.0,
        vmax=None if vmax is None else float(vmax),
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("RMSD (Å)")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Frame")
    ax.set_title(f"Pairwise RMSD ({n} frames)")

    # Tick every ~8–12 labels so the axis stays readable.
    step = max(1, n // 10)
    ticks = list(range(0, n, step))
    if ticks[-1] != n - 1:
        ticks.append(n - 1)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([str(frames[i]) for i in ticks], rotation=45, ha="right")
    ax.set_yticklabels([str(frames[i]) for i in ticks])

    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)
    print(f"TrajCluster: wrote {outfile}")
    return outfile


def write_clusters(u, frames, labels, matrix, outdir, sel="protein and name CA"):
    """Write summary, membership, matrix, and one PDB per frame under clusterN/."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for old in outdir.glob("cluster*"):
        if old.is_dir() and old.name.startswith("cluster"):
            for child in old.iterdir():
                child.unlink()
            old.rmdir()

    write_matrix(outdir / "rmsd_matrix.dat", frames, matrix)

    # Rank clusters by size: cluster1 = largest.
    groups = {}
    for idx, label in enumerate(labels):
        groups.setdefault(int(label), []).append(idx)
    ranked = sorted(groups.values(), key=len, reverse=True)

    summary = outdir / "summary.dat"
    membership = outdir / "frames.dat"
    with summary.open("w") as sum_h, membership.open("w") as mem_h:
        sum_h.write("# TrajCluster RMSD hierarchical clusters\n")
        sum_h.write(f"# sel: {sel}\n")
        sum_h.write(f"# ncluster: {len(ranked)}\n")
        sum_h.write("# pandas: read_csv(path, sep=r'\\s+', comment='#')\n")
        sum_h.write("cluster n_frames center_frame mean_rmsd_to_center_A\n")
        mem_h.write("# TrajCluster membership\n")
        mem_h.write("# pandas: read_csv(path, sep=r'\\s+', comment='#')\n")
        mem_h.write("frame cluster rmsd_to_center_A\n")

        atoms = u.atoms
        for cid, members in enumerate(ranked, start=1):
            folder = outdir / f"cluster{cid}"
            folder.mkdir(exist_ok=True)
            center = min(
                members,
                key=lambda i: float(np.mean([matrix[i, j] for j in members])),
            )
            mean_rmsd = float(np.mean([matrix[center, j] for j in members]))
            sum_h.write(
                f"{cid} {len(members)} {frames[center]} {mean_rmsd:.4f}\n"
            )
            for member in sorted(members):
                frame = frames[member]
                mem_h.write(
                    f"{frame} {cid} {matrix[member, center]:.4f}\n"
                )
                u.trajectory[frame]
                atoms.write(str(folder / f"frame_{frame}.pdb"))
            print(
                f"TrajCluster: cluster{cid} has {len(members)} frames, "
                f"center frame {frames[center]}"
            )
    print(f"TrajCluster: wrote PDBs under {outdir}/cluster1 .. cluster{len(ranked)}")
    print(f"TrajCluster: wrote {summary} and {membership}")
    return ranked


def cluster_trajectory(
    u,
    sel="protein and name CA",
    ncluster=None,
    cutoff=None,
    method="average",
    first=0,
    last=None,
    stride=1,
    outdir=None,
    write_pdbs=True,
):
    """
    Cluster frames by pairwise fitted RMSD of one selection.

    Returns
    -------
    labels : ndarray of cluster ids (1-based, ranked by size if outdir is set)
    matrix : pairwise RMSD matrix aligned with the analyzed frames
    frames : list of trajectory frame indices
    """
    frames = frame_indices(u, first=first, last=last, stride=stride)
    print(f"TrajCluster: building RMSD matrix for {len(frames)} frames")
    matrix = cal_rmsd_matrix(u, sel=sel, frames=frames)
    labels, _tree = hierarchical_cluster(
        matrix, ncluster=ncluster, cutoff=cutoff, method=method
    )
    print(f"TrajCluster: identified {len(set(labels))} clusters")
    if outdir and write_pdbs:
        write_clusters(u, frames, labels, matrix, outdir, sel=sel)
    return labels, matrix, frames


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="trajClustering.py",
        description="RMSD clustering of one trajectory selection",
        add_help=False,
    )
    parser.add_argument("command", nargs="?", default="help")
    parser.add_argument("-top", default=None)
    parser.add_argument("-traj", default=None)
    parser.add_argument("-sel", default="protein and name CA")
    parser.add_argument("-first", default=0)
    parser.add_argument("-last", default="end")
    parser.add_argument("-stride", type=int, default=1)
    parser.add_argument("-out", default=None)
    parser.add_argument("-outdir", default="traj_clusters")
    parser.add_argument("-ncluster", type=int, default=None)
    parser.add_argument("-cutoff", type=float, default=None)
    parser.add_argument("-method", default="average")
    parser.add_argument("-matrix", default=None)
    parser.add_argument("-vmax", type=float, default=None)
    parser.add_argument("-cmap", default="viridis")
    parser.add_argument("-h", "--help", action="store_true")
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    if args.command in {"help", "-h", "--help"} or args.help:
        print(HELP)
        return 0

    if args.command == "plot":
        if args.matrix is None:
            print("TrajCluster plot: -matrix <rmsd_matrix.dat> is required", file=sys.stderr)
            return 1
        plot_rmsd_heatmap(
            args.matrix,
            outfile=args.out or "rmsd_heatmap.png",
            vmax=args.vmax,
            cmap=args.cmap,
        )
        return 0

    if args.top is None:
        print("TrajCluster: -top is required", file=sys.stderr)
        print(HELP)
        return 1

    last = None if args.last in (None, "end") else int(args.last)
    u = load_universe(args.top, args.traj)
    frames = frame_indices(u, first=args.first, last=last, stride=args.stride)

    if args.command == "rmsd":
        values = cal_rmsd(u, sel=args.sel, ref_frame=frames[0], frames=frames)
        out = args.out or "traj_rmsd.dat"
        write_rmsd_table(out, frames, values)
    elif args.command == "matrix":
        matrix = cal_rmsd_matrix(u, sel=args.sel, frames=frames)
        out = args.out or "traj_rmsd_matrix.dat"
        write_matrix(out, frames, matrix)
    elif args.command == "cluster":
        if (args.ncluster is None) == (args.cutoff is None):
            print("TrajCluster: set -ncluster or -cutoff", file=sys.stderr)
            return 1
        cluster_trajectory(
            u,
            sel=args.sel,
            ncluster=args.ncluster,
            cutoff=args.cutoff,
            method=args.method,
            first=args.first,
            last=last,
            stride=args.stride,
            outdir=args.outdir,
        )
    else:
        print(f"TrajCluster: unknown command '{args.command}'", file=sys.stderr)
        print(HELP)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
