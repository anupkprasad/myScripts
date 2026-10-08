#!/usr/bin/env python3
"""
Plot tables written by InterDomain.tcl.

Each data file has '#' comment lines and one whitespace-separated header.
Import the functions, or run this file from the command line.

    python3 interDomain_plot.py com-angle geometry.dat -o com_vs_angle.png
    python3 interDomain_plot.py clusters geomet_clusters/frames.dat -o clusters.png

    from interDomain_plot import plot_clusters
    plot_clusters("geomet_clusters/frames.dat", outfile="clusters.png")
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


LABELS = {
    "frame": "Frame",
    "cluster": "Cluster",
    "distance_A": "COM–COM distance (Å)",
    "axis_angle_deg": "Principal-axis angle (deg)",
    "dihedral_deg": "Dihedral (deg)",
    "com1_x": "Domain 1 COM x (Å)",
    "com1_y": "Domain 1 COM y (Å)",
    "com1_z": "Domain 1 COM z (Å)",
    "com2_x": "Domain 2 COM x (Å)",
    "com2_y": "Domain 2 COM y (Å)",
    "com2_z": "Domain 2 COM z (Å)",
    "occupancy": "Occupancy",
    "min_distance_A": "Minimum contact distance (Å)",
    "n_frames": "Frames in contact",
    "n_atom_contacts": "Atom contacts",
    "n_residue_pairs": "Residue pairs",
    "n_res_sel1": "sel1 residues at the interface",
    "n_res_sel2": "sel2 residues at the interface",
    "rmsd_to_center_A": "RMSD to center (Å)",
}

COM_COLUMNS = ("com1_x", "com1_y", "com1_z", "com2_x", "com2_y", "com2_z")

FEATURE_ALIASES = {
    "com": "distance_A",
    "distance": "distance_A",
    "distance_a": "distance_A",
    "angle": "axis_angle_deg",
    "axis_angle": "axis_angle_deg",
    "axis_angle_deg": "axis_angle_deg",
    "dihedral": "dihedral_deg",
    "dihedral_deg": "dihedral_deg",
}

DEFAULT_CLUSTER_AXES = ("distance_A", "axis_angle_deg", "dihedral_deg")


def read_table(path):
    """Read one InterDomain .dat file. Comment lines start with '#'."""
    table = pd.read_csv(path, sep=r"\s+", comment="#")
    if "frame" in table.columns:
        table = table.sort_values("frame").reset_index(drop=True)
    return table


def read_comment_meta(path):
    """Return key/value pairs from '# key: value' comment lines."""
    meta = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line.startswith("#"):
                break
            body = line[1:].strip()
            if ":" not in body:
                continue
            key, value = body.split(":", 1)
            meta[key.strip().lower()] = value.strip()
    return meta


def label_of(column):
    return LABELS.get(column, column.replace("_", " "))


def resolve_feature_names(names):
    """Map aliases like com/angle/dihedral to frames.dat column names."""
    if isinstance(names, str):
        names = [tok.strip() for tok in names.replace(",", " ").split() if tok.strip()]
    columns = []
    for name in names:
        key = name.lower()
        if key not in FEATURE_ALIASES:
            raise ValueError(
                f"unknown feature '{name}'. Use com, angle, dihedral "
                "or distance_A, axis_angle_deg, dihedral_deg"
            )
        column = FEATURE_ALIASES[key]
        if column not in columns:
            columns.append(column)
    return columns


def add_com_distance(table):
    """Add distance_A from the two center-of-mass columns when it is absent."""
    missing = [name for name in COM_COLUMNS if name not in table.columns]
    if missing:
        raise ValueError(
            "COM distance needs distance_A, or columns " + ", ".join(COM_COLUMNS)
        )
    out = table.copy()
    dx = out["com2_x"] - out["com1_x"]
    dy = out["com2_y"] - out["com1_y"]
    dz = out["com2_z"] - out["com1_z"]
    out["distance_A"] = (dx * dx + dy * dy + dz * dz) ** 0.5
    return out


def _as_table(data):
    if isinstance(data, pd.DataFrame):
        return data.copy()
    return read_table(data)


def combine_on_frame(left, right):
    """Merge two InterDomain tables on the frame column."""
    left_table = _as_table(left)
    right_table = _as_table(right)
    if "frame" not in left_table.columns or "frame" not in right_table.columns:
        raise ValueError("both tables need a frame column")
    overlap = [
        name for name in right_table.columns
        if name != "frame" and name in left_table.columns
    ]
    right_table = right_table.drop(columns=overlap)
    return left_table.merge(right_table, on="frame", how="inner")


# Distinct colors for cluster groups (easy to tell apart in 2D/3D plots).
CLUSTER_PALETTE = [
    "#e41a1c",  # red
    "#377eb8",  # blue
    "#ff7f00",  # orange
    "#984ea3",  # purple
    "#00bfc4",  # cyan
    "#4daf4a",  # green
    "#a65628",  # brown
    "#f781bf",  # magenta
    "#ffff33",  # yellow
    "#999999",  # gray
]


def _cluster_colors(clusters):
    ids = sorted(pd.unique(clusters))
    colors = {}
    for i, cid in enumerate(ids):
        colors[cid] = CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)]
    return colors, ids


def _scatter_by_cluster_2d(ax, data, x, y):
    colors, ids = _cluster_colors(data["cluster"])
    for cid in ids:
        part = data[data["cluster"] == cid]
        ax.scatter(
            part[x],
            part[y],
            s=28,
            linewidths=0.3,
            edgecolors="k",
            color=colors[cid],
            label=f"cluster {int(cid)} (n={len(part)})",
        )
    ax.set_xlabel(label_of(x))
    ax.set_ylabel(label_of(y))
    ax.legend(frameon=False, fontsize=8)


def _scatter_by_cluster_3d(ax, data, x, y, z):
    colors, ids = _cluster_colors(data["cluster"])
    for cid in ids:
        part = data[data["cluster"] == cid]
        ax.scatter(
            part[x],
            part[y],
            part[z],
            s=28,
            linewidths=0.3,
            edgecolors="k",
            color=colors[cid],
            label=f"cluster {int(cid)} (n={len(part)})",
            depthshade=True,
        )
    ax.set_xlabel(label_of(x))
    ax.set_ylabel(label_of(y))
    ax.set_zlabel(label_of(z))
    ax.legend(frameon=False, fontsize=8, loc="upper left")


def choose_cluster_axes(data, features=None, path=None):
    """
    Decide which columns to plot for a geometcluster frames.dat file.

    Priority:
      1. explicit features / column names
      2. '# features:' comment in the file
      3. all available among distance_A, axis_angle_deg, dihedral_deg
    """
    if features:
        return resolve_feature_names(features)

    if path is not None:
        meta = read_comment_meta(path)
        if "features" in meta and meta["features"]:
            try:
                return resolve_feature_names(meta["features"])
            except ValueError:
                pass

    axes = [name for name in DEFAULT_CLUSTER_AXES if name in data.columns]
    if not axes:
        raise ValueError(
            "no geometry axes found; need distance_A, axis_angle_deg, and/or dihedral_deg"
        )
    return axes


def plot_clusters(
    data,
    features=None,
    outfile=None,
    title=None,
    also_timeline=False,
):
    """
    Visualize geometcluster frames.dat.

    Points are colored by the cluster column.
    1 feature  -> strip / jitter plot
    2 features -> 2D scatter
    3 features -> 3D scatter
    Set also_timeline=True to add a frame-vs-cluster panel.
    """
    path = None if isinstance(data, pd.DataFrame) else data
    table = _as_table(data)
    if "cluster" not in table.columns:
        raise ValueError("frames.dat needs a cluster column")
    axes = choose_cluster_axes(table, features=features, path=path)
    missing = [name for name in axes if name not in table.columns]
    if missing:
        raise ValueError(f"missing columns for clustering plot: {missing}")

    n_axes = len(axes)
    if also_timeline:
        fig = plt.figure(figsize=(11.5, 4.8) if n_axes < 3 else (12.5, 5.2))
        if n_axes == 3:
            ax = fig.add_subplot(1, 2, 1, projection="3d")
        else:
            ax = fig.add_subplot(1, 2, 1)
        ax_time = fig.add_subplot(1, 2, 2)
    else:
        fig = plt.figure(figsize=(6.5, 5.2) if n_axes < 3 else (7.2, 5.8))
        ax = fig.add_subplot(111, projection="3d" if n_axes == 3 else None)
        ax_time = None

    if n_axes == 1:
        x = axes[0]
        colors, ids = _cluster_colors(table["cluster"])
        rng = np.random.default_rng(0)
        for cid in ids:
            part = table[table["cluster"] == cid]
            jitter = rng.normal(0.0, 0.08, size=len(part))
            ax.scatter(
                part["cluster"].to_numpy(dtype=float) + jitter,
                part[x],
                s=28,
                linewidths=0.3,
                edgecolors="k",
                color=colors[cid],
                label=f"cluster {int(cid)} (n={len(part)})",
            )
        ax.set_xlabel(label_of("cluster"))
        ax.set_ylabel(label_of(x))
        ax.set_xticks(ids)
        ax.legend(frameon=False, fontsize=8)
    elif n_axes == 2:
        _scatter_by_cluster_2d(ax, table, axes[0], axes[1])
    else:
        _scatter_by_cluster_3d(ax, table, axes[0], axes[1], axes[2])

    ax.set_title(title or f"Geometry clusters ({', '.join(axes)})")

    if ax_time is not None:
        colors, ids = _cluster_colors(table["cluster"])
        for cid in ids:
            part = table[table["cluster"] == cid]
            ax_time.scatter(
                part["frame"],
                part["cluster"],
                s=22,
                color=colors[cid],
                label=f"cluster {int(cid)}",
            )
        ax_time.set_xlabel(label_of("frame"))
        ax_time.set_ylabel(label_of("cluster"))
        ax_time.set_yticks(ids)
        ax_time.set_title("Cluster membership over frames")
        ax_time.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    _finish(fig, outfile)
    return ax


def cluster_sizes(summary, outfile=None, title=None):
    """Bar chart of cluster populations from geometcluster summary.dat."""
    data = _as_table(summary)
    if "cluster" not in data.columns or "n_frames" not in data.columns:
        raise ValueError("summary.dat needs cluster and n_frames columns")
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    colors, ids = _cluster_colors(data["cluster"])
    ax.bar(
        [str(int(c)) for c in data["cluster"]],
        data["n_frames"],
        color=[colors[c] for c in data["cluster"]],
        edgecolor="k",
        linewidth=0.4,
    )
    ax.set_xlabel(label_of("cluster"))
    ax.set_ylabel("Number of frames")
    ax.set_title(title or "Cluster sizes")
    fig.tight_layout()
    _finish(fig, outfile)
    return ax


def scatter(table, x, y, outfile=None, color="frame", title=None, ax=None):
    """
    Scatter plot of two columns.

    color="frame" colors points by time. color="cluster" uses discrete groups.
    color=None uses one color.
    """
    data = _as_table(table)
    for name in (x, y):
        if name not in data.columns:
            raise ValueError(f"column '{name}' not in table: {list(data.columns)}")
    created = ax is None
    if created:
        _, ax = plt.subplots(figsize=(6.2, 5.0))
    if color == "cluster" and "cluster" in data.columns:
        _scatter_by_cluster_2d(ax, data, x, y)
    elif color and color in data.columns:
        points = ax.scatter(
            data[x], data[y], c=data[color], cmap="viridis", s=18, linewidths=0
        )
        bar = ax.figure.colorbar(points, ax=ax)
        bar.set_label(label_of(color))
    else:
        ax.scatter(data[x], data[y], s=18, linewidths=0, color="#1f4e79")
    ax.set_xlabel(label_of(x))
    ax.set_ylabel(label_of(y))
    ax.set_title(title or f"{label_of(y)} vs {label_of(x)}")
    if created:
        ax.figure.tight_layout()
        _finish(ax.figure, outfile)
    return ax


def scatter_com_angle(data, angle=None, outfile=None, color="frame", title=None):
    """
    Scatter of COM–COM distance against the angle between principal axes.

    data can be a geometry file (it already has both columns), a COM file plus
    an angle file, or a DataFrame.
    """
    table = _as_table(data)
    if angle is not None:
        table = combine_on_frame(table, angle)
    if "distance_A" not in table.columns:
        table = add_com_distance(table)
    if "axis_angle_deg" not in table.columns:
        raise ValueError("axis_angle_deg is missing; pass a geometry file or --angle")
    return scatter(
        table,
        x="distance_A",
        y="axis_angle_deg",
        outfile=outfile,
        color=color,
        title=title or "COM distance vs inter-domain angle",
    )


def timeseries(table, columns, outfile=None, title=None):
    """Plot one or more columns against frame. columns is a name or a list of names."""
    data = _as_table(table)
    if isinstance(columns, str):
        columns = [columns]
    if "frame" not in data.columns:
        raise ValueError("table has no frame column")
    missing = [name for name in columns if name not in data.columns]
    if missing:
        raise ValueError(f"missing columns: {missing}")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for name in columns:
        ax.plot(data["frame"], data[name], label=label_of(name), linewidth=1.2)
    ax.set_xlabel(label_of("frame"))
    ax.set_ylabel(label_of(columns[0]) if len(columns) == 1 else "Value")
    if len(columns) > 1:
        ax.legend()
    ax.set_title(title or "Inter-domain time series")
    fig.tight_layout()
    _finish(fig, outfile)
    return ax


def interface_occupancy(table, outfile=None, top=30, min_occupancy=0.0):
    """Horizontal bar chart of residue-pair occupancy."""
    data = _as_table(table)
    required = ("resname1", "resid1", "resname2", "resid2", "occupancy")
    missing = [name for name in required if name not in data.columns]
    if missing:
        raise ValueError(f"interface table is missing {missing}")
    data = data[data["occupancy"] >= min_occupancy].copy()
    data = data.sort_values("occupancy", ascending=False).head(top)
    data = data.iloc[::-1]
    labels = [
        f"{rn1}{int(ri1)}–{rn2}{int(ri2)}"
        for rn1, ri1, rn2, ri2 in zip(
            data["resname1"], data["resid1"], data["resname2"], data["resid2"]
        )
    ]
    height = max(3.5, 0.28 * len(data) + 1.0)
    fig, ax = plt.subplots(figsize=(7.0, height))
    ax.barh(labels, data["occupancy"], color="#1f4e79")
    ax.set_xlabel(label_of("occupancy"))
    ax.set_xlim(0, 1)
    ax.set_title(f"Interface residue pairs (top {len(data)})")
    fig.tight_layout()
    _finish(fig, outfile)
    return ax


def interface_frames(table, column="n_residue_pairs", outfile=None):
    """Time series of one interface-per-frame column."""
    return timeseries(table, column, outfile=outfile, title="Interface over time")


def _finish(fig, outfile):
    if outfile:
        path = Path(outfile)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=160)
        print(f"wrote {path}")
        plt.close(fig)
    else:
        plt.show()


def _build_parser():
    import argparse

    parser = argparse.ArgumentParser(
        description="Plot InterDomain.tcl output tables."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    com_angle = sub.add_parser(
        "com-angle",
        help="scatter COM–COM distance against the principal-axis angle",
    )
    com_angle.add_argument(
        "data",
        help="geometry file, or a COM/distance file when --angle is set",
    )
    com_angle.add_argument("--angle", help="angle file, if data is not a geometry file")
    com_angle.add_argument("-o", "--outfile", help="png path; omit to open a window")

    scatter_cmd = sub.add_parser("scatter", help="scatter any two columns")
    scatter_cmd.add_argument("data")
    scatter_cmd.add_argument("x")
    scatter_cmd.add_argument("y")
    scatter_cmd.add_argument(
        "--color",
        default="frame",
        help="column used for color (frame, cluster, none)",
    )
    scatter_cmd.add_argument("-o", "--outfile")

    clusters = sub.add_parser(
        "clusters",
        help="plot geometcluster frames.dat colored by cluster group",
    )
    clusters.add_argument(
        "data",
        help="geomet_clusters/frames.dat from InterDomain geometcluster",
    )
    clusters.add_argument(
        "--features",
        default=None,
        help='axes to plot, e.g. "com,angle" or "com,angle,dihedral"',
    )
    clusters.add_argument(
        "--timeline",
        action="store_true",
        help="also show a frame-vs-cluster membership panel",
    )
    clusters.add_argument(
        "--summary",
        default=None,
        help="optional summary.dat; also writes a cluster-size bar chart",
    )
    clusters.add_argument("-o", "--outfile")

    trace = sub.add_parser("trace", help="columns vs frame")
    trace.add_argument("data")
    trace.add_argument("columns", nargs="+")
    trace.add_argument("-o", "--outfile")

    interface = sub.add_parser("interface", help="occupancy bars for residue pairs")
    interface.add_argument("data")
    interface.add_argument("--top", type=int, default=30)
    interface.add_argument("--min-occupancy", type=float, default=0.0)
    interface.add_argument("-o", "--outfile")

    frames = sub.add_parser("interface-frames", help="interface counts vs frame")
    frames.add_argument("data")
    frames.add_argument("--column", default="n_residue_pairs")
    frames.add_argument("-o", "--outfile")
    return parser


def main(argv=None):
    args = _build_parser().parse_args(argv)
    if args.command == "com-angle":
        scatter_com_angle(args.data, angle=args.angle, outfile=args.outfile)
    elif args.command == "scatter":
        color = None if args.color in ("none", "None", "") else args.color
        scatter(args.data, args.x, args.y, outfile=args.outfile, color=color)
    elif args.command == "clusters":
        plot_clusters(
            args.data,
            features=args.features,
            outfile=args.outfile,
            also_timeline=args.timeline,
        )
        if args.summary:
            size_out = None
            if args.outfile:
                out = Path(args.outfile)
                size_out = out.with_name(out.stem + "_sizes" + out.suffix)
            cluster_sizes(args.summary, outfile=size_out)
    elif args.command == "trace":
        timeseries(args.data, args.columns, outfile=args.outfile)
    elif args.command == "interface":
        interface_occupancy(
            args.data,
            outfile=args.outfile,
            top=args.top,
            min_occupancy=args.min_occupancy,
        )
    elif args.command == "interface-frames":
        interface_frames(args.data, column=args.column, outfile=args.outfile)


if __name__ == "__main__":
    main()
