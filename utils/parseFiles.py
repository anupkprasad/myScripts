#!/usr/bin/env python3
"""
Load structures / trajectories into an MDAnalysis Universe.

Same role as utils/parseFiles.tcl (LoadStruct / LoadAf3 in VMD).
Same base name on purpose: parseFiles.py ↔ parseFiles.tcl.

    from utils.parseFiles import LoadUniverse, LoadAf3, LoadPdbs, LoadCifs

    u = LoadUniverse("top.pdb", "traj.xtc")
    u = LoadAf3("/path/to/af3_job")
    u = LoadCifs(directory=".")

CLI:

    python parseFiles.py help
    python parseFiles.py af3 -dir ./af3_job
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Union

import numpy as np

PathLike = Union[str, Path]


HELP = """
parseFiles — MDAnalysis loaders (companion to parseFiles.tcl)

Commands:
  universe   Load topology + optional trajectory (pdb/gro/psf + xtc/dcd/…)
  pdb        Load one or more PDBs (multi-file → multi-frame)
  cif        Load one or more CIFs as frames
  af3        Load AF3 seed-*_sample-* CIFs as one multi-frame Universe
  help       Show this message

Examples:
  python parseFiles.py universe -top top.pdb -traj traj.xtc
  python parseFiles.py af3 -dir ./af3_job
  python -c "from utils.parseFiles import LoadAf3; print(len(LoadAf3('.').trajectory))"
"""


def _requireMda():
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise SystemExit(
            "MDAnalysis is required. Activate ComBio, e.g.\n"
            "  conda activate ComBio"
        ) from exc
    return mda


def LoadUniverse(
    top: PathLike,
    traj: Optional[Union[PathLike, Sequence[PathLike]]] = None,
    **kwargs,
):
    """Load a topology and optional trajectory into an MDAnalysis Universe."""
    mda = _requireMda()
    topPath = Path(top).expanduser().resolve()
    if not topPath.is_file():
        raise FileNotFoundError(f"topology not found: {topPath}")
    if traj is None:
        return mda.Universe(str(topPath), **kwargs)
    if isinstance(traj, (str, Path)):
        trajPaths = [Path(traj).expanduser().resolve()]
    else:
        trajPaths = [Path(t).expanduser().resolve() for t in traj]
    for trajPath in trajPaths:
        if not trajPath.is_file():
            raise FileNotFoundError(f"trajectory not found: {trajPath}")
    return mda.Universe(str(topPath), *[str(t) for t in trajPaths], **kwargs)


def _globFiles(directory: PathLike, pattern: str) -> List[Path]:
    directory = Path(directory).expanduser().resolve()
    files = sorted(directory.glob(pattern), key=lambda p: p.name)
    return [p for p in files if p.is_file()]


def _isCif(path: Path) -> bool:
    return path.suffix.lower() in {".cif", ".mmcif"}


def _requireBioPdb():
    try:
        from Bio.PDB import MMCIFParser, PDBIO, PDBParser
    except ImportError as exc:
        raise SystemExit(
            "Biopython is required to read mmCIF (MDAnalysis has no CIF topology reader).\n"
            "  conda activate ComBio && pip install biopython"
        ) from exc
    return MMCIFParser, PDBIO, PDBParser


def _cifToTempPdb(cifPath: Path, pdbPath: Path) -> Path:
    """Convert mmCIF → PDB for MDAnalysis (which cannot parse CIF topologies)."""
    MMCIFParser, PDBIO, _ = _requireBioPdb()
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(cifPath.stem, str(cifPath))
    io = PDBIO()
    io.set_structure(structure)
    pdbPath.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(pdbPath))
    return pdbPath


def _coordsFromPath(path: Path) -> np.ndarray:
    """Atomic coordinates (nAtoms, 3) from PDB or CIF."""
    if _isCif(path):
        MMCIFParser, _, _ = _requireBioPdb()
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(path.stem, str(path))
        return np.asarray([a.coord for a in structure.get_atoms()], dtype=np.float32)
    mda = _requireMda()
    u = mda.Universe(str(path))
    return np.asarray(u.atoms.positions, dtype=np.float32)


def _bfactorsFromPath(path: Path) -> np.ndarray:
    """Per-atom B-factors / pLDDT (nAtoms,) from PDB or CIF."""
    if _isCif(path):
        MMCIFParser, _, _ = _requireBioPdb()
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(path.stem, str(path))
        return np.asarray([a.bfactor for a in structure.get_atoms()], dtype=np.float32)
    mda = _requireMda()
    u = mda.Universe(str(path))
    return np.asarray(u.atoms.tempfactors, dtype=np.float32)


def _universeFromStructureFile(path: Path):
    """Build an MDAnalysis Universe from PDB or CIF (CIF via temp PDB)."""
    import tempfile

    mda = _requireMda()
    if not _isCif(path):
        return mda.Universe(str(path))
    tmp = Path(tempfile.mkdtemp(prefix="parseFiles_cif_")) / f"{path.stem}.pdb"
    _cifToTempPdb(path, tmp)
    return mda.Universe(str(tmp))


def _coordsFromStructures(paths: Sequence[Path]) -> tuple:
    """Return (templateUniverse, coords array shape (nFrames, nAtoms, 3))."""
    if not paths:
        raise ValueError("no structure files to load")
    template = _universeFromStructureFile(paths[0])
    nAtoms = template.atoms.n_atoms
    frames = np.empty((len(paths), nAtoms, 3), dtype=np.float32)
    frames[0] = _coordsFromPath(paths[0])
    if frames[0].shape[0] != nAtoms:
        raise ValueError(
            f"atom count mismatch between topology and coords for {paths[0]}"
        )
    for i, path in enumerate(paths[1:], start=1):
        xyz = _coordsFromPath(path)
        if xyz.shape[0] != nAtoms:
            raise ValueError(
                f"atom count mismatch: {paths[0].name} has {nAtoms}, "
                f"{path.name} has {xyz.shape[0]}"
            )
        frames[i] = xyz
    return template, frames


def _universeFromFrames(template, frames: np.ndarray):
    """Attach in-memory coordinates as trajectory frames."""
    template.load_new(frames, order="fac")
    return template


def LoadPdbs(
    files: Optional[Iterable[PathLike]] = None,
    directory: PathLike = ".",
):
    """Load PDB(s). Multiple files → multi-frame Universe (same topology)."""
    if files is None:
        paths = _globFiles(directory, "*.pdb")
    else:
        paths = [Path(f).expanduser().resolve() for f in files]
    if not paths:
        raise FileNotFoundError(f"no PDB files in {directory}")
    if len(paths) == 1:
        return LoadUniverse(paths[0])
    template, frames = _coordsFromStructures(paths)
    return _universeFromFrames(template, frames)


def LoadCifs(
    files: Optional[Iterable[PathLike]] = None,
    directory: PathLike = ".",
):
    """Load CIF(s). Multiple files → multi-frame Universe (via Biopython→PDB)."""
    if files is None:
        paths = _globFiles(directory, "*.cif")
    else:
        paths = [Path(f).expanduser().resolve() for f in files]
    if not paths:
        raise FileNotFoundError(f"no CIF files in {directory}")
    template, frames = _coordsFromStructures(paths)
    if len(paths) == 1:
        return template
    return _universeFromFrames(template, frames)


def Af3CifPaths(directory: PathLike = ".", pattern: str = "seed") -> List[Path]:
    """Collect AF3 sample CIFs under seed-* / seed-*_sample-* directories."""
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"not a directory: {directory}")
    seedDirs = sorted(
        [p for p in directory.glob(f"{pattern}*") if p.is_dir()],
        key=lambda p: p.name,
    )
    if not seedDirs:
        raise FileNotFoundError(f"no directories matching '{pattern}*' in {directory}")
    cifs: List[Path] = []
    for seedDir in seedDirs:
        found = sorted(seedDir.glob("*.cif"), key=lambda p: p.name)
        if not found:
            raise FileNotFoundError(f"no .cif in {seedDir}")
        cifs.append(found[0])
    return cifs


def LoadAf3(directory: PathLike = ".", pattern: str = "seed"):
    """
    Load an AlphaFold3 job folder as one multi-frame Universe.

    AF3 mmCIF is converted with Biopython (MDAnalysis cannot read CIF topologies).
    Same behaviour as ``LoadStruct af3`` / ``LoadAf3 cif`` in VMD.
    Per-residue pLDDT is stored as B-factor / tempfactor on each CIF atom.
    """
    paths = Af3CifPaths(directory, pattern=pattern)
    print(f"parseFiles: loading {len(paths)} AF3 CIF models (CIF→PDB for MDAnalysis)…")
    template, frames = _coordsFromStructures(paths)
    u = _universeFromFrames(template, frames)
    u._af3CifPaths = paths  # type: ignore[attr-defined]
    print(f"parseFiles: AF3 ensemble {len(paths)} frames from {directory}")
    return u


def LoadXtc(top: PathLike, traj: PathLike, **kwargs):
    """Load GRO/PDB topology + XTC trajectory."""
    return LoadUniverse(top, traj, **kwargs)


def LoadDcd(top: PathLike, traj: Union[PathLike, Sequence[PathLike]], **kwargs):
    """Load PSF/PDB topology + DCD trajectory (one or more files)."""
    return LoadUniverse(top, traj, **kwargs)


# Backward-compatible aliases (older camelCase / snake_case names)
loadUniverse = LoadUniverse
loadAf3 = LoadAf3
loadPdbs = LoadPdbs
loadCifs = LoadCifs
loadXtc = LoadXtc
loadDcd = LoadDcd
af3CifPaths = Af3CifPaths


def load_universe(*args, **kwargs):
    return LoadUniverse(*args, **kwargs)


def load_af3(*args, **kwargs):
    return LoadAf3(*args, **kwargs)


def load_pdbs(*args, **kwargs):
    return LoadPdbs(*args, **kwargs)


def load_cifs(*args, **kwargs):
    return LoadCifs(*args, **kwargs)


def main(argv=None):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("command", nargs="?", default="help")
    parser.add_argument("-top", default=None)
    parser.add_argument("-traj", default=None)
    parser.add_argument("-dir", default=".")
    parser.add_argument("-pattern", default="seed")
    parser.add_argument("-h", "--help", action="store_true")
    args = parser.parse_args(argv)

    if args.command in {"help", "-h", "--help"} or args.help:
        print(HELP)
        return 0

    if args.command == "universe":
        if args.top is None:
            print("parseFiles universe: -top is required", file=sys.stderr)
            return 1
        u = LoadUniverse(args.top, args.traj)
        print(f"Universe: {u.atoms.n_atoms} atoms, {len(u.trajectory)} frames")
        return 0

    if args.command == "pdb":
        u = LoadPdbs(directory=args.dir)
        print(f"Universe: {u.atoms.n_atoms} atoms, {len(u.trajectory)} frames")
        return 0

    if args.command == "cif":
        u = LoadCifs(directory=args.dir)
        print(f"Universe: {u.atoms.n_atoms} atoms, {len(u.trajectory)} frames")
        return 0

    if args.command == "af3":
        u = LoadAf3(args.dir, pattern=args.pattern)
        print(f"Universe: {u.atoms.n_atoms} atoms, {len(u.trajectory)} frames")
        return 0

    print(f"parseFiles: unknown command '{args.command}'", file=sys.stderr)
    print(HELP)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
