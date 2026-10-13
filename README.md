# myScripts

Personal computational-biology toolkit for **MD** (GROMACS, NAMD, VMD), **bioinformatics** (sequence / structure), and **machine learning** (AlphaFold3, PyRosetta).

Code is organized by **research field**, not by language. Python, Tcl, and bash live side by side under the same folder.

## Layout

```text
src/
  md/                 # setup, analysis, visualization, enhanced sampling
  bioinformatics/     # sequence + structure (PyMOL, HexDock)
  machine_learning/   # AlphaFold3 + PyRosetta
  myscripts/          # shared paths (ROOT, DATA, REF)
utils/                # plotting, image→movie helpers
scripts/              # thin CLI wrappers
docs/tutorials/       # how-to guides
data/
  reference/          # kinase PDBs and other reference structures
  forcefields/        # GMX + NAMD parameter sets
backup/               # local dumps (gitignored)
trash/                # unused vendors (e.g. old pconpy)
experimentals/        # one-off lab scripts
callScripts.tcl       # VMD entry → auto-loads src/md/**/*.tcl
```

## Install

```bash
conda activate ComBio   # or your env with MDAnalysis / matplotlib
cd ~/myScripts
pip install -e .
```

Optional: add `~/myScripts/scripts` to `PATH`.

```bash
export PATH="$HOME/myScripts/scripts:$PATH"
```

Override repo root if needed:

```bash
export MYSCRIPTS_ROOT=/path/to/myScripts
```

## VMD

In `~/.vmdrc` (already set in this machine’s [`.vmdrc`](.vmdrc)):

```tcl
source ~/myScripts/callScripts.tcl
```

That loads [`src/md/vmd_autoload.tcl`](src/md/vmd_autoload.tcl), which:

- sets `::myscripts::ROOT`, `::myscripts::DATA`, `refStructures`, `pathSetup`, …
- sources core Tcl (`basicFunctions`, `trajFunctions`, `utils/parseFiles.tcl`) first
- then auto-sources `src/md/`, `src/bioinformatics/`, and other `utils/*.tcl`

**Adding a new VMD command:** put `MyTool.tcl` in `src/md/analysis/` (or extend `utils/parseFiles.tcl` for loaders) and restart VMD. Naming: [docs/convention.md](docs/convention.md).

```tcl
LoadStruct help       ;# pdb / cif / dcd / xtc / AF3
LoadAf3 help          ;# AF3 shortcuts
TrajCluster help
InterDomain help
```

## Quick commands

```bash
# Pairwise RMSD clustering + heatmap
traj-cluster help
python src/md/analysis/trajClustering.py cluster -top top.pdb -traj traj.xtc \
  -sel "protein and name CA" -ncluster 4 -outdir traj_clusters
python src/md/analysis/trajClustering.py plot \
  -matrix traj_clusters/rmsd_matrix.dat -out heatmap.png

# Secondary-structure heatmap
secstr-module from-traj -top top.pdb -traj traj.xtc -sel "protein" -out ss.png
secstr-module from-vmd -file ss_traj.dat -resNum 17 -out ss_vmd.png

# From Python
python -c "from myscripts.paths import ROOT, DATA, REF; print(ROOT, REF)"
```

In VMD Tk console:

```tcl
InterDomain help
TrajCluster help
```

## Docs

| Guide | Topic |
|-------|--------|
| [docs/tutorials/01_install_and_vmd.md](docs/tutorials/01_install_and_vmd.md) | Install + VMD autoload |
| [docs/tutorials/02_md_analysis.md](docs/tutorials/02_md_analysis.md) | InterDomain, TrajCluster, secStr |
| [docs/tutorials/03_af3_and_sequence.md](docs/tutorials/03_af3_and_sequence.md) | AF3 prep + sequence tools |
| [src/md/README.md](src/md/README.md) | MD field map |
| [src/bioinformatics/README.md](src/bioinformatics/README.md) | Sequence / PyMOL / docking |
| [src/machine_learning/README.md](src/machine_learning/README.md) | AF3 / PyRosetta |
| [utils/README.md](utils/README.md) | Shared utilities |

## License

MIT — see repository history / project convention.
