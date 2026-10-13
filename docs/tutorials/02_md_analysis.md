# Tutorial 02 — MD analysis

## InterDomain (VMD)

Two-domain geometry (COM distance, axis angle, dihedral), contacts, and clustering.

```tcl
# After callScripts.tcl is loaded
InterDomain help
InterDomain -sel1 "resid 1 to 100" -sel2 "resid 200 to 300" -out geometry.dat
```

Plot with:

```bash
conda activate ComBio
python ~/myScripts/src/md/analysis/interDomain_plot.py --help
```

## TrajCluster (VMD + Python)

Pairwise fit + RMSD hierarchical clustering of one selection.

**VMD:**

```tcl
TrajCluster help
TrajCluster cluster -sel "protein and name CA" -ncluster 4 -outdir traj_clusters
# or
TrajCluster cluster -sel "name CA" -cutoff 2.5 -stride 5
```

Writes `outdir/rmsd_matrix.dat`, `summary.dat`, `frames.dat`, and `clusterN/*.pdb`.

**Heatmap of the matrix:**

```bash
python ~/myScripts/src/md/analysis/trajClustering.py plot \
  -matrix traj_clusters/rmsd_matrix.dat -out heatmap.png
```

**Python-only** (no VMD; needs MDAnalysis + scipy):

```bash
traj-cluster cluster -top top.pdb -traj traj.xtc \
  -sel "protein and name CA" -ncluster 4 -outdir traj_clusters
```

## Secondary structure (`secStrModule.py`)

From a trajectory (MDAnalysis DSSP):

```bash
python ~/myScripts/src/md/analysis/secStrModule.py from-traj \
  -top top.pdb -traj traj.xtc \
  -sel "protein and resid 1-200" -stride 5 \
  -matrix ss.dat -out ss_heatmap.png
```

From a VMD letter file (`ss_traj.dat`):

```bash
python ~/myScripts/src/md/analysis/secStrModule.py from-vmd \
  -file ss_traj.dat -resNum 17 -out ss_vmd.png
```

With an existing Universe in Python:

```python
from md.analysis.secStrModule import computeSecStr, plotSecStrHeatmap
frames, resids, codes = computeSecStr(universe=u, sel="protein")
plotSecStrHeatmap(codes, frames, resids, outFile="ss.png")
```

## GROMACS quick setup

```bash
bash ~/myScripts/src/md/setup/gmx/simSetup.sh protein.pdb
# Force field copied from data/forcefields/gmx/charmm36-jul2022.ff
```

Place MDP templates in `src/md/setup/gmx/mdp/` if your pipeline expects them.
