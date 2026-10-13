# MD (`src/md`)

GROMACS + NAMD setup, VMD visualization/analysis, and Python trajectory tools.

| Subdir | Contents |
|--------|----------|
| `setup/` | Topology helpers, GMX shell pipelines (`setup/gmx/`), NAMD/VMD system build (`setup/namd/`) |
| `analysis/` | InterDomain, TrajCluster, secStr heatmap, RMSF, gyration, waters, … (Tcl + Python) |
| `visualization/` | VMD representations, render, play/rotate |
| `enhanced_sampling/` | ABF / PMF Python post-processing |
| `vmd_autoload.tcl` | Recursive Tcl loader used by `~/myScripts/callScripts.tcl` |

Force fields: `~/myScripts/data/forcefields/{gmx,namd}/`  
Reference PDBs: `~/myScripts/data/reference/KinaseStructures/`

## Most used

```tcl
# VMD
InterDomain help
TrajCluster cluster -sel "protein and name CA" -ncluster 4 -outdir traj_clusters
```

```bash
python src/md/analysis/trajClustering.py plot -matrix traj_clusters/rmsd_matrix.dat -out hm.png
python src/md/analysis/secStrModule.py from-traj -top top.pdb -traj traj.xtc -sel "protein" -out ss.png
python src/md/analysis/secStrModule.py from-vmd -file ss_traj.dat -resNum 17 -out ss_vmd.png
bash src/md/setup/gmx/simSetup.sh my.pdb
```
