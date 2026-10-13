# Tutorial 01 — Install and VMD autoload

## 1. Python environment

```bash
conda activate ComBio
cd ~/myScripts
pip install -e .
```

Check paths:

```bash
python -c "from myscripts.paths import ROOT, DATA, REF; print(ROOT); print(REF)"
```

## 2. VMD hook

Ensure `~/.vmdrc` contains:

```tcl
source ~/myScripts/callScripts.tcl
```

Start VMD. You should see something like:

```text
myScripts VMD: loaded N Tcl modules from src/md (+ bioinformatics)
  ROOT=/home/you/myScripts
  DATA=.../data  REF=.../KinaseStructures
Custom myScripts functions are integrated
```

Globals available: `pathVmdScript`, `pathAnalysis`, `pathSetup`, `pathGraphic`, `refStructures`, `::myscripts::ROOT`, `::myscripts::DATA`.

## 3. Add your own VMD tool

1. Create `~/myScripts/src/md/analysis/myTool.tcl` with a proc, e.g. `proc MyTool {args} { ... }`.
2. Restart VMD (or `source ~/myScripts/callScripts.tcl` again).
3. Call `MyTool` in the Tk console.

No need to edit `callScripts.tcl`.

## 4. CLI without install

```bash
python ~/myScripts/src/md/analysis/trajClustering.py help
~/myScripts/scripts/secstr_heatmap help
```
