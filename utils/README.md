# Utils

Shared helpers used across research fields.

Naming: see [docs/convention.md](../docs/convention.md) — **lowerCamelCase** files, **PascalCase** public functions; same stem across languages.

| Path | Role |
|------|------|
| `parseFiles.tcl` | VMD — **`LoadStruct help`** / **`LoadAf3 help`** |
| `parseFiles.py` | MDAnalysis — `LoadUniverse`, `LoadAf3`, … (same job as the Tcl file) |
| `plotStyle.py`, `basicFunction.py` | Plotting / small helpers |
| `imageMovie/` | PNG stacks → GIF/video |
| `arialFonts/` | Bundled fonts |

## VMD

```tcl
LoadStruct help
LoadStruct pdb
LoadStruct xtc -top ../equilibration/min.gro -traj md.xtc
LoadStruct af3 -rename af3_ensemble
LoadAf3 cif
```

## Python

```python
from utils.parseFiles import LoadUniverse, LoadAf3, LoadCifs

u = LoadUniverse("top.pdb", "traj.xtc")
u = LoadAf3("/path/to/af3_job")
```

```bash
python ~/myScripts/utils/parseFiles.py help
python ~/myScripts/utils/parseFiles.py af3 -dir ./af3_job
```
