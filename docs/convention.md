# Naming & coding conventions (myScripts)

One style for **Python, Tcl, and bash**. Use this for all new code; rename older files when you touch them.

---

## 1. File names — lowerCamelCase

- Start with a **lowercase** letter.
- Capitalize the start of each following word (no underscores, no kebab-case).

| Good | Bad |
|------|-----|
| `parseFiles.tcl` | `parse_files.tcl` |
| `parseFiles.py` | `load_struct.py` |
| `trajClustering.py` | `traj_clustering.py` |
| `secStrModule.py` | `sec_str_module.py` |
| `basicFunctions.tcl` | `BasicFunctions.tcl` |
| `inputAF3.py` | `AF3_input_prep.py` |

### Same job → same base name

If Tcl / Python / bash do the **same task**, use the **same stem**; only the extension changes:

```text
utils/parseFiles.tcl     ← VMD
utils/parseFiles.py      ← MDAnalysis
utils/parseFiles.sh      ← optional shell wrapper (if needed)
```

You should see the pair and know they belong together.

---

## 2. Functions / procs — PascalCase for public APIs

| Kind | Style | Examples |
|------|--------|----------|
| **Public API** (import / VMD / CLI entry users call) | **PascalCase** | `LoadStruct`, `LoadAf3`, `LoadUniverse`, `ExtractMsa`, `TrajCluster` |
| **Private helpers** | `_` + lowerCamelCase | `_lsArgs`, `_requireMda`, `_af3CifPaths` |

```python
# Python — public
def LoadUniverse(top, traj=None): ...
def LoadAf3(directory=".", pattern="seed"): ...

# Python — private
def _af3CifPaths(directory, pattern="seed"): ...
```

```tcl
# Tcl — public tool
proc LoadStruct {args} { ... }
proc LoadAf3 {args} { ... }
# Tcl — helper
proc _lsAf3Frames {dir pattern} { ... }
```

```bash
# bash — public entry
LoadAf3() { ... }
```

---

## 3. Arguments / variables

Use **lowerCamelCase** (not `snake_case`):

| Good | Bad |
|------|-----|
| `topPath`, `trajFile`, `dirPath` | `top_path`, `traj_file` |
| `nAtoms`, `seedDirs` | `n_atoms`, `seed_dirs` |
| `-dir`, `-pattern` (CLI flags stay short) | |

Python:

```python
def LoadAf3(directory=".", pattern="seed"):
    seedDirs = ...
```

Tcl options stay `-key value` for the user; internal variables use lowerCamelCase:

```tcl
set seedDirs [_lsSeedDirs $dir $pattern]
```

---

## 4. Where code lives

Organize by **research field**, not by language.

| Location | Use for |
|----------|---------|
| `src/md/` | MD setup, analysis, visualization, enhanced sampling |
| `src/bioinformatics/` | Sequence / MSA / PyMOL / docking |
| `src/machine_learning/` | AlphaFold3, PyRosetta |
| `utils/` | Shared loaders & helpers (`parseFiles.*`, plotting, …) |
| `scripts/` | Thin CLI wrappers only |
| `data/` | Reference structures and force fields |
| `docs/` | Tutorials and this file |

Immediate-run VMD snippets → `src/md/analysis/oneshot/` (not auto-loaded).

---

## 5. Public tool pattern (Tcl)

```tcl
proc ToolName {args} {
    if {[llength $args] == 0 || [lindex $args 0] in {help -h --help}} {
        _toolHelp
        return
    }
    ...
}
puts "ToolName loaded. Type: ToolName help"
```

Rules:

1. `ToolName help` always works.
2. Options: `-key value`.
3. No `mol new` / side effects at `source` time.
4. Paths via `::myscripts::ROOT` / `refStructures`, not hardcoded `/home/...`.

---

## 6. Autoload

- `src/md/`, `src/bioinformatics/`, and top-level `utils/*.tcl` load via `vmd_autoload.tcl`.
- New analysis tool → `src/md/analysis/myTool.tcl`.
- Structure loading → extend `utils/parseFiles.tcl` **and** `utils/parseFiles.py` together.
- Do **not** edit `callScripts.tcl` for every new file.

---

## 7. Quick checklist

- [ ] File name is lowerCamelCase (`parseFiles.py`)
- [ ] Matching stems if multi-language (`parseFiles.tcl` + `parseFiles.py`)
- [ ] Public functions / VMD entries are PascalCase + `help` where applicable
- [ ] Args / locals are lowerCamelCase
- [ ] Private helpers start with `_`
- [ ] No side effects on source / import
