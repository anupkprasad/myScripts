# Tutorial 03 — AlphaFold3 and sequence tools

## Naming

Files: **lowerCamelCase**. Public functions: **PascalCase** (see `docs/convention.md`).

## Load AF3 models (VMD + Python)

AF3 writes folders like `seed-<seed>_sample-<0..4>/` with a `.cif` in each. pLDDT is in the CIF B-factor column.

**VMD** (`utils/parseFiles.tcl`):

```tcl
LoadStruct help
LoadAf3 cif -rename af3_ensemble
TrajCluster rmsd -sel "protein and name CA" -out af3_rmsd.dat
LoadAf3 multimol -prefix DCLK3_ -ref DCLK3_KD.pdb
```

**Python / MDAnalysis** (`utils/parseFiles.py` — same stem as the Tcl file):

```python
from utils.parseFiles import LoadAf3
u = LoadAf3("/path/to/af3_job")
print(len(u.trajectory), "frames")
```

## Extract MSA / templates / pLDDT

MSA and templates live in the job’s `*_data.json`. Models are under `seed-*_sample-*/model.cif`.

```python
from machine_learning.alphafold3 import ExtractMsaTemplate, PlotCAPlddt
from utils.parseFiles import LoadAf3

job = "/path/to/af3_job"
ExtractMsaTemplate(job, outDir="af3_features", chainId="A")  # MSA + templates, one dir

u = LoadAf3(job)
PlotCAPlddt(u, sel="protein and name CA and chainID A", outFile="plddt.png")
```

CLI:

```bash
python -m machine_learning.alphafold3.extractAF3 extract -dir ./af3_job -out ./af3_features -chain A
python -m machine_learning.alphafold3.plddtAF3 -dir ./af3_job \
  -sel "protein and name CA and chainID A" -out plddt.png
```

## AlphaFold3 JSON prep

Module: `src/machine_learning/alphafold3/inputAF3.py`

```python
from machine_learning.alphafold3 import GenerateAf3InputJson

GenerateAf3InputJson(
    name="my_job",
    modelSeeds=[1, 2, 3, 4, 5],
    sequence="MKTFF...",
    outputFile="af3_input.json",
    chainId="A",
    ptms=[{"type": "SEP", "position": 12}],  # 1-based
)
```

Related helpers:

- `msaAF3.py` — `CustomMsa`, jackhmmer, MMA→A3M
- `extractAF3.py` — `ExtractMsaTemplate`
- `plddtAF3.py` — CA pLDDT

## Sequence / MSA

Under `src/bioinformatics/sequence/` (see module top docstrings):

| File | Role |
|------|------|
| `msaIO.py` | `LoadMsa` / format convert (FASTA, A3M, CMA) |
| `msaStats.py` | conservation / entropy / gap profiles |
| `msaClean.py` | `CleanMsa` + QC plots |
| `msaView.py` | query-centric display / Dash viewer |
| `seqCluster.py` | HMMER domtblout → UMAP |
| `uniProt.py` | UniProt REST helpers |

Full examples: `src/bioinformatics/sequence/README.md`.

```python
from bioinformatics.sequence import LoadMsa, PlotConservation
records = LoadMsa("af3_features/chainA_unpaired.a3m", fmt="a3m")
PlotConservation(records[:200], outFile="conservation.png")  # subset if huge
```

## Reference structures

Kinase KD (± tail) PDBs:

```text
~/myScripts/data/reference/KinaseStructures/
```

In Python:

```python
from myscripts.paths import REF
print(list(REF.glob("*.pdb"))[:5])
```
