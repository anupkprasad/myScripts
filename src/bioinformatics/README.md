# Bioinformatics (`src/bioinformatics`)

| Subdir | Role | Guide |
|--------|------|--------|
| `sequence/` | MSA I/O, cleaning, conservation, viewing, clustering, UniProt | [`sequence/README.md`](sequence/README.md) |
| `structure/` | Structure helpers | [`structure/README.md`](structure/README.md) |
| `structure/pymol/` | PyMOL figure recipes | [`structure/pymol/README.md`](structure/pymol/README.md) |
| `structure/hexDock/` | Hex docking | [`structure/hexDock/README.md`](structure/hexDock/README.md) |

Naming: **lowerCamelCase** files, **PascalCase** public APIs — see `docs/convention.md`.
Each module’s top docstring also has purpose + examples.

## Quick start (sequence)

```python
from bioinformatics.sequence import LoadMsa, PlotMsa, PlotPhyloTree

records = LoadMsa("clean.fasta")
PlotMsa(records, outFile="msa.png", maxSeqs=40, queryId="query")
PlotPhyloTree(records, outFile="tree.png", method="nj", maxSeqs=80)
```

```bash
python -m bioinformatics.sequence.msaStats aln.fasta -plot conservation -out stats_out
```

AF3 job MSA/templates → `machine_learning.alphafold3.ExtractMsaTemplate`, then
`LoadMsa(..., fmt="a3m")` here. Details: [`sequence/README.md`](sequence/README.md).
