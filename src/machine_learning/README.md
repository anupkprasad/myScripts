# Machine learning (`src/machine_learning`)

| Subdir | Contents |
|--------|----------|
| `alphafold3/` | AF3 input JSON, extract MSA/templates, CA pLDDT, MSA helpers |
| `pyrosetta/` | Minimization / binding-energy helpers, ligand params |

## AlphaFold3

```python
from machine_learning.alphafold3 import (
    GenerateAf3InputJson, ExtractMsaTemplate, PlotCAPlddt,
)
from utils.parseFiles import LoadAf3

GenerateAf3InputJson(name="job", modelSeeds=[0], sequence="ACDE", outputFile="job.json")
ExtractMsaTemplate("./af3_job", outDir="./af3_features", chainId="A")
u = LoadAf3("./af3_job")
PlotCAPlddt(u, sel="protein and name CA and chainID A", outFile="plddt.png")
```

```bash
python -m machine_learning.alphafold3.extractAF3 extract -dir ./af3_job -out ./af3_features -chain A
python -m machine_learning.alphafold3.plddtAF3 -dir ./af3_job \
  -sel "protein and name CA and chainID A" -out plddt.png
```

See `docs/tutorials/03_af3_and_sequence.md` and `alphafold3/README.md`.

PyRosetta requires a licensed PyRosetta install in the active environment.
