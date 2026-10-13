# alphafold3

Open any `.py` file — the top docstring has purpose, when to use, and examples.

| File | Role | Main entry |
|------|------|------------|
| `inputAF3.py` | Build AF3 **input** JSON | `GenerateAf3InputJson` |
| `extractAF3.py` | Pull MSA + templates from finished job `*_data.json` | `ExtractMsaTemplate` |
| `plddtAF3.py` | Mean CA pLDDT over predicted CIFs | `PlotCAPlddt` |
| `msaAF3.py` | Pack / convert / jackhmmer MSA for custom inputs | `CustomMsa` |

Structures: `utils.parseFiles.LoadAf3`.

```python
from machine_learning.alphafold3 import (
    GenerateAf3InputJson, ExtractMsaTemplate, PlotCAPlddt,
)

ExtractMsaTemplate("./af3_job", outDir="./af3_features", chainId="A")
```

```bash
python -m machine_learning.alphafold3.extractAF3 extract -dir ./af3_job -out ./af3_features -chain A
python -m machine_learning.alphafold3.plddtAF3 -dir ./af3_job -sel "protein and name CA" -out plddt.png
```
