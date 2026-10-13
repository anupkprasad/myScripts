"""
alphafold3 — AlphaFold3 helpers (input → extract → pLDDT / MSA tools).

Modules (open the file for full examples)
    inputAF3.py    Build AF3 input JSON          GenerateAf3InputJson
    extractAF3.py  MSA + templates from job      ExtractMsaTemplate
    plddtAF3.py    Mean CA pLDDT over models     PlotCAPlddt / MeanCAPlddt
    msaAF3.py      Custom MSA pack / convert      CustomMsa, MmaToA3m, BuildAf3Msa

Typical flow
    1) GenerateAf3InputJson(...)     → run AF3 externally
    2) ExtractMsaTemplate(job, ...)  → MSA + template CIFs from *_data.json
    3) LoadAf3(job) + PlotCAPlddt    → ensemble pLDDT (utils.parseFiles)

Quick import
    from machine_learning.alphafold3 import (
        GenerateAf3InputJson, ExtractMsaTemplate, PlotCAPlddt, CustomMsa,
    )
"""

from .extractAF3 import (
    ExtractMsa,
    ExtractMsaTemplate,
    ExtractTemplates,
    FindAf3DataJson,
)
from .inputAF3 import GenerateAf3InputJson, GenerateMultichainAf3InputJson
from .msaAF3 import (
    BuildAf3Msa,
    CustomMsa,
    FastaToA3m,
    MmaToA3m,
    SeqModification,
)
from .plddtAF3 import MeanCAPlddt, PlotCAPlddt

__all__ = [
    "GenerateAf3InputJson",
    "GenerateMultichainAf3InputJson",
    "ExtractMsaTemplate",
    "ExtractMsa",
    "ExtractTemplates",
    "FindAf3DataJson",
    "MeanCAPlddt",
    "PlotCAPlddt",
    "CustomMsa",
    "SeqModification",
    "MmaToA3m",
    "BuildAf3Msa",
    "FastaToA3m",
]
