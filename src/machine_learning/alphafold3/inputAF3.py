#!/usr/bin/env python3
"""
inputAF3 — build AlphaFold3 *input* JSON before running a prediction.

When to use
    You have a sequence (and optional PTMs / ligand / custom MSA) and need a
    job JSON that AF3 can read. Official fields:
    https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md

Public API
    GenerateAf3InputJson(...)           single chain (± ligand)
    GenerateMultichainAf3InputJson(...) multi-chain complex

Related modules
    msaAF3.CustomMsa     → pack a local A3M into unpairedMsa string
    extractAF3           → pull MSA/templates *after* a finished AF3 job
    utils.parseFiles.LoadAf3 → load predicted CIFs

Example
    from machine_learning.alphafold3 import GenerateAf3InputJson, CustomMsa

    # Auto MSA + template search (leave MSA/templates as None)
    GenerateAf3InputJson(
        name="job",
        modelSeeds=[0, 100, 200],
        sequence="MSSIGTGYDLSASTF",
        outputFile="job.json",
        ptms=[{"type": "SEP", "position": 11}],  # 1-based
        chainId="A",
    )

    # Custom MSA, template-free
    seq, oneLineMsa = CustomMsa("my.a3m")
    GenerateAf3InputJson(
        name="job_custom_msa",
        modelSeeds=[0],
        sequence=seq,
        outputFile="job_custom.json",
        unpairedMsa=oneLineMsa,
        pairedMsa="",
        templates=[],
    )
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

PathLike = Union[str, Path]


def _formatPtms(ptmList: Optional[Sequence[dict]]) -> List[dict]:
    if not ptmList:
        return []
    return [{"ptmType": ptm["type"], "ptmPosition": ptm["position"]} for ptm in ptmList]


def GenerateAf3InputJson(
    name: str,
    modelSeeds: List[int],
    sequence: str,
    outputFile: PathLike,
    ptms: Optional[List[dict]] = None,
    unpairedMsa: Optional[str] = None,
    pairedMsa: Optional[str] = None,
    templates: Optional[List] = None,
    chainId: str = "A",
    ligand: Optional[str] = None,
):
    """
    Write a single-chain (optional ligand) AlphaFold3 input JSON.

    PTMs: ``[{"type": "SEP", "position": 11}, ...]`` (1-based).
    Templates: ``None`` = AF3 search; ``[]`` = template-free; or custom list.
    """
    ptmBlock = _formatPtms(ptms or [])
    protein = {
        "id": chainId,
        "sequence": sequence,
        "modifications": ptmBlock,
        "unpairedMsa": unpairedMsa,
        "pairedMsa": pairedMsa,
        "templates": templates,
    }
    sequences = [{"protein": protein}]
    if ligand:
        sequences.extend(
            [
                {"ligand": {"id": ["L"], "ccdCodes": [ligand]}},
                {"ligand": {"id": ["I"], "ccdCodes": ["MG"]}},
                {"ligand": {"id": ["J"], "ccdCodes": ["MG"]}},
            ]
        )

    af3Input = {
        "dialect": "alphafold3",
        "version": 1,
        "name": name,
        "sequences": sequences,
        "modelSeeds": modelSeeds,
    }
    outPath = Path(outputFile).expanduser().resolve()
    outPath.parent.mkdir(parents=True, exist_ok=True)
    with open(outPath, "w") as fh:
        json.dump(af3Input, fh, indent=2)
    print(f"inputAF3: wrote {outPath}")
    return outPath


def GenerateMultichainAf3InputJson(
    name: str,
    modelSeeds: List[int],
    sequence: List[str],
    outputFile: PathLike,
    ptms: Optional[List] = None,
    unpairedMsa: Optional[str] = None,
    pairedMsa: Optional[str] = None,
    templates: Optional[List] = None,
    chainId: Optional[List[str]] = None,
):
    """
    Write a multi-chain AlphaFold3 input JSON.

    ptms: per-chain list-of-lists, or a flat list applied to every chain (legacy).
    MSA fields are attached to the first chain only (AF3 convention for complexes).
    """
    if chainId is None:
        chainId = [chr(ord("A") + i) for i in range(len(sequence))]
    if len(chainId) != len(sequence):
        raise ValueError("chainId length must match sequence length")

    if ptms is None:
        ptmBlocks = [[] for _ in sequence]
    elif ptms and isinstance(ptms[0], dict):
        formatted = _formatPtms(ptms)
        ptmBlocks = [formatted for _ in sequence]
    else:
        ptmBlocks = [_formatPtms(chainPtms) for chainPtms in ptms]

    if len(ptmBlocks) != len(sequence):
        raise ValueError(
            f"ptms length ({len(ptmBlocks)}) must match sequence length ({len(sequence)})"
        )

    sequencesBlock = []
    for i, (cid, seq) in enumerate(zip(chainId, sequence)):
        protein: Dict = {
            "id": cid,
            "sequence": seq,
            "modifications": ptmBlocks[i],
            "templates": templates,
        }
        if i == 0:
            protein["unpairedMsa"] = unpairedMsa
            protein["pairedMsa"] = pairedMsa
        sequencesBlock.append({"protein": protein})

    af3Input = {
        "dialect": "alphafold3",
        "version": 1,
        "name": name,
        "sequences": sequencesBlock,
        "modelSeeds": modelSeeds,
    }
    outPath = Path(outputFile).expanduser().resolve()
    outPath.parent.mkdir(parents=True, exist_ok=True)
    with open(outPath, "w") as fh:
        json.dump(af3Input, fh, indent=2)
    print(f"inputAF3: wrote {outPath}")
    return outPath


def _mapLegacyKwargs(kwargs: dict) -> dict:
    alias = {
        "model_seeds": "modelSeeds",
        "output_file": "outputFile",
        "PTMs": "ptms",
        "unpaired_msa": "unpairedMsa",
        "paired_msa": "pairedMsa",
        "chainID": "chainId",
    }
    out = dict(kwargs)
    for old, new in alias.items():
        if old in out and new not in out:
            out[new] = out.pop(old)
    return out


def generateAF3InputJson(*args, **kwargs):
    return GenerateAf3InputJson(*args, **_mapLegacyKwargs(kwargs))


def generateMultichainAF3InputJson(*args, **kwargs):
    return GenerateMultichainAf3InputJson(*args, **_mapLegacyKwargs(kwargs))


USAGE = """
inputAF3 — build AlphaFold3 input JSON

from machine_learning.alphafold3 import GenerateAf3InputJson, CustomMsa

GenerateAf3InputJson(
    name="job_single_chain",
    modelSeeds=[0, 100, 200],
    sequence="MSSIGTGYDLSASTF",
    outputFile="single_chain.json",
    ptms=[{"type": "SEP", "position": 11}],
    chainId="A",
)

sequence, oneLineMsa = CustomMsa("my_alignment.a3m")
"""


if __name__ == "__main__":
    print(USAGE)
