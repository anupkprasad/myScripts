#!/usr/bin/env python3
"""
extractAF3 — pull MSA + structural templates from a finished AF3 job.

When to use
    AF3 already ran. You want the searched MSA and/or template CIFs that were
    used as *inputs* to the network (not the predicted models).

Where data lives
    <job>/<name>_data.json
        sequences[].protein.unpairedMsa / pairedMsa
        sequences[].protein.templates[].mmcif + query/templateIndices

Public API
    ExtractMsaTemplate(job, outDir, chainId="A")  ← usual one-call entry
    ExtractMsa(...)       MSA only (.a3m)
    ExtractTemplates(...) template_XX.cif + chainX_templates.json
    FindAf3DataJson(job)  locate *_data.json

Outputs (same outDir for ExtractMsaTemplate)
    chainA_unpaired.a3m / chainA_paired.a3m
    chainA_template_00.cif …
    chainA_templates.json   (index arrays kept on one line each)

Related modules
    plddtAF3 / LoadAf3 → analyze *predicted* seed-*_sample-*/model.cif
    inputAF3           → build a new input JSON

Example (Python)
    from machine_learning.alphafold3 import ExtractMsaTemplate
    ExtractMsaTemplate("/path/to/af3_job", outDir="af3_features", chainId="A")

Example (CLI)
    python -m machine_learning.alphafold3.extractAF3 extract \\
        -dir ./af3_job -out ./af3_features -chain A
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

PathLike = Union[str, Path]

HELP = """
extractAF3 — pull MSA + templates from AF3 *_data.json

Commands:
  extract     MSA and templates into one directory (default)
  msa         MSA only
  templates   Templates only
  help        Show this message

Examples:
  python -m machine_learning.alphafold3.extractAF3 extract -dir ./af3_job -out ./af3_features -chain A
  python -m machine_learning.alphafold3.extractAF3 msa -dir ./af3_job -out ./af3_features -chain A
"""


def FindAf3DataJson(directory: PathLike) -> Path:
    """Locate ``*_data.json`` (or a bare JSON path) for an AF3 job."""
    path = Path(directory).expanduser().resolve()
    if path.is_file():
        return path
    if not path.is_dir():
        raise FileNotFoundError(f"not found: {path}")
    matches = sorted(path.glob("*_data.json"))
    if not matches:
        matches = sorted(path.glob("*data.json"))
    if not matches:
        raise FileNotFoundError(f"no *_data.json in {path}")
    return matches[0]


def _loadAf3Json(dataJson: PathLike) -> dict:
    path = FindAf3DataJson(dataJson)
    with open(path, "r") as fh:
        return json.load(fh)


def _proteinByChain(jdata: dict, chainId: str) -> dict:
    for block in jdata.get("sequences", []):
        protein = block.get("protein")
        if not protein:
            continue
        cid = protein.get("id")
        ids = [cid] if isinstance(cid, str) else list(cid or [])
        if chainId in ids or (len(ids) == 1 and ids[0] == chainId):
            return protein
    raise KeyError(f"no protein chain '{chainId}' in AF3 JSON sequences")


def _cleanMsaLines(
    msaText: str,
    noLowercase: bool = False,
    noGap: bool = False,
) -> List[str]:
    lines: List[str] = []
    for line in (msaText or "").split("\n"):
        if not line:
            continue
        if line.startswith(">"):
            lines.append(line)
            continue
        if noLowercase:
            line = "".join(c for c in line if not c.islower())
        if noGap:
            line = line.replace("-", "")
        lines.append(line)
    return lines


def _compactJsonArray(values) -> str:
    """JSON array on one line, e.g. [0,1,2]."""
    return json.dumps(values, separators=(",", ":"))


def _writeTemplatesMeta(metaPath: Path, meta: List[dict]) -> None:
    """Write templates.json with query/template index arrays each on one line."""
    parts = ["["]
    for i, entry in enumerate(meta):
        comma = "," if i < len(meta) - 1 else ""
        block = (
            "  {\n"
            f'    "index": {int(entry["index"])},\n'
            f'    "cif": {json.dumps(entry["cif"])},\n'
            f'    "queryIndices": {_compactJsonArray(entry.get("queryIndices") or [])},\n'
            f'    "templateIndices": {_compactJsonArray(entry.get("templateIndices") or [])}\n'
            f"  }}{comma}"
        )
        parts.append(block)
    parts.append("]")
    metaPath.write_text("\n".join(parts) + "\n")


def ExtractMsa(
    dataJsonOrDir: PathLike,
    outDir: PathLike,
    chainId: str = "A",
    kind: str = "both",
    noLowercase: bool = False,
    noGap: bool = False,
) -> Dict[str, Path]:
    """
    Extract MSA strings from AF3 ``*_data.json`` into ``outDir``.

    kind: ``unpaired``, ``paired``, or ``both``.
    Returns mapping of kind → written A3M path.
    """
    jdata = _loadAf3Json(dataJsonOrDir)
    protein = _proteinByChain(jdata, chainId)
    outPath = Path(outDir).expanduser().resolve()
    outPath.mkdir(parents=True, exist_ok=True)

    kinds = ["unpaired", "paired"] if kind == "both" else [kind]
    written: Dict[str, Path] = {}
    keyMap = {"unpaired": "unpairedMsa", "paired": "pairedMsa"}
    for k in kinds:
        if k not in keyMap:
            raise ValueError("kind must be unpaired, paired, or both")
        text = protein.get(keyMap[k]) or ""
        if not text:
            print(f"extractAF3: no {k} MSA for chain {chainId}")
            continue
        lines = _cleanMsaLines(text, noLowercase=noLowercase, noGap=noGap)
        outFile = outPath / f"chain{chainId}_{k}.a3m"
        with open(outFile, "w") as fh:
            fh.write("\n".join(lines) + ("\n" if lines else ""))
        written[k] = outFile
        nSeq = sum(1 for ln in lines if ln.startswith(">"))
        print(f"extractAF3: wrote {k} MSA ({nSeq} sequences) → {outFile}")
    if not written:
        raise ValueError(f"no MSA content for chain {chainId} (kind={kind})")
    return written


def ExtractTemplates(
    dataJsonOrDir: PathLike,
    outDir: PathLike,
    chainId: str = "A",
) -> List[Path]:
    """
    Write each AF3 template mmCIF from ``*_data.json`` into ``outDir``.

    Also writes ``chainX_templates.json`` (index arrays kept on one line each).
    """
    jdata = _loadAf3Json(dataJsonOrDir)
    protein = _proteinByChain(jdata, chainId)
    templates = protein.get("templates") or []
    outPath = Path(outDir).expanduser().resolve()
    outPath.mkdir(parents=True, exist_ok=True)
    if not templates:
        print(f"extractAF3: no templates for chain {chainId}")
        return []

    written: List[Path] = []
    meta: List[dict] = []
    for i, tmpl in enumerate(templates):
        mmcif = tmpl.get("mmcif") if isinstance(tmpl, dict) else None
        if not mmcif:
            continue
        cifPath = outPath / f"chain{chainId}_template_{i:02d}.cif"
        with open(cifPath, "w") as fh:
            fh.write(mmcif if mmcif.endswith("\n") else mmcif + "\n")
        written.append(cifPath)
        meta.append(
            {
                "index": i,
                "cif": cifPath.name,
                "queryIndices": tmpl.get("queryIndices"),
                "templateIndices": tmpl.get("templateIndices"),
            }
        )
        print(f"extractAF3: wrote template {i} → {cifPath}")

    metaPath = outPath / f"chain{chainId}_templates.json"
    _writeTemplatesMeta(metaPath, meta)
    print(f"extractAF3: wrote template index map → {metaPath}")
    return written


def ExtractMsaTemplate(
    dataJsonOrDir: PathLike,
    outDir: PathLike,
    chainId: str = "A",
    kind: str = "both",
    noLowercase: bool = False,
    noGap: bool = False,
) -> Dict[str, object]:
    """
    Extract MSA and templates into the **same** ``outDir`` in one call.

    Returns ``{"msa": {...}, "templates": [...], "outDir": Path}``.
    """
    outPath = Path(outDir).expanduser().resolve()
    outPath.mkdir(parents=True, exist_ok=True)
    msaWritten: Dict[str, Path] = {}
    try:
        msaWritten = ExtractMsa(
            dataJsonOrDir,
            outDir=outPath,
            chainId=chainId,
            kind=kind,
            noLowercase=noLowercase,
            noGap=noGap,
        )
    except ValueError as exc:
        print(f"extractAF3: MSA skipped ({exc})")
    tmplWritten = ExtractTemplates(dataJsonOrDir, outDir=outPath, chainId=chainId)
    if not msaWritten and not tmplWritten:
        raise ValueError(f"nothing extracted for chain {chainId} from {dataJsonOrDir}")
    print(f"extractAF3: MSA + templates → {outPath}")
    return {"msa": msaWritten, "templates": tmplWritten, "outDir": outPath}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("command", nargs="?", default="help")
    parser.add_argument("-dir", default=".")
    parser.add_argument("-out", default=None)
    parser.add_argument("-chain", default="A")
    parser.add_argument("-kind", default="both", choices=["unpaired", "paired", "both"])
    parser.add_argument("-nolowcase", action="store_true")
    parser.add_argument("-nogap", action="store_true")
    parser.add_argument("-h", "--help", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.command in {"help", "-h", "--help"} or args.help:
        print(HELP)
        return 0

    if args.command in {"extract", "msa-template", "msa_template", "both"}:
        outDir = args.out or "af3_features"
        ExtractMsaTemplate(
            args.dir,
            outDir=outDir,
            chainId=args.chain,
            kind=args.kind,
            noLowercase=args.nolowcase,
            noGap=args.nogap,
        )
        return 0

    if args.command == "msa":
        outDir = args.out or "af3_features"
        ExtractMsa(
            args.dir,
            outDir=outDir,
            chainId=args.chain,
            kind=args.kind,
            noLowercase=args.nolowcase,
            noGap=args.nogap,
        )
        return 0

    if args.command == "templates":
        outDir = args.out or "af3_features"
        ExtractTemplates(args.dir, outDir=outDir, chainId=args.chain)
        return 0

    print(f"extractAF3: unknown command '{args.command}'", file=sys.stderr)
    print(HELP)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
