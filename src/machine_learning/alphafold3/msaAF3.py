#!/usr/bin/env python3
"""
msaAF3 — MSA utilities used when *preparing* custom inputs for AF3.

When to use
    You are building an MSA yourself (jackhmmer / MMA / local A3M) and need
    to convert or pack it for GenerateAf3InputJson(... unpairedMsa=...).

    Not for pulling MSA out of a finished job → use extractAF3 instead.

Public API
    CustomMsa(file)              A3M/FASTA → (querySeq, oneLineMsa for JSON)
    SeqModification(file, ...)   strip lowercase / gaps → *Mod.fasta
    MmaToA3m(in, out)            MMA/CMA → A3M
    BuildAf3Msa(query, db, ...)  jackhmmer → Stockholm → FASTA → A3M
    RunJackhmmer / StoToFasta / FastaToA3m   lower-level steps

Related modules
    inputAF3.GenerateAf3InputJson → write job JSON with custom MSA
    extractAF3                    → MSA already inside *_data.json

Example
    from machine_learning.alphafold3 import CustomMsa, MmaToA3m, BuildAf3Msa
    from machine_learning.alphafold3 import GenerateAf3InputJson

    MmaToA3m("aln.mma", "aln.a3m")
    # or: BuildAf3Msa("query.fasta", "uniref.fasta", outDir="/tmp/msa")

    seq, oneLineMsa = CustomMsa("aln.a3m")
    GenerateAf3InputJson(
        name="job",
        modelSeeds=[0],
        sequence=seq,
        outputFile="job.json",
        unpairedMsa=oneLineMsa,
        pairedMsa="",
        templates=[],
    )

CLI
    python -m machine_learning.alphafold3.msaAF3 mma in.mma out.a3m
    python -m machine_learning.alphafold3.msaAF3 jackhmmer -query q.fa -db db.fa -out /tmp/msa
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path
from typing import List, Optional, Union

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Pack / clean MSA for AF3 JSON
# ---------------------------------------------------------------------------

def CustomMsa(file: PathLike):
    """
    Read an A3M/FASTA MSA and return ``(querySequence, oneLineMsa)`` for AF3 JSON.

    Lowercase inserts are stripped; query gaps become ``X``.
    """
    lines = Path(file).expanduser().resolve().read_text().splitlines(keepends=True)
    oneLineMsa = ""
    sequence = ""
    seqIndex = 0
    for line in lines:
        if len(line) <= 1:
            continue
        if line[0] == ">":
            oneLineMsa += line
            continue
        noLowerSeq = "".join(c for c in line if not c.islower())
        if seqIndex == 0:
            gapped = noLowerSeq.replace("-", "X")
            oneLineMsa += gapped
            sequence = gapped.rstrip("\n")
            print(
                f"msaAF3: query length {len(sequence)} "
                f"({sequence.count('X')} gaps as X)"
            )
        else:
            oneLineMsa += noLowerSeq
        seqIndex += 1
    return sequence, oneLineMsa


def SeqModification(file: PathLike, modRemove: Optional[List[str]] = None):
    """Strip characters from a FASTA/MSA (default: lowercase inserts) → ``*Mod.fasta``."""
    if modRemove is None:
        modRemove = ["lowercase"]
    path = Path(file).expanduser().resolve()
    lines = path.read_text().splitlines(keepends=True)
    modified: List[str] = []
    for line in lines:
        if len(line) > 1 and not line.startswith(">"):
            for token in modRemove:
                if token == "lowercase":
                    line = "".join(c for c in line if not c.islower() or c in "\n\r")
                else:
                    line = line.replace(token, "")
        modified.append(line)

    stem = path.name
    if stem.endswith(".fasta"):
        outName = stem[:-6] + "Mod.fasta"
    else:
        outName = path.stem + "Mod.fasta"
    outPath = path.with_name(outName)
    outPath.write_text("".join(modified))
    print(f"msaAF3: wrote {outPath}")
    return outPath


# ---------------------------------------------------------------------------
# MMA / CMA → A3M
# ---------------------------------------------------------------------------

def _processMmaSequence(line: str) -> str:
    newLine = []
    insideParentheses = False
    for char in line:
        if char == "(":
            insideParentheses = True
            newLine.append(char)
        elif char == ")":
            insideParentheses = False
            newLine.append(char)
        elif insideParentheses:
            newLine.append(char.lower())
        else:
            newLine.append(char)
    processed = "".join(newLine)
    return re.sub(r"[\}\{\]\[\)\(\*]", "", processed)


def MmaToA3m(inputFile: PathLike, outputFile: PathLike) -> Path:
    """Convert an MMA/CMA file to A3M."""
    inPath = Path(inputFile).expanduser().resolve()
    outPath = Path(outputFile).expanduser().resolve()
    outPath.parent.mkdir(parents=True, exist_ok=True)
    with open(inPath, "r") as infile, open(outPath, "w") as outfile:
        for line in infile:
            if line.startswith("$"):
                continue
            if line.startswith(">"):
                outfile.write(line)
                continue
            if line.startswith("{"):
                outfile.write(_processMmaSequence(line) + "\n")
    print(f"msaAF3: wrote {outPath}")
    return outPath


# ---------------------------------------------------------------------------
# jackhmmer → A3M
# ---------------------------------------------------------------------------

def RunJackhmmer(
    queryFasta: PathLike,
    databaseFasta: PathLike,
    outputSto: PathLike,
    iterations: int = 3,
    jackhmmerPath: str = "jackhmmer",
):
    """Run jackhmmer and write a Stockholm alignment."""
    queryFasta = Path(queryFasta).expanduser().resolve()
    databaseFasta = Path(databaseFasta).expanduser().resolve()
    outputSto = Path(outputSto).expanduser().resolve()
    outputSto.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        jackhmmerPath,
        "-N",
        str(iterations),
        "-A",
        str(outputSto),
        str(queryFasta),
        str(databaseFasta),
    ]
    print(f"msaAF3: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return outputSto


def StoToFasta(stockholmFile: PathLike, fastaOutput: PathLike):
    """Convert Stockholm → FASTA."""
    from Bio import AlignIO

    stockholmFile = Path(stockholmFile).expanduser().resolve()
    fastaOutput = Path(fastaOutput).expanduser().resolve()
    aln = AlignIO.read(str(stockholmFile), "stockholm")
    with open(fastaOutput, "w") as outHandle:
        AlignIO.write(aln, outHandle, "fasta")
    print(f"msaAF3: wrote {fastaOutput}")
    return fastaOutput


def FastaToA3m(inputFasta: PathLike, outputFile: PathLike):
    """Convert FASTA MSA → A3M by removing columns that are gaps in the query."""
    from Bio import SeqIO

    inputFasta = Path(inputFasta).expanduser().resolve()
    outputFile = Path(outputFile).expanduser().resolve()
    records = list(SeqIO.parse(str(inputFasta), "fasta"))
    if not records:
        raise ValueError(f"no sequences in {inputFasta}")
    querySeq = records[0].seq
    with open(outputFile, "w") as out:
        for record in records:
            a3mSeq = "".join(
                res if q != "-" else "" for q, res in zip(querySeq, record.seq)
            )
            out.write(f">{record.id}\n{a3mSeq}\n")
    print(f"msaAF3: wrote {outputFile}")
    return outputFile


def BuildAf3Msa(
    queryFasta: PathLike,
    databaseFasta: PathLike,
    outDir: PathLike = "/tmp/jackhmmer",
    iterations: int = 3,
    jackhmmerPath: str = "jackhmmer",
) -> Path:
    """jackhmmer → Stockholm → FASTA → A3M; returns path to ``.a3m``."""
    outPath = Path(outDir).expanduser().resolve()
    outPath.mkdir(parents=True, exist_ok=True)
    stoFile = outPath / "jackhmmer_output.sto"
    fastaFile = outPath / "msa.fasta"
    a3mFile = outPath / "msa.a3m"
    RunJackhmmer(
        queryFasta,
        databaseFasta,
        stoFile,
        iterations=iterations,
        jackhmmerPath=jackhmmerPath,
    )
    StoToFasta(stoFile, fastaFile)
    return FastaToA3m(fastaFile, a3mFile)


# Aliases
customMSA = CustomMsa
seqModification = SeqModification
remove_mma_patterns = MmaToA3m
run_jackhmmer = RunJackhmmer
sto_to_fasta = StoToFasta
fasta_to_GapRemoved_a3m = FastaToA3m


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="MSA helpers for AF3 (jackhmmer / mma / custom pack)"
    )
    sub = parser.add_subparsers(dest="command")

    pJack = sub.add_parser("jackhmmer", help="Build A3M via jackhmmer")
    pJack.add_argument("-query", required=True)
    pJack.add_argument("-db", required=True)
    pJack.add_argument("-out", default="/tmp/jackhmmer")
    pJack.add_argument("-N", type=int, default=3)
    pJack.add_argument("-jackhmmer", default="jackhmmer")

    pMma = sub.add_parser("mma", help="Convert MMA/CMA → A3M")
    pMma.add_argument("inputFile")
    pMma.add_argument("outputFile")

    args = parser.parse_args(argv)
    if args.command == "jackhmmer":
        BuildAf3Msa(
            args.query,
            args.db,
            outDir=args.out,
            iterations=args.N,
            jackhmmerPath=args.jackhmmer,
        )
    elif args.command == "mma":
        MmaToA3m(args.inputFile, args.outputFile)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
