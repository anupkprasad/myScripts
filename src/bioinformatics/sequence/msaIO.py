#!/usr/bin/env python3
"""
msaIO — read / write / convert sequence alignments (FASTA, A3M, CMA).

When to use
    Load an MSA for cleaning, stats, clustering, or later conservation plots.
    Convert between formats used by AF3 / jackhmmer / CMA tools.

Public API
    LoadMsa(path, fmt="fasta"|"a3m"|"clustal") → list[SeqRecord]
    WriteMsa(records, path, fmt="fasta")
    MultiToSingleLineFasta(inFile, outFile)
    A3mToFasta(inFile, outFile, stripInserts=True)
    CmaToFasta(inFile, outFile)

Related
    msaStats / msaClean / msaView — analysis on loaded records
    machine_learning.alphafold3.extractAF3 — pull A3M from AF3 *_data.json

Example
    from bioinformatics.sequence.msaIO import LoadMsa, A3mToFasta, WriteMsa

    records = LoadMsa("aln.a3m", fmt="a3m")
    WriteMsa(records, "aln.fasta")
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from textwrap import fill
from typing import List, Optional, Sequence, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PathLike = Union[str, Path]


def LoadMsa(path: PathLike, fmt: str = "fasta") -> List[SeqRecord]:
    """
    Load an MSA as a list of SeqRecords (equal length for fasta/clustal).

    fmt:
      - ``fasta`` / ``clustal`` — Biopython readers
      - ``a3m`` — strips lowercase inserts so columns align
    """
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    fmt = fmt.lower()
    if fmt == "a3m":
        return _loadA3m(path)
    records = list(SeqIO.parse(str(path), fmt))
    if not records:
        raise ValueError(f"no sequences in {path}")
    return records


def WriteMsa(
    records: Sequence[SeqRecord],
    path: PathLike,
    fmt: str = "fasta",
    singleLine: bool = True,
) -> Path:
    """
    Write SeqRecords to FASTA (or other Biopython-supported format).

    For ``fmt="fasta"``, ``singleLine=True`` (default) writes one sequence
    line per record (no 60-char wrapping).
    """
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    if fmt == "fasta" and singleLine:
        with open(path, "w") as fh:
            for rec in records:
                header = rec.id
                desc = (rec.description or "").strip()
                if (
                    desc
                    and desc != rec.id
                    and desc != "<unknown description>"
                    and not desc.startswith(rec.id)
                ):
                    header = f"{rec.id} {desc}"
                elif desc.startswith(rec.id) and desc != rec.id:
                    header = desc
                fh.write(f">{header}\n{str(rec.seq)}\n")
    else:
        SeqIO.write(list(records), str(path), fmt)
    print(f"msaIO: wrote {len(records)} sequences → {path}")
    return path


def _loadA3m(path: Path) -> List[SeqRecord]:
    records: List[SeqRecord] = []
    header: Optional[str] = None
    chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(c for c in "".join(chunks) if not c.islower())
                    records.append(SeqRecord(Seq(seq), id=header, description=""))
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seq = "".join(c for c in "".join(chunks) if not c.islower())
            records.append(SeqRecord(Seq(seq), id=header, description=""))
    if not records:
        raise ValueError(f"no sequences in {path}")
    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1:
        raise ValueError(
            f"A3M columns unequal after stripping inserts: lengths={sorted(lengths)[:5]}…"
        )
    return records


def MultiToSingleLineFasta(inputFile: PathLike, outputFile: PathLike) -> Path:
    """Convert wrapped FASTA → one sequence line per record."""
    inputFile = Path(inputFile).expanduser().resolve()
    outputFile = Path(outputFile).expanduser().resolve()
    data: List[str] = []
    firstHeader = True
    with open(inputFile) as infile:
        for line in infile:
            if not line.startswith(">"):
                cleaned = line.strip().replace(" ", "")
                if cleaned:
                    data.append(cleaned)
            else:
                if firstHeader:
                    data.append(line.strip() + "\n")
                    firstHeader = False
                else:
                    data.append("\n" + line.strip() + "\n")
    outputFile.parent.mkdir(parents=True, exist_ok=True)
    with open(outputFile, "w") as outfile:
        for line in data:
            outfile.write(line)
        if data and not str(data[-1]).endswith("\n"):
            outfile.write("\n")
    print(f"msaIO: single-line FASTA → {outputFile}")
    return outputFile


def A3mToFasta(
    inputFile: PathLike,
    outputFile: PathLike,
    stripInserts: bool = True,
) -> Path:
    """A3M → FASTA; by default drop lowercase inserts."""
    inputFile = Path(inputFile).expanduser().resolve()
    outputFile = Path(outputFile).expanduser().resolve()
    outputFile.parent.mkdir(parents=True, exist_ok=True)
    with open(inputFile) as infile, open(outputFile, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                outfile.write(line)
            elif stripInserts:
                outfile.write("".join(c for c in line if not c.islower()))
            else:
                outfile.write(line)
    print(f"msaIO: A3M→FASTA → {outputFile}")
    return outputFile


def CmaToFasta(
    inputFile: PathLike,
    outputFile: PathLike,
    width: int = 60,
    headerTokenFirst: bool = False,
) -> Path:
    """Extract letter sequences from CMA blocks → FASTA."""
    inputFile = Path(inputFile).expanduser().resolve()
    outputFile = Path(outputFile).expanduser().resolve()
    text = inputFile.read_text(encoding="utf-8", errors="ignore")
    seqs = []
    for match in re.finditer(r">([^\n\r]+)", text):
        header = match.group(1).strip()
        startPos = match.end()
        openBrace = text.find("{", startPos)
        if openBrace == -1:
            seqs.append((header, ""))
            continue
        closeMarker = text.find("()}", openBrace)
        if closeMarker != -1:
            seqBlock = text[openBrace + 1 : closeMarker]
        else:
            closeBrace = text.find("}", openBrace)
            seqBlock = text[openBrace + 1 :] if closeBrace == -1 else text[openBrace + 1 : closeBrace]
        seqLetters = re.sub(r"[^A-Za-z]", "", seqBlock).upper()
        seqs.append((header, seqLetters))

    outputFile.parent.mkdir(parents=True, exist_ok=True)
    nEmpty = 0
    with open(outputFile, "w") as out:
        for header, seq in seqs:
            fastaHeader = header.split()[0] if headerTokenFirst else header.replace(">", "").strip()
            out.write(f">{fastaHeader}\n")
            if seq:
                out.write(fill(seq, width=width) + "\n")
            else:
                nEmpty += 1
                out.write("\n")
    if nEmpty:
        print(f"msaIO: warning — {nEmpty} empty sequence(s)")
    print(f"msaIO: CMA→FASTA ({len(seqs)} entries) → {outputFile}")
    return outputFile


# Aliases
multi2SingleLineFasta = MultiToSingleLineFasta
a3mtoFasta = A3mToFasta
cmatoFasta = CmaToFasta


def main(argv=None):
    parser = argparse.ArgumentParser(description="MSA I/O and format conversion")
    parser.add_argument("inputFile")
    parser.add_argument("outputFile")
    parser.add_argument(
        "-task",
        default="a3m2fasta",
        choices=["a3m2fasta", "multi2single", "cma2fasta"],
    )
    args = parser.parse_args(argv)
    if args.task == "cma2fasta":
        CmaToFasta(args.inputFile, args.outputFile)
    elif args.task == "multi2single":
        MultiToSingleLineFasta(args.inputFile, args.outputFile)
    else:
        A3mToFasta(args.inputFile, args.outputFile)


if __name__ == "__main__":
    main()
