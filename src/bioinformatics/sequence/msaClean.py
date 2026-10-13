#!/usr/bin/env python3
"""
msaClean — clean an MSA for coevolution / DCA-style analysis + QC plots.

When to use
    Raw jackhmmer / AF3 MSAs are gappy and redundant. CleanMsa filters by
    coverage, optional paralogs, trims gappy columns, computes Neff, and
    writes QC plots via msaStats.

Public API
    CleanMsa(inputFile, outputFile, outDir=..., minCov=0.7, maxColGap=0.3, ...)
    FilterByCoverage(records, minCov)
    RemoveParalogs(records)
    TrimGappyColumns(records, maxColGap)

Defaults (common DCA practice)
    minCov=0.70, maxColGap=0.30, theta=0.80

Example
    from bioinformatics.sequence import CleanMsa
    CleanMsa("raw.fasta", "clean.fasta", outDir="msa_qc")

CLI
    python -m bioinformatics.sequence.msaClean raw.fasta clean.fasta -outdir msa_qc
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
from Bio.Seq import Seq

from .msaIO import LoadMsa, WriteMsa
from .msaStats import (
    EffectiveSeqCount,
    GapFraction,
    PlotColumnProfile,
    PlotEntropy,
    PlotSeqIdMatrix,
    SeqIdentityMatrix,
)

PathLike = Union[str, Path]


def GetSpeciesId(header: str) -> str:
    """Heuristic species key from a FASTA header (OS=… or trailing taxon token)."""
    h = header.strip()
    match = re.search(r"OS=([A-Za-z_ ]+)", h)
    if match:
        return match.group(1).strip().replace(" ", "_")
    match = re.search(r"([A-Z]{3,6})$", h)
    if match:
        return match.group(1)
    return h.split()[0]


def CoverageFraction(seq: str) -> float:
    if not seq:
        return 0.0
    return 1.0 - (seq.count("-") / len(seq))


def FilterByCoverage(records: Sequence, minCov: float = 0.7) -> List:
    return [r for r in records if CoverageFraction(str(r.seq)) >= minCov]


def RemoveParalogs(records: Sequence) -> List:
    """Keep one representative per species (best coverage)."""
    speciesMap = defaultdict(list)
    for rec in records:
        speciesMap[GetSpeciesId(rec.id)].append(rec)
    reps = []
    for recs in speciesMap.values():
        rep = max(recs, key=lambda x: len(str(x.seq)) - str(x.seq).count("-"))
        reps.append(rep)
    return reps


def TrimGappyColumns(
    records: Sequence,
    maxColGap: float = 0.3,
) -> Tuple[List, np.ndarray]:
    """Drop columns with gap fraction > maxColGap. Returns (newRecords, gapFracBefore)."""
    arr = np.asarray([list(str(r.seq)) for r in records])
    gapFrac = (arr == "-").mean(axis=0)
    keepMask = gapFrac <= maxColGap
    if keepMask.sum() == 0:
        raise ValueError("no columns remain after trimming; relax maxColGap")
    from Bio.SeqRecord import SeqRecord

    out = []
    for i, rec in enumerate(records):
        out.append(
            SeqRecord(
                Seq("".join(arr[i, keepMask])),
                id=rec.id,
                name=getattr(rec, "name", rec.id),
                description=getattr(rec, "description", ""),
            )
        )
    return out, gapFrac


def CleanMsa(
    inputFile: PathLike,
    outputFile: PathLike,
    outDir: PathLike = "msa_cleaning_results",
    minCov: float = 0.70,
    maxColGap: float = 0.30,
    theta: float = 0.80,
    removeParalogs: bool = False,
    fmt: str = "fasta",
) -> dict:
    """
    Full cleaning pipeline. Writes cleaned MSA, JSON report, summary, QC plots.

    Returns the report dict.
    """
    outPath = Path(outDir).expanduser().resolve()
    outPath.mkdir(parents=True, exist_ok=True)
    outputFile = Path(outputFile).expanduser().resolve()

    report = {
        "inputMsa": str(Path(inputFile).expanduser().resolve()),
        "cleanedMsa": str(outputFile),
        "parameters": {
            "minCov": minCov,
            "maxColGap": maxColGap,
            "theta": theta,
            "removeParalogs": removeParalogs,
            "fmt": fmt,
        },
        "startTime": datetime.now(timezone.utc).isoformat(),
    }

    records = LoadMsa(inputFile, fmt=fmt)
    report["nInitialSequences"] = len(records)
    report["initialAlignmentLength"] = len(records[0].seq)

    records = FilterByCoverage(records, minCov=minCov)
    report["nAfterCoverage"] = len(records)
    if not records:
        raise ValueError("no sequences survive coverage filter; relax minCov")

    if removeParalogs:
        records = RemoveParalogs(records)
    report["nAfterParalogFilter"] = len(records)

    origGapFrac = GapFraction(records)
    records, _ = TrimGappyColumns(records, maxColGap=maxColGap)

    weights, neff = EffectiveSeqCount(records, theta=theta)
    length = len(records[0].seq)
    report["nFinalSequences"] = len(records)
    report["finalAlignmentLength"] = length
    report["Neff"] = neff
    report["NeffOverL"] = float(neff / length) if length else 0.0

    # QC plots
    PlotColumnProfile(
        origGapFrac,
        outPath / "gap_profile_original.png",
        ylabel="gap fraction",
        title="Gap fraction per column (before trim)",
        ylim=(0, 1),
    )
    PlotEntropy(records, outPath / "entropy_profile.png")
    PlotSeqIdMatrix(records, outPath / "seqid_matrix.png")
    mat = SeqIdentityMatrix(records)
    triu = mat[np.triu_indices_from(mat, k=1)]
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.hist(triu, bins=40)
    ax.set_xlabel("pairwise sequence identity")
    ax.set_ylabel("count")
    ax.set_title("Pairwise sequence identity distribution")
    fig.tight_layout()
    fig.savefig(outPath / "seqid_hist.png", dpi=150)
    plt.close(fig)
    print(f"msaClean: wrote {outPath / 'seqid_hist.png'}")

    WriteMsa(records, outputFile, singleLine=True)

    report["endTime"] = datetime.now(timezone.utc).isoformat()
    report["plots"] = {
        "seqidMatrix": "seqid_matrix.png",
        "seqidHistogram": "seqid_hist.png",
        "gapProfileOriginal": "gap_profile_original.png",
        "entropyProfile": "entropy_profile.png",
    }
    reportPath = outPath / "cleaning_report.json"
    reportPath.write_text(json.dumps(report, indent=2))
    summaryPath = outPath / "summary.txt"
    summaryPath.write_text(
        "\n".join(
            [
                "MSA cleaning report",
                f"Input MSA: {report['inputMsa']}",
                f"Cleaned MSA: {report['cleanedMsa']}",
                f"Initial sequences: {report['nInitialSequences']}",
                f"After coverage filter: {report['nAfterCoverage']}",
                f"After paralog filter: {report['nAfterParalogFilter']}",
                f"Final sequences: {report['nFinalSequences']}",
                f"Initial alignment length: {report['initialAlignmentLength']}",
                f"Final alignment length: {report['finalAlignmentLength']}",
                f"Neff: {report['Neff']:.2f}",
                f"Neff/L: {report['NeffOverL']:.3f}",
                f"Parameters: minCov={minCov}, maxColGap={maxColGap}, theta={theta}",
                f"Plots directory: {outPath}",
                "",
            ]
        )
    )
    print(f"msaClean: done → {outPath}")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description="Clean MSA FASTA/A3M + QC plots")
    parser.add_argument("inputFile")
    parser.add_argument("outputFile")
    parser.add_argument("-outdir", default="msa_cleaning_results")
    parser.add_argument("-fmt", default="fasta", choices=["fasta", "a3m", "clustal"])
    parser.add_argument("-removeParalogs", action="store_true")
    parser.add_argument("-minCov", type=float, default=0.70)
    parser.add_argument("-maxColGap", type=float, default=0.30)
    parser.add_argument("-theta", type=float, default=0.80)
    args = parser.parse_args(argv)
    CleanMsa(
        args.inputFile,
        args.outputFile,
        outDir=args.outdir,
        minCov=args.minCov,
        maxColGap=args.maxColGap,
        theta=args.theta,
        removeParalogs=args.removeParalogs,
        fmt=args.fmt,
    )


if __name__ == "__main__":
    main()
