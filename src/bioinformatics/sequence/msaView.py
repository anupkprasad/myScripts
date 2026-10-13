#!/usr/bin/env python3
"""
msaView — inspect an MSA: sort by query similarity, highlight diffs, Dash viewer.

When to use
    Quick visual check of an alignment before / after cleaning. For publication
    conservation strips, prefer msaStats.PlotConservation (and later richer plots).

Public API
    SortBySimilarity(alignment|records, queryId)
    HighlightDifferences(querySeq, alignedSeq)
    DisplayAlignment(records, queryId)          print to terminal
    ShowConservationDigits(records)             crude 0–9 conservation bar
    ViewAlignment(fastaFile, port=8050)         interactive Dash Bio chart

Example
    from bioinformatics.sequence import LoadMsa, DisplayAlignment

    records = LoadMsa("aln.fasta")
    DisplayAlignment(records)                 # first sequence as query
    DisplayAlignment(records, queryId="query")  # AF3-style query header

CLI
    python -m bioinformatics.sequence.msaView aln.fasta
    python -m bioinformatics.sequence.msaView aln.fasta query
    python -m bioinformatics.sequence.msaView aln.fasta -dash
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Sequence, Union

from Bio.SeqRecord import SeqRecord

from .msaIO import LoadMsa
from .msaStats import ColumnConservation

PathLike = Union[str, Path]


def HammingDistance(seq1, seq2) -> int:
    """Count mismatches ignoring gap–gap and gap–residue pairs."""
    return sum(
        (a != b) and (a != "-") and (b != "-") for a, b in zip(str(seq1), str(seq2))
    )


def _resolveQuery(
    records: Sequence[SeqRecord],
    queryId: Optional[str] = None,
) -> SeqRecord:
    """Pick query record; default = first sequence. ``QUERY_ID`` is treated as unset."""
    if not records:
        raise ValueError("empty alignment")
    if queryId is None or queryId in {"", "QUERY_ID", "queryId"}:
        return records[0]
    query = next((r for r in records if queryId in r.id), None)
    if query is not None:
        return query
    sample = ", ".join(repr(r.id) for r in records[:8])
    more = "" if len(records) <= 8 else f", … ({len(records)} total)"
    raise KeyError(
        f"queryId {queryId!r} not found in alignment.\n"
        f"  First sequence id: {records[0].id!r}\n"
        f"  Examples: {sample}{more}\n"
        f"  Tip: use queryId={records[0].id!r} or omit queryId to use the first sequence."
    )


def SortBySimilarity(
    records: Sequence[SeqRecord],
    queryId: Optional[str] = None,
) -> List[SeqRecord]:
    """Reorder with query on top; others by increasing Hamming distance to query."""
    query = _resolveQuery(records, queryId)
    others = [r for r in records if r.id != query.id]
    others.sort(key=lambda r: HammingDistance(query.seq, r.seq))
    return [query] + others


def HighlightDifferences(querySeq, alignedSeq) -> str:
    """Matches uppercase; mismatches lowercase (relative to query)."""
    q = str(querySeq)
    a = str(alignedSeq)
    return "".join(b if a_i == b else b.lower() for a_i, b in zip(q, a))


def DisplayAlignment(
    records: Sequence[SeqRecord],
    queryId: Optional[str] = None,
    maxSeqs: Optional[int] = 20,
) -> None:
    """
    Print query, then others with mismatches lowercased.

    queryId: substring match on record.id; default / ``QUERY_ID`` → first sequence.
    maxSeqs: limit how many sequences to print (None = all). Large MSAs are noisy.
    """
    ordered = SortBySimilarity(records, queryId)
    if maxSeqs is not None:
        ordered = ordered[: max(1, maxSeqs)]
    query = ordered[0]
    print(f"msaView: query = {query.id!r} ({len(ordered)} sequences shown)")
    print(f"\n> {query.id}\n{query.seq}")
    for rec in ordered[1:]:
        print(f"\n> {rec.id}\n{HighlightDifferences(query.seq, rec.seq)}")


def ShowConservationDigits(records: Sequence[SeqRecord]) -> str:
    """Crude 0–9 conservation string (consensus frequency × 10)."""
    scores = ColumnConservation(records)
    text = "".join(str(min(9, int(c * 10))) for c in scores)
    print("\nPer-position conservation (0–9):")
    print(text)
    return text


def ViewAlignment(fastaFile: PathLike, port: int = 8050):
    """Launch interactive MSA viewer (requires dash + dash-bio)."""
    try:
        from dash import Dash, html
        import dash_bio as dashbio
    except ImportError as exc:
        raise SystemExit(
            "ViewAlignment needs dash and dash-bio:\n"
            "  conda activate ComBio && pip install dash dash-bio"
        ) from exc

    path = Path(fastaFile).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    data = path.read_text()
    app = Dash(__name__)
    app.layout = html.Div(
        [
            html.H3(f"MSA Viewer — {path.name}"),
            dashbio.AlignmentChart(id="alignment-viewer", data=data),
        ]
    )
    print(f"msaView: http://127.0.0.1:{port}/")
    app.run(port=port, debug=True)


# Aliases
viewAlignment = ViewAlignment
sort_by_similarity = SortBySimilarity


def main(argv=None):
    parser = argparse.ArgumentParser(description="Query-centric MSA viewer")
    parser.add_argument("msaFile")
    parser.add_argument(
        "queryId",
        nargs="?",
        default=None,
        help="substring of query id (default: first sequence)",
    )
    parser.add_argument("-fmt", default="fasta", choices=["fasta", "a3m", "clustal"])
    parser.add_argument("-dash", action="store_true", help="open Dash Bio viewer")
    parser.add_argument("-port", type=int, default=8050)
    parser.add_argument("-maxSeqs", type=int, default=20)
    args = parser.parse_args(argv)
    if args.dash:
        ViewAlignment(args.msaFile, port=args.port)
        return
    records = LoadMsa(args.msaFile, fmt=args.fmt)
    DisplayAlignment(records, args.queryId, maxSeqs=args.maxSeqs)
    ShowConservationDigits(records)


if __name__ == "__main__":
    main()
