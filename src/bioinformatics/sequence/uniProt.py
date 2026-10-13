#!/usr/bin/env python3
"""
uniProt — fetch sequences / domain annotations from the UniProt REST API.

When to use
    Look up accessions, pull full sequences or annotated domains (e.g. kinase),
    write FASTA for MSA / AF3 input prep.

Public API
    SearchUniprot(query, size=500) → list[dict] or DataFrame-like rows
    GetUniprotIds(query) → list[str]
    GetUniprotSequence(uniprotId, domain="kinase") → dict
    SaveFasta(entries, filename)
    SaveFastaDomains(entries, filePrefix)

Example
    from bioinformatics.sequence import GetUniprotSequence, SaveFasta

    entry = GetUniprotSequence("P00533", domain="kinase")
    SaveFasta([entry], "egfr.fasta")

CLI
    python -m bioinformatics.sequence.uniProt P00533
    python -m bioinformatics.sequence.uniProt -search 'gene:DCLK3 AND organism_id:9606'
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Sequence, Tuple, Union

PathLike = Union[str, Path]


def _requireRequests():
    try:
        import requests
    except ImportError as exc:
        raise SystemExit(
            "uniProt needs requests:\n  conda activate ComBio && pip install requests"
        ) from exc
    return requests


def SearchUniprot(query: str, size: int = 500) -> List[dict]:
    """
    UniProt text search → list of {accession, gene, protein, organism}.

    Example query: ``keyword:\"Pseudokinase [KW-1028]\" AND organism_id:9606``
    """
    requests = _requireRequests()
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "fields": "accession,gene_names,protein_name,organism_name",
        "format": "json",
        "size": size,
    }
    response = requests.get(url, params=params, timeout=60)
    response.raise_for_status()
    data = response.json()
    results = []
    for entry in data.get("results", []):
        gene = ""
        genes = entry.get("genes") or []
        if genes:
            gene = (genes[0].get("geneName") or {}).get("value", "")
        name = ""
        try:
            name = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
        except (KeyError, TypeError):
            pass
        results.append(
            {
                "accession": entry.get("primaryAccession"),
                "gene": gene,
                "protein": name,
                "organism": (entry.get("organism") or {}).get("scientificName", ""),
            }
        )
    print(f"uniProt: {len(results)} hits")
    return results


def GetUniprotIds(query: str) -> List[str]:
    """Return primary accessions for a UniProt search query."""
    return [r["accession"] for r in SearchUniprot(query) if r.get("accession")]


def GetUniprotSequence(uniprotId: str, domain: str = "kinase") -> dict:
    """
    Fetch one entry: full sequence + domains whose description contains ``domain``.

    Returns dict with keys: uniprotId, fullSequence, domains [(start,end,seq)], header
    """
    requests = _requireRequests()
    url = f"https://rest.uniprot.org/uniprotkb/{uniprotId}.json"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    record = response.json()
    sequence = (record.get("sequence") or {}).get("value", "")
    domains = []
    for feature in record.get("features") or []:
        if feature.get("type") == "Domain" and domain.lower() in (
            feature.get("description") or ""
        ).lower():
            start = int(feature["location"]["start"]["value"]) - 1
            end = int(feature["location"]["end"]["value"])
            domains.append((start + 1, end, sequence[start:end]))

    organism = (record.get("organism") or {}).get("scientificName", "")
    header = f"{record.get('primaryAccession')}|{record.get('uniProtkbId')}|OS={organism}"
    return {
        "uniprotId": uniprotId,
        "fullSequence": sequence,
        "domains": domains,
        "header": header,
        # backward-compatible keys
        "uniprot_id": uniprotId,
        "full_sequence": sequence,
        "kinase_domains": domains,
    }


def SaveFasta(entries: Sequence[dict], filename: PathLike = "sequences.fasta") -> Path:
    """Write full sequences + numbered domains to one FASTA."""
    path = Path(filename).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for entry in entries:
            uid = entry.get("uniprotId") or entry.get("uniprot_id")
            full = entry.get("fullSequence") or entry.get("full_sequence") or ""
            fh.write(f">{uid}_full\n{full}\n")
            domains = entry.get("domains") or entry.get("kinase_domains") or []
            for i, (start, end, seq) in enumerate(domains, 1):
                fh.write(f">{uid}_domain{i}_{start}-{end}\n{seq}\n")
    print(f"uniProt: wrote {path}")
    return path


def SaveFastaDomains(
    entries: Sequence[dict],
    filePrefix: PathLike = "uniprotseqs",
) -> Tuple[Path, Path]:
    """Write separate full-length and domain FASTA files."""
    prefix = Path(filePrefix)
    fullPath = prefix.parent / f"{prefix.name}_full_length.fasta"
    domPath = prefix.parent / f"{prefix.name}_kinase_domains.fasta"
    fullPath.parent.mkdir(parents=True, exist_ok=True)
    with open(fullPath, "w") as fh:
        for entry in entries:
            fh.write(
                f">FullLength|{entry.get('header', entry.get('uniprotId'))}\n"
                f"{entry.get('fullSequence') or entry.get('full_sequence')}\n"
            )
    with open(domPath, "w") as fh:
        for entry in entries:
            domains = entry.get("domains") or entry.get("kinase_domains") or []
            for i, (start, end, seq) in enumerate(domains, 1):
                fh.write(
                    f">Dom{i}_{start}-{end}|{entry.get('header', '')}\n{seq}\n"
                )
    print(f"uniProt: wrote {fullPath} and {domPath}")
    return fullPath, domPath


# Aliases
getUniprotIds = GetUniprotIds
getUniprotSequence = GetUniprotSequence
saveFasta = SaveFasta


def main(argv=None):
    parser = argparse.ArgumentParser(description="UniProt sequence / domain fetch")
    parser.add_argument("uniprotId", nargs="?", default=None)
    parser.add_argument("-search", default=None, help="UniProt search query")
    parser.add_argument("-domain", default="kinase")
    parser.add_argument("-out", default=None)
    args = parser.parse_args(argv)

    if args.search:
        hits = SearchUniprot(args.search)
        for h in hits[:20]:
            print(h["accession"], h.get("gene"), h.get("protein"))
        if len(hits) > 20:
            print(f"... {len(hits)} total")
        return 0

    if not args.uniprotId:
        parser.print_help()
        return 1

    result = GetUniprotSequence(args.uniprotId, domain=args.domain)
    print(f"uniProt: {args.uniprotId}")
    print(f"  length: {len(result['fullSequence'])}")
    print(f"  domains ({args.domain}): {len(result['domains'])}")
    for start, end, seq in result["domains"]:
        print(f"    {start}-{end} ({len(seq)} aa)")
    print(f"  header: {result['header']}")
    if args.out:
        SaveFasta([result], args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
