#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 23:38:12 2024

@author: anupkumar
"""
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import argparse

def hamming_distance(seq1, seq2):
    """Count mismatches between two sequences (ignoring gaps)."""
    return sum((a != b) and (a != '-' and b != '-') for a, b in zip(seq1, seq2))

def compute_conservation(alignment):
    """Calculate conservation score for each column (as frequency of most common residue)."""
    from collections import Counter
    scores = []
    for i in range(alignment.get_alignment_length()):
        column = [record.seq[i] for record in alignment]
        column = [res for res in column if res != '-']
        if not column:
            scores.append(0.0)
        else:
            common = Counter(column).most_common(1)[0][1]
            scores.append(common / len(column))
    return scores

def sort_by_similarity(alignment, query_id):
    """Reorder alignment with query on top, and others sorted by similarity to query."""
    query_seq = next(record for record in alignment if query_id in record.id)
    others = [rec for rec in alignment if query_id not in rec.id]
    others.sort(key=lambda rec: hamming_distance(query_seq.seq, rec.seq))
    return MultipleSeqAlignment([query_seq] + others)

def highlight_differences(query_seq, aligned_seq):
    """Highlight mismatches with lowercase; matches kept as is."""
    return ''.join(a if a == b else b.lower() for a, b in zip(query_seq, aligned_seq))

def display_alignment(alignment, query_id):
    query_record = next(record for record in alignment if query_id in record.id)
    print(f"\n> {query_record.id}\n{query_record.seq}")
    for rec in alignment:
        if rec.id != query_record.id:
            highlighted = highlight_differences(query_record.seq, rec.seq)
            print(f"\n> {rec.id}\n{highlighted}")

def main(aln_file, query_id, format="fasta"):
    alignment = AlignIO.read(aln_file, format)
    sorted_alignment = sort_by_similarity(alignment, query_id)
    display_alignment(sorted_alignment, query_id)
    
    print("\nPer-position conservation (0–1):")
    conservation = compute_conservation(sorted_alignment)
    print("".join(f"{int(c*10)}" for c in conservation))  # crude visual, 0–10 scale


msa_file = "/mnt/mydrive/alpfaFold3/dynamics/DCLK_with_tail/DCLK1/aligned.fasta"
query_id = "NP_001317000"  # Replace with your actual query sequence ID


main(msa_file, query_id)
dsdf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Query-centric MSA viewer with conservation and difference highlighting.")
    parser.add_argument("aln_file", help="Path to alignment file (FASTA or Clustal format).")
    parser.add_argument("query_id", help="ID of query sequence to fix on top.")
    parser.add_argument("--format", default="fasta", help="Alignment file format (default: fasta)")
    args = parser.parse_args()
    
    main(args.aln_file, args.query_id, args.format)

