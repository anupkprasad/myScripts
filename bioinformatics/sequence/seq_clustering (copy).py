#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-12-05 (Y/M/D) at 14:00
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

#!/usr/bin/env python3

"""
Extract protein family assignment and confidence metrics from HMMER domtblout.
"""

import csv

def parse_domtblout(file):
    """
    Parse domtblout file.
    Returns a list of dicts.
    """

    results = []

    with open(file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split()

            # Handle minimal valid length
            if len(fields) < 22:
                continue

            entry = {
                "sequence_id": fields[0],
                "pfam_id": fields[1],
                "i_evalue": float(fields[12]),
                "domain_score": float(fields[13]),
                "hmm_from": int(fields[15]),
                "hmm_to": int(fields[16]),
                "ali_from": int(fields[17]),
                "ali_to": int(fields[18]),
                "family_description": " ".join(fields[22:])
            }

            results.append(entry)

    return results


def write_output(rows, outfile):
    """
    Write results table for visualization.
    """

    # define output columns
    columns = [
        "sequence_id",
        "pfam_id",
        "family_description",
        "i_evalue",
        "domain_score",
        "hmm_from",
        "hmm_to",
        "ali_from",
        "ali_to"
    ]

    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved: {outfile}")


















if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract domain hits from domtblout.")
    parser.add_argument("-i", "--input", required=True, help="domtblout file")
    parser.add_argument("-o", "--output", required=True, help="output CSV")

    args = parser.parse_args()

    hits = parse_domtblout(args.input)
    write_output(hits, args.output)
