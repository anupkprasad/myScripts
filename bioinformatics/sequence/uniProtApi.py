#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-08-05 (Y/M/D) at 13:16
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

import requests
import time
import argparse

def getUniprotIds(query:str=None):
    """Only search by text string, UniProt REST API supports simple keyword search
    query = 'catalytically inactive AND organism_id:9606 AND keyword:Kinase'
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "fields": "accession",
        "format": "json",
        "size": 500  # max UniProt allows
    }
    response = requests.get(url, params=params)
    response.raise_for_status()
    data = response.json()
    ids = [entry["primaryAccession"] for entry in data.get("results", [])]
    return ids


def getUniprotSequence(uniprot_id, domain="kinase"):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    response.raise_for_status()
    record = response.json()

    sequence = record.get("sequence", {}).get("value", "")
    domainResult = []
    for feature in record.get("features", []):
        if feature["type"] == "Domain" and domain in feature.get("description", "").lower():
            start = int(feature["location"]["start"]["value"]) - 1
            end = int(feature["location"]["end"]["value"])
            domain_seq = sequence[start:end]
            domainResult.append((start + 1, end, domain_seq)) ## has three elements: start, end, sequence

  
    data = record
    info = {
        "Accession": data.get("primaryAccession"),
        "Entry Name": data.get("uniProtkbId"),
        #"Protein Name": data["proteinDescription"]["recommendedName"]["fullName"]["value"],
        "Organism": data["organism"]["scientificName"],
        # "Taxonomy ID": data["organism"]["taxonId"],
        #"Gene Name": data["genes"][0]["geneName"]["value"] if data.get("genes") else None,
        #"Protein Existence (PE)": data["proteinExistence"]["evidenceCode"],
        # "Sequence Version (SV)": data["sequence"]["entryVersion"],
        "Sequence": data["sequence"]["value"]
    }
    header = f">{info['Accession']}|{info['Entry Name']}|OS={info['Organism']}"
    return {
        "uniprot_id": uniprot_id,
        "full_sequence": sequence,
        "kinase_domains": domainResult, # List of tuples (start, end, sequence)
        "header": header
    }

def saveFasta(data, filename="pseudokinases.fasta"):
    with open(filename, "w") as f:
        for entry in data:
            f.write(f">{entry['uniprot_id']}_full\n{entry['full_sequence']}\n")
            for i, (start, end, seq) in enumerate(entry["kinase_domains"], 1):
                f.write(f">{entry['uniprot_id']}_domain{i}_{start}-{end}\n{seq}\n")


def saveFastaKDnfullLength(data, filePrefix="uniprotseqs"):
    ### data is a list of dictionaries with keys: uniprot_id, full_sequence, kinase_domains, header
    with open(filePrefix + "_full_length.fasta", "w") as f:
        for entry in data:
            f.write(f">FullLength|{entry['header']}\n{entry['full_sequence']}\n")

    with open(filePrefix + "_kinase_domains.fasta", "w") as f:
        for entry in data:
            for i, (start, end, seq) in enumerate(entry["kinase_domains"], 1):
                f.write(f">KinDom{i}_{start}-{end}|{entry['header']}\n{seq}\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch human pseudokinase sequences from UniProt")
    parser.add_argument("uniprot_ids", help="Enter a uniprot ID" )
    args = parser.parse_args()
    result = getUniprotSequence(args.uniprot_ids)
    print("\n")
    print(f"Fetched sequence for {args.uniprot_ids}")
    print(f"Full sequence length: {result['full_sequence']}")
    print("\n")
    print(f"Kinase domains found: {result['kinase_domains']}")
    print(f"Header: {result['header']}")

