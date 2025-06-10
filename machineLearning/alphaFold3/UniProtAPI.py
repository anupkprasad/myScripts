import requests
from Bio import SeqIO
from io import StringIO

# Step 1: Query UniProt for human kinase sequences (limited example)
query_url = "https://rest.uniprot.org/uniprotkb/stream"
params = {
    "query": 'family:"kinase family" AND organism_id:9606',
    "format": "fasta",
    "size": 100  # adjust as needed
}
response = requests.get(query_url, params=params)
fasta_data = response.text

# Save to a .fasta file
with open("kinases_human.fasta", "w") as f:
    f.write(fasta_data)

# Step 2: Parse sequences and extract C-terminal tails
sequences = list(SeqIO.parse(StringIO(fasta_data), "fasta"))
for record in sequences:
    tail = record.seq[-50:]  # last 50 residues as the C-terminal tail
    print(f">{record.id}_tail\n{tail}")

