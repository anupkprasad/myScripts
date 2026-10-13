# sequence/

MSA and sequence utilities: load/convert alignments, clean MSAs, column stats
(conservation / entropy), view alignments, cluster sequences, fetch UniProt.

Naming: **lowerCamelCase** files, **PascalCase** public APIs. Each `.py` file has a
top docstring with purpose and examples.

| File | Role | Main API |
|------|------|----------|
| `msaIO.py` | Read/write/convert FASTA, A3M, CMA | `LoadMsa`, `WriteMsa`, `A3mToFasta`, `CmaToFasta` |
| `msaStats.py` | Column stats + profile plots | `ColumnConservation`, `PlotConservation`, `PlotEntropy` |
| `msaPlot.py` | Colored MSA figure + conservation strip | `PlotMsa` |
| `phyloTree.py` | Distance tree (NJ / UPGMA) | `PlotPhyloTree`, `BuildDistanceTree` |
| `msaClean.py` | Filter/trim MSA + QC report | `CleanMsa` |
| `msaView.py` | Terminal print / Dash viewer | `DisplayAlignment`, `ViewAlignment` |
| `seqCluster.py` | HMMER domtblout → k-mer → UMAP | `ClusterByDomtblout` |
| `uniProt.py` | UniProt REST fetch | `GetUniprotSequence`, `SearchUniprot` |

Import from the package (after `pip install -e .` or `PYTHONPATH=.../src`):

```python
from bioinformatics.sequence import LoadMsa, PlotMsa, PlotPhyloTree, CleanMsa
```

---

## Typical workflows

### 1) AF3 MSA → conservation plot

```python
from bioinformatics.sequence import LoadMsa, PlotConservation

# After: ExtractMsaTemplate(job, outDir="af3_features", chainId="A")
records = LoadMsa("af3_features/chainA_unpaired.a3m", fmt="a3m")
PlotConservation(records[:200], outFile="conservation.png")  # subset if huge
```

### 2) Clean a raw MSA

```python
from bioinformatics.sequence import CleanMsa

CleanMsa("raw.fasta", "clean.fasta", outDir="msa_qc", minCov=0.7, maxColGap=0.3)
```

```bash
python -m bioinformatics.sequence.msaClean raw.fasta clean.fasta -outdir msa_qc -fmt fasta
```

### 3) Format conversion

```python
from bioinformatics.sequence import A3mToFasta, CmaToFasta, MultiToSingleLineFasta

A3mToFasta("aln.a3m", "aln.fasta")
CmaToFasta("aln.cma", "aln.fasta")
MultiToSingleLineFasta("wrapped.fasta", "single.fasta")
```

```bash
python -m bioinformatics.sequence.msaIO aln.a3m aln.fasta -task a3m2fasta
```

### 4) Stats / profiles (CLI)

```bash
python -m bioinformatics.sequence.msaStats aln.fasta -plot all -out stats_out
python -m bioinformatics.sequence.msaStats aln.a3m -fmt a3m -plot conservation -out stats_out
```

### 5) Publication MSA plot (+ conservation strip)

```python
from bioinformatics.sequence import LoadMsa, PlotMsa

records = LoadMsa("clean.fasta")
PlotMsa(records, outFile="msa.png", maxSeqs=40, queryId="query")
```

```bash
python -m bioinformatics.sequence.msaPlot clean.fasta -out msa.png -maxSeqs 40 -queryId query
```

### 6) Phylogenetic distance tree (NJ / UPGMA)

Built from \(d = 1 -\) sequence identity (quick QC tree; not a substitute for IQ-TREE).

```python
from bioinformatics.sequence import LoadMsa, PlotPhyloTree, BuildDistanceTree

records = LoadMsa("clean.fasta")
PlotPhyloTree(records, outFile="tree.png", method="nj", maxSeqs=80)  # also writes .nwk
tree = BuildDistanceTree(records[:80], method="nj")
```

```bash
python -m bioinformatics.sequence.phyloTree clean.fasta -out tree.png -method nj -maxSeqs 80
```

### 7) View alignment (terminal or Dash)

```python
from bioinformatics.sequence import LoadMsa, DisplayAlignment

records = LoadMsa("aln.fasta")
DisplayAlignment(records)                  # first sequence as query (often "query")
DisplayAlignment(records, queryId="query") # AF3-style header
# ViewAlignment("aln.fasta", port=8050)    # needs dash + dash-bio
```

```bash
python -m bioinformatics.sequence.msaView aln.fasta
python -m bioinformatics.sequence.msaView aln.fasta query
```

### 8) Cluster by HMMER family

```bash
python -m bioinformatics.sequence.seqCluster \
  -domtblout hits.domtblout -fasta seqs.fasta -outdir cluster_out
```

Needs: `pandas`, `scikit-learn`, `umap-learn`.

### 9) UniProt

```python
from bioinformatics.sequence import GetUniprotSequence, SaveFasta

entry = GetUniprotSequence("P00533", domain="kinase")
SaveFasta([entry], "egfr.fasta")
```

```bash
python -m bioinformatics.sequence.uniProt P00533 -out egfr.fasta
python -m bioinformatics.sequence.uniProt -search 'gene:DCLK3 AND organism_id:9606'
```

Needs: `requests`.

---

## Interpreting `CleanMsa` QC plots

Both plots use **pairwise sequence identity** (fraction of matching columns) among sequences left after cleaning.

### `seqid_hist.png`

Histogram of all unique pairs \((i < j)\).

| Pattern | Meaning |
|---------|---------|
| Peak near ~0.3 | Most pairs are **distant homologs** (diverse MSA — typical for AF3) |
| Hump near ~0.6–0.7 | Subset of **closer relatives** |
| Bump near 1.0 | Near-duplicates still present |
| Single peak near 1.0 | MSA is highly **redundant** |

Compare to report **Neff**: if Neff ≪ \(N\), many rows are similar even when the hist peak is low (reweighting collapses redundancy).

### `seqid_matrix.png`

Heatmap: cell \((i,j)\) = identity of sequence \(i\) vs \(j\) (yellow = high, dark = low). Diagonal is always 1.0.

| Pattern | Meaning |
|---------|---------|
| Bright blocks on diagonal | **Groups / subfamilies** similar within the block |
| Dark rectangles off-diagonal | Those groups are **distant from each other** |
| Bright patches inside a block | **Redundancy** (near-identical sequences) |

Order follows the MSA file order (not reclustered). Use the matrix to spot clades; use the hist for overall diversity. For less redundancy: `removeParalogs=True` or a tighter `theta`.

---

## Related

- Parent overview: [`../README.md`](../README.md)
- AF3 extract MSA/templates: `machine_learning.alphafold3.ExtractMsaTemplate`
- Load AF3 structures: `utils.parseFiles.LoadAf3`
