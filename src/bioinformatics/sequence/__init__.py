"""
sequence — MSA I/O, cleaning, stats/conservation, viewing, clustering, UniProt.

Modules (open file for full examples)
    msaIO.py       LoadMsa / WriteMsa / A3mToFasta / CmaToFasta
    msaStats.py    ColumnConservation / PlotConservation / PlotEntropy
    msaPlot.py     PlotMsa (colored MSA + conservation strip)
    phyloTree.py   BuildDistanceTree / PlotPhyloTree
    msaClean.py    CleanMsa (+ QC plots)
    msaView.py     DisplayAlignment / ViewAlignment (Dash)
    seqCluster.py  ClusterByDomtblout / PlotUmap
    uniProt.py     GetUniprotSequence / SearchUniprot

See README.md in this directory for usage examples.

Typical flow for MSA / conservation
    records = LoadMsa("aln.a3m", fmt="a3m")
    CleanMsa(...)                     # optional
    PlotMsa(records, "msa.png")
    PlotPhyloTree(records[:80], "tree.png")
"""

from .msaClean import CleanMsa, FilterByCoverage, RemoveParalogs, TrimGappyColumns
from .msaIO import A3mToFasta, CmaToFasta, LoadMsa, MultiToSingleLineFasta, WriteMsa
from .msaPlot import PlotMsa
from .msaStats import (
    ColumnConservation,
    ColumnEntropy,
    EffectiveSeqCount,
    GapFraction,
    PlotColumnProfile,
    PlotConservation,
    PlotEntropy,
    PlotGapProfile,
    PlotSeqIdMatrix,
    SeqIdentityMatrix,
)
from .msaView import DisplayAlignment, SortBySimilarity, ViewAlignment
from .phyloTree import BuildDistanceTree, DistanceFromMsa, PlotPhyloTree, WriteNewick
from .seqCluster import ClusterByDomtblout, KmerEmbedding, ParseDomtblout, PlotUmap
from .uniProt import GetUniprotIds, GetUniprotSequence, SaveFasta, SearchUniprot

__all__ = [
    "LoadMsa",
    "WriteMsa",
    "A3mToFasta",
    "CmaToFasta",
    "MultiToSingleLineFasta",
    "ColumnConservation",
    "ColumnEntropy",
    "GapFraction",
    "SeqIdentityMatrix",
    "EffectiveSeqCount",
    "PlotConservation",
    "PlotEntropy",
    "PlotGapProfile",
    "PlotSeqIdMatrix",
    "PlotColumnProfile",
    "PlotMsa",
    "BuildDistanceTree",
    "DistanceFromMsa",
    "PlotPhyloTree",
    "WriteNewick",
    "CleanMsa",
    "FilterByCoverage",
    "RemoveParalogs",
    "TrimGappyColumns",
    "DisplayAlignment",
    "SortBySimilarity",
    "ViewAlignment",
    "ClusterByDomtblout",
    "ParseDomtblout",
    "KmerEmbedding",
    "PlotUmap",
    "GetUniprotSequence",
    "GetUniprotIds",
    "SearchUniprot",
    "SaveFasta",
]
