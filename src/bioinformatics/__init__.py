"""
bioinformatics — sequence analysis and structure visualization helpers.

Subpackages
    sequence/   MSA I/O, cleaning, conservation stats, clustering, UniProt
    structure/  PyMOL figure recipes, Hex docking

Quick import
    from bioinformatics.sequence import LoadMsa, PlotConservation, CleanMsa
"""

from . import sequence

__all__ = ["sequence"]
