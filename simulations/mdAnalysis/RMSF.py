#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-07-08 (Y/M/D) at 09:56
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import numpy as np

def rmsfCal(u, outpath = None, sel = "protein and name CA", avg_str = None):

    average = align.AverageStructure(u, u, select=sel, ref_frame=0).run()
    ref = average.universe
    aligner = align.AlignTraj(u, ref, select=sel, in_memory=True).run() 
    c_alphas = u.select_atoms(sel)
    
    R = rms.RMSF(c_alphas).run()
    rmsf_results = R.results.rmsf
    print (f"CA: {len(c_alphas), len(rmsf_results)}")
    if outpath:
        u.trajectory[-1] # set the trajectory to the last frame, this will write the bfactor to the last frame
        u.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
        protein = u.select_atoms(sel) # select protein atoms
        for residue, r_value in zip(protein.residues, rmsf_results):
            residue.atoms.tempfactors = r_value

        np.savetxt(outpath + '/rmsf.csv', rmsf_results, delimiter=',')
        u.atoms.write(outpath + '/bfactor_rmsf.pdb')
        #protein.write(Outfile[:-4] + '.xtc', frames='all')
        #protein.write(Outfile[:-4] + '.gro', reindex=False)

        avg_universe = mda.Merge(u.select_atoms("protein or resname THP1 or resname SP2 "))  # Create a new Universe for the sel
        avg_universe.add_TopologyAttr('tempfactors')
        for residue, r_value in zip(avg_universe.residues, rmsf_results):
            residue.atoms.tempfactors = r_value
        avg_universe.atoms.write(outpath + '/AvgProtStrBeta.pdb')
    return rmsf_results, u




path = "/mnt/mydrive/alpfaFold3/dynamics/DCLK/DCLK_KD_tail/MD/DCLK_KD_tail/dclk1_amsa_atemp"
sel = "protein and backbone and resid 390:690"
u = mda.Universe(path + "/equilibration/" + 'min.gro', path + "/rep1/" + 'md.xtc', in_memory=True, in_memory_step=10)
rmsf, u1 = rmsfCal(u, outpath = path , sel = sel)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(3.5, 2.2), dpi = 600, constrained_layout=True)
c_alphas = u.select_atoms(sel)
ax.plot(c_alphas.resids, rmsf, color = "black", label =path.split("/")[-2])
ax.legend(fontsize = 6, ncol=2);
from matplotlib.ticker import (MultipleLocator)
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(MultipleLocator(25))
#ax.yaxis.set_major_locator(MultipleLocator(5))
#ax.yaxis.set_minor_locator(MultipleLocator(2.5))
ax.tick_params(axis='both', which='major', labelsize=8)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.set_xlabel("Residue index", fontsize=10)
ax.set_ylabel("RMSF (\u00c5)", fontsize=10)
fig.savefig(path +  "/RMSF.png", dpi = 600)

