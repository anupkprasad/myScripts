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
from MDAnalysis.transformations.positionaveraging import PositionAverager

def CalRmsf(u, outpath = None, selRMSF = "protein and backbone", selAlign = "protein and backbone and resid 1:258", ref:str = 0):    
    if ref==0:
        align.AlignTraj(u, reference=u, ref_frame=0, select=selAlign, in_memory=True).run()
    elif ref=="average":
        average = align.AverageStructure(u, u, select=selRMSF, ref_frame=0).run()
        alignuni = average.universe
        align.AlignTraj(u, reference=alignuni, select=selAlign, in_memory=True).run()
    elif ref:
        u = mda.Universe(ref, u.trajectory)
        align.AlignTraj(u, reference=u, select=selAlign, in_memory=True).run()


    rmsf_atoms = u.select_atoms(selRMSF)
    R = rms.RMSF(rmsf_atoms).run()
    rmsf_results = R.results.rmsf
    print (f"CA: {len(rmsf_atoms), len(rmsf_results)}")
    if outpath:
        np.savetxt(outpath + '/rmsf.csv', rmsf_results, delimiter=',')
        u.trajectory[0]
        u.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
        protein = u.select_atoms(selRMSF)
        for residue, r_value in zip(protein.residues, rmsf_results):
            residue.atoms.tempfactors = r_value
            print(f"Residue {residue.resname} {residue.resid} RMSF: {r_value:.2f} A")
        
        u.atoms.write(outpath + '/AvgProtStrBeta.pdb') #write the bfactor to the first frame of trajectory
    return rmsf_results





path = "/mnt/mydrive/alpfaFold3/dynamics/DCLK/DCLK_KD_tail/dclk3_amsa_atemp"
u = mda.Universe(path + '/frm0.pdb', path + '/traj.trr', in_memory=True, in_memory_step=1)

rmsf = CalRmsf(u, outpath = path)


fdfsfd


align.AlignTraj(u, u, select="protein and backbone and resid 1:258", ref_frame=0).run()
#u.trajectory[-1]
#u.atoms.write(path + '/lastframe.pdb')
#rmsf = CalRmsf(u, outpath = path)
#average = align.AverageStructure(u, u, select="protein and backbone and resid 1:258", ref_frame=0).run()
# ref = average.universe
# ref.atoms.write(path + '/test.pdb')

# Select backbone atoms
backbone = u.select_atoms("protein and backbone")

# Initialize array to accumulate positions
n_atoms = backbone.n_atoms
sum_positions = np.zeros((n_atoms, 3))
n_frames = 0

# Loop over trajectory
for ts in u.trajectory:
    sum_positions += backbone.positions
    n_frames += 1

# Compute average position
mean_positions = sum_positions / n_frames

print("Mean position of backbone atoms (shape: Nx3):")
print(mean_positions)


# Create a new universe with backbone atoms
avg_universe = mda.Merge(backbone)
avg_universe.atoms.positions = mean_positions
avg_universe.atoms.write(path + "/avg_backbone.pdb")













lfdlf
from MDAnalysis.analysis import rms, align

average = align.AverageStructure(u, u, select='protein and name CA',
                                 ref_frame=0).run()
ref = average.universe
aligner = align.AlignTraj(u, ref,
                          select='protein and name CA',
                          in_memory=True).run()

c_alphas = u.select_atoms('protein and name CA')
R = rms.RMSF(c_alphas).run()








jkdfbksdf




def CalRmsf(u, outpath = None, sel = "protein and backbone", ref:str = None):
    if ref:
        u = mda.Universe(ref, u.trajectory)
    else:   
        average = align.AverageStructure(u, u, select=sel, ref_frame=0).run()
        ref = average.universe
        print (f"Average structure saved to {outpath + '/average.pdb'}")

    #aligner = align.AlignTraj(u, ref, select=sel, in_memory=True).run() 
    rmsf_atoms = u.select_atoms(sel)

    R = rms.RMSF(rmsf_atoms).run()
    rmsf_results = R.results.rmsf
    print (f"CA: {len(rmsf_atoms), len(rmsf_results)}")
    if outpath:
        u.trajectory[-1] # set the trajectory to the last frame, this will write the bfactor to the last frame
        u.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
        protein = u.select_atoms(sel) # select protein atoms
        for residue, r_value in zip(protein.residues, rmsf_results):
            residue.atoms.tempfactors = r_value

        np.savetxt(outpath + '/rmsf.csv', rmsf_results, delimiter=',')
        u.atoms.write(outpath + '/bfactor_rmsf.pdb')

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

