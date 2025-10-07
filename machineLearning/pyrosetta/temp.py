#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-08-13 (Y/M/D) at 14:12
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

import sys
import os
import numpy as np
sys.path.append('/home/anup/myScripts/machineLearning/pyrosetta')

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.core.pose import Pose, append_pose_to_pose
from freeEnePyrosetta import minimizePose, computeBindingEnergy
init()

pose = pose_from_pdb("/mnt/mydrive/pseudokinase/James_murphy/AF3_ions/atp_f2ygr7_kd1__amsa_atemp/pdbs/0_min.pdb")

protein_pose = pose.split_by_chain(1)  # protein chain A
ligand_pose = pose.split_by_chain(2)   # ligand chain L
ion_pose = pose.split_by_chain(3)     # ion chain I

ligand_ion_pose = Pose()
ligand_ion_pose.assign(ligand_pose)
append_pose_to_pose(ligand_ion_pose, ion_pose)


res_names = [ligand_ion_pose.residue(i).name() for i in range(1, ligand_ion_pose.total_residue() + 1)]
print(res_names)


if pose.num_chains() < 3:
    raise ValueError("Pose must have at least 3 chains (protein, ligand, ion) for binding energy calculation.")

from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
iam = InterfaceAnalyzerMover("L")  # replace "X" with the chain ID of the ligand block; or define a jump if protein–ligand via docking jump
iam.set_pack_separated(True)       # repack both bound and separated states
iam.apply(pose)
print("Interface ΔG (REU):", iam.get_interface_dG())

scorefxn = get_fa_scorefxn()  # ref2015 by default
scorefxn(pose)
scorefxn.show(pose)  # look at fa