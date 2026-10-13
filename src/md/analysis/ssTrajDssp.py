#!/usr/bin/env python

# Author: Guillaume Bouvier -- guillaume.bouvier@pasteur.fr
# https://research.pasteur.fr/en/member/guillaume-bouvier/
# 2016-08-26 13:20:35 (UTC+0200)

# DSSP_SECONCARY.PY

# Trajectories of Secondary Structure by DSSP
#    The keys  [" ","B","E","G","H","I","T","S"]
#    change to   0   1   2   3   4   5   6   7 
# B	Beta bridge
# E	Strand
# G	Helix-3
# H	Alpha helix
# I	Helix-5
# T	Turn
# S	Bend
# 
# python dssp_secondary.py pdbfile dcdfile
#
# Adapted by Guillaume Bouvier
# from Naoyuki Miyashita original script (see: https://gist.github.com/yukimya/b774e104dc8308c19ce3)

import subprocess
import numpy as np
from MDAnalysis import Universe, Writer
import tempfile
import sys
import os

def main_file(uni, outfile = "ss_traj.dat"):
    u = uni
    protein = u.select_atoms("protein")
    output=np.array([])
    line1=np.array([])
    sam=np.array([0.0])
    #sam = 0.0
    lsw=0
    icnt=0
    list_x = ["B","E","G","H","I","T","S"] # list of secondary structure types
    f=open(outfile,"w")
    for ts in u.trajectory:
        temp_file = tempfile.NamedTemporaryFile(dir='/dev/shm', suffix='.pdb',
                                                delete=False)
        pdb = Writer(temp_file.name, multiframe=False)
        pdb.write(protein)
        pdb.close()
        output = subprocess.check_output(['dssp', temp_file.name])
        line1=output.split("\n")
        ss_dict = {e:0 for e in list_x} # Dicionnary of secondary structure count
        for line in line1:
           if line.find("#  RESIDUE") > 0:
               lsw=1
               continue
           if lsw == 1:
               resn=line[0:5]
               seco=line[16:17]
               if ss_dict.has_key(seco):
                   ss_dict[seco] += 1
        for i in list_x:
            f.write(str(ss_dict[i]))
            f.write(" ")
        f.write("\n")
        lsw = 0
        icnt=icnt+1
        os.remove(temp_file.name)
    f.close()
    #for i in arange(sam):
    #    sam[i]=float(sam[i])/float(icnt)
    #print sam

if __name__ == '__main__':
    psf = sys.argv[1]
    dcd = sys.argv[2]
    outfile = "dssp.out"
    main_file(psf, dcd, outfile)
