# Python script to calculate oligomerization state for each monomer chain
import numpy as np
import copy as cp
import subprocess as sp
import os as os
import shutil as sh
import MDAnalysis as mdana
import sys
from MDAnalysis.analysis.distances import distance_array
import networkx as nx
import pandas as pd
import mdtraj as md
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from natsort import natsorted # pip install natsort
import glob
from datetime import datetime

# import MDAnalysis
# from MDAnalysis.tests.datafiles import PDB, XTC
#### writing pdb
# protein = uni1.select_atoms("protein")
# with MDAnalysis.Writer("protein_test.pdb", protein.n_atoms) as W:
#     for ts in uni1.trajectory[0:2]:
#         W.write(protein)
#################  Defination
def aggregate_analysis(uni, sel_all_atoms, Min_Distance, outpath):
    #####  sel_all_atoms : of which we have to analyze  
    uni = uni
    n,t = list(enumerate(uni.trajectory))[0]
    box = t.dimensions[:6]
    frames =  len(uni.trajectory)
    #### all atoms selection for which aggre analyz, (protein/resname TFE) as input
    all_atoms = uni.select_atoms(sel_all_atoms)
    num_molecule = len(all_atoms.fragments)  ### number of total molecules
    #atoms = len(all_atoms.residues[0].atoms)
    #AminoAcids = len(protein.segments[0].residues)
    # isum=1
    # atoms_list=[]
    # atomsperAminoAcid=[]
    # residue_list=[]
    
    atom_Groups = list(all_atoms.fragments)
    fileout1 =  open(outpath +'oligomer-groups.dat','w')
    fileout2 =  open(outpath + 'oligomer-states.dat','w')
    
    for tt in uni.trajectory:  #[-110::10]:
        fileout1.write (str(tt.frame) + '\t')
        fileout2.write (str(tt.frame) + '\t')
        
        mySet = set([])
        graph = []
        # m_start=0
        # count = 0
        for i in range(num_molecule-1):
            for j in range(i+1,num_molecule):
                dist = distance_array(atom_Groups[i].positions,atom_Groups[j].positions,box).min()
                if(dist <= Min_Distance):
                    my_tuple = (i,j)
                    mySet.add(my_tuple)                
        graph = nx.from_edgelist(mySet)   
        for i in range(num_molecule):
            if i not in list(graph):
                fileout1.write ('['+ str(i)+']' + '\t')
                fileout2.write ('1' + '\t')
            else:
                fileout1.write (str(list(nx.node_connected_component(graph, i))) + '\t')
                fileout2.write (str(len(list(nx.node_connected_component(graph, i)))) + '\t')
        fileout1.write ('\n')
        fileout2.write ('\n')
    fileout1.close()
    fileout2.close()
    
    # Get oligomerization data
    
    OligoStates = [[0 for z in range(num_molecule)] for x in range(frames+1)]
    file = open(outpath + "oligomer-groups.dat",'r')
    line = file.readline()
    j = 0
    while line:
        temp = line.split('\t')
        for k in range(num_molecule):
            OligoStates[j][k] = temp[k + 1][1:-1].split(',')
        j += 1
        line = file.readline()
    file.close
    
    #Highest Oligomer size in each frame
    states=open(outpath + 'oligomer-states.dat', 'r')
    ter=states.readlines()[0:frames+1]
    
    result=[]
    for freq in (ter):
        	result.append([int(hist) for hist in freq.strip().split('\t')[1:]])
    
    fileout3 = open (outpath + 'oligo-highest-size.dat', 'w')
    for oli_count in range(len(ter)):
    	fileout3.write("{} {} {}\n".format(oli_count, '\t', np.max(result[oli_count])))
    fileout3.close()
    
    # Block Average
    size_data = np.loadtxt(outpath + 'oligomer-states.dat')
    window = 25                                         	# specify over how many frames the running average is to be calculated
    weights = np.repeat(1.0,window)/window
    size_data_m = np.convolve(size_data[:,1],weights,'valid')
    
    fileout4 = open(outpath + 'oligo-block-average.dat', 'w')
    for t,b in enumerate(size_data_m):
        	fileout4.write("{} {} {}\n".format(t, '\t', b))
    fileout4.close()



import subprocess
import numpy as np
from MDAnalysis import Universe, Writer
import tempfile
import sys
import os
def ss_traj(uni, outfile = "ss_traj.dat", outpath = os.getcwd()):
    u = uni
    protein = u.select_atoms("protein")
    output=np.array([])
    line1=np.array([])
    sam=np.array([0.0])
    #sam = 0.0
    lsw=0
    icnt=0
    list_x = ["B","E","G","H","I","T","S"] # list of secondary structure types
    ss_allinfo_allframe = []
    for ts in u.trajectory:
        temp_file = tempfile.NamedTemporaryFile(dir='/dev/shm', suffix='.pdb',
                                                delete=False)
        pdb = Writer(temp_file.name, multiframe=False)
        pdb.write(protein)
        pdb.close()
        output = subprocess.check_output(['dssp', temp_file.name])
        line1=output.split(b"\n")
        ss_dict = {e:0 for e in list_x} # Dicionnary of secondary structure count
        ss_allinfo = []
        for line in line1:
           if line.find(b"#  RESIDUE") > 0:
               lsw=1
               continue
           if (lsw == 1 and line[13:14].decode("utf-8") != "!"):
                resn=line[0:5]
                ss_allinfo.append(line[16:24])
                seco=line[16:17]
                
                if ss_dict.__contains__(seco):
                    ss_dict[seco.decode("utf-8")] += 1  #### binary to string seco.decode("utf-8")
        lsw = 0
        icnt=icnt+1
        os.remove(temp_file.name)
        ss_allinfo_allframe.append(ss_allinfo)
        
        
    ss_allres_allframes = []
    for frame in ss_allinfo_allframe:
        ss_allres = ''
        for res in frame:
            if res.decode("utf-8").split(" ")[0] in ["S", ""]:
                ss_allres = ss_allres + "C"
            else:
                ss_allres = ss_allres + res.decode("utf-8").split(" ")[0]
        ss_allres_allframes.append(ss_allres[:-1])
   
    file = open(outpath+outfile,"w")
    for line in ss_allres_allframes:
        file.write(line)
        file.write("\n")
    file.close()
    #for i in arange(sam):
    #    sam[i]=float(sam[i])/float(icnt)
    #print sam
    return ss_allres_allframes



def ss_percent(filename = "ss_traj.dat", pep_name = "all"):
    """
    file contain the tml of all peptide in line(each line is frame)
    
    pep_name = 'all' or [1, 2 ,3 , 5] ###{list of peptides}
    
    return = array3 ### containg percent of ss
    """

    f = open(filename, "r")
    ### read and make list of line(\n) and remove space by replace()###########
    
    data = f.read().replace(" ", "").split("\n")[:-1]
    
    if pep_name == 'all':
        data = data
    else:    
        ############### array for selected peptide  ########
        #array1 = [[] for i in range(len(data))]
        array1 = []
        for i in range(len(data)):
            data1 = ''
            for j in pep_name:
                c1 = 0 + 17*(j-1)  ##### j-1 is index of jth peptide
                c2 = 17 + 17*(j-1)
                data1 = data1 + data[i][c1:c2]
            
            array1.append(data1)
        data = array1
    
    
    array1 = [[] for i in range(len(data))]
    for i in range(len(data)):
        for line in data[i]:
            data1 = line
            array1[i].append(data1)
            
    ############ array2 = ({row = all residue}, {column = frame}) ##########
    Tot_resid = len(array1[0])
    array2 = np.arange(len(data)*Tot_resid).reshape(Tot_resid, len(data))
    array2 = array2.astype(str)
    
    for i in range(len(data)):  # iteration of frames
        for j in range(Tot_resid):  # iteration of total resid
            line = data[i][j]
            (array2[:, i]) = line
    
    ############ array3 for percent of secondary structure ########
    
    array3 = np.zeros(len(array1)*8).reshape(len(array1), 8)
    
    for j in range(len(array1)):
        array3[j, 1] = (array1[j].count("H")/Tot_resid)*100
        array3[j, 2] = (array1[j].count("G")/Tot_resid)*100
        array3[j, 3] = (array1[j].count("I")/Tot_resid)*100
        array3[j, 4] = (array1[j].count("E")/Tot_resid)*100
        array3[j, 5] = (array1[j].count("B")/Tot_resid)*100
        array3[j, 6] = (array1[j].count("T")/Tot_resid)*100
        array3[j, 7] = (array1[j].count("C")/Tot_resid)*100
        array3[j, 0] = indexs = j
    
    ########## saving data for gnuplot of ss percent
    # savetxt(filepath + filename + "gnuploting_prcnt_ss.dat",
    #         array3, delimiter=' ', fmt='%s')
    return(array3)





"""
initialization of code
"""

curdir = os.getcwd()
datapath = "/media/anupkumar/Backup Plus/simulation/Reeju/simulations_spacetime/Reeju/4uperin/"
job_dir_list = ["10_nosalt", "15_nosalt", "20_nosalt","40_nosalt"]
job_dir_list = ['0_nosalt']
#job_dir_list = ["0_nosalt"]

set_list =["set1", "set2", "set3"]
set_list =["set2"]
current_date_and_time = datetime.now()
logfile =  open(str(current_date_and_time)[0:10] + 'aggregate_log.dat','w')

sel_all_atoms = "protein"
Min_Distance = 5

logfile.write ("#############    Aggregate analysis:"+   str(current_date_and_time) + "    ########### \n"
              +"Aggregate analysis:    " + sel_all_atoms + '\t' +"Cou-off:    " +str(Min_Distance))


for job in job_dir_list:
    for cur_set in set_list:
        print("\n\n\n\n" +"##############" + job + "    " + cur_set  +  "##############")
        cur_datapath = datapath + job
        cur_outpath = curdir + "/" + job + "/" + cur_set  + "/aggregate_analysis/"+ sel_all_atoms+"/"
        logfile.write ("\n\n\n\n" +"##############    " + job + "   " + cur_set  +  "     ##############"
                       + "\n\n\n" + "datapath:"+ cur_datapath
                       +"\n" + "Output:" + cur_outpath + "\n")
        if not os.path.exists(cur_outpath):
            os.makedirs(cur_outpath)
        for root, dirs, files in os.walk(cur_datapath, topdown=False):
            for name in files:
                if name.endswith(".psf"):
                    psf = os.path.join(root, name)
            matches = []
            if cur_set == "set1":
                for name in dirs:
                    if name.startswith(("prod1", "Prod1")):
                        matches.append(name)
            if cur_set == "set2":
                for name in dirs:
                    if name.startswith(("prod2", "Prod2")):
                        matches.append(name)
            if cur_set == "set3":
                for name in dirs:
                    if name.startswith(("prod3", "Prod3")):
                        matches.append(name)
            if cur_set == "set4":
                for name in dirs:
                    if name.startswith(("prod4", "Prod4")):
                        matches.append(name)
        matches = natsorted(matches) ### sort the value according to number not by char
        #print(matches)
        traj = []
        for dcd_dir in matches:
            dcd_dir = os.path.join(root,dcd_dir)
            traj.append(glob.glob(dcd_dir + "/*dcd")[0])
        #print(traj)
        logfile.write (psf + "\n")
        for dcdfile in traj:
            logfile.write (dcdfile + "\n")
        
        #### Create universe ######
        uni = mdana.Universe(psf,traj, in_memory = True, in_memory_step=10)
        aggregate_analysis(uni = uni, sel_all_atoms = sel_all_atoms, Min_Distance = Min_Distance, outpath= cur_outpath)
        #ss_allinfo_allframe1 = ss_traj(uni, outfile="ss_traj.dat", outpath= "/media/anupkumar/Backup Plus/simulation/Reeju/analysis_anup/4uperin/10_nosalt/set3/")
                      
logfile.close()

      