#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 10:16:36 2022
@author: anupkumar
"""
import matplotlib.pyplot as plt
import os
import numpy as np

class ABFTraj (object):
    def __init__(self, datapath, dirPrefix = "prod"):
        self.datapath = datapath
        self.dirPrefix = dirPrefix
        self.data = self.getData()
        
    def getData (self):
        datapath = self.datapath
        dirPrefix = self.dirPrefix
        listDirs = next(os.walk(datapath))[1] # list of all dir
        data = []
        prodlist = []
        for d in listDirs:
            if d.startswith(dirPrefix):
                prodlist.append(d)
        print("Total prods = " + str(len(prodlist)) + "\n")
        prodlist.sort()
        for prod in prodlist:
            curr_path = datapath + prod + '/'
            for file in os.listdir(curr_path):
                if file.endswith(".traj"):
                    print(file)
                    cur_data = open(curr_path + file, 'r')
                    cur_data = cur_data.readlines()[1:]
                    for line in cur_data:
                        line = line.split(' ')
                        line_val = []
                        for val in line:
                            if val != '':
                                line_val.append(val)
                        if line_val[0] != '#':                          
                            frame = float(line_val[0])/500000
                            alpha_rmsd = float(line_val[1])
                            beta_rmsd = float(line_val[2])
                            line_val_all = [frame,alpha_rmsd,beta_rmsd]
                            data.append(line_val_all)
                            
        data_array = np.asarray(data)
        return (data_array)
    
    def getTimeofRCs (self, RC1_range, RC2_range,  grid = 0.5):
        data_array = self.data
        step = grid                                  
        ###### reaction cordinates points for which frames will be searched
        match_1 = np.arange(RC1_range[0], RC1_range[1]+step, step)
        match_2 = np.arange(RC2_range[0], RC2_range[1]+step, step)
        match_points = []
        for i in match_1:
            for j in match_2:
                point = [i,j]
                match_points.append(point)
        ###### finding frames corresponding to reaction cordinates points
        matched_value = []
        for p in match_points:   
            for i in range(len(data_array[:,1])):
                if round(data_array[:,1][i],0) == p[0]:
                    if round(data_array[:,2][i],0) == p[1]:
                        matched_value.append(i)
           
        time_of_rcs = []
        for i in matched_value:
            time_of_rcs_line = [data_array[i,0], data_array[i,1], data_array[i,2]]
            time_of_rcs.append(time_of_rcs_line)
            
        time_of_rcs = np.asarray(time_of_rcs)
        return (time_of_rcs)
    
    def plotTraj (self, RC1_range, RC2_range):
        data_array = self.data
        time_of_rcs = self.getTimeofRCs(RC1_range, RC2_range)
        datapath = self.datapath
        fig, axs = plt.subplots(1, 1, figsize=(3.5,2.5), constrained_layout = True)
        axs.plot(data_array[:,0], data_array[:,1], linewidth = 0.1, color = 'black', label = 'alpha-rmsd', alpha = 0.5)
        axs.plot(data_array[:,0], data_array[:,2], linewidth = 0.1, color = 'red', label = 'beta-rmsd', alpha = 0.5)
        axs.scatter(time_of_rcs[:,0], time_of_rcs[:,1], s= 20, c ='b')
        axs.scatter(time_of_rcs[:,0], time_of_rcs[:,2], s= 20, c ='r')
        axs.legend(loc="upper right", fontsize=5)
        plt.savefig(datapath + 'traj.png', dpi=300)
        
    def plotSampling (self):
        data_array = self.data
        datapath = self.datapath
        fig, axs = plt.subplots(1, 1, figsize=(3.5,2.5), constrained_layout = True)
        axs.scatter(data_array[:,1], data_array[:,2], s= 2, c ='black', marker = ".")
        plt.savefig(datapath + 'sampled.png', dpi=300)
    
    
    
 ###################### Code initialisation ###########################
   
                  
# datapath = '/media/anupkumar/Backup Plus/simulation/anup_data/u3.x_.15M/concen_effect/abf/dimer_u3.5R7AK8A/x_0_10_y_0_10/'
# data = ABFTraj(datapath).data
# RC1_range = [11.5, 12.5]
# RC2_range = [6,7.5]
# ABFTraj(datapath).plotSampling()

def plotSamplingAllQuad (pathlist):
    fig, axs = plt.subplots(1, 1, figsize=(3.5,2.5), constrained_layout = True)
    color = ["black", "red", "blue", "cyan"]
    for i, datapath in enumerate(pathlist):
        data_array = ABFTraj(datapath).data
        print(len(data_array))
        axs.scatter(data_array[:,1], data_array[:,2], s= 2, c = color[i], marker = ".", alpha = 0.5)
        axs.plot([10,10], [0,20], color = "black")
        axs.plot([0,20], [10,10], color = "black")
    axs.set_xlim(0,20)
    axs.set_ylim(0,20)
    plt.savefig(path + 'sampledAllQuad.png', dpi=300)

path = "/media/anupkumar/Backup Plus/simulation/anup_data/u3.x_.15M/concen_effect/abf/dimer_u3.5D4A/"
pathlist = [path + "x_0_10_y_0_10/", path + "x_0_10_y_10_20/", path + "x_10_20_y_0_10/", path + "x_10_20_y_10_20/"]

plotSamplingAllQuad(pathlist)