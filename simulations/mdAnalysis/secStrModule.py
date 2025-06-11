"""
Created on Mon Dec 21 13:53:15 2020
@author: anupkumar
anupkprasad121@gmail.com
PhD @IITB-Monash Research Academy
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from tsmoothie.smoother import *
import sys

def running_mean(y, N):
    smoother1 = ConvolutionSmoother(window_len=N, window_type='ones')
    smoother1.smooth(y)
    return smoother1.smooth_data[0]


class Sec_Str(object):
    """
    Works for input files having status of secondry structure of each residue of proteins (all proteins are same).
    """
    def __init__(self, filename = 'ss_traj.dat', pep_idx = "all", res_num = 17):
        self.filename = filename
        self.pep_idx = pep_idx
        self.res_num = res_num
        self.ss_order = {r"$\alpha$-helix":"H",r"$3_{10}$-helix":"G",r"$\pi$-helix":"I",r"$\beta$-sheet":"E",r"$\beta$-bridge":"B","turn":"T","coil":"C",}
        self.data_all = open(filename, "r").read().replace(" ", "").split("\n")[:-1]
        
    def data_pep_idx(self):
        """
        Returns : array of sec_structure of selected pep_idx
        -------
       """
        pep_idx = self.pep_idx
        data_all = self.data_all
        res_num = self.res_num
        
        if pep_idx == 'all':
            data = data_all
        else:    
            ############### array for selected peptide  ########
            arr = []
            for i in range(len(data_all)):
                data1 = ''
                for j in pep_idx:
                    c1 = 0 + res_num*(j)  ##### j is index of jth peptide
                    if c1 >= len(data_all[0]):
                        sys.exit("pep_idx index out of range")
                    c2 = res_num + res_num*(j)
                    data1 = data1 + data_all[i][c1:c2]
                arr.append(data1)
            data = arr  
        array = []
        for row in data:
            array.append(list(row)) ## list(row) is easy in array conversion
        array = np.array(array) 
        return(array)
    
    
    def ss_percent(self):
        """
        Returns: np array of sec_structure % over trajectories for pep_idx
       
        """
        array = self.data_pep_idx()
        array_sp = np.zeros(len(array)*8).reshape(len(array), 8)
        sso = list(self.ss_order.values())
        all_res = len(array[0])
        for j in range(len(array)):
            array_sp[j, 0] = j + 1  #### time frame
            for i, st in enumerate(sso):
                array_sp[j, i+1] = (np.count_nonzero(array[j] == st)/all_res)*100            
        return(array_sp)   

    def ss_freq(self, time_i = 0, time_f =-1):
        """/media/anupkumar/Backup Plus/simulation/ML_sim/LK_e_40k_num_77
        Returns: np array of sec_structure frequency for trajectories of pep_idx
       
        """
        array = self.data_pep_idx()
        self.time_i = time_i
        self.time_f = time_f
        res_num = self.res_num
        if not (time_i == 0 and time_f == -1):
            array = array[time_i:time_f]
        array_ssf = np.zeros((res_num, 8))
        all_res = len(array[0])
        divisor_factor = len(array[0:,0])*(all_res/res_num) ### 
        array = array.T
        sso = list(self.ss_order.values())
        for j in range(res_num):
            array_ssf[j, 0] = j + 1  #### time frame
            for i, st in enumerate(sso):
                array_ssf[j, i+1] = (np.count_nonzero(array[j:all_res:res_num] == st)/divisor_factor)*100            
        return(array_ssf)
    
    def ss_timeavg (self,time_i = 0, time_f = -1):
        """array of ss_percent array of each peptide,
        Ini = initial value , Fin = final value , Ran = range for which taking avarge for both Ini Fin
        """
        array_sp = self.ss_percent()
        self.time_i = time_i
        self.time_f = time_f
        
        if (time_i == 0 and time_f == -1):
            array_sp = array_sp
        else:
            array_sp = array_sp[time_i:time_f]
            
        div = len(array_sp)
        array_sp = np.sum(array_sp, axis = 0)/div
        array_sp[0] = time_f - time_i
        return(array_sp)
    
    def ss_heatmap (self,time_i = 0, time_f = -1):
        """
        data conversion for heat map
        """
        array = self.data_pep_idx()
        self.time_i = time_i
        self.time_f = time_f
        
        if not (time_i == 0 and time_f == -1):
            array = array[time_i : time_f]
        sso = list(self.ss_order.values())
        for i, st in enumerate(sso):
            array[array==st] = i+1
        array = array.astype(int).T
        return array
    
    def plot_basic_ss(self, out_file, time_i = 0, time_f = -1):
        self.time_i = time_i
        self.time_f = time_f 
        
        fig, axs = plt.subplots(3,1, figsize=(3.5,6), constrained_layout = True)
        ssf = self.ss_freq(time_i = time_i, time_f = time_f)
        ssp = self.ss_percent()
        ss_hmp = self.ss_heatmap(time_i = time_i, time_f = time_f)
        axs[1].plot(ssp[:,1] , color = "red", alpha = 0.3)
        axs[1].plot(running_mean(ssp[:,1], 10) , label = r"$\alpha$-helix", color = "red")
        axs[1].plot(running_mean(ssp[:,2], 10), label = r"$3_{10}$-helix", color = "green")
        axs[1].plot(ssp[:,2], color = "green", alpha = 0.3)
        axs[1].plot(ssp[:,4], label = r"$\beta$-sheet", color = "blue")
        axs[2].plot(ssf[:,1], label = r"$\alpha$-helix", color = "red",  marker = 'o')
        axs[2].plot(ssf[:,2], label = r"$3_{10}$-helix", color = "green", marker = 'o')
        axs[2].plot(ssf[:,4], label = r"$\beta$-sheet", color = "blue", marker = 'o')
                
        #ss_hmp[0:7,-1] = [0,1,2,3,4,5,6]  #### this line is intentionally added to maintain the cmap for all ss range
        cmap = (matplotlib.colors.ListedColormap(['red','limegreen','magenta','blue','tan','cyan','white'], 7))
        c = axs[0].pcolormesh(ss_hmp, cmap = cmap)
        evenly_dist = (ss_hmp.max()-ss_hmp.min())/7
        t = np.arange(ss_hmp.min()+(evenly_dist/2), ss_hmp.max()+(evenly_dist/2), evenly_dist)
        cbar = fig.colorbar(c, ax = axs[0], ticks= t, location = "top",fraction=1, aspect=50)
        cbar.ax.set_xticklabels( list(self.ss_order.keys()), fontsize = 6)  # horizontal colorbar
        y_label = ["Residue index", "Secondary structure (%)", "Frequency (%)"]
        x_label = ["Time (ns)", "Time (ns)", "Residue index"]
        for i in range(3):
            axs[i].set_ylabel(y_label[i])
            axs[i].set_xlabel(x_label[i])
            if i != 0:
                axs[i].legend(fontsize = 6)
        from matplotlib.ticker import MaxNLocator
        axs[0].yaxis.set_major_locator(MaxNLocator(integer=True))
        fig.savefig(out_file, dpi = 600)
        

