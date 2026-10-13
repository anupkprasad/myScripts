"""
Created on Mon Dec 21 13:53:15 2020

@author: anup
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import plotly
import plotly.graph_objects as go
import plotly.express as px


class PMF_analysis:
    def __init__(self, filename, x_range=None, y_range =None, z_range = None):
        self.filename = filename
        self.x_range = x_range
        self.y_range = y_range
        self.z_range = z_range
   
    def get_dataXYZ(self):
        filename = self.filename
        x_range, y_range, z_range =self.x_range, self.y_range, self.z_range
        data = [[], [], []]
        with open(filename, 'r') as file : 
            file_   =  file.readlines()[4:] ### [4:] skip 4 lines of file
            for lines in file_ :
                if lines  != '\n' :                    
                    data[0].append(float(lines.split()[0]))
                    data[1].append(float(lines.split()[1]))
                    data[2].append(float(lines.split()[2]))
                    
        data  =  np.array(data, dtype = float).T
        row_num = len(np.unique(data[:,0]))
        column_num = len(np.unique(data[:,1]))
        data = data.reshape(row_num, column_num, -1)
        X           =    data[:,:,0]
        Y           =    data[:,:,1]
        Z           =    data[:,:,2]

        if x_range != None:
            xi_range = int(x_range[0]*2)
            xf_range = int(x_range[1]*2)
            [X,Y,Z] = [X[xi_range:xf_range], Y[xi_range:xf_range], Z[xi_range:xf_range]]
        if y_range != None:
            yi_range = int(y_range[0]*2)
            yf_range = int(y_range[1]*2)
            [X,Y,Z] = [X[:, yi_range:yf_range], Y[:, yi_range:yf_range], Z[:, yi_range:yf_range]]
        if z_range != None:
            zi = z_range[0]
            zf = z_range[1]
            Z = np.where(Z < zi, np.nan, Z)
            Z = np.where(Z > zf, np.nan, Z)
        data_final = [X,Y,Z]
        return(data_final)
    
    def get_XYZ_of_interpolated_line(self, line_points:list, bin_size = 0.5):
        data = self.get_dataXYZ()
        [X,Y,Z] =  data
        XY = []
        for i in range(len(line_points)-1):
            P1 = line_points[i] # first point 
            P2 = line_points[i+1] # second point
            print(P1,P2)
            if P1[0] < P2[0]:
                X_cur = list(np.arange(P1[0],P2[0],bin_size))
            else:
                X_cur = list(np.arange(P2[0],P1[0],bin_size))
            for x in X_cur:
                print(x)
                m = (P2[1]-P1[1])/(P2[0]-P1[0])
                y = m*(x-P1[0]) + P1[1]
                d = round(y)-y
                if -0.25 < d < 0.25:
                    y = round(y)
                elif d < 0:
                    y = round(y) + 0.5
                else:
                    y = round(y) - 0.5
                    
                XY.append([x,y])
        XYZ = []
        for [x,y] in XY:
            index_x = np.where(X==x)
            index_x = index_x[0][0]
            index_y = np.where(Y.T==y)
            index_y = index_y[0][0]
            z = Z[index_x, index_y]
            XYZ.append([x,y,z])
        XYZ = np.array(XYZ)
        return(XYZ)
        
    def annotate_track (self, points:list, XZ_output = True):
        data = self.get_dataXYZ()
        [X,Y,Z] =  data
        XY = points
        XYZ= []
        for [x,y] in XY:
            index_x = np.where(X==x)
            index_x = index_x[0][0]
            index_y = np.where(Y.T==y)
            index_y = index_y[0][0]
            z = Z[index_x, index_y]
            XYZ.append([x,y,z])
        XYZ = np.array(XYZ)
        if XZ_output == True:
            output = XYZ[:,[0,2]]
        else:
            output = XYZ[:,[0,1]]
            
        return output
        
    def get_XYZ (self, points:list, XZ_output = True):
        data = self.get_dataXYZ()
        [X,Y,Z] =  data
        XY = points
        XYZ= []
        for [x,y] in XY:
            index_x = np.where(X==x)
            index_x = index_x[0][0]
            index_y = np.where(Y.T==y)
            index_y = index_y[0][0]
            z = Z[index_x, index_y]
            XYZ.append([x,y,z])
        XYZ = np.array(XYZ)
        return XYZ
    
    def plot_contour(self, save_file, minima = True, maxima = False, saddle = False, cband = 15):
        data = self.get_dataXYZ()
        fig, ax = plt.subplots(figsize=(3.5, 2.5), constrained_layout = True)
        cmap = cm.get_cmap("bwr")
        plt.gca().set_facecolor("lightgray") 
        cs = ax.contourf(data[0],data[1],data[2], cband, cmap=cmap)
        ax.contour(data[0],data[1],data[2], cband, colors = 'black', linewidths= 0.2)
        local_minmax = self.get_minmax()
        if minima == True:
            minima_points =  local_minmax[0]
            x = [i[0] for i in minima_points]
            y = [i[1] for i in minima_points]
            z =  [i[2] for i in minima_points]
            ax.scatter(x,y, c ="black", s = 2)
            for (a,b,c) in zip(x,y,z):
                ax.text(a + 0.2, b, round(c, 2), fontsize = 4, c ="black")
        
        if maxima == True:
            maxima_points =  local_minmax[1]
            x = [i[0] for i in maxima_points]
            y = [i[1] for i in maxima_points]
            z =  [i[2] for i in maxima_points]
            ax.scatter(x,y, c ="blue", s = 2)
            for (a,b,c) in zip(x,y,z):
                ax.text(a + 0.2, b, round(c, 2), fontsize = 4, c ="blue")
                
        if saddle == True:
            saddle_points =  self.get_saddles()
            x = [i[0] for i in saddle_points]
            y = [i[1] for i in saddle_points]
            z =  [i[2] for i in saddle_points]
            print("Total saddles = " + str(len(z)))
            ax.scatter(x,y, c ="blue", s = 2)
            for (a,b,c) in zip(x,y,z):
                ax.text(a + 0.2, b, round(c, 2), fontsize = 4, color = "blue") 
        #plt.clabel(contours, inline=True, fontsize=5, colors = 'black')
        fig.colorbar(cs, ax=ax, label = "PMF (Kcal/mol)")
        
        ########## Minor tick + Major ##############
        from matplotlib.ticker import (MultipleLocator)
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.xaxis.set_tick_params(which='both', top = 'on', direction= 'out')
        ax.yaxis.set_minor_locator(MultipleLocator(1.5))
        ax.yaxis.set_major_locator(MultipleLocator(3))
        ax.yaxis.set_tick_params(which='both', right = 'on', direction= 'out')
        plt.xlabel('RMSD alpha (x)')        
        plt.ylabel('RMSD beta (y)')
        plt.savefig(save_file, dpi=600)
        
    def get_minmax(self):
        data = self.get_dataXYZ()
        local_minima = []
        local_maxima = []
        data3 = data[2]
        # Iterate over the 2D grid
        for i in range(data3.shape[0]):
            for j in range(data3.shape[1]):
                current_value = data3[i, j]
                neighbors = []
                for x in [-1, 0, 1]:
                    for y in [-1, 0, 1]:
                        ni, nj = i + x, j + y

                        if (0 <= ni < data3.shape[0]) and (0 <= nj < data3.shape[1]) and (x, y) != (0, 0):
                            neighbors.append((ni, nj))        
                # Check if the current point is less than all its neighbors
                if all(current_value < data3[n[0], n[1]] for n in neighbors):
                    local_minima.append((data[0][i, j], data[1][i, j], current_value))
                # Check if the current point is less than all its neighbors
                if all(current_value > data3[n[0], n[1]] for n in neighbors):
                    local_maxima.append((data[0][i, j], data[1][i, j], current_value))
        return [local_minima, local_maxima]
    
    def get_saddles(self):
        data3D = self.get_dataXYZ()
        data3 = data3D[2] ### 2D energy matrix
        xdiff = np.diff(data3, axis=0) # Derivative along 1st dim
        ydiff = np.diff(data3, axis=1) # Derivative along 2nd dim
    
    # Saddle points are not in the "boundary" elements
    
        xdiff1 = xdiff[1:  , 1:-1]  # Shift matrix skipping 1st row
        xdiff2 = xdiff[ :-1, 1:-1]  # Cut the last row
        ydiff1 = ydiff[1:-1, 1:  ]  # Same, for columns
        ydiff2 = ydiff[1:-1,  :-1]
    
    # A saddle points is a local maximum/minimum in one dimension (xdiff1*xdiff2 < 0) and in the other dimension (ydiff1*ydiff2 < 0), 
    # but only in the combination max1&min2 or min1&max2 (xdiff1*ydiff1 < 0)
        ind = np.where((xdiff1*xdiff2 < 0) & (ydiff1*ydiff2 < 0) &
    (xdiff1*ydiff1 < 0))
        saddle_idx = []
        saddle_points = []
        if len(ind) > 0:
            for j, x in enumerate(ind[0]):
                xid = x + 1
                yid = ind[1][j] +1
                saddle_idx.append([xid, yid])
                saddle_points.append([data3D[0][xid,yid], data3D[1][xid,yid], data3D[2][xid,yid]])
    
        return saddle_points
                    
        
    def interactive_contour (self, save_file, minima = True, maxima = False, saddle = False,msize = 10):
        data = self.get_dataXYZ()
        [X,Y,Z] =  data
        #fig = go.Figure(data=[go.Surface(z=z_data.values)])
        terrain_cmap = cm.get_cmap('terrain')
        def matplotlib_to_plotly(cmap, pl_entries):
            h = 1.0/(pl_entries-1)
            pl_colorscale = []
        
            for k in range(pl_entries):
                C = list(map(np.uint8, np.array(cmap(k*h)[:3])*255))
                pl_colorscale.append([k*h, 'rgb'+str((C[0], C[1], C[2]))])
        
            return pl_colorscale
        
        terrain = matplotlib_to_plotly(terrain_cmap, 255)

        plotly.offline.init_notebook_mode(connected=True)      
        fig = go.Figure(data=[go.Surface(colorscale=terrain,z=Z, x=X, y=Y)])
        fig.update_layout(#title = "U3.4 R7A",
            #autosize=False, width=2000, height=500,
            font_color="blue",
            font_size = 12,
            title_font_color="red",
            legend_title_font_color="green"
        )
        
        fig.update_layout(scene = dict(
                            xaxis_title='RMSD_alpha (x)',
                            yaxis_title='RMSD_beta (y)',
                            zaxis_title='PMF'))
        
        fig.update_layout(scene = dict(
                            xaxis = dict(
                                 backgroundcolor="rgb(200, 200, 230)",
                                 gridcolor="white",
                                 showbackground=True,
                                 zerolinecolor="white",),
                            yaxis = dict(
                                backgroundcolor="rgb(230, 200,230)",
                                gridcolor="white",
                                showbackground=True,
                                zerolinecolor="white"),
                            zaxis = dict(
                                backgroundcolor="rgb(230, 230,200)",
                                gridcolor="white",
                                showbackground=True,
                                zerolinecolor="white",)))       
        
        local_minmax = self.get_minmax()
        if minima == True:
            minima_points =  local_minmax[0]
            x = [i[0] for i in minima_points]
            y = [i[1] for i in minima_points]
            z =  [i[2] for i in minima_points]
            print("Total minima = " + str(len(z)))
            fig.add_scatter3d(x=x, y=y, z = z, mode='markers', 
                  marker=dict(size=msize, color='black'), name='Legendary')
        
        
        if saddle == True:
            saddle_points =  self.get_saddles()
            x = [i[0] for i in saddle_points]
            y = [i[1] for i in saddle_points]
            z =  [i[2] for i in saddle_points]
            print("Total saddles = " + str(len(z)))
            fig.add_scatter3d(x=x, y=y, z = z, mode='markers', 
                  marker=dict(size=msize, color='orange'))
            
        fig.update_layout(showlegend=True)
        
        # Note that include_plotlyjs is used as cdn so that the static site generator can read it and present it on the browser. This is not typically required.
        plotly.offline.plot(fig, filename= save_file,include_plotlyjs='cdn')
        
        





# Function to get the neighbors of a given point (i,j)
def get_neighbors(i, j, shape):
    neighbors = []

    for x in [-1, 0, 1]:
        for y in [-1, 0, 1]:
            ni, nj = i + x, j + y

            if (0 <= ni < shape[0]) and (0 <= nj < shape[1]) and (x, y) != (0, 0):
                neighbors.append((ni, nj))
    return neighbors


def get_minmax (data):
    local_minima = []
    local_maxima = []
    data3 = data[2]
    # Iterate over the 2D grid
    for i in range(data3.shape[0]):
        for j in range(data3.shape[1]):
            current_value = data3[i, j]
            neighbors = get_neighbors(i, j, data3.shape)
            
            # Check if the current point is less than all its neighbors
            if all(current_value < data3[n[0], n[1]] for n in neighbors):
                local_minima.append((data[0][i, j], data[1][i, j], current_value))
            # Check if the current point is less than all its neighbors
            if all(current_value > data3[n[0], n[1]] for n in neighbors):
                local_maxima.append((data[0][i, j], data[1][i, j], current_value))
    return [local_minima, local_maxima]
    
    
def allSaddles(data3D):
    
    data3 = data3D[2] ### 2D energy matrix
    rowmins = []
    rowmaxs = []
    colmins = []
    colmaxs = []

    for i,row in enumerate(data3):
        m = min(row)
        M = max(row)
        for j,x in enumerate(row):
            if x == m: rowmins.append((i,j))
            if x == M: rowmaxs.append((i,j))

    t = [list(column) for column in zip(*data3)] #transpose of matrix

    for j,col in enumerate(t): 
        m = min(col)
        M = max(col)
        for i,x in enumerate(col):
            if x == m: colmins.append((i,j))
            if x == M: colmaxs.append((i,j))
    ### perform set theory unione (|) and intersectiona (&)
    saddles_indices = (set(rowmins) & set(colmaxs)) | (set(rowmaxs) & set(colmins))
    saddle_points = []
    for  index in saddles_indices:
        i, j = index
        saddle_x = data3D[0][i,j]
        saddle_y = data3D[1][i,j]
        saddle_z = data3D[2][i,j]
        saddle_points.append([saddle_x,saddle_y,saddle_z])

    return saddle_points


def find_saddle(data3D):
    data3 = data3D[2] ### 2D energy matrix
    xdiff = np.diff(data3, axis=0) # Derivative along 1st dim
    ydiff = np.diff(data3, axis=1) # Derivative along 2nd dim

# Saddle points are not in the "boundary" elements

    xdiff1 = xdiff[1:  , 1:-1]  # Shift matrix skipping 1st row
    xdiff2 = xdiff[ :-1, 1:-1]  # Cut the last row
    ydiff1 = ydiff[1:-1, 1:  ]  # Same, for columns
    ydiff2 = ydiff[1:-1,  :-1]

# A saddle points is a local maximum/minimum in one dimension (xdiff1*xdiff2 < 0) and in the other dimension (ydiff1*ydiff2 < 0), 
# but only in the combination max1&min2 or min1&max2 (xdiff1*ydiff1 < 0)
    ind = np.where((xdiff1*xdiff2 < 0) & (ydiff1*ydiff2 < 0) &
(xdiff1*ydiff1 < 0))
    saddle_idx = []
    saddle_points = []
    if len(ind) > 0:
        for j, x in enumerate(ind[0]):
            xid = x + 1
            yid = ind[1][j] +1
            saddle_idx.append([xid, yid])
            saddle_points.append([data3D[0][xid,yid], data3D[1][xid,yid], data3D[2][xid,yid]])
    
    return saddle_points
    

"""
code initialization
"""
## variables

path = "/media/anupkumar/Backup Plus/simulation/anup_data/u3.x_.15M/concen_effect/abf/dimer_r7au3.5/randomcoil/prod10/"

filename    =   path + 'prod10.pmf'
z_range =  [0, 20]

save_path = path

# data = PMF_analysis(filename,x_range= [3, 18],y_range= [0, 20]).get_dataXYZ()
# get_minimas = get_minmax(data)
# test = PMF_analysis(filename,x_range= [0, 20],y_range= [0, 20]).get_saddles()
# data = PMF_analysis(filename).get_dataXYZ()

PMF_analysis(filename, x_range= [0, 20], y_range= [0, 20] , z_range= z_range).plot_contour(save_path+"contour.png", saddle= True, cband = 20)
PMF_analysis(filename,x_range= [0, 20], y_range= [0, 20], z_range= z_range).interactive_contour(save_path+"contour.html", saddle= True)

