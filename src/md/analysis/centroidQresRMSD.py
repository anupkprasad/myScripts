#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 23:38:12 2024
@author: anupkumar
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import AgglomerativeClustering
import shutil
import os
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
import seaborn as sns
import matplotlib.cm as cm
from sklearn.manifold import TSNE


class centroidQresRMSD (object):
    def __init__(self, filepath, num_clusters, fileqres = 'qres.dat',filermsd = 'rmsd.dat'):
        self.filepath = filepath
        self.num_clusters = num_clusters
        self.fileqres = fileqres
        self.filermsd = filermsd
        self.outpath = filepath + "centroids/"+ "centroid_" + str(num_clusters)
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)
        
    def getData(self):
        filename_qres= self.filepath + self.fileqres
        filename_rmsd= self.filepath + self.filermsd
        
        qres = open(filename_qres, "r").readlines()
        data1 = []
        for line in qres:
            line =  line.split(" ")
            line = [float(qres_rmsd) for qres_rmsd in line]
            data1.append(line)
        data1 = np.array(data1).astype(float)

        rmsd = open(filename_rmsd, "r").readlines()
        data2 = []
        for line in rmsd:
            line =  line.split(" ")
            line = [float(qres_rmsd) for qres_rmsd in line]
            data2.append(line)
        data2 = np.array(data2).astype(float)
        frame = data2 [:, 0]
        normframes = (frame-np.min(frame))/(np.max(frame)-np.min(frame)).reshape(-1, 1)
        frame_qres_rmsd = np.hstack((normframes.T, data1[:, 1:], data2[:, 1:]))       
        return frame_qres_rmsd
    
    def hierarchiacalClustering (self):
        outpath = self.outpath
        frame_qres_rmsd = self.getData()
        qres_rmsd = frame_qres_rmsd[:, 1:]
        fig, ax = plt.subplots(figsize = (3.5,2), constrained_layout = True)
        linkage_data = linkage(qres_rmsd, method='ward', metric='euclidean')
        dendrogram(linkage_data, ax=ax, color_threshold=5)
        ax.set_xticks([])
        fig.savefig(outpath + '/1_dandogram.png', format='png',dpi = 600)
        
    def elbow_method (self):
        """
        Elbow Method 
        Finding the ideal number of groups to divide the data into is a basic stage
        in any unsupervised algorithm. One of the most common techniques for figuring 
        out this ideal value of k is the elbow approach.
        """
        outpath = self.outpath
        frame_qres_rmsd = self.getData()
        qres_rmsd = frame_qres_rmsd[:, 1:]
        
        #Find optimum number of cluster
        sse = [] #SUM OF SQUARED ERROR
        for k in range(1,11):
            km = KMeans(n_clusters=k, random_state=2)
            km.fit(qres_rmsd)
            sse.append(km.inertia_)
            
        fig, ax = plt.subplots(1,1, figsize = (3.5,2), constrained_layout = True)
        ax.plot(sse)
        ax.set(xlabel ="Number of cluster (k)", 
              ylabel = "Sum Squared Error")

        fig.savefig(outpath + '/2_elbow.png', format='png', dpi = 600)
        
    def kmeansClustering (self):
        num_clusters = self.num_clusters
        frame_qres_rmsd = self.getData()
        qres_rmsd = frame_qres_rmsd[:, 1:]
        
        ### Clustering process
        kmeans = KMeans(n_clusters = num_clusters, random_state = 2)
        km = kmeans.fit(qres_rmsd)
        closest, _ = pairwise_distances_argmin_min(km.cluster_centers_, qres_rmsd)
        pred = kmeans.fit_predict(qres_rmsd)
        
        ##### finding closet data points from centroids
        tf_matrix = np.array(qres_rmsd)
        m_clusters = km.labels_.tolist()   #### clusters tag
        centers = np.array(km.cluster_centers_)
        
        centroid = {} ### data point near centers of clusters
        for i in range(len(centers)):
            center_vec = centers[i]
            data_idx_within_i_cluster = [ idx for idx, clu_num in enumerate(m_clusters) if clu_num == i ]
        
            one_cluster_tf_matrix = np.zeros((len(data_idx_within_i_cluster) , centers.shape[1] ) )
            for row_num, data_idx in enumerate(data_idx_within_i_cluster):
                one_row = tf_matrix[data_idx]
                one_cluster_tf_matrix[row_num] = one_row
        
            closest, _ = pairwise_distances_argmin_min(center_vec.reshape(1, -1), one_cluster_tf_matrix)
            closest_idx_in_one_cluster_tf_matrix = closest[0]
            closest_data_row_num = data_idx_within_i_cluster[closest_idx_in_one_cluster_tf_matrix]
            data_id = qres_rmsd[closest_data_row_num]
            centroid[closest_data_row_num] = data_id  
        assert len(centroid) == num_clusters
        
        return( [pred, centroid, kmeans] )
    
    def plotTsne2DN3Dscatter(self):
        ## applied tSNE on multidimentional data and centroids to represent them in 2D and 3D
        outpath = self.outpath
        [pred, centroid, kmeans] = self.kmeansClustering()
        cntarray = np.array(list(centroid.values()))
        frame_qres_rmsd = self.getData()
        qres_rmsd = frame_qres_rmsd[:, 1:]
        qres_rmsd_cnteroids  = np.vstack((qres_rmsd,cntarray))
        num_aa = int(len(qres_rmsd[0])/2)
        tsne = TSNE(n_components=2, random_state=0, learning_rate=2, n_iter=5000, n_iter_without_progress=400, verbose=1)
        tsne2D = tsne.fit_transform(qres_rmsd_cnteroids)
        centtsne2D = tsne2D[-len(cntarray):]
        tsne2D = tsne2D[:-len(cntarray)]
       
        ### 2D Ploting of clusters
        s = 6
        fig = plt.figure(figsize = (3.5, 2.5), constrained_layout = True)
        ax = plt.axes()
        ax.scatter(tsne2D[:,0],tsne2D[:,1],s = s, linewidth=0.5, c = pred, cmap=cm.tab20)
        xyc = []
        for f,c in centroid.items():
            center_mean_qres = np.mean(c[:num_aa])
            center_mean_rmsd = np.mean(c[num_aa:])
            c = pred[f]
            xyc.append([center_mean_qres,center_mean_rmsd,c, f])
        xyc = np.array(xyc)
        ### plotting centroids only
        scat = ax.scatter(centtsne2D[:,0], centtsne2D[:,1], marker = '*', edgecolors = 'r' , linewidth=0.5, s = s +40, c = xyc[:,2], cmap=cm.tab20)
        handles, labels = scat.legend_elements(prop="colors", num = 20)
        labels = xyc[:,3].astype(int)
        ax.legend(handles, labels, loc="upper left", ncol = 5, fontsize = 5, bbox_to_anchor=(0.1, 1.3))
        ax.set(ylabel ="Dimension-2", 
              xlabel = "Dimension-1")
        fig.savefig(outpath + '/3_scattered_cluster.png', format='png', dpi = 600)
        
        ### 3D Plotting of Clusters
        tsne = TSNE(n_components=3, random_state=0, learning_rate=2, n_iter=5000, n_iter_without_progress=400, verbose=1)
        tsne3D = tsne.fit_transform(qres_rmsd_cnteroids)
        centtsne2D = tsne3D[-len(cntarray):]
        tsne3D = tsne3D[:-len(cntarray)]
        fig = plt.figure(figsize = (3.5, 2.5))
        ax = plt.axes(projection ="3d")
        #ax.view_init(30, 70, 20)
        ax.scatter3D(tsne3D[:,0], tsne3D[:,1], tsne3D[:,2],s = s,c = pred, cmap=cm.tab20)
        scat = ax.scatter3D(centtsne2D[:,0],centtsne2D[:,1],centtsne2D[:,2] ,s = s + 40, marker = '*', edgecolors = 'r', linewidth=1, c = xyc[:,2], cmap=cm.tab20)
        handles, labels = scat.legend_elements(prop="colors", num = 20)
        labels = xyc[:,3].astype(int)
        ax.legend(handles, labels, loc="upper left", ncol = 5, fontsize = 5)
        ax.set_xlabel("Dimension-1", fontsize = 8)
        ax.set_ylabel("Dimension-2", fontsize = 8)
        ax.set_zlabel("Dimension-3", fontsize = 8)
        #ax.set_ylim(0,50)
        fig.savefig(outpath + '/3_3Dscattered_cluster.png', format='png', dpi = 600, bbox_inches = 'tight')
        return handles

        
        
    def plot2DN3Dscatter(self):
        outpath = self.outpath
        [pred, centroid, kmeans] = self.kmeansClustering()
        frame_qres_rmsd = self.getData()
        qres_rmsd = frame_qres_rmsd[:, 1:]
        frame = frame_qres_rmsd[:, 0]
        num_aa = int(len(qres_rmsd[0])/2)
        qres_avg = np.mean(qres_rmsd[:, :num_aa], axis= 1)
        rmsd_avg = np.mean(qres_rmsd[:, num_aa:], axis= 1)
        
        ### 2D Ploting of clusters
        s = 6
        fig = plt.figure(figsize = (3.5, 2.5), constrained_layout = True)
        ax = plt.axes()
        ### plotting all data including cenroids
        ax.scatter(qres_avg,rmsd_avg,s = s, linewidth=0.5, c = pred, cmap=cm.tab20)
        xyc = []
        for f,c in centroid.items():
            center_mean_qres = np.mean(c[:num_aa])
            center_mean_rmsd = np.mean(c[num_aa:])
            c = pred[f]
            xyc.append([center_mean_qres,center_mean_rmsd,c, f])
        xyc = np.array(xyc)
        ### plotting centroids only
        scat = ax.scatter(xyc[:,0], xyc[:,1], marker = '*', edgecolors = 'r' , linewidth=0.5, s = s +40, c = xyc[:,2], cmap=cm.tab20)
        handles, labels = scat.legend_elements(prop="colors", alpha=0.6)
        labels = xyc[:,3].astype(int)
        legend2 = ax.legend(handles, labels, loc="upper left", ncol = 3, fontsize = 6)
        ax.set(ylabel ="RMSD ($A$)", 
              xlabel = "$Q_{avg}$")
        fig.savefig(outpath + '/3_scattered_cluster.png', format='png', dpi = 600)
        
        ### 3D Plotting of Clusters
        fig = plt.figure(figsize = (3.5, 2.5), constrained_layout = True)
        ax = plt.axes(projection ="3d")
        #ax.view_init(30, 70, 20)
        # Creating plot
        ax.scatter3D(frame, qres_avg, rmsd_avg,s = s, c = pred, cmap=cm.Accent)
        for center in kmeans.cluster_centers_:
            center = np.array(center)
            center_mean_qres = np.mean(center[:num_aa])
            center_mean_rmsd = np.mean(center[num_aa:])
            ax.scatter3D(center[0],center_mean_qres,center_mean_rmsd ,s = s +2,marker = '*')
    
        ax.set(xlabel ="Time (ns)", 
              ylabel = "$Q_{avg}$", zlabel = "RMSD (A)")
        
        fig.savefig(outpath + '/3_3Dscattered_cluster.png', format='png', dpi = 600, bbox_inches = 'tight')

    
    def moveCentroidsPdbs(self):
        outpath = self.outpath
        filepath = self.filepath
        [pred, centroid, kmeans] = self.kmeansClustering()
        ###### Copy pdb of centroid list in corresponding directories
        outpath = qres_rmsd_instance.outpath
        dock = outpath + "/dock"
        if not os.path.exists(dock):
            os.makedirs(dock)

        centroid_pdbs = [str(x) + ".pdb" for x in list(centroid.keys())]
        for file_name in centroid_pdbs:
            full_file_name = os.path.join(filepath, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, dock)
    


#%%
filepath = "/media/anupkumar/Backup Plus/project/ML_project/pytorch_AMP_BiGAN/W_BiGAN_GP_monash/model_saved_231206/simulation_model2/mono_simulations/analysis/mdl2_P18/centroids/"
num_clusters = 30

qres_rmsd_instance= centroidQresRMSD(filepath, num_clusters)
qres_rmsd = qres_rmsd_instance.getData()
qres_rmsd_instance.hierarchiacalClustering()
qres_rmsd_instance.elbow_method()
qres_rmsd_instance.moveCentroidsPdbs()
test = qres_rmsd_instance.plotTsne2DN3Dscatter()
[pred, centroid, kmeans] = qres_rmsd_instance.kmeansClustering()

