#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 14:23:17 2019

@author: cunchengzhu
"""
import matplotlib.pyplot as plt
#import ipyvolume as ipv
from mpl_toolkits.mplot3d import Axes3D  
import numpy as np
import pymesh 
from matplotlib import cm

# ipv.figure()
# mesh = ipv.plot_trisurf(v[:,0],v[:,1],v[:,2],triangles = f, color = 'orange');
# ipv.scatter(v[:,0],v[:,1], v[:,2], marker='sphere', color='blue');
# ipv.xyzlim(-2, 2);
# ipv.show()

# fig, axs = plt.subplots(2)
# fig.suptitle('Vertically stacked subplots')
# axs[0].plot(x, y)
# axs[1].plot(x, -y)

# fig = plt.figure()
# ax = plt.axes(projection='3d')

class visualization(object):
    def __init__(self,v1,v2,f):
        self.v1 = v1;
        self.v2 = v2;
        self.f = f;
        
    def H_on_face(self,vertices,faces):
        mesh = pymesh.form_mesh(vertices, faces)
        mesh.add_attribute("vertex_mean_curvature")
        H = mesh.get_attribute("vertex_mean_curvature")
        H_on_face = np.zeros(mesh.faces.shape[0])
        #print(H_on_face)
        for i in range(mesh.faces.shape[0]):
            H_on_face[i] = np.mean([H[mesh.faces[i]]])
            #print([mesh.faces[i]])
        #print(H_on_face)
        
        # norm = plt.Normalize()
        # colors = plt.cm.jet(norm(H_on_face))
            
        scamap = plt.cm.ScalarMappable(cmap='inferno')
        colors = scamap.to_rgba(H_on_face)
        return colors 

    def plotPoly(self):
        
        # fig,(ax1, ax2) = plt.subplots(1, 2)
        
        # ax1 = plt.axes(projection='3d')

 
        
        # ax1.scatter3D(self.v1[n][0],self.v1[n][1],self.v1[n][2]);
    
        # #ax2 = plt.axes(projection='3d')
        
        # ax2.plot_trisurf(self.v2[:,0], self.v2[:,1], self.v2[:,2], triangles=self.f,
        #                 cmap='viridis', linewidths=0.2);
    
        # ax2.set_xlim(-1, 1); 
        # ax2.set_ylim(-1, 1); 
        # ax2.set_zlim(-1, 1);
        
        # ax2.scatter3D(self.v2[n][0],self.v2[n][1],self.v2[n][2]);
        
        
#         scamap = plt.cm.ScalarMappable(cmap='inferno')
# fcolors = scamap.to_rgba(C)
# ax.plot_surface(X, Y, Z, facecolors=fcolors, cmap='inferno')
# fig.colorbar(scamap)
# plt.show()
        
        
        
        color = self.H_on_face(self.v1,self.f)
        fig = plt.figure(figsize = plt.figaspect(1))
        Axes3D(fig);
        ax1 = plt.axes(projection='3d')
        xm = np.mean(self.v1[:,0])
        ym = np.mean(self.v1[:,1])
        zm = np.mean(self.v1[:,2])
        ax1.plot_trisurf(self.v1[:,0]-xm, self.v1[:,1]-ym, self.v1[:,2]-zm, triangles=self.f,cmap='viridis',linewidths=0.2);
        #ax1.scatter3D(self.v1[n][0],self.v1[n][1],self.v1[n][2],s = 100,marker = 'o');
        
     
        
        ax1.set_xlim(-1, 1); 
        ax1.set_ylim(-1, 1); 
        ax1.set_zlim(-1, 1);
        
        ax1.set_title('V = 0.85');
        # ax2 = fig.add_subplot(1, 2, 2,projection='3d')
        # ax2.plot_trisurf(self.v2[:,0], self.v2[:,1], self.v2[:,2], triangles=self.f,
        #                 cmap='viridis', linewidths=0.2);
        # #ax2.scatter3D(self.v2[n][0],self.v2[n][1],self.v2[n][2],s = 100, marker='o');
        # ax2.set_xlim(-1, 1); 
        # ax2.set_ylim(-1, 1); 
        # ax2.set_zlim(-1, 1);
        
        plt.show()

