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


    def plotPoly(self,n):
        
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

        
        fig = plt.figure(figsize = plt.figaspect(0.5))
        #Axes3D(fig);
        ax1 = fig.add_subplot(1, 2, 1,projection='3d')
        ax1.plot_trisurf(self.v1[:,0], self.v1[:,1], self.v1[:,2], triangles=self.f,
                        cmap='viridis', linewidths=0.2);
        ax1.scatter3D(self.v1[n][0],self.v1[n][1],self.v1[n][2],s = 100,marker = 'o');
        
        v_new = np.load('vertices.npy')
        xm = np.mean(v_new[:,0])
        ym = np.mean(v_new[:,1])
        zm = np.mean(v_new[:,2])
        ax1.set_xlim(xm-1, xm+1)
        ax1.set_ylim(ym-1, ym+1)
        ax1.set_zlim(zm-1, zm+1)
        
        ax2 = fig.add_subplot(1, 2, 2,projection='3d')
        ax2.plot_trisurf(self.v2[:,0], self.v2[:,1], self.v2[:,2], triangles=self.f,
                        cmap='viridis', linewidths=0.2);
        ax2.scatter3D(self.v2[n][0],self.v2[n][1],self.v2[n][2],s = 100, marker='o');
        ax2.set_xlim(-1, 1); 
        ax2.set_ylim(-1, 1); 
        ax2.set_zlim(-1, 1);
        
        plt.show()
        

