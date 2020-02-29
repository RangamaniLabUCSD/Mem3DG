#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 14:23:17 2019

@author: cunchengzhu
"""
import matplotlib.pyplot as plt
import ipyvolume as ipv
from mpl_toolkits.mplot3d import Axes3D  

# ipv.figure()
# mesh = ipv.plot_trisurf(v[:,0],v[:,1],v[:,2],triangles = f, color = 'orange');
# ipv.scatter(v[:,0],v[:,1], v[:,2], marker='sphere', color='blue');
# ipv.xyzlim(-2, 2);
# ipv.show()


def plotPoly(v,f):
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(v[:,0], v[:,1], v[:,2], triangles=f,
                    cmap='viridis', linewidths=0.2);

    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1);