#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 14:27:02 2019

@author: cunchengzhu
"""

from numpy import linalg as LA
import numpy as np;
import scipy

def curvature(M,L,v):
    curvature = LA.norm(scipy.sparse.linalg.inv(M)*L*v, axis=1);
    return curvature 

def bending_energy(curvature,M):
    bending_energy = 1/2 * curvature.transpose() * M * curvature;
    return bending_energy

def surface_area(M):
    surface_area = np.sum(M)
    return surface_area
