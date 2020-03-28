#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 16:23:12 2020

@author: cunchengzhu
"""
import numpy as np
import pymesh
from Optimizer import optimizer
from Visualization import visualization
#import trimesh;


center = np.array([0,0,0]);
mesh = pymesh.generate_icosphere(1.0,center,1);
vertices = mesh.vertices;

shape_vertices = np.shape(mesh.vertices);
#vertices = np.reshape(x,shape_vertices,order='F');
#mesh = pymesh.form_mesh(vertices, self.mesh.faces);
mesh.add_attribute("face_normal")
face_normal = mesh.get_attribute("face_normal")
face_normal = np.reshape(face_normal,np.shape(mesh.faces))
mesh.add_attribute("face_area")
At = mesh.get_attribute("face_area")
a = np.zeros(np.size(mesh.vertices))
print(face_normal)

