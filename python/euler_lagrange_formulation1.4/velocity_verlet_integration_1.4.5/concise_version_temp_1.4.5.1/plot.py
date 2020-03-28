import numpy as np
import pymesh
#from scipy.sparse.linalg import inv as sp_inv
from scipy.sparse.linalg import spsolve
#import scipy
#import matplotlib.pyplot as plt
from Visualization import visualization



R = 1
center = np.array([0,0,0]);
init_mesh = pymesh.generate_icosphere(R,center,2);
v_new = np.load('vertices.npy')

mesh = pymesh.form_mesh(v_new, init_mesh.faces)
assembler = pymesh.Assembler(mesh)
M = assembler.assemble("mass")
L = assembler.assemble("laplacian")

Hn = spsolve(M,L @ v_new) 
H = np.linalg.norm(Hn,axis = 1)

fig = visualization(v_new,init_mesh.vertices,init_mesh.faces);
fig.plotPoly(2)
    