from Visualization import visualization
import numpy as np
import pymesh


R = 1
center = np.array([0,0,0]);
init_mesh = pymesh.generate_icosphere(R,center,2);
v_new = np.load('vertices.npy')
fig = visualization(v_new,init_mesh.vertices,init_mesh.faces);
fig.plotPoly(2)
    