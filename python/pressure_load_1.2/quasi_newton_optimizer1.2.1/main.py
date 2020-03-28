
import numpy as np;
import pymesh;
from Optimizer import optimizer
from Visualization import visualization
import trimesh;


center = np.array([0,0,0]);
mesh = pymesh.generate_icosphere(1.0,center,1);
v = mesh.vertices;
# n = 34;
# tri_mesh = trimesh.Trimesh(vertices=v,
#                        faces=mesh.faces);
# print(tri_mesh.volume)
for i in range(50):
    pressure = -0.5*i;
    opt = optimizer(50,pressure,mesh,1e-15);
    v_new = np.reshape(opt.minimize_newton(),np.shape(v),order='F');
    fig = visualization(v_new,v,mesh.faces);
    fig.plotPoly(2)
    mesh = pymesh.form_mesh(v_new, mesh.faces)
    #print(v_new-v);



