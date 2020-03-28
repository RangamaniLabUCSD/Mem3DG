
import numpy as np
import pymesh
from Optimizer import optimizer
from Visualization import visualization
#import trimesh;


center = np.array([0,0,0]);
mesh = pymesh.generate_icosphere(1.0,center,3);
v = mesh.vertices;
# mesh.enable_connectivity();
# adj_faces = mesh.get_vertex_adjacent_faces(1)
# adj_vertices = mesh.get_vertex_adjacent_vertices(1)
# mesh.add_attribute("vertex_normal")
# vertex_normal = mesh.get_vertex_attribute("vertex_normal")
# face_normal = mesh.get_attribute("face_normal")
# n = 34;
# tri_mesh = trimesh.Trimesh(vertices=v,
#                        faces=mesh.faces);
# print(tri_mesh.volume)


#for i in range(5):
pressure = -1
opt = optimizer(500,pressure,mesh,1e-15)
v_new = np.reshape(opt.minimize_newton(),np.shape(v),order='F');
fig = visualization(v_new,v,mesh.faces);
fig.plotPoly(2)
mesh = pymesh.form_mesh(v_new, mesh.faces)
    #print(v_new-v);



