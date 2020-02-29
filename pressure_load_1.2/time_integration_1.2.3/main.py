import numpy as np
import pymesh
from Integrator import integrator

center = np.array([0,0,0]);
mesh_perfect = pymesh.generate_icosphere(1.0,center,3);
mesh = pymesh.form_mesh(mesh_perfect.vertices + np.random.uniform(-0.01,0.01,np.shape(mesh_perfect.vertices)),mesh_perfect.faces);
v = mesh.vertices;
#x = np.reshape(v,np.size(v),order='F');
#np.shape(v)
pressure = -5
opt = integrator(50000,pressure,mesh)
opt.verlet_int(10,0.00001)