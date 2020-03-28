import numpy as np
import pymesh
from Integrator import integrator

center = np.array([0,0,0]);
mesh_perfect = pymesh.generate_icosphere(1.0,center,1);
mesh = pymesh.form_mesh(mesh_perfect.vertices + np.random.uniform(-1e-1,1e-1,np.shape(mesh_perfect.vertices)),mesh_perfect.faces);
v = mesh.vertices;
#x = np.reshape(v,np.size(v),order='F');
#np.shape(v)
pressure = -5
opt = integrator(50000,pressure,mesh_perfect)
opt.verlet_int(0.0001,0.0000001)