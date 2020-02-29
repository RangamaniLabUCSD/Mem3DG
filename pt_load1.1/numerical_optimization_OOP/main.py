
import numpy as np;
import pymesh;
from Optimizer import optimizer
from Visualization import visualization

center = np.array([0,0,0]);
mesh = pymesh.generate_icosphere(1.0,center,1);
v = mesh.vertices;
n = 34;


force = np.zeros(np.size(v));
force[3*n:3*n+3] = v[n]*(3);

opt = optimizer(10,force,mesh,1e-15);
v_new = np.reshape(opt.minimize_newton(),np.shape(v),order='F');
print(v_new-v);


fig = visualization(v_new,v,mesh.faces);
fig.plotPoly(n)
