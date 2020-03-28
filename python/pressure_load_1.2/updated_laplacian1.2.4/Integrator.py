import numpy as np
import pymesh
from Problem import problem 
from Visualization import visualization

class integrator(object):
    
    def __init__(self,penalty,ex_force,mesh):
        self.prob = problem(penalty,ex_force,mesh) 
        self.mesh = mesh

    def verlet_int(self,T,h):
        x_past = np.reshape(self.mesh.vertices,np.size(self.mesh.vertices),order='F') 
        x = np.reshape(self.mesh.vertices,np.size(self.mesh.vertices),order='F') + np.random.uniform(-0.1,0.1,np.size(self.mesh.vertices))
        for i in range(int(T/h)):    
            an = -self.prob.potential_gradient(x)
            x_next = 2*x - x_past + an*h*h
            
            v_new = np.reshape(x_next,np.shape(self.mesh.vertices),order='F');
            fig = visualization(v_new,self.mesh.vertices,self.mesh.faces);
            fig.plotPoly(2)
            
            x_temp = x
            x = x_next
            x_past = x_temp
            
