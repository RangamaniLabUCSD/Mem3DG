import pymesh;
import visualization
import geo_attribute
import numpy as np;
import minimization;

center = np.array([0,0,0]);
mesh = pymesh.generate_icosphere(1.0,center,1);
v = mesh.vertices;
f = mesh.faces;
mesh.add_attribute("vertex_mean_curvature");
Km = mesh.get_attribute("vertex_mean_curvature") ;

force = np.zeros(np.size(v));
force[1] = 3;
#print('force =',force[10]);
v_new = minimization.minimize_newton(mesh,10,force)
visualization.plotPoly(v_new,f);
#mesh = pymesh.form_mesh(v_new, f);
    
    
# assembler = pymesh.Assembler(mesh);
# L = assembler.assemble("laplacian");
# M = assembler.assemble("mass");
# curvature = geo_attribute.curvature(M, L, v);
# bending_energy = geo_attribute.bending_energy(curvature, M);
# surface_area = geo_attribute.surface_area(M);
    
