#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 16:19:59 2019

@author: cunchengzhu
"""

from scipy import optimize
import geo_attribute
import pymesh
import numpy as np;
import scipy

def objective_func_bfgs(x,eps,force,penalty,faces,shape_vertices,A0):
    
    vertices = np.reshape(x,shape_vertices,order='F');
    mesh = pymesh.form_mesh(vertices, faces);
    assembler = pymesh.Assembler(mesh);
    L = assembler.assemble("laplacian");
    M = assembler.assemble("lumped_mass");
    curvature = geo_attribute.curvature(M, L, vertices);
    bending_energy = geo_attribute.bending_energy(curvature, M);
    #surface_area = geo_attribute.surface_area(M);
    objective_func = sum(bending_energy) + penalty * scipy.sparse.linalg.norm(M-A0,'fro');
    
    return objective_func 

def gradient_numerical(x,eps,force,penalty,faces,shape_vertices,A0):
    grad_int = optimize.approx_fprime(x, objective_func_bfgs, eps, eps,force,penalty,faces,shape_vertices,A0);
    grad = grad_int + force;
    print(np.linalg.norm(grad));
    return grad


def minimize_bfgs(mesh,penalty,force):
    v = mesh.vertices;
    faces = mesh.faces;
    x0 = np.reshape(v,np.size(v),order='F');
    shape_vertices = np.shape(v);
    assembler = pymesh.Assembler(mesh);
    A0 = assembler.assemble("lumped_mass");
    eps = 1e-8;
    #x_new = optimize.minimize(objective_func, x0, args=(eps,force,penalty,faces,shape_vertices,A0), method='L-BFGS-B', jac=gradient, bounds=None, tol=None, callback=None, options={'disp': None, 'maxcor': 10, 'ftol': 2.220446049250313e-09, 'gtol': 0.1, 'eps': 1e-08, 'maxfun': 15000, 'maxiter': 15000, 'iprint': -1, 'maxls': 20})
    # x_new = optimize.minimize(objective_func, x0, args= (eps,force,penalty,faces,shape_vertices,A0), \
    #                         method='CG', jac=gradient, \
    #                             tol=None, callback=None, options={'gtol': 1e-05, 'norm': 2, 'eps': 1.4901161193847656e-08, 'maxiter': None, 'disp': False, 'return_all': False})
    x_new = optimize.fmin_bfgs(f= objective_func_bfgs,x0 = x0, fprime = gradient_numerical\
        ,args = (eps,force,penalty,faces,shape_vertices,A0)\
            ,gtol=1e-02, norm=2, epsilon=None, maxiter=None, full_output=1, disp=1, retall=1, callback=None);
    v_new = np.reshape(x_new.x,np.shape(v),order='F');
    return v_new

#############################################################################################################
    
def objective_func_newton(x,force,penalty,faces,hess,shape_vertices,A0):
    
    vertices = np.reshape(x,shape_vertices,order='F');
    mesh = pymesh.form_mesh(vertices, faces);
    #assembler = pymesh.Assembler(mesh);
#    L = assembler.assemble("laplacian");
    M = mesh.assembler.assemble("lumped_mass");
    # curvature = geo_attribute.curvature(M, L, vertices);
    # bending_energy = geo_attribute.bending_energy(curvature, M);
    #surface_area = geo_attribute.surface_area(M);
    # objective_func = sum(bending_energy) + penalty * scipy.sparse.linalg.norm(M-A0,'fro');
    # print(objective_func);
    objective_func = np.asscalar(np.matrix(x)*hess.todense()*np.matrix(x).T/2) + penalty * scipy.sparse.linalg.norm(M-A0,'fro');
    
    return objective_func 

def gradient(x,force,penalty,faces,hess,shape_vertices,A0):
    # grad = hess*np.reshape(x,shape_vertices,order='F');
    # grad = np.reshape(grad,np.size(x),order='F') + force;
    def penalty_energy(x,force,penalty,faces,hess,shape_vertices,A0):
        vertices = np.reshape(x,shape_vertices,order='F');
        mesh = pymesh.form_mesh(vertices, faces);
        assembler = pymesh.Assembler(mesh);
        M = assembler.assemble("lumped_mass");
        penalty_energy = penalty * scipy.sparse.linalg.norm(M-A0,'fro');
        return penalty_energy;
    eps=1e-9;     
    grad_penalty = optimize.approx_fprime(x, penalty_energy, eps, force,penalty,faces,hess,shape_vertices,A0);
    grad = hess*x+force+grad_penalty;
    print(np.linalg.norm(grad));
    return grad

def hessian(x,force,penalty,faces,hess,shape_vertices,A0):
    return hess.todense();

def minimize_newton(mesh,penalty,force):
    v = mesh.vertices;
    faces = mesh.faces;
    x0 = np.reshape(v,np.size(v),order='F');
    shape_vertices = np.shape(v);
    assembler = pymesh.Assembler(mesh);
    A0 = assembler.assemble("lumped_mass");
    L = assembler.assemble('laplacian');
    hess_block = scipy.sparse.spmatrix.transpose(L)\
        *scipy.sparse.linalg.inv(A0)* L;
    hess = scipy.sparse.block_diag((hess_block,hess_block,hess_block)); 
    #scipy.optimize.newton(objective_func_newton, x0, fprime=gradient,\
    #                      args=(force,penalty,faces,hess,shape_vertices,A0), tol=1.48e-08, maxiter=50, fprime2=hessian, x1=None, rtol=0.1, full_output=False, disp=False)
    #x_new = scipy.optimize.minimize(objective_func_newton, x0, args=(force,penalty,faces,hess,shape_vertices,A0), method='Newton-CG', jac=gradient, hess=None, hessp=None, tol=None, callback=None, options={'xtol': 1e-05, 'eps': 1.4901161193847656e-08, 'maxiter': None, 'disp': True, 'return_all': True});
    x_new = optimize.fmin_bfgs(f= objective_func_newton,x0 = x0, fprime = gradient\
        ,args=(force,penalty,faces,hess,shape_vertices,A0)\
            ,gtol=1e-02, norm=2, epsilon=None, maxiter=None, full_output=1, disp=1, retall=1, callback=None);
    v_new = np.reshape(x_new.x,np.shape(v),order='F');
    return v_new