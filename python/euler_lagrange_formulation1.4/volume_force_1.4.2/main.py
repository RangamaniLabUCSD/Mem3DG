import numpy as np
import pymesh
#from scipy.sparse.linalg import inv as sp_inv
from scipy.sparse.linalg import spsolve
#import scipy
#import matplotlib.pyplot as plt
from Visualization import visualization
import trimesh

Kb = 0.01
R = 1
center = np.array([0,0,0]);
init_mesh = pymesh.generate_icosphere(R,center,2);
H0 = 0*np.ones(init_mesh.vertices.shape[0])
Ksl = 2
Ksg = 2
Kv = 1
gamma = 1
#P = 0*np.ones(init_mesh.vertices.shape[0])
vol0 = 0.6 * 4/3 * np.pi* R**3

def bending_force(init_mesh,x,H0,Kb):
    v = np.reshape(x,init_mesh.vertices.shape,order='F')
    mesh = pymesh.form_mesh(v, init_mesh.faces)
    assembler = pymesh.Assembler(mesh)
    M = assembler.assemble("mass")
    L = assembler.assemble("laplacian")
    
    Hn = spsolve(M,L @ v) 
    H = np.linalg.norm(Hn,axis = 1)
    n = Hn / H[:,None]
    
    mesh.add_attribute("vertex_gaussian_curvature")
    G = mesh.get_attribute("vertex_gaussian_curvature")
    
    lap_H =  spsolve(M,L @ H) 
    f_mag = M @ (-2 * Kb * (2 * (H-H0) * (H**2 + H0 * H - G) 
                   + lap_H))
    f = n * f_mag[:,None]
    f = np.reshape(f,np.size(v),order='F') 
    return f,H

def stretching_force(init_mesh,x,Ksl,Ksg):
    v = np.reshape(x,init_mesh.vertices.shape,order="F")
    mesh = pymesh.form_mesh(v, init_mesh.faces)
    
    mesh.add_attribute("face_normal")
    face_normal = mesh.get_attribute("face_normal")
    face_normal = np.reshape(face_normal,init_mesh.faces.shape,order='C')
    
    mesh.add_attribute("face_area")
    At = mesh.get_attribute("face_area")
    a = np.zeros(np.size(x))
    b = np.zeros(np.size(x))
    
    init_mesh.add_attribute("face_area")
    At0 = init_mesh.get_attribute("face_area")
    
    mesh.enable_connectivity()
    for i in range(init_mesh.vertices.shape[0]):
        
        adj_faces = mesh.get_vertex_adjacent_faces(i)
        ai = np.zeros(3)
        bi = np.zeros(3)
        for j in range(adj_faces.size):
            face = adj_faces[j]
            tri_v = np.setdiff1d(mesh.faces[face],i)
            edge = v[tri_v[0]] - v[tri_v[1]]
            gradient = np.cross(face_normal[face],edge)
            
            if np.dot(gradient,v[i] - v[tri_v[0]])<0:
                gradient = -gradient 
            
            ai = ai + 2*gradient * (At[face] -At0[face])/At0[face] 
            bi = bi + gradient
            
        a[[i,i+mesh.vertices.shape[0], i + 2*mesh.vertices.shape[0]]] = ai
        b[[i,i+mesh.vertices.shape[0], i + 2*mesh.vertices.shape[0]]] = bi \
            * 2 * (At.sum()-At0.sum())/At.sum()
    return -Ksl * a - Ksg * b


def pressure_force(init_mesh,x,Kv,cente,vol0):
    v = np.reshape(x,init_mesh.vertices.shape,order='F')
    mesh = pymesh.form_mesh(v, init_mesh.faces)
    mesh.add_attribute("face_normal")
    face_normal = mesh.get_attribute("face_normal")
    face_normal = np.reshape(face_normal,init_mesh.faces.shape,order='C')
    
    tri_mesh = trimesh.Trimesh(vertices=v,
                        faces=mesh.faces);
    vol = tri_mesh.volume
    
    a = np.zeros(x.size)
    mesh.enable_connectivity()
    for i in range(init_mesh.vertices.shape[0]):
        
        adj_faces = mesh.get_vertex_adjacent_faces(i)
        ai = np.zeros(3)
        for j in range(adj_faces.size):
            face = adj_faces[j]
            tri_v = np.setdiff1d(mesh.faces[face],i)
            edge1 = v[tri_v[0]]  - center 
            edge2 = v[tri_v[1]]  - center
            A_n = 0.5*np.cross(edge1,edge2)
            
            if np.dot(A_n,v[i] - v[tri_v[0]])<0:
                A_n = -A_n
                
            if np.dot(face_normal[face],edge1)<0:
                A_n = -A_n
              
            ai = ai + 1/3*A_n
            
        a[[i,i+mesh.vertices.shape[0], i + 2*mesh.vertices.shape[0]]] = ai
    a = 2*(vol - vol0)/vol0 *a
        
    return -Kv*a 
               


def damping_force(init_mesh,gamma,x_past,x,h):
    velo = np.reshape((x-x_past)/h,init_mesh.vertices.shape,order="F")
    
    v = np.reshape(x,init_mesh.vertices.shape,order="F")
    mesh = pymesh.form_mesh(v, init_mesh.faces)
    
    mesh.enable_connectivity()
    a = np.zeros(np.size(x))
    
    for i in range(init_mesh.vertices.shape[0]):
        adj_vertex = mesh.get_vertex_adjacent_vertices(i)
        ai = np.zeros(3)
        for j in range(adj_vertex.size):
            velo_diff = velo[i] - velo[adj_vertex[j]]
            posi_diff = (v[i] - v[adj_vertex[j]])\
                /np.linalg.norm((v[i] - v[adj_vertex[j]]))
            ai = ai + np.dot(velo_diff,posi_diff)*posi_diff
            
        a[[i,i+mesh.vertices.shape[0], i + 2*mesh.vertices.shape[0]]] = ai
    return -gamma*a

T = 500
h = 0.01

x_past = np.reshape(init_mesh.vertices,np.size(init_mesh.vertices),order='F') 
x = x_past+np.random.uniform(-0.001,0.001,np.size(init_mesh.vertices))

grad_norm_past = 0.0
grad_norm =0.0

for i in range(int(T/h)):    
    an1,H = bending_force(init_mesh,x,np.zeros(np.size(init_mesh.vertices.shape[0])),Kb)
    an2 = stretching_force(init_mesh,x,Ksl,Ksg)
    an3 = pressure_force(init_mesh,x,Kv,center,vol0)
    an4 = damping_force(init_mesh,gamma,x_past,x,h)
    print("fb: ",np.linalg.norm(an1),
          " fs: ",np.linalg.norm(an2),
          " fp: ",np.linalg.norm(an2),
          " fd: ",np.linalg.norm(an3))
    
    an = an1+an2 + an3 + an4
    
    grad_norm = (np.linalg.norm(an))
    print('i = ',i,'',grad_norm)
    
#   if abs(grad_norm-grad_norm_past)<1e-4:
#        gamma = 0.5
#        print('increase damping')
        
    if abs(grad_norm)<1e-4:
        print('converged')
        break 
    grad_norm_past = grad_norm 
    
    x_next = 2*x - x_past + an*h*h
    v_new = np.reshape(x_next,init_mesh.vertices.shape,order='F');
    fig = visualization(v_new,init_mesh.vertices,init_mesh.faces);
    fig.plotPoly(2)
    np.save('vertices',v_new)
    np.save('faces',init_mesh.faces)
    x_temp = x
    x = x_next
    x_past = x_temp
