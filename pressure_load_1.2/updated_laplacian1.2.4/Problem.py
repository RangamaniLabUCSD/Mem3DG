import pymesh
import numpy as np
from scipy.sparse.linalg import inv as sp_inv
import scipy
import trimesh

class problem(object):
    
    def __init__(self,penalty,ex_force,mesh):
        
        self.mesh = mesh;
        assembler = pymesh.Assembler(mesh);
        self.init_M = assembler.assemble("lumped_mass");
        mesh.add_attribute("face_area")
        self.At0 = mesh.get_attribute("face_area")
        L = assembler.assemble("laplacian");
        hess_block = scipy.sparse.spmatrix.transpose(L)@sp_inv(self.init_M)@ L;
        self.hess = scipy.sparse.block_diag((hess_block,hess_block,hess_block));
        self.laplacian = sp_inv(self.init_M)@L
        self.penalty = penalty;
        #self.x0 = np.reshape(self.mesh.vertices,np.size(self.mesh.vertices),order='F');
        self.pressure = ex_force;
       
      
    ######################################################################
 
    def gradient_bending(self,x):
        
        shape_vertices = np.shape(self.mesh.vertices)
        vertices = np.reshape(x,shape_vertices,order='F')
        mesh = pymesh.form_mesh(vertices, self.mesh.faces)
        assembler = pymesh.Assembler(mesh)
        L = assembler.assemble("laplacian");
        M = assembler.assemble("lumped_mass");
        hess_block = scipy.sparse.spmatrix.transpose(L)@sp_inv(M)@L
        hess = scipy.sparse.block_diag((hess_block,hess_block,hess_block))
        
        return np.zeros(np.size(x)) #hess @ x
    
    def gradient_penalty(self,x):
        
        
        shape_vertices = np.shape(self.mesh.vertices);
        vertices = np.reshape(x,shape_vertices,order='F');
        mesh = pymesh.form_mesh(vertices, self.mesh.faces);
        mesh.add_attribute("face_normal")
        face_normal = mesh.get_attribute("face_normal")
        face_normal = np.reshape(face_normal,np.shape(mesh.faces))
        mesh.add_attribute("face_area")
        At = mesh.get_attribute("face_area")
        a = np.zeros(np.size(x))
        
        for v in range(shape_vertices[0]):
            mesh.enable_connectivity()
            adj_faces = mesh.get_face_adjacent_faces(v)
            av = 0
            for f in range(np.size(adj_faces)):
                tri_v = np.setdiff1d(mesh.faces[f],v)
                edge = vertices[tri_v[0]] - vertices[tri_v[1]]
                gradient = np.cross(face_normal[f],edge)
                
                if np.dot(gradient,vertices[v] - vertices[tri_v[0]])<0:
                    gradient = - gradient
                    
                av = av + 2*gradient * (At[f] - self.At0[f])/self.At0[f] 
            a[[v,v+shape_vertices[0],v+ 2 * shape_vertices[0]]] = av
        #a = scipy.optimize.approx_fprime(x, self.obj_func_penalty,self.nu_grad_size);
        #print(a)
        return self.penalty*a
    
    def force(self,x):
        #mesh = trimesh.Trimesh(vertices=self.mesh.vertices;,
        #               faces=self.mesh.faces);
        shape_vertices = np.shape(self.mesh.vertices);
        vertices = np.reshape(x,shape_vertices,order='F');
        mesh = pymesh.form_mesh(vertices, self.mesh.faces)
        mesh.add_attribute("vertex_normal")
        vector = mesh.get_vertex_attribute("vertex_normal")
        # assembler = pymesh.Assembler(mesh)
        # L = assembler.assemble("laplacian");
        # M = assembler.assemble("mass");
        # laplacian = sp_inv(M)@ L
        # vector = laplacian*vertices;
        # norm = np.linalg.norm(vector,2,1)
        # norm = norm[:,None]
        # vector = vector/norm;
        #print((np.diagonal(vertices@np.transpose(vector)))>0.1)
        vector = np.reshape(vector,np.size(self.mesh.vertices),order='F');
        force = self.pressure*vector
        return np.zeros(np.size(x)) #force
        
        
        
    def potential_gradient(self,x):
        a = self.gradient_bending(x)
        b = self.gradient_penalty(x) 
        c = self.force(x);
        print(np.linalg.norm(a),np.linalg.norm(b),np.linalg.norm(c))
        print(np.linalg.norm(a+b+c))
        return a+b+c
    
