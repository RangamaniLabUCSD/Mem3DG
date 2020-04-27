#pragma once

#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <Eigen/Core>
#include <Eigen/SparseLU>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;


class Force{
public:
    gcs::HalfedgeMesh& mesh;
    gcs::VertexPositionGeometry& vpg;
    
    double time_step;

    Eigen::Matrix<double, Eigen::Dynamic, 3> bf;
    Eigen::Matrix<double, Eigen::Dynamic, 3> sf;
    Eigen::Matrix<double, Eigen::Dynamic, 3> pf;
    Eigen::Matrix<double, Eigen::Dynamic, 3> df;
    Eigen::Matrix<double, Eigen::Dynamic, 3> xf;

    Eigen::SparseMatrix<double> M;
    Eigen::SparseMatrix<double> M_inv;

    Eigen::SparseMatrix<double> L;

    gcs::FaceData<double> face_area_init;
    double total_face_area_init = 0.0;

    double volume_init;
    
    gcs::VertexData<gc::Vector3> vertex_position_past;

    Force(gcs::HalfedgeMesh& mesh_, gcs::VertexPositionGeometry& vpg_, double time_step_): mesh(mesh_), vpg(vpg_), time_step(time_step_) {
        vpg.requireFaceNormals();
        vpg.requireVertexGalerkinMassMatrix();
        vpg.requireCotanLaplacian();
        vpg.requireFaceAreas();
        vpg.requireVertexIndices();
        vpg.requireVertexGaussianCurvatures();



        //initialize force matrix
        bf.setZero(mesh.nVertices(), 3);
        sf.setZero(mesh.nVertices(), 3);
        pf.setZero(mesh.nVertices(), 3);
        df.setZero(mesh.nVertices(), 3);
        xf.setZero(mesh.nVertices(), 3);

        // find the mass matrix 
        M = vpg.vertexGalerkinMassMatrix;
        // find the inverse of mass matrix 
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
        solver.compute(M);
        std::size_t n = mesh.nVertices();
        Eigen::SparseMatrix<double> I(n, n);
        I.setIdentity();
        M_inv = solver.solve(I);
        // find the confomal laplacian matrix 

        L = vpg.cotanLaplacian;
        // find the initial faceArea 
        // !!!!!!!!!!!
        // !!!!!!!!!!! this might raise problem?? if the underlying vpg got changes, will it be affected? 
        // Or not, since the face in rooted in intrinsic geometry, change the vertex positon will not affect. 
        // anyway, worth coming back and revise.
        face_area_init = vpg.faceAreas;
        for (gcs::Face f : mesh.faces()) {
            total_face_area_init += face_area_init[f];
        }   
        
        for (gcs::Face f : mesh.faces()) {
            double face_volume = signed_volume_from_face(f, vpg);
            volume_init += face_volume;
        }

        // initialize the vertex position of the last iteration
        vertex_position_past = vpg.inputVertexPositions;
    }

    ~Force(){
        vpg.unrequireFaceAreas();
        vpg.unrequireVertexGalerkinMassMatrix();
        vpg.unrequireCotanLaplacian();
        vpg.unrequireFaceNormals();
        vpg.unrequireVertexIndices();
        vpg.unrequireVertexGaussianCurvatures();
    }

    void update_cached_values();
    
    void bending_force(double Kb, double H0);
    void bending_force(double Kb, Eigen::Matrix<double, Eigen::Dynamic, 1> H0);
    void stretching_force(double Ksl, double Ksg);
    void pressure_force(double Kv, double Vt);
    void damping_force(double gamma);
    void stochastic_force(double sigma);

    double signed_volume_from_face(gcs::Face& f, gcs::VertexPositionGeometry& vpg);
    gc::Vector3 vec_from_halfedge(gcs::Halfedge& he, gcs::VertexPositionGeometry& vpg);

    void update_Vertex_positions() {
        vpg.refreshQuantities();
    }
};


