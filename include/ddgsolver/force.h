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

    Eigen::SparseMatrix<double> M;
    Eigen::SparseMatrix<double> M_inv;

    Eigen::SparseMatrix<double> L;

    gcs::FaceData<double> face_area_init;

    Force(gcs::HalfedgeMesh& mesh_, gcs::VertexPositionGeometry& vpg_): mesh(mesh_), vpg(vpg_) {
        // find the mass matrix 
        vpg.requireVertexGalerkinMassMatrix();
        M = vpg.vertexGalerkinMassMatrix;
        // find the inverse of mass matrix 
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
        solver.compute(M);
        std::size_t n = mesh.nVertices();
        Eigen::SparseMatrix<double> I(n, n);
        I.setIdentity();
        M_inv = solver.solve(I);
        // find the confomal laplacian matrix 
        vpg.requireCotanLaplacian();
        L = vpg.cotanLaplacian;
        // find the initial faceArea 
        vpg.requireFaceAreas();
        gcs::FaceData<double> face_area_init = vpg.faceAreas;
    }

    ~Force(){vpg.unrequireVertexGalerkinMassMatrix();}

    Eigen::Matrix<double, Eigen::Dynamic, 3> bending_force(double Kb, double H0);

    Eigen::Matrix<double, Eigen::Dynamic, 3> bending_force(double Kb, Eigen::Matrix<double, Eigen::Dynamic, 1> H0);

    Eigen::Matrix<double, Eigen::Dynamic, 3> stretching_force(double Ksl, double Ksg);
};

/*
Eigen::Matrix<double, Eigen::Dynamic, 1> stretching_force(std::unique_ptr<HalfedgeMesh> mesh,
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> vpg);

Eigen::Matrix<double, Eigen::Dynamic, 1> pressure_force(std::unique_ptr<HalfedgeMesh> mesh,
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> vpg);
*/
