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

    Force(gcs::HalfedgeMesh& mesh_, gcs::VertexPositionGeometry& vpg_): mesh(mesh_), vpg(vpg_) {
        vpg.requireVertexGalerkinMassMatrix();
        M = vpg.vertexGalerkinMassMatrix;

        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
        solver.compute(M);
        std::size_t n = mesh.nVertices();
        Eigen::SparseMatrix<double> I(n, n);
        I.setIdentity();
        M_inv = solver.solve(I);
    }

    ~Force(){vpg.unrequireVertexGalerkinMassMatrix();}

    Eigen::Matrix<double, Eigen::Dynamic, 3> bending_force(double Kb, double H0);

    Eigen::Matrix<double, Eigen::Dynamic, 3> bending_force(double Kb, Eigen::Matrix<double, Eigen::Dynamic, 1> H0);
};

/*
Eigen::Matrix<double, Eigen::Dynamic, 1> stretching_force(std::unique_ptr<HalfedgeMesh> mesh,
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> vpg);

Eigen::Matrix<double, Eigen::Dynamic, 1> pressure_force(std::unique_ptr<HalfedgeMesh> mesh,
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> vpg);
*/
