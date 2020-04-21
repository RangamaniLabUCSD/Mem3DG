#include "ddgsolver/force.h"
#include "geometrycentral/utilities/vector3.h"
#include <math.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include "geometrycentral/utilities/vector3.h"
#include <iostream>
#include <Eigen/IterativeLinearSolvers>
#include "geometrycentral/numerical/linear_solvers.h"


namespace gc  = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

Eigen::Matrix<double, Eigen::Dynamic, 3> force::bending_force(double Kb, double H0) {

    //gcs::IntrinsicGeometryInterface& geometry = vpg;

    vpg.requireCotanLaplacian();
    //vpg.requireVertexGalerkinMassMatrix();
    vpg.requireVertexGaussianCurvatures();

    // weak(conformal) laplacian operator
    Eigen::SparseMatrix<double>& L = vpg.cotanLaplacian;
    //std::cout << Eigen::MatrixXd(L) << std::endl;

    // Mass matrix (Galerkin approximation)
    //Eigen::SparseMatrix<double>& M = vpg.vertexGalerkinMassMatrix;
    //std::cout << "mass" << Eigen::MatrixXd(M).inverse() << std::endl;

    // Gaussian curvature
    gcs::VertexData<double>& gaussian = vpg.vertexGaussianCurvatures;
    Eigen::Matrix<double, Eigen::Dynamic, 1> KG = gaussian.toVector();
    //std::cout << "Gaussian" << KG << std::endl;

    ////std::cout << "force cpp, KG size:  " << KG.size() << std::endl;
    //std::vector<std::vector<std::size_t>>& face_list = mesh->getFaceVertexList();
    //gcs::VertexData<Vector3> vertex_list = vpg->inputVertexPositions;
    //size_t n_vertices = (mesh->nVertices());
    ////std::cout << "no vertices" << std::endl;

    //Eigen::Matrix<double, Eigen::Dynamic, 3> vl;
    //vl.resize(n_vertices,3)
    ////std::cout << "1,2,: "<< vertex_list[5][2] << std::endl;

    //for (int col = 0; col < 3; ++col) {
    //    for (int row = 0; row < vertex_list.size(); ++row) {
    //        //std::cout << "1,2,: " << vertex_list[5][2] << std::endl;
    //        vl(row,col) = vertex_list[row][col];
    //    }
    //}

    size_t n_vertices = (mesh.nVertices());
    gc::Vector3* d = vpg.inputVertexPositions.rawdata().data();

    Eigen::Map<Eigen::Matrix<double, 3, Eigen::Dynamic>> positions(reinterpret_cast<double*>(d), 3, n_vertices);
    //std::cout << "evecdouble" << evecdouble.cols() << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, 3> Hn = M_inv * L * positions.transpose() / 2;
    Eigen::Matrix<double, Eigen::Dynamic, 1> H;
    Eigen::Matrix<double, Eigen::Dynamic, 3> n;
    H.resize(n_vertices, 1);
    n.resize(n_vertices, 3);


    for (int row = 0; row < n_vertices; ++row) {
        H(row) = Hn.row(row).norm();
        std::cout << "curvature" << H(row) << std::endl;
        n.row(row) = Hn.row(row) / H(row);
    }

    Eigen::Matrix<double, Eigen::Dynamic, 1> lap_H = M_inv * L * H;
    //std::cout << (H - MatrixXd::Constant(n_vertices, 1, H0)).array() * square(H.array()) << std::endl;
    auto f_mag = M * (-2 * Kb * (2 * ((H - Eigen::MatrixXd::Constant(n_vertices, 1, H0)).array() * square(H.array())).matrix() + H0 * H - KG) + lap_H);
    Eigen::Matrix<double, Eigen::Dynamic, 3> f;
    f.resize(n_vertices, 3);

    for (size_t row = 0; row < mesh.nVertices(); ++row) {
        f.row(row) = f_mag(row) * n.row(row);
    }

    return f;
}
