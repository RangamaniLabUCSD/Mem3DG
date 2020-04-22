

#include <iostream>
#include <math.h>

#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include "geometrycentral/numerical/linear_solvers.h"

#include <Eigen/Core>

#include "ddgsolver/force.h"


namespace gc  = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

Eigen::Matrix<double, Eigen::Dynamic, 3> Force::bending_force(double Kb, double H0) {

    //vpg.requireCotanLaplacian();
    //vpg.requireVertexGalerkinMassMatrix();
    vpg.requireVertexGaussianCurvatures();

    // weak(conformal) laplacian operator
    //Eigen::SparseMatrix<double>& L = vpg.cotanLaplacian;
    //std::cout << Eigen::MatrixXd(L) << std::endl;

    // Gaussian curvature
    // gcs::VertexData<double>& gaussian = vpg.vertexGaussianCurvatures;
    // Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::AlignedMax>
    auto KG = vpg.vertexGaussianCurvatures.toMappedVector();
    //std::cout << "Gaussian" << KG << std::endl;

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
        // std::cout << "curvature" << H(row) << std::endl;
        n.row(row) = Hn.row(row) / H(row);
    }

    Eigen::Matrix<double, Eigen::Dynamic, 1> lap_H = M_inv * L * H;
    auto f_mag = M * (-2 * Kb * (2 * (H.array() - Eigen::ArrayXd::Constant(n_vertices, 1, H0))
        * (square(H.array()) + H0 * H.array() - (M_inv * KG).array()) + lap_H.array()).matrix());
    //std::cout << "force Magnitude" << f_mag << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, 3> f;
    f.resize(n_vertices, 3);

    for (size_t row = 0; row < mesh.nVertices(); ++row) {
        f.row(row) = f_mag(row) * n.row(row);
    }

    vpg.unrequireCotanLaplacian();
    vpg.unrequireVertexGaussianCurvatures();
    return f;
}
