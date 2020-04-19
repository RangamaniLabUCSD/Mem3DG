#include "ddgsolver/force.h"
#include "geometrycentral/utilities/vector3.h"
#include <math.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include "geometrycentral/utilities/vector3.h"
#include <iostream>
#include<Eigen/IterativeLinearSolvers>
#include "geometrycentral/numerical/linear_solvers.h"

using namespace Eigen;
using Vector3 = geometrycentral::Vector3;

Eigen::Matrix<double, Eigen::Dynamic, 3> bending_force(std::unique_ptr<HalfedgeMesh>& mesh,
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry>& vpg, double Kb, double H0) {

    namespace geosurf = geometrycentral::surface;
    geosurf::IntrinsicGeometryInterface& geometry = *vpg;

    geometry.requireCotanLaplacian();
    geometry.requireVertexGalerkinMassMatrix();
    geometry.requireVertexGaussianCurvatures();

	// weak(conformal) laplacian operator
	Eigen::SparseMatrix<double> L = geometry.cotanLaplacian;
	//std::cout << Eigen::MatrixXd(L) << std::endl;

	// Mass matrix (Galerkin approximation)
	Eigen::SparseMatrix<double> M = geometry.vertexGalerkinMassMatrix;
	//std::cout << "mass" << Eigen::MatrixXd(M).inverse() << std::endl;

	// Gaussian curvature
	geosurf::VertexData<double> gaussian = geometry.vertexGaussianCurvatures;
	Eigen::Matrix<double, Eigen::Dynamic, 1> KG = gaussian.toVector();
	//std::cout << "Gaussian" << KG << std::endl;

    //std::cout << "force cpp, KG size:  " << KG.size() << std::endl;
    std::vector<std::vector<std::size_t>> face_list = mesh->getFaceVertexList();
    geosurf::VertexData<Vector3> vertex_list = vpg->inputVertexPositions;
    static size_t n_vertices = (mesh->nVertices());
    //std::cout << "no vertices" << std::endl;

    Eigen::Matrix<double, 162, 3> vl;

    //std::cout << "1,2,: "<< vertex_list[5][2] << std::endl;

    for (int col = 0; col < 3; ++col) {
        for (int row = 0; row < vertex_list.size(); ++row) {
            //std::cout << "1,2,: " << vertex_list[5][2] << std::endl;
            vl(row,col) = vertex_list[row][col];
        }
    }
    //std::cout << MatrixXd(M).inverse() << std::endl;
    Matrix<double, 162, 3> Hn = MatrixXd(M).inverse() * MatrixXd(L) * vl/2;
    Matrix<double, 162, 1> H;
    for (int row = 0; row < n_vertices; ++row) {
        H(row) = Hn.row(row).norm();
    }
    Matrix<double, 162, 3> n;
    for (size_t row = 0; row < n_vertices; ++row) {
        n.row(row) = Hn.row(row) / H(row);
    }

    //std::cout << n << std::endl;

    Matrix<double, 162, 1> lap_H = MatrixXd(M).inverse() * MatrixXd(L) * H;
    //std::cout << (H - MatrixXd::Constant(162, 1, H0)).array() * square(H.array()) << std::endl;
    Matrix<double, 162, 1> f_mag = M * (-2 * Kb * (2 * ((H - MatrixXd::Constant(162,1,H0)).array() * square(H.array())).matrix() + H0 * H - KG) + lap_H);
    Matrix<double, 162, 3> f;
    for (size_t row = 0; row < n_vertices; ++row) {
        f.row(row) = f_mag(row) * n.row(row);
    }

    return f;

}

