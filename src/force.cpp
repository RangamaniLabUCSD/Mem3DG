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
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

Eigen::Matrix<double, Eigen::Dynamic, 3> bending_force(std::unique_ptr<HalfedgeMesh>& mesh,
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry>& vpg, double Kb, double H0) {


    gcs::IntrinsicGeometryInterface& geometry = *vpg;

    geometry.requireCotanLaplacian();
    geometry.requireVertexGalerkinMassMatrix();
    geometry.requireVertexGaussianCurvatures();

	// weak(conformal) laplacian operator
	Eigen::SparseMatrix<double>& L = geometry.cotanLaplacian;
	//std::cout << Eigen::MatrixXd(L) << std::endl;

	// Mass matrix (Galerkin approximation)
	Eigen::SparseMatrix<double>& M = geometry.vertexGalerkinMassMatrix;
	//std::cout << "mass" << Eigen::MatrixXd(M).inverse() << std::endl;

	// Gaussian curvature
	gcs::VertexData<double>& gaussian = geometry.vertexGaussianCurvatures;
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

    size_t n_vertices = (mesh->nVertices());
    gc::Vector3* d = vpg->inputVertexPositions.rawdata().data();
    Eigen::Map<Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1>> evec3(d, n_vertices);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3>> evecdouble(reinterpret_cast<double*>(d), n_vertices, 3);
    //std::cout << MatrixXd(M).inverse() << std::endl;
    auto Hn = MatrixXd(M).inverse() * MatrixXd(L) * evecdouble/2;
    Matrix<double, Dynamic, 1> H;                     
    Matrix<double, Dynamic, 3> n;
    H.resize(n_vertices, 1);
    n.resize(n_vertices, 3);


    for (int row = 0; row < n_vertices; ++row) {
        H(row) = Hn.row(row).norm();
        n.row(row) = Hn.row(row) / H(row);
    }


    //std::cout << n << std::endl;

    auto lap_H = MatrixXd(M).inverse() * MatrixXd(L) * H;
    std::cout << (H - MatrixXd::Constant(n_vertices, 1, H0)).array() * square(H.array()) << std::endl;
    auto f_mag = M * (-2 * Kb * (2 * ((H - MatrixXd::Constant(162,1,H0)).array() * square(H.array())).matrix() + H0 * H - KG) + lap_H);
    Matrix<double, Dynamic, 3> f;
    f.resize(n_vertices, 3);

    for (size_t row = 0; row < mesh->nVertices(); ++row) {
        f.row(row) = f_mag(row) * n.row(row);
    }

    return n;

}

