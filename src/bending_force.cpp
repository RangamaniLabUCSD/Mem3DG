
#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include "ddgsolver/force.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getBendingForces(double &Kb, double &H0) {

  //// update L and M_inv if necessary 

  //// Initialize the mass matrix
  //M = vpg.vertexLumpedMassMatrix;

  //// Initialize the inverted Mass matrix
  //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  //solver.compute(M);
  //Eigen::SparseMatrix<double> I(mesh.nVertices(), mesh.nVertices());
  //I.setIdentity();
  //M_inv = solver.solve(I);

  //// Initialize the conformal Laplacian matrix
  // L = vpg.cotanLaplacian;

  bendingForces.fill({ 0,0,0 });

  // Gaussian curvature per vertex Area
  Eigen::Matrix<double, Eigen::Dynamic, 1> KG = M_inv * (vpg.vertexGaussianCurvatures.toMappedVector());

  // number of vertices for convenience 
  size_t n_vertices = (mesh.nVertices());

  // map ivp to eigen matrix position 
  auto positions = ddgsolver::EigenMap<double, 3>(vpg.inputVertexPositions);

  // map the VertexData bendingForces to eigen matrix bendingForces_e
  auto bendingForces_e = ddgsolver::EigenMap<double, 3>(bendingForces);
  
  // the build-in angle weight vertex normal
  auto vertexAngleNormal_e = ddgsolver::EigenMap<double, 3>(vpg.vertexNormals);

  // calcuate mean curvature per vertex area by laplacian matrix
  Hn = M_inv * L * positions / 2;

  //// project the mean curvature to angle weight vertex normal
  //for (size_t row = 0; row < mesh.nVertices(); ++row) {
  //  Hn.row(row) = Hn.row(row).dot(vertexNormal_e.row(row)) * vertexNormal_e.row(row);
  //}

  // initialize area gradient vertex normal 
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexAreaGradientNormal_e;
  vertexAreaGradientNormal_e.resize(n_vertices, 3);

  // calculate the area gradient vertex normal
  for (size_t row = 0; row < n_vertices; ++row) {
    vertexAreaGradientNormal_e.row(row) = Hn.row(row).normalized();
    vertexAreaGradientNormal_e.row(row) *= 
      copysign(1.0, vertexAreaGradientNormal_e.row(row).dot(vertexAngleNormal_e.row(row)));
  }

  // initialize the spontaneous curvature matrix 
  H0n = H0 * vertexAreaGradientNormal_e;

  // calculate laplacian H 
  Eigen::Matrix<double, Eigen::Dynamic, 3> lap_H = M_inv * L * Hn;

  // initialize and calculate intermediary result productTerms
  Eigen::Matrix<double, Eigen::Dynamic, 3> productTerms;
  productTerms.resize(n_vertices, 3);

  for (size_t row = 0; row < mesh.nVertices(); ++row) {
    if (Hn.row(row).dot(Hn.row(row)) + H0n.row(row).dot(Hn.row(row)) - KG(row) == 0) {
      Eigen::Matrix<double, 1, 3> zeros;
      zeros << 0.0, 0.0, 0.0;
      productTerms.row(row) = zeros;
      //std::cout << "smaller than 0 !!" << std::endl;
    }
    else {
      productTerms.row(row) = 2 * (Hn.row(row) - H0n.row(row))
        * (Hn.row(row).dot(Hn.row(row)) + H0n.row(row).dot(Hn.row(row)) - KG(row));
    }
  }

  // calculate bendingForce 
  bendingForces_e = M * (-2 * Kb * (productTerms + lap_H));

}
} // end namespace ddgsolver
