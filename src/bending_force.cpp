
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
#include "ddgsolver/meshops.h"

namespace ddgsolver {

  namespace gc = ::geometrycentral;
  namespace gcs = ::geometrycentral::surface;

void Force::getBendingForces() {
  // Gaussian curvature per vertex Area
  Eigen::Matrix<double, Eigen::Dynamic, 1> KG =
      M_inv * (vpg.vertexGaussianCurvatures.toMappedVector());

  // number of vertices for convenience
  std::size_t n_vertices = (mesh.nVertices());

  // map ivp to eigen matrix position
  auto positions = ddgsolver::EigenMap<double, 3>(vpg.inputVertexPositions);

  // map the VertexData bendingForces to eigen matrix bendingForces_e
  auto bendingForces_e = ddgsolver::EigenMap<double, 3>(bendingForces);
  bendingForces_e.setZero();

  // the build-in angle-weighted vertex normal
  auto vertexAngleNormal_e = ddgsolver::EigenMap<double, 3>(vpg.vertexNormals);

  // calculate mean curvature per vertex area and map it to angle-weighted normal
  Hn = rowwiseScaling(rowwiseDotProduct(vertexAngleNormal_e, 
            M_inv * L * positions / 2.0), vertexAngleNormal_e);
 
  // calculate the Laplacian of mean curvature H 
  Eigen::Matrix<double, Eigen::Dynamic, 3> lap_H = M_inv * L * Hn;

  // initialize the spontaneous curvature matrix
  H0n = H0 * vertexAngleNormal_e;

  // initialize and calculate intermediary result scalerTerms, set to zero if negative
  Eigen::Matrix<double, Eigen::Dynamic, 1> scalerTerms =
    rowwiseDotProduct(Hn, Hn) + rowwiseDotProduct(H0n, H0n) - KG;
  Eigen::Matrix<double, Eigen::Dynamic, 1> zeroMatrix;
  zeroMatrix.resize(n_vertices, 1);
  zeroMatrix.setZero();
  scalerTerms = scalerTerms.array().max(zeroMatrix.array());

  // initialize and calculate intermediary result productTerms
  Eigen::Matrix<double, Eigen::Dynamic, 3> productTerms;
  productTerms.resize(n_vertices, 3);
  productTerms = 2 * rowwiseScaling(scalerTerms, Hn - H0n);

  // calculate bendingForce
  bendingForces_e = M * (-2.0 * Kb * (productTerms + lap_H));
}
} // end namespace ddgsolver
