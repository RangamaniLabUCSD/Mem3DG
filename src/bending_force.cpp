
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

  // the build-in angle weight vertex normal
  auto vertexAngleNormal_e = ddgsolver::EigenMap<double, 3>(vpg.vertexNormals);

  // calculate mean curvature per vertex area by Laplacian matrix
  Hn = M_inv * L * positions / 2.0;
  vertexAreaGradientNormal = Hn.rowwise().normalized();
  auto projection = (vertexAreaGradientNormal.array()
                * vertexAngleNormal_e.array()).rowwise().sum();
  vertexAreaGradientNormal = (vertexAreaGradientNormal.array().colwise()
              *((projection > 0) - (projection < 0)).cast<double>()).matrix();

  // calculate laplacian H
  Eigen::Matrix<double, Eigen::Dynamic, 3> lap_H = M_inv * L * Hn;

  // initialize the spontaneous curvature matrix
  H0n = H0 * vertexAreaGradientNormal;

  // initialize and calculate intermediary result productTerms
  Eigen::Matrix<double, Eigen::Dynamic, 3> productTerms;
  productTerms.resize(n_vertices, 3);
  productTerms = 2 * ((Hn - H0n).array().colwise()
                * ((Hn.array() * Hn.array()).rowwise().sum()
                  + (H0n.array() * H0n.array()).rowwise().sum()
                 - KG.array())).matrix();

  // calculate bendingForce
  bendingForces_e = M * (-2.0 * Kb * (productTerms + lap_H));
}
} // end namespace ddgsolver
