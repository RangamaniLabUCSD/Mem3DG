// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2020:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include "mem3dg/solver/force.h"
#include "mem3dg/solver/meshops.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getBendingForces() {

  // update the (tufted) mass and conformal Laplacian matrix
  if (isTuftedLaplacian) {
    getTuftedLaplacianAndMass(M, L, mesh, vpg, mollifyFactor);
  } else {
    M = vpg.vertexLumpedMassMatrix;
    L = vpg.cotanLaplacian;
  }

  // Cache the inverse mass matrix
  M_inv = (1 / (M.diagonal().array())).matrix().asDiagonal();

  // Gaussian curvature per vertex Area
  Eigen::Matrix<double, Eigen::Dynamic, 1> KG =
      M_inv * vpg.vertexGaussianCurvatures.raw();

  // number of vertices for convenience
  std::size_t n_vertices = (mesh.nVertices());

  // map ivp to eigen matrix position
  auto positions = EigenMap<double, 3>(vpg.inputVertexPositions);

  // map the VertexData bendingForces to eigen matrix bendingForces_e
  auto bendingForces_e = EigenMap<double, 3>(bendingForces);
  bendingForces_e.setZero();

  // the build-in angle-weighted vertex normal
  auto vertexAngleNormal_e = EigenMap<double, 3>(vpg.vertexNormals);

  // calculate mean curvature and map it to angle-weighted normal
  Hn = rowwiseScaling(
      rowwiseDotProduct(vertexAngleNormal_e, M_inv * L * positions / 2.0),
      vertexAngleNormal_e);

  // calculate the Laplacian of mean curvature H
  Eigen::Matrix<double, Eigen::Dynamic, 3> lap_H = M_inv * L * Hn;

  // initialize the spontaneous curvature matrix
  H0n = P.H0 * vertexAngleNormal_e;

  // initialize and calculate intermediary result scalerTerms, set to zero if
  // negative
  Eigen::Matrix<double, Eigen::Dynamic, 1> scalerTerms =
      rowwiseDotProduct(Hn, Hn) + rowwiseDotProduct(Hn, H0n) - KG;
  Eigen::Matrix<double, Eigen::Dynamic, 1> zeroMatrix;
  zeroMatrix.resize(n_vertices, 1);
  zeroMatrix.setZero();
  scalerTerms = scalerTerms.array().max(zeroMatrix.array());

  // initialize and calculate intermediary result productTerms
  Eigen::Matrix<double, Eigen::Dynamic, 3> productTerms;
  productTerms.resize(n_vertices, 3);
  productTerms = 2 * rowwiseScaling(scalerTerms, Hn - H0n);

  // calculate bendingForce
  bendingForces_e = M * (-2.0 * P.Kb * (productTerms + lap_H));
}
} // end namespace ddgsolver
