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
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>
#include "mem3dg/solver/force.h"
#include "mem3dg/solver/meshops.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getBendingForces() {

  // map the MeshData to eigen matrix XXX_e
  auto bendingPressure_e = gc::EigenMap<double, 3>(bendingPressure);
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg.vertexNormals);
  auto positions = gc::EigenMap<double, 3>(vpg.inputVertexPositions);

  // Alias
  std::size_t n_vertices = (mesh.nVertices());

  // update the (tufted) mass and conformal Laplacian matrix
  if (isTuftedLaplacian) {
    getTuftedLaplacianAndMass(M, L, mesh, vpg, mollifyFactor);
  } else {
    M = vpg.vertexLumpedMassMatrix;
    L = vpg.cotanLaplacian;
  }
  // Cache the inverse mass matrix
  M_inv = (1 / (M.diagonal().array())).matrix().asDiagonal();

  // calculate mean curvature
  H = rowwiseDotProduct(M_inv * L * positions / 2.0, vertexAngleNormal_e);

  // Gaussian curvature per vertex Area
  Eigen::Matrix<double, Eigen::Dynamic, 1> KG =
      M_inv * vpg.vertexGaussianCurvatures.raw();

  // calculate the Laplacian of mean curvature H
  Eigen::Matrix<double, Eigen::Dynamic, 1> lap_H = M_inv * L * (H - H0);

  // initialize and calculate intermediary result scalerTerms
  Eigen::Matrix<double, Eigen::Dynamic, 1> scalerTerms =
      rowwiseProduct(H, H) + rowwiseProduct(H, H0) - KG;
  Eigen::Matrix<double, Eigen::Dynamic, 1> zeroMatrix;
  zeroMatrix.resize(n_vertices, 1);
  zeroMatrix.setZero();
  scalerTerms = scalerTerms.array().max(zeroMatrix.array());

  // initialize and calculate intermediary result productTerms
  Eigen::Matrix<double, Eigen::Dynamic, 1> productTerms;
  productTerms.resize(n_vertices, 1);
  productTerms = 2 * rowwiseProduct(scalerTerms, H - H0);

  // calculate bendingForce
  bendingPressure_e =
      rowwiseScaling(-2.0 * P.Kb * (productTerms + lap_H), vertexAngleNormal_e);
}

void Force::getChemicalPotential() { 
  
  Eigen::Matrix<double, Eigen::Dynamic, 1> proteinDensitySq =
      (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();

  H0 = (P.H0 * proteinDensitySq.array() / (1 + proteinDensitySq.array()))
           .matrix();

  Eigen::Matrix<double, Eigen::Dynamic, 1> dH0dphi =
      (2 * P.H0 * proteinDensity.raw().array() /
            ((1 + proteinDensitySq.array()) * (1 + proteinDensitySq.array()))).matrix();

  chemicalPotential.raw() =
      (P.epsilon - (2 * P.Kb * (H - H0)).array() * dH0dphi.array()).matrix();
}

} // end namespace ddgsolver
