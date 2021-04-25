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
#include <iostream>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void System::computeBendingEnergy() {

  // no particular reason, just experiment
  E.BE = 0;
  for (gcs::Vertex v : mesh->vertices()) {
    E.BE += Kb[v] *
            (vpg->vertexMeanCurvature(v) / vpg->vertexDualArea(v) - H0[v]) *
            (vpg->vertexMeanCurvature(v) / vpg->vertexDualArea(v) - H0[v]) *
            vpg->vertexDualArea(v);
  }

  // Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = H.raw() - H0.raw();
  // E.BE = P.Kb * H_difference.transpose() * M * H_difference;

  // Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference =
  //     abs(vpg->vertexMeanCurvatures.raw().array() /
  //             vpg->vertexDualAreas.raw().array() -
  //         H0.raw().array());
  // E.BE = (Kb.raw().array() * vpg->vertexDualAreas.raw().array() *
  //         H_difference.array().square())
  //            .sum();

  // when considering topological changes, additional term of gauss curvature
  // E.BE = P.Kb * H_difference.transpose() * M * H_difference + P.KG * (M *
  // K).sum();
}

void System::computeSurfaceEnergy() {
  // cotan laplacian normal is exact for area variation

  double A_difference = surfaceArea - refSurfaceArea;
  if (mesh->hasBoundary()) {
    // non moving boundary
    E.sE = P.Ksg * A_difference;
  } else {
    E.sE = P.Ksg * A_difference * A_difference / refSurfaceArea / 2 +
           P.lambdaSG * A_difference;
  }
}

void System::computePressureEnergy() {
  // Note: angled weight normal is exact volume variation

  double V_difference = volume - refVolume * P.Vt;
  if (mesh->hasBoundary()) {
    E.pE = -P.Kv * V_difference;
  } else if (O.isReducedVolume) {
    E.pE = P.Kv * V_difference * V_difference / (refVolume * P.Vt) / 2 +
           P.lambdaV * V_difference;
  } else {
    E.pE = P.Kv * (P.cam * volume - log(P.cam * volume) - 1);
  }
}

void System::computeChemicalEnergy() {
  E.cE = (vpg->vertexLumpedMassMatrix * P.epsilon * proteinDensity.raw()).sum();
}

void System::computeLineTensionEnergy() {
  // scale the dH0 such that it is integrated over the edge
  // this is under the case where the resolution is low, WIP
  // auto dH0 = vpg->edgeLengths.raw().array() *  ((vpg->d0 *
  // H0.raw()).cwiseAbs()).array(); auto dH0 = (vpg->d0 * H0.raw()).cwiseAbs();
  E.lE = (vpg->hodge1Inverse * F.lineTension.raw()).sum();
}

void System::computeExternalForceEnergy() {
  E.exE = -rowwiseDotProduct(
               rowwiseScaling(F.externalForce.raw(),
                              gc::EigenMap<double, 3>(vpg->vertexNormals)),
               gc::EigenMap<double, 3>(vpg->inputVertexPositions))
               .sum();
}

void System::computeKineticEnergy() {
  auto velocity = gc::EigenMap<double, 3>(vel);
  // auto velocity =
  //     rowwiseDotProduct(gc::EigenMap<double, 3>(vel),
  //                       gc::EigenMap<double, 3>(vpg->inputVertexPositions));
  E.kE = 0.5 * (vpg->vertexLumpedMassMatrix *
                (velocity.array() * velocity.array()).matrix())
                   .sum();
}

void System::computePotentialEnergy() {
  computeBendingEnergy();
  if (P.Ksg != 0) {
    computeSurfaceEnergy();
  }
  if (P.Kv != 0) {
    computePressureEnergy();
  }
  if (O.isProtein) {
    computeChemicalEnergy();
  }
  if (P.eta != 0) {
    computeLineTensionEnergy();
  }
  if (P.Kf != 0) {
    computeExternalForceEnergy();
  }
  E.potE = E.BE + E.sE + E.pE + E.cE + E.lE + E.exE;
}

void System::computeFreeEnergy() {
  // zero all energy
  E = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  computeKineticEnergy();
  computePotentialEnergy();
  E.totalE = E.kE + E.potE;
}

// double
// System::computeL1Norm(Eigen::Matrix<double, Eigen::Dynamic, 1> force) const {
//   // return sqrt((rowwiseDotProduct(force, force)).sum() / surfaceArea);

//   // L2 Norm
//   // return sqrt(rowwiseDotProduct(force, force).sum()) / surfaceArea;

//   // auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg.vertexNormals);
//   // return (M * rowwiseDotProduct(pressure,
//   // vertexAngleNormal_e).cwiseAbs()).sum() / surfaceArea;

//   // L1 Norm
//   std::cout << "in compute error: " << force.cwiseAbs().sum() << std::endl;
//   return force.cwiseAbs().sum() / surfaceArea;
// }

double
System::computeL1Norm(Eigen::Matrix<double, Eigen::Dynamic, 1> &force) const {

  return force.cwiseAbs().sum() / surfaceArea;
}

double
System::computeL1Norm(Eigen::Matrix<double, Eigen::Dynamic, 1> &&force) const {
  // L1 Norm
  return force.cwiseAbs().sum() / surfaceArea;

  // Mask the mutated (unsmooth) set of data
  // return ((smoothingMask.raw().array() == false).cast<double>() *
  // force.array())
  //            .cwiseAbs()
  //            .sum() /
  //        (smoothingMask.raw().array() == false).cast<double>().sum();
}

} // namespace mem3dg
