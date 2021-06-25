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
  // E.BE = 0;
  // for (gcs::Vertex v : mesh->vertices()) {
  //   E.BE += Kb[v] *
  //           (vpg->vertexMeanCurvature(v) / vpg->vertexDualArea(v) - H0[v]) *
  //           (vpg->vertexMeanCurvature(v) / vpg->vertexDualArea(v) - H0[v]) *
  //           vpg->vertexDualArea(v);
  // }

  // Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference = H.raw() - H0.raw();
  // E.BE = P.Kb * H_difference.transpose() * M * H_difference;

  Eigen::Matrix<double, Eigen::Dynamic, 1> H_difference =
      abs(vpg->vertexMeanCurvatures.raw().array() /
              vpg->vertexDualAreas.raw().array() -
          H0.raw().array());
  E.BE = (Kb.raw().array() * vpg->vertexDualAreas.raw().array() *
          H_difference.array().square())
             .sum();

  // when considering topological changes, additional term of gauss curvature
  // E.BE = P.Kb * H_difference.transpose() * M * H_difference + P.KG * (M *
  // K).sum();
}

void System::computeSurfaceEnergy() {
  // cotan laplacian normal is exact for area variation
  double A_difference = surfaceArea - refSurfaceArea;
  E.sE = O.isConstantSurfaceTension ? F.surfaceTension * surfaceArea
                                    : F.surfaceTension * A_difference / 2 +
                                          P.lambdaSG * A_difference / 2;
}

void System::computePressureEnergy() {
  // Note: area weighted normal is exact volume variation
  if (O.isReducedVolume) {
    double V_difference = volume - refVolume * P.Vt;
    E.pE = -F.osmoticPressure * V_difference / 2 + P.lambdaV * V_difference / 2;
  } else if (O.isConstantOsmoticPressure) {
    E.pE = -F.osmoticPressure * volume;
  } else {
    E.pE = P.Kv * (P.cam * volume - log(P.cam * volume) - 1);
  }
}

void System::computeAdsorptionEnergy() {
  E.aE =
      P.epsilon *
      (vpg->vertexDualAreas.raw().array() * proteinDensity.raw().array()).sum();
}

void System::computeProteinInteriorPenaltyEnergy() {
  // interior method to constrain protein density to remain from 0 to 1
  E.inE =
      -P.lambdaPhi * ((proteinDensity.raw().array()).log().sum() +
                      (1 - proteinDensity.raw().array()).log().sum());
}

void System::computeDirichletEnergy() {
  if (false) {
    throw std::runtime_error(
        "computeDirichletEnergy: out of date implementation, "
        "shouldn't be called!");
    // scale the dH0 such that it is integrated over the edge
    // this is under the case where the resolution is low, WIP
    // auto dH0 = vpg->edgeLengths.raw().array() *  ((vpg->d0 *
    // H0.raw()).cwiseAbs()).array(); auto dH0 = (vpg->d0 *
    // H0.raw()).cwiseAbs();
    // E.dE = (vpg->hodge1Inverse * F.lineTension.raw()).sum();
  }

  // explicit dirichlet energy
  E.dE = 0;
  for (gcs::Face f : mesh->faces()) {
    E.dE += 0.5 * P.eta * proteinDensityGradient[f].norm2() * vpg->faceAreas[f];
  }

  // alternative dirichlet energy after integration by part
  // E.dE = 0.5 * P.eta * proteinDensity.raw().transpose() * vpg->cotanLaplacian
  // *
  //        proteinDensity.raw();
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
  if (P.Kb != 0 || P.Kbc != 0) {
    computeBendingEnergy();
  }
  if (P.Ksg != 0) {
    computeSurfaceEnergy();
  }
  if (P.Kv != 0) {
    computePressureEnergy();
  }
  if (P.epsilon != 0) {
    computeAdsorptionEnergy();
  }
  if (P.eta != 0) {
    computeDirichletEnergy();
  }
  if (P.Kf != 0) {
    computeExternalForceEnergy();
  }
  if (O.isProteinVariation) {
    computeProteinInteriorPenaltyEnergy();
  }
  E.potE = E.BE + E.sE + E.pE + E.aE + E.dE + E.exE + E.inE;
}

void System::computeFreeEnergy() {
  // zero all energy
  E = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  computeKineticEnergy();
  computePotentialEnergy();
  E.totalE = E.kE + E.potE;
}

void System::computeGradient(gcs::VertexData<double> &quantities,
                             gcs::FaceData<gc::Vector3> &gradient) {
  if ((quantities.raw().array() == quantities.raw()[0]).all()) {
    gradient.fill({0, 0, 0});
  } else {
    for (gcs::Face f : mesh->faces()) {
      gc::Vector3 normal = vpg->faceNormals[f];
      gc::Vector3 gradientVec{0, 0, 0};
      for (gcs::Halfedge he : f.adjacentHalfedges()) {
        gradientVec += quantities[he.next().tipVertex()] *
                       gc::cross(normal, vecFromHalfedge(he, *vpg));
      }
      gradient[f] = gradientVec / 2 / vpg->faceAreas[f];
    }
  }
}

} // namespace mem3dg
