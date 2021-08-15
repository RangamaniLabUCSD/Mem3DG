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
#include "mem3dg/constants.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {

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
  energy.BE = (Kb.raw().array() * vpg->vertexDualAreas.raw().array() *
               H_difference.array().square())
                  .sum();

  // when considering topological changes, additional term of gauss curvature
  // E.BE = P.Kb * H_difference.transpose() * M * H_difference + P.KG * (M *
  // K).sum();
}

void System::computeSurfaceEnergy() {
  // cotan laplacian normal is exact for area variation
  double A_difference = surfaceArea - parameters.tension.At;
  energy.sE = parameters.tension.isConstantSurfaceTension
                  ? forces.surfaceTension * surfaceArea
                  : forces.surfaceTension * A_difference / 2 +
                        parameters.tension.lambdaSG * A_difference / 2;
}

void System::computePressureEnergy() {
  // Note: area weighted normal is exact volume variation
  if (parameters.osmotic.isPreferredVolume) {
    double V_difference = volume - parameters.osmotic.Vt;
    energy.pE = -forces.osmoticPressure * V_difference / 2 +
                parameters.osmotic.lambdaV * V_difference / 2;
  } else if (parameters.osmotic.isConstantOsmoticPressure) {
    energy.pE = -forces.osmoticPressure * volume;
  } else {
    double ratio = parameters.osmotic.cam * volume / parameters.osmotic.n;
    energy.pE = mem3dg::constants::i * mem3dg::constants::R * parameters.temperature *
                parameters.osmotic.n * (ratio - log(ratio) - 1);
  }
}

void System::computeAdsorptionEnergy() {
  energy.aE =
      parameters.adsorption.epsilon *
      (vpg->vertexDualAreas.raw().array() * proteinDensity.raw().array()).sum();
}

void System::computeProteinInteriorPenaltyEnergy() {
  // interior method to constrain protein density to remain from 0 to 1
  energy.inE = -parameters.proteinDistribution.lambdaPhi *
               ((proteinDensity.raw().array()).log().sum() +
                (1 - proteinDensity.raw().array()).log().sum());
}

void System::computeDirichletEnergy() {
  if (false) {
    mem3dg_runtime_error("computeDirichletEnergy: out of date implementation, "
                         "shouldn't be called!");
    // scale the dH0 such that it is integrated over the edge
    // this is under the case where the resolution is low, WIP
    // auto dH0 = vpg->edgeLengths.raw().array() *  ((vpg->d0 *
    // H0.raw()).cwiseAbs()).array(); auto dH0 = (vpg->d0 *
    // H0.raw()).cwiseAbs();
    // E.dE = (vpg->hodge1Inverse * F.lineTension.raw()).sum();
  }

  // explicit dirichlet energy
  energy.dE = 0;
  for (gcs::Face f : mesh->faces()) {
    energy.dE += 0.5 * parameters.dirichlet.eta *
                 proteinDensityGradient[f].norm2() * vpg->faceAreas[f];
  }

  // alternative dirichlet energy after integration by part
  // E.dE = 0.5 * P.eta * proteinDensity.raw().transpose() * vpg->cotanLaplacian
  // *
  //        proteinDensity.raw();
}

void System::computeExternalForceEnergy() {
  energy.exE =
      -rowwiseDotProduct(
           rowwiseScalarProduct(forces.externalForce.raw(),
                                gc::EigenMap<double, 3>(vpg->vertexNormals)),
           gc::EigenMap<double, 3>(vpg->inputVertexPositions))
           .sum();
}

void System::computeKineticEnergy() {
  auto vel = gc::EigenMap<double, 3>(velocity);
  // auto velocity =
  //     rowwiseDotProduct(gc::EigenMap<double, 3>(vel),
  //                       gc::EigenMap<double, 3>(vpg->inputVertexPositions));
  energy.kE =
      0.5 * (vpg->vertexLumpedMassMatrix * (vel.array() * vel.array()).matrix())
                .sum();
}

void System::computePotentialEnergy() {
  if (parameters.bending.Kb != 0 || parameters.bending.Kbc != 0) {
    computeBendingEnergy();
  }
  if (parameters.tension.Ksg != 0) {
    computeSurfaceEnergy();
  }
  if (parameters.osmotic.Kv != 0) {
    computePressureEnergy();
  }
  if (parameters.adsorption.epsilon != 0) {
    computeAdsorptionEnergy();
  }
  if (parameters.dirichlet.eta != 0) {
    computeDirichletEnergy();
  }
  if (parameters.external.Kf != 0) {
    computeExternalForceEnergy();
  }
  if (parameters.variation.isProteinVariation) {
    computeProteinInteriorPenaltyEnergy();
  }
  energy.potE = energy.BE + energy.sE + energy.pE + energy.aE + energy.dE +
                energy.exE + energy.inE;
}

void System::computeFreeEnergy() {
  // zero all energy
  energy = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  computeKineticEnergy();
  computePotentialEnergy();
  energy.totalE = energy.kE + energy.potE;
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

} // namespace solver
} // namespace mem3dg
