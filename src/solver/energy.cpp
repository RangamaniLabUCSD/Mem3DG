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

#include "Eigen/src/Core/GlobalFunctions.h"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "mem3dg/constants.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/type_utilities.h"

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
  energy.bendingEnergy =
      (Kb.raw().array() * vpg->vertexDualAreas.raw().array() *
       H_difference.array().square())
          .sum();

  // when considering topological changes, additional term of gauss curvature
  // E.BE = P.Kb * H_difference.transpose() * M * H_difference + P.KG * (M *
  // K).sum();
}

void System::computeDeviatoricEnergy() {
  energy.deviatoricEnergy =
      (Kd.raw().array() * (vpg->vertexMeanCurvatures.raw().array().square() /
                               vpg->vertexDualAreas.raw().array() -
                           vpg->vertexGaussianCurvatures.raw().array()))
          .sum();
  // energy.deviatoricEnergy =
  //     (Kd.raw().array() * (vpg->vertexMeanCurvatures.raw().array().square() /
  //                          vpg->vertexDualAreas.raw().array()))
  //         .sum();
  // energy.deviatoricEnergy =
  //     (Kd.raw().array() *
  //     (-vpg->vertexGaussianCurvatures.raw().array())).sum();
}

void System::computeSurfaceEnergy() {
  // cotan laplacian normal is exact for area variation
  double A_difference = surfaceArea - parameters.tension.At;
  energy.surfaceEnergy =
      parameters.tension.isConstantSurfaceTension
          ? forces.surfaceTension * surfaceArea
          : forces.surfaceTension * A_difference / 2 +
                parameters.tension.lambdaSG * A_difference / 2;
}

void System::computePressureEnergy() {
  // Note: area weighted normal is exact volume variation
  if (parameters.osmotic.isPreferredVolume) {
    double V_difference = volume - parameters.osmotic.Vt;
    energy.pressureEnergy = -forces.osmoticPressure * V_difference / 2 +
                            parameters.osmotic.lambdaV * V_difference / 2;
  } else if (parameters.osmotic.isConstantOsmoticPressure) {
    energy.pressureEnergy = -forces.osmoticPressure * volume;
  } else {
    double ratio = parameters.osmotic.cam * volume / parameters.osmotic.n;
    energy.pressureEnergy = mem3dg::constants::i * mem3dg::constants::R *
                            parameters.temperature * parameters.osmotic.n *
                            (ratio - log(ratio) - 1);
  }
}

// void System::computeAdsorptionEnergy() {
//   energy.adsorptionEnergy =
//       parameters.adsorption.epsilon * (proteinDensity.raw().array()).sum();
// }

// void System::computeAggregationEnergy() {
//   energy.aggregationEnergy =
//       parameters.aggregation.chi *
//       (proteinDensity.raw().array() * proteinDensity.raw().array()).sum();
// }

void System::computeAdsorptionEnergy() {
  energy.adsorptionEnergy =
      parameters.adsorption.epsilon *
      (vpg->vertexDualAreas.raw().array() * proteinDensity.raw().array()).sum();
}

void System::computeAggregationEnergy() {
  energy.aggregationEnergy =
      parameters.aggregation.chi *
      (vpg->vertexDualAreas.raw().array() * proteinDensity.raw().array() *
       proteinDensity.raw().array())
          .sum();
}

void System::computeProteinInteriorPenalty() {
  // interior method to constrain protein density to remain from 0 to 1
  energy.proteinInteriorPenalty =
      -parameters.proteinDistribution.lambdaPhi *
      ((proteinDensity.raw().array()).log().sum() +
       (1 - proteinDensity.raw().array()).log().sum());
}

void System::computeSelfAvoidanceEnergy() {
  const double d0 = parameters.selfAvoidance.d;
  const double mu = parameters.selfAvoidance.mu;
  const double n = parameters.selfAvoidance.n;
  double e = 0.0;
  projectedCollideTime = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex vi{mesh->vertex(i)};
    gc::VertexData<bool> neighborList(*mesh, false);
    meshProcessor.meshMutator.markVertices(neighborList, vi, n);
    for (std::size_t j = i + 1; j < mesh->nVertices(); ++j) {
      if (neighborList[j])
        continue;
      gc::Vertex vj{mesh->vertex(j)};

      // double penalty = mu * vpg->vertexDualAreas[vi] * proteinDensity[vi] *
      //                  vpg->vertexDualAreas[vj] * proteinDensity[vj];
      double penalty = mu * proteinDensity[vi] * proteinDensity[vj];
      // double penalty = mu;
      // double penalty = mu * vpg->vertexDualAreas[vi] *
      // vpg->vertexDualAreas[vj];

      gc::Vector3 r =
          vpg->inputVertexPositions[vj] - vpg->inputVertexPositions[vi];
      double distance = gc::norm(r) - d0;
      double collideTime = distance / gc::dot(velocity[vi] - velocity[vj], r);
      if (collideTime < projectedCollideTime &&
          gc::dot(velocity[vi] - velocity[vj], r) > 0)
        projectedCollideTime = collideTime;
      // e -= penalty * log(distance);
      e += penalty / distance;
    }
  }
  if (projectedCollideTime == std::numeric_limits<double>::max())
    projectedCollideTime = 0;
  energy.selfAvoidancePenalty = e;
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
  energy.dirichletEnergy = 0;
  for (gcs::Face f : mesh->faces()) {
    energy.dirichletEnergy += 0.5 * parameters.dirichlet.eta *
                              proteinDensityGradient[f].norm2() *
                              vpg->faceAreas[f];
  }

  // alternative dirichlet energy after integration by part
  // E.dE = 0.5 * P.eta * proteinDensity.raw().transpose() *
  // vpg->cotanLaplacian
  // *
  //        proteinDensity.raw();
}

double System::computePotentialEnergy() {
  // fundamental internal potential energy
  energy.bendingEnergy = 0;
  energy.deviatoricEnergy = 0;
  energy.surfaceEnergy = 0;
  energy.pressureEnergy = 0;
  energy.adsorptionEnergy = 0;
  energy.aggregationEnergy = 0;
  energy.dirichletEnergy = 0;
  energy.selfAvoidancePenalty = 0;
  energy.proteinInteriorPenalty = 0;

  computeBendingEnergy();
  computeSurfaceEnergy();
  computePressureEnergy();

  // optional internal potential energy
  computeDeviatoricEnergy();
  if (parameters.adsorption.epsilon != 0) {
    computeAdsorptionEnergy();
  }
  if (parameters.aggregation.chi != 0) {
    computeAggregationEnergy();
  }
  if (parameters.dirichlet.eta != 0) {
    computeDirichletEnergy();
  }
  if (parameters.selfAvoidance.mu != 0) {
    computeSelfAvoidanceEnergy();
  }
  if (parameters.variation.isProteinVariation &&
      parameters.proteinDistribution.lambdaPhi != 0) {
    computeProteinInteriorPenalty();
  }

  // summerize internal potential energy
  energy.potentialEnergy =
      energy.bendingEnergy + energy.deviatoricEnergy + energy.surfaceEnergy +
      energy.pressureEnergy + energy.adsorptionEnergy + energy.dirichletEnergy +
      energy.aggregationEnergy + energy.selfAvoidancePenalty +
      energy.proteinInteriorPenalty;
  return energy.potentialEnergy;
}

double System::computeIntegratedPower(double dt) {
  prescribeExternalForce();
  return dt * rowwiseDotProduct(toMatrix(forces.externalForceVec),
                                toMatrix(velocity))
                  .sum();
}
double System::computeIntegratedPower(double dt, EigenVectorX3dr &&velocity) {
  prescribeExternalForce();
  return dt *
         rowwiseDotProduct(toMatrix(forces.externalForceVec), velocity).sum();
}

double System::computeExternalWork(double currentTime, double dt) {
  energy.time = currentTime;
  energy.externalWork += computeIntegratedPower(dt);
  return energy.externalWork;
}

double System::computeKineticEnergy() {
  energy.kineticEnergy =
      0.5 * Eigen::square(toMatrix(velocity).array()).matrix().sum();
  // energy.kineticEnergy =
  //     0.5 * rowwiseScalarProduct(toMatrix(vpg->vertexDualAreas).array(),
  //                                Eigen::square(toMatrix(velocity).array()))
  //               .sum();
  return energy.kineticEnergy;
}

double System::computeTotalEnergy() {
  computePotentialEnergy();
  computeKineticEnergy();
  if (time == energy.time) {
    energy.totalEnergy =
        energy.kineticEnergy + energy.potentialEnergy - energy.externalWork;
  } else if (parameters.external.Kf == 0) {
    energy.totalEnergy = energy.kineticEnergy + energy.potentialEnergy;
  } else {
    mem3dg_runtime_error("energy.externalWork not updated!")
  }
  return energy.totalEnergy;
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
