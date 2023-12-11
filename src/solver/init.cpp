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

#include "mem3dg/solver/system.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {
namespace solver {

void System::initialize(bool ifMutateMesh) {
  checkConfiguration();

  pcg_extras::seed_seq_from<std::random_device> seed_source;
  rng = pcg32(seed_source);

  bool ifUpdateNotableVertex = true, ifUpdateGeodesics = true,
       ifUpdateProteinDensityDistribution = true, ifUpdateMask = true;
  updatePrescription(ifMutateMesh, ifUpdateNotableVertex, ifUpdateGeodesics,
                     ifUpdateProteinDensityDistribution, ifUpdateMask);

  updateConfigurations();
  geometry.updateReferenceConfigurations();

  if (geometry.mesh->hasBoundary()) {
    boundaryForceMask(*geometry.mesh, forces.forceMask,
                      parameters.boundary.shapeBoundaryCondition);
    boundaryProteinMask(*geometry.mesh, forces.proteinMask,
                        parameters.boundary.proteinBoundaryCondition);
  }
  computeConservativeForcing();
  computeTotalEnergy();
  mechErrorNorm = toMatrix(forces.mechanicalForceVec).norm();
  chemErrorNorm = forces.chemicalPotential.raw().norm();
}

void System::checkConfiguration() {
  parameters.checkParameters(geometry.mesh->hasBoundary(),
                             geometry.mesh->nVertices());
  meshProcessor.summarizeStatus();
  if (meshProcessor.isMeshMutate && !parameters.variation.isShapeVariation) {
    mem3dg_runtime_error("Mesh mutation operation not allowed for non shape "
                         "variation simulation");
  }
  if (parameters.selfAvoidance.mu != 0) {
    for (std::size_t i = 0; i < geometry.mesh->nVertices(); ++i) {
      gc::Vertex vi{geometry.mesh->vertex(i)};
      gc::VertexData<bool> neighborList(*geometry.mesh, false);
      meshProcessor.meshMutator.markVertices(neighborList, vi,
                                             parameters.selfAvoidance.n);
      for (std::size_t j = i + 1; j < geometry.mesh->nVertices(); ++j) {
        if (neighborList[j])
          continue;
        gc::Vertex vj{geometry.mesh->vertex(j)};
        gc::Vector3 r = geometry.vpg->inputVertexPositions[vj] -
                        geometry.vpg->inputVertexPositions[vi];
        double distance = gc::norm(r);
        if (distance < parameters.selfAvoidance.d)
          mem3dg_runtime_error(
              "Input mesh violates the self avoidance constraint!");
      }
    }
  }

  if (((proteinDensity - proteinDensity[0]).raw().norm() == 0) &&
      (parameters.protein.prescribeProteinDensityDistribution ==
       NULL)) { // homogeneous distribution
    if (parameters.variation.isProteinVariation) {
      if (proteinDensity[0] < 0 || proteinDensity[0] > 1)
        mem3dg_runtime_error("{0<=phi<=1}");
    } else {
      if (proteinDensity[0] != 1 || parameters.bending.Kb != 0 ||
          parameters.dirichlet.eta != 0 || parameters.adsorption.epsilon != 0 ||
          parameters.aggregation.chi != 0)
        mem3dg_runtime_warning(
            "For homogeneous membrane simulation, good practice is to set "
            "proteinDensity = 1, Kb = 0, eta  = 0, "
            "epsilon = 0, chi = "
            "0 to "
            "avoid ambiguity & save computation!");
    }
  }
}

void System::updateConfigurations() {
  geometry.updateConfigurations();

  // compute face gradient of protein density
  if (parameters.dirichlet.eta != 0) {
    geometry.computeFaceTangentialDerivative(proteinDensity,
                                             proteinDensityGradient);
  }

  // Update protein density dependent quantities
  if (parameters.bending.relation == "linear") {
    H0.raw() = proteinDensity.raw() * parameters.bending.H0c;
    Kb.raw() = parameters.bending.Kb +
               parameters.bending.Kbc * proteinDensity.raw().array();
    Kd.raw() = parameters.bending.Kd +
               parameters.bending.Kdc * proteinDensity.raw().array();
  } else if (parameters.bending.relation == "hill") {
    EigenVectorX1d proteinDensitySq =
        (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();
    H0.raw() = parameters.bending.H0c * proteinDensitySq.array() /
               (1 + proteinDensitySq.array());
    Kb.raw() = parameters.bending.Kb + parameters.bending.Kbc *
                                           proteinDensitySq.array() /
                                           (1 + proteinDensitySq.array());
    Kd.raw() = parameters.bending.Kd + parameters.bending.Kdc *
                                           proteinDensitySq.array() /
                                           (1 + proteinDensitySq.array());
  } else {
    mem3dg_runtime_error("updateVertexPosition: P.relation is invalid option!");
  }

  /// volume and osmotic pressure
  if (parameters.osmotic.form != NULL)
    std::tie(forces.osmoticPressure, energy.pressureEnergy) =
        parameters.osmotic.form(geometry.volume);

  // area and surface tension
  if (parameters.tension.form != NULL)
    std::tie(forces.surfaceTension, energy.surfaceEnergy) =
        parameters.tension.form(geometry.surfaceArea);
}

bool System::updatePrescription(std::map<std::string, double> &lastUpdateTime,
                                double timeStep) {
  bool ifMutateMesh = (time - lastUpdateTime["mutateMesh"] >
                       (meshProcessor.meshMutator.mutateMeshPeriod * timeStep)),
       ifUpdateNotableVertex =
           (time - lastUpdateTime["notableVertex"] >
            (parameters.point.updateNotableVertexPeriod * timeStep)),
       ifUpdateGeodesics =
           (time - lastUpdateTime["geodesics"] >
            (parameters.point.updateGeodesicsPeriod * timeStep)),
       ifUpdateProteinDensityDistribution =
           (time - lastUpdateTime["protein"] >
            (parameters.protein.updateProteinDensityDistributionPeriod *
             timeStep)),
       ifUpdateMask = (time - lastUpdateTime["mask"] >
                       (parameters.variation.updateMaskPeriod * timeStep));

  bool updated =
      updatePrescription(ifMutateMesh, ifUpdateNotableVertex, ifUpdateGeodesics,
                         ifUpdateProteinDensityDistribution, ifUpdateMask);
  if (ifMutateMesh)
    lastUpdateTime["mutateMesh"] = time;
  if (ifUpdateNotableVertex)
    lastUpdateTime["notableVertex"] = time;
  if (ifUpdateGeodesics)
    lastUpdateTime["geodesics"] = time;
  if (ifUpdateProteinDensityDistribution)
    lastUpdateTime["protein"] = time;
  if (ifUpdateMask)
    lastUpdateTime["mask"] = time;

  return updated;
}

bool System::updatePrescription(bool &ifMutateMesh, bool &ifUpdateNotableVertex,
                                bool &ifUpdateGeodesics,
                                bool &ifUpdateProteinDensityDistribution,
                                bool &ifUpdateMask) {

  if (ifMutateMesh) {
    if (meshProcessor.isMeshMutate) {
      mutateMesh();
      updateConfigurations();
      geometry.updateReferenceMeshFromVPG();
      geometry.updateReferenceConfigurations();
    } else {
      ifMutateMesh = false;
      // mem3dg_runtime_warning("Meshmutator is not activated!");
    }
  }

  if (ifUpdateNotableVertex) {
    if (parameters.point.prescribeNotableVertex != NULL) {
      geometry.notableVertex.raw() = parameters.point.prescribeNotableVertex(
          geometry.mesh->getFaceVertexMatrix<std::size_t>(),
          toMatrix(geometry.vpg->vertexPositions),
          geometry.geodesicDistance.raw());
    } else {
      ifUpdateNotableVertex = false;
      // mem3dg_runtime_warning("Parameter.point.prescribeNotableVertex is
      // NULL!");
    }
  }

  if (ifUpdateGeodesics) {
    geometry.computeGeodesicDistance();
  }

  if (ifUpdateProteinDensityDistribution) {
    if (parameters.protein.prescribeProteinDensityDistribution != NULL) {
      proteinDensity.raw() =
          parameters.protein.prescribeProteinDensityDistribution(
              time, geometry.vpg->vertexMeanCurvatures.raw(),
              geometry.geodesicDistance.raw());
    } else {
      ifUpdateProteinDensityDistribution = false;
      // mem3dg_runtime_warning("Parameter protein form is NULL!")
    }
  }

  if (ifUpdateMask) {
    if (parameters.variation.geodesicMask != -1) {
      prescribeGeodesicMasks();
    } else {
      ifUpdateMask = false;
      // mem3dg_runtime_warning("geodesicMask not activated!")
    }
  }

  return ifMutateMesh || ifUpdateNotableVertex || ifUpdateGeodesics ||
         ifUpdateProteinDensityDistribution || ifUpdateMask;
}

} // namespace solver
} // namespace mem3dg
