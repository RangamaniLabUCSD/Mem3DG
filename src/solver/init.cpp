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

void System::initialize(std::size_t nMutation, bool ifMute) {
  checkConfiguration();
  initializeConstants(ifMute);
  meshProcessor.summarizeStatus();
  if (!meshProcessor.isMeshMutate && nMutation != 0) {
    mem3dg_runtime_message("mesh mutator not activated!");
  } else {
    updateConfigurations();
    updateReferenceConfigurations();
    mutateMesh(nMutation);
  }
  if (nMutation != 0) {
    updateConfigurations();
    refVpg = vpg->copy();
    updateReferenceConfigurations();
  }
}

void System::checkConfiguration() {

  isOpenMesh = mesh->hasBoundary();
  parameters.checkParameters(isOpenMesh, mesh->nVertices());
  meshProcessor.summarizeStatus();
  if (meshProcessor.isMeshMutate && !parameters.variation.isShapeVariation) {
    mem3dg_runtime_error("Mesh mutation operation not allowed for non shape "
                         "variation simulation");
  }
  if (!isOpenMesh && mesh->genus() != 0) {
    mem3dg_runtime_error(
        "Do not support closed mesh with nonzero number of genus!")
  }
  if (parameters.point.pt.rows() == 2 && !isOpenMesh) {
    mem3dg_runtime_message(
        "specifying x-y coordinate on closed surface may "
        "lead to ambiguity! Please check by visualizing it first!");
  }
  if (parameters.selfAvoidance.mu != 0) {
    for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
      gc::Vertex vi{mesh->vertex(i)};
      gc::VertexData<bool> neighborList(*mesh, false);
      meshProcessor.meshMutator.markVertices(neighborList, vi,
                                             parameters.selfAvoidance.n);
      for (std::size_t j = i + 1; j < mesh->nVertices(); ++j) {
        if (neighborList[j])
          continue;
        gc::Vertex vj{mesh->vertex(j)};
        gc::Vector3 r =
            vpg->inputVertexPositions[vj] - vpg->inputVertexPositions[vi];
        double distance = gc::norm(r);
        if (distance < parameters.selfAvoidance.d)
          mem3dg_runtime_error(
              "Input mesh violates the self avoidance constraint!");
      }
    }
  }

  if (((proteinDensity - proteinDensity[0]).raw().norm() == 0) &&
      (!parameters.protein.ifPrescribe)) { // homogeneous distribution
    if (parameters.variation.isProteinVariation) {
      if (proteinDensity[0] < 0 || proteinDensity[0] > 1)
        mem3dg_runtime_error("{0<=phi<=1}");
    } else {
      if (proteinDensity[0] != 1 || parameters.bending.Kb != 0 ||
          parameters.dirichlet.eta != 0 || parameters.adsorption.epsilon != 0 ||
          parameters.aggregation.chi != 0)
        mem3dg_runtime_message(
            "For homogenous membrane simulation, good practice is to set "
            "proteinDensity = 1, Kb = 0, eta  = 0, "
            "epsilon = 0, chi = "
            "0 to "
            "avoid ambiguity & save computation!");
    }
  }
}

void System::initializeConstants(bool ifMute) {
  pcg_extras::seed_seq_from<std::random_device> seed_source;
  rng = pcg32(seed_source);

  if (parameters.point.isFloatVertex) {
    findFloatCenter();
  } else {
    findVertexCenter();
  }
  geodesicDistance.raw() = computeGeodesicDistance();
  prescribeGeodesicMasks();

  if (mesh->hasBoundary()) {
    boundaryForceMask(*mesh, forces.forceMask,
                      parameters.boundary.shapeBoundaryCondition);
    boundaryProteinMask(*mesh, forces.proteinMask,
                        parameters.boundary.proteinBoundaryCondition);
  }

  surfaceArea = vpg->faceAreas.raw().sum() + parameters.tension.A_res;
  volume = getMeshVolume(*mesh, *vpg, true) + parameters.osmotic.V_res;
  if (!ifMute) {
    std::cout << "area_init = " << surfaceArea << std::endl;
    std::cout << "vol_init = " << volume << std::endl;
    std::cout << "Characteristic volume wrt to At = "
              << (isOpenMesh
                      ? parameters.osmotic.V_res
                      : std::pow(parameters.tension.At / constants::PI / 4,
                                 1.5) *
                            (4 * constants::PI / 3))
              << std::endl;
  }
}

void System::updateReferenceConfigurations() {
  refVpg->requireEdgeLengths();
  refVpg->requireFaceAreas();
  for (std::size_t i = 0; i < mesh->nHalfedges(); ++i) {
    gcs::Halfedge he{mesh->halfedge(i)};
    refLcrs[he] = computeLengthCrossRatio(*refVpg, he);
  }
}

void System::updateConfigurations() {

  // refresh cached quantities after regularization
  vpg->refreshQuantities();

  // compute face gradient of protein density
  if (parameters.dirichlet.eta != 0) {
    computeFaceTangentialDerivative(proteinDensity, proteinDensityGradient);
  }

  // Update protein density dependent quantities
  if (parameters.bending.relation == "linear") {
    H0.raw() = proteinDensity.raw() * parameters.bending.H0c;
    Kb.raw() = parameters.bending.Kb +
               parameters.bending.Kbc * proteinDensity.raw().array();
    Kd.raw() = parameters.bending.Kd +
               parameters.bending.Kdc * proteinDensity.raw().array();
  } else if (parameters.bending.relation == "hill") {
    Eigen::Matrix<double, Eigen::Dynamic, 1> proteinDensitySq =
        (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();
    H0.raw() = (parameters.bending.H0c * proteinDensitySq.array() /
                (1 + proteinDensitySq.array()))
                   .matrix();
    Kb.raw() = (parameters.bending.Kb + parameters.bending.Kbc *
                                            proteinDensitySq.array() /
                                            (1 + proteinDensitySq.array()))
                   .matrix();
    Kd.raw() = (parameters.bending.Kd + parameters.bending.Kdc *
                                            proteinDensitySq.array() /
                                            (1 + proteinDensitySq.array()))
                   .matrix();
  } else {
    mem3dg_runtime_error("updateVertexPosition: P.relation is invalid option!");
  }

  /// initialize/update enclosed volume
  volume = getMeshVolume(*mesh, *vpg, true) + parameters.osmotic.V_res;

  // update global osmotic pressure
  if (parameters.osmotic.isPreferredVolume) {
    forces.osmoticPressure =
        -(parameters.osmotic.Kv * (volume - parameters.osmotic.Vt) /
              parameters.osmotic.Vt / parameters.osmotic.Vt +
          parameters.osmotic.lambdaV);
  } else if (parameters.osmotic.isConstantOsmoticPressure) {
    forces.osmoticPressure = parameters.osmotic.Kv;
  } else {
    forces.osmoticPressure =
        mem3dg::constants::i * mem3dg::constants::R * parameters.temperature *
        (parameters.osmotic.n / volume - parameters.osmotic.cam);
  }

  // initialize/update total surface area
  surfaceArea = vpg->faceAreas.raw().sum() + parameters.tension.A_res;

  // update global surface tension
  forces.surfaceTension = parameters.tension.isConstantSurfaceTension
                              ? parameters.tension.Ksg
                              : parameters.tension.Ksg *
                                        (surfaceArea - parameters.tension.At) /
                                        parameters.tension.At +
                                    parameters.tension.lambdaSG;
}
} // namespace solver
} // namespace mem3dg
