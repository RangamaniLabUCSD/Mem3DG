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
#include "mem3dg/solver/parameters.h"
#include "mem3dg/macros.h"
#include <cstddef>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {
namespace solver {

void Parameters::Bending::checkParameters() {
  if (alpha != 0) {
    if (D == 0)
      mem3dg_runtime_error("Membrane thickness D has to be larger than 0!");
  }
}

void Parameters::Boundary::checkParameters() {
  if (shapeBoundaryCondition != "roller" && shapeBoundaryCondition != "pin" &&
      shapeBoundaryCondition != "fixed" && shapeBoundaryCondition != "none") {
    mem3dg_runtime_error("Invalid option for shapeBoundaryCondition!");
  }
  if (proteinBoundaryCondition != "pin" && proteinBoundaryCondition != "none") {
    mem3dg_runtime_error("Invalid option for proteinBoundaryCondition!");
  }
}

void Parameters::Variation::checkParameters() {
  if (geodesicMask <= 0 && geodesicMask != -1) {
    mem3dg_runtime_error("Radius > 0 or radius = -1 to disable!");
  }
}

void Parameters::checkParameters(bool hasBoundary, size_t nVertex) {
  bending.checkParameters();
  variation.checkParameters();

  // variation
  if (!variation.isShapeVariation) {
    if (tension.form != NULL) {
      mem3dg_runtime_warning("tension.form is not NULL!");
    }
    if (osmotic.form != NULL) {
      mem3dg_runtime_warning("osmotic.form is not NULL!");
    }
    if (boundary.shapeBoundaryCondition != "none") {
      mem3dg_runtime_warning("shape boundary condition is not none for non "
                             "shape variation simulation");
    }
  }

  if (variation.isProteinVariation != (proteinMobility > 0)) {
    mem3dg_runtime_warning("proteinMobility value is not consistent with the "
                           "protein variation option!");
  }

  if (variation.isProteinConservation) {
    if (adsorption.epsilon != 0)
      mem3dg_runtime_warning(
          "protein adsorption has no effect when conserve protein!")
  }

  // boundary
  if (hasBoundary) {
    if (boundary.shapeBoundaryCondition == "none" &&
        variation.isShapeVariation) {
      mem3dg_runtime_warning(
          "Shape boundary condition type (roller, pin or fixed) "
          "has not been specified for open boundary mesh! May result in "
          "unexpected behavior(e.g. osmotic force). ");
    }
    if (boundary.proteinBoundaryCondition != "pin" &&
        variation.isProteinVariation) {
      mem3dg_runtime_warning("Protein boundary condition type (pin) "
                             "has not been specified for open boundary mesh!");
    }
  } else {
    if (boundary.shapeBoundaryCondition != "none") {
      mem3dg_runtime_warning(
          "Shape boundary condition type should be disable (= \"none\") "
          "for closed boundary mesh!");
    }
    if (boundary.proteinBoundaryCondition != "none") {
      mem3dg_runtime_warning(
          "Protein boundary condition type should be disable (= \"none\") "
          "for closed boundary mesh!");
    }
  }
}

} // namespace solver
} // namespace mem3dg
