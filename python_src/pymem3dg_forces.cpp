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
#include <cstdarg>
#include <cstddef>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Eigen/src/Core/util/Constants.h"

#include "mem3dg/solver/geometry.h"
#include "mem3dg/solver/integrator/conjugate_gradient.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/version.h"

#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include "mem3dg/mem3dg"
#include "pybind11/cast.h"

namespace gc = ::geometrycentral;
namespace mem3dg {
namespace solver {
namespace integrator {
namespace py = pybind11;

void init_forces(py::module_ &pymem3dg) {
  // ==========================================================
  // =============     Forces             ===============
  // ==========================================================
  py::class_<Forces> forces(pymem3dg, "Forces",
                            R"delim(
        The forces object
    )delim");
  forces.def(
      "getSurfaceTension", [](Forces &s) { return s.surfaceTension; },
      py::return_value_policy::copy,
      R"delim(
          get the Surface tension
      )delim");
  forces.def(
      "getOsmoticPressure", [](Forces &s) { return s.osmoticPressure; },
      py::return_value_policy::copy,
      R"delim(
          get the osmotic pressure
      )delim");

  /**
   * @brief Mechanical force
   */
  forces.def(
      "getDeviatoricCurvatureForceVec",
      [](Forces &s) { return toMatrix(s.deviatoricCurvatureForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the deviatoric curvature force of the system
      )delim");
  forces.def(
      "getAreaDifferenceForceVec",
      [](Forces &s) { return toMatrix(s.areaDifferenceForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the area difference force of the system
      )delim");
  forces.def(
      "getSpontaneousCurvatureForceVec",
      [](Forces &s) { return toMatrix(s.spontaneousCurvatureForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the spontaneous curvature force of the system
      )delim");
  forces.def(
      "getSpontaneousCurvatureForceVec_areaGrad",
      [](Forces &s) {
        return toMatrix(s.spontaneousCurvatureForceVec_areaGrad);
      },
      py::return_value_policy::copy,
      R"delim(
          get the area gradient component of the spontaneous curvature force of the system
      )delim");
  forces.def(
      "getSpontaneousCurvatureForceVec_gaussVec",
      [](Forces &s) {
        return toMatrix(s.spontaneousCurvatureForceVec_gaussVec);
      },
      py::return_value_policy::copy,
      R"delim(
          get the the gaussian curvature vector component of the spontaneous curvature force of the system
      )delim");
  forces.def(
      "getSpontaneousCurvatureForceVec_schlafliVec",
      [](Forces &s) {
        return toMatrix(s.spontaneousCurvatureForceVec_schlafliVec);
      },
      py::return_value_policy::copy,
      R"delim(
          get the Schlaflic (smoothing) component of the spontaneous curvature force of the system
      )delim");
  forces.def(
      "getCapillaryForceVec",
      [](Forces &s) { return toMatrix(s.capillaryForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the tension-induced capillary Force
      )delim");
  forces.def(
      "getLineCapillaryForceVec",
      [](Forces &s) { return toMatrix(s.lineCapillaryForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the interfacial line tension
      )delim");
  forces.def(
      "getExternalForceVec",
      [](Forces &s) { return toMatrix(s.externalForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the externally-applied Force
      )delim");
  forces.def(
      "getOsmoticForceVec",
      [](Forces &s) { return toMatrix(s.osmoticForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the osmotic force
      )delim");
  forces.def(
      "getAdsorptionForceVec",
      [](Forces &s) { return toMatrix(s.adsorptionForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the adsorption force
      )delim");
  forces.def(
      "getAggregationForceVec",
      [](Forces &s) { return toMatrix(s.aggregationForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the aggregation force
      )delim");
  forces.def(
      "getEntropyForceVec",
      [](Forces &s) { return toMatrix(s.entropyForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the entropy force
      )delim");
  forces.def(
      "getSpringForceVec", [](Forces &s) { return toMatrix(s.springForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the spring force
      )delim");
  forces.def(
      "getMechanicalForceVec",
      [](Forces &s) { return toMatrix(s.mechanicalForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the the total mechanical force
      )delim");

  /**
   * @brief Chemical Potential
   */
  forces.def(
      "getSpontaneousCurvaturePotential",
      [](Forces &s) { return s.spontaneousCurvaturePotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the spontaneous curvature potential
      )delim");
  forces.def(
      "getDeviatoricCurvaturePotential",
      [](Forces &s) { return s.deviatoricCurvaturePotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the deviatoric curvature potential
      )delim");
  forces.def(
      "getInteriorPenaltyPotential",
      [](Forces &s) { return s.interiorPenaltyPotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the interior point potential
      )delim");
  forces.def(
      "getAdsorptionPotential",
      [](Forces &s) { return s.adsorptionPotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the adsorption potential
      )delim");
  forces.def(
      "getAggregationPotential",
      [](Forces &s) { return s.aggregationPotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the aggregation potential
      )delim");
  forces.def(
      "getDirichletPotential",
      [](Forces &s) { return s.dirichletPotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the dirichlet Potential
      )delim");
  forces.def(
      "getEntropyPotential", [](Forces &s) { return s.entropyPotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the entropy Potential
      )delim");
  forces.def(
      "getChemicalPotential",
      [](Forces &s) { return s.chemicalPotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the chemical Potential
      )delim");
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
