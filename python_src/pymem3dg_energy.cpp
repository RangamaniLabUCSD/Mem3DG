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

void init_energy(py::module_ &pymem3dg) {
  // ==========================================================
  // =============   Free energy                ===============
  // ==========================================================
  py::class_<Energy> energy(pymem3dg, "Energy", R"delim(
        The energy
    )delim");
  energy.def(py::init<>());
  energy.def_readonly("totalEnergy", &Energy::totalEnergy,
                      R"delim(
          get total Energy of the system
      )delim");
  energy.def_readonly("kineticEnergy", &Energy::kineticEnergy,
                      R"delim(
          get kinetic energy of the membrane
      )delim");
  energy.def_readonly("potentialEnergy", &Energy::potentialEnergy,
                      R"delim(
          get potential energy of the membrane
      )delim");
  energy.def_readonly("spontaneousCurvatureEnergy",
                      &Energy::spontaneousCurvatureEnergy,
                      R"delim(
          get bending energy of the membrane
      )delim");
  energy.def_readonly("deviatoricCurvatureEnergy",
                      &Energy::deviatoricCurvatureEnergy,
                      R"delim(
          get deviatoric energy of the membrane
      )delim");
  energy.def_readonly("areaDifferenceEnergy", &Energy::areaDifferenceEnergy,
                      R"delim(
          get area difference energy of the membrane
      )delim");
  energy.def_readonly("surfaceEnergy", &Energy::surfaceEnergy,
                      R"delim(
          get stretching energy of the membrane
      )delim");
  energy.def_readonly("pressureEnergy", &Energy::pressureEnergy,
                      R"delim(
          get work of pressure within membrane
      )delim");
  energy.def_readonly("adsorptionEnergy", &Energy::adsorptionEnergy,
                      R"delim(
          get adsorption energy of the membrane protein
      )delim");
  energy.def_readonly("aggregationEnergy", &Energy::aggregationEnergy,
                      R"delim(
          get aggregation energy of the membrane protein
      )delim");
  energy.def_readonly("entropyEnergy", &Energy::entropyEnergy,
                      R"delim(
          get entropy energy of the membrane protein
      )delim");
  energy.def_readonly("edgeSpringEnergy", &Energy::edgeSpringEnergy,
                      R"delim(
          get edgeSpring energy of the membrane
      )delim");
  energy.def_readonly("faceSpringEnergy", &Energy::faceSpringEnergy,
                      R"delim(
          get faceSpring energy of the membrane
      )delim");
  energy.def_readonly("lcrSpringEnergy", &Energy::lcrSpringEnergy,
                      R"delim(
          get lcrSpring energy of the membrane
      )delim");
  energy.def_readonly("dirichletEnergy", &Energy::dirichletEnergy,
                      R"delim(
          get  line tension (dirichlet) energy of interface energy
      )delim");
  energy.def_readonly("externalWork", &Energy::externalWork,
                      R"delim(
          get work of external force
      )delim");
  energy.def_readonly("proteinInteriorPenalty", &Energy::proteinInteriorPenalty,
                      R"delim(
          get protein interior penalty energy (numerical energy)
      )delim");
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
