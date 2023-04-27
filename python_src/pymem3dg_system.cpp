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

void init_system(py::module_ &pymem3dg) {
  // ==========================================================
  // =============          System              ===============
  // ==========================================================
  py::class_<System> system(pymem3dg, "System",
                            R"delim(
        The system
    )delim");

  /**
   * @brief Constructors
   */
  system.def(py::init<Geometry &, EigenVectorX1d &, EigenVectorX3dr &,
                      Parameters &, double>(),
             py::arg("geometry"), py::arg("proteinDensity"),
             py::arg("velocity"), py::arg("parameters"), py::arg("time") = 0,
             R"delim(
        System constructor with Matrices
      )delim");
  system.def(py::init<Geometry &, Parameters &, double>(), py::arg("geometry"),
             py::arg("parameters"), py::arg("time") = 0,
             R"delim(
        System constructor with Matrices
      )delim");
#ifdef MEM3DG_WITH_NETCDF
  system.def(py::init<Geometry &, std::string, int, Parameters &>(),
             py::arg("geometry"), py::arg("trajFile"), py::arg("startingFrame"),
             py::arg("parameters"),
             R"delim(
        System constructor with NetCDF trajectory file
      )delim");
#endif

  /**
   * @brief Read/write attribute
   */
  system.def_readwrite("parameters", &System::parameters,
                       "Get the Parameters struct");
  system.def_readwrite("meshProcessor", &System::meshProcessor,
                       "Get the mesh processor object");
  system.def_readwrite("time", &System::time, "Get the time");
  /**
   * @brief Getters
   */
  system.def(
      "getGeometry", [](System &s) -> Geometry & { return s.geometry; },
      py::return_value_policy::automatic_reference,
      R"delim(
            Get the geometry

            Args:
                System (:py:class:`System`): System of interest

            Returns:
                :py:class:`pymem3dg.Geometry`: Geometry
      )delim");
  system.def(
      "getForces", [](System &s) { return s.forces; },
      py::return_value_policy::copy,
      R"delim(
            Get the pointwise forces

            Args:
                System (:py:class:`System`): System of interest

            Returns:
                :py:class:`pymem3dg.Forces`: Forces
      )delim");
  system.def(
      "getEnergy", [](System &s) { return s.energy; },
      py::return_value_policy::copy,
      R"delim(
            Get the energy

            Args:
                System (:py:class:`System`): System of interest

            Returns:
                :py:class:`float`: Energy
      )delim");
  system.def(
      "getSpontaneousCurvature", [](System &s) { return s.H0.raw(); },
      py::return_value_policy::copy, "get the pointwise spontaneous curvature");
  system.def(
      "getVelocity",
      [](System &s) { return gc::EigenMap<double, 3>(s.velocity); },
      py::return_value_policy::copy, "get the vertex velocity matrix");
  system.def(
      "getProteinDensity", [](System &s) { return s.proteinDensity.raw(); },
      py::return_value_policy::copy, "get the protein Density");
  system.def(
      "getProteinRateOfChange",
      [](System &s) { return s.proteinRateOfChange.raw(); },
      py::return_value_policy::copy, "get the protein rate of change");

  /**
   * @brief Method: force computation
   */
  system.def("computeInPlaneFluxForm", &System::computeInPlaneFluxForm,
             py::arg("chemicalPotential"),
             R"delim(
            Compute in plane flux form from chemical potential

            Args:
                chemicalPotential (npt.NDArray[np.float64]): components (or sum) of chemical potential

            Returns
                np.NDArray[np.float64]: in plane flux form on edges
        )delim");
  system.def("computeConservativeForcing", &System::computeConservativeForcing,
             R"delim(
            Compute and update all conservative forces, update mechanicalForce(Vec) with conservativeForce(Vec)
        )delim");
  system.def("addNonconservativeForcing", &System::addNonconservativeForcing,
             py::arg("timeStep") = 0,
             R"delim(
            Compute and append all non-conservative forces, update mechanicalForce(Vec) and mechErrorNorm
        )delim");
  system.def("prescribeExternalForce", &System::prescribeExternalForce,
             py::return_value_policy::copy,
             R"delim(
              prescribe the External Force
          )delim");

  /**
   * @brief Method: Energy computation
   */
  system.def("computeTotalEnergy", &System::computeTotalEnergy,
             R"delim(
          compute the total energy, where total energy = kinetic energy + potential energy - external work
      )delim");
  system.def(
      "computeIntegratedPower",
      static_cast<double (System::*)(double)>(&System::computeIntegratedPower),
      py::arg("dt"),
      R"delim(
            Intermediate function to integrate the power
        )delim");

  /**
   * @brief Method: initialize System
   */
  system.def("initialize", &System::initialize, py::arg("ifMutateMesh") = false,
             R"delim(
          initialize the system
      )delim");

  /**
   * @brief Method: updateVertexPosition
   */
  system.def("updateConfigurations", &System::updateConfigurations,
             R"delim(
          update the system configuration due to changes in state variables (e.g vertex positions or protein density)
      )delim");

  /**
   * @brief Method: I/O
   */
  system.def("saveRichData", &System::saveRichData, py::arg("pathToSave"),
             py::arg("isJustGeometry") = false,
             R"delim(
          save snapshot data to directory
      )delim");

  /**
   * @brief Method: mutate the mesh
   */
  system.def("mutateMesh", &System::mutateMesh, py::arg("nMutation") = 1,
             R"delim(
          mutate the mesh
      )delim");

  /**
   * @brief Method: test force computation
   */
  system.def("testConservativeForcing", &System::testConservativeForcing,
             py::arg("timeStep"),
             R"delim(
          test conservative force computation by validating energy decrease
      )delim");

  /**
   * @brief Method: smoothen the mesh
   */
  system.def("smoothenMesh", &System::smoothenMesh, py::arg("initStep"),
             py::arg("target"), py::arg("maxIteration"),
             R"delim(
          smoothen the mesh using spontaneous curvature force
      )delim");
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
