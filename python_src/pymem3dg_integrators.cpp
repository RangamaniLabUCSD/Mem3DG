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

void init_integrators(py::module_ &pymem3dg) {

  // ==========================================================
  // =============     Integrator Template      ===============
  // ==========================================================
  //   py::class_<Integrator> integrator(pymem3dg, "Integrator", R"delim(
  //         The integrator
  //     )delim");
  //   integrator.def(py::init<System &, double, double, double, double,
  //   std::string,
  //                           bool, std::string, std::size_t>(),
  //                  py::arg("f"), py::arg("dt"), py::arg("total_time"),
  //                  py::arg("tSave"), py::arg("tolerance"),
  //                  py::arg("outputDir"), py::arg("ifAdaptiveStep") = true,
  //                  py::arg("trajFileName") = "traj.nc", py::arg("verbosity")
  //                  = 3, R"delim(
  //         Integrator constructor
  //       )delim");
  //   integrator.def("saveData", &Integrator::saveData,
  //                  R"delim(
  //           save data to output directory
  //       )delim");
  //   integrator.def("saveRichData", &Integrator::saveRichData,
  //   py::arg("plyName"),
  //                  "save to richData and output .ply file to output
  //                  directory", R"delim(
  //           save data to output directory
  //       )delim");
  // #ifdef MEM3DG_WITH_NETCDF
  //   integrator.def("saveNetcdfData", &Integrator::saveNetcdfData,
  //                  "save to netcdf file in output directory",
  //                  R"delim(
  //           save data to output directory
  //       )delim");
  // #endif

  // ==========================================================
  // =============     Velocity Verlet          ===============
  // ==========================================================
  py::class_<VelocityVerlet> velocityverlet(pymem3dg, "VelocityVerlet",
                                            R"delim(
        Velocity Verlet integration
    )delim");

  velocityverlet.def(py::init<System &, double, double, double, double,
                              std::string, std::size_t>(),
                     py::arg("system"), py::arg("characteristicTimeStep"),
                     py::arg("totalTime"), py::arg("savePeriod"),
                     py::arg("tolerance"), py::arg("outputDirectory"),
                     py::arg("frame") = 0,
                     R"delim(
        Velocity Verlet integrator constructor
      )delim");

  velocityverlet.def_readwrite("isBacktrack", &VelocityVerlet::isBacktrack,
                               R"delim(
         whether do backtracking line search
      )delim");
  velocityverlet.def_readwrite("rho", &VelocityVerlet::rho,
                               R"delim(
          backtracking coefficient
      )delim");
  velocityverlet.def_readwrite("c1", &VelocityVerlet::c1,
                               R"delim(
          Wolfe condition parameter
      )delim");
  velocityverlet.def_readwrite("ifOutputTrajFile",
                               &VelocityVerlet::ifOutputTrajFile, R"delim(
          if output trajectory file
      )delim");
  velocityverlet.def_readwrite("ifOutputMeshFile",
                               &VelocityVerlet::ifOutputMeshFile, R"delim(
          if output mesh file
      )delim");
  velocityverlet.def_readwrite("ifPrintToConsole",
                               &VelocityVerlet::ifPrintToConsole, R"delim(
          if print to console
      )delim");
  velocityverlet.def_readwrite("trajFileName", &VelocityVerlet::trajFileName,
                               R"delim(
          name of the trajectory file
      )delim");
  velocityverlet.def_readwrite("ifAdaptiveStep",
                               &VelocityVerlet::ifAdaptiveStep,
                               R"delim(
          option to scale time step according to mesh size
      )delim");
  velocityverlet.def_readwrite("isCapEnergy", &VelocityVerlet::isCapEnergy,
                               R"delim(
          option to exit if exceed initial energy cap
      )delim");
  velocityverlet.def_readwrite("outputDirectory",
                               &VelocityVerlet::outputDirectory,
                               R"delim(
        collapse small triangles
      )delim");

  velocityverlet.def("integrate", &VelocityVerlet::integrate,
                     R"delim(
          integrate
      )delim");
  velocityverlet.def("status", &VelocityVerlet::status,
                     R"delim(
          status computation and thresholding
      )delim");
  velocityverlet.def("march", &VelocityVerlet::march,
                     R"delim(
          stepping forward
      )delim");
  velocityverlet.def("saveData", &VelocityVerlet::saveData,
                     py::arg("ifOutputTrajFile"), py::arg("ifOutputMeshFile"),
                     py::arg("ifPrintToConsole"),
                     R"delim(
          save data to output directory
      )delim");
  velocityverlet.def("step", &VelocityVerlet::step, py::arg("n"),
                     R"delim(
          step for n iterations
      )delim");

  // ==========================================================
  // =============     Forward Euelr            ===============
  // ==========================================================
  py::class_<Euler> euler(pymem3dg, "Euler", R"delim(
        forward euler (gradient descent) integration
    )delim");

  euler.def(py::init<System &, double, double, double, double, std::string,
                     std::size_t>(),
            py::arg("system"), py::arg("characteristicTimeStep"),
            py::arg("totalTime"), py::arg("savePeriod"), py::arg("tolerance"),
            py::arg("outputDirectory"), py::arg("frame") = 0,
            R"delim(
        Euler integrator (steepest descent) constructor
      )delim");
  euler.def(py::init<System &, double, double, std::string>(),
            py::arg("system"), py::arg("characteristicTimeStep"),
            py::arg("tolerance"), py::arg("outputDirectory"),
            R"delim(
        Euler integrator (steepest descent) constructor
      )delim");
  /**
   * @brief state variable
   */
  euler.def_readonly("EXIT", &Euler::EXIT,
                     R"delim(
          Flag for terminating the simulation
      )delim");
  euler.def_readonly("timeStep", &Euler::timeStep,
                     R"delim(
          time step
      )delim");

  /**
   * @brief attributes, integration optionsP
   */
  euler.def_readonly("characteristicTimeStep", &Euler::characteristicTimeStep,
                     R"delim(
          characteristic time step
      )delim");
  euler.def_readonly("totalTime", &Euler::totalTime,
                     R"delim(
          time limit
      )delim");
  euler.def_readonly("savePeriod", &Euler::savePeriod,
                     R"delim(
         period of saving output data
      )delim");
  euler.def_readonly("tolerance", &Euler::tolerance,
                     R"delim(
          tolerance for termination
      )delim");
  euler.def_readwrite("ifOutputTrajFile", &Euler::ifOutputTrajFile, R"delim(
          if output trajectory file
      )delim");
  euler.def_readwrite("ifOutputMeshFile", &Euler::ifOutputMeshFile, R"delim(
          if output mesh file
      )delim");
  euler.def_readwrite("ifPrintToConsole", &Euler::ifPrintToConsole, R"delim(
          if print to console
      )delim");
  euler.def_readwrite("trajFileName", &Euler::trajFileName,
                      R"delim(
          name of the trajectory file
      )delim");
  euler.def_readwrite("ifAdaptiveStep", &Euler::ifAdaptiveStep,
                      R"delim(
          option to scale time step according to mesh size
      )delim");
  euler.def_readwrite("outputDirectory", &Euler::outputDirectory,
                      R"delim(
        collapse small triangles
      )delim");
  euler.def_readwrite("ifJustGeometryPly", &Euler::ifJustGeometryPly,
                      R"delim(
           save .ply with just geometry
      )delim");
  euler.def_readwrite("isBacktrack", &Euler::isBacktrack,
                      R"delim(
         whether do backtracking line search
      )delim");
  euler.def_readwrite("rho", &Euler::rho,
                      R"delim(
          backtracking coefficient
      )delim");
  euler.def_readwrite("c1", &Euler::c1,
                      R"delim(
          Wolfe condition parameter
      )delim");

  /**
   * @brief methods
   */
  euler.def(
      "integrate",
      [](Euler &self) {
        py::scoped_ostream_redirect stream(
            std::cout,                                // std::ostream&
            py::module_::import("sys").attr("stdout") // Python output
        );
        self.integrate();
      },
      R"delim(
            integrate
        )delim");

  euler.def("status", &Euler::status,
            R"delim(
          status computation and thresholding
      )delim");
  euler.def("march", &Euler::march,
            R"delim(
          stepping forward
      )delim");
  euler.def("saveData", &Euler::saveData, py::arg("ifOutputTrajFile"),
            py::arg("ifOutputMeshFile"), py::arg("ifPrintToConsole"),
            R"delim(
          save data to output directory
      )delim");
  euler.def("step", &Euler::step, py::arg("n"),
            R"delim(
          step for n iterations
      )delim");

#ifdef MEM3DG_WITH_NETCDF
  euler.def("createMutableNetcdfFile", &Euler::createMutableNetcdfFile,
            py::arg("isContinue"),
            R"delim(
          create netcdf file
      )delim");
  euler.def("saveMutableNetcdfData", &Euler::saveMutableNetcdfData,
            R"delim(
          write to netcdf file
      )delim");
  euler.def("closeMutableNetcdfFile", &Euler::closeMutableNetcdfFile,
            R"delim(
          close netcdf file
      )delim");

#endif

  // ==========================================================
  // =============     Conjugate Gradient       ===============
  // ==========================================================
  py::class_<ConjugateGradient> conjugategradient(pymem3dg, "ConjugateGradient",
                                                  R"delim(
        conjugate Gradient propagator
    )delim");

  conjugategradient.def(py::init<System &, double, double, double, double,
                                 std::string, std::size_t>(),
                        py::arg("system"), py::arg("characteristicTimeStep"),
                        py::arg("totalTime"), py::arg("savePeriod"),
                        py::arg("tolerance"), py::arg("outputDirectory"),
                        py::arg("frame") = 0,
                        R"delim(
        Conjugate Gradient optimizer constructor
      )delim");

  /**
   * @brief attributes, integration options
   */
  conjugategradient.def_readonly("characteristicTimeStep",
                                 &ConjugateGradient::characteristicTimeStep,
                                 R"delim(
          characteristic time step
      )delim");
  conjugategradient.def_readonly("totalTime", &ConjugateGradient::totalTime,
                                 R"delim(
          time limit
      )delim");
  conjugategradient.def_readonly("savePeriod", &ConjugateGradient::savePeriod,
                                 R"delim(
         period of saving output data
      )delim");
  conjugategradient.def_readonly("tolerance", &ConjugateGradient::tolerance,
                                 R"delim(
          tolerance for termination
      )delim");
  conjugategradient.def_readwrite("ifOutputTrajFile",
                                  &ConjugateGradient::ifOutputTrajFile, R"delim(
          if output trajectory file
      )delim");
  conjugategradient.def_readwrite("ifOutputMeshFile",
                                  &ConjugateGradient::ifOutputMeshFile, R"delim(
          if output mesh file
      )delim");
  conjugategradient.def_readwrite("ifPrintToConsole",
                                  &ConjugateGradient::ifPrintToConsole, R"delim(
          if print to console
      )delim");
  conjugategradient.def_readwrite("trajFileName",
                                  &ConjugateGradient::trajFileName,
                                  R"delim(
          name of the trajectory file
      )delim");
  conjugategradient.def_readwrite("ifAdaptiveStep",
                                  &ConjugateGradient::ifAdaptiveStep,
                                  R"delim(
          option to scale time step according to mesh size
      )delim");
  conjugategradient.def_readwrite("outputDirectory",
                                  &ConjugateGradient::outputDirectory,
                                  R"delim(
        collapse small triangles
      )delim");
  conjugategradient.def_readwrite("ifJustGeometryPly",
                                  &ConjugateGradient::ifJustGeometryPly,
                                  R"delim(
           save .ply with just geometry
      )delim");
  conjugategradient.def_readwrite("isBacktrack",
                                  &ConjugateGradient::isBacktrack,
                                  R"delim(
         whether do backtracking line search
      )delim");
  conjugategradient.def_readwrite("rho", &ConjugateGradient::rho,
                                  R"delim(
          backtracking coefficient
      )delim");
  conjugategradient.def_readwrite("c1", &ConjugateGradient::c1,
                                  R"delim(
          Wolfe condition parameter
      )delim");
  conjugategradient.def_readwrite("restartPeriod",
                                  &ConjugateGradient::restartPeriod,
                                  R"delim(
          option to restart conjugate gradient using gradient descent
      )delim");

  conjugategradient.def_readwrite("constraintTolerance",
                                  &ConjugateGradient::constraintTolerance,
                                  R"delim(
            tolerance for constraints
      )delim");
  conjugategradient.def_readwrite("isAugementedLagrangian",
                                  &ConjugateGradient::isAugmentedLagrangian,
                                  R"delim(
            whether use augmented lagrangian method
      )delim");

  /**
   * @brief methods
   */
  conjugategradient.def("integrate", &ConjugateGradient::integrate,
                        R"delim(
          integrate
      )delim");
  conjugategradient.def("status", &ConjugateGradient::status,
                        R"delim(
          status computation and thresholding
      )delim");
  conjugategradient.def("march", &ConjugateGradient::march,
                        R"delim(
          stepping forward
      )delim");
  conjugategradient.def(
      "saveData", &ConjugateGradient::saveData, py::arg("ifOutputTrajFile"),
      py::arg("ifOutputMeshFile"), py::arg("ifPrintToConsole"),
      R"delim(
          save data to output directory
      )delim");
  conjugategradient.def("step", &ConjugateGradient::step, py::arg("n"),
                        R"delim(
          step for n iterations
      )delim");
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
