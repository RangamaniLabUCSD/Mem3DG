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

#include "mem3dg/solver/integrator/conjugate_gradient.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/solver/system.h"

#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>

#include "mem3dg/mem3dg"
#include "pybind11/cast.h"

namespace gc = ::geometrycentral;
namespace mem3dg {
namespace solver {
namespace integrator {
namespace py = pybind11;

// Initialize the `pymem3dg` module
PYBIND11_MODULE(_core, pymem3dg) {
  pymem3dg.doc() = "Python wrapper around the DDG solver C++ library.";

#pragma region integrators
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
  velocityverlet.def_readwrite("updateGeodesicsPeriod",
                               &VelocityVerlet::updateGeodesicsPeriod,
                               R"delim(
          period of update geodesics
      )delim");
  velocityverlet.def_readwrite("processMeshPeriod",
                               &VelocityVerlet::processMeshPeriod,
                               R"delim(
          period of processing mesh
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
  euler.def_readwrite("updateGeodesicsPeriod", &Euler::updateGeodesicsPeriod,
                      R"delim(
          period of update geodesics
      )delim");
  euler.def_readwrite("processMeshPeriod", &Euler::processMeshPeriod,
                      R"delim(
          period of processing mesh
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
  euler.def("integrate", &Euler::integrate,
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
  conjugategradient.def_readwrite(
      "ifOutputTrajFile", &ConjugateGradient::ifOutputTrajFile, R"delim(
          if output trajectory file
      )delim");
  conjugategradient.def_readwrite(
      "ifOutputMeshFile", &ConjugateGradient::ifOutputMeshFile, R"delim(
          if output mesh file
      )delim");
  conjugategradient.def_readwrite(
      "ifPrintToConsole", &ConjugateGradient::ifPrintToConsole, R"delim(
          if print to console
      )delim");
  conjugategradient.def_readwrite("updateGeodesicsPeriod",
                                  &ConjugateGradient::updateGeodesicsPeriod,
                                  R"delim(
          period of update geodesics
      )delim");
  conjugategradient.def_readwrite("processMeshPeriod",
                                  &ConjugateGradient::processMeshPeriod,
                                  R"delim(
          period of processing mesh
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

#pragma endregion integrators

#pragma region forces
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
      "getDeviatoricForceVec",
      [](Forces &s) { return toMatrix(s.deviatoricForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the deviatoric force of the system
      )delim");
  forces.def(
      "getBendingForceVec",
      [](Forces &s) { return toMatrix(s.bendingForceVec); },
      py::return_value_policy::copy,
      R"delim(
          get the bending force of the system
      )delim");
  forces.def(
      "getBendingForceVec_areaGrad",
      [](Forces &s) { return toMatrix(s.bendingForceVec_areaGrad); },
      py::return_value_policy::copy,
      R"delim(
          get the area gradient component of the bending force of the system
      )delim");
  forces.def(
      "getBendingForceVec_gaussVec",
      [](Forces &s) { return toMatrix(s.bendingForceVec_gaussVec); },
      py::return_value_policy::copy,
      R"delim(
          get the the gaussian curvature vector component of the bending force of the system
      )delim");
  forces.def(
      "getBendingForceVec_schlafliVec",
      [](Forces &s) { return toMatrix(s.bendingForceVec_schlafliVec); },
      py::return_value_policy::copy,
      R"delim(
          get the Schlaflic (smoothing) component of the bending force of the system
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
      "getBendingPotential", [](Forces &s) { return s.bendingPotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the bending potential
      )delim");
  forces.def(
      "getDeviatoricPotential",
      [](Forces &s) { return s.deviatoricPotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the deviatoric potential
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
      "getInPlaneFluxForm",
      [](Forces &s) { return s.inPlaneFluxForm.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the inPlaneFluxForm on edge (one-form)
      )delim");
  forces.def(
      "getChemicalPotential",
      [](Forces &s) { return s.chemicalPotential.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the chemical Potential
      )delim");

#pragma endregion forces

#pragma region mesh_mutator
  // ==========================================================
  // =============     MeshMutator              ===============
  // ==========================================================
  py::class_<MeshProcessor::MeshMutator> meshmutator(pymem3dg, "MeshMutator",
                                                     R"delim(
        The mesh mutator settings 
    )delim");
  meshmutator.def(py::init<>(),
                  R"delim(
       meshmutator constructor
      )delim");
  meshmutator.def_readonly("isFlipEdge",
                           &MeshProcessor::MeshMutator::isFlipEdge,
                           R"delim(
          get the option of whether do edge flip
      )delim");
  meshmutator.def_readonly("isSplitEdge",
                           &MeshProcessor::MeshMutator::isSplitEdge,
                           R"delim(
          get the option of whether split edge to grow mesh
      )delim");
  meshmutator.def_readonly("isCollapseEdge",
                           &MeshProcessor::MeshMutator::isCollapseEdge,
                           R"delim(
          get the option of whether Collapse edge to grow mesh
      )delim");
  meshmutator.def_readonly("isChangeTopology",
                           &MeshProcessor::MeshMutator::isChangeTopology,
                           R"delim(
          get the option of change topology
      )delim");

  /**
   * @brief vertex shifting
   */
  meshmutator.def_readwrite("isShiftVertex",
                            &MeshProcessor::MeshMutator::isShiftVertex,
                            R"delim(
          get the option of whether do vertex shift  
      )delim");

  /**
   * @brief mesh smoothing
   */
  meshmutator.def_readwrite("isSmoothenMesh",
                            &MeshProcessor::MeshMutator::isSmoothenMesh,
                            R"delim(
          whether conduct mesh smoothing operation
      )delim");

  /**
   * @brief flipping criterion
   */
  meshmutator.def_readwrite("flipNonDelaunay",
                            &MeshProcessor::MeshMutator::flipNonDelaunay,
                            R"delim(
          whether flip non-Delaunay edge
      )delim");
  meshmutator.def_readwrite(
      "flipNonDelaunayRequireFlat",
      &MeshProcessor::MeshMutator::flipNonDelaunayRequireFlat,
      R"delim(
          whether require flatness condition when flipping non-Delaunay edge
      )delim");

  /**
   * @brief splitting criterion
   */
  meshmutator.def_readwrite("splitLarge",
                            &MeshProcessor::MeshMutator::splitLarge,
                            R"delim(
          split edge with large faces
      )delim");
  meshmutator.def_readwrite("splitLong", &MeshProcessor::MeshMutator::splitLong,
                            R"delim(
          split long edge
      )delim");
  meshmutator.def_readwrite("splitCurved",
                            &MeshProcessor::MeshMutator::splitCurved,
                            R"delim(
          split edge on high curvature domain
      )delim");
  meshmutator.def_readwrite("splitSharp",
                            &MeshProcessor::MeshMutator::splitSharp,
                            R"delim(
          split edge with sharp membrane property change
      )delim");
  meshmutator.def_readwrite("splitFat", &MeshProcessor::MeshMutator::splitFat,
                            R"delim(
          split obtuse triangle
      )delim");
  meshmutator.def_readwrite("splitSkinnyDelaunay",
                            &MeshProcessor::MeshMutator::splitSkinnyDelaunay,
                            R"delim(
          split poor aspected triangle that is still Delaunay
      )delim");
  meshmutator.def_readwrite("minimumEdgeLength",
                            &MeshProcessor::MeshMutator::minimumEdgeLength,
                            R"delim(
          minimum edge length 
      )delim");

  /**
   * @brief collapsing criterion
   */
  meshmutator.def_readwrite("collapseSkinny",
                            &MeshProcessor::MeshMutator::collapseSkinny,
                            R"delim(
          collapse skinny triangles
      )delim");
  meshmutator.def_readwrite("collapseSmall",
                            &MeshProcessor::MeshMutator::collapseSmall,
                            R"delim(
          collapse small triangles
      )delim");
  meshmutator.def_readwrite("collapseFlat",
                            &MeshProcessor::MeshMutator::collapseFlat,
                            R"delim(
         collapse flat edge
      )delim");
  meshmutator.def_readwrite("curvTol", &MeshProcessor::MeshMutator::curvTol,
                            R"delim(
          tolerance for curvature approximation
      )delim");
  meshmutator.def_readwrite("targetFaceArea",
                            &MeshProcessor::MeshMutator::targetFaceArea,
                            R"delim(
          target face area 
      )delim");

  py::class_<MeshProcessor> meshprocessor(pymem3dg, "MeshProcessor",
                                          R"delim(
        The mesh processor settings 
    )delim");
  meshprocessor.def(py::init<>(),
                    R"delim(
       meshprocessor constructor
      )delim");
  meshprocessor.def_readwrite("meshMutator", &MeshProcessor::meshMutator,
                              R"delim(
          meshMutator struct
      )delim");
  meshprocessor.def_readonly("isMeshMutate", &MeshProcessor::isMeshMutate,
                             R"delim(
          get the option of whether do mesh mutation
      )delim");
#pragma endregion mesh_mutator

#pragma region system
  // ==========================================================
  // =============          System              ===============
  // ==========================================================
  py::class_<System> system(pymem3dg, "System",
                            R"delim(
        The system
    )delim");

  /**
   * @brief Constructors by .ply file
   */
  system.def(py::init<std::string, std::string, EigenVectorX1d &,
                      EigenVectorX3dr &, Parameters &, double>(),
             py::arg("inputMesh"), py::arg("referenceMesh"),
             py::arg("proteinDensity"), py::arg("velocity"),
             py::arg("parameters"), py::arg("time") = 0,
             R"delim(
        System constructor with .ply files. 
      )delim");
  system.def(py::init<std::string, EigenVectorX1d &, EigenVectorX3dr &,
                      Parameters &, double>(),
             py::arg("inputMesh"), py::arg("proteinDensity"),
             py::arg("velocity"), py::arg("parameters"), py::arg("time") = 0,
             R"delim(
        System constructor with .ply files. 
      )delim");
  system.def(py::init<std::string, Parameters &, double>(),
             py::arg("inputMesh"), py::arg("parameters"), py::arg("time") = 0,

             R"delim(
        System constructor with .ply files. 
      )delim");
  system.def(
      py::init<std::string, EigenVectorX1d &, EigenVectorX3dr &, double>(),
      py::arg("inputMesh"), py::arg("proteinDensity"), py::arg("velocity"),
      py::arg("time") = 0,
      R"delim(
        System constructor with .ply files. 
      )delim");
  system.def(py::init<std::string, double>(), py::arg("inputMesh"),
             py::arg("time") = 0,
             R"delim(
        System constructor with .ply files. 
      )delim");

  /**
   * @brief Constructors by matrices
   */
  system.def(
      py::init<EigenVectorX3sr &, EigenVectorX3dr &, EigenVectorX3dr &,
               EigenVectorX1d &, EigenVectorX3dr &, Parameters &, double>(),
      py::arg("topologyMatrix"), py::arg("vertexMatrix"),
      py::arg("referenceVertexMatrix"), py::arg("proteinDensity"),
      py::arg("velocity"), py::arg("parameters"), py::arg("time") = 0,
      R"delim(
        System constructor with Matrices 
      )delim");
  system.def(py::init<EigenVectorX3sr &, EigenVectorX3dr &, EigenVectorX1d &,
                      EigenVectorX3dr &, Parameters &, double>(),
             py::arg("topologyMatrix"), py::arg("vertexMatrix"),
             py::arg("proteinDensity"), py::arg("velocity"),
             py::arg("parameters"), py::arg("time") = 0,
             R"delim(
        System constructor with Matrices 
      )delim");
  system.def(
      py::init<EigenVectorX3sr &, EigenVectorX3dr &, Parameters &, double>(),
      py::arg("topologyMatrix"), py::arg("vertexMatrix"), py::arg("parameters"),
      py::arg("time") = 0,
      R"delim(
        System constructor with Matrices 
      )delim");
  system.def(py::init<EigenVectorX3sr &, EigenVectorX3dr &, EigenVectorX1d &,
                      EigenVectorX3dr &, double>(),
             py::arg("topologyMatrix"), py::arg("vertexMatrix"),
             py::arg("proteinDensity"), py::arg("velocity"),
             py::arg("time") = 0,
             R"delim(
        System constructor with Matrices 
      )delim");
  system.def(py::init<EigenVectorX3sr &, EigenVectorX3dr &, double>(),
             py::arg("topologyMatrix"), py::arg("vertexMatrix"),
             py::arg("time") = 0,
             R"delim(
        System constructor with Matrices 
      )delim");

  /**
   * @brief Constructors by NetCDF trajectory file
   */
#ifdef MEM3DG_WITH_NETCDF
  system.def(py::init<std::string, int, Parameters &>(), py::arg("trajFile"),
             py::arg("startingFrame"), py::arg("parameters"),
             R"delim(
        System constructor with NetCDF trajectory file
      )delim");
  system.def(py::init<std::string, int>(), py::arg("trajFile"),
             py::arg("startingFrame"),
             R"delim(
        System constructor with NetCDF trajectory file
      )delim");
#endif

  /**
   * @brief Initializing arguments
   */
  system.def_readwrite("parameters", &System::parameters,
                       R"delim(
          get the Parameters struct
      )delim");
  system.def_readwrite("meshProcessor", &System::meshProcessor,
                       R"delim(
          get the mesh processor object
      )delim");
  system.def_readwrite("time", &System::time,
                       R"delim(
          get the time
      )delim");

  /**
   * @brief    Geometric properties (Geometry central) getter
   */
  system.def_readonly("surfaceArea", &System::surfaceArea,
                      R"delim(
          get the surface area of the mesh
      )delim");
  system.def_readonly("volume", &System::volume,
                      R"delim(
          get the enclosed volume of the mesh
      )delim");
  system.def(
      "getLumpedMassMatrix",
      [](System &s) { return s.vpg->vertexLumpedMassMatrix; },
      py::return_value_policy::copy,
      R"delim(
          get the lumped mass matrix of the mesh
      )delim");
  system.def(
      "getCotanLaplacian", [](System &s) { return s.vpg->cotanLaplacian; },
      py::return_value_policy::copy,
      R"delim(
          get the Cotan Laplacian matrix of the mesh
      )delim");
  system.def(
      "getVertexNormals",
      [](System &s) { return toMatrix(s.vpg->vertexNormals); },
      py::return_value_policy::copy,
      R"delim(
          get angle-weighted normal on vertices
      )delim");
  system.def(
      "getInputVertexPositions",
      [](System &s) {
        return gc::EigenMap<double, 3>(s.vpg->inputVertexPositions);
      },
      py::return_value_policy::copy,
      R"delim(
          get the vertex position matrix
      )delim");
  system.def(
      "getFaceVertexMatrix",
      [](System &s) { return s.mesh->getFaceVertexMatrix<std::size_t>(); },
      py::return_value_policy::copy,
      R"delim(
          get the face vertex matrix
      )delim");
  system.def(
      "getVertexAdjacencyMatrix", [](System &s) { return s.vpg->d0; },
      py::return_value_policy::copy,
      R"delim(
          get the signed E-V vertex adjacency matrix, equivalent of d0 operator
      )delim");
  system.def(
      "getEdgeAdjacencyMatrix", [](System &s) { return s.vpg->d1; },
      py::return_value_policy::copy,
      R"delim(
          get the signed F-E edge adjacency matrix, equivalent of d1 operator
      )delim");
  system.def(
      "getVertexDualAreas",
      [](System &s) { return s.vpg->vertexDualAreas.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get vertex dual area
      )delim");
  system.def(
      "getVertexMeanCurvatures",
      [](System &s) {
        s.vpg->requireVertexMeanCurvatures();
        return s.vpg->vertexMeanCurvatures.raw();
      },
      py::return_value_policy::copy,
      R"delim(
          get the integrated scalar mean curvature
      )delim");
  system.def(
      "getVertexGaussianCurvatures",
      [](System &s) {
        s.vpg->requireVertexGaussianCurvatures();
        return s.vpg->vertexGaussianCurvatures.raw();
      },
      py::return_value_policy::copy,
      R"delim(
          get the integrated scalar Gaussian Curvature
      )delim");
  system.def(
      "getVertexGaussianCurvatureVectors",
      [](System &s) {
        auto vector = s.computeVertexGaussianCurvatureVectors();
        return toMatrix(vector);
      },
      py::return_value_policy::copy,
      R"delim(
          get the integrated vector Gaussian Curvature
      )delim");
  system.def(
      "getVertexVolumeVariationVectors",
      [](System &s) {
        auto vector = s.computeVertexVolumeVariationVectors();
        return toMatrix(vector);
      },
      py::return_value_policy::copy,
      R"delim(
          get the integrated vector Volume Variation (dual area)
      )delim");
  system.def(
      "getVertexMeanCurvatureVectors",
      [](System &s) {
        auto vector = s.computeVertexMeanCurvatureVectors();
        return toMatrix(vector);
      },
      py::return_value_policy::copy,
      R"delim(
          get the integrated vector Mean Curvature
      )delim");
  system.def(
      "getVertexSchlafliVectors",
      [](System &s) {
        auto vector = s.computeVertexSchlafliVectors();
        return toMatrix(vector);
      },
      py::return_value_policy::copy,
      R"delim(
          get the vertex Schlafli Vectors
      )delim");

  /**
   * @brief    geometry setter
   */
  system.def(
      "setInputVertexPositions",
      [](System &s, EigenVectorX3dr newGeo) {
        gc::EigenMap<double, 3>(s.vpg->inputVertexPositions) = newGeo;
      },
      R"delim(
          set the vertex position matrix
      )delim");

  /**
   * @brief Membrane dynamics properties (Mem3DG)
   */
  system.def(
      "getForces", [](System &s) { return s.forces; },
      py::return_value_policy::copy,
      R"delim(
          get the pointwise spontaneous curvature
      )delim");
  system.def(
      "getEnergy", [](System &s) { return s.energy; },
      py::return_value_policy::copy,
      R"delim(
          get the pointwise spontaneous curvature
      )delim");
  system.def(
      "getSpontaneousCurvature", [](System &s) { return s.H0.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the pointwise spontaneous curvature
      )delim");
  system.def(
      "getVelocity",
      [](System &s) { return gc::EigenMap<double, 3>(s.velocity); },
      py::return_value_policy::copy,
      R"delim(
          get the vertex velocity matrix
      )delim");
  system.def(
      "getProteinDensity", [](System &s) { return s.proteinDensity.raw(); },
      py::return_value_policy::copy,
      R"delim(
          get the protein Density
      )delim");

  /**
   * @brief Method: force computation
   */
  system.def(
      "computePhysicalForcing",
      static_cast<void (System::*)(double)>(&System::computePhysicalForcing),
      py::arg("timeStep") = 0,
      R"delim(
            compute all the forces
        )delim");
  //   system.def("computeBendingForce", &System::computeBendingForce,
  //              py::return_value_policy::copy,
  //              R"delim(
  //           compute the bending force
  //       )delim");
  //   system.def("computeChemicalPotential", &System::computeChemicalPotential,
  //              py::return_value_policy::copy,
  //              R"delim(
  //           compute the chemical potential
  //       )delim");
  //   system.def("computeCapillaryForce", &System::computeCapillaryForce,
  //              py::return_value_policy::copy,
  //              R"delim(
  //           compute the capillary force
  //       )delim");
  //   system.def("computeOsmoticForce", &System::computeOsmoticForce,
  //              py::return_value_policy::copy,
  //              R"delim(
  //           compute the osmotic force
  //       )delim");
  //   system.def("computeLineCapillaryForce",
  //   &System::computeLineCapillaryForce,
  //              py::return_value_policy::copy,
  //              R"delim(
  //           compute the LineTensionForce
  //       )delim");
  system.def("prescribeExternalForce", &System::prescribeExternalForce,
             py::return_value_policy::copy,
             R"delim(
            prescribe the External Force
        )delim");
  system.def("computeDPDForces", &System::computeDPDForces, py::arg("dt"),
             R"delim(
            compute the DPDForces
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
  //   system.def("computeL1Norm", &System::computeL1Norm,
  //              R"delim(
  //                    compute error norm
  //             Args:
  //                   force (:py:class:`list`): mesh vertex force
  //         )delim");

  /**
   * @brief Method: initialize System
   */
  system.def("initialize", &System::initialize, py::arg("nMutation") = 0,
             py::arg("ifMute") = false,
             R"delim(
          initialize the system
      )delim");

  /**
   * @brief Method: updateGeodesics
   */
  system.def("updateGeodesicsDistance", &System::updateGeodesicsDistance,
             R"delim(
          update the geodesics and related configurations
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
   * @brief Method: I/O
   */
  system.def(
      "getPolyscopePermutations",
      [](System &s) { return gcs::polyscopePermutations(*s.mesh); },
      R"delim(
          get polyscope permutation 
      )delim");

  /**
   * @brief Method: I/O
   */
  system.def(
      "getPolyscopeEdgeOrientations",
      [](System &s) { return gcs::polyscopeEdgeOrientations(*s.mesh).raw(); },
      R"delim(
          get polyscope edge orientation 
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
  system.def("testForceComputation",
             py::overload_cast<const double, const EigenVectorX3dr,
                               const EigenVectorX1d, const Energy>(
                 &System::testForceComputation),
             py::arg("timeStep"), py::arg("previousPosition"),
             py::arg("previousProteinDensity"), py::arg("previousEnergy"),
             R"delim(
          test force computation by validating energy decrease (reverse mode)
      )delim");
  system.def("testForceComputation",
             py::overload_cast<const double>(&System::testForceComputation),
             py::arg("timeStep"),
             R"delim(
          test force computation by validating energy decrease (forward mode)
      )delim");

  /**
   * @brief Method: smoothen the mesh
   */
  system.def("smoothenMesh", &System::smoothenMesh, py::arg("initStep"),
             py::arg("target"), py::arg("maxIteration"),
             R"delim(
          smoothen the mesh using bending force
      )delim");
#pragma endregion system

#pragma region parameters
  // ==========================================================
  // =============   Simulation parameters      ===============
  // ==========================================================
  py::class_<Parameters::Boundary> boundary(pymem3dg, "Boundary", R"delim(
        The boundary conditions
    )delim");
  boundary.def(py::init<>());
  boundary.def_readwrite("shapeBoundaryCondition",
                         &Parameters::Boundary::shapeBoundaryCondition,
                         R"delim(
          get the option of "roller", "pin", "fixed", or "none" to specify shape boundary condition.
      )delim");
  boundary.def_readwrite("proteinBoundaryCondition",
                         &Parameters::Boundary::proteinBoundaryCondition,
                         R"delim(
          get the option of "pin", or "none" to specify protein boundary condition.
      )delim");

  py::class_<Parameters::Variation> variation(pymem3dg, "Variation", R"delim(
        Variation
    )delim");
  variation.def(py::init<>());
  variation.def_readwrite("isProteinVariation",
                          &Parameters::Variation::isProteinVariation,
                          R"delim(
          get the option of whether simulate protein variation
      )delim");
  variation.def_readwrite("isProteinConservation",
                          &Parameters::Variation::isProteinConservation,
                          R"delim(
          get the option of whether conserve protein mass
      )delim");
  variation.def_readwrite("isShapeVariation",
                          &Parameters::Variation::isShapeVariation,
                          R"delim(
          get the option of whether simulate shape variation
      )delim");
  variation.def_readwrite("geodesicMask", &Parameters::Variation::geodesicMask,
                          R"delim(
          get domain of shape variation
      )delim");

  py::class_<Parameters::Bending> bending(pymem3dg, "Bending", R"delim(
        The bending parameters
    )delim");
  bending.def(py::init<>());
  bending.def_readwrite("Kd", &Parameters::Bending::Kd,
                        R"delim(
          get deviatoric rigidity of the membrane 
      )delim");
  bending.def_readwrite("Kdc", &Parameters::Bending::Kdc,
                        R"delim(
          get constant of deviatoric modulus vs protein density
      )delim");
  bending.def_readwrite("Kb", &Parameters::Bending::Kb,
                        R"delim(
          get Bending rigidity of the bare membrane 
      )delim");
  bending.def_readwrite("Kbc", &Parameters::Bending::Kbc,
                        R"delim(
          get constant of bending modulus vs protein density
      )delim");
  bending.def_readwrite("H0c", &Parameters::Bending::H0c,
                        R"delim(
          get constant of spontaneous curvature vs protein density
      )delim");
  bending.def_readwrite("relation", &Parameters::Bending::relation,
                        R"delim(
          get relation between H0 and protein densit, "linear" or "hill"
      )delim");

  py::class_<Parameters::Tension> tension(pymem3dg, "Tension", R"delim(
        The surface tension parameters
    )delim");
  tension.def(py::init<>());
  tension.def_readwrite("isConstantSurfaceTension",
                        &Parameters::Tension::isConstantSurfaceTension,
                        R"delim(
          get the option of whether adopt constant surface tension 
      )delim");
  tension.def_readwrite("Ksg", &Parameters::Tension::Ksg,
                        R"delim(
          get Global stretching modulus 
      )delim");
  tension.def_readwrite("A_res", &Parameters::Tension::A_res,
                        R"delim(
          get area reservoir
      )delim");
  tension.def_readwrite("At", &Parameters::Tension::At,
                        R"delim(
          get preferred surface area
      )delim");
  tension.def_readwrite("lambdaSG", &Parameters::Tension::lambdaSG,
                        R"delim(
          get augmented Lagrangian parameter for area
      )delim");

  py::class_<Parameters::Osmotic> osmotic(pymem3dg, "Osmotic", R"delim(
        The osmotic pressure parameters
    )delim");
  osmotic.def_readwrite("Kv", &Parameters::Osmotic::Kv,
                        R"delim(
          get Volume regularization 
      )delim");
  osmotic.def_readwrite("V_res", &Parameters::Osmotic::V_res,
                        R"delim(
          get volume reservoir
      )delim");
  osmotic.def_readwrite("Vt", &Parameters::Osmotic::Vt,
                        R"delim(
          get the preferred volume
      )delim");
  osmotic.def_readwrite("cam", &Parameters::Osmotic::cam,
                        R"delim(
          get the ambient concentration 
      )delim");
  osmotic.def_readwrite("n", &Parameters::Osmotic::n,
                        R"delim(
          get the enclosed solute amount
      )delim");
  osmotic.def_readwrite("isPreferredVolume",
                        &Parameters::Osmotic::isPreferredVolume,
                        R"delim(
          get the option of whether adopt the preferred volume parametrization  
      )delim");
  osmotic.def_readwrite("isConstantOsmoticPressure",
                        &Parameters::Osmotic::isConstantOsmoticPressure,
                        R"delim(
          get the option of whether adopt constant osmotic pressure
      )delim");
  osmotic.def_readwrite("lambdaV", &Parameters::Osmotic::lambdaV,
                        R"delim(
          get augmented Lagrangian parameter for volume
      )delim");

  py::class_<Parameters::Adsorption> adsorption(pymem3dg, "Adsorption",
                                                R"delim(
        The adsorption parameters
    )delim");
  adsorption.def_readwrite("epsilon", &Parameters::Adsorption::epsilon,
                           R"delim(
          get adsorption energy per protein
      )delim");

  py::class_<Parameters::Aggregation> aggregation(pymem3dg, "Aggregation",
                                                  R"delim(
        The aggregation parameters
    )delim");
  aggregation.def_readwrite("chi", &Parameters::Aggregation::chi,
                            R"delim(
          get aggregation energy 
      )delim");

  py::class_<Parameters::Entropy> entropy(pymem3dg, "Entropy",
                                          R"delim(
        The entropy parameters
    )delim");
  entropy.def_readwrite("xi", &Parameters::Entropy::xi,
                        R"delim(
          get entropy parameters
      )delim");

  py::class_<Parameters::External> external(pymem3dg, "External",
                                            R"delim(
        The external force parameters
    )delim");
  external.def(
      "setForm",
      [](Parameters::External &external,
         std::function<EigenVectorX3dr(EigenVectorX3dr, EigenVectorX1d, double,
                                       EigenVectorX1d)>
             &externalForceFunction) { external.form = externalForceFunction; },
      R"delim(
          get the vertex position matrix
      )delim");

  py::class_<Parameters::DPD> dpd(pymem3dg, "DPD",
                                  R"delim(
        The DPD parameters
    )delim");
  dpd.def_readwrite("gamma", &Parameters::DPD::gamma,
                    R"delim(
          get Dissipation coefficient 
      )delim");

  py::class_<Parameters::Dirichlet> dirichlet(pymem3dg, "Dirichlet",
                                              R"delim(
        The Dirichlet energy parameters
    )delim");
  dirichlet.def_readwrite("eta", &Parameters::Dirichlet::eta,
                          R"delim(
          get coefficient
      )delim");

  py::class_<Parameters::SelfAvoidance> selfAvoidance(pymem3dg, "SelfAvoidance",
                                                      R"delim(
        The SelfAvoidance energy parameters
    )delim");
  selfAvoidance.def_readwrite("d", &Parameters::SelfAvoidance::d,
                              R"delim(
          get coefficient of limit distance
      )delim");
  selfAvoidance.def_readwrite("mu", &Parameters::SelfAvoidance::mu,
                              R"delim(
          get coefficient of penalty coefficient
      )delim");
  selfAvoidance.def_readwrite("n", &Parameters::SelfAvoidance::n,
                              R"delim(
          get the number excluding neighborhood layers 
      )delim");
  selfAvoidance.def_readwrite("p", &Parameters::SelfAvoidance::p,
                              R"delim(
          get the period factor of self-avoidance computation
      )delim");

  py::class_<Parameters::Point> point(pymem3dg, "Point",
                                      R"delim(
        The Point energy parameters
    )delim");
  point.def_readwrite("pt", &Parameters::Point::pt,
                      R"delim(
          get the point
      )delim");
  point.def_readwrite("isFloatVertex", &Parameters::Point::isFloatVertex,
                      R"delim(
          whether use floating vertex option
      )delim");

  py::class_<Parameters::Protein> protein(pymem3dg, "Protein",
                                          R"delim(
        The protein distribution parameters
    )delim");
  protein.def_readwrite("profile", &Parameters::Protein::profile, R"delim(
          get the profile type: 'gaussian' or 'tanh'
      )delim");
  protein.def_readwrite(
      "geodesicProteinDensityDistribution",
      &Parameters::Protein::geodesicProteinDensityDistribution, R"delim(
          get (initial) protein density
      )delim");
  protein.def_readwrite("tanhSharpness", &Parameters::Protein::tanhSharpness,
                        R"delim(
          get protein density sharpness of tanh transition
      )delim");
  protein.def_readwrite("proteinInteriorPenalty",
                        &Parameters::Protein::proteinInteriorPenalty,
                        R"delim(
          get interior point parameter for protein density
      )delim");

  py::class_<Parameters::Spring> spring(pymem3dg, "spring",
                                        R"delim(
        The mesh mutator settings 
    )delim");
  spring.def_readwrite("Kst", &Parameters::Spring::Kst,
                       R"delim(
          get Vertex shifting constant 
      )delim");
  spring.def_readwrite("Ksl", &Parameters::Spring::Ksl,
                       R"delim(
          get Local stretching modulus 
      )delim");
  spring.def_readwrite("Kse", &Parameters::Spring::Kse,
                       R"delim(
          get Edge spring constant 
      )delim");

  py::class_<Parameters> parameters(pymem3dg, "Parameters", R"delim(
        The parameters
    )delim");
  parameters.def(py::init<>());
  //   parameters.def(
  //       py::init<double, double, double, EigenVectorX1d, double, double,
  //       double,
  //                double, double, double, double, double, double, double,
  //                double, double, double, double, EigenVectorX1d, double,
  //                double, double, double, double, double, double, double,
  //                std::string>(),
  //       py::arg("Kb") = 0, py::arg("Kbc") = 0, py::arg("H0c") = 0,
  //       py::arg("protein0") = Eigen::MatrixXd::Constant(1, 1, 1),
  //       py::arg("Ksg") = 0, py::arg("A_res") = 0, py::arg("Kst") = 0,
  //       py::arg("Ksl") = 0, py::arg("Kse") = 0, py::arg("Kv") = 0,
  //       py::arg("V_re") = 0, py::arg("eta") = 0, py::arg("epsilon") = 0,
  //       py::arg("Bc") = 0, py::arg("gamma") = 0, py::arg("Vt") = -1,
  //       py::arg("cam") = 0, py::arg("temp") = 0,
  //       py::arg("pt") = Eigen::MatrixXd::Constant(1, 1, 1), py::arg("Kf") =
  //       0, py::arg("conc") = -1, py::arg("height") = 0, py::arg("radius") =
  //       0, py::arg("lambdaSG") = 0, py::arg("lambdaV") = 0,
  //       py::arg("proteinInteriorPenalty") = 1e-7, py::arg("sharpness") =
  //       20, py::arg("relation") = "linear");
  parameters.def_readwrite("bending", &Parameters::bending,
                           R"delim(
          bending parameters
      )delim");
  parameters.def_readwrite("tension", &Parameters::tension,
                           R"delim(
          tension parameters
      )delim");
  parameters.def_readwrite("osmotic", &Parameters::osmotic,
                           R"delim(
        "osmotic parameters
      )delim");
  parameters.def_readwrite("adsorption", &Parameters::adsorption,
                           R"delim(
          adsorption parameters
      )delim");
  parameters.def_readwrite("aggregation", &Parameters::aggregation,
                           R"delim(
          aggregation parameters
      )delim");
  parameters.def_readwrite("entropy", &Parameters::entropy,
                           R"delim(
          entropy parameters
      )delim");
  parameters.def_readwrite("dirichlet", &Parameters::dirichlet,
                           R"delim(
          dirichlet parameters
      )delim");
  parameters.def_readwrite("selfAvoidance", &Parameters::selfAvoidance,
                           R"delim(
          selfAvoidance parameters
      )delim");
  parameters.def_readwrite("dpd", &Parameters::dpd,
                           R"delim(
          dpd parameters
      )delim");
  parameters.def_readwrite("external", &Parameters::external,
                           R"delim(
          external parameters
      )delim");
  parameters.def_readwrite("boundary", &Parameters::boundary,
                           R"delim(
          boundary parameters
      )delim");
  parameters.def_readwrite("point", &Parameters::point,
                           R"delim(
          point parameters
      )delim");
  parameters.def_readwrite("protein", &Parameters::protein,
                           R"delim(
          protein parameters
      )delim");
  parameters.def_readwrite("variation", &Parameters::variation,
                           R"delim(
          variation parameters
      )delim");
  parameters.def_readwrite("temperature", &Parameters::temperature,
                           R"delim(
          get Temperature 
      )delim");
  parameters.def_readwrite("proteinMobility", &Parameters::proteinMobility,
                           R"delim(
          get protein mobility constant 
      )delim");
  parameters.def_readwrite("damping", &Parameters::damping,
                           R"delim(
          get damping constant 
      )delim");
  parameters.def_readwrite("spring", &Parameters::spring,
                           R"delim(
          get spring parameters
      )delim");
#pragma endregion parameters

#pragma region energy
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
  energy.def_readonly("bendingEnergy", &Energy::bendingEnergy,
                      R"delim(
          get bending energy of the membrane  
      )delim");
  energy.def_readonly("deviatoricEnergy", &Energy::deviatoricEnergy,
                      R"delim(
          get deviatoric energy of the membrane  
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

#pragma endregion energy

#pragma region mesh_io
  // ==========================================================
  // =============      mesh generation    ===============
  // ==========================================================
  pymem3dg.def(
      "getCylinder", &getCylinderMatrix,
      "get topology and vertex position matrix of a non-capped cylinder",
      py::arg("radius"), py::arg("radialSubdivision"),
      py::arg("axialSubdivision"), py::arg("frequency") = 1,
      py::arg("amplitude") = 0);

  pymem3dg.def("getIcosphere", &getIcosphereMatrix,
               "get topology and vertex position matrix of icosphere",
               py::arg("radius"), py::arg("subdivision") = 0);

  pymem3dg.def("getTetrahedron", &getTetrahedronMatrix,
               "get topology and vertex position matrix of tetrahedron");

  pymem3dg.def("getDiamond", &getDiamondMatrix,
               "get topology and vertex position matrix of diamond",
               py::arg("dihedral"));

  pymem3dg.def("getHexagon", &getHexagonMatrix,
               "get topology and vertex position matrix of Hexagon",
               py::arg("radius"), py::arg("subdivision") = 0);

  pymem3dg.def(
      "subdivide",
      py::overload_cast<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &,
                        Eigen::Matrix<double, Eigen::Dynamic, 3> &,
                        std::size_t>(&subdivide),
      "subdivide the mesh", py::arg("face"), py::arg("vertex"),
      py::arg("nSub"));

  pymem3dg.def(
      "loopSubdivide",
      py::overload_cast<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &,
                        Eigen::Matrix<double, Eigen::Dynamic, 3> &,
                        std::size_t>(&loopSubdivide),
      "subdivide the mesh in Loop scheme", py::arg("face"), py::arg("vertex"),
      py::arg("nSub"));

  pymem3dg.def("readMesh", &readMesh,
               "read vertex and face matrix from .ply file",
               py::arg("plyName"));

  pymem3dg.def("readData", py::overload_cast<std::string &>(&readData),
               "read data in the format of matrix from .ply file",
               py::arg("plyName"));
  pymem3dg.def("readData",
               py::overload_cast<std::string &, std::string &>(&readData),
               "read data in the format of matrix from .ply file",
               py::arg("plyName"), py::arg("elementName"));
  pymem3dg.def(
      "readData",
      py::overload_cast<std::string &, std::string &, std::string &>(&readData),
      "read data in the format of matrix from .ply file", py::arg("plyName"),
      py::arg("elementName"), py::arg("propertyName"));

  pymem3dg.def("processSoup", &processSoup, "process polygon soup",
               py::arg("meshName"));

#pragma endregion mesh_io
};
} // namespace integrator
} // namespace solver
} // namespace mem3dg
