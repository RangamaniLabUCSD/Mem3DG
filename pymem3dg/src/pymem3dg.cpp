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
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Eigen/src/Core/util/Constants.h"
#include "mem3dg/solver/mem3dg.h"
#include "mem3dg/solver/mesh.h"

#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>

#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/system.h"
#include "pybind11/cast.h"

namespace gc = ::geometrycentral;
namespace mem3dg {
namespace py = pybind11;

// Initialize the `pymem3dg` module
PYBIND11_MODULE(pymem3dg, pymem3dg) {
  pymem3dg.doc() = "Python wrapper around the DDG solver C++ library.";

  /// Integrator object
  py::class_<Integrator> integrator(pymem3dg, "Integrator", R"delim(
        The integrator
    )delim");
  integrator.def(py::init<System &, double, bool, double, double, double,
                          std::string, std::string, size_t>(),
                 py::arg("f"), py::arg("dt_"), py::arg("isAdaptiveStep_"),
                 py::arg("total_time_"), py::arg("tSave_"),
                 py::arg("tolerance_"), py::arg("outputDir_"),
                 py::arg("trajFileName_"), py::arg("verbosity_"),
                 R"delim(
        Integrator constructor
      )delim");
  integrator.def("saveData", &Integrator::saveData,
                 R"delim(
          save data to output directory
      )delim");
  integrator.def("saveRichData", (&Integrator::saveRichData),
                 "save to richData and output .ply file to output directory",
                 R"delim(
          save data to output directory
      )delim");
#ifdef MEM3DG_WITH_NETCDF
  integrator.def("saveNetcdfData", &Integrator::saveNetcdfData,
                 "save to netcdf file in output directory",
                 R"delim(
          save data to output directory
      )delim");
#endif
  /// Integrator-velocity verlet object
  py::class_<VelocityVerlet> velocityverlet(pymem3dg, "VelocityVerlet",
                                            R"delim(
        Velocity Verlet integration
    )delim");

  velocityverlet.def(py::init<System &, double, bool, double, double, double,
                              std::string, std::string, size_t>(),
                     py::arg("f"), py::arg("dt_"), py::arg("isAdaptiveStep_"),
                     py::arg("total_time_"), py::arg("tSave_"),
                     py::arg("tolerance_"), py::arg("outputDir_"),
                     py::arg("trajFileName_"), py::arg("verbosity_"),
                     R"delim(
        Velocity Verlet integrator constructor
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
                     R"delim(
          save data to output directory
      )delim");
  velocityverlet.def("step", &VelocityVerlet::step, py::arg("n"),
                     R"delim(
          step for n iterations
      )delim");

  /// Integrator-euler object
  py::class_<Euler> euler(pymem3dg, "Euler", R"delim(
        forward euler (gradient descent) integration
    )delim");

  euler.def(py::init<System &, double, bool, double, double, double,
                     std::string, std::string, size_t, bool, double, double>(),
            py::arg("f"), py::arg("dt_"), py::arg("isAdaptiveStep_"),
            py::arg("total_time_"), py::arg("tSave_"), py::arg("tolerance_"),
            py::arg("outputDir_"), py::arg("trajFileName_"),
            py::arg("verbosity_"), py::arg("isBacktrack_"), py::arg("rho_"),
            py::arg("c1_"),
            R"delim(
        Euler integrator (steepest descent) constructor
      )delim");
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
  euler.def("saveData", &Euler::saveData,
            R"delim(
          save data to output directory
      )delim");
  euler.def("step", &Euler::step, py::arg("n"),
            R"delim(
          step for n iterations
      )delim");

  /// Integrator-conjugate gradient object
  py::class_<ConjugateGradient> conjugategradient(pymem3dg, "ConjugateGradient",
                                                  R"delim(
        conjugate Gradient propagator
    )delim");

  conjugategradient.def(
      py::init<System &, double, bool, double, double, double, std::string,
               std::string, size_t, bool, double, double, double, bool,
               size_t>(),
      py::arg("f"), py::arg("dt_"), py::arg("isAdaptiveStep_"),
      py::arg("total_time_"), py::arg("tSave_"), py::arg("tolerance_"),
      py::arg("outputDir_"), py::arg("trajFileName_"), py::arg("verbosity_"),
      py::arg("isBacktrack_"), py::arg("rho_"), py::arg("c1_"), py::arg("ctol"),
      py::arg("isAugmentedLagrangian_"), py::arg("restartNum_"),
      R"delim(
        Conjugate Gradient optimizer constructor
      )delim");
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
  conjugategradient.def("saveData", &ConjugateGradient::saveData,
                        R"delim(
          save data to output directory
      )delim");
  conjugategradient.def("step", &ConjugateGradient::step, py::arg("n"),
                        R"delim(
          step for n iterations
      )delim");

  /// Integrator-BFGS object
  py::class_<BFGS> bfgs(pymem3dg, "BFGS",
                        R"delim(
        conjugate Gradient propagator
    )delim");

  bfgs.def(py::init<System &, double, bool, double, double, double, std::string,
                    std::string, size_t, bool, double, double, double, bool>(),
           py::arg("f"), py::arg("dt_"), py::arg("isAdaptiveStep_"),
           py::arg("total_time_"), py::arg("tSave_"), py::arg("tolerance_"),
           py::arg("outputDir_"), py::arg("trajFileName_"),
           py::arg("verbosity_"), py::arg("isBacktrack_"), py::arg("rho_"),
           py::arg("c1_"), py::arg("ctol"), py::arg("isAugmentedLagrangian_"),
           R"delim(
        BFGS optimizer constructor
      )delim");
  bfgs.def("integrate", &BFGS::integrate,
           R"delim(
          integrate 
      )delim");
  bfgs.def("status", &BFGS::status,
           R"delim(
          status computation and thresholding
      )delim");
  bfgs.def("march", &BFGS::march,
           R"delim(
          stepping forward 
      )delim");
  bfgs.def("saveData", &BFGS::saveData,
           R"delim(
          save data to output directory
      )delim");
  bfgs.def("step", &BFGS::step, py::arg("n"),
           R"delim(
          step for n iterations
      )delim");

  /// Mesh mutator object
  py::class_<MeshMutator> meshmutator(pymem3dg, "MeshMutator",
                                      R"delim(
        The mesh mutator settings 
    )delim");

  meshmutator.def(py::init<>(),
                  R"delim(
       meshmutator constructor
      )delim");

  meshmutator.def_readwrite("flipNonDelaunay", &MeshMutator::flipNonDelaunay,
                            R"delim(
          whether flip non-Delaunay edge
      )delim");
  meshmutator.def_readwrite("flipNonDelaunayRequireFlat",
                            &MeshMutator::flipNonDelaunayRequireFlat,
                            R"delim(
          whether require flatness condition when flipping non-Delaunay edge
      )delim");
  meshmutator.def_readwrite("splitLarge", &MeshMutator::splitLarge,
                            R"delim(
          split edge with large faces
      )delim");
  meshmutator.def_readwrite("splitLong", &MeshMutator::splitLong,
                            R"delim(
          split long edge
      )delim");
  meshmutator.def_readwrite("splitCurved", &MeshMutator::splitCurved,
                            R"delim(
          split edge on high curvature domain
      )delim");
  meshmutator.def_readwrite("splitSharp", &MeshMutator::splitSharp,
                            R"delim(
          split edge with sharp membrane property change
      )delim");
  meshmutator.def_readwrite("splitFat", &MeshMutator::splitFat,
                            R"delim(
          split obtuse triangle
      )delim");
  meshmutator.def_readwrite("splitSkinnyDelaunay",
                            &MeshMutator::splitSkinnyDelaunay,
                            R"delim(
          split poor aspected triangle that is still Delaunay
      )delim");
  meshmutator.def_readwrite("collapseSkinny", &MeshMutator::collapseSkinny,
                            R"delim(
          collapse skinny triangles
      )delim");
  meshmutator.def_readwrite("collapseSmall", &MeshMutator::collapseSmall,
                            R"delim(
          collapse small triangles
      )delim");
  meshmutator.def_readwrite("collapseSmallNeedFlat",
                            &MeshMutator::collapseSmallNeedFlat,
                            R"delim(
         whether require flatness condition when collapsing small edge
      )delim");
  meshmutator.def_readwrite("curvTol", &MeshMutator::curvTol,
                            R"delim(
          tolerance for curvature approximation
      )delim");
  meshmutator.def_readwrite("targetFaceArea", &MeshMutator::targetFaceArea,
                            R"delim(
          target face area 
      )delim");

  /// System object
  py::class_<System> system(pymem3dg, "System",
                            R"delim(
        The system
    )delim");
  system.def(py::init<std::string, std::string, size_t, bool, Parameters &,
                      Options &>(),
             py::arg("inputMesh"), py::arg("refMesh"), py::arg("nSub"),
             py::arg("isContinue"), py::arg("p"), py::arg("o"),
             R"delim(
        System constructor
      )delim");
  system.def(
      py::init<std::string, int, size_t, bool, Parameters &, Options &>(),
      py::arg("trajFile"), py::arg("startingFrame"), py::arg("nSub"),
      py::arg("isContinue"), py::arg("p"), py::arg("o"),
      R"delim(
        System constructor
      )delim");
  system.def(py::init<Eigen::Matrix<double, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>, size_t,
                      Parameters &, Options &>(),
             py::arg("topologyMatrix"), py::arg("vertexMatrix"),
             py::arg("refVertexMatrix"), py::arg("nSub"), py::arg("p"),
             py::arg("o"),
             R"delim(
        System constructor
      )delim");
  system.def_readwrite("E", &System::E,
                       R"delim(
          get the Energy components struct
      )delim");
  system.def_readwrite("P", &System::P,
                       R"delim(
          get the Parameters struct
      )delim");
  system.def_readwrite("meshMutator", &System::meshMutator,
                       R"delim(
          get the mesh mutator object
      )delim");
  system.def_readonly("O", &System::O,
                      R"delim(
          get the Options struct
      )delim");
  system.def_readwrite("time", &System::time,
                       R"delim(
          get the time
      )delim");
  system.def_readwrite("surfaceArea", &System::surfaceArea,
                       R"delim(
          get the surface area of the mesh
      )delim");
  system.def_readwrite("volume", &System::volume,
                       R"delim(
          get the enclosed volume of the mesh
      )delim");
  system.def_readwrite("refVolume", &System::refVolume,
                       R"delim(
          get the reference enclosed volume of the mesh
      )delim");

  system.def(
      "getVertexPositionMatrix",
      [](System &s) {
        return gc::EigenMap<double, 3>(s.vpg->inputVertexPositions);
      },
      py::return_value_policy::reference_internal,
      R"delim(
          get the vertex position matrix
      )delim");
  system.def(
      "getVertexVelocityMatrix",
      [](System &s) { return gc::EigenMap<double, 3>(s.vel); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the vertex velocity matrix
      )delim");
  system.def(
      "getReferenceVertexPositionMatrix",
      [](System &s) {
        return gc::EigenMap<double, 3>(s.refVpg->inputVertexPositions);
      },
      py::return_value_policy::reference_internal,
      R"delim(
          get the face vertex matrix
      )delim");
  system.def(
      "getFaceVertexMatrix",
      [](System &s) { return s.mesh->getFaceVertexMatrix<size_t>(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the face vertex matrix
      )delim");
  system.def(
      "getVertexAdjacencyMatrix", [](System &s) { return s.vpg->d0; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the signed E-V vertex adjacency matrix, equivalent of d0 operator
      )delim");
  system.def(
      "getEdgeAdjacencyMatrix", [](System &s) { return s.vpg->d1; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the signed F-E edge adjacency matrix, equivalent of d1 operator
      )delim");
  system.def(
      "getBendingForce", [](System &s) { return s.F.bendingForce.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the bending Force
      )delim");
  system.def(
      "getCapillaryForce", [](System &s) { return s.F.capillaryForce.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the tension-induced capillary Force
      )delim");
  system.def(
      "getLineCapillaryForce",
      [](System &s) { return s.F.lineCapillaryForce.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the interfacial line tension
      )delim");
  system.def(
      "getLineTension", [](System &s) { return s.F.lineTension.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the interfacial line tension
      )delim");
  system.def(
      "getExternalForce", [](System &s) { return s.F.externalForce.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the externally-applied Force
      )delim");
  system.def(
      "getOsmoticForce", [](System &s) { return s.F.osmoticForce.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the osmotic force
      )delim");
  system.def(
      "getSurfaceTension", [](System &s) { return s.F.surfaceTension; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the Surface tension
      )delim");
  system.def(
      "getOsmoticPressure", [](System &s) { return s.F.osmoticPressure; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the osmotic pressure
      )delim");
  system.def(
      "getProteinDensity", [](System &s) { return s.proteinDensity.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the protein Density
      )delim");
  system.def(
      "getMeanCurvature",
      [](System &s) {
        return s.vpg->vertexMeanCurvatures.raw().array() /
               s.vpg->vertexDualAreas.raw().array();
      },
      py::return_value_policy::reference_internal,
      R"delim(
          get the mean curvature
      )delim");
  system.def(
      "getGaussianCurvature",
      [](System &s) {
        return s.vpg->vertexGaussianCurvatures.raw().array() /
               s.vpg->vertexDualAreas.raw().array();
      },
      py::return_value_policy::reference_internal,
      R"delim(
          get the Gaussian Curvature
      )delim");
  system.def(
      "getSpontaneousCurvature", [](System &s) { return s.H0.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the spontaneous curvature
      )delim");
  system.def(
      "getBendingRigidity", [](System &s) { return s.Kb.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the bending rigidity
      )delim");
  system.def(
      "getSurfaceArea", [](System &s) { return s.surfaceArea; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the surface area
      )delim");
  system.def(
      "getVolume", [](System &s) { return s.volume; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the enclosed volume
      )delim");
  system.def(
      "getLumpedMassMatrix",
      [](System &s) { return s.vpg->vertexLumpedMassMatrix; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the lumped mass matrix of the mesh
      )delim");
  system.def(
      "getCotanLaplacian", [](System &s) { return s.vpg->cotanLaplacian; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the Cotan Laplacian matrix of the mesh
      )delim");

  system.def("computeBendingForce", &System::computeBendingForce,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the bending force
      )delim");
  system.def("computeChemicalPotential", &System::computeChemicalPotential,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the chemical potential
      )delim");
  system.def("computeCapillaryForce", &System::computeCapillaryForce,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the capillary force
      )delim");
  system.def("computeOsmoticForce", &System::computeOsmoticForce,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the osmotic force
      )delim");
  system.def("computeLineCapillaryForce", &System::computeLineCapillaryForce,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the LineTensionForce
      )delim");
  system.def("computeExternalForce", &System::computeExternalForce,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the External Force
      )delim");
  system.def("computePhysicalForces", &System::computePhysicalForces,
             R"delim(
          compute all the forces
      )delim");
  system.def("computeDPDForces", &System::computeDPDForces,
             py::return_value_policy::reference_internal, py::arg("dt"),
             R"delim(
          compute the DPDForces
      )delim");
  system.def("computeFreeEnergy", &System::computeFreeEnergy,
             R"delim(
          compute the free energy of the system
      )delim");
  //   system.def("computeL1Norm", &System::computeL1Norm,
  //              R"delim(
  //                    compute error norm
  //             Args:
  //                   force (:py:class:`list`): mesh vertex force
  //         )delim");

  /// Options struct
  py::class_<Options> options(pymem3dg, "Options", R"delim(
        The options
    )delim");
  options.def(py::init<>());
  options.def_readwrite("isVertexShift", &Options::isVertexShift,
                        R"delim(
          get the option of whether do vertex shift  
      )delim");
  options.def_readwrite("isProteinVariation", &Options::isProteinVariation,
                        R"delim(
          get the option of whether simulate protein variation
      )delim");
  options.def_readwrite("isReducedVolume", &Options::isReducedVolume,
                        R"delim(
          get the option of whether adopt reduced volume  
      )delim");
  options.def_readwrite("isConstantOsmoticPressure",
                        &Options::isConstantOsmoticPressure,
                        R"delim(
          get the option of whether adopt constant osmotic pressure
      )delim");
  options.def_readwrite("isConstantSurfaceTension",
                        &Options::isConstantSurfaceTension,
                        R"delim(
          get the option of whether adopt constant surface tension 
      )delim");
  options.def_readwrite("isEdgeFlip", &Options::isEdgeFlip,
                        R"delim(
          get the option of whether do edge flip
      )delim");
  options.def_readwrite("isSplitEdge", &Options::isSplitEdge,
                        R"delim(
          get the option of whether split edge to grow mesh
      )delim");
  options.def_readwrite("isCollapseEdge", &Options::isCollapseEdge,
                        R"delim(
          get the option of whether Collapse edge to grow mesh
      )delim");
  options.def_readwrite("isRefMesh", &Options::isRefMesh,
                        R"delim(
          get the option of whether have a reference mesh
      )delim");
  options.def_readwrite("isFloatVertex", &Options::isFloatVertex,
                        R"delim(
          get the option of whether have "the" vertex floating in embedded space
      )delim");
  options.def_readwrite("boundaryConditionType",
                        &Options::boundaryConditionType,
                        R"delim(
          get the option of "roller", "pin", "fixed", or "none" to specify boundary condition.
      )delim");

  /// Parameter struct
  py::class_<Parameters> parameters(pymem3dg, "Parameters", R"delim(
        The parameters
    )delim");
  parameters.def(py::init<>());
  parameters.def(
      py::init<double, double, double, EigenVectorX1D, double, double, double,
               double, double, double, double, double, double, double, double,
               double, double, double, EigenVectorX1D, double, double, double,
               double, double, double>());
  parameters.def_readwrite("Kb", &Parameters::Kb,
                           R"delim(
          get Bending rigidity of the bare membrane 
      )delim");
  parameters.def_readwrite("Kbc", &Parameters::Kbc,
                           R"delim(
          get constant of bending modulus vs protein density
      )delim");
  parameters.def_readwrite("H0c", &Parameters::H0c,
                           R"delim(
          get constant of spontaneous curvature vs protein density
      )delim");
  parameters.def_readwrite("protein0", &Parameters::protein0,
                           R"delim(
          get setting of initial protein density
      )delim");
  parameters.def_readwrite("Ksg", &Parameters::Ksg,
                           R"delim(
          get Global stretching modulus 
      )delim");
  parameters.def_readwrite("A_res", &Parameters::A_res,
                           R"delim(
          get area reservoir
      )delim");
  parameters.def_readwrite("Kst", &Parameters::Kst,
                           R"delim(
          get Vertex shifting constant 
      )delim");
  parameters.def_readwrite("Ksl", &Parameters::Ksl,
                           R"delim(
          get Local stretching modulus 
      )delim");
  parameters.def_readwrite("Kse", &Parameters::Kse,
                           R"delim(
          get Edge spring constant 
      )delim");
  parameters.def_readwrite("Kv", &Parameters::Kv,
                           R"delim(
          get Volume regularization 
      )delim");
  parameters.def_readwrite("V_res", &Parameters::V_res,
                           R"delim(
          get volume reservoir
      )delim");
  parameters.def_readwrite("eta", &Parameters::eta,
                           R"delim(
          get Line tension 
      )delim");
  parameters.def_readwrite("epsilon", &Parameters::epsilon,
                           R"delim(
          get binding energy per protein
      )delim");
  parameters.def_readwrite("Bc", &Parameters::Bc,
                           R"delim(
          get binding constant 
      )delim");
  parameters.def_readwrite("gamma", &Parameters::gamma,
                           R"delim(
          get Dissipation coefficient 
      )delim");
  parameters.def_readwrite("Vt", &Parameters::Vt,
                           R"delim(
          get Reduced volume 
      )delim");
  parameters.def_readwrite("cam", &Parameters::cam,
                           R"delim(
          get Ambient Pressure 
      )delim");
  parameters.def_readwrite("temp", &Parameters::temp,
                           R"delim(
          get Temperature 
      )delim");
  parameters.def_readwrite("pt", &Parameters::pt,
                           R"delim(
          get specification for the point
      )delim");
  parameters.def_readwrite("Kf", &Parameters::Kf,
                           R"delim(
          get Magnitude of external force 
      )delim");
  parameters.def_readwrite("conc", &Parameters::conc,
                           R"delim(
          get level of concentration of the external force 
      )delim");
  parameters.def_readwrite("height", &Parameters::height,
                           R"delim(
          get target height 
      )delim");
  parameters.def_readwrite("radius", &Parameters::radius,
                           R"delim(
          get domain of integration 
      )delim");
  parameters.def_readwrite("lambdaSG", &Parameters::lambdaSG,
                           R"delim(
          get augmented Lagrangian parameter for area 
      )delim");
  parameters.def_readwrite("lambdaV", &Parameters::lambdaV,
                           R"delim(
          get augmented Lagrangian parameter for volume 
      )delim");

  /// Energy struct
  py::class_<Energy> energy(pymem3dg, "Energy", R"delim(
        The energy
    )delim");
  energy.def(py::init<double, double, double, double, double, double, double,
                      double, double>());
  energy.def_readwrite("totalE", &Energy::totalE,
                       R"delim(
          get total Energy of the system  
      )delim");
  energy.def_readwrite("kE", &Energy::kE,
                       R"delim(
          get kinetic energy of the membrane  
      )delim");
  energy.def_readwrite("potE", &Energy::potE,
                       R"delim(
          get potential energy of the membrane  
      )delim");
  energy.def_readwrite("BE", &Energy::BE,
                       R"delim(
          get bending energy of the membrane  
      )delim");
  energy.def_readwrite("sE", &Energy::sE,
                       R"delim(
          get stretching energy of the membrane  
      )delim");
  energy.def_readwrite("pE", &Energy::pE,
                       R"delim(
          get work of pressure within membrane  
      )delim");
  energy.def_readwrite("aE", &Energy::aE,
                       R"delim(
          get adsorption energy of the membrane protein  
      )delim");
  energy.def_readwrite("dE", &Energy::dE,
                       R"delim(
          get  line tension (dirichlet) energy of interface energy
      )delim");
  energy.def_readwrite("exE", &Energy::exE,
                       R"delim(
          get work of external force  
      )delim");
  energy.def_readwrite("inE", &Energy::inE,
                       R"delim(
          get protein interior penalty energy (numerical energy)
      )delim");

  pymem3dg.def("visualize", &visualize, py::arg("f"),
               R"delim(
          visualization of the system object
      )delim");

  pymem3dg.def(
      "snapshot_ply", &snapshot_ply,
      " Visualize .ply file in polysope with options of additional quantities",
      py::arg("fileName"), py::arg("options"), py::arg("transparency") = 1,
      py::arg("fov") = 50, py::arg("edgeWidth") = 1);

  pymem3dg.def(
      "animate_ply", &animate_ply,
      " Visualize .ply files in polysope with options of additional quantities",
      py::arg("framesDir"), py::arg("options"), py::arg("frameNum"),
      py::arg("transparency") = 1, py::arg("fov") = 50,
      py::arg("edgeWidth") = 1);

#ifdef MEM3DG_WITH_NETCDF
  pymem3dg.def(
      "snapshot_nc", &snapshot_nc, "Visualize netcdf file in single frame",
      py::arg("fileName"), py::arg("options"), py::arg("frame"),
      py::arg("transparency") = 1, py::arg("angle") = 0, py::arg("fov") = 50,
      py::arg("edgeWidth") = 1, py::arg("isShow") = true,
      py::arg("isSave") = false, py::arg("screenshotName") = "screenshot.png");

  pymem3dg.def("animate_nc", &animate_nc,
               "Animate netcdf file with options of additional quantities",
               py::arg("fileName"), py::arg("options"),
               py::arg("transparency") = 1, py::arg("fov") = 50,
               py::arg("edgeWidth") = 1);
#endif
  /// visualization quantities struct
  py::class_<Quantities> quantities(pymem3dg, "Quantities", R"delim(
        The quantities for visualization
    )delim");
  quantities.def(py::init<>());
  quantities.def(py::init<bool, bool, bool, bool, bool, bool, bool, bool, bool,
                          bool, bool, bool, bool>());
  quantities.def_readwrite("ref_coord", &Quantities::ref_coord,
                           R"delim(
        visualize reference coordinate 
      )delim");
  quantities.def_readwrite("velocity", &Quantities::velocity,
                           R"delim(
        visualize velocity
      )delim");
  quantities.def_readwrite("mean_curvature", &Quantities::mean_curvature,
                           R"delim(
        visualize mean_curvature
      )delim");
  quantities.def_readwrite("gauss_curvature", &Quantities::gauss_curvature,
                           R"delim(
        visualize gaussian curvature
      )delim");
  quantities.def_readwrite("spon_curvature", &Quantities::spon_curvature,
                           R"delim(
       visualize spontaneous curvature
      )delim");
  quantities.def_readwrite("ext_force", &Quantities::ext_force,
                           R"delim(
        visualize external force
      )delim");
  quantities.def_readwrite("physical_force", &Quantities::physical_force,
                           R"delim(
        visualize (total) physical force
      )delim");
  quantities.def_readwrite("capillary_force", &Quantities::capillary_force,
                           R"delim(
        visualize capillary_force
      )delim");
  quantities.def_readwrite("osmotic_force", &Quantities::osmotic_force,
                           R"delim(
        visualize osmotic_force
      )delim");
  quantities.def_readwrite("bending_force", &Quantities::bending_force,
                           R"delim(
        visualize bending force
      )delim");
  quantities.def_readwrite("line_force", &Quantities::line_force,
                           R"delim(
       visualize line tension force
      )delim");
  quantities.def_readwrite("mask", &Quantities::mask,
                           R"delim(
        visualize mask for integration
      )delim");
  quantities.def_readwrite("H_H0", &Quantities::H_H0,
                           R"delim(
        visualize H - H0
      )delim");
  quantities.def_readwrite("the_point", &Quantities::the_point,
                           R"delim(
        visualize "the" point
      )delim");
  quantities.def_readwrite("smoothing_mask", &Quantities::smoothing_mask,
                           R"delim(
        visualize the smoothing mask
      )delim");
  quantities.def_readwrite("chemical_potential",
                           &Quantities::chemical_potential,
                           R"delim(
        visualize total chemical potential
      )delim");
  quantities.def_readwrite("bending_potential", &Quantities::bending_potential,
                           R"delim(
        visualize bending component of chemical potential
      )delim");
  quantities.def_readwrite("diffusion_potential",
                           &Quantities::diffusion_potential,
                           R"delim(
        visualize diffusion component of chemical potential
      )delim");
  quantities.def_readwrite("adsorption_potential",
                           &Quantities::adsorption_potential,
                           R"delim(
        visualize adsorption component of chemical potential
      )delim");

  pymem3dg.def("driver_ply", &driver_ply,
               "Run single simulation starting with .ply files",
               py::arg("verbosity"), py::arg("inputMesh"), py::arg("refMesh"),
               py::arg("nSub"), py::arg("isReducedVolume"),
               py::arg("isProtein"), py::arg("isLocalCurvature"),
               py::arg("isVertexShift"), py::arg("isEdgeFlip"),
               py::arg("isGrowMesh"), py::arg("isRefMesh"),
               py::arg("isFloatVertex"), py::arg("isLaplacianMeanCurvature"),
               py::arg("Kb"), py::arg("Kbc"), py::arg("H0"), py::arg("r_H0"),
               py::arg("Kse"), py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"),
               py::arg("Kv"), py::arg("eta"), py::arg("epsilon"), py::arg("Bc"),
               py::arg("Vt"), py::arg("cam"), py::arg("gamma"), py::arg("temp"),
               py::arg("pt"), py::arg("Kf"), py::arg("conc"), py::arg("height"),
               py::arg("radius"), py::arg("h"), py::arg("T"), py::arg("eps"),
               py::arg("tSave"), py::arg("outputDir"), py::arg("integration"),
               py::arg("isBacktrack"), py::arg("rho"), py::arg("c1"),
               py::arg("ctol"), py::arg("isAugmentedLagrangian"),
               py::arg("restartNum"), py::arg("isAdaptiveStep"),
               R"delim(
                    Run single simulation starting with .ply files
               Args:
                   verbosity (:py:class:`int`): verbosity of output data
                   inputMesh (:py:class:`str`): input mesh path
                   refMesh (:py:class:`str`): reference mesh path
                   nSub (:py:class:`int`): number of subdivision 
                   isReducedVolume (:py:class:`bool`): whether adopt reduced volume parametrization
                   isProtein (:py:class:`bool`): whether consider protein binding
                   isLocalCurvature (:py:class:`bool`): whether has local spontaneous curvature profile
                   isVertexShfit (:py:class:`bool`): whether conduct vertex shift during integration
                   isEdgeFlip (:py:class:`bool`): whether conduct edge flip during integration
                   isGrowMesh (:py:class:`bool`): whether conduct mesh growth during integration
                   isRefMesh (:py:class:`bool`): whether whether have a reference mesh
                   isFloatVertex (:py:class:`bool`): whether have "the" vertex floating in embedded space
                   isLaplacianMeanCurvature (:py:class:`bool`): whether adopt Laplacian mean curvature definition, otherwise dihedral angle definition
                   Kb (:py:class:`double`): bending modulus of the bare membrane 
                   Kbc (:py:class:`double`): bending modulus of the coated membrane 
                   H0 (:py:class:`double`): spontaneous curvature of the membrane
                   r_H0 (:py:class:`list`): principal axis of elliptical domain of H0
                   Kse (:py:class:`double`): edge modulus for mesh regularization
                   Kst (:py:class:`double`): modulus for conformal regularization 
                   Ksl (:py:class:`double`): local area mesh regularization
                   Ksg (:py:class:`double`): global stretching modulus of membrane
                   Kv (:py:class:`double`): pressure modulus 
                   eta (:py:class:`double`): interfacial linetension
                   epsilon (:py:class:`double`): adsorption energy per unit of protein 
                   Bc (:py:class:`double`): binding constant of the protein
                   Vt (:py:class:`double`): targeted reduced volume of closed membrane 
                   cam (:py:class:`double`): anbient "concentration" outside closed membrane
                   gamma (:py:class:`double`): dissipation coefficient of DPD force
                   kt (:py:class:`double`): stochastic coeffcient of DPD force 
                   pt (:py:class:`list`): 3-D spatial coordinate of a point 
                   Kf (:py:class:`double`): "spring" constant of external force 
                   conc (:py:class:`double`): extent of local concentration of external force 
                   height (:py:class:`double`): targeted height of the external force
                   radius (:py:class:`double`): radius of integration 
                   h (:py:class:`double`): (time) step size
                   T (:py:class:`double`): maximum duration of integration
                   eps (:py:class:`double`): tolerance of convergence 
                   tSave (:py:class:`double`): time period for saving data 
                   outputDir (:py:class:`str`): output directory path
                   integration (:py:class:`str`): method of integration (optimization)
                   isBacktrack (:py:class:`bool`): whether conduct backtracking 
                   rho (:py:class:`double`): discount of step size when backtracking
                   c1 (:py:class:`double`): const of Wolfe condition 0 < c1 < 1, usually ~ 1e-4
                   ctol (:py:class:`double`): tolerance of constraints
                   isAugmentedLagrangian (:py:class:`bool`): whether use augmented lagrangian method
                   restartNum (:py::class:`int`): number of iteration to restart conjugate gradient with steepest descent
               Returns:
                   :py:class:`int`: success.
            )delim");

  pymem3dg.def("forwardsweep_ply", &forwardsweep_ply,
               "Run forward sweep simulation starting with .ply files",
               py::arg("inputMesh"), py::arg("refMesh"), py::arg("nSub"),
               py::arg("isReducedVolume"), py::arg("isProtein"),
               py::arg("isLocalCurvature"), py::arg("isVertexShift"),
               py::arg("isEdgeFlip"), py::arg("isGrowMesh"),
               py::arg("isRefMesh"), py::arg("isFloatVertex"),
               py::arg("isLaplacianMeanCurvature"), py::arg("Kb"),
               py::arg("Kbc"), py::arg("H0"), py::arg("r_H0"), py::arg("Kse"),
               py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"), py::arg("Kv"),
               py::arg("eta"), py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"),
               py::arg("cam"), py::arg("gamma"), py::arg("temp"), py::arg("pt"),
               py::arg("Kf"), py::arg("conc"), py::arg("height"),
               py::arg("radius"), py::arg("h"), py::arg("T"), py::arg("eps"),
               py::arg("tSave"), py::arg("outputDir"), py::arg("isBacktrack"),
               py::arg("rho"), py::arg("c1"), py::arg("ctol"),
               py::arg("isAugmentedLagrangian"), py::arg("restartNum"),
               py::arg("isAdaptiveStep"),
               R"delim(
                    Run forward sweep simulation starting with .ply files
               Args:
                   inputMesh (:py:class:`str`): input mesh path
                   refMesh (:py:class:`str`): reference mesh path
                   nSub (:py:class:`int`): number of subdivision 
                   isReducedVolume (:py:class:`bool`): whether adopt reduced volume parametrization
                   isProtein (:py:class:`bool`): whether consider protein binding
                   isLocalCurvature (:py:class:`bool`): whether has local spontaneous curvature profile
                   isVertexShfit (:py:class:`bool`): whether conduct vertex shift during integration
                   isEdgeFlip (:py:class:`bool`): whether conduct edge flip during integration
                   isGrowMesh (:py:class:`bool`): whether conduct mesh growth during integration
                   isRefMesh (:py:class:`bool`): whether whether have a reference mesh
                   isFloatVertex (:py:class:`bool`): whether have "the" vertex floating in embedded space
                   isLaplacianMeanCurvature (:py:class:`bool`): whether adopt Laplacian mean curvature definition, otherwise dihedral angle definition
                   Kb (:py:class:`double`): bending modulus of the bare membrane 
                   Kbc (:py:class:`double`): bending modulus of the coated membrane 
                   H0 (:py:class:`double`): spontaneous curvature of the membrane
                   r_H0 (:py:class:`list`): principal axis of elliptical domain of H0
                   Kse (:py:class:`double`): edge modulus for mesh regularization
                   Kst (:py:class:`double`): modulus for conformal regularization 
                   Ksl (:py:class:`double`): local area mesh regularization
                   Ksg (:py:class:`double`): global stretching modulus of membrane
                   Kv (:py:class:`double`): pressure modulus 
                   eta (:py:class:`double`): interfacial linetension
                   epsilon (:py:class:`double`): adsorption energy per unit of protein 
                   Bc (:py:class:`double`): binding constant of the protein
                   Vt (:py:class:`double`): targeted reduced volume of closed membrane 
                   cam (:py:class:`double`): anbient "concentration" outside closed membrane
                   gamma (:py:class:`double`): dissipation coefficient of DPD force
                   temp (:py:class:`double`): stochastic coeffcient of DPD force 
                   pt (:py:class:`list`): 3-D spatial coordinate of a point 
                   Kf (:py:class:`double`): "spring" constant of external force 
                   conc (:py:class:`double`): extent of local concentration of external force 
                   height (:py:class:`double`): targeted height of the external force
                   radius (:py:class:`double`): radius of integration 
                   h (:py:class:`double`): (time) step size
                   T (:py:class:`double`): maximum duration of integration
                   eps (:py:class:`double`): tolerance of convergence 
                   tSave (:py:class:`double`): time period for saving data 
                   outputDir (:py:class:`str`): output directory path
                   isBacktrack (:py:class:`bool`): whether conduct backtracking 
                   rho (:py:class:`double`): discount of step size when backtracking
                   c1 (:py:class:`double`): const of Wolfe condition 0 < c1 < 1, usually ~ 1e-4
                   ctol (:py:class:`double`): tolerance of constraints
                   isAugmentedLagrangian (:py:class:`bool`): whether use augmented lagrangian method
                   restartNum (:py::class:`int`): number of iteration to restart conjugate gradient with steepest descent
               Returns:
                   :py:class:`int`: success.
            )delim");

  pymem3dg.def(
      "getCylinder", &getCylinderMatrix,
      "get topology and vertex position matrix of a non-capped cylinder",
      py::arg("R"), py::arg("nR"), py::arg("nh"), py::arg("freq") = 1,
      py::arg("amp") = 0);

  pymem3dg.def("getIcosphere", &getIcosphereMatrix,
               "get topology and vertex position matrix of icosphere",
               py::arg("n"), py::arg("R"));

  pymem3dg.def("getTetrahedron", &getTetrahedronMatrix,
               "get topology and vertex position matrix of tetrahedron");

  pymem3dg.def("getDiamond", &getDiamondMatrix,
               "get topology and vertex position matrix of diamond",
               py::arg("dihedral"));

#ifdef MEM3DG_WITH_NETCDF

  pymem3dg.def(
      "driver_nc", &driver_nc,
      "Run single simulation starting with netcdf files", py::arg("verbosity"),
      py::arg("trajFile"), py::arg("startingFrame"), py::arg("nSub"),
      py::arg("isContinue"), py::arg("isReducedVolume"), py::arg("isProtein"),
      py::arg("isLocalCurvature"), py::arg("isVertexShift"),
      py::arg("isEdgeFlip"), py::arg("isGrowMesh"), py::arg("isRefMesh"),
      py::arg("isFloatVertex"), py::arg("isLaplacianMeanCurvature"),
      py::arg("Kb"), py::arg("Kbc"), py::arg("H0"), py::arg("r_H0"),
      py::arg("Kse"), py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"),
      py::arg("Kv"), py::arg("eta"), py::arg("epsilon"), py::arg("Bc"),
      py::arg("Vt"), py::arg("cam"), py::arg("gamma"), py::arg("temp"),
      py::arg("pt"), py::arg("Kf"), py::arg("conc"), py::arg("height"),
      py::arg("radius"), py::arg("h"), py::arg("T"), py::arg("eps"),
      py::arg("tSave"), py::arg("outputDir"), py::arg("integration"),
      py::arg("isBacktrack"), py::arg("rho"), py::arg("c1"), py::arg("ctol"),
      py::arg("isAugmentedLagrangian"), py::arg("restartNum"),
      py::arg("isAdaptiveStep"),
      R"delim(
                   Run single simulation starting with netcdf files
               Args:
                   verbosity (:py:class:`int`): verbosity of output data
                   trajFile (:py:class:`str`): input trajectory file path
                   startingFrame (:py:class:`int`): starting frame of continuation
                   nSub (:py:class:`int`): number of loop subdivision
                   isContinue (:py:class:`bool`): whether continue the simulation from trajectory
                   isReducedVolume (:py:class:`bool`): whether adopt reduced volume parametrization
                   isProtein (:py:class:`bool`): whether consider protein binding
                   isLocalCurvature (:py:class:`bool`): whether has local spontaneous curvature profile
                   isVertexShfit (:py:class:`bool`): whether conduct vertex shift during integration
                   isEdgeFlip (:py:class:`bool`): whether conduct edge flip during integration
                   isGrowMesh (:py:class:`bool`): whether conduct mesh growth during integration
                   isRefMesh (:py:class:`bool`): whether whether have a reference mesh
                   isFloatVertex (:py:class:`bool`): whether have "the" vertex floating in embedded space
                   isLaplacianMeanCurvature (:py:class:`bool`): whether adopt Laplacian mean curvature definition, otherwise dihedral angle definition
                   Kb (:py:class:`double`): bending modulus of the bare membrane 
                   Kbc (:py:class:`double`): bending modulus of the coated membrane 
                   H0 (:py:class:`double`): spontaneous curvature of the membrane
                   r_H0 (:py:class:`list`): principal axis of elliptical domain of H0
                   Kse (:py:class:`double`): edge modulus for mesh regularization
                   Kst (:py:class:`double`): modulus for conformal regularization 
                   Ksl (:py:class:`double`): local area mesh regularization
                   Ksg (:py:class:`double`): global stretching modulus of membrane
                   Kv (:py:class:`double`): pressure modulus 
                   eta (:py:class:`double`): interfacial linetension
                   epsilon (:py:class:`double`): adsorption energy per unit of protein 
                   Bc (:py:class:`double`): binding constant of the protein
                   Vt (:py:class:`double`): targeted reduced volume of closed membrane 
                   cam (:py:class:`double`): anbient "concentration" outside closed membrane
                   gamma (:py:class:`double`): dissipation coefficient of DPD force
                   temp (:py:class:`double`): stochastic coeffcient of DPD force 
                   pt (:py:class:`list`): 3-D spatial coordinate of a point 
                   Kf (:py:class:`double`): "spring" constant of external force 
                   conc (:py:class:`double`): extent of local concentration of external force 
                   height (:py:class:`double`): targeted height of the external force
                   radius (:py:class:`double`): radius of integration 
                   h (:py:class:`double`): (time) step size
                   T (:py:class:`double`): maximum duration of integration
                   eps (:py:class:`double`): tolerance of convergence 
                   tSave (:py:class:`double`): time period for saving data 
                   outputDir (:py:class:`str`): output directory path
                   integration (:py:class:`str`): method of integration (optimization)
                   isBacktrack (:py:class:`bool`): whether conduct backtracking 
                   rho (:py:class:`double`): discount of step size when backtracking
                   c1 (:py:class:`double`): const of Wolfe condition 0 < c1 < 1, usually ~ 1e-4
                   ctol (:py:class:`double`): tolerance of constraints
                   isAugmentedLagrangian (:py:class:`bool`): whether use augmented lagrangian method
                   restartNum (:py::class:`int`): number of iteration to restart conjugate gradient with steepest descent
               Returns:
                   :py:class:`int`: success.
            )delim");

  pymem3dg.def("forwardsweep_nc", &forwardsweep_nc,
               "Run forward sweep simulation starting with netcdf files",
               py::arg("trajFile"), py::arg("startingFrame"), py::arg("nSub"),
               py::arg("isContinue"), py::arg("isReducedVolume"),
               py::arg("isProtein"), py::arg("isLocalCurvature"),
               py::arg("isVertexShift"), py::arg("isEdgeFlip"),
               py::arg("isGrowMesh"), py::arg("isRefMesh"),
               py::arg("isFloatVertex"), py::arg("isLaplacianMeanCurvature"),
               py::arg("Kb"), py::arg("Kbc"), py::arg("H0"), py::arg("r_H0"),
               py::arg("Kse"), py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"),
               py::arg("Kv"), py::arg("eta"), py::arg("epsilon"), py::arg("Bc"),
               py::arg("Vt"), py::arg("cam"), py::arg("gamma"), py::arg("temp"),
               py::arg("pt"), py::arg("Kf"), py::arg("conc"), py::arg("height"),
               py::arg("radius"), py::arg("h"), py::arg("T"), py::arg("eps"),
               py::arg("tSave"), py::arg("outputDir"), py::arg("isBacktrack"),
               py::arg("rho"), py::arg("c1"), py::arg("ctol"),
               py::arg("isAugmentedLagrangian"), py::arg("restartNum"),
               py::arg("isAdaptiveStep"),
               R"delim(
                   Run forward sweep simulation starting with netcdf files
               Args:
                   trajFile (:py:class:`str`): input trajectory file path
                   startingFrame (:py:class:`int`): starting frame of continuation
                   nSub (:py:class:`int`): number of loop subdivision
                   isContinue (:py:class:`bool`): whether continue the simulation from trajectory
                   isReducedVolume (:py:class:`bool`): whether adopt reduced volume parametrization
                   isProtein (:py:class:`bool`): whether consider protein binding
                   isLocalCurvature (:py:class:`bool`): whether has local spontaneous curvature profile
                   isVertexShfit (:py:class:`bool`): whether conduct vertex shift during integration
                   isEdgeFlip (:py:class:`bool`): whether conduct edge flip during integration
                   isGrowMesh (:py:class:`bool`): whether conduct mesh growth during integration
                   isRefMesh (:py:class:`bool`): whether whether have a reference mesh
                   isFloatVertex (:py:class:`bool`): whether have "the" vertex floating in embedded space
                   isLaplacianMeanCurvature (:py:class:`bool`): whether adopt Laplacian mean curvature definition, otherwise dihedral angle definition
                   Kb (:py:class:`double`): bending modulus of the bare membrane 
                   Kbc (:py:class:`double`): bending modulus of the coated membrane 
                   H0 (:py:class:`double`): spontaneous curvature of the membrane
                   r_H0 (:py:class:`list`): principal axis of elliptical domain of H0
                   Kse (:py:class:`double`): edge modulus for mesh regularization
                   Kst (:py:class:`double`): modulus for conformal regularization 
                   Ksl (:py:class:`double`): local area mesh regularization
                   Ksg (:py:class:`double`): global stretching modulus of membrane
                   Kv (:py:class:`double`): pressure modulus 
                   eta (:py:class:`double`): interfacial linetension
                   epsilon (:py:class:`double`): adsorption energy per unit of protein 
                   Bc (:py:class:`double`): binding constant of the protein
                   Vt (:py:class:`double`): targeted reduced volume of closed membrane 
                   cam (:py:class:`double`): anbient "concentration" outside closed membrane
                   gamma (:py:class:`double`): dissipation coefficient of DPD force
                   temp (:py:class:`double`): stochastic coeffcient of DPD force 
                   pt (:py:class:`list`): 3-D spatial coordinate of a point 
                   Kf (:py:class:`double`): "spring" constant of external force 
                   conc (:py:class:`double`): extent of local concentration of external force 
                   height (:py:class:`double`): targeted height of the external force
                   radius (:py:class:`double`): radius of integration 
                   h (:py:class:`double`): (time) step size
                   T (:py:class:`double`): maximum duration of integration
                   eps (:py:class:`double`): tolerance of convergence 
                   tSave (:py:class:`double`): time period for saving data 
                   outputDir (:py:class:`str`): output directory path
                   isBacktrack (:py:class:`bool`): whether conduct backtracking 
                   rho (:py:class:`double`): discount of step size when backtracking
                   c1 (:py:class:`double`): const of Wolfe condition 0 < c1 < 1, usually ~ 1e-4
                   ctol (:py:class:`double`): tolerance of constraints
                   isAugmentedLagrangian (:py:class:`bool`): whether use augmented lagrangian method
                   restartNum (:py::class:`int`): number of iteration to restart conjugate gradient with steepest descent
               Returns:
                   :py:class:`int`: success.
            )delim");
#endif
};
} // namespace mem3dg
