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

  /// System object
  py::class_<System> system(pymem3dg, "System",
                            R"delim(
        The system
    )delim");
  system.def(
      py::init<std::string, std::string, size_t, Parameters &, Options &>(),
      py::arg("inputMesh"), py::arg("refMesh"), py::arg("nSub"), py::arg("p"),
      py::arg("o"),
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
      "getBendingPressure", [](System &s) { return s.bendingPressure.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the bending Pressure
      )delim");
  system.def(
      "getCapillaryPressure",
      [](System &s) { return s.capillaryPressure.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the tension-induced capillary pressure
      )delim");
  system.def(
      "getLineCapillaryForce",
      [](System &s) { return s.lineCapillaryForce.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the interfacial line tension
      )delim");
  system.def(
      "getLineTension", [](System &s) { return s.lineTension.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the interfacial line tension
      )delim");
  system.def(
      "getExternalPressure", [](System &s) { return s.externalPressure.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the externally-applied pressure
      )delim");
  system.def(
      "getInsidePressure", [](System &s) { return s.insidePressure.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the relative inside pressure
      )delim");
  system.def(
      "getSurfaceTension", [](System &s) { return s.surfaceTension; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the Surface tension
      )delim");
  system.def(
      "getProteinDensity", [](System &s) { return s.proteinDensity.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the protein Density
      )delim");
  system.def(
      "getMeanCurvature", [](System &s) { return s.H.raw(); },
      py::return_value_policy::reference_internal,
      R"delim(
          get the mean curvature
      )delim");
  system.def(
      "getGaussianCurvature", [](System &s) { return s.K.raw(); },
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
      "getLumpedMassMatrix", [](System &s) { return s.M; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the lumped mass matrix of the mesh
      )delim");
  system.def(
      "getCotanLaplacian", [](System &s) { return s.L; },
      py::return_value_policy::reference_internal,
      R"delim(
          get the Cotan Laplacian matrix of the mesh
      )delim");

  system.def("computeBendingPressure", &System::computeBendingPressure,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the bending pressures
      )delim");
  system.def("computeChemicalPotential", &System::computeChemicalPotential,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the chemical potential
      )delim");
  system.def("computeCapillaryPressure", &System::computeCapillaryPressure,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the capillary Pressure
      )delim");
  system.def("computeInsidePressure", &System::computeInsidePressure,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the InsidePressure
      )delim");
  system.def("computeLineCapillaryForce", &System::computeLineCapillaryForce,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the LineTensionPressure
      )delim");
  system.def("computeDPDForces", &System::computeDPDForces,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the DPDForces
      )delim");
  system.def("computeExternalPressure", &System::computeExternalPressure,
             py::return_value_policy::reference_internal,
             R"delim(
          compute the ExternalPressure
      )delim");
  system.def("computePhysicalForces", &System::computePhysicalForces,
             R"delim(
          compute all the forces
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
  options.def(py::init<bool, bool, bool, bool, bool, bool, bool, bool, bool>());
  options.def_readwrite("isVertexShift", &Options::isVertexShift,
                        R"delim(
          get the option of whether do vertex shift  
      )delim");
  options.def_readwrite("isProtein", &Options::isProtein,
                        R"delim(
          get the option of whether simulate protein  
      )delim");
  options.def_readwrite("isReducedVolume", &Options::isReducedVolume,
                        R"delim(
          get the option of whether adopt reduced volume  
      )delim");
  options.def_readwrite("isLocalCurvature", &Options::isLocalCurvature,
                        R"delim(
          get the option of whether consider local curvature  
      )delim");
  options.def_readwrite("isEdgeFlip", &Options::isEdgeFlip,
                        R"delim(
          get the option of whether do edge flip
      )delim");
  options.def_readwrite("isGrowMesh", &Options::isGrowMesh,
                        R"delim(
          get the option of whether grow mesh 
      )delim");
  options.def_readwrite("isRefMesh", &Options::isRefMesh,
                        R"delim(
          get the option of whether have a reference mesh
      )delim");
  options.def_readwrite("isFloatVertex", &Options::isFloatVertex,
                        R"delim(
          get the option of whether have "the" vertex floating in embedded space
      )delim");
  options.def_readwrite("isLaplacianMeanCurvature",
                        &Options::isLaplacianMeanCurvature,
                        R"delim(
          get the option of whether adopt Laplacian mean curvature definition, 
          otherwise dihedral angle definition.
      )delim");

  /// Parameter struct
  py::class_<Parameters> parameters(pymem3dg, "Parameters", R"delim(
        The parameters
    )delim");
  parameters.def(py::init<>());
  parameters.def(
      py::init<double, double, double, std::vector<double>, double, double,
               double, double, double, double, double, double, double, double,
               double, double, double, std::vector<double>, double, double,
               double, double, double, double>());
  parameters.def_readwrite("Kb", &Parameters::Kb,
                           R"delim(
          get Bending rigidity of the bare membrane 
      )delim");
  parameters.def_readwrite("Kbc", &Parameters::Kbc,
                           R"delim(
          get Bending rigidity of the coated membrane
      )delim");
  parameters.def_readwrite("H0", &Parameters::H0,
                           R"delim(
          get Spontaneous curvature 
      )delim");
  parameters.def_readwrite("r_H0", &Parameters::r_H0,
                           R"delim(
          get radius of non-zero spontaneous curvature 
      )delim");
  parameters.def_readwrite("Ksg", &Parameters::Ksg,
                           R"delim(
          get Global stretching modulus 
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
  parameters.def_readwrite("sigma", &Parameters::sigma,
                           R"delim(
          get Noise 
      )delim");
  parameters.def_readwrite("pt", &Parameters::pt,
                           R"delim(
          get index of node with applied external force 
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
  energy.def_readwrite("cE", &Energy::cE,
                       R"delim(
          get chemical energy of the membrane protein  
      )delim");
  energy.def_readwrite("lE", &Energy::lE,
                       R"delim(
          get  line tension energy of interface   energy
      )delim");
  energy.def_readwrite("exE", &Energy::exE,
                       R"delim(
          get work of external force  
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
  quantities.def_readwrite("ext_pressure", &Quantities::ext_pressure,
                           R"delim(
        visualize external pressure
      )delim");
  quantities.def_readwrite("physical_pressure", &Quantities::physical_pressure,
                           R"delim(
        visualize (total) physical pressure
      )delim");
  quantities.def_readwrite("capillary_pressure",
                           &Quantities::capillary_pressure,
                           R"delim(
        visualize capillary_pressure
      )delim");
  quantities.def_readwrite("inside_pressure", &Quantities::inside_pressure,
                           R"delim(
        visualize inside_pressure
      )delim");
  quantities.def_readwrite("bending_pressure", &Quantities::bending_pressure,
                           R"delim(
        visualize bending pressure
      )delim");
  quantities.def_readwrite("line_pressure", &Quantities::line_pressure,
                           R"delim(
       visualize line tension pressure
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

  pymem3dg.def("getIcosphere", &getIcosphereMatrix,
               "get topology and vertex position matrix of icosphere",
               py::arg("n"), py::arg("R"));

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
