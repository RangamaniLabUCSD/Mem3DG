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

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mem3dg/solver/mem3dg.h"
#include "mem3dg/solver/mesh.h"

#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>

#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/system.h"

namespace gc = ::geometrycentral;
namespace mem3dg {
namespace py = pybind11;

// Initialize the `pymem3dg` module
PYBIND11_MODULE(pymem3dg, pymem3dg) {
  pymem3dg.doc() = "Python wrapper around the DDG solver C++ library.";

  /// System object
  py::class_<System> system(pymem3dg, "System",
                            R"delim(
        The system
    )delim");
  system.def(py::init<std::string, std::string, size_t, Parameters, bool, bool,
                      bool, bool>());
  system.def_readwrite("proteinDensity", &System::proteinDensity,
                       R"delim(
          get the protein Density
      )delim");
  system.def_readwrite("insidePressure", &System::insidePressure,
                       R"delim(
          get the inside pressures
      )delim");
  system.def_readwrite("E", &System::E,
                       R"delim(
          get the Energy components struct
      )delim");
  system.def_readwrite("P", &System::P,
                       R"delim(
          get the Parameters struct
      )delim");
  system.def("getBendingPressure", &System::getBendingPressure,
             R"delim(
          get the bending pressures
      )delim");
  system.def("getChemicalPotential", &System::getChemicalPotential,
             R"delim(
          get the chemical potential
      )delim");
  system.def("getCapillaryPressure", &System::getCapillaryPressure,
             R"delim(
          get the capillary Pressure
      )delim");
  system.def("getInsidePressure", &System::getInsidePressure,
             R"delim(
          get the getInsidePressure
      )delim");
  system.def("getLineTensionPressure", &System::getLineTensionPressure,
             R"delim(
          get the getLineTensionPressure
      )delim");
  system.def("getDPDForces", &System::getDPDForces,
             R"delim(
          get the getDPDForces
      )delim");
  system.def("getExternalPressure", &System::getExternalPressure,
             R"delim(
          get the getExternalPressure
      )delim");
  system.def("getAllForces", &System::getAllForces,
             R"delim(
          get all forces
      )delim");
  system.def("getFreeEnergy", &System::getFreeEnergy,
             R"delim(
          get the free energy of the system
      )delim");

  /// Parameter struct
  py::class_<Parameters> parameters(pymem3dg, "Parameters", R"delim(
        The parameters
    )delim");
  parameters.def(
      py::init<double, double, double, std::vector<double>, double, double,
               double, double, double, double, double, double, double, double,
               double, double, double, std::vector<double>, double, double,
               double, double, double, double>());
  parameters.def_readwrite("Kb", &Parameters::Kb,
                           R"delim(
          get Bending modulus 
      )delim");
  parameters.def_readwrite("H0", &Parameters::H0,
                           R"delim(
          get Spontaneous curvature 
      )delim");
  parameters.def_readwrite("sharpness", &Parameters::sharpness,
                           R"delim(
          get Sharpness of the spontaneous curvature hetergeneity 
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
        The parameters
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

  /// Velocity Verlet
  pymem3dg.def("velocityVerlet", &mem3dg::integration::velocityVerlet,
               "Run velocity verlet time integration", py::arg("f"),
               py::arg("dt"), py::arg("total_time"), py::arg("tSave"),
               py::arg("tolerance"), py::arg("verbosity"),
               py::arg("isAdaptiveStep"), py::arg("outputDir"),
               R"delim(
        )delim");

  /// Euler integration
  pymem3dg.def("euler", &mem3dg::integration::euler,
               "Run forward euler time integration", py::arg("f"),
               py::arg("dt"), py::arg("total_time"), py::arg("tSave"),
               py::arg("tolerance"), py::arg("verbosity"), py::arg("outputDir"),
               py::arg("isBacktrack"), py::arg("rho"), py::arg("c1"),
               py::arg("isAdaptiveStep"),
               R"delim(
        )delim");

  /// CG integration
  pymem3dg.def("conjugateGradient", &mem3dg::integration::conjugateGradient,
               "Run conjugate gradient time integration", py::arg("f"),
               py::arg("dt"), py::arg("total_time"), py::arg("tSave"),
               py::arg("tol"), py::arg("ctol"), py::arg("verbosity"),
               py::arg("outputDir"), py::arg("isBacktrack"), py::arg("rho"),
               py::arg("c1"), py::arg("isAugmentedLagrangian"),
               py::arg("isAdaptiveStep"), py::arg("trajFileName"),
               R"delim(
        )delim");

  /// Driver function for system generation
  pymem3dg.def("system_ply", &system_ply,
               "Run single simulation starting with .ply files",
               py::arg("verbosity"), py::arg("inputMesh"), py::arg("refMesh"),
               py::arg("nSub"), py::arg("isReducedVolume"),
               py::arg("isProtein"), py::arg("isLocalCurvature"),
               py::arg("isVertexShift"), py::arg("Kb"), py::arg("H0"),
               py::arg("sharpness"), py::arg("r_H0"), py::arg("Kse"),
               py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"), py::arg("Kv"),
               py::arg("eta"), py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"),
               py::arg("cam"), py::arg("gamma"), py::arg("temp"), py::arg("pt"),
               py::arg("Kf"), py::arg("conc"), py::arg("height"),
               py::arg("radius"), py::arg("h"), py::arg("T"), py::arg("eps"),
               py::arg("tSave"), py::arg("outputDir"), py::arg("integration"),
               py::arg("isBacktrack"), py::arg("rho"), py::arg("c1"),
               py::arg("ctol"), py::arg("isAugmentedLagrangian"),
               R"delim(
        )delim");

  pymem3dg.def(
      "viewer_ply", &viewer_ply,
      " Visualize .ply file in polysope with options of additional quantities",
      py::arg("fileName"), py::arg("mean_curvature"),
      py::arg("gauss_curvature"), py::arg("spon_curvature"),
      py::arg("ext_pressure"), py::arg("physical_pressure"),
      py::arg("capillary_pressure"), py::arg("bending_pressure"),
      py::arg("line_pressure"));

  pymem3dg.def(
      "driver_ply", &driver_ply,
      "Run single simulation starting with .ply files", py::arg("verbosity"),
      py::arg("inputMesh"), py::arg("refMesh"), py::arg("nSub"),
      py::arg("isReducedVolume"), py::arg("isProtein"),
      py::arg("isLocalCurvature"), py::arg("isVertexShift"), py::arg("Kb"),
      py::arg("H0"), py::arg("sharpness"), py::arg("r_H0"), py::arg("Kse"),
      py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"), py::arg("Kv"),
      py::arg("eta"), py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"),
      py::arg("cam"), py::arg("gamma"), py::arg("temp"), py::arg("pt"),
      py::arg("Kf"), py::arg("conc"), py::arg("height"), py::arg("radius"),
      py::arg("h"), py::arg("T"), py::arg("eps"), py::arg("tSave"),
      py::arg("outputDir"), py::arg("integration"), py::arg("isBacktrack"),
      py::arg("rho"), py::arg("c1"), py::arg("ctol"),
      py::arg("isAugmentedLagrangian"), py::arg("isAdaptiveStep"),
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
                   Kb (:py:class:`double`): bending modulus of the membrane 
                   H0 (:py:class:`double`): spontaneous curvature of the membrane
                   sharpness (:py:class:`double`): sharpness of the interfacial transition of H0
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
               Returns:
                   :py:class:`int`: success.
            )delim");

  pymem3dg.def("forwardsweep_ply", &forwardsweep_ply,
               "Run forward sweep simulation starting with .ply files",
               py::arg("inputMesh"), py::arg("refMesh"), py::arg("nSub"),
               py::arg("isReducedVolume"), py::arg("isProtein"),
               py::arg("isLocalCurvature"), py::arg("isVertexShift"),
               py::arg("Kb"), py::arg("H0"), py::arg("sharpness"),
               py::arg("r_H0"), py::arg("Kse"), py::arg("Kst"), py::arg("Ksl"),
               py::arg("Ksg"), py::arg("Kv"), py::arg("eta"),
               py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"), py::arg("cam"),
               py::arg("gamma"), py::arg("temp"), py::arg("pt"), py::arg("Kf"),
               py::arg("conc"), py::arg("height"), py::arg("radius"),
               py::arg("h"), py::arg("T"), py::arg("eps"), py::arg("tSave"),
               py::arg("outputDir"), py::arg("isBacktrack"), py::arg("rho"),
               py::arg("c1"), py::arg("ctol"), py::arg("isAugmentedLagrangian"),
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
                   Kb (:py:class:`double`): bending modulus of the membrane 
                   H0 (:py:class:`double`): spontaneous curvature of the membrane
                   sharpness (:py:class:`double`): sharpness of the interfacial transition of H0
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
               Returns:
                   :py:class:`int`: success.
            )delim");

#ifdef MEM3DG_WITH_NETCDF

  pymem3dg.def("snapshot_nc", &snapshot_nc,
               "Visualize netcdf file in single frame", py::arg("fileName"),
               py::arg("frame"), py::arg("transparency"), py::arg("angle"),
               py::arg("fov"), py::arg("edgeWidth"), py::arg("isShow"),
               py::arg("isSave"), py::arg("screenshotName"),
               py::arg("ref_coord"), py::arg("velocity"),
               py::arg("mean_curvature"), py::arg("gauss_curvature"),
               py::arg("spon_curvature"), py::arg("ext_pressure"),
               py::arg("physical_pressure"), py::arg("capillary_pressure"),
               py::arg("inside_pressure"), py::arg("bending_pressure"),
               py::arg("line_pressure"), py::arg("mask"), py::arg("H_H0"));

  pymem3dg.def("animation_nc", &animation_nc,
               "Animate netcdf file with options of additional quantities",
               py::arg("fileName"), py::arg("transparency"), py::arg("angle"),
               py::arg("fov"), py::arg("edgeWidth"), py::arg("ref_coord"),
               py::arg("velocity"), py::arg("mean_curvature"),
               py::arg("gauss_curvature"), py::arg("spon_curvature"),
               py::arg("ext_pressure"), py::arg("physical_pressure"),
               py::arg("capillary_pressure"), py::arg("inside_pressure"),
               py::arg("bending_pressure"), py::arg("line_pressure"),
               py::arg("mask"), py::arg("H_H0"));

  pymem3dg.def(
      "driver_nc", &driver_nc,
      "Run single simulation starting with netcdf files", py::arg("verbosity"),
      py::arg("trajFile"), py::arg("startingFrame"), py::arg("nSub"),
      py::arg("isContinue"), py::arg("isReducedVolume"), py::arg("isProtein"),
      py::arg("isLocalCurvature"), py::arg("isVertexShift"), py::arg("Kb"),
      py::arg("H0"), py::arg("sharpness"), py::arg("r_H0"), py::arg("Kse"),
      py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"), py::arg("Kv"),
      py::arg("eta"), py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"),
      py::arg("cam"), py::arg("gamma"), py::arg("temp"), py::arg("pt"),
      py::arg("Kf"), py::arg("conc"), py::arg("height"), py::arg("radius"),
      py::arg("h"), py::arg("T"), py::arg("eps"), py::arg("tSave"),
      py::arg("outputDir"), py::arg("integration"), py::arg("isBacktrack"),
      py::arg("rho"), py::arg("c1"), py::arg("ctol"),
      py::arg("isAugmentedLagrangian"), py::arg("isAdaptiveStep"),
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
                   Kb (:py:class:`double`): bending modulus of the membrane 
                   H0 (:py:class:`double`): spontaneous curvature of the membrane
                   sharpness (:py:class:`double`): sharpness of the interfacial transition of H0
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
               Returns:
                   :py:class:`int`: success.
            )delim");

  pymem3dg.def("forwardsweep_nc", &forwardsweep_nc,
               "Run forward sweep simulation starting with netcdf files",
               py::arg("trajFile"), py::arg("startingFrame"), py::arg("nSub"),
               py::arg("isContinue"), py::arg("isReducedVolume"),
               py::arg("isProtein"), py::arg("isLocalCurvature"),
               py::arg("isVertexShift"), py::arg("Kb"), py::arg("H0"),
               py::arg("sharpness"), py::arg("r_H0"), py::arg("Kse"),
               py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"), py::arg("Kv"),
               py::arg("eta"), py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"),
               py::arg("cam"), py::arg("gamma"), py::arg("temp"), py::arg("pt"),
               py::arg("Kf"), py::arg("conc"), py::arg("height"),
               py::arg("radius"), py::arg("h"), py::arg("T"), py::arg("eps"),
               py::arg("tSave"), py::arg("outputDir"), py::arg("isBacktrack"),
               py::arg("rho"), py::arg("c1"), py::arg("ctol"),
               py::arg("isAugmentedLagrangian"), py::arg("isAdaptiveStep"),
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
                   Kb (:py:class:`double`): bending modulus of the membrane 
                   H0 (:py:class:`double`): spontaneous curvature of the membrane
                   sharpness (:py:class:`double`): sharpness of the interfacial transition of H0
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
               Returns:
                   :py:class:`int`: success.
            )delim");
#endif
};
} // namespace mem3dg
