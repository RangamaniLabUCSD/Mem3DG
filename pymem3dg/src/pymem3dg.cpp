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

#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mem3dg/solver/mem3dg.h"
#include "mem3dg/solver/mesh.h"

#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>

#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace py = pybind11;

// Initialize the `pymem3dg` module
PYBIND11_MODULE(pymem3dg, pymem3dg) {
  pymem3dg.doc() = "Python wrapper around the DDG solver C++ library.";

  pymem3dg.def("genIcosphere", &genIcosphere, "Generate a icosphere .ply file",
               py::arg("nSub"), py::arg("path"), py::arg("R"));

  pymem3dg.def(
      "viewer_ply", &viewer_ply,
      " Visualize .ply file in polysope with options of additional quantities",
      py::arg("fileName"), py::arg("mean_curvature"), py::arg("spon_curvature"),
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

  pymem3dg.def(
      "snapshot_nc", &snapshot_nc, "Visualize netcdf file in single frame",
      py::arg("fileName"), py::arg("frame"), py::arg("transparency"),
      py::arg("angle"), py::arg("fov"), py::arg("edgeWidth"), py::arg("isShow"),
      py::arg("isSave"), py::arg("screenshotName"), py::arg("ref_coord"),
      py::arg("velocity"), py::arg("mean_curvature"), py::arg("spon_curvature"),
      py::arg("ext_pressure"), py::arg("physical_pressure"),
      py::arg("capillary_pressure"), py::arg("inside_pressure"),
      py::arg("bending_pressure"), py::arg("line_pressure"), py::arg("mask"),
      py::arg("H_H0"));

  pymem3dg.def("animation_nc", &animation_nc,
               "Animate netcdf file with options of additional quantities",
               py::arg("fileName"), py::arg("transparency"), py::arg("angle"),
               py::arg("fov"), py::arg("edgeWidth"), py::arg("ref_coord"),
               py::arg("velocity"), py::arg("mean_curvature"),
               py::arg("spon_curvature"), py::arg("ext_pressure"),
               py::arg("physical_pressure"), py::arg("capillary_pressure"),
               py::arg("inside_pressure"), py::arg("bending_pressure"),
               py::arg("line_pressure"), py::arg("mask"), py::arg("H_H0"));

  pymem3dg.def(
      "driver_nc", &driver_nc,
      "Run single simulation starting with netcdf files", py::arg("verbosity"),
      py::arg("trajFile"), py::arg("startingFrame"), py::arg("isReducedVolume"),
      py::arg("isProtein"), py::arg("isLocalCurvature"),
      py::arg("isVertexShift"), py::arg("Kb"), py::arg("H0"),
      py::arg("sharpness"), py::arg("r_H0"), py::arg("Kse"), py::arg("Kst"),
      py::arg("Ksl"), py::arg("Ksg"), py::arg("Kv"), py::arg("eta"),
      py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"), py::arg("cam"),
      py::arg("gamma"), py::arg("temp"), py::arg("pt"), py::arg("Kf"),
      py::arg("conc"), py::arg("height"), py::arg("radius"), py::arg("h"),
      py::arg("T"), py::arg("eps"), py::arg("tSave"), py::arg("outputDir"),
      py::arg("integration"), py::arg("isBacktrack"), py::arg("rho"),
      py::arg("c1"), py::arg("ctol"), py::arg("isAugmentedLagrangian"),
      py::arg("isAdaptiveStep"),
      R"delim(
                   Run single simulation starting with netcdf files
               Args:
                   verbosity (:py:class:`int`): verbosity of output data
                   trajFile (:py:class:`str`): input trajectory file path
                   startingFrame (:py:class:`int`): starting frame of continuation
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

  pymem3dg.def(
      "forwardsweep_nc", &forwardsweep_nc,
      "Run forward sweep simulation starting with netcdf files",
      py::arg("trajFile"), py::arg("startingFrame"), py::arg("isReducedVolume"),
      py::arg("isProtein"), py::arg("isLocalCurvature"),
      py::arg("isVertexShift"), py::arg("Kb"), py::arg("H0"),
      py::arg("sharpness"), py::arg("r_H0"), py::arg("Kse"), py::arg("Kst"),
      py::arg("Ksl"), py::arg("Ksg"), py::arg("Kv"), py::arg("eta"),
      py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"), py::arg("cam"),
      py::arg("gamma"), py::arg("temp"), py::arg("pt"), py::arg("Kf"),
      py::arg("conc"), py::arg("height"), py::arg("radius"), py::arg("h"),
      py::arg("T"), py::arg("eps"), py::arg("tSave"), py::arg("outputDir"),
      py::arg("isBacktrack"), py::arg("rho"), py::arg("c1"), py::arg("ctol"),
      py::arg("isAugmentedLagrangian"), py::arg("isAdaptiveStep"),
      R"delim(
                   Run forward sweep simulation starting with netcdf files
               Args:
                   trajFile (:py:class:`str`): input trajectory file path
                   startingFrame (:py:class:`int`): starting frame of continuation
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
