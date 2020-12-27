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

#include "mem3dg/solver/ddgsolver.h"
#include "mem3dg/solver/mesh.h"

#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>

#include "mem3dg/solver/system.h"

namespace ddgsolver {
namespace py = pybind11;

// Initialize the `pymem3dg` module
PYBIND11_MODULE(pymem3dg, pymem3dg) {
  pymem3dg.doc() = "Python wrapper around the DDG solver C++ library.";

  pymem3dg.def(
      "driver_ply", &driver_ply,
      "the driver function for input .ply mesh file ", py::arg("verbosity"),
      py::arg("inputMesh"), py::arg("refMesh"), py::arg("nSub"),
      py::arg("isTuftedLaplacian"), py::arg("isProtein"),
      py::arg("mollifyFactor"), py::arg("isVertexShift"), py::arg("Kb"),
      py::arg("H0"), py::arg("sharpness"), py::arg("r_H0"), py::arg("Kse"),
      py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"), py::arg("Kv"),
      py::arg("eta"), py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"),
      py::arg("gamma"), py::arg("kt"), py::arg("pt"), py::arg("Kf"),
      py::arg("conc"), py::arg("height"), py::arg("radius"), py::arg("h"),
      py::arg("T"), py::arg("eps"), py::arg("tSave"), py::arg("outputDir"),
      py::arg("integration"), py::arg("isBacktrack"), py::arg("rho"),
      py::arg("c1"), py::arg("ctol"), py::arg("isAugmentedLagrangian"),
      R"delim(
                    Run the driver for .ply file input

               Args:
                   verbosity (:py:class:`int`): verbosity of output data
                   inputMesh (:py:class:`str`): input mesh path
                   refMesh (:py:class:`str`): reference mesh path
                   nSub (:py:class:`int`): number of subdivision 
                   isTuftedLaplacian (:py:class:`bool`): whether adopt tufted laplacian
                   mollifyFactor (:py:class:`double`): mollify factor for tufted laplacian
                   isProtein (:py:class:`bool`): whether consider protein binding
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

  pymem3dg.def(
      "driver_ply_sweep", &driver_ply_sweep,
      "the driver function for serially sweeping H0", py::arg("inputMesh"),
      py::arg("refMesh"), py::arg("nSub"), py::arg("isTuftedLaplacian"),
      py::arg("isProtein"), py::arg("mollifyFactor"), py::arg("isVertexShift"),
      py::arg("Kb"), py::arg("H0"), py::arg("sharpness"), py::arg("r_H0"),
      py::arg("Kse"), py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"),
      py::arg("Kv"), py::arg("eta"), py::arg("epsilon"), py::arg("Bc"),
      py::arg("Vt"), py::arg("gamma"), py::arg("kt"), py::arg("pt"),
      py::arg("Kf"), py::arg("conc"), py::arg("height"), py::arg("radius"),
      py::arg("h"), py::arg("T"), py::arg("eps"), py::arg("tSave"),
      py::arg("outputDir"), py::arg("isBacktrack"), py::arg("rho"),
      py::arg("c1"), py::arg("ctol"), py::arg("isAugmentedLagrangian"),
      R"delim(
                    Run the driver for .ply file input

               Args:
                   inputMesh (:py:class:`str`): input mesh path
                   refMesh (:py:class:`str`): reference mesh path
                   nSub (:py:class:`int`): number of subdivision 
                   isTuftedLaplacian (:py:class:`bool`): whether adopt tufted laplacian
                   mollifyFactor (:py:class:`double`): mollify factor for tufted laplacian
                   isProtein (:py:class:`bool`): whether consider protein binding
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
                   isBacktrack (:py:class:`bool`): whether conduct backtracking 
                   rho (:py:class:`double`): discount of step size when backtracking
                   c1 (:py:class:`double`): const of Wolfe condition 0 < c1 < 1, usually ~ 1e-4
                   ctol (:py:class:`double`): tolerance of constraints
                   isAugmentedLagrangian (:py:class:`bool`): whether use augmented lagrangian method
               Returns:
                   :py:class:`int`: success.
            )delim");

  pymem3dg.def("viewer", &viewer, " a visualization function",
               py::arg("fileName"), py::arg("mean_curvature"),
               py::arg("spon_curvature"), py::arg("ext_pressure"),
               py::arg("physical_pressure"), py::arg("capillary_pressure"),
               py::arg("bending_pressure"), py::arg("line_pressure"));

  pymem3dg.def("viewer", &viewPly, " a visualization function",
               py::arg("fileName"));

  pymem3dg.def("genIcosphere", &genIcosphere, "Generate a icosphere .ply file",
               py::arg("nSub"), py::arg("path"), py::arg("R"));

#ifdef MEM3DG_WITH_NETCDF
  pymem3dg.def(
      "driver_nc", &driver_nc,
      " a driver function for input .nc trajectory file", py::arg("verbosity"),
      py::arg("trajFile"), py::arg("startingFrame"),
      py::arg("isTuftedLaplacian"), py::arg("isProtein"),
      py::arg("mollifyFactor"), py::arg("isVertexShift"), py::arg("Kb"),
      py::arg("H0"), py::arg("sharpness"), py::arg("r_H0"), py::arg("Kse"),
      py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"), py::arg("Kv"),
      py::arg("eta"), py::arg("epsilon"), py::arg("Bc"), py::arg("Vt"),
      py::arg("gamma"), py::arg("kt"), py::arg("pt"), py::arg("Kf"),
      py::arg("conc"), py::arg("height"), py::arg("radius"), py::arg("h"),
      py::arg("T"), py::arg("eps"), py::arg("tSave"), py::arg("outputDir"),
      py::arg("integration"), py::arg("isBacktrack"), py::arg("rho"),
      py::arg("c1"), py::arg("ctol"), py::arg("isAugmentedLagrangian"),
      R"delim(
                   Run the driver for netcdf input file (continuation)

               Args:
                   verbosity (:py:class:`int`): verbosity of output data
                   trajFile (:py:class:`str`): input trajectory file path
                   startingFrame (:py:class:`int`): starting frame of continuation
                   isTuftedLaplacian (:py:class:`bool`): whether adopt tufted laplacian
                   mollifyFactor (:py:class:`double`): mollify factor for tufted laplacian
                   isProtein (:py:class:`bool`): whether consider protein binding
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

  pymem3dg.def(
      "driver_nc_sweep", &driver_nc_sweep,
      " a driver function for input .nc trajectory file", py::arg("trajFile"),
      py::arg("startingFrame"), py::arg("isTuftedLaplacian"),
      py::arg("isProtein"), py::arg("mollifyFactor"), py::arg("isVertexShift"),
      py::arg("Kb"), py::arg("H0"), py::arg("sharpness"), py::arg("r_H0"),
      py::arg("Kse"), py::arg("Kst"), py::arg("Ksl"), py::arg("Ksg"),
      py::arg("Kv"), py::arg("eta"), py::arg("epsilon"), py::arg("Bc"),
      py::arg("Vt"), py::arg("gamma"), py::arg("kt"), py::arg("pt"),
      py::arg("Kf"), py::arg("conc"), py::arg("height"), py::arg("radius"),
      py::arg("h"), py::arg("T"), py::arg("eps"), py::arg("tSave"),
      py::arg("outputDir"), py::arg("isBacktrack"), py::arg("rho"),
      py::arg("c1"), py::arg("ctol"), py::arg("isAugmentedLagrangian"),
      R"delim(
                   Run the driver for netcdf input file (continuation)

               Args:
                   trajFile (:py:class:`str`): input trajectory file path
                   startingFrame (:py:class:`int`): starting frame of continuation
                   isTuftedLaplacian (:py:class:`bool`): whether adopt tufted laplacian
                   mollifyFactor (:py:class:`double`): mollify factor for tufted laplacian
                   isProtein (:py:class:`bool`): whether consider protein binding
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
                   isBacktrack (:py:class:`bool`): whether conduct backtracking 
                   rho (:py:class:`double`): discount of step size when backtracking
                   c1 (:py:class:`double`): const of Wolfe condition 0 < c1 < 1, usually ~ 1e-4
                   ctol (:py:class:`double`): tolerance of constraints
                   isAugmentedLagrangian (:py:class:`bool`): whether use augmented lagrangian method
            
               Returns:
                   :py:class:`int`: success.
            )delim");

  pymem3dg.def("view_animation", &view_animation, " a visualization function",
               py::arg("fileName"), py::arg("ref_coord"), py::arg("velocity"),
               py::arg("mean_curvature"), py::arg("spon_curvature"),
               py::arg("ext_pressure"), py::arg("physical_pressure"),
               py::arg("capillary_pressure"), py::arg("bending_pressure"),
               py::arg("line_pressure"), py::arg("mask"), py::arg("H_H0"));
#endif
};
} // namespace ddgsolver
