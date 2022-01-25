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

#pragma once
#include "mem3dg/mem3dg"
#include <geometrycentral/surface/surface_mesh.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

/**
 * @brief Quantities option for visualization
 *
 */
struct DLL_PUBLIC Quantities {
  /// reference coordinate
  bool ref_coord = false;
  /// vertex velocity
  bool velocity = false;
  /// vertex mean curvature
  bool mean_curvature = false;
  /// vertex gaussian curvature
  bool gauss_curvature = false;
  /// vertex spontaneous curvature
  bool spon_curvature = false;
  /// vertex external force
  bool ext_force = false;
  /// vertex self avoidance force
  bool avoidance_force = false;
  /// vertex (total) physical force
  bool physical_force = false;
  /// vertex capillary force
  bool capillary_force = false;
  /// constant osmotic force
  bool osmotic_force = false;
  /// vertex bending force
  bool bending_force = false;
  /// vertex deviatoric force
  bool deviatoric_force = false;
  /// vertex adsorption force
  bool adsorption_force = false;
  /// vertex aggregation force
  bool aggregation_force = false;
  /// vertex line tension force
  bool line_force = false;
  /// vertex mask for time integration
  bool mask = false;
  /// vertex mean - spontaneous curvature difference
  bool H_H0 = false;
  /// "the" point
  bool the_point = false;
  /// smoothing mask
  bool smoothing_mask = false;
  /// chemical potential
  bool chemical_potential = false;
  /// bending chemical potential
  bool bending_potential = false;
  /// deviatoric chemical potential
  bool deviatoric_potential = false;
  /// diffusion chemical potential
  bool diffusion_potential = false;
  /// adsorption chemical potential
  bool adsorption_potential = false;
  /// aggregation chemical potential
  bool aggregation_potential = false;
};

// ==========================================================
// =============        Viewers                ==============
// ==========================================================

/**
 * @brief visualize system object
 */
DLL_PUBLIC void visualize(mem3dg::solver::System &f);

/**
 * @brief Visualize .ply file in polysope with options of additional quantities
 */
DLL_PUBLIC int snapshot_ply(std::string fileName, const Quantities &options,
                            float transparency = 1, float fov = 50,
                            float edgeWidth = 1);

/**
 * @brief Visualize .ply files in polysope with options of additional quantities
 */
DLL_PUBLIC int animate_ply(std::string framesDir, const Quantities &options,
                           std::vector<std::size_t> frameNum,
                           double mapMinLim, double mapMaxLim,
                           float transparency = 1, float fov = 50,
                           float edgeWidth = 1);

#ifdef MEM3DG_WITH_NETCDF
/**
 * @brief Visualize netcdf file in single frame
 */
DLL_PUBLIC int snapshot_nc(std::string &filename, const Quantities &options,
                           int frame, float transparency = 1, float angle = 0,
                           float fov = 50, float edgeWidth = 1,
                           bool isShow = true, bool isSave = false,
                           std::string screenshotName = "screenshot.png");

/**
 * @brief Animate netcdf file with options of additional quantities
 */
DLL_PUBLIC int animate_nc(std::string &filename, const Quantities &options,
                          float transparency = 1, float fov = 50,
                          float edgeWidth = 1);
#endif

// ==========================================================
// =============        Helper functions       ==============
// ==========================================================
/**
 * @brief Initialize GUI window with customized options
 */
void initGui();

/**
 * @brief Wait function for slow motion animation
 */
void wait(unsigned timeout);

/**
 * @brief Play the next frame of .ply frame files
 */
void play(polyscope::SurfaceMesh *&polyscopeMesh, std::string framesDir,
          int &idx, int &waitTime, const Quantities options, bool &toggle,
          std::vector<std::size_t> frameNum, double mapMinLim = 0,
          double mapMaxLim = 0);

/**
 * @brief Register Polyscope surface mesh from .ply file with options of data
 * quantities
 */
polyscope::SurfaceMesh *registerSurfaceMesh(std::string plyName,
                                            const Quantities &options,
                                            double mapMinLim = 0,
                                            double mapMaxLim = 0);
/**
 * @brief get discrete count from every nonzero entries
 */
std::vector<std::pair<std::size_t, int>>
getCountQuantities(gc::VertexData<int> &&meshData);

#ifdef MEM3DG_WITH_NETCDF
/**
 * @brief Play the next frame of the NetCDF trajectory file
 */
void play(polyscope::SurfaceMesh *&polyscopeMesh, mem3dg::solver::TrajFile &fd,
          int &idx, int &waitTime, const Quantities options, bool &toggle);
void play(polyscope::SurfaceMesh *&polyscopeMesh,
          mem3dg::solver::MutableTrajFile &fd, int &idx, int &waitTime,
          const Quantities options, bool &toggle);

/**
 * @brief Register Polyscope surface mesh from certain frame of the NetCDF
 * trajectory file with options of data quantities
 */
polyscope::SurfaceMesh *registerSurfaceMesh(mem3dg::solver::TrajFile &fd,
                                            int idx, const Quantities &options);
polyscope::SurfaceMesh *registerSurfaceMesh(mem3dg::solver::MutableTrajFile &fd,
                                            int idx, const Quantities &options);
#endif
