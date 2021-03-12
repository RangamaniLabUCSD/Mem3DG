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
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <geometrycentral/surface/surface_mesh.h>

/**
 * @brief Quantities option for visualization
 *
 */
struct Quantities {
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
  /// vertex external pressure
  bool ext_pressure = false;
  /// vertex (total) physical pressure
  bool physical_pressure = false;
  /// vertex capillary pressure
  bool capillary_pressure = false;
  /// constant global inside pressure
  bool inside_pressure = false;
  /// vertex bending pressure
  bool bending_pressure = false;
  /// vertex line tension pressure
  bool line_pressure = false;
  /// vertex mask for time integration
  bool mask = false;
  /// vertex mean - spontaneous curvature difference
  bool H_H0 = false;
};

// ==========================================================
// =============        Viewers                ==============
// ==========================================================

/**
 * @brief visualize system object
 */
void visualize(mem3dg::System &f);

/**
 * @brief Visualize .ply file in polysope with options of additional quantities
 */
int snapshot_ply(std::string fileName, const Quantities &options,
                 float transparency = 1, float fov = 50, float edgeWidth = 1);

/**
 * @brief Visualize .ply files in polysope with options of additional quantities
 */
int animate_ply(std::string framesDir, const Quantities &options,
                std::vector<size_t> frameNum, float transparency = 1,
                float fov = 50, float edgeWidth = 1);

#ifdef MEM3DG_WITH_NETCDF
/**
 * @brief Visualize netcdf file in single frame
 */
int snapshot_nc(std::string &filename, const Quantities &options, int frame,
                float transparency = 1, float angle = 0, float fov = 50,
                float edgeWidth = 1, bool isShow = true, bool isSave = false,
                std::string screenshotName = "screenshot.png");

/**
 * @brief Animate netcdf file with options of additional quantities
 */
int animate_nc(std::string &filename, const Quantities &options,
               float transparency = 1, float fov = 50, float edgeWidth = 1);
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
          std::vector<size_t> frameNum);

/**
 * @brief Register Polyscope surface mesh from .ply file with options of data
 * quantities
 */
polyscope::SurfaceMesh *registerSurfaceMesh(std::string plyName,
                                            const Quantities &options);

#ifdef MEM3DG_WITH_NETCDF
/**
 * @brief Play the next frame of the NetCDF trajectory file
 */
void play(polyscope::SurfaceMesh *&polyscopeMesh, mem3dg::TrajFile &fd,
          int &idx, int &waitTime, const Quantities options, bool &toggle);

/**
 * @brief Register Polyscope surface mesh from certain frame of the NetCDF
 * trajectory file with options of data quantities
 */
polyscope::SurfaceMesh *registerSurfaceMesh(mem3dg::TrajFile &fd, int idx,
                                            const Quantities &options);
#endif