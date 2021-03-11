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
#include "mem3dg/solver/system.h"
#include "polyscope/polyscope.h"

#include "polyscope/surface_mesh.h"
#include <vector>

struct Quantities {
  bool ref_coord;
  bool velocity;
  bool mean_curvature;
  bool gauss_curvature;
  bool spon_curvature;
  bool ext_pressure;
  bool physical_pressure;
  bool capillary_pressure;
  bool inside_pressure;
  bool bending_pressure;
  bool line_pressure;
  bool mask;
  bool H_H0;
};

/**
 * @brief visualize system object
 *
 */
void visualize(mem3dg::System &f);

/**
 * @brief Visualize .ply file in polysope with options of additional quantities
 *
 */
int snapshot_ply(std::string fileName, const Quantities &options);

/**
 * @brief Visualize .ply files in polysope with options of additional quantities
 *
 */
int animate_ply(std::string framesDir, std::vector<size_t> frameNum,
                const Quantities &options);

#ifdef MEM3DG_WITH_NETCDF
/**
 * @brief Visualize netcdf file in single frame
 *
 */
int snapshot_nc(std::string &filename, int frame, const Quantities &options,
                float transparency, float angle, float fov, float edgeWidth,
                bool isShow, bool isSave, std::string screenshotName);

/**
 * @brief Animate netcdf file with options of additional quantities
 *
 */
int animate_nc(std::string &filename, const Quantities &options,
                 float transparency = 1, float angle = 0, float fov = 50,
                 float edgeWidth = 1);
#endif

void initGui();
void wait(unsigned timeout);
void play(mem3dg::TrajFile &fd, int &idx, int &waitTime, Quantities options,
          bool &toggle);
void play(std::string framesDir, int &idx, int &waitTime, Quantities options,
          bool &toggle, int maxFrame);
polyscope::SurfaceMesh *registerSurfaceMesh(mem3dg::TrajFile &fd, int &idx,
                                            const Quantities &options);
polyscope::SurfaceMesh *registerSurfaceMesh(mem3dg::TrajFile &fd, int &&idx,
                                            const Quantities &options);
polyscope::SurfaceMesh *registerSurfaceMesh(std::string plyName,
                                            const Quantities &options);
