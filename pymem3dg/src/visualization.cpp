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

#include <iostream>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "mem3dg/solver/ddgsolver.h"
#include "mem3dg/solver/force.h"
#include "mem3dg/solver/icosphere.h"
#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/solver/typetraits.h"
#include "mem3dg/solver/util.h"

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;



void otherCallback(ddgsolver::TrajFile &fd) {
  double time;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> coords;

  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor> top =
      fd.getTopology();

  Eigen::Matrix<double, Eigen::Dynamic, 1> H = fd.getMeanCurvature(anim_index);
  std::tie(time, coords) = fd.getTimeAndCoords(anim_index++);

  if (anim_index >= fd.getNextFrameIndex()) {
    anim_index = 0;
  }

  polyscope::registerSurfaceMesh("Vesicle surface", coords, top);
  polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexScalarQuantity("mean_curvature", H);
}

int view_animation(std::string &filename) {
  ddgsolver::TrajFile fd = ddgsolver::TrajFile::openReadOnly(filename);

  polyscope::init();

  //   polyscope::getSurfaceMesh("Vesicle surface")
  //   ->addVertexScalarQuantity("mean_curvature", meanCurvature_e);

  auto myCallback = [&fd]() {
    // Since options::openImGuiWindowForUserCallback == true by default,
    // we can immediately start using ImGui commands to build a UI
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
                               // instead of full width. Must have
                               // matching PopItemWidth() below.

    ImGui::InputInt("index", &anim_index); // set a float variable

    if (ImGui::Button("Play/Pause")) {
      play = !play;
    }

    if (ImGui::Button("Rerender")) {
      otherCallback(fd);
    }

    mySubroutine(fd);
    ImGui::PopItemWidth();
  };

  polyscope::state::userCallback = myCallback;

  polyscope::show();

  return 0;
}
