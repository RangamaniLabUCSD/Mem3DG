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

#include "mem3dg/mem3dg"

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

using EigenVectorX1D = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using EigenVectorX3D =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenTopVec =
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

void updateSurfaceMesh(ddgsolver::TrajFile &fd, int &idx) {
  if (idx >= fd.getNextFrameIndex()) {
    idx = 0;
  }

  double time;
  EigenVectorX3D coords;
  EigenTopVec top = fd.getTopology();
  EigenVectorX1D H = fd.getMeanCurvature(idx);
  std::tie(time, coords) = fd.getTimeAndCoords(idx);

  polyscope::registerSurfaceMesh("Mesh", coords, top);
  polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("mean_curvature",
                                                             H);
}

void animate(ddgsolver::TrajFile &fd, int &idx, bool &play) {
  if (play) {
    updateSurfaceMesh(fd, idx);
    idx++;
    if (idx >= fd.getNextFrameIndex()) {
      idx = 0;
    }
  }
}

int view_animation(std::string &filename) {

  ddgsolver::TrajFile fd = ddgsolver::TrajFile::openReadOnly(filename);

  int idx = 0;
  bool play = false;

  // Some settings for polyscope
  polyscope::options::programName = "Mem3DG Visualization";
  polyscope::options::verbosity = 0;
  polyscope::options::usePrefsFile = false;
  polyscope::options::autocenterStructures = false;
  polyscope::options::autoscaleStructures = false;

  polyscope::view::upDir = polyscope::view::UpDir::ZUp;
  polyscope::view::style = polyscope::view::NavigateStyle::Free;

  polyscope::init();

  auto myCallback = [&fd, &play, &idx]() {
    // Since options::openImGuiWindowForUserCallback == true by default,
    // we can immediately start using ImGui commands to build a UI
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
                               // instead of full width. Must have
                               // matching PopItemWidth() below.

    ImGui::InputInt("index", &idx); // set a float variable

    if (ImGui::Button("Play/Pause")) {
      play = !play;
    }

    if (ImGui::Button("Update Frame")) {
      updateSurfaceMesh(fd, idx);
    }

    animate(fd, idx, play);
    ImGui::PopItemWidth();
  };

  polyscope::state::userCallback = myCallback;

  polyscope::show();

  return 0;
}
