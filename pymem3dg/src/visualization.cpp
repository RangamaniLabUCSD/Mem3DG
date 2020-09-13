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

#ifdef MEM3DG_WITH_NETCDF
#include <iostream>
#include <time.h>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "mem3dg/mem3dg"

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

void wait(unsigned timeout)
{
  timeout += std::clock(); 
  while(std::clock() < timeout) continue;
}

using EigenVectorX1D = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using EigenVectorX3D =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenTopVec =
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

void updateSurfaceMesh(polyscope::SurfaceMesh *mesh, ddgsolver::TrajFile &fd,
                       int &idx) {
  if (idx >= fd.getNextFrameIndex()) {
    idx = 0;
  }

  double time;
  EigenVectorX3D coords;
  EigenVectorX3D refcoords = fd.getRefcoordinate();
  EigenVectorX3D vel = fd.getVelocity(idx);
  EigenVectorX1D H = fd.getMeanCurvature(idx);
  EigenVectorX1D H0 = fd.getSponCurvature(idx);
  EigenVectorX1D f_ext = fd.getExternalPressure(idx);
  EigenVectorX1D fn = fd.getPhysicalPressure(idx);
  EigenVectorX1D ft = fd.getCapillaryPressure(idx);
  EigenVectorX1D fb = fd.getBendingPressure(idx);
  std::tie(time, coords) = fd.getTimeAndCoords(idx);

  // polyscope::registerSurfaceMesh("Mesh", coords, top);
  mesh->updateVertexPositions(coords);
  mesh->addVertexVectorQuantity("velocity", vel);
  mesh->addVertexVectorQuantity("ref_coordinate", refcoords);
  mesh->addVertexScalarQuantity("mean_curvature", H);
  mesh->addVertexScalarQuantity("spon_curvature", H0);
  mesh->addVertexScalarQuantity("external_pressure", f_ext);
  mesh->addVertexScalarQuantity("physical_pressure", fn);
  mesh->addVertexScalarQuantity("capillary_pressure", ft);
  mesh->addVertexScalarQuantity("bending_pressure", fb);
}

polyscope::SurfaceMesh *registerSurfaceMesh(ddgsolver::TrajFile &fd) {
  double time;
  EigenVectorX3D coords;
  EigenVectorX3D refcoords = fd.getRefcoordinate();
  EigenVectorX3D vel = fd.getVelocity(0);
  EigenTopVec top = fd.getTopology();
  EigenVectorX1D H = fd.getMeanCurvature(0);
  EigenVectorX1D H0 = fd.getSponCurvature(0);
  EigenVectorX1D f_ext = fd.getExternalPressure(0);
  EigenVectorX1D fn = fd.getPhysicalPressure(0);
  EigenVectorX1D ft = fd.getCapillaryPressure(0);
  EigenVectorX1D fb = fd.getBendingPressure(0);
  std::tie(time, coords) = fd.getTimeAndCoords(0);

  polyscope::SurfaceMesh *mesh =
      polyscope::registerSurfaceMesh("Mesh", coords, top);
  polyscope::getSurfaceMesh("Mesh")->addVertexVectorQuantity("velocity",
                                                             vel);
  polyscope::getSurfaceMesh("Mesh")->addVertexVectorQuantity("ref_coordinate", refcoords);
  polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("mean_curvature",
                                                             H);
  polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("spon_curvature",
                                                             H0);
  polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("external_pressure",
                                                             f_ext);
  polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("physical_pressure",
                                                             fn);
  polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("capillary_pressure",
                                                             ft);
  polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("bending_pressure",
                                                             fb);
  return mesh;
}

void animate(polyscope::SurfaceMesh *mesh, ddgsolver::TrajFile &fd, int &idx, int &waitTime) {
  updateSurfaceMesh(mesh, fd, idx);
  idx++;
  if (idx >= fd.getNextFrameIndex()) {
    idx = 0;
  }
  wait(waitTime);
}

int view_animation(std::string &filename) {

  ddgsolver::TrajFile fd = ddgsolver::TrajFile::openReadOnly(filename);

  // Visualization state variables
  int prevFrame = 0;
  int currFrame = 0;
  bool play = false;
  int maxFrame = fd.getNextFrameIndex() - 1;
  int maxWaitTime = 500;
  int waitTime = 0;

  // Some settings for polyscope
  polyscope::options::programName = "Mem3DG Visualization";
  polyscope::options::verbosity = 0;
  polyscope::options::usePrefsFile = false;
  polyscope::options::autocenterStructures = false;
  polyscope::options::autoscaleStructures = false;

  polyscope::view::upDir = polyscope::view::UpDir::ZUp;
  polyscope::view::style = polyscope::view::NavigateStyle::Free;

  polyscope::init();

  auto mesh = registerSurfaceMesh(fd);

  auto myCallback = [&]() {
    // Since options::openImGuiWindowForUserCallback == true by default,
    // we can immediately start using ImGui commands to build a UI
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
                               // instead of full width. Must have
                               // matching PopItemWidth() below.

    ImGui::SliderInt("index", &currFrame, 0, maxFrame); // set a float variable

    ImGui::SliderInt("slow-mo", &waitTime, 0, maxWaitTime); // set a float variable

    if (ImGui::Button("Play/Pause")) {
      play = !play;
    }

    if (prevFrame != currFrame) {
      updateSurfaceMesh(mesh, fd, currFrame);
      prevFrame = currFrame;
    }

    if (play) {
      animate(mesh, fd, currFrame, waitTime);
      prevFrame = currFrame;
    }
    ImGui::PopItemWidth();
  };

  polyscope::state::userCallback = myCallback;

  polyscope::show();

  return 0;
}
#endif
