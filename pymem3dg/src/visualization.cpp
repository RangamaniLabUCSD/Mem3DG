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

#include <csignal>
#include <iostream>
#include <time.h>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"

#include <geometrycentral/surface/meshio.h>

#include "mem3dg/mem3dg"

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

using EigenVectorX1D = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using EigenVectorX1D_i = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using EigenVectorX3D =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenTopVec =
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

void signalHandler(int signum);

int viewer_ply(std::string fileName, const bool mean_curvature,
               const bool spon_curvature, const bool ext_pressure,
               const bool physical_pressure, const bool capillary_pressure,
               const bool bending_pressure, const bool line_pressure) {

  signal(SIGINT, signalHandler);

  /// alias eigen matrix
  using EigenVectorX1D = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using EigenVectorX3D =
      Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
  using EigenTopVec =
      Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

  /// Declare pointers to mesh, geometry and richdata objects
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;

  /// Load input mesh
  if (!mean_curvature && !spon_curvature && !ext_pressure &&
      !physical_pressure && !capillary_pressure && !bending_pressure &&
      !line_pressure) {
    std::tie(ptrMesh, ptrVpg) = gcs::readManifoldSurfaceMesh(fileName);
  } else {
    std::tie(ptrMesh, ptrRichData) =
        gcs::RichSurfaceMeshData::readMeshAndData(fileName);
    ptrVpg = ptrRichData->getGeometry();
  }

  // Initialize visualization variables
  float transparency = 1;

  // Set preference for polyscope
  polyscope::options::programName = "Mem3DG Visualization";
  polyscope::options::verbosity = 0;
  polyscope::options::usePrefsFile = false;
  polyscope::options::autocenterStructures = false;
  polyscope::options::autoscaleStructures = false;
  polyscope::options::groundPlaneEnabled = false;
  polyscope::options::transparencyMode = polyscope::TransparencyMode::Pretty;
  polyscope::view::upDir = polyscope::view::UpDir::ZUp;
  polyscope::view::style = polyscope::view::NavigateStyle::Turntable;

  /// Initiate in polyscope
  polyscope::init();
  polyscope::SurfaceMesh *mesh = polyscope::registerSurfaceMesh(
      "Vesicle surface", ptrVpg->inputVertexPositions,
      ptrMesh->getFaceVertexList());
  mesh->setSmoothShade(true);
  mesh->setEnabled(true);
  mesh->setEdgeWidth(1);

  /// Read element data
  if (mean_curvature) {
    gcs::VertexData<double> meanCurvature =
        ptrRichData->getVertexProperty<double>("mean_curvature");
    EigenVectorX1D meanCurvature_e = meanCurvature.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("mean_curvature", meanCurvature_e);
  }
  if (spon_curvature) {
    gcs::VertexData<double> sponCurvature =
        ptrRichData->getVertexProperty<double>("spon_curvature");
    EigenVectorX1D sponCurvature_e = sponCurvature.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("spon_curvature", sponCurvature_e);
  }
  if (ext_pressure) {
    gcs::VertexData<double> extPressure =
        ptrRichData->getVertexProperty<double>("external_pressure");
    EigenVectorX1D extPressure_e = extPressure.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("applied_pressure", extPressure_e);
  }
  if (physical_pressure) {
    gcs::VertexData<double> physicalPressure =
        ptrRichData->getVertexProperty<double>("physical_pressure");
    EigenVectorX1D physicalPressure_e = physicalPressure.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("physical_pressure", physicalPressure_e);
  }
  if (capillary_pressure) {
    gcs::VertexData<double> capillaryPressure =
        ptrRichData->getVertexProperty<double>("capillary_pressure");
    EigenVectorX1D capillaryPressure_e = capillaryPressure.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("surface_tension", capillaryPressure_e);
  }
  if (bending_pressure) {
    gcs::VertexData<double> bendingPressure =
        ptrRichData->getVertexProperty<double>("bending_pressure");
    EigenVectorX1D bendingPressure_e = bendingPressure.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("bending_pressure", bendingPressure_e);
  }
  if (line_pressure) {
    gcs::VertexData<double> linePressure =
        ptrRichData->getVertexProperty<double>("line_tension_pressure");
    EigenVectorX1D linePressure_e = linePressure.raw();
    polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexScalarQuantity("line_tension_pressure", linePressure_e);
  }

  /*gcs::VertexData<gc::Vector3> vertexVelocity =
      ptrRichData->getVertexProperty<gc::Vector3>("vertex_velocity");*/
  /*gcs::VertexData<gc::Vector3> normalForce =
  ptrRichData->getVertexProperty<gc::Vector3>("normal_force");
  gcs::VertexData<gc::Vector3> tangentialForce =
  ptrRichData->getVertexProperty<gc::Vector3>("tangential_force");*/
  // EigenVectorX3D vertexVelocity_e =
  //    mem3dg::EigenMap<double, 3>(vertexVelocity);
  /*EigenVectorX3D normalForce_e =
  gc::EigenMap<double, 3>(normalForce); Eigen::Matrix<double,
  Eigen::Dynamic, 3> tangentialForce_e = gc::EigenMap<double,
  3>(tangentialForce);*/
  /*polyscope::getSurfaceMesh("Vesicle surface")
      ->addVertexVectorQuantity("vertexVelocity", vertexVelocity_e);*/
  /*polyscope::getSurfaceMesh("Vesicle
  surface")->addVertexVectorQuantity("tangential_force", tangentialForce_e);
  polyscope::getSurfaceMesh("Vesicle
  surface")->addVertexVectorQuantity("normal_force", normalForce_e);*/

  // Callback function for interactive GUI
  auto myCallback = [&]() {
    // Since options::openImGuiWindowForUserCallback == true by default,
    // we can immediately start using ImGui commands to build a UI
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
                               // instead of full width. Must have
                               // matching PopItemWidth() below.

    // Initialize sliders
    ImGui::SliderFloat("transparency", &transparency, 0, 1);

    // Execute transparency slider
    mesh->setTransparency(transparency);

    // Define buttons
    if (ImGui::Button("Screenshot")) {
      // char buff[50];
      // snprintf(buff, 50, "screenshot_frame%06d.png", frame);
      // std::string defaultName(buff);
      polyscope::screenshot("screenshot.png", true);
    }

    ImGui::PopItemWidth();
  };
  polyscope::state::userCallback = myCallback;
  polyscope::show();

  return 0;
}

#ifdef MEM3DG_WITH_NETCDF
void wait(unsigned timeout) {
  timeout += std::clock();
  while (std::clock() < timeout)
    continue;
}

struct checkBox {
  const bool ref_coord;
  const bool velocity;
  const bool mean_curvature;
  const bool spon_curvature;
  const bool ext_pressure;
  const bool physical_pressure;
  const bool capillary_pressure;
  const bool bending_pressure;
  const bool line_pressure;
  const bool mask;
  const bool H_H0;
};

void updateSurfaceMesh(polyscope::SurfaceMesh *mesh, mem3dg::TrajFile &fd,
                       int &idx, checkBox options) {

  if (idx >= fd.getNextFrameIndex()) {
    idx = 0;
  }

  double time;
  EigenVectorX3D coords;
  std::tie(time, coords) = fd.getTimeAndCoords(idx);
  mesh->updateVertexPositions(coords);

  if (options.ref_coord) {
    EigenVectorX3D refcoords = fd.getRefcoordinate();
    mesh->addVertexVectorQuantity("ref_coordinate", refcoords);
  }

  if (options.velocity) {
    EigenVectorX3D vel = fd.getVelocity(idx);
    mesh->addVertexVectorQuantity("velocity", vel);
  }
  if (options.mean_curvature) {
    EigenVectorX1D H = fd.getMeanCurvature(idx);
    mesh->addVertexScalarQuantity("mean_curvature", H);
  }
  if (options.spon_curvature) {
    EigenVectorX1D H0 = fd.getSponCurvature(idx);
    mesh->addVertexScalarQuantity("spon_curvature", H0);
  }
  if (options.ext_pressure) {
    EigenVectorX1D f_ext = fd.getExternalPressure(idx);
    mesh->addVertexScalarQuantity("external_pressure", f_ext);
  }
  if (options.physical_pressure) {
    EigenVectorX1D fn = fd.getPhysicalPressure(idx);
    mesh->addVertexScalarQuantity("physical_pressure", fn);
  }
  if (options.capillary_pressure) {
    EigenVectorX1D ft = fd.getCapillaryPressure(idx);
    mesh->addVertexScalarQuantity("capillary_pressure", ft);
  }
  if (options.bending_pressure) {
    EigenVectorX1D fb = fd.getBendingPressure(idx);
    mesh->addVertexScalarQuantity("bending_pressure", fb);
  }
  if (options.line_pressure) {
    EigenVectorX1D fl = fd.getLinePressure(idx);
    mesh->addVertexScalarQuantity("line_tension_pressure", fl);
  }
  if (options.mask) {
    EigenVectorX1D_i msk = fd.getMask();
    mesh->addVertexScalarQuantity("mask", msk);
  }
  if (options.H_H0) {
    EigenVectorX1D h_h0 = fd.getH_H0_diff(idx);
    mesh->addVertexScalarQuantity("curvature_diff", h_h0);
  }
  // polyscope::registerSurfaceMesh("Mesh", coords, top);
}

polyscope::SurfaceMesh *registerSurfaceMesh(mem3dg::TrajFile &fd,
                                            checkBox options) {
  double time;
  EigenVectorX3D coords;
  EigenTopVec top = fd.getTopology();
  std::tie(time, coords) = fd.getTimeAndCoords(0);
  polyscope::SurfaceMesh *mesh =
      polyscope::registerSurfaceMesh("Mesh", coords, top);
  mesh->setSmoothShade(true);
  mesh->setEnabled(true);

  if (options.ref_coord) {
    EigenVectorX3D refcoords = fd.getRefcoordinate();
    polyscope::getSurfaceMesh("Mesh")->addVertexVectorQuantity("ref_coordinate",
                                                               refcoords);
  }

  if (options.velocity) {
    EigenVectorX3D vel = fd.getVelocity(0);
    polyscope::getSurfaceMesh("Mesh")->addVertexVectorQuantity("velocity", vel);
  }
  if (options.mean_curvature) {
    EigenVectorX1D H = fd.getMeanCurvature(0);
    polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("mean_curvature",
                                                               H);
  }
  if (options.spon_curvature) {
    EigenVectorX1D H0 = fd.getSponCurvature(0);
    polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("spon_curvature",
                                                               H0);
  }
  if (options.ext_pressure) {
    EigenVectorX1D f_ext = fd.getExternalPressure(0);
    polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity(
        "external_pressure", f_ext);
  }
  if (options.physical_pressure) {
    EigenVectorX1D fn = fd.getPhysicalPressure(0);
    polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity(
        "physical_pressure", fn);
  }
  if (options.capillary_pressure) {
    EigenVectorX1D ft = fd.getCapillaryPressure(0);
    polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity(
        "capillary_pressure", ft);
  }
  if (options.bending_pressure) {
    EigenVectorX1D fb = fd.getBendingPressure(0);
    polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity(
        "bending_pressure", fb);
  }
  if (options.line_pressure) {
    EigenVectorX1D fb = fd.getLinePressure(0);
    polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity(
        "line_tension_pressure", fb);
  }
  if (options.mask) {
    EigenVectorX1D_i msk = fd.getMask();
    polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("mask", msk);
  }
  if (options.H_H0) {
    EigenVectorX1D h_h0 = fd.getH_H0_diff(0);
    polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("curvature_diff",
                                                               h_h0);
  }
  return mesh;
}

void animate(polyscope::SurfaceMesh *mesh, mem3dg::TrajFile &fd, int &idx,
             int &waitTime, checkBox options, bool &toggle) {

  updateSurfaceMesh(mesh, fd, idx, options);
  idx++;
  if (idx >= fd.getNextFrameIndex()) {
    idx = 0;
    toggle = !toggle;
  }
  wait(waitTime);
}

int animation_nc(std::string &filename, const bool ref_coord,
                 const bool velocity, const bool mean_curvature,
                 const bool spon_curvature, const bool ext_pressure,
                 const bool physical_pressure, const bool capillary_pressure,
                 const bool bending_pressure, const bool line_pressure,
                 const bool mask, const bool H_H0) {

  // Activate signal handling
  signal(SIGINT, signalHandler);

  // Read netcdf trajectory file
  mem3dg::TrajFile fd = mem3dg::TrajFile::openReadOnly(filename);

  // Initialize visualization variables
  int prevFrame = 0;
  int currFrame = 0;
  bool play = false;
  bool record = false;
  int maxFrame = fd.getNextFrameIndex() - 1;
  int maxWaitTime = 500;
  int waitTime = 0;
  float transparency = 1;
  checkBox options({ref_coord, velocity, mean_curvature, spon_curvature,
                    ext_pressure, physical_pressure, capillary_pressure,
                    bending_pressure, line_pressure, mask, H_H0});

  // Set preference for polyscope
  polyscope::options::programName = "Mem3DG Visualization";
  polyscope::options::verbosity = 0;
  polyscope::options::usePrefsFile = false;
  polyscope::options::autocenterStructures = false;
  polyscope::options::autoscaleStructures = false;
  polyscope::options::groundPlaneEnabled = false;
  polyscope::options::transparencyMode = polyscope::TransparencyMode::Pretty;
  polyscope::view::upDir = polyscope::view::UpDir::ZUp;
  polyscope::view::style = polyscope::view::NavigateStyle::Turntable;

  // Initialize polyscope
  polyscope::init();

  // Initialize surface mesh
  auto mesh = registerSurfaceMesh(fd, options);

  // Callback function for interactive GUI
  auto myCallback = [&]() {
    // Since options::openImGuiWindowForUserCallback == true by default,
    // we can immediately start using ImGui commands to build a UI
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
                               // instead of full width. Must have
                               // matching PopItemWidth() below.

    // Initialize sliders
    ImGui::SliderInt("index", &currFrame, 0, maxFrame); // set a float variable
    ImGui::SliderFloat("transparency", &transparency, 0, 1);
    ImGui::SliderInt("slow-mo", &waitTime, 0,
                     maxWaitTime); // set a float variable

    // Execute transparency slider
    mesh->setTransparency(transparency);

    // Define buttons
    if (ImGui::Button("Play/Pause")) {
      play = !play;
    }
    if (ImGui::Button("Screenshot")) {
      char buff[50];
      snprintf(buff, 50, "screenshot_frame%06d.png", currFrame);
      std::string defaultName(buff);
      polyscope::screenshot(defaultName, true);
      // polyscope::screenshot("screenshot.png", true);
    }
    if (ImGui::Button("Record")) {
      record = !record;
    }
    if (prevFrame != currFrame) {
      updateSurfaceMesh(mesh, fd, currFrame, options);
      prevFrame = currFrame;
    }
    if (record) {
      char buff[50];
      snprintf(buff, 50, "video/frame%06d.png", currFrame);
      std::string defaultName(buff);
      polyscope::screenshot(defaultName, true);
      animate(mesh, fd, currFrame, waitTime, options, record);
      prevFrame = currFrame;
    }
    if (play) {
      animate(mesh, fd, currFrame, waitTime, options, play);
      prevFrame = currFrame;
    }

    ImGui::PopItemWidth();
  };

  polyscope::state::userCallback = myCallback;
  polyscope::show();

  return 0;
}

int snapshot_nc(std::string &filename, int frame, float angle, bool isShow,
                bool isSave, std::string screenshotName, const bool ref_coord,
                const bool velocity, const bool mean_curvature,
                const bool spon_curvature, const bool ext_pressure,
                const bool physical_pressure, const bool capillary_pressure,
                const bool bending_pressure, const bool line_pressure,
                const bool mask, const bool H_H0) {

  static int ENTRY = 0;

  // Activate signal handling
  signal(SIGINT, signalHandler);

  // Read netcdf trajectory file
  mem3dg::TrajFile fd = mem3dg::TrajFile::openReadOnly(filename);

  // Initialize visualization variables
  int maxFrame = fd.getNextFrameIndex() - 1;
  if (frame > maxFrame) {
    throw std::runtime_error("Snapshot frame exceed maximum frame index!");
  }
  float transparency = 1;
  checkBox options({ref_coord, velocity, mean_curvature, spon_curvature,
                    ext_pressure, physical_pressure, capillary_pressure,
                    bending_pressure, line_pressure, mask, H_H0});

  // Set preference for polyscope
  polyscope::options::programName = "Mem3DG Visualization";
  polyscope::options::verbosity = 0;
  polyscope::options::usePrefsFile = false;
  polyscope::options::autocenterStructures = false;
  polyscope::options::autoscaleStructures = false;
  polyscope::options::groundPlaneEnabled = false;
  polyscope::options::transparencyMode = polyscope::TransparencyMode::Pretty;
  polyscope::view::upDir = polyscope::view::UpDir::ZUp;
  polyscope::view::style = polyscope::view::NavigateStyle::Turntable;
  // polyscope::view::fov = 50;

  // Initialize polyscope
  if (!polyscope::state::initialized) {
    polyscope::init();
  }

  // Initialize surface mesh and switch to specific frame
  auto mesh = registerSurfaceMesh(fd, options);
  updateSurfaceMesh(mesh, fd, frame, options);
  mesh->setEdgeWidth(1);

  if (options.ref_coord) {
    EigenVectorX3D refcoords = fd.getRefcoordinate();
    mesh->addVertexVectorQuantity("ref_coordinate", refcoords)
        ->setEnabled(true);
  }
  if (options.velocity) {
    EigenVectorX3D vel = fd.getVelocity(frame);
    mesh->addVertexVectorQuantity("velocity", vel)->setEnabled(true);
  }
  if (options.mean_curvature) {
    EigenVectorX1D H = fd.getMeanCurvature(frame);
    mesh->addVertexScalarQuantity("mean_curvature", H)->setEnabled(true);
  }
  if (options.spon_curvature) {
    EigenVectorX1D H0 = fd.getSponCurvature(frame);
    mesh->addVertexScalarQuantity("spon_curvature", H0)->setEnabled(true);
  }
  if (options.ext_pressure) {
    EigenVectorX1D f_ext = fd.getExternalPressure(frame);
    mesh->addVertexScalarQuantity("external_pressure", f_ext)->setEnabled(true);
  }
  if (options.physical_pressure) {
    EigenVectorX1D fn = fd.getPhysicalPressure(frame);
    mesh->addVertexScalarQuantity("physical_pressure", fn)->setEnabled(true);
  }
  if (options.capillary_pressure) {
    EigenVectorX1D ft = fd.getCapillaryPressure(frame);
    mesh->addVertexScalarQuantity("capillary_pressure", ft)->setEnabled(true);
  }
  if (options.bending_pressure) {
    EigenVectorX1D fb = fd.getBendingPressure(frame);
    mesh->addVertexScalarQuantity("bending_pressure", fb)->setEnabled(true);
  }
  if (options.line_pressure) {
    EigenVectorX1D fl = fd.getLinePressure(frame);
    mesh->addVertexScalarQuantity("line_tension_pressure", fl)
        ->setEnabled(true);
  }
  if (options.mask) {
    EigenVectorX1D_i msk = fd.getMask();
    mesh->addVertexScalarQuantity("mask", msk)->setEnabled(true);
  }
  if (options.H_H0) {
    EigenVectorX1D h_h0 = fd.getH_H0_diff(frame);
    mesh->addVertexScalarQuantity("curvature_diff", h_h0)->setEnabled(true);
  }

  // Callback function for interactive GUI
  auto myCallback = [&]() {
    // Since options::openImGuiWindowForUserCallback == true by default,
    // we can immediately start using ImGui commands to build a UI
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
                               // instead of full width. Must have
                               // matching PopItemWidth() below.

    // Initialize sliders
    ImGui::SliderFloat("transparency", &transparency, 0, 1);

    // Execute transparency slider
    mesh->setTransparency(transparency);

    // Define buttons
    if (ImGui::Button("Screenshot")) {
      // char buff[50];
      // snprintf(buff, 50, "screenshot_frame%06d.png", frame);
      // std::string defaultName(buff);
      polyscope::screenshot(screenshotName, true);
    }

    ImGui::PopItemWidth();
  };

  // For some reasons we have to run screenshot twice for fov (field of view)
  // setting to work, not ideal
  if (isSave) {
    polyscope::screenshot(screenshotName, true);
  }
  polyscope::state::userCallback = myCallback;

  polyscope::view::fov = 50;
  if (ENTRY == 0) {
    // glm::vec3 frameLookDir, frameUpDir, frameRightDir;
    // polyscope::view::getCameraFrame(frameLookDir, frameUpDir, frameRightDir);
    glm::vec3 frameRightDir = glm::vec3(1., 0, 0.);
    glm::mat4x4 phiCamR = glm::rotate(glm::mat4x4(1.0), angle, frameRightDir);
    polyscope::view::viewMat = polyscope::view::viewMat * phiCamR;
  }

  if (isSave) {
    polyscope::screenshot(screenshotName, true);
  }
  if (isShow) {
    polyscope::show();
  }

  ENTRY++;
  return 0;
}
#endif
