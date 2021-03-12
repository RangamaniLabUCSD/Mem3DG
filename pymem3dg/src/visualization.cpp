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
#include <pybind11/embed.h>

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

void visualize(mem3dg::System &f) {
  signal(SIGINT, mem3dg::signalHandler);
  // Initialize visualization variables
  float transparency = 1;

  // Set preference for polyscope
  initGui();

  /// Initiate in polyscope
  polyscope::init();

  polyscope::SurfaceMesh *polyMesh = polyscope::registerSurfaceMesh(
      "Membrane", f.vpg->inputVertexPositions, f.mesh->getFaceVertexList());
  polyMesh->setSmoothShade(true);
  polyMesh->setEnabled(true);
  polyMesh->setEdgeWidth(1);

  // Process attributes
  Eigen::Matrix<double, Eigen::Dynamic, 1> fn;

  fn = f.bendingPressure.raw() + f.capillaryPressure.raw() +
       f.insidePressure.raw() + f.externalPressure.raw() +
       f.M_inv * f.lineCapillaryForce.raw();

  /// Read element data
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("mean_curvature", f.H.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("gauss_curvature", f.K.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("spon_curvature", f.H0.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("external_pressure", f.externalPressure.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("bending_pressure", f.bendingPressure.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("line_tension_pressure",
                                f.M_inv * f.lineCapillaryForce.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("capillary_pressure",
                                f.capillaryPressure.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("inside_pressure", f.insidePressure.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("physical_pressure", fn);
  polyscope::getSurfaceMesh("Membrane")
      ->addEdgeScalarQuantity("line_tension", f.lineTension.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addEdgeScalarQuantity("edge_dihedral", f.vpg->edgeDihedralAngles.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addEdgeScalarQuantity(
          "edge_line_capillary",
          f.vpg->hodge1Inverse *
              ((f.vpg->hodge1 *
                (f.lineTension.raw().array() / f.vpg->edgeLengths.raw().array())
                    .matrix())
                   .array() *
               f.vpg->edgeDihedralAngles.raw().array().max(0))
                  .matrix());
  // polyscope::getSurfaceMesh("Membrane")
  //     ->addEdgeScalarQuantity("isFlip", isFlip.raw().cast<double>());

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
    polyMesh->setTransparency(transparency);

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
}

int snapshot_ply(std::string fileName, const Quantities &option) {

  signal(SIGINT, mem3dg::signalHandler);

  auto polyscopeMesh = registerSurfaceMesh(fileName, option);

  // Initialize visualization variables
  float transparency = 1;

  // Set preference for polyscope
  initGui();

  /// Initiate in polyscope
  polyscope::init();

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
    polyscopeMesh->setTransparency(transparency);

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

int animate_ply(std::string frameDir, std::vector<size_t> frameNum,
                const Quantities &options) {

  // Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);

  // Initialize visualization variables
  int prevFrame = frameNum[0];
  int currFrame = frameNum[0];
  bool isStart = false;
  bool isRecord = false;
  int maxFrame = frameNum[1];
  int maxWaitTime = 500;
  int waitTime = 0;

  // Set preference for polyscope
  initGui();

  // Initialize polyscope
  polyscope::init();

  // Initialize surface mesh
  char buffer[50];
  sprintf(buffer, "/frame%d.ply", (int)frameNum[0]);
  std::string plyName(buffer);
  plyName = frameDir + plyName;
  auto polyscopeMesh = registerSurfaceMesh(plyName, options);
  polyscopeMesh->setSmoothShade(true);

  // Callback function for interactive GUI
  auto myCallback = [&]() {
    // Since options::openImGuiWindowForUserCallback == true by default,
    // we can immediately start using ImGui commands to build a UI
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
                               // instead of full width. Must have
                               // matching PopItemWidth() below.

    // Initialize sliders
    ImGui::SliderInt("index", &currFrame, frameNum[0],
                     maxFrame); // set a float variable
    ImGui::SliderInt("slow-mo", &waitTime, 0,
                     maxWaitTime); // set a float variable

    // Define buttons
    if (ImGui::Button("Play/Pause")) {
      isStart = !isStart;
    }
    if (ImGui::Button("Screenshot")) {
      char buff[50];
      snprintf(buff, 50, "screenshot_frame%06d.png", currFrame);
      std::string defaultName(buff);
      polyscope::screenshot(defaultName, true);
      // polyscope::screenshot("screenshot.png", true);
    }
    if (ImGui::Button("Record")) {
      isRecord = !isRecord;
    }
    if (prevFrame != currFrame) {
      char buffer[50];
      sprintf(buffer, "/frame%d.ply", (int)currFrame);
      std::string plyName(buffer);
      plyName = frameDir + plyName;
      polyscopeMesh = registerSurfaceMesh(plyName, options);
      prevFrame = currFrame;
    }
    if (isRecord) {
      char buff[50];
      snprintf(buff, 50, "video/frame%06d.png", currFrame);
      std::string defaultName(buff);
      polyscope::screenshot(defaultName, true);
      play(frameDir, currFrame, waitTime, options, isRecord, frameNum);
      prevFrame = currFrame;
    }
    if (isStart) {
      play(frameDir, currFrame, waitTime, options, isStart, frameNum);
      prevFrame = currFrame;
    }

    ImGui::PopItemWidth();
  };

  polyscope::state::userCallback = myCallback;
  polyscope::show();

  return 0;
}

#ifdef MEM3DG_WITH_NETCDF

int animate_nc(std::string &filename, const Quantities &options,
               float transparency, float angle, float fov, float edgeWidth) {

  // Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);

  // Read netcdf trajectory file
  mem3dg::TrajFile fd = mem3dg::TrajFile::openReadOnly(filename);

  // Initialize visualization variables
  int prevFrame = 0;
  int currFrame = 0;
  bool isStart = false;
  bool isRecord = false;
  int maxFrame = fd.getNextFrameIndex() - 1;
  int maxWaitTime = 500;
  int waitTime = 0;

  // Set preference for polyscope
  initGui();
  polyscope::view::fov = fov;

  // Initialize polyscope
  polyscope::init();

  // Initialize surface mesh
  auto polyscopeMesh = registerSurfaceMesh(fd, 0, options);
  polyscopeMesh->setSmoothShade(true);
  // polyscopeMesh->setEdgeWidth(edgeWidth);
  // polyscopeMesh->setTransparency(transparency);

  // Callback function for interactive GUI
  auto myCallback = [&]() {
    // Since options::openImGuiWindowForUserCallback == true by default,
    // we can immediately start using ImGui commands to build a UI
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
                               // instead of full width. Must have
                               // matching PopItemWidth() below.

    // Initialize sliders
    ImGui::SliderInt("index", &currFrame, 0,
                     maxFrame); // set a float variable
    ImGui::SliderFloat("transparency", &transparency, 0, 1);
    ImGui::SliderInt("slow-mo", &waitTime, 0,
                     maxWaitTime); // set a float variable

    // // Execute transparency slider
    // polyscopeMesh->setTransparency(transparency);

    // Define buttons
    if (ImGui::Button("Play/Pause")) {
      isStart = !isStart;
    }
    if (ImGui::Button("Screenshot")) {
      char buff[50];
      snprintf(buff, 50, "screenshot_frame%06d.png", currFrame);
      std::string defaultName(buff);
      polyscope::screenshot(defaultName, true);
      // polyscope::screenshot("screenshot.png", true);
    }
    if (ImGui::Button("Record")) {
      isRecord = !isRecord;
    }
    if (prevFrame != currFrame) {
      polyscopeMesh = registerSurfaceMesh(fd, currFrame, options);
      prevFrame = currFrame;
    }
    if (isRecord) {
      char buff[50];
      snprintf(buff, 50, "video/frame%06d.png", currFrame);
      std::string defaultName(buff);
      polyscope::screenshot(defaultName, true);
      play(fd, currFrame, waitTime, options, isRecord);
      prevFrame = currFrame;
    }
    if (isStart) {
      play(fd, currFrame, waitTime, options, isStart);
      prevFrame = currFrame;
    }

    ImGui::PopItemWidth();
  };

  polyscope::state::userCallback = myCallback;
  polyscope::show();

  return 0;
}

int snapshot_nc(std::string &filename, const Quantities &options, int frame,
                float transparency, float angle, float fov, float edgeWidth,
                bool isShow, bool isSave, std::string screenshotName) {

  static int ENTRY = 0;

  // Activate signal handling
  signal(SIGINT, mem3dg::signalHandler);

  // Read netcdf trajectory file
  mem3dg::TrajFile fd = mem3dg::TrajFile::openReadOnly(filename);
  fd.getNcFrame(frame);

  // Set preference for polyscope
  initGui();
  // polyscope::view::fov = 50;

  // Initialize polyscope
  if (!polyscope::state::initialized) {
    polyscope::init();
  }

  // Initialize surface mesh and switch to specific frame
  auto polyscopeMesh = registerSurfaceMesh(fd, frame, options);
  polyscopeMesh->setEdgeWidth(edgeWidth);
  polyscopeMesh->setTransparency(transparency);
  polyscopeMesh->setSmoothShade(true);

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
    polyscopeMesh->setTransparency(transparency);

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

  polyscope::view::fov = fov;
  if (ENTRY == 0) {
    // glm::vec3 frameLookDir, frameUpDir, frameRightDir;
    // polyscope::view::getCameraFrame(frameLookDir, frameUpDir,
    // frameRightDir);
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

#ifdef MEM3DG_WITH_NETCDF
polyscope::SurfaceMesh *registerSurfaceMesh(mem3dg::TrajFile &fd, int idx,
                                            const Quantities &options) {
  if (idx >= fd.getNextFrameIndex()) {
    idx = 0;
  }

  EigenTopVec topo_frame = fd.getTopoFrame(idx);
  EigenVectorX3D coords = fd.getCoords(idx);
  // mesh->updateVertexPositions(coords);
  polyscope::SurfaceMesh *polyscopeMesh =
      polyscope::registerSurfaceMesh("Mesh", coords, topo_frame);
  polyscopeMesh->setEnabled(true);

  if (options.ref_coord) {
    EigenVectorX3D refcoords = fd.getRefcoordinate();
    polyscopeMesh->addVertexVectorQuantity("ref_coordinate", refcoords);

    // Show quantities at the opening
    // polyscopeMesh->addVertexVectorQuantity("ref_coordinate", refcoords)
    //     ->setEnabled(true);
  }
  if (options.mean_curvature) {
    EigenVectorX1D H = fd.getMeanCurvature(idx);
    polyscopeMesh->addVertexScalarQuantity("mean_curvature", H);
  }
  if (options.gauss_curvature) {
    EigenVectorX1D K = fd.getGaussCurvature(idx);
    polyscopeMesh->addVertexScalarQuantity("gauss_curvature", K);
  }
  if (options.spon_curvature) {
    EigenVectorX1D H0 = fd.getSponCurvature(idx);
    polyscopeMesh->addVertexScalarQuantity("spon_curvature", H0);
  }
  if (options.mask) {
    EigenVectorX1D_i msk = fd.getMask();
    polyscopeMesh->addVertexScalarQuantity("mask", msk);
  }
  if (options.H_H0) {
    EigenVectorX1D h_h0 = fd.getH_H0_diff(idx);
    polyscopeMesh->addVertexScalarQuantity("curvature_diff", h_h0);
  }
  if (options.velocity) {
    EigenVectorX3D vel = fd.getVelocity(idx);
    polyscopeMesh->addVertexVectorQuantity("velocity", vel);
  }
  if (options.bending_pressure) {
    EigenVectorX1D fb = fd.getBendingPressure(idx);
    polyscopeMesh->addVertexScalarQuantity("bending_pressure", fb);
  }
  if (options.capillary_pressure) {
    EigenVectorX1D ft = fd.getCapillaryPressure(idx);
    polyscopeMesh->addVertexScalarQuantity("capillary_pressure", ft);
  }
  if (options.inside_pressure) {
    EigenVectorX1D ft = fd.getInsidePressure(idx);
    polyscopeMesh->addVertexScalarQuantity("inside_pressure", ft);
  }
  if (options.ext_pressure) {
    EigenVectorX1D f_ext = fd.getExternalPressure(idx);
    polyscopeMesh->addVertexScalarQuantity("external_pressure", f_ext);
  }
  if (options.line_pressure) {
    EigenVectorX1D fl = fd.getLinePressure(idx);
    polyscopeMesh->addVertexScalarQuantity("line_tension_pressure", fl);
  }
  if (options.physical_pressure) {
    EigenVectorX1D fn = fd.getPhysicalPressure(idx);
    polyscopeMesh->addVertexScalarQuantity("physical_pressure", fn);
  }
  // polyscope::registerSurfaceMesh("Mesh", coords, top);

  return polyscopeMesh;
}

void play(mem3dg::TrajFile &fd, int &idx, int &waitTime, Quantities options,
          bool &toggle) {

  registerSurfaceMesh(fd, idx, options);
  idx++;
  if (idx >= fd.getNextFrameIndex()) {
    idx = 0;
    toggle = !toggle;
  }
  wait(waitTime);
}
#endif

void wait(unsigned timeout) {
  timeout += std::clock();
  while (std::clock() < timeout)
    continue;
}

void initGui() {
  polyscope::options::programName = "Mem3DG Visualization";
  polyscope::options::verbosity = 0;
  polyscope::options::usePrefsFile = false;
  polyscope::options::autocenterStructures = false;
  polyscope::options::autoscaleStructures = false;
  polyscope::options::groundPlaneEnabled = false;
  polyscope::options::transparencyMode = polyscope::TransparencyMode::Pretty;
  polyscope::view::upDir = polyscope::view::UpDir::ZUp;
  polyscope::view::style = polyscope::view::NavigateStyle::Turntable;
}

polyscope::SurfaceMesh *registerSurfaceMesh(std::string plyName,
                                            const Quantities &options) {

  /// Declare pointers to mesh, geometry and richdata objects
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;

  /// Load input mesh
  if (!options.mean_curvature && !options.gauss_curvature &&
      !options.spon_curvature && !options.ext_pressure &&
      !options.physical_pressure && !options.capillary_pressure &&
      !options.bending_pressure && !options.line_pressure) {
    std::tie(ptrMesh, ptrVpg) = gcs::readManifoldSurfaceMesh(plyName);
  } else {
    std::tie(ptrMesh, ptrRichData) =
        gcs::RichSurfaceMeshData::readMeshAndData(plyName);
    ptrVpg = ptrRichData->getGeometry();
  }

  polyscope::SurfaceMesh *polyscopeMesh = polyscope::registerSurfaceMesh(
      "Vesicle surface", ptrVpg->inputVertexPositions,
      ptrMesh->getFaceVertexList());
  polyscopeMesh->setSmoothShade(true);
  polyscopeMesh->setEnabled(true);
  polyscopeMesh->setEdgeWidth(1);

  /// Read element data
  if (options.mean_curvature) {
    gcs::VertexData<double> meanCurvature =
        ptrRichData->getVertexProperty<double>("mean_curvature");
    EigenVectorX1D meanCurvature_e = meanCurvature.raw();
    polyscopeMesh->addVertexScalarQuantity("mean_curvature", meanCurvature_e);
  }
  if (options.gauss_curvature) {
    gcs::VertexData<double> gaussCurvature =
        ptrRichData->getVertexProperty<double>("gauss_curvature");
    EigenVectorX1D gaussCurvature_e = gaussCurvature.raw();
    polyscopeMesh->addVertexScalarQuantity("gauss_curvature", gaussCurvature_e);
  }
  if (options.spon_curvature) {
    gcs::VertexData<double> sponCurvature =
        ptrRichData->getVertexProperty<double>("spon_curvature");
    EigenVectorX1D sponCurvature_e = sponCurvature.raw();
    polyscopeMesh->addVertexScalarQuantity("spon_curvature", sponCurvature_e);
  }
  if (options.ext_pressure) {
    gcs::VertexData<double> extPressure =
        ptrRichData->getVertexProperty<double>("external_pressure");
    EigenVectorX1D extPressure_e = extPressure.raw();
    polyscopeMesh->addVertexScalarQuantity("applied_pressure", extPressure_e);
  }
  if (options.physical_pressure) {
    gcs::VertexData<double> physicalPressure =
        ptrRichData->getVertexProperty<double>("physical_pressure");
    EigenVectorX1D physicalPressure_e = physicalPressure.raw();
    polyscopeMesh->addVertexScalarQuantity("physical_pressure",
                                           physicalPressure_e);
  }
  if (options.capillary_pressure) {
    gcs::VertexData<double> capillaryPressure =
        ptrRichData->getVertexProperty<double>("capillary_pressure");
    EigenVectorX1D capillaryPressure_e = capillaryPressure.raw();
    polyscopeMesh->addVertexScalarQuantity("surface_tension",
                                           capillaryPressure_e);
  }
  if (options.bending_pressure) {
    gcs::VertexData<double> bendingPressure =
        ptrRichData->getVertexProperty<double>("bending_pressure");
    EigenVectorX1D bendingPressure_e = bendingPressure.raw();
    polyscopeMesh->addVertexScalarQuantity("bending_pressure",
                                           bendingPressure_e);
  }
  if (options.line_pressure) {
    gcs::VertexData<double> linePressure =
        ptrRichData->getVertexProperty<double>("line_tension_pressure");
    EigenVectorX1D linePressure_e = linePressure.raw();
    polyscopeMesh->addVertexScalarQuantity("line_tension_pressure",
                                           linePressure_e);
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

  return polyscopeMesh;
}

void play(std::string framesDir, int &idx, int &waitTime, Quantities options,
          bool &toggle, std::vector<size_t> frameNum) {
  char buffer[50];
  sprintf(buffer, "/frame%d.ply", (int)idx);
  std::string plyName(buffer);
  plyName = framesDir + plyName;
  registerSurfaceMesh(plyName, options);

  idx++;
  if (idx > frameNum[1]) {
    idx = frameNum[0];
    toggle = !toggle;
  }
  wait(waitTime);
}
