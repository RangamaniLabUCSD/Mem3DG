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
#include <memory>
#include <time.h>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"

#include <geometrycentral/surface/meshio.h>

#include "mem3dg/mem3dg"
//#include <pybind11/embed.h>

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
      ->addVertexScalarQuantity("mean_curvature", f.vpg->vertexMeanCurvatures);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("gauss_curvature", f.K);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("spon_curvature", f.H0);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("external_pressure", f.externalPressure);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("bending_pressure", f.bendingPressure);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("line_tension_pressure",
                                f.M_inv * f.lineCapillaryForce.raw());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("capillary_pressure", f.capillaryPressure);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("inside_pressure", f.insidePressure);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("physical_pressure", fn);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("bending_rigidity", f.Kb);
  polyscope::getSurfaceMesh("Membrane")
      ->addEdgeScalarQuantity("line_tension", f.lineTension);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity(
          "-lapH(smoothing)",
          -f.Kb.raw().array() *
              (f.M_inv * f.L * (f.H.raw() - f.H0.raw())).array());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity(
          "spon part)",
          f.bendingPressure.raw().array() +
              f.Kb.raw().array() *
                  (f.M_inv * f.L * (f.H.raw() - f.H0.raw())).array());
  polyscope::getSurfaceMesh("Membrane")
      ->addEdgeScalarQuantity("edge_dihedral", f.vpg->edgeDihedralAngles);
  std::cout << "no of edges: " << f.mesh->nEdges() << std::endl;
  std::cout << "no of edges: " << f.vpg->edgeLengths.raw().rows() << std::endl;
  polyscope::getSurfaceMesh("Membrane")
      ->addEdgeScalarQuantity("edge_length", f.vpg->edgeLengths);
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
  std::vector<std::pair<size_t, int>> values;
  values.push_back(
      std::make_pair(f.thePoint.nearestVertex().getIndex(), 1));
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexCountQuantity("the point", values);
  std::vector<std::pair<size_t, int>> values2;
  gcs::Vertex theVertex;
  double shorestDistance = 1e18;
  double distance;
  for (gcs::Vertex v : f.mesh->vertices()) {
    distance = (gc::Vector2{f.vpg->inputVertexPositions[v].x,
                            f.vpg->inputVertexPositions[v].y})
                   .norm();
    if (distance < shorestDistance) {
      shorestDistance = distance;
      theVertex = v;
    }
  }
  std::cout << "the vertex is: " << theVertex.getIndex()
            << "and the distance is: " << shorestDistance << std::endl;
  values2.push_back(std::make_pair(theVertex.getIndex(), 2));
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexCountQuantity("the point2", values2);

  gc::VertexData<double> xydistance(*f.mesh);
  for (gcs::Vertex v : f.mesh->vertices()) {
    xydistance[v] = std::sqrt(
        f.vpg->inputVertexPositions[v].x * f.vpg->inputVertexPositions[v].x +
        f.vpg->inputVertexPositions[v].y * f.vpg->inputVertexPositions[v].y);
  }
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("x y distance", xydistance);

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
      polyscope::screenshot("screenshot.png", true);
    }

    ImGui::PopItemWidth();
  };
  polyscope::state::userCallback = myCallback;
  polyscope::show();
}

int snapshot_ply(std::string fileName, const Quantities &option,
                 float transparency, float fov, float edgeWidth) {

  signal(SIGINT, mem3dg::signalHandler);

  // Set preference for polyscope
  initGui();
  polyscope::view::fov = fov;

  /// Initiate in polyscope
  polyscope::init();

  // Initialize surface mesh
  auto polyscopeMesh = registerSurfaceMesh(fileName, option);
  polyscopeMesh->setSmoothShade(true);
  polyscopeMesh->setEdgeWidth(edgeWidth);
  polyscopeMesh->setTransparency(transparency);

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

int animate_ply(std::string frameDir, const Quantities &options,
                std::vector<size_t> frameNum, float transparency, float fov,
                float edgeWidth) {

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
  polyscope::view::fov = fov;

  // Initialize polyscope
  polyscope::init();

  // Initialize surface mesh
  char buffer[50];
  sprintf(buffer, "/frame%d.ply", (int)frameNum[0]);
  std::string plyName(buffer);
  plyName = frameDir + plyName;
  auto polyscopeMesh = registerSurfaceMesh(plyName, options);
  polyscopeMesh->setSmoothShade(true);
  polyscopeMesh->setEdgeWidth(edgeWidth);
  polyscopeMesh->setTransparency(transparency);

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
    ImGui::SliderFloat("transparency", &transparency, 0, 1);
    ImGui::SliderInt("slow-mo", &waitTime, 0,
                     maxWaitTime); // set a float variable

    // Execute transparency slider
    polyscopeMesh->setTransparency(transparency);

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
      play(polyscopeMesh, frameDir, currFrame, waitTime, options, isRecord,
           frameNum);
      prevFrame = currFrame;
    }
    if (isStart) {
      play(polyscopeMesh, frameDir, currFrame, waitTime, options, isStart,
           frameNum);
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
               float transparency, float fov, float edgeWidth) {

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
  polyscopeMesh->setEdgeWidth(edgeWidth);
  polyscopeMesh->setTransparency(transparency);

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

    // Execute transparency slider
    polyscopeMesh->setTransparency(transparency);

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
      play(polyscopeMesh, fd, currFrame, waitTime, options, isRecord);
      prevFrame = currFrame;
    }
    if (isStart) {
      play(polyscopeMesh, fd, currFrame, waitTime, options, isStart);
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
    polyscopeMesh->addVertexScalarQuantity("mean_curvature",
                                           fd.getMeanCurvature(idx));
  }
  if (options.gauss_curvature) {
    polyscopeMesh->addVertexScalarQuantity("gauss_curvature",
                                           fd.getGaussCurvature(idx));
  }
  if (options.spon_curvature) {
    polyscopeMesh->addVertexScalarQuantity("spon_curvature",
                                           fd.getSponCurvature(idx));
  }
  if (options.mask) {
    polyscopeMesh->addVertexScalarQuantity("mask", fd.getMask());
  }
  if (options.H_H0) {
    polyscopeMesh->addVertexScalarQuantity("curvature_diff",
                                           fd.getH_H0_diff(idx));
  }
  if (options.velocity) {
    polyscopeMesh->addVertexVectorQuantity("velocity", fd.getVelocity(idx));
  }
  if (options.bending_pressure) {
    polyscopeMesh->addVertexScalarQuantity("bending_pressure",
                                           fd.getBendingPressure(idx));
  }
  if (options.capillary_pressure) {
    polyscopeMesh->addVertexScalarQuantity("capillary_pressure",
                                           fd.getCapillaryPressure(idx));
  }
  if (options.inside_pressure) {
    polyscopeMesh->addVertexScalarQuantity("inside_pressure",
                                           fd.getInsidePressure(idx));
  }
  if (options.ext_pressure) {
    polyscopeMesh->addVertexScalarQuantity("external_pressure",
                                           fd.getExternalPressure(idx));
  }
  if (options.line_pressure) {
    polyscopeMesh->addVertexScalarQuantity("line_tension_pressure",
                                           fd.getLinePressure(idx));
  }
  if (options.physical_pressure) {
    polyscopeMesh->addVertexScalarQuantity("physical_pressure",
                                           fd.getPhysicalPressure(idx));
  }
  // polyscope::registerSurfaceMesh("Mesh", coords, top);

  return polyscopeMesh;
}

void play(polyscope::SurfaceMesh *&polyscopeMesh, mem3dg::TrajFile &fd,
          int &idx, int &waitTime, Quantities options, bool &toggle) {

  polyscopeMesh = registerSurfaceMesh(fd, idx, options);
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

  // Declare pointers to mesh, geometry and richdata objects
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;

  // Load input mesh
  // Case of viewing raw .ply file
  // if (!options.mean_curvature && !options.gauss_curvature &&
  //     !options.spon_curvature && !options.ext_pressure &&
  //     !options.physical_pressure && !options.capillary_pressure &&
  //     !options.bending_pressure && !options.line_pressure) {
  //   std::tie(ptrMesh, ptrVpg) = gcs::readManifoldSurfaceMesh(plyName);
  // } else {
  std::tie(ptrMesh, ptrRichData) =
      gcs::RichSurfaceMeshData::readMeshAndData(plyName);
  ptrVpg = ptrRichData->getGeometry();
  // }

  polyscope::SurfaceMesh *polyscopeMesh = polyscope::registerSurfaceMesh(
      "Vesicle surface", ptrVpg->inputVertexPositions,
      ptrMesh->getFaceVertexList());
  polyscopeMesh->setSmoothShade(true);
  polyscopeMesh->setEnabled(true);
  polyscopeMesh->setEdgeWidth(1);

  /// Read element data
  if (options.mean_curvature) {
    polyscopeMesh->addVertexScalarQuantity(
        "mean_curvature",
        ptrRichData->getVertexProperty<double>("mean_curvature"));
  }
  if (options.gauss_curvature) {
    polyscopeMesh->addVertexScalarQuantity(
        "gauss_curvature",
        ptrRichData->getVertexProperty<double>("gauss_curvature"));
  }
  if (options.spon_curvature) {
    polyscopeMesh->addVertexScalarQuantity(
        "spon_curvature",
        ptrRichData->getVertexProperty<double>("spon_curvature"));
  }

  if (options.ext_pressure) {
    polyscopeMesh->addVertexScalarQuantity(
        "external_pressure",
        ptrRichData->getVertexProperty<double>("external_pressure"));
  }
  if (options.physical_pressure) {
    polyscopeMesh->addVertexScalarQuantity(
        "physical_pressure",
        ptrRichData->getVertexProperty<double>("physical_pressure"));
  }
  if (options.capillary_pressure) {
    polyscopeMesh->addVertexScalarQuantity(
        "capillary_pressure",
        ptrRichData->getVertexProperty<double>("capillary_pressure"));
  }
  if (options.bending_pressure) {
    polyscopeMesh->addVertexScalarQuantity(
        "bending_pressure",
        ptrRichData->getVertexProperty<double>("bending_pressure"));
  }
  if (options.line_pressure) {
    polyscopeMesh->addVertexScalarQuantity(
        "line_tension_pressure",
        ptrRichData->getVertexProperty<double>("line_tension_pressure"));
  }
  if (options.inside_pressure) {
    polyscopeMesh->addVertexScalarQuantity(
        "inside_pressure",
        ptrRichData->getVertexProperty<double>("inside_pressure"));
  }

  if (options.mask) {
    polyscopeMesh->addVertexScalarQuantity(
        "mask", ptrRichData->getVertexProperty<int>("mask"));
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

void play(polyscope::SurfaceMesh *&polyscopeMesh, std::string framesDir,
          int &idx, int &waitTime, Quantities options, bool &toggle,
          std::vector<size_t> frameNum) {
  char buffer[50];
  sprintf(buffer, "/frame%d.ply", (int)idx);
  std::string plyName(buffer);
  plyName = framesDir + plyName;
  polyscopeMesh = registerSurfaceMesh(plyName, options);
  idx++;
  if (idx > frameNum[1]) {
    idx = frameNum[0];
    toggle = !toggle;
  }
  wait(waitTime);
}
