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

#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/solver/mutable_trajfile.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"

#include <geometrycentral/surface/meshio.h>

#include "mem3dg/mem3dg"
#include "visualization.h"
//#include <pybind11/embed.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

// ==========================================================
// =============        Viewers                ==============
// ==========================================================
void visualize(mem3dg::solver::System &f) {
  signal(SIGINT, mem3dg::signalHandler);
  // Initialize visualization variables
  float transparency = 1;

  // Set preference for polyscope
  initGui();

  /// Initiate in polyscope
  polyscope::init();

  polyscope::SurfaceMesh *polyMesh = polyscope::registerSurfaceMesh(
      "Membrane", f.vpg->inputVertexPositions, f.mesh->getFaceVertexList(),
      gcs::polyscopePermutations(*f.mesh));
  polyMesh->setSmoothShade(true);
  polyMesh->setEnabled(true);
  polyMesh->setEdgeWidth(1);

  // Process attributes
  Eigen::Matrix<double, Eigen::Dynamic, 1> fn;

  fn = f.forces.bendingForce.raw() + f.forces.capillaryForce.raw() +
       f.forces.osmoticForce.raw() + f.forces.externalForce.raw() +
       f.forces.lineCapillaryForce.raw();

  /// Read element data
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("mean_curvature",
                                f.vpg->vertexMeanCurvatures.raw().array() /
                                    f.vpg->vertexDualAreas.raw().array());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("gauss_curvature",
                                f.vpg->vertexGaussianCurvatures.raw().array() /
                                    f.vpg->vertexDualAreas.raw().array());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("spon_curvature", f.H0);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("external_Force", f.forces.externalForce);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("line_tension_pressure",
                                f.forces.lineCapillaryForce.raw().array() /
                                    f.vpg->vertexDualAreas.raw().array());
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("physical_force", fn);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity("bending_rigidity", f.Kb);
  polyscope::getSurfaceMesh("Membrane")
      ->addVertexScalarQuantity(
          "-lapH(smoothing)",
          -f.Kb.raw().array() * (f.vpg->vertexLumpedMassMatrix.cwiseInverse() *
                                 f.vpg->cotanLaplacian *
                                 (f.vpg->vertexLumpedMassMatrix.cwiseInverse() *
                                      f.vpg->vertexMeanCurvatures.raw() -
                                  f.H0.raw()))
                                    .array());
  polyscope::getSurfaceMesh("Membrane")
      ->addEdgeScalarQuantity("edge_dihedral", f.vpg->edgeDihedralAngles);
  polyscope::getSurfaceMesh("Membrane")
      ->addEdgeScalarQuantity("cotan weight", f.vpg->edgeCotanWeights);
  polyscope::getSurfaceMesh("Membrane")
      ->addEdgeScalarQuantity("edge_length", f.vpg->edgeLengths);
  polyscope::getSurfaceMesh("Membrane")
      ->addFaceCountQuantity(
          "the point", std::vector<std::pair<std::size_t, int>>{std::make_pair(
                           f.thePoint.inSomeFace().face.getIndex(), 1)})
      ->setPointRadius(sqrt(f.surfaceArea / f.mesh->nFaces() * 4 / sqrt(3)) / 2,
                       false);

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
                std::vector<std::size_t> frameNum, double mapMinLim,
                double mapMaxLim, float transparency, float fov,
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
  auto polyscopeMesh =
      registerSurfaceMesh(plyName, options, mapMinLim, mapMaxLim);
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
      polyscopeMesh =
          registerSurfaceMesh(plyName, options, mapMinLim, mapMaxLim);
      prevFrame = currFrame;
    }
    if (isRecord) {
      char buff[50];
      snprintf(buff, 50, "video/frame%06d.png", currFrame);
      std::string defaultName(buff);
      polyscope::screenshot(defaultName, true);
      play(polyscopeMesh, frameDir, currFrame, waitTime, options, isRecord,
           frameNum, mapMinLim, mapMaxLim);
      prevFrame = currFrame;
    }
    if (isStart) {
      play(polyscopeMesh, frameDir, currFrame, waitTime, options, isStart,
           frameNum, mapMinLim, mapMaxLim);
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
  // mem3dg::solver::TrajFile fd =
  //     mem3dg::solver::TrajFile::openReadOnly(filename);
  mem3dg::solver::MutableTrajFile fd =
      mem3dg::solver::MutableTrajFile::openReadOnly(filename);

  // Initialize visualization variables
  int prevFrame = 0;
  int currFrame = 0;
  bool isStart = false;
  bool isRecord = false;
  int maxFrame = fd.nFrames() - 1;
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
  // mem3dg::solver::TrajFile fd =
  // mem3dg::solver::TrajFile::openReadOnly(filename);
  mem3dg::solver::MutableTrajFile fd =
      mem3dg::solver::MutableTrajFile::openReadOnly(filename);
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

// ==========================================================
// =============        Helper functions       ==============
// ==========================================================

#ifdef MEM3DG_WITH_NETCDF
polyscope::SurfaceMesh *registerSurfaceMesh(mem3dg::solver::TrajFile &fd,
                                            int idx,
                                            const Quantities &options) {
  if (idx >= fd.nFrames()) {
    idx = 0;
  }

  mem3dg::EigenVectorX3ur topo_frame = fd.getTopoFrame(idx);
  mem3dg::EigenVectorX3dr coords = fd.getCoords(idx);
  // mesh->updateVertexPositions(coords);
  polyscope::SurfaceMesh *polyscopeMesh =
      polyscope::registerSurfaceMesh("Mesh", coords, topo_frame);
  // polyscopeMesh->setEnabled(true);

  if (options.ref_coord) {
    mem3dg::EigenVectorX3dr refcoords = fd.getRefcoordinate();
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
  if (options.bending_force) {
    polyscopeMesh->addVertexScalarQuantity("bending_force",
                                           fd.getBendingForce(idx));
  }
  if (options.capillary_force) {
    polyscopeMesh->addVertexScalarQuantity("capillary_force",
                                           fd.getCapillaryForce(idx));
  }
  if (options.osmotic_force) {
    polyscopeMesh->addVertexScalarQuantity("inside_force",
                                           fd.getOsmoticForce(idx));
  }
  if (options.ext_force) {
    polyscopeMesh->addVertexScalarQuantity("external_force",
                                           fd.getExternalForce(idx));
  }
  if (options.line_force) {
    polyscopeMesh->addVertexScalarQuantity("line_tension_force",
                                           fd.getLineForce(idx));
  }
  if (options.physical_force) {
    polyscopeMesh->addVertexScalarQuantity("physical_force",
                                           fd.getPhysicalForce(idx));
  }
  // polyscope::registerSurfaceMesh("Mesh", coords, top);

  return polyscopeMesh;
}

polyscope::SurfaceMesh *registerSurfaceMesh(mem3dg::solver::MutableTrajFile &fd,
                                            int idx,
                                            const Quantities &options) {
  if (idx >= fd.nFrames()) {
    idx = 0;
  }

  mem3dg::EigenVectorX3ur topo_frame = fd.getTopology(idx);
  mem3dg::EigenVectorX3dr coords = fd.getCoords(idx);
  // mesh->updateVertexPositions(coords);
  polyscope::SurfaceMesh *polyscopeMesh =
      polyscope::registerSurfaceMesh("Mesh", coords, topo_frame);
  // polyscopeMesh->setEnabled(true);

  return polyscopeMesh;
}

void play(polyscope::SurfaceMesh *&polyscopeMesh, mem3dg::solver::TrajFile &fd,
          int &idx, int &waitTime, Quantities options, bool &toggle) {

  polyscopeMesh = registerSurfaceMesh(fd, idx, options);
  idx++;
  if (idx >= fd.nFrames()) {
    idx = 0;
    toggle = !toggle;
  }
  wait(waitTime);
}

void play(polyscope::SurfaceMesh *&polyscopeMesh,
          mem3dg::solver::MutableTrajFile &fd, int &idx, int &waitTime,
          Quantities options, bool &toggle) {

  polyscopeMesh = registerSurfaceMesh(fd, idx, options);
  idx++;
  if (idx >= fd.nFrames()) {
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
                                            const Quantities &options,
                                            double mapMinLim,
                                            double mapMaxLim) {

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
      ptrMesh->getFaceVertexList(), gcs::polyscopePermutations(*ptrMesh));
  polyscopeMesh->setSmoothShade(true);
  polyscopeMesh->setEnabled(true);
  // polyscopeMesh->setEdgeWidth(1);

  /// Read element data
  if (mapMinLim == 0 && mapMaxLim == 0) {
    if (options.mean_curvature) {
      polyscopeMesh->addVertexScalarQuantity(
          "mean_curvature",
          ptrRichData->getVertexProperty<double>("mean_curvature"));
      // polyscopeMesh
      //     ->addVertexScalarQuantity(
      //         "mean_curvature",
      //         ptrRichData->getVertexProperty<double>("mean_curvature"))
      //     ->setMapRange(std::make_pair(-2.54, 10.29));
    }
    if (options.velocity) {
      polyscopeMesh->addVertexScalarQuantity(
          "velocity", ptrRichData->getVertexProperty<double>("velocity"));
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

    if (options.ext_force) {
      polyscopeMesh->addVertexScalarQuantity(
          "external_force",
          ptrRichData->getVertexProperty<double>("external_force"));
    }
    if (options.avoidance_force) {
      polyscopeMesh->addVertexScalarQuantity(
          "avoidance_force",
          ptrRichData->getVertexProperty<double>("avoidance_force"));
    }
    if (options.physical_force) {
      polyscopeMesh->addVertexScalarQuantity(
          "physical_force",
          ptrRichData->getVertexProperty<double>("physical_force"));
    }
    if (options.capillary_force) {
      polyscopeMesh->addVertexScalarQuantity(
          "capillary_force",
          ptrRichData->getVertexProperty<double>("capillary_force"));
    }
    if (options.bending_force) {
      polyscopeMesh->addVertexScalarQuantity(
          "bending_force",
          ptrRichData->getVertexProperty<double>("bending_force"));
    }
    if (options.deviatoric_force) {
      polyscopeMesh->addVertexScalarQuantity(
          "deviatoric_force",
          ptrRichData->getVertexProperty<double>("deviatoric_force"));
    }
    if (options.line_force) {
      polyscopeMesh->addVertexScalarQuantity(
          "line_tension_force",
          ptrRichData->getVertexProperty<double>("line_tension_force"));
    }
    if (options.osmotic_force) {
      polyscopeMesh->addVertexScalarQuantity(
          "osmotic_force",
          ptrRichData->getVertexProperty<double>("osmotic_force"));
    }
    if (options.adsorption_force) {
      polyscopeMesh->addVertexScalarQuantity(
          "adsorption_force",
          ptrRichData->getVertexProperty<double>("adsorption_force"));
    }
    if (options.aggregation_force) {
      polyscopeMesh->addVertexScalarQuantity(
          "aggregation_force",
          ptrRichData->getVertexProperty<double>("aggregation_force"));
    }
    if (options.mask) {
      polyscopeMesh->addVertexScalarQuantity(
          "force_mask", ptrRichData->getVertexProperty<double>("force_mask"));
    }
    if (options.the_point) {
      polyscopeMesh
          ->addVertexCountQuantity(
              "the_point",
              getCountQuantities(
                  ptrRichData->getVertexProperty<int>("the_point")))
          ->setPointRadius(0.01, true);
    }
    if (options.smoothing_mask) {
      polyscopeMesh->addVertexScalarQuantity(
          "smoothing_mask",
          ptrRichData->getVertexProperty<int>("smoothing_mask"));
    }
    if (options.chemical_potential) {
      polyscopeMesh->addVertexScalarQuantity(
          "chemical_potentialx1000",
          ptrRichData->getVertexProperty<double>("chemical_potential") * 1000);
    }
    if (options.bending_potential) {
      polyscopeMesh->addVertexScalarQuantity(
          "bending_potentialx1000",
          ptrRichData->getVertexProperty<double>("bending_potential") * 1000);
    }
    if (options.deviatoric_potential) {
      polyscopeMesh->addVertexScalarQuantity(
          "deviatoric_potentialx1000",
          ptrRichData->getVertexProperty<double>("deviatoric_potential") *
              1000);
    }
    if (options.diffusion_potential) {
      polyscopeMesh->addVertexScalarQuantity(
          "diffusion_potentialx1000",
          ptrRichData->getVertexProperty<double>("diffusion_potential") * 1000);
    }
    if (options.adsorption_potential) {
      polyscopeMesh->addVertexScalarQuantity(
          "adsorption_potentialx1000",
          ptrRichData->getVertexProperty<double>("adsorption_potential") *
              1000);
    }
    if (options.aggregation_potential) {
      polyscopeMesh->addVertexScalarQuantity(
          "aggregation_potentialx1000",
          ptrRichData->getVertexProperty<double>("aggregation_potential") *
              1000);
    }
    /*gcs::VertexData<gc::Vector3> vertexVelocity =
        ptrRichData->getVertexProperty<gc::Vector3>("vertex_velocity");*/
    /*gcs::VertexData<gc::Vector3> normalForce =
    ptrRichData->getVertexProperty<gc::Vector3>("normal_force");
    gcs::VertexData<gc::Vector3> tangentialForce =
    ptrRichData->getVertexProperty<gc::Vector3>("tangential_force");*/
    // mem3dg::EigenVectorX3dr vertexVelocity_e =
    //    mem3dg::EigenMap<double, 3>(vertexVelocity);
    /*mem3dg::EigenVectorX3dr normalForce_e =
    gc::EigenMap<double, 3>(normalForce); Eigen::Matrix<double,
    Eigen::Dynamic, 3> tangentialForce_e = gc::EigenMap<double,
    3>(tangentialForce);*/
    /*polyscope::getSurfaceMesh("Vesicle surface")
        ->addVertexVectorQuantity("vertexVelocity", vertexVelocity_e);*/
    /*polyscope::getSurfaceMesh("Vesicle
    surface")->addVertexVectorQuantity("tangential_force", tangentialForce_e);
    polyscope::getSurfaceMesh("Vesicle
    surface")->addVertexVectorQuantity("normal_force", normalForce_e);*/
  } else {
    if (options.mean_curvature) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "mean_curvature",
              ptrRichData->getVertexProperty<double>("mean_curvature"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.velocity) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "velocity", ptrRichData->getVertexProperty<double>("velocity"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.gauss_curvature) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "gauss_curvature",
              ptrRichData->getVertexProperty<double>("gauss_curvature"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.spon_curvature) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "spon_curvature",
              ptrRichData->getVertexProperty<double>("spon_curvature"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }

    if (options.ext_force) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "external_force",
              ptrRichData->getVertexProperty<double>("external_force"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.physical_force) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "physical_force",
              ptrRichData->getVertexProperty<double>("physical_force"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.capillary_force) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "capillary_force",
              ptrRichData->getVertexProperty<double>("capillary_force"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.bending_force) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "bending_force",
              ptrRichData->getVertexProperty<double>("bending_force"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.deviatoric_force) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "deviatoric_force",
              ptrRichData->getVertexProperty<double>("deviatoric_force"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.line_force) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "line_tension_force",
              ptrRichData->getVertexProperty<double>("line_tension_force"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.osmotic_force) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "osmotic_force",
              ptrRichData->getVertexProperty<double>("osmotic_force"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.adsorption_force) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "adsorption_force",
              ptrRichData->getVertexProperty<double>("adsorption_force"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.aggregation_force) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "aggregation_force",
              ptrRichData->getVertexProperty<double>("aggregation_force"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.avoidance_force) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "avoidance_force",
              ptrRichData->getVertexProperty<double>("avoidance_force"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.mask) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "force_mask",
              ptrRichData->getVertexProperty<double>("force_mask"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.smoothing_mask) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "smoothing_mask",
              ptrRichData->getVertexProperty<int>("smoothing_mask"))
          ->setMapRange(std::make_pair(mapMinLim, mapMaxLim));
    }
    if (options.chemical_potential) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "chemical_potentialx1000",
              ptrRichData->getVertexProperty<double>("chemical_potential") *
                  1000)
          ->setMapRange(std::make_pair(mapMinLim * 1000, mapMaxLim * 1000));
    }
    if (options.bending_potential) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "bending_potentialx1000",
              ptrRichData->getVertexProperty<double>("bending_potential") *
                  1000)
          ->setMapRange(std::make_pair(mapMinLim * 1000, mapMaxLim * 1000));
    }
    if (options.deviatoric_potential) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "deviatoric_potentialx1000",
              ptrRichData->getVertexProperty<double>("deviatoric_potential") *
                  1000)
          ->setMapRange(std::make_pair(mapMinLim * 1000, mapMaxLim * 1000));
    }
    if (options.diffusion_potential) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "diffusion_potentialx1000",
              ptrRichData->getVertexProperty<double>("diffusion_potential") *
                  1000)
          ->setMapRange(std::make_pair(mapMinLim * 1000, mapMaxLim * 1000));
    }
    if (options.adsorption_potential) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "adsorption_potentialx1000",
              ptrRichData->getVertexProperty<double>("adsorption_potential") *
                  1000)
          ->setMapRange(std::make_pair(mapMinLim * 1000, mapMaxLim * 1000));
    }
    if (options.aggregation_potential) {
      polyscopeMesh
          ->addVertexScalarQuantity(
              "aggregation_potentialx1000",
              ptrRichData->getVertexProperty<double>("aggregation_potential") *
                  1000)
          ->setMapRange(std::make_pair(mapMinLim * 1000, mapMaxLim * 1000));
    }
    if (options.the_point) {
      polyscopeMesh
          ->addVertexCountQuantity(
              "the_point",
              getCountQuantities(
                  ptrRichData->getVertexProperty<int>("the_point")))
          ->setPointRadius(0.01, true);
    }
  }
  return polyscopeMesh;
}

void play(polyscope::SurfaceMesh *&polyscopeMesh, std::string framesDir,
          int &idx, int &waitTime, Quantities options, bool &toggle,
          std::vector<std::size_t> frameNum, double mapMinLim,
          double mapMaxLim) {
  char buffer[50];
  sprintf(buffer, "/frame%d.ply", (int)idx);
  std::string plyName(buffer);
  plyName = framesDir + plyName;
  polyscopeMesh = registerSurfaceMesh(plyName, options, mapMinLim, mapMaxLim);
  idx++;
  if (idx > frameNum[1]) {
    idx = frameNum[0];
    toggle = !toggle;
  }
  wait(waitTime);
}

std::vector<std::pair<std::size_t, int>>
getCountQuantities(gc::VertexData<int> &&meshData) {
  std::vector<std::pair<std::size_t, int>> values;
  for (gc::Vertex v : meshData.getMesh()->vertices()) {
    if (meshData[v] != 0) {
      values.push_back(std::make_pair(v.getIndex(), meshData[v]));
    }
  }
  return values;
}
