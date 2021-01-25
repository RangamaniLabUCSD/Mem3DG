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

#include "mem3dg/solver/system.h"
#include <iomanip>

namespace mem3dg {

void System::checkParameters() {
  // check validity of parameters / options
  if (mesh.hasBoundary()) {
    if (isReducedVolume || P.Vt != -1.0 || P.cam != 0.0) {
      throw std::logic_error(
          "For open boundary simulation, isReducedVolume has to be false, Vt "
          "has to be -1, and cam has to be 0 ");
    }
  }

  if (!isLocalCurvature) {
    if (P.eta != 0) {
      throw std::logic_error(
          "line tension eta has to be 0 for nonlocal curvature!");
    }
    if (P.r_H0 != std::vector<double>({-1, -1})) {
      throw std::logic_error("r_H0 has to be {-1, -1} for nonlocal curvature!");
    }
    if (P.sharpness != 0) {
      throw std::logic_error(
          "sharpness of transition has to be 0 for nonlocal curvature!");
    }
  }

  if (isReducedVolume) {
    if (P.cam != -1) {
      throw std::logic_error("ambient concentration cam has to be -1 for "
                             "reduced volume parametrized simulation!");
    }
  } else {
    if (P.Vt != -1) {
      throw std::logic_error("reduced volume Vt has to be -1 for "
                             "ambient pressure parametrized simulation!");
    }
  }

  if (!isProtein) {
    if (P.epsilon != -1 || P.Bc != -1) {
      throw std::logic_error("Binding constant Bc and binding energy "
                             "epsilon has to be both -1 for "
                             "protein binding disabled simulation!");
    }
  } else {
    if (isLocalCurvature) {
      throw std::logic_error("Local curvature should be deactivated with "
                             "protein binding activated!");
    }
  }

  if (P.Kf == 0) {
    if (P.conc != -1 || P.height != 0) {
      throw std::logic_error(
          "With no external force, its concentration should be disabled (=-1) "
          "and prescribed height should be set to 0!");
    }
  }
}

void System::pcg_test() {
  // Generate a normal distribution around that mean
  std::normal_distribution<> normal_dist(0, 2);

  // Make a copy of the RNG state to use later
  pcg32 rng_checkpoint = rng;

  // Produce histogram
  std::map<int, int> hist;
  for (int n = 0; n < 10000; ++n) {
    ++hist[std::round(normal_dist(rng))];
  }
  std::cout << "Normal distribution around " << 0 << ":\n";
  for (auto p : hist) {
    std::cout << std::fixed << std::setprecision(1) << std::setw(2) << p.first
              << ' ' << std::string(p.second / 30, '*') << '\n';
  }

  // Produce information about RNG usage
  std::cout << "Required " << (rng - rng_checkpoint) << " random numbers.\n";
}

void System::initConstants() {
  // Initialize RNG
  pcg_extras::seed_seq_from<std::random_device> seed_source;
  rng = pcg32(seed_source);

  /// compute constant values during simulation
  // Find the closest point index to P.pt in refVpg
  closestPtIndToPt(mesh, refVpg, P.pt, ptInd);

  // Initialize the constant mask based on distance from the point specified
  mask =
      (heatMethodDistance(refVpg, mesh.vertex(ptInd)).raw().array() < P.radius)
          .matrix();
  // Mask boundary element
  if (mesh.hasBoundary()) {
    boundaryMask(mesh, mask);
  }

  // Initialize the constant target face/surface areas
  targetFaceAreas = refVpg.faceAreas;
  targetSurfaceArea = targetFaceAreas.raw().sum();

  // Initialize the constant target edge length
  targetEdgeLengths = refVpg.edgeLengths.reinterpretTo(mesh);

  // Initialize the target constant cross length ration
  getCrossLengthRatio(mesh, refVpg, targetLcr);

  // Initialize the constant reference volume
  if (mesh.hasBoundary()) {
    refVolume = 0.0;
  } else {
    refVolume = std::pow(targetSurfaceArea / constants::PI / 4, 1.5) *
                (4 * constants::PI / 3);
  }

  // Initialize the constant spontaneous curvature
  H0.setConstant(mesh.nVertices(), 1, P.H0);
}

void System::updateVertexPositions() {
  vpg.refreshQuantities();

  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(vpg.vertexNormals);
  auto positions = gc::EigenMap<double, 3>(vpg.inputVertexPositions);

  // initialize/update Laplacian matrix
  M_inv = (1 / (M.diagonal().array())).matrix().asDiagonal();

  // initialize/update distance from the point specified
  // if (isLocalCurvature) {
  //   geodesicDistanceFromPtInd = heatSolver.computeDistance(mesh.vertex(ptInd));
  // }

  // initialize/update spontaneous curvature (protein binding)
  if (isProtein) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> proteinDensitySq =
        (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();
    H0 = (P.H0 * proteinDensitySq.array() / (1 + proteinDensitySq.array()))
             .matrix();
  }

  // initialize/update spontaneous curvature (local spontaneous curvature)
  if (isLocalCurvature) {
    tanhDistribution(vpg, H0, geodesicDistanceFromPtInd.raw(), P.sharpness,
                     P.r_H0);
    H0 *= P.H0;
  }

  // initialize/update mean curvature
  H = rowwiseDotProduct(M_inv * L * positions / 2.0, vertexAngleNormal_e);

  /// initialize/update enclosed volume
  volume = 0;
  for (gcs::Face f : mesh.faces()) {
    volume += signedVolumeFromFace(
        f, vpg, refVpg.inputVertexPositions[mesh.vertex(ptInd)]);
  }

  // initialize/update total surface area
  surfaceArea = vpg.faceAreas.raw().sum();

  // initialize/update intersection area
  interArea = 0.0;
  for (gcs::Vertex v : mesh.vertices()) {
    if ((H0[v.getIndex()] > (0.1 * P.H0)) &&
        (H0[v.getIndex()] < (0.9 * P.H0)) && (H[v.getIndex()] != 0)) {
      interArea += vpg.vertexDualAreas[v];
    }
  }

  // initialize/update external force
  getExternalPressure();

  // initialize/update the vertex position of the last iteration
  pastPositions = vpg.inputVertexPositions;
}
} // namespace mem3dg