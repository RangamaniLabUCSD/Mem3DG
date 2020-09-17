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

#include <cassert>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/heat_method_distance.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <pcg_random.hpp>
#include <random>

#include <math.h>

#include "mem3dg/solver/macros.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/util.h"

#include <vector>

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

struct Parameters {
  /// Bending modulus
  double Kb;
  /// Spontaneous curvature
  double H0;
  /// Sharpness of the spontaneous curvature hetergeneity 
  double sharpness;
  /// radius of non-zero spontaneous curvature
  double r_H0;
  /// Global stretching modulus
  double Ksg;
  /// Vertex shifting constant
  double Kst;
  /// Local stretching modulus
  double Ksl;
  /// Edge spring constant
  double Kse;
  /// Volume regularization
  double Kv;
  /// binding energy per protein
  double epsilon;
  /// binding constant
  double Bc;
  /// Dissipation coefficient
  double gamma;
  /// Reduced volume
  double Vt;
  /// Boltzmann constant*Temperature
  double kt;
  /// Noise
  double sigma;
  /// index of node with applied external force
  size_t ptInd;
  /// Magnitude of external force
  double Kf;
  /// level of concentration of the external force
  double conc;
  /// target height
  double height;
  /// domain of integration
  double radius;
};

class DLL_PUBLIC Force {
public:
  /// Parameters
  Parameters P;
  /// Cached mesh of interest
  gcs::SurfaceMesh &mesh;
  /// Cached mesh data
  gcs::RichSurfaceMeshData &richData;
  /// Embedding and other geometric details
  gcs::VertexPositionGeometry &vpg;
  /// reference embedding geometry
  gcs::VertexPositionGeometry &refVpg;

  /// Cached bending stress
  gcs::VertexData<gc::Vector3> bendingPressure;
  /// Cached tension-induced capillary pressure
  gcs::VertexData<gc::Vector3> capillaryPressure;
  /// Cached relative inside pressure
  gcs::VertexData<gc::Vector3> insidePressure;
  /// Cached externally-applied pressure
  gcs::VertexData<gc::Vector3> externalPressure;

  /// Cached local stretching forces (in-plane regularization)
  gcs::VertexData<gc::Vector3> regularizationForce;
  /// Cached damping forces
  gcs::VertexData<gc::Vector3> dampingForce;
  /// Cached stochastic forces
  gcs::VertexData<gc::Vector3> stochasticForce;

  /// Cached protein surface density 
  gcs::VertexData<double> proteinDensity;
  /// Cached chemical potential
  gcs::VertexData<double> chemicalPotential;

  /// Whether or not use tufted laplacian matrix
  bool isTuftedLaplacian;
  /// Mollify Factor in constructing tufted laplacian matrix
  double mollifyFactor;
  /// Whether or not do vertex shift
  bool isVertexShift;
  /// Whether or not consider protein binding
  bool isProtein;

  /// Cached galerkin mass matrix
  Eigen::SparseMatrix<double> M;
  /// Inverted galerkin mass matrix
  Eigen::SparseMatrix<double> M_inv;
  /// Cotangent Laplacian
  Eigen::SparseMatrix<double> L;
  /// Cached geodesic distance
  gcs::VertexData<double> geodesicDistanceFromAppliedForce;

  /// Target area per face
  gcs::FaceData<double> targetFaceAreas;
  /// Target total face area
  double targetSurfaceArea = 0.0;
  /// surface area
  double surfaceArea = 0.0;
  /// Maximal volume
  double refVolume = 0.0;
  /// Volume
  double volume = 0.0;
  /// Target length per edge
  gcs::EdgeData<double> targetEdgeLengths;
  /// Cached vertex positions from the previous step
  gcs::VertexData<gc::Vector3> pastPositions;
  /// Cached vertex velocity by finite differencing past and current position
  gcs::VertexData<gc::Vector3> vel;
  // Mean curvature of the mesh
  Eigen::Matrix<double, Eigen::Dynamic, 1> H;
  // Spontaneous curvature of the mesh
  Eigen::Matrix<double, Eigen::Dynamic, 1> H0;
  /// Random number engine
  pcg32 rng;
  std::normal_distribution<double> normal_dist;
  /// magnitude of externally-applied pressure
  Eigen::Matrix<double, Eigen::Dynamic, 1> externalPressureMagnitude;
  /// indices for vertices chosen for integration
  Eigen::Matrix<bool, Eigen::Dynamic, 1> mask;

  /*
   * @brief Construct a new Force object
   *
   * @param mesh_         Mesh connectivity
   * @param vpg_          Embedding and geometry information
   * @param time_step_    Numerical timestep
   */

  Force(gcs::SurfaceMesh &mesh_, gcs::VertexPositionGeometry &vpg_,
        gcs::VertexPositionGeometry &refVpg_,
        gcs::RichSurfaceMeshData &richData_, Parameters &p, bool isProtein_ = false,
        bool isTuftedLaplacian_ = false, double mollifyFactor_ = 1e-6,
        bool isVertexShift_ = false)
      : mesh(mesh_), vpg(vpg_), richData(richData_), refVpg(refVpg_),
        P(p), isTuftedLaplacian(isTuftedLaplacian_), isProtein(isProtein_), mollifyFactor(mollifyFactor_),
        isVertexShift(isVertexShift_), bendingPressure(mesh_, {0, 0, 0}),
        insidePressure(mesh_, {0, 0, 0}), capillaryPressure(mesh_, {0, 0, 0}),
        externalPressure(mesh_, {0, 0, 0}), regularizationForce(mesh_, {0,0,0}), 
        stochasticForce(mesh_, {0, 0, 0}), dampingForce(mesh_, {0, 0, 0}),
        proteinDensity(mesh_, 0), vel(mesh_, {0, 0, 0}) {

    // Initialize RNG
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    rng = pcg32(seed_source);

    // GC computed properties
    vpg.requireFaceNormals();
    vpg.requireVertexGalerkinMassMatrix();
    vpg.requireVertexLumpedMassMatrix();
    vpg.requireCotanLaplacian();
    vpg.requireFaceAreas();
    vpg.requireVertexIndices();
    vpg.requireVertexGaussianCurvatures();
    vpg.requireFaceIndices();
    vpg.requireEdgeLengths();
    vpg.requireVertexNormals();
    vpg.requireVertexDualAreas();
    vpg.requireCornerAngles();

    refVpg.requireFaceAreas();
    refVpg.requireEdgeLengths();

    // Initialize the geodesic distance from ptInd
    geodesicDistanceFromAppliedForce =
        heatMethodDistance(vpg, mesh.vertex(P.ptInd));
    auto &dist_e = geodesicDistanceFromAppliedForce.raw();

    // Initialize the external pressure magnitude distribution
    gaussianDistribution(externalPressureMagnitude, dist_e,
                         dist_e.maxCoeff() / P.conc);
    externalPressureMagnitude *= P.Kf;

    // Initialize the spontaneous curvature distribution
    if (isProtein) {
      proteinDensity.raw().setZero();
      H0.setZero(mesh.nVertices(), 1);
    } else {
      tanhDistribution(H0, dist_e, P.sharpness, P.r_H0);
      H0 *= P.H0;
    }

    // Initialize the mask on choosing integration vertices based on geodesic
    // distance from the local external force location on the reference geometry
    mask = (heatMethodDistance(refVpg, mesh.vertex(P.ptInd)).raw().array() <
            P.radius)
               .matrix();

    // Regularize the vetex position geometry if needed
    if (isVertexShift) {
      vertexShift(mesh, vpg, mask);
      update_Vertex_positions();
    }

    // Initialize the mass and conformal Laplacian matrix
    if (isTuftedLaplacian) {
      getTuftedLaplacianAndMass(M, L, mesh, vpg, mollifyFactor);
    } else {
      M = vpg.vertexLumpedMassMatrix;
      L = vpg.cotanLaplacian;
    }

    // Initialize the inverse mass matrix
    M_inv = (1 / (M.diagonal().array())).matrix().asDiagonal();
    // Alternatively, use the Galerkin mass matrix
    // M = vpg.vertexGalerkinMassMatrix;
    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    // solver.compute(vpg.vertexGalerkinMassMatrix);
    // std::size_t n = mesh.nVertices();
    // Eigen::SparseMatrix<double> I(n, n);
    // I.setIdentity();
    // M_inv = solver.solve(I);

    // Initialize target face/surface areas
    targetFaceAreas = refVpg.faceAreas.reinterpretTo(mesh);
    targetSurfaceArea = targetFaceAreas.raw().sum();

    // Initialize edge length
    targetEdgeLengths = refVpg.edgeLengths.reinterpretTo(mesh);

    // Initialize reference volume
    if (mesh.hasBoundary()) {
      refVolume = 0.0;
    } else {
      refVolume = std::pow(targetSurfaceArea / M_PI / 4, 1.5) * (4 * M_PI / 3); 
    }

    // Initialize surface area
    surfaceArea = vpg.faceAreas.raw().sum();

    // Initialize volume
    for (gcs::Face f : mesh.faces()) {
      volume += signedVolumeFromFace(
          f, vpg, vpg.inputVertexPositions[mesh.vertex(P.ptInd)]);
    }

    // Initialize the vertex position of the last iteration
    pastPositions = vpg.inputVertexPositions;
  }

  /**
   * @brief Destroy the Force object
   *
   * Explicitly unrequire values required by the constructor. In case, there
   * is another pointer to the HalfEdgeMesh and VertexPositionGeometry
   * elsewhere, calculation of dependent quantities should be respected.
   */
  ~Force() {
    vpg.unrequireFaceNormals();
    vpg.unrequireVertexGalerkinMassMatrix();
    vpg.unrequireVertexLumpedMassMatrix();
    vpg.unrequireCotanLaplacian();
    vpg.unrequireFaceAreas();
    vpg.unrequireVertexIndices();
    vpg.unrequireVertexGaussianCurvatures();
    vpg.unrequireFaceIndices();
    vpg.unrequireEdgeLengths();
    vpg.unrequireVertexNormals();
    vpg.unrequireVertexDualAreas();
    vpg.unrequireCornerAngles();
  }

  void getBendingForces();

  void getChemicalPotential();

  void getStretchingForces();

  void getPressureForces();

  void getVesicleForces();

  void getTubeForces();

  void getDPDForces();

  void getExternalForces();

  /**
   * @brief Get velocity from the position of the last iteration
   *
   * @param timeStep
   */
  void getVelocityFromPastPosition(double dt);

  /**
   * @brief Update the vertex position and recompute cached values
   *
   */
  void update_Vertex_positions() { vpg.refreshQuantities(); }

  void pcg_test();
};
} // end namespace ddgsolver
