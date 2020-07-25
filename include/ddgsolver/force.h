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

#include "ddgsolver/macros.h"
#include "ddgsolver/meshops.h"
#include "ddgsolver/util.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

struct Parameters {
  /// Bending modulus
  double Kb;
  /// Spontaneous curvature
  double H0;
  /// Local stretching modulus
  double Ksl;
  /// Global stretching modulus
  double Ksg;
  /// Edge spring constant
  double Kse;
  /// Volume regularization
  double Kv;
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
  double extF;
  /// level of concentration of the external force
  double conc;
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
  /// Cached bending forces
  gcs::VertexData<gc::Vector3> bendingForces;
  /// Cached stretching forces
  gcs::VertexData<gc::Vector3> stretchingForces;
  /// Cached pressure induced forces
  gcs::VertexData<gc::Vector3> pressureForces;
  /// Cached damping forces
  gcs::VertexData<gc::Vector3> dampingForces;
  /// Cached stochastic forces
  gcs::VertexData<gc::Vector3> stochasticForces;
  /// Cached external forces
  gcs::VertexData<gc::Vector3> externalForces;

  /// Cached galerkin mass matrix
  Eigen::SparseMatrix<double> M;
  /// Inverted galerkin mass matrix
  Eigen::SparseMatrix<double> M_inv;
  /// Cotangent Laplacian
  Eigen::SparseMatrix<double> L;
  /// Target area per face
  gcs::FaceData<double> targetFaceAreas;
  /// Target total face area
  double targetSurfaceArea = 0.0;
  /// surface area
  double surfaceArea = 0.0;
  /// Target length per edge
  gcs::EdgeData<double> targetEdgeLength;
  /// Maximal volume
  double maxVolume = 0.0;
  /// Volume
  double volume = 0.0;
  /// Cached vertex positions from the previous step
  gcs::VertexData<gc::Vector3> pastPositions;
  /// Cached vertex velocity by finite differencing past and current position
  gcs::VertexData<gc::Vector3> vel;
  // Mean curvature of the mesh
  Eigen::Matrix<double, Eigen::Dynamic, 3> Hn;
  // Mean curvature of the mesh
  Eigen::Matrix<double, Eigen::Dynamic, 1> H;
  // Spontaneous curvature of the mesh
  Eigen::Matrix<double, Eigen::Dynamic, 3> H0n;
  // Spontaneous curvature of the mesh
  Eigen::Matrix<double, Eigen::Dynamic, 1> H0;
  /// Random number engine
  pcg32 rng;
  std::normal_distribution<double> normal_dist;
  /// magnitude of externally applied force
  Eigen::Matrix<double, Eigen::Dynamic, 1> appliedForceMagnitude;

  /**
   * @brief Construct a new Force object
   *
   * @param mesh_         Mesh connectivity
   * @param vpg_          Embedding and geometry information
   * @param time_step_    Numerical timestep
   */

  Force(gcs::SurfaceMesh &mesh_, gcs::VertexPositionGeometry &vpg_, gcs::VertexPositionGeometry &refVpg_,
        gcs::RichSurfaceMeshData &richData_, Parameters &p)
      : mesh(mesh_), vpg(vpg_), richData(richData_), refVpg(refVpg_), bendingForces(mesh_, {0, 0, 0}), P(p),
        stretchingForces(mesh_, {0, 0, 0}), dampingForces(mesh_, {0, 0, 0}),
        pressureForces(mesh_, {0, 0, 0}), stochasticForces(mesh_, {0, 0, 0}),
        externalForces(mesh_, {0, 0, 0}), vel(mesh_, {0, 0, 0}) {

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

    refVpg.requireFaceAreas();
    refVpg.requireEdgeLengths();

    /// Initialize the mass matrix
    M = vpg.vertexLumpedMassMatrix;
    M_inv = (1 / (M.diagonal().array())).matrix().asDiagonal();
    // M = vpg.vertexGalerkinMassMatrix;
    //// Initialize the inverted Mass matrix
    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    // solver.compute(vpg.vertexGalerkinMassMatrix);
    // std::size_t n = mesh.nVertices();
    // Eigen::SparseMatrix<double> I(n, n);
    // I.setIdentity();
    // M_inv = solver.solve(I);

    // Initialize the conformal Laplacian matrix
    L = vpg.cotanLaplacian;

    // Initialize target face/surface areas
    targetFaceAreas = refVpg.faceAreas.reinterpretTo(mesh);
    targetSurfaceArea = targetFaceAreas.raw().sum();

    // Initialize edge length
    targetEdgeLength = refVpg.edgeLengths.reinterpretTo(mesh);

    // Initialize maximal volume
    double pi = 2 * std::acos(0.0);
    maxVolume = std::pow(targetSurfaceArea / pi / 4, 1.5) * (4 * pi / 3);

    // Initialize surface area
    surfaceArea = vpg.faceAreas.raw().sum();

    // Initialize volume
    for (gcs::Face f : mesh.faces()) {
      volume += signedVolumeFromFace(f, vpg);
    }

    // Initialize the vertex position of the last iteration
    pastPositions = vpg.inputVertexPositions;

    // Initialize mean curvature vector
    Hn.resize(mesh.nVertices(), 3);

    // Initialize spontaneous curvature vector
    H0n.resize(mesh.nVertices(), 3);

    // Initialize the magnitude of externally applied force
    gcs::VertexData<double> geodesicDistanceFromAppliedForce =
        heatMethodDistance(vpg, mesh.vertex(P.ptInd));
    auto &dist_e = geodesicDistanceFromAppliedForce.raw();
    double stdDev = dist_e.maxCoeff() / P.conc;
    appliedForceMagnitude =
        P.extF / (stdDev * pow(pi * 2, 0.5)) *
        (-dist_e.array() * dist_e.array() / (2 * stdDev * stdDev)).exp();
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
  }

  void getBendingForces();

  void getStretchingForces();

  void getPressureForces();

  void getConservativeForces();

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
