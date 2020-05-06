#pragma once

#include <cassert>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <pcg_random.hpp>
#include <random>

#include "util.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class Force {
public:
  /// Cached mesh of interest
  gcs::HalfedgeMesh& mesh;
  /// Embedding and other geometric details
  gcs::VertexPositionGeometry &vpg;
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
  /// Target volume
  double targetVolume = 0.0;
  /// Cached vertex positions from the previous step
  gcs::VertexData<gc::Vector3> pastPositions;
  /// Cached vertex velocity by finite differecing past and current position
  gcs::VertexData<gc::Vector3> vertexVelocity;
  /// Random numer engine
  pcg32 rng;
  std::normal_distribution<double> normal_dist;

  /**
   * @brief Construct a new Force object
   *
   * @param mesh_         Mesh connectivity
   * @param vpg_          Embedding and goemetry information
   * @param time_step_    Numerical timestep
   */
  Force(gcs::HalfedgeMesh &mesh_, gcs::VertexPositionGeometry &vpg_)
      : mesh(mesh_), vpg(vpg_),
        bendingForces(mesh_, {0, 0, 0}), stretchingForces(mesh_, {0, 0, 0}),
        dampingForces(mesh_, {0, 0, 0}), pressureForces(mesh_, {0, 0, 0}),
    stochasticForces(mesh_, { 0, 0, 0 }), vertexVelocity(mesh_, {0, 0, 0 }) {

    // Initialize RNG
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    rng = pcg32(seed_source);

    // GC computed properties
    vpg.requireFaceNormals();
    vpg.requireVertexGalerkinMassMatrix();
    vpg.requireCotanLaplacian();
    vpg.requireFaceAreas();
    vpg.requireVertexIndices();
    vpg.requireVertexGaussianCurvatures();
    vpg.requireFaceIndices();

    // Initialize the mass matrix
    M = vpg.vertexGalerkinMassMatrix;
    // Initialize the inverted Mass matrix
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M);
    std::size_t n = mesh.nVertices();
    Eigen::SparseMatrix<double> I(n, n);
    I.setIdentity();
    M_inv = solver.solve(I);

    // Initialize the conformal Laplacian matrix
    L = vpg.cotanLaplacian;

    // Initialize face areas
    targetFaceAreas = vpg.faceAreas;
    auto faceAreas_e = EigenMap(targetFaceAreas);
    targetSurfaceArea = faceAreas_e.sum();

    // Initialize initial volume
    for (gcs::Face f : mesh.faces()) {
      targetVolume += signedVolumeFromFace(f, vpg);
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
    vpg.unrequireFaceAreas();
    vpg.unrequireVertexGalerkinMassMatrix();
    vpg.unrequireCotanLaplacian();
    vpg.unrequireFaceNormals();
    vpg.unrequireVertexIndices();
    vpg.unrequireVertexGaussianCurvatures();
  }

  /**
   * @brief Compute the bending force
   *
   * @param Kb
   * @param H0
   */
  void getBendingForces(double &Kb, double &H0);
  // void bending_force(double Kb, Eigen::Matrix<double, Eigen::Dynamic, 1>
  // H0);

  /**
   * @brief Compute force from stretching
   *
   * @param Ksl
   * @param Ksg
   */
  void getStretchingForces(double &Ksl, double &Ksg);

  /**
   * @brief Compute forces from pressure
   *
   * @param Kv
   * @param Vt
   */
  void getPressureForces(double &Kv, double &Vt);

  /**
   * @brief Compute forces from damping
   *
   * @param gamma
   */
  void getDampingForces(double &gamma);

  /**
  * @brief Get velocity from the position of the last iteration
  *
  * @param timeStep
  */
  void getVelocityFromPastPosition(double& timeStep);
  /**
   * @brief Compute forces from random noise
   *
   * @param sigma
   */
  void getStochasticForces(double &sigma);

  /**
   * @brief Get volume from a face
   *
   * @param f
   * @param vpg
   * @return double
   */
  double signedVolumeFromFace(gcs::Face &f, gcs::VertexPositionGeometry &vpg);

  /**
   * @brief Get the vector from halfedge vertices
   *
   * @param he
   * @param vpg
   * @return gc::Vector3
   */
  inline gc::Vector3 vecFromHalfedge(gcs::Halfedge &he,
                                     gcs::VertexPositionGeometry &vpg) {
    return vpg.inputVertexPositions[he.next().vertex()] -
           vpg.inputVertexPositions[he.vertex()];
  }

  /**
   * @brief Update the vertex position and recompute cached values
   *
   */
  void update_Vertex_positions() { vpg.refreshQuantities(); }

  void pcg_test();
};
} // end namespace ddgsolver
