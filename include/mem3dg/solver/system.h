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
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <pcg_random.hpp>
#include <random>

#include <math.h>

#include "mem3dg/solver/constants.h"
#include "mem3dg/solver/macros.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/util.h"

#include <vector>

namespace mem3dg {

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
  std::vector<double> r_H0;
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
  /// Line tension
  double eta;
  /// binding energy per protein
  double epsilon;
  /// binding constant
  double Bc;
  /// Dissipation coefficient
  double gamma;
  /// Reduced volume
  double Vt;
  /// Boltzmann constant*Temperature
  double temp;
  /// Noise
  double sigma;
  /// index of node with applied external force
  std::vector<double> pt;
  /// Magnitude of external force
  double Kf;
  /// level of concentration of the external force
  double conc;
  /// target height
  double height;
  /// domain of integration
  double radius;
  /// augmented Lagrangian parameter for area
  double lambdaSG = 0;
  /// augmented Lagrangian parameter for volume
  double lambdaV = 0;
};

struct Energy {
  /// total Energy of the system
  double totalE;
  /// kinetic energy of the membrane
  double kE;
  /// potential energy of the membrane
  double potE;
  /// bending energy of the membrane
  double BE;
  /// stretching energy of the membrane
  double sE;
  /// work of pressure within membrane
  double pE;
  /// chemical energy of the membrane protein
  double cE;
  /// line tension energy of interface
  double lE;
  /// work of external force
  double exE;
};

class DLL_PUBLIC System {
public:
  /// Parameters
  Parameters P;
  /// Cached mesh of interest
  gcs::ManifoldSurfaceMesh &mesh;
  /// Cached mesh data
  gcs::RichSurfaceMeshData &richData;
  /// Embedding and other geometric details
  gcs::VertexPositionGeometry &vpg;
  /// reference embedding geometry
  gcs::VertexPositionGeometry &refVpg;
  /// Energy
  Energy E;

  /// Cached bending stress
  gcs::VertexData<gc::Vector3> bendingPressure;
  /// Cached tension-induced capillary pressure
  gcs::VertexData<gc::Vector3> capillaryPressure;
  /// Cached interfacial line tension
  gcs::VertexData<gc::Vector3> lineTensionPressure;
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
  const bool isTuftedLaplacian;
  /// Whether or not do vertex shift
  const bool isVertexShift;
  /// Whether or not consider protein binding
  const bool isProtein;
  /// Whether adopt reduced volume parametrization
  const bool isReducedVolume;

  /// Target area per face
  gcs::FaceData<double> targetFaceAreas;
  /// Target total face area
  double targetSurfaceArea;
  /// Maximal volume
  double refVolume;
  /// Target length per edge
  gcs::EdgeData<double> targetEdgeLengths;
  /// Target edge cross length ratio
  gcs::EdgeData<double> targetLcr;
  /// Distance solver
  gcs::HeatMethodDistanceSolver heatSolver;

  /// Cached galerkin mass matrix
  Eigen::SparseMatrix<double> M;
  /// Inverted galerkin mass matrix
  Eigen::SparseMatrix<double> M_inv;
  /// Cotangent Laplacian
  Eigen::SparseMatrix<double> L;
  /// Cached geodesic distance
  gcs::VertexData<double> geodesicDistanceFromPtInd;

  /// L2 error norm
  double L2ErrorNorm;
  /// surface area
  double surfaceArea;
  /// Volume
  double volume;
  /// Interface Area;
  double interArea;
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
  /// indices for vertices chosen for integration
  Eigen::Matrix<bool, Eigen::Dynamic, 1> mask;
  /// "the point" index
  size_t ptInd;

  /*
   * @brief Construct a new Force object
   *
   * @param mesh_         Mesh connectivity
   * @param vpg_          Embedding and geometry information
   * @param time_step_    Numerical timestep
   */

  System(gcs::ManifoldSurfaceMesh &mesh_, gcs::VertexPositionGeometry &vpg_,
         gcs::VertexPositionGeometry &refVpg_,
         gcs::RichSurfaceMeshData &richData_, Parameters &p,
         bool isReducedVolume_, bool isProtein_, bool isVertexShift_,
         bool isTuftedLaplacian_)
      : mesh(mesh_), vpg(vpg_), richData(richData_), refVpg(refVpg_), P(p),
        isTuftedLaplacian(isTuftedLaplacian_),
        isReducedVolume(isReducedVolume_), isProtein(isProtein_),
        isVertexShift(isVertexShift_), bendingPressure(mesh_, {0, 0, 0}),
        insidePressure(mesh_, {0, 0, 0}), capillaryPressure(mesh_, {0, 0, 0}),
        lineTensionPressure(mesh_, {0, 0, 0}), chemicalPotential(mesh_, 0.0),
        externalPressure(mesh_, {0, 0, 0}),
        regularizationForce(mesh_, {0, 0, 0}), targetLcr(mesh_),
        stochasticForce(mesh_, {0, 0, 0}), dampingForce(mesh_, {0, 0, 0}),
        proteinDensity(mesh_, 0), vel(mesh_, {0, 0, 0}),
        E({0, 0, 0, 0, 0, 0, 0, 0, 0}), heatSolver(vpg) {

    // Initialize RNG
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    rng = pcg32(seed_source);

    // GC computed properties
    vpg.requireFaceNormals();
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
    vpg.requireCornerScaledAngles();
    // vpg.requireVertexTangentBasis();

    /// compute constant values during simulation
    // Find the closest point index to P.pt in refVpg
    closestPtIndToPt(mesh, refVpg, P.pt, ptInd);

    // Initialize the constant mask based on distance from the point specified
    // Or mask boundary element
    mask = (heatMethodDistance(refVpg, mesh.vertex(ptInd)).raw().array() <
            P.radius)
               .matrix();
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

    // Regularize the vetex position geometry if needed
    if (isVertexShift) {
      vertexShift(mesh, vpg, mask);
      update_Vertex_positions();
    }

    /// compute nonconstant values during simulation
    update_Vertex_positions();
  }

  /**
   * @brief Destroy the Force object
   *
   * Explicitly unrequire values required by the constructor. In case, there
   * is another pointer to the HalfEdgeMesh and VertexPositionGeometry
   * elsewhere, calculation of dependent quantities should be respected.
   */
  ~System() {
    vpg.unrequireFaceNormals();
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
    vpg.unrequireCornerScaledAngles();
  }


  // ==========================================================
  // ================        Update          ==================
  // ==========================================================
  /**
   * @brief Update the vertex position and recompute cached values
   * (all quantities that characterizes the current energy state)
   *
   */
  void update_Vertex_positions();


  // ==========================================================
  // ================        Pressure        ==================
  // ==========================================================
  /**
   * @brief Get bending pressure component of the system
   */
  void getBendingPressure();

  /**
   * @brief Get chemical potential of the system
   */
  void getChemicalPotential();

  /**
   * @brief Get capillary pressure component of the system
   */
  void getCapillaryPressure();

  /**
   * @brief Get inside pressure component of the system
   */
  void getInsidePressure();

  /**
   * @brief Get regularization pressure component of the system
   */
  void getRegularizationForce();

  /**
   * @brief Get line tension pressure component of the system
   */
  void getLineTensionPressure();

  /**
   * @brief Get DPD forces of the system
   */
  void getDPDForces();

  /**
   * @brief testing of random number generator pcg
   *
   */
  void pcg_test();

  /**
   * @brief Get velocity from the position of the last iteration
   *
   * @param dt timestep
   */
  void getVelocityFromPastPosition(double dt);

  /**
   * @brief Get external pressure component of the system
   */
  void getExternalPressure();

  /**
   * @brief Get all forces of the system
   */
  void getAllForces();

  // ==========================================================
  // ================        Energy          ==================
  // ==========================================================
  /**
   * @brief Get bending energy
   */
  void getBendingEnergy();

  /**
   * @brief Get surface energy
   */
  void getSurfaceEnergy();

  /**
   * @brief Get pressure work
   */
  void getPressureEnergy();

  /**
   * @brief Get chemical energy
   */
  void getChemicalEnergy();

  /**
   * @brief Get line tension energy
   */
  void getLineTensionEnergy();

  /**
   * @brief Get external force energy
   */
  void getExternalForceEnergy();

  /**
   * @brief Get kinetic energy
   */
  void getKineticEnergy();

  /**
   * @brief Get potential energy
   */
  void getPotentialEnergy();

  /**
   * @brief Get all components of energy (free energy)
   */
  void getFreeEnergy();

  /**
   * @brief Get the L2 norm of the force (pressure), which is the residual of
   * the PDE
   */
  void
  getL2ErrorNorm(Eigen::Matrix<double, Eigen::Dynamic, 3> physicalPressure);
};
} // namespace mem3dg
