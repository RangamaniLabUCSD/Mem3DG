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
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <pcg_random.hpp>
#include <random>

#include <math.h>
#include <vector>

#include "mem3dg/solver/constants.h"
#include "mem3dg/solver/macros.h"
#include "mem3dg/solver/mesh.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/util.h"

using EigenVectorX1D = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using EigenVectorX1D_i = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using EigenVectorX3D =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenTopVec =
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {

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
  /// Ambient Pressure
  double cam;
  /// Temperature
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
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  /// Embedding and other geometric details
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  /// reference embedding geometry
  std::unique_ptr<gcs::VertexPositionGeometry> refVpg;
  /// Cached mesh data
  gcs::RichSurfaceMeshData richData;
  /// Energy
  Energy E;
  /// Time
  double time;

  /// Cached bending stress
  gcs::VertexData<double> bendingPressure;
  /// Cached tension-induced capillary pressure
  gcs::VertexData<double> capillaryPressure;
  /// Cached interfacial line tension
  gcs::VertexData<double> lineTensionPressure;
  gcs::EdgeData<double> lineTension;
  /// Cached externally-applied pressure
  gcs::VertexData<double> externalPressure;
  /// Cached relative inside pressure
  gcs::VertexData<double> insidePressure;
  /// Cached Surface tension
  double surfaceTension;

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

  /// Whether or not do vertex shift
  const bool isVertexShift;
  /// Whether or not consider protein binding
  const bool isProtein;
  /// Whether adopt reduced volume parametrization
  const bool isReducedVolume;
  /// Whether calculate geodesic distance
  const bool isLocalCurvature;

  /// Target area per face
  gcs::FaceData<double> &targetFaceAreas;
  /// Target total face area
  double targetSurfaceArea;
  /// Maximal volume
  double refVolume;
  /// Target length per edge
  gcs::EdgeData<double> &targetEdgeLengths;
  /// Target edge cross length ratio
  gcs::EdgeData<double> targetLcr;
  /// Distance solver
  gcs::HeatMethodDistanceSolver heatSolver;

  /// Cached galerkin mass matrix
  Eigen::SparseMatrix<double> &M;
  /// Inverted galerkin mass matrix
  Eigen::SparseMatrix<double> M_inv;
  /// Cotangent Laplacian
  Eigen::SparseMatrix<double> &L;
  /// Cached geodesic distance
  gcs::VertexData<double> geodesicDistanceFromPtInd;
  /// V-E distribution matrix
  Eigen::SparseMatrix<double> D;

  /// L1 error norm
  double L1ErrorNorm;
  /// surface area
  double surfaceArea;
  /// Volume
  double volume;
  /// Cached vertex positions from the previous step
  gcs::VertexData<gc::Vector3> pastPositions;
  /// Cached vertex velocity by finite differencing past and current position
  gcs::VertexData<gc::Vector3> vel;
  /// Mean curvature of the mesh
  gcs::VertexData<double> H;
  /// Gaussian curvature of the mesh
  gcs::VertexData<double> K;
  /// Spontaneous curvature of the mesh
  gcs::VertexData<double> H0;
  /// Random number engine
  pcg32 rng;
  std::normal_distribution<double> normal_dist;
  /// Indices for vertices chosen for integration
  gcs::VertexData<bool> mask;
  /// "the vertex"
  gcs::Vertex theVertex;

  /// If flipping the edge
  gcs::EdgeData<bool> isFlip;
  /// If splitting the edge
  gcs::EdgeData<bool> isSplit;
  /// If collapsing the edge
  gcs::EdgeData<bool> isCollapse;

  // ==========================================================
  // =============        Constructors           +=============
  // ==========================================================
  /**
   * @brief Construct a new Force object by reading topology and vertex matrices
   *
   * @param topologyMatrix,  topology matrix, F x 3
   * @param vertexMatrix,    input Mesh coordinate matrix, V x 3
   * @param refVertexMatrix, reference mesh coordinate matrix V x 3
   * @param nSub          Number of subdivision
   * @param p             Parameter of simulation
   * @param isReducedVolume Option of whether adopting reduced volume
   * parametrization
   * @param isProtein     Option of considering protein adsorption
   * @param isLocalCurvature Option of whether membrane has local curvature
   * @param isVertexShift Option of whether conducting vertex shift
   * regularization
   */
  System(Eigen::Matrix<double, Eigen::Dynamic, 3> topologyMatrix,
         Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix,
         Eigen::Matrix<double, Eigen::Dynamic, 3> refVertexMatrix, size_t nSub,
         Parameters &p, bool isReducedVolume_, bool isProtein_,
         bool isLocalCurvature_, bool isVertexShift_)
      : System(readMeshes(topologyMatrix, vertexMatrix, refVertexMatrix, nSub),
               p, isReducedVolume_, isProtein_, isLocalCurvature_,
               isVertexShift_){};

  /**
   * @brief Construct a new Force object by reading mesh file path
   *
   * @param inputMesh     Input Mesh
   * @param refMesh       Reference Mesh
   * @param nSub          Number of subdivision
   * @param p             Parameter of simulation
   * @param isReducedVolume Option of whether adopting reduced volume
   * parametrization
   * @param isProtein     Option of considering protein adsorption
   * @param isLocalCurvature Option of whether membrane has local curvature
   * @param isVertexShift Option of whether conducting vertex shift
   * regularization
   */
  System(std::string inputMesh, std::string refMesh, size_t nSub, Parameters &p,
         bool isReducedVolume_, bool isProtein_, bool isLocalCurvature_,
         bool isVertexShift_)
      : System(readMeshes(inputMesh, refMesh, nSub), p, isReducedVolume_,
               isProtein_, isLocalCurvature_, isVertexShift_){};

  /**
   * @brief Construct a new System object by reading netcdf trajectory file path
   *
   * @param trajFile      Netcdf trajectory file
   * @param startingFrame Starting frame for the input mesh
   * @param nSub          Number of subdivision
   * @param p             Parameter of simulation
   * @param isReducedVolume Option of whether adopting reduced volume
   * parametrization
   * @param isProtein     Option of considering protein adsorption
   * @param isLocalCurvature Option of whether membrane has local curvature
   * @param isVertexShift Option of whether conducting vertex shift
   * regularization
   */
  System(std::string trajFile, int startingFrame, size_t nSub, bool isContinue,
         Parameters &p, bool isReducedVolume_, bool isProtein_,
         bool isLocalCurvature_, bool isVertexShift_)
      : System(readTrajFile(trajFile, startingFrame, nSub), p, isReducedVolume_,
               isProtein_, isLocalCurvature_, isVertexShift_) {
    if (isContinue) {
      mapContinuationVariables(trajFile, startingFrame);
    }
  };

  /**
   * @brief Construct a new System object by reading tuple of unique_ptrs
   *
   * @param tuple        Mesh connectivity, Embedding and geometry
   * information, Mesh rich data
   * @param p             Parameter of simulation
   * @param isReducedVolume Option of whether adopting reduced volume
   * parametrization
   * @param isProtein     Option of considering protein adsorption
   * @param isLocalCurvature Option of whether membrane has local
   * curvature
   * @param isVertexShift Option of whether conducting vertex shift
   * regularization
   */
  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>,
                    std::unique_ptr<gcs::VertexPositionGeometry>>
             meshVpgTuple,
         Parameters &p, bool isReducedVolume_, bool isProtein_,
         bool isLocalCurvature_, bool isVertexShift_)
      : System(std::move(std::get<0>(meshVpgTuple)),
               std::move(std::get<1>(meshVpgTuple)),
               std::move(std::get<2>(meshVpgTuple)), p, isReducedVolume_,
               isProtein_, isLocalCurvature_, isVertexShift_){};

  /**
   * @brief Construct a new System object by reading unique_ptrs to mesh and
   * geometry objects
   * @param ptrmesh_         Mesh connectivity
   * @param ptrvpg_          Embedding and geometry information
   * @param ptrRefvpg_       Embedding and geometry information
   * @param p             Parameter of simulation
   * @param isReducedVolume Option of whether adopting reduced volume
   * parametrization
   * @param isProtein     Option of considering protein adsorption
   * @param isLocalCurvature Option of whether membrane has local curvature
   * @param isVertexShift Option of whether conducting vertex shift
   * regularization
   */
  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrrefVpg_, Parameters &p,
         bool isReducedVolume_, bool isProtein_, bool isLocalCurvature_,
         bool isVertexShift_)
      : mesh(std::move(ptrmesh_)), vpg(std::move(ptrvpg_)),
        refVpg(std::move(ptrrefVpg_)), richData(*mesh), P(p), time(0),
        isReducedVolume(isReducedVolume_), isProtein(isProtein_),
        isLocalCurvature(isLocalCurvature_), isVertexShift(isVertexShift_),
        M(vpg->vertexLumpedMassMatrix), L(vpg->cotanLaplacian),
        bendingPressure(*mesh, 0), insidePressure(*mesh, 0), D(),
        capillaryPressure(*mesh, 0), lineTension(*mesh, 0),
        lineTensionPressure(*mesh, 0), chemicalPotential(*mesh, 0),
        externalPressure(*mesh, 0), regularizationForce(*mesh, {0, 0, 0}),
        targetLcr(*mesh), targetEdgeLengths(refVpg->edgeLengths),
        targetFaceAreas(refVpg->faceAreas), stochasticForce(*mesh, {0, 0, 0}),
        dampingForce(*mesh, {0, 0, 0}), proteinDensity(*mesh, 0),
        vel(*mesh, {0, 0, 0}), isFlip(*mesh, false), isSplit(*mesh, false),
        isCollapse(*mesh, false), E({0, 0, 0, 0, 0, 0, 0, 0, 0}),
        heatSolver(*vpg) {

    // GC computed properties
    vpg->requireFaceNormals();
    vpg->requireVertexLumpedMassMatrix();
    vpg->requireCotanLaplacian();
    vpg->requireFaceAreas();
    vpg->requireVertexIndices();
    vpg->requireVertexGaussianCurvatures();
    vpg->requireFaceIndices();
    vpg->requireEdgeLengths();
    vpg->requireVertexNormals();
    vpg->requireVertexDualAreas();
    vpg->requireCornerAngles();
    vpg->requireCornerScaledAngles();
    vpg->requireDECOperators();
    vpg->requireEdgeDihedralAngles();
    // vpg->requireVertexTangentBasis();

    // Check confliciting parameters and options
    checkParameters();

    // Initialize reference values
    initConstants();

    // Regularize the vetex position geometry if needed
    if (isVertexShift) {
      vertexShift();
    }

    /// compute nonconstant values during simulation
    updateVertexPositions();
  }

  /**
   * @brief Destroy the Force object
   *
   * Explicitly unrequire values required by the constructor. In case, there
   * is another pointer to the HalfEdgeMesh and VertexPositionGeometry
   * elsewhere, calculation of dependent quantities should be respected.
   */
  ~System() {
    vpg->unrequireFaceNormals();
    vpg->unrequireVertexLumpedMassMatrix();
    vpg->unrequireCotanLaplacian();
    vpg->unrequireFaceAreas();
    vpg->unrequireVertexIndices();
    vpg->unrequireVertexGaussianCurvatures();
    vpg->unrequireFaceIndices();
    vpg->unrequireEdgeLengths();
    vpg->unrequireVertexNormals();
    vpg->unrequireVertexDualAreas();
    vpg->unrequireCornerAngles();
    vpg->unrequireCornerScaledAngles();
  }

  // ==========================================================
  // ================     Initialization     ==================
  // ==========================================================

  /**
   * @brief Construct a tuple of unique_ptrs from topology matrix and vertex
   * position matrix
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMeshes(Eigen::Matrix<double, Eigen::Dynamic, 3> faceVertexMatrix,
             Eigen::Matrix<double, Eigen::Dynamic, 3> vertexPositionMatrix,
             Eigen::Matrix<double, Eigen::Dynamic, 3> refVertexPositionMatrix,
             size_t nSub);

  /**
   * @brief Construct a tuple of unique_ptrs from mesh and refMesh path
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMeshes(std::string inputMesh, std::string refMesh, size_t nSub);

#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Construct a tuple of unique_ptrs from netcdf path
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readTrajFile(std::string trajFile, int startingFrame, size_t nSub);
  /**
   * @brief Map the continuation variables
   *
   */
  void mapContinuationVariables(std::string trajFile, int startingFrame);
#endif
  /**
   * @brief Check all conflicting parameters and options
   *
   */
  void checkParameters();

  /**
   * @brief testing of random number generator pcg
   *
   */
  void pcg_test();

  /**
   * @brief Initialize all constant values (on refVpg) needed for computation
   *
   */
  void initConstants();

  /**
   * @brief Update the vertex position and recompute cached values
   * (all quantities that characterizes the current energy state)
   * Careful when using eigenMap: memory address may change after update!!
   */
  void updateVertexPositions();

  /**
   * @brief Visualization using Polyscope
   *
   */
  void visualize();

  // ==========================================================
  // ================        Pressure        ==================
  // ==========================================================
  /**
   * @brief Compute bending pressure component of the system
   */
  EigenVectorX1D computeBendingPressure();

  /**
   * @brief Compute chemical potential of the system
   */
  EigenVectorX1D computeChemicalPotential();

  /**
   * @brief Compute capillary pressure component of the system
   */
  EigenVectorX1D computeCapillaryPressure();

  /**
   * @brief Compute inside pressure component of the system
   */
  EigenVectorX1D computeInsidePressure();

  /**
   * @brief Compute line tension pressure component of the system
   */
  EigenVectorX1D computeLineTensionPressure();

  /**
   * @brief Compute DPD forces of the system
   */
  std::tuple<EigenVectorX3D, EigenVectorX3D> computeDPDForces();

  /**
   * @brief Compute external pressure component of the system
   */
  EigenVectorX1D computeExternalPressure();

  /**
   * @brief Compute all forces of the system
   */
  void computePhysicalForces();

  // ==========================================================
  // ================        Energy          ==================
  // ==========================================================
  /**
   * @brief Compute bending energy
   */
  void computeBendingEnergy();

  /**
   * @brief Compute surface energy
   */
  void computeSurfaceEnergy();

  /**
   * @brief Compute pressure work
   */
  void computePressureEnergy();

  /**
   * @brief Compute chemical energy
   */
  void computeChemicalEnergy();

  /**
   * @brief Compute line tension energy
   */
  void computeLineTensionEnergy();

  /**
   * @brief Compute external force energy
   */
  void computeExternalForceEnergy();

  /**
   * @brief Compute kinetic energy
   */
  void computeKineticEnergy();

  /**
   * @brief Compute potential energy
   */
  void computePotentialEnergy();

  /**
   * @brief Compute all components of energy (free energy)
   */
  void computeFreeEnergy();

  /**
   * @brief Compute the L1 norm of the pressure
   */
  double computeL1Norm(Eigen::Matrix<double, Eigen::Dynamic, 1> &force) const;
  double computeL1Norm(Eigen::Matrix<double, Eigen::Dynamic, 1> &&force) const;

  // ==========================================================
  // =============        Regularization        ===============
  // ==========================================================
  /**
   * @brief Apply vertex shift by moving the vertices chosen for integration to
   * the Barycenter of the it neighbors
   */
  void vertexShift();

  /**
   * @brief Compute regularization pressure component of the system
   */
  void computeRegularizationForce();

  /**
   * @brief Edge flip if not Delaunay
   */
  void edgeFlip();

  /**
   * @brief Get regularization pressure component of the system
   */
  void growMesh();

  /**
   * @brief Get length cross ratio of the mesh
   */
  gcs::EdgeData<double>
  computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg) const;
};
} // namespace mem3dg
