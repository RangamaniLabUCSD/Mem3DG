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

// #include <cassert>

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/heat_method_distance.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <pcg_random.hpp>
#include <random>

#include <functional>
#include <iomanip>
#include <math.h>
#include <vector>

#include "mem3dg/constants.h"
#include "mem3dg/macros.h"
#include "mem3dg/mesh_io.h"
#include "mem3dg/meshops.h"
#include "mem3dg/type_utilities.h"

#include "mem3dg/solver/forces.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/solver/parameters.h"
#ifdef MEM3DG_WITH_NETCDF
#include "mem3dg/solver/mutable_trajfile.h"
#include "mem3dg/solver/trajfile.h"
#endif

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {

namespace solver {

struct Energy {
  /// time
  double time = 0;
  /// total Energy of the system
  double totalEnergy = 0;
  /// kinetic energy of the membrane
  double kineticEnergy = 0;
  /// potential energy of the membrane
  double potentialEnergy = 0;
  /// spontaneous curvature energy of the membrane
  double spontaneousCurvatureEnergy = 0;
  /// deviatoric curvature energy of the membrane
  double deviatoricCurvatureEnergy = 0;
  /// area difference energy of the membrane
  double areaDifferenceEnergy = 0;
  /// stretching energy of the membrane
  double surfaceEnergy = 0;
  /// work of pressure within membrane
  double pressureEnergy = 0;
  /// adsorption energy of the membrane protein
  double adsorptionEnergy = 0;
  /// aggregation energy of the membrane protein
  double aggregationEnergy = 0;
  /// entropy energy of the membrane protein
  double entropyEnergy = 0;
  /// line tension energy of interface
  double dirichletEnergy = 0;
  /// work of external force
  double externalWork = 0;
  /// protein interior penalty energy
  double proteinInteriorPenalty = 0;
  /// membrane self-avoidance penalty energy
  double selfAvoidancePenalty = 0;
  /// mesh edge spring energy
  double edgeSpringEnergy = 0;
  /// mesh face spring energy
  double faceSpringEnergy = 0;
  /// mesh LCR spring energy
  double lcrSpringEnergy = 0;
};

class DLL_PUBLIC System {
protected:
  /// Random number engine
  pcg32 rng;
  std::normal_distribution<double> normal_dist;

public:
  /// Parameters
  Parameters parameters;
  /// Mesh processor
  MeshProcessor meshProcessor;

  /// Cached mesh of interest
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  /// Embedding and other geometric details
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  /// Embedding and other geometric details of reference mesh
  std::unique_ptr<gcs::VertexPositionGeometry> refVpg;
  /// reference length cross ratio
  gcs::HalfedgeData<double> refLcrs;
  /// Energy
  Energy energy;
  /// Time
  double time;

  /// Forces of the system
  Forces forces;

  /// mechanical error norm
  double mechErrorNorm;
  /// chemical error norm
  double chemErrorNorm;
  /// surface area
  double surfaceArea;
  /// Volume
  double volume;
  /// Cached protein surface density
  gcs::VertexData<double> proteinDensity;
  /// Spontaneous curvature gradient of the mesh
  gcs::FaceData<gc::Vector3> proteinDensityGradient;
  /// Cached vertex velocity
  gcs::VertexData<gc::Vector3> velocity;
  /// Cached vertex protein rate of change
  gcs::VertexData<double> proteinRateOfChange;
  /// Spontaneous curvature of the mesh
  gcs::VertexData<double> H0;
  /// Bending rigidity of the membrane
  gcs::VertexData<double> Kb;
  /// deviatoric rigidity of the membrane
  gcs::VertexData<double> Kd;

  /// Cached geodesic distance
  gcs::VertexData<double> geodesicDistance;
  /// is Smooth
  bool isSmooth;
  /// if being mutated
  gcs::VertexData<bool> mutationMarker;
  /// if has boundary
  bool isOpenMesh;
  /// "the vertex"
  gcs::SurfacePoint center;
  gcs::VertexData<bool> centerTracker;
  /// projected time of collision
  double projectedCollideTime;

  // ==========================================================
  // =============        Constructors           ==============
  // ==========================================================

  // =======================================
  // =======       Matrices         ========
  // =======================================
  System(EigenVectorX3sr &topologyMatrix, EigenVectorX3dr &vertexMatrix,
         EigenVectorX3dr &referenceVertexMatrx, EigenVectorX1d &proteinDensity_,
         EigenVectorX3dr &velocity_, Parameters &p, double time_ = 0)
      : System(readMatrices(topologyMatrix, vertexMatrix, referenceVertexMatrx),
               proteinDensity_, velocity_, p, time_){};

  System(EigenVectorX3sr &topologyMatrix, EigenVectorX3dr &vertexMatrix,
         EigenVectorX1d &proteinDensity_, EigenVectorX3dr &velocity_,
         Parameters &p, double time_ = 0)
      : System(readMatrices(topologyMatrix, vertexMatrix), proteinDensity_,
               velocity_, p, time_){};

  System(EigenVectorX3sr &topologyMatrix, EigenVectorX3dr &vertexMatrix,
         Parameters &p, double time_ = 0)
      : System(readMatrices(topologyMatrix, vertexMatrix), p, time_){};

  System(EigenVectorX3sr &topologyMatrix, EigenVectorX3dr &vertexMatrix,
         EigenVectorX1d &proteinDensity_, EigenVectorX3dr &velocity_,
         double time_ = 0)
      : System(readMatrices(topologyMatrix, vertexMatrix), proteinDensity_,
               velocity_, time_){};

  System(EigenVectorX3sr &topologyMatrix, EigenVectorX3dr &vertexMatrix,
         double time_ = 0)
      : System(readMatrices(topologyMatrix, vertexMatrix), time_){};

  // =======================================
  // =======       Mesh Files       ========
  // =======================================
  System(std::string inputMesh, std::string referenceMesh,
         EigenVectorX1d &proteinDensity_, EigenVectorX3dr &velocity_,
         Parameters &p, double time_ = 0)
      : System(readMeshFile(inputMesh, referenceMesh), proteinDensity_,
               velocity_, p, time_){};

  System(std::string inputMesh, EigenVectorX1d &proteinDensity_,
         EigenVectorX3dr &velocity_, Parameters &p, double time_ = 0)
      : System(readMeshFile(inputMesh), proteinDensity_, velocity_, p, time_){};

  System(std::string inputMesh, Parameters &p, double time_ = 0)
      : System(readMeshFile(inputMesh), p, time_){};

  System(std::string inputMesh, EigenVectorX1d &proteinDensity_,
         EigenVectorX3dr &velocity_, double time_ = 0)
      : System(readMeshFile(inputMesh), proteinDensity_, velocity_, time_){};

  System(std::string inputMesh, double time_ = 0)
      : System(readMeshFile(inputMesh), time_){};

  // =======================================
  // =======       NetCDF Files     ========
  // =======================================
#ifdef MEM3DG_WITH_NETCDF
  System(std::string trajFile, int startingFrame, Parameters &p)
      : System(readTrajFile(trajFile, startingFrame), p){};

  System(std::string trajFile, int startingFrame)
      : System(readTrajFile(trajFile, startingFrame)){};

#endif

  // =======================================
  // =======       Tuple            ========
  // =======================================
  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>,
                    std::unique_ptr<gcs::VertexPositionGeometry>,
                    EigenVectorX1d, EigenVectorX3dr, double>
             initialConditionsTuple,
         Parameters &p)
      : System(std::move(std::get<0>(initialConditionsTuple)),
               std::move(std::get<1>(initialConditionsTuple)),
               std::move(std::get<2>(initialConditionsTuple)),
               std::get<3>(initialConditionsTuple),
               std::get<4>(initialConditionsTuple), p,
               std::get<5>(initialConditionsTuple)){};

  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>,
                    EigenVectorX1d, EigenVectorX3dr, double>
             initialConditionsTuple,
         Parameters &p)
      : System(std::move(std::get<0>(initialConditionsTuple)),
               std::move(std::get<1>(initialConditionsTuple)),
               std::get<2>(initialConditionsTuple),
               std::get<3>(initialConditionsTuple), p,
               std::get<4>(initialConditionsTuple)){};

  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>,
                    std::unique_ptr<gcs::VertexPositionGeometry>,
                    EigenVectorX1d, EigenVectorX3dr, double>
             initialConditionsTuple)
      : System(std::move(std::get<0>(initialConditionsTuple)),
               std::move(std::get<1>(initialConditionsTuple)),
               std::move(std::get<2>(initialConditionsTuple)),
               std::get<3>(initialConditionsTuple),
               std::get<4>(initialConditionsTuple),
               std::get<5>(initialConditionsTuple)){};

  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>,
                    EigenVectorX1d, EigenVectorX3dr, double>
             initialConditionsTuple)
      : System(std::move(std::get<0>(initialConditionsTuple)),
               std::move(std::get<1>(initialConditionsTuple)),
               std::get<2>(initialConditionsTuple),
               std::get<3>(initialConditionsTuple),
               std::get<4>(initialConditionsTuple)){};

  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>,
                    std::unique_ptr<gcs::VertexPositionGeometry>>
             meshVpgTuple,
         EigenVectorX1d &proteinDensity_, EigenVectorX3dr &velocity_,
         Parameters &p, double time_ = 0)
      : System(std::move(std::get<0>(meshVpgTuple)),
               std::move(std::get<1>(meshVpgTuple)),
               std::move(std::get<2>(meshVpgTuple)), proteinDensity_, velocity_,
               p, time_){};

  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>>
             meshVpgTuple,
         EigenVectorX1d &proteinDensity_, EigenVectorX3dr &velocity_,
         Parameters &p, double time_ = 0)
      : System(std::move(std::get<0>(meshVpgTuple)),
               std::move(std::get<1>(meshVpgTuple)), proteinDensity_, velocity_,
               p, time_){};

  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>>
             meshVpgTuple,
         Parameters &p, double time_ = 0)
      : System(std::move(std::get<0>(meshVpgTuple)),
               std::move(std::get<1>(meshVpgTuple)), p, time_){};

  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>>
             meshVpgTuple,
         EigenVectorX1d &proteinDensity_, EigenVectorX3dr &velocity_,
         double time_ = 0)
      : System(std::move(std::get<0>(meshVpgTuple)),
               std::move(std::get<1>(meshVpgTuple)), proteinDensity_, velocity_,
               time_){};

  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>>
             meshVpgTuple,
         double time_ = 0)
      : System(std::move(std::get<0>(meshVpgTuple)),
               std::move(std::get<1>(meshVpgTuple)), time_){};

  // =======================================
  // =======    Geometry Central    ========
  // =======================================
  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrrefvpg_,
         EigenVectorX1d &proteinDensity_, EigenVectorX3dr &velocity_,
         Parameters &p, double time_ = 0)
      : System(std::move(ptrmesh_), std::move(ptrvpg_), std::move(ptrrefvpg_),
               proteinDensity_, velocity_, time_) {
    parameters = p;
  }

  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrrefvpg_,
         EigenVectorX1d &proteinDensity_, EigenVectorX3dr &velocity_,
         double time_ = 0)
      : System(std::move(ptrmesh_), std::move(ptrvpg_), proteinDensity_,
               velocity_, time_) {
    refVpg = std::move(ptrrefvpg_);
  }

  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_,
         EigenVectorX1d &proteinDensity_, EigenVectorX3dr &velocity_,
         Parameters &p, double time_ = 0)
      : System(std::move(ptrmesh_), std::move(ptrvpg_), proteinDensity_,
               velocity_, time_) {
    parameters = p;
  }

  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_,
         EigenVectorX1d &proteinDensity_, EigenVectorX3dr &velocity_,
         double time_ = 0)
      : System(std::move(ptrmesh_), std::move(ptrvpg_), time_) {
    proteinDensity.raw() = proteinDensity_;
    toMatrix(velocity) = velocity_;
  }

  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_, Parameters &p,
         double time_ = 0)
      : System(std::move(ptrmesh_), std::move(ptrvpg_), time_) {
    parameters = p;
  }

  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_, double time_ = 0)
      : mesh(std::move(ptrmesh_)), vpg(std::move(ptrvpg_)), forces(*mesh, *vpg),
        time(time_) {
    refLcrs = gcs::HalfedgeData<double>(*mesh);
    refVpg = vpg->copy();

    energy = Energy({time, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    proteinDensity = gc::VertexData<double>(*mesh, 1);
    proteinDensityGradient = gcs::FaceData<gc::Vector3>(*mesh, {0, 0, 0});
    velocity = gcs::VertexData<gc::Vector3>(*mesh, {0, 0, 0});
    proteinRateOfChange = gcs::VertexData<double>(*mesh, 0);
    H0 = gcs::VertexData<double>(*mesh);
    Kb = gcs::VertexData<double>(*mesh);
    Kd = gcs::VertexData<double>(*mesh);

    chemErrorNorm = 0;
    mechErrorNorm = 0;

    geodesicDistance = gcs::VertexData<double>(*mesh, 0);

    isSmooth = true;
    mutationMarker = gc::VertexData<bool>(*mesh, false);
    centerTracker = gc::VertexData<bool>(*mesh, false);

    // GC computed properties
    vpg->requireFaceNormals();
    vpg->requireVertexLumpedMassMatrix();
    vpg->requireCotanLaplacian();
    vpg->requireFaceAreas();
    vpg->requireVertexIndices();
    vpg->requireVertexGaussianCurvatures();
    vpg->requireVertexMeanCurvatures();
    vpg->requireFaceIndices();
    vpg->requireEdgeLengths();
    vpg->requireVertexNormals();
    vpg->requireVertexDualAreas();
    vpg->requireCornerAngles();
    vpg->requireCornerScaledAngles();
    vpg->requireDECOperators();
    vpg->requireEdgeDihedralAngles();
    vpg->requireHalfedgeCotanWeights();
    vpg->requireEdgeCotanWeights();
    // vpg->requireVertexTangentBasis();
  }

public:
  /**
   * @brief Destroy the System
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
    vpg->unrequireVertexMeanCurvatures();
    vpg->unrequireFaceIndices();
    vpg->unrequireEdgeLengths();
    vpg->unrequireVertexNormals();
    vpg->unrequireVertexDualAreas();
    vpg->unrequireCornerAngles();
    vpg->unrequireCornerScaledAngles();
    vpg->unrequireDECOperators();
    vpg->unrequireEdgeDihedralAngles();
    vpg->unrequireHalfedgeCotanWeights();
    vpg->unrequireEdgeCotanWeights();
  }

  // ==========================================================
  // ================          io.cpp        ==================
  // ==========================================================

  /**
   * @brief Construct a tuple of unique_ptrs from topology matrix and vertex
   * position matrix
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMatrices(EigenVectorX3sr &faceVertexMatrix,
               EigenVectorX3dr &vertexPositionMatrix);

  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMatrices(EigenVectorX3sr &faceVertexMatrix,
               EigenVectorX3dr &vertexPositionMatrix,
               EigenVectorX3dr &refVertexPositionMatrix);

  /**
   * @brief Construct a tuple of unique_ptrs from mesh and refMesh path
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMeshFile(std::string inputMesh);

  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMeshFile(std::string inputMesh, std::string refMesh);

  /**
   * @brief Map the continuation variables
   *
   */
  void mapContinuationVariables(std::string plyFile);

  /**
   * @brief Save RichData to .ply file
   *
   */
  void saveRichData(std::string PathToSave, bool isJustGeometry = false);

#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Construct a tuple of unique_ptrs from netcdf path
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>,
             std::unique_ptr<gcs::VertexPositionGeometry>, EigenVectorX1d,
             EigenVectorX3dr, double>
  readTrajFile(std::string trajFile, int startingFrame);
#endif

  // ==========================================================
  // ================     init.cpp           ==================
  // ==========================================================
  /**
   * @brief Check all conflicting parameters and options
   *
   */
  void checkConfiguration();

  /**
   * @brief Initialize system
   *
   */
  void initialize(std::size_t nMutation = 0, bool ifMute = false);

  /**
   * @brief Initialize all constant values (on refVpg) needed for computation
   *
   */
  void initializeConstants(bool ifMute);

  /**
   * @brief Update the vertex position and recompute cached values
   * (all quantities that characterizes the current energy state)
   * Careful: 1. when using eigenMap: memory address may change after update!!
   */
  void updateConfigurations();

  /**
   * @brief Update the reference edge length, face area, and length cross ratio
   */
  void updateReferenceConfigurations();

  // ==========================================================
  // ================  variation_vector.cpp  ==================
  // ==========================================================
  /**
   * @brief template code for populate verttexwise using halfedge vector
   * computation
   */
  static gcs::VertexData<gc::Vector3> halfedgeVectorToVertexVector(
      gcs::ManifoldSurfaceMesh &mesh, gcs::VertexPositionGeometry &vpg,
      std::function<gc::Vector3(gcs::VertexPositionGeometry &, gc::Halfedge &)>
          computeHalfedgeVariationalVector);

  /**
   * @brief Compute vertex volume variation vector
   */
  gcs::VertexData<gc::Vector3> computeVertexVolumeVariationVectors();

  /**
   * @brief Compute halfedge volume variation vector
   */
  static gc::Vector3
  computeHalfedgeVolumeVariationVector(gcs::VertexPositionGeometry &vpg,
                                       gc::Halfedge &he);

  /**
   * @brief Compute vertex mean curvature vector using cotan
   */
  gcs::VertexData<gc::Vector3> computeVertexMeanCurvatureVectors();

  /**
   * @brief Compute halfedge mean curvature vector using cotan
   */
  static gc::Vector3
  computeHalfedgeMeanCurvatureVector(gcs::VertexPositionGeometry &vpg,
                                     gc::Halfedge &he);

  /**
   * @brief Compute vertex Gaussian curvature vector
   */
  gcs::VertexData<gc::Vector3> computeVertexGaussianCurvatureVectors();

  /**
   * @brief Compute halfedge Gaussian curvature vector
   */
  static gc::Vector3
  computeHalfedgeGaussianCurvatureVector(gcs::VertexPositionGeometry &vpg,
                                         gc::Halfedge &he);

  /**
   * @brief Compute vertex Schlafli-based Laplacian of mean curvature vector
   */
  gcs::VertexData<gc::Vector3>
  computeVertexSchlafliLaplacianMeanCurvatureVectors();

  /**
   * @brief Compute halfedge Schlafli vector
   */
  static std::tuple<gc::Vector3, gc::Vector3>
  computeHalfedgeSchlafliVector(gcs::VertexPositionGeometry &vpg,
                                gc::Halfedge &he);

  /**
   * @brief Helper functions to compute shape variation of corner angles
   */
  gc::Vector3 computeCornerAngleVariation(gcs::Corner c, gcs::Vertex v);

  /**
   * @brief Helper functions to compute shape variation of dihedral angles
   */
  gc::Vector3 computeDihedralAngleVariation(gcs::Halfedge he, gcs::Vertex v);

  /**
   * @brief Compute halfedge |\int grad phi|^2 variation vector
   */
  gc::Vector3 computeHalfedgeSquaredIntegratedDerivativeNormVariationVector(
      const gcs::VertexData<double> &quantities, const gcs::Halfedge &he);

  // ==========================================================
  // ================        Force.cpp       ==================
  // ==========================================================
  /**
   * @brief Compute and update all conservative forces, update
   * mechanicalForce(Vec) with conservativeForce(Vec)
   */
  void computeConservativeForcing();

  /**
   * @brief Compute and append all non-conservative forces, update
   * mechanicalForce(Vec) and mechErrorNorm
   */
  void addNonconservativeForcing(double timeStep);

  /**
   * @brief Compute geometric forces, including
   * - spontaneous curvature force
   * - deviatoric curvature force
   * - area difference force
   * - capillary (surface tension) force
   * - osmotic force
   * - line capillary (line tension) force
   * - adsorption (area expansion) force
   * - aggregation (area expansion) force
   * - entropy (area expansion) force
   */
  void computeGeometricForces();
  void computeGeometricForces(size_t i);
  void computeGeometricForces(gcs::Vertex &v);

  /**
   * @brief Compute regularization pressure component of the system
   */
  void computeSpringForces();

  /**
   * @brief Compute Self Avoidance force
   */
  void computeSelfAvoidanceForce();

  /**
   * @brief Compute external force component of the system
   */
  EigenVectorX3dr prescribeExternalForce();

  /**
   * @brief Compute DPD forces of the system
   */
  void computeDPDForces(double dt);

  /**
   * @brief Compute damping forces of the system
   */
  gc::VertexData<gc::Vector3> computeDampingForce();

  /**
   * @brief Compute chemical potential of the system, including
   * - spontaneous curvature potential
   * - adsorption potential
   * - aggregation potential
   * - entropy potential
   * - deviatoric curvature potential
   * - dirichlet potential
   * - interior penalty potential
   */
  void computeChemicalPotentials();

  /**
   * @brief Compute in plane flux form on edge
   */
  EigenVectorX1d computeInPlaneFluxForm(EigenVectorX1d &chemicalPotential);

  // ==========================================================
  // ================        energy.cpp      ==================
  // ==========================================================
  /**
   * @brief Compute spontaneous curvature energy
   */
  void computeSpontaneousCurvatureEnergy();

  /**
   * @brief Compute deviatoric curvature energy
   */
  void computeDeviatoricCurvatureEnergy();

  /**
   * @brief Compute area difference energy
   */
  void computeAreaDifferenceEnergy();

  /**
   * @brief Compute surface energy
   */
  void computeSurfaceEnergy();

  /**
   * @brief Compute pressure work
   */
  void computePressureEnergy();

  /**
   * @brief Compute adsorption energy
   */
  void computeAdsorptionEnergy();

  /**
   * @brief Compute aggregation energy
   */
  void computeAggregationEnergy();

  /**
   * @brief Compute entropy penalty
   */
  void computeEntropyEnergy();

  /**
   * @brief Compute protein interior penalty
   */
  void computeProteinInteriorPenalty();

  /**
   * @brief Compute Dirichlet energy
   */
  void computeDirichletEnergy();

  /**
   * @brief Compute self-avoidance energy
   */
  void computeSelfAvoidanceEnergy();

  /**
   * @brief Compute edge spring energy
   */
  void computeEdgeSpringEnergy();

  /**
   * @brief Compute face spring energy
   */
  void computeFaceSpringEnergy();

  /**
   * @brief Compute LCR spring energy
   */
  void computeLcrSpringEnergy();

  /**
   * @brief Compute external work
   */
  double computeExternalWork(double currentTime, double dt);

  /**
   * @brief Compute kinetic energy
   */
  double computeKineticEnergy();

  /**
   * @brief Compute potential energy
   */
  double computePotentialEnergy();

  /**
   * @brief compute total energy
   */
  double computeTotalEnergy();

  /**
   * @brief Intermediate function to integrate the power
   */
  double computeIntegratedPower(double dt);
  double computeIntegratedPower(double dt, EigenVectorX3dr &&velocity);

  /**
   * @brief Get tangential derivative of quantities on face
   */
  void computeFaceTangentialDerivative(gcs::VertexData<double> &quantities,
                                       gcs::FaceData<gc::Vector3> &gradient);

  /**
   * @brief helper function to compute LCR
   */
  double computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg,
                                 gcs::Halfedge &he) const;

  // ==========================================================
  // =============        regularization.cpp    ===============
  // ==========================================================
  /**
   * @brief Mesh mutation
   */
  void mutateMesh(size_t nMutation = 1);

  /**
   * @brief Apply vertex shift by moving the vertices chosen for integration to
   * the Barycenter of the it neighbors
   */
  void vertexShift();

  /**
   * @brief Edge flip if not Delaunay
   */
  bool edgeFlip();
  void flipEdge();

  /**
   * @brief Get regularization pressure component of the system
   */
  bool meshGrowth();
  bool growMesh();

  /**
   * @brief global smoothing after mutation of the mesh
   * @param initStep init guess of time step
   * @param target target reduce of force norm
   * @param maxIteration maximum number of iteration
   */
  Eigen::Matrix<bool, Eigen::Dynamic, 1>
  smoothenMesh(double initStep = 0.01, double target = 0.7,
               size_t maxIteration = 1000);
  /**
   * @brief pointwise smoothing after mutation of the mesh
   */
  void localSmoothing(const gcs::Vertex &v, std::size_t num = 10,
                      double stepSize = 0.01);
  void localSmoothing(const gcs::Halfedge &he, std::size_t num = 10,
                      double stepSize = 0.01);

  /**
   * @brief global update of quantities after mutation of the mesh
   */
  void globalUpdateAfterMutation();

  // ==========================================================
  // =============          misc.cpp            ===============
  // ==========================================================
  /**
   * @brief test conservative force computation by validating energy decrease
   * @return
   */
  bool testConservativeForcing(const double timeStep);

  /**
   * @brief backtrace energy increase from the system perturbation
   * @return
   */
  void backtraceEnergyGrowth(const double timeStep,
                             const Energy previousEnergy);

  /**
   * @brief infer the target surface area of the system
   */
  double inferTargetSurfaceArea();

  /**
   * @brief Check finiteness of forcing and energy
   * @return whether the system is finite
   */
  bool checkFiniteness();

  /**
   * @brief testing of random number generator pcg
   *
   */
  void check_pcg();

  /**
   * @brief Find "the" vertex
   */
  void findFloatCenter(double range = std::numeric_limits<double>::max());
  void findVertexCenter(double range = std::numeric_limits<double>::max());

  /**
   * @brief update cache of geodesicDistance
   */
  EigenVectorX1d computeGeodesicDistance();

  /**
   * @brief update cache of proteinDensity based on geodesicDistance and
   * parameters on profile type
   */
  void prescribeProteinDensityDistribution();

  /**
   * @brief prescribe mask based on geodesic disk
   */
  void prescribeGeodesicMasks();
};
} // namespace solver
} // namespace mem3dg
