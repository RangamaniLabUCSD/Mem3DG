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

#include <functional>
#include <math.h>
#include <vector>

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include "mem3dg/constants.h"
#include "mem3dg/macros.h"
#include "mem3dg/mesh_io.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/forces.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/solver/parameters.h"
#include "mem3dg/type_utilities.h"

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
  /// bending energy of the membrane
  double bendingEnergy = 0;
  /// deviatoric energy of the membrane
  double deviatoricEnergy = 0;
  /// stretching energy of the membrane
  double surfaceEnergy = 0;
  /// work of pressure within membrane
  double pressureEnergy = 0;
  /// adsorption energy of the membrane protein
  double adsorptionEnergy = 0;
  /// aggregation energy of the membrane protein
  double aggregationEnergy = 0;
  /// line tension energy of interface
  double dirichletEnergy = 0;
  /// work of external force
  double externalWork = 0;
  /// protein interior penalty energy
  double proteinInteriorPenalty = 0;
  /// membrane self-avoidance penalty energy
  double selfAvoidancePenalty = 0;
};

class DLL_PUBLIC System {
protected:
  /// Cached geodesic distance
  gcs::VertexData<double> geodesicDistanceFromPtInd;
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
  /// Cached vertex protein velocity
  gcs::VertexData<double> proteinVelocity;
  /// Spontaneous curvature of the mesh
  gcs::VertexData<double> H0;
  /// Bending rigidity of the membrane
  gcs::VertexData<double> Kb;
  /// deviatoric rigidity of the membrane
  gcs::VertexData<double> Kd;

  /// is Smooth
  bool isSmooth;
  /// if being mutated
  gcs::VertexData<bool> mutationMarker;
  /// if has boundary
  bool isOpenMesh;
  /// "the vertex"
  gcs::SurfacePoint thePoint;
  gcs::VertexData<bool> thePointTracker;
  /// projected time of collision
  double projectedCollideTime;

  // ==========================================================
  // =============        Constructors           ==============
  // ==========================================================
  /**
   * @brief Construct a new (geometry) System by reading topology and vertex
   * matrices
   *
   * @param topologyMatrix,  topology matrix, F x 3
   * @param vertexMatrix,    input Mesh coordinate matrix, V x 3
   * @param nSub          Number of subdivision
   */
  System(Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &topologyMatrix,
         Eigen::Matrix<double, Eigen::Dynamic, 3> &vertexMatrix,
         std::size_t nSub)
      : System(readMeshes(topologyMatrix, vertexMatrix, nSub)) {

    // Initialize reference values
    initConstants();

    // compute nonconstant values during simulation
    updateConfigurations();
  };

  /**
   * @brief Construct a new System by reading topology and vertex matrices
   *
   * @param topologyMatrix,  topology matrix, F x 3
   * @param vertexMatrix,    input Mesh coordinate matrix, V x 3
   * @param p             Parameter of simulation
   * @param nSub          Number of subdivision
   */
  System(Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &topologyMatrix,
         Eigen::Matrix<double, Eigen::Dynamic, 3> &vertexMatrix, Parameters &p,
         std::size_t nSub)
      : System(readMeshes(topologyMatrix, vertexMatrix, nSub), p) {
    // Check incompatible configuration
    checkConfiguration();

    // Initialize reference values
    initConstants();

    // compute nonconstant values during simulation
    updateConfigurations();

    // Smoothen the mesh using bending force
    smoothenMesh(0.01, 0.1, 1000);
  };

  /**
   * @brief Construct a new System by reading topology and vertex matrices
   *
   * @param topologyMatrix,  topology matrix, F x 3
   * @param vertexMatrix,    input Mesh coordinate matrix, V x 3
   * @param p             Parameter of simulation
   * @param mp         Setting for mesh processing
   * @param nSub          Number of subdivision
   * @param nMutation     Number of mutation
   */
  System(Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &topologyMatrix,
         Eigen::Matrix<double, Eigen::Dynamic, 3> &vertexMatrix, Parameters &p,
         MeshProcessor &mp, std::size_t nSub, std::size_t nMutation)
      : System(readMeshes(topologyMatrix, vertexMatrix, nSub), p, mp) {
    // Check incompatible configuration
    checkConfiguration();

    // Initialize reference values
    initConstants();

    // Process the mesh by regularization and mutation
    mutateMesh(nMutation);

    // compute nonconstant values during simulation
    updateConfigurations();

    // Smoothen the mesh using bending force
    smoothenMesh(0.01, 0.1, 1000);
  };

  /**
   * @brief Construct a new (geometry) System by reading mesh file path
   *
   * @param inputMesh     Input Mesh
   * @param nSub          Number of subdivision
   */
  System(std::string inputMesh, std::size_t nSub)
      : System(readMeshes(inputMesh, nSub)) {
    // Initialize reference values
    initConstants();

    // compute nonconstant values during simulation
    updateConfigurations();
  };

  /**
   * @brief Construct a new System by reading mesh file path
   *
   * @param inputMesh     Input Mesh
   * @param p             Parameter of simulation
   * @param nSub          Number of subdivision
   * @param isContinue    Wether continue simulation
   */
  System(std::string inputMesh, Parameters &p, std::size_t nSub,
         bool isContinue)
      : System(readMeshes(inputMesh, nSub), p) {

    // Check incompatible configuration
    checkConfiguration();

    // Initialize reference values
    initConstants();

    // Map continuation variables
    if (isContinue) {
      std::cout << "\nWARNING: isContinue is on and make sure mesh file "
                   "supports richData!"
                << std::endl;
      mapContinuationVariables(inputMesh);
    }

    // compute nonconstant values during simulation
    updateConfigurations();

    // Smoothen the mesh using bending force
    smoothenMesh(0.01, 0.1, 1000);
  };

  /**
   * @brief Construct a new System by reading mesh file path
   *
   * @param inputMesh     Input Mesh
   * @param p             Parameter of simulation
   * @param mp         Setting for mesh processing
   * @param nSub          Number of subdivision
   * @param nMutation     Number of mutation
   * @param isContinue    Wether continue simulation
   */
  System(std::string inputMesh, Parameters &p, MeshProcessor &mp,
         std::size_t nSub, std::size_t nMutation, bool isContinue)
      : System(readMeshes(inputMesh, nSub), p, mp) {

    // Check incompatible configuration
    checkConfiguration();

    // Initialize reference values
    initConstants();

    // Map continuation variables
    if (isContinue) {
      std::cout << "\nWARNING: isContinue is on and make sure mesh file "
                   "supports richData!"
                << std::endl;
      mapContinuationVariables(inputMesh);
    }

    // Process the mesh by regularization and mutation
    mutateMesh(nMutation);

    // compute nonconstant values during simulation
    updateConfigurations();

    // Smoothen the mesh using bending force
    smoothenMesh(0.01, 0.1, 1000);
  };

#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Construct a new System object by reading netcdf trajectory file path
   *
   * @param trajFile      Netcdf trajectory file
   * @param startingFrame Starting frame for the input mesh
   * @param nSub          Number of subdivision
   */
  System(std::string trajFile, int startingFrame, std::size_t nSub)
      : System(readTrajFile(trajFile, startingFrame, nSub)) {

    // Initialize reference values
    initConstants();

    // compute nonconstant values during simulation
    updateConfigurations();
  };

  /**
   * @brief Construct a new System object by reading netcdf trajectory file path
   *
   * @param trajFile      Netcdf trajectory file
   * @param startingFrame Starting frame for the input mesh
   * @param p             Parameter of simulation
   * @param nSub          Number of subdivision
   * @param isContinue    Wether continue simulation
   */
  System(std::string trajFile, int startingFrame, Parameters &p,
         std::size_t nSub, bool isContinue)
      : System(readTrajFile(trajFile, startingFrame, nSub), p) {

    // Check incompatible configuration
    checkConfiguration();

    // Initialize reference values
    initConstants();

    // Map continuation variables
    if (isContinue) {
      mapContinuationVariables(trajFile, startingFrame);
    }

    // compute nonconstant values during simulation
    updateConfigurations();

    // Smoothen the mesh using bending force
    smoothenMesh(0.01, 0.1, 1000);
  };

  /**
   * @brief Construct a new System object by reading netcdf trajectory file path
   *
   * @param trajFile      Netcdf trajectory file
   * @param startingFrame Starting frame for the input mesh
   * @param p             Parameter of simulation
   * @param mp         Setting for mesh processing
   * @param nSub          Number of subdivision
   * @param nMutation     Number of mutation
   * @param isContinue    Wether continue simulation
   */
  System(std::string trajFile, int startingFrame, Parameters &p,
         MeshProcessor &mp, std::size_t nSub, std::size_t nMutation,
         bool isContinue)
      : System(readTrajFile(trajFile, startingFrame, nSub), p, mp) {

    // Check incompatible configuration
    checkConfiguration();

    // Initialize reference values
    initConstants();

    // Map continuation variables
    if (isContinue) {
      mapContinuationVariables(trajFile, startingFrame);
    }

    // Process the mesh by regularization and mutation
    mutateMesh(nMutation);

    // compute nonconstant values during simulation
    updateConfigurations();

    // Smoothen the mesh using bending force
    smoothenMesh(0.01, 0.1, 1000);
  };
#endif

private:
  /**
   * @brief Construct a new System object by reading tuple of unique_ptrs
   *
   * @param tuple        Mesh connectivity, Embedding and geometry
   * information, Mesh rich data
   * @param p             Parameter of simulation
   * @param mp         Setting for mesh processing
   */
  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>>
             meshVpgTuple,
         Parameters &p, MeshProcessor &mp)
      : System(std::move(std::get<0>(meshVpgTuple)),
               std::move(std::get<1>(meshVpgTuple)), p, mp){};

  /**
   * @brief Construct a new System object by reading tuple of unique_ptrs
   *
   * @param tuple        Mesh connectivity, Embedding and geometry
   * information, Mesh rich data
   * @param p             Parameter of simulation
   */
  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>>
             meshVpgTuple,
         Parameters &p)
      : System(std::move(std::get<0>(meshVpgTuple)),
               std::move(std::get<1>(meshVpgTuple)), p){};

  /**
   * @brief Construct a new System object by reading tuple of unique_ptrs
   *
   * @param tuple        Mesh connectivity, Embedding and geometry
   * information, Mesh rich data
   */
  System(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                    std::unique_ptr<gcs::VertexPositionGeometry>>
             meshVpgTuple)
      : System(std::move(std::get<0>(meshVpgTuple)),
               std::move(std::get<1>(meshVpgTuple))){};

  /**
   * @brief Construct a new System object by reading unique_ptrs to mesh and
   * geometry objects
   * @param ptrmesh_         Mesh connectivity
   * @param ptrvpg_          Embedding and geometry information
   * @param p             Parameter of simulation
   * @param mp         Setting for mesh processing
   */
  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_, Parameters &p,
         MeshProcessor &mp)
      : System(std::move(ptrmesh_), std::move(ptrvpg_)) {
    parameters = p;
    meshProcessor = mp;
  }

  /**
   * @brief Construct a new System object by reading unique_ptrs to mesh and
   * geometry objects
   * @param ptrmesh_         Mesh connectivity
   * @param ptrvpg_          Embedding and geometry information
   * @param p             Parameter of simulation
   */
  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_, Parameters &p)
      : System(std::move(ptrmesh_), std::move(ptrvpg_)) {
    parameters = p;
  }

  /**
   * @brief Construct a new System object by reading unique_ptrs to mesh and
   * geometry objects
   * @param ptrmesh_         Mesh connectivity
   * @param ptrvpg_          Embedding and geometry information
   */
  System(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
         std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_)
      : mesh(std::move(ptrmesh_)), vpg(std::move(ptrvpg_)),
        forces(*mesh, *vpg) {

    time = 0;
    energy = Energy({time, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    proteinDensity = gc::VertexData<double>(*mesh, 0);
    proteinDensityGradient = gcs::FaceData<gc::Vector3>(*mesh, {0, 0, 0});
    velocity = gcs::VertexData<gc::Vector3>(*mesh, {0, 0, 0});
    proteinVelocity = gcs::VertexData<double>(*mesh, 0);
    H0 = gcs::VertexData<double>(*mesh);
    Kb = gcs::VertexData<double>(*mesh);
    Kd = gcs::VertexData<double>(*mesh);

    geodesicDistanceFromPtInd = gcs::VertexData<double>(*mesh, 0);

    isSmooth = true;
    mutationMarker = gc::VertexData<bool>(*mesh, false);
    thePointTracker = gc::VertexData<bool>(*mesh, false);

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
  // ================          I/O           ==================
  // ==========================================================

  /**
   * @brief Construct a tuple of unique_ptrs from topology matrix and vertex
   * position matrix
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMeshes(Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &faceVertexMatrix,
             Eigen::Matrix<double, Eigen::Dynamic, 3> &vertexPositionMatrix,
             std::size_t nSub);

  /**
   * @brief Construct a tuple of unique_ptrs from mesh and refMesh path
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMeshes(std::string inputMesh, std::size_t nSub);

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
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readTrajFile(std::string trajFile, int startingFrame, std::size_t nSub);
  /**
   * @brief Map the continuation variables
   *
   */
  void mapContinuationVariables(std::string trajFile, int startingFrame);
#endif

  // ==========================================================
  // ================     Initialization     ==================
  // ==========================================================
  /**
   * @brief Check all conflicting parameters and options
   *
   */
  void checkConfiguration();

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
   * @brief Mesh mutation
   */
  void mutateMesh(size_t nRepetition = 1);

  /**
   * @brief Update the vertex position and recompute cached values
   * (all quantities that characterizes the current energy state)
   * Careful: 1. when using eigenMap: memory address may change after update!!
   * Careful: 2. choosing to update geodesics and spatial properties may lead to
   * failing in backtrack!!
   */
  void updateConfigurations(bool isUpdateGeodesics = false);

  // ==========================================================
  // ================   Variational vectors  ==================
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
  gcs::VertexData<gc::Vector3> computeVertexVolumeVariationVector();

  /**
   * @brief Compute halfedge volume variation vector
   */
  static gc::Vector3
  computeHalfedgeVolumeVariationVector(gcs::VertexPositionGeometry &vpg,
                                       gc::Halfedge &he);

  /**
   * @brief Compute vertex mean curvature vector using cotan
   */
  gcs::VertexData<gc::Vector3> computeVertexMeanCurvatureVector();

  /**
   * @brief Compute halfedge mean curvature vector using cotan
   */
  static gc::Vector3
  computeHalfedgeMeanCurvatureVector(gcs::VertexPositionGeometry &vpg,
                                     gc::Halfedge &he);

  /**
   * @brief Compute vertex Gaussian curvature vector
   */
  gcs::VertexData<gc::Vector3> computeVertexGaussianCurvatureVector();

  /**
   * @brief Compute halfedge Gaussian curvature vector
   */
  static gc::Vector3
  computeHalfedgeGaussianCurvatureVector(gcs::VertexPositionGeometry &vpg,
                                         gc::Halfedge &he);

  /**
   * @brief Compute vertex Schlafli vector
   */
  gcs::VertexData<gc::Vector3> computeVertexSchlafliVector();

  /**
   * @brief Compute halfedge Schlafli vector
   */
  static std::tuple<gc::Vector3, gc::Vector3>
  computeHalfedgeSchlafliVector(gcs::VertexPositionGeometry &vpg,
                                gc::Halfedge &he);

  /**
   * @brief Helper functions to compute geometric derivatives
   */
  gc::Vector3 cornerAngleGradient(gcs::Corner c, gcs::Vertex v);
  gc::Vector3 dihedralAngleGradient(gcs::Halfedge he, gcs::Vertex v);

  // ==========================================================
  // ================        Pressure        ==================
  // ==========================================================
  /**
   * @brief Compute all forcing of the system, include DPD if given time step
   */
  void computePhysicalForcing();
  void computePhysicalForcing(double timeStep);

  /**
   * @brief Compute chemical potential of the system
   */
  void computeChemicalPotentials();

  /**
   * @brief Compute Self Avoidance force
   */
  void computeSelfAvoidanceForce();

  /**
   * @brief Compute mechanical forces
   */
  void computeMechanicalForces();
  void computeMechanicalForces(size_t i);
  void computeMechanicalForces(gcs::Vertex &v);

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
  // ==========================================================
  // ================        Energy          ==================
  // ==========================================================
  /**
   * @brief Compute bending energy
   */
  void computeBendingEnergy();

  /**
   * @brief Compute deviatoric energy
   */
  void computeDeviatoricEnergy();

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
   * @brief Compute the L1 norm of the pressure
   */
  double computeNorm(
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &force) const;
  double computeNorm(
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &&force) const;
  /**
   * @brief Intermediate function to integrate the power
   */
  double computeIntegratedPower(double dt);
  double computeIntegratedPower(double dt, EigenVectorX3dr &&velocity);

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
  bool edgeFlip();

  /**
   * @brief Get regularization pressure component of the system
   */
  bool growMesh();

  // ==========================================================
  // =============          Helpers             ===============
  // ==========================================================

  /**
   * @brief Get gradient of quantities on face
   */
  void computeGradient(gcs::VertexData<double> &quantities,
                       gcs::FaceData<gc::Vector3> &gradient);

  /**
   * @brief Get gradient of quantities on face
   */
  gc::Vector3
  computeGradientNorm2Gradient(const gcs::Halfedge &he,
                               const gcs::VertexData<double> &quantities);

  /**
   * @brief Find "the" vertex
   */
  void findThePoint(gcs::VertexPositionGeometry &vpg,
                    gcs::VertexData<double> &geodesicDistance,
                    double range = 1e10);

  /**
   * @brief global smoothing after mutation of the mesh
   * @param initStep init guess of time step
   * @param target target reduce of force norm
   * @param maxIteration maximum number of iteration
   */
  Eigen::Matrix<bool, Eigen::Dynamic, 1>
  smoothenMesh(double initStep, double target = 0.7,
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

  /**
   * @brief infer the target surface area of the system
   */
  double inferTargetSurfaceArea();
};
} // namespace solver
} // namespace mem3dg
