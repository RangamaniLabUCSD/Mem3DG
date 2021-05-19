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

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"
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

struct DLL_PUBLIC MeshMutator {
  /// flip non-delaunay edge
  bool flipNonDelaunay = false;
  /// whether require flatness condition when flipping non-Delaunay edge
  bool flipNonDelaunayRequireFlat = false;

  /// split edge with large faces
  bool splitLarge = false;
  /// split long edge
  bool splitLong = false;
  /// split edge on high curvature domain
  bool splitCurved = false;
  /// split edge with sharp membrane property change
  bool splitSharp = false;
  /// split obtuse triangle
  bool splitFat = false;
  /// split poor aspected triangle that is still Delaunay
  bool splitSkinnyDelaunay = false;

  /// collapse skinny triangles
  bool collapseSkinny = false;
  /// collapse small triangles
  bool collapseSmall = false;
  /// whether require flatness condition when collapsing small edge
  bool collapseSmallNeedFlat = false;

  /// tolerance for curvature approximation
  double curvTol = 0.0012;
  /// target face area
  double targetFaceArea = 0.001;

  MeshMutator() {}
  ~MeshMutator() {}

  /**
   * @brief return condition of edge flip
   */
  bool ifFlip(const gcs::Edge e, const gcs::VertexPositionGeometry &vpg);

  /**
   * @brief return condition of edge split
   */
  bool ifSplit(const gc::Edge e, const gcs::VertexPositionGeometry &vpg);

  /**
   * @brief return condition of edge collapse
   */
  bool ifCollapse(const gc::Edge e, const gcs::VertexPositionGeometry &vpg);

  void maskAllNeighboring(gcs::VertexData<bool> &smoothingMask,
                          const gcs::Vertex v);

  void neighborAreaSum(const gcs::Edge e,
                       const gcs::VertexPositionGeometry &vpg, double &area,
                       size_t &num_neighbor);

  double
  computeCurvatureThresholdLength(const gcs::Edge e,
                                  const gcs::VertexPositionGeometry &vpg);
};

struct Forces {
  /// Cached mesh of interest
  gcs::ManifoldSurfaceMesh &mesh;
  /// Embedding and other geometric details
  gcs::VertexPositionGeometry &vpg;

  /// Cached bending force
  gcs::VertexData<double> bendingForce;
  /// Cached tension-induced capillary force
  gcs::VertexData<double> capillaryForce;
  /// Cached Surface tension
  double surfaceTension;
  /// Cached interfacial line tension force
  gcs::VertexData<double> lineCapillaryForce;
  gcs::EdgeData<double> lineTension;
  /// Cached externally-applied force
  gcs::VertexData<double> externalForce;
  /// Cached osmotic force
  gcs::VertexData<double> osmoticForce;
  double osmoticPressure;
  /// Cached three fundamentals
  gcs::VertexData<gc::Vector3> vectorForces;
  /// Cached bending force
  gcs::VertexData<gc::Vector3> bendingForceVec;
  /// Cached tension-induced capillary force
  gcs::VertexData<gc::Vector3> capillaryForceVec;
  /// Cached osmotic force
  gcs::VertexData<gc::Vector3> osmoticForceVec;
  /// Cached osmotic force
  gcs::VertexData<gc::Vector3> lineCapillaryForceVec;

  /// Cached local stretching forces (in-plane regularization)
  gcs::VertexData<gc::Vector3> regularizationForce;
  /// Cached damping forces
  gcs::VertexData<gc::Vector3> dampingForce;
  /// Cached stochastic forces
  gcs::VertexData<gc::Vector3> stochasticForce;

  /// Cached bending related chemical potential
  gcs::VertexData<double> bendingPotential;
  /// Cached adsorption related chemical potential
  gcs::VertexData<double> adsorptionPotential;
  /// Cached dirichlet energy related chemical potential
  gcs::VertexData<double> diffusionPotential;
  /// Cached chemical potential
  gcs::VertexData<double> chemicalPotential;

  /// force mask
  gcs::VertexData<gc::Vector3> forceMask;

  Forces(gcs::ManifoldSurfaceMesh &mesh_, gcs::VertexPositionGeometry &vpg_)
      : mesh(mesh_), vpg(vpg_), vectorForces(mesh, {0, 0, 0}),
        bendingForceVec(mesh, {0, 0, 0}), capillaryForceVec(mesh, {0, 0, 0}),
        osmoticForceVec(mesh, {0, 0, 0}),
        lineCapillaryForceVec(mesh, {0, 0, 0}), bendingForce(mesh, 0),
        capillaryForce(mesh, 0), surfaceTension(0), lineTension(mesh, 0),
        lineCapillaryForce(mesh, 0), externalForce(mesh, 0),
        osmoticForce(mesh, 0), osmoticPressure(0),
        regularizationForce(mesh, {0, 0, 0}), stochasticForce(mesh, {0, 0, 0}),
        dampingForce(mesh, {0, 0, 0}), bendingPotential(mesh, 0),
        adsorptionPotential(mesh, 0), diffusionPotential(mesh, 0),
        chemicalPotential(mesh, 0), forceMask(mesh, {1.0, 1.0, 1.0}) {}

  ~Forces() {}

  // ==========================================================
  // =============      Data interop helpers    ===============
  // ==========================================================

  /**
   * @brief Return raw buffer of a vertexData that contains gc::Vector3 or
   * scaler values
   */
  inline auto toMatrix(gcs::VertexData<gc::Vector3> &vector) {
    return gc::EigenMap<double, 3>(vector);
  }
  inline auto toMatrix(gcs::VertexData<gc::Vector3> &&vector) {
    return gc::EigenMap<double, 3>(vector);
  }

  inline auto toMatrix(gcs::VertexData<double> &vector) { return vector.raw(); }
  inline auto toMatrix(gcs::VertexData<double> &&vector) {
    return vector.raw();
  }

  /**
   * @brief Return constructed vertexData that from Eigen matrix or vector
   */
  inline gcs::VertexData<double> toVertexData(EigenVectorX1D &vector) {
    return gcs::VertexData<double>(mesh, vector);
  }
  inline gcs::VertexData<double> toVertexData(EigenVectorX1D &&vector) {
    return gcs::VertexData<double>(mesh, vector);
  }

  inline gcs::VertexData<gc::Vector3> toVertexData(EigenVectorX3D &vector) {
    gcs::VertexData<gc::Vector3> vertexData(mesh);
    gc::EigenMap<double, 3>(vertexData) = vector;
    return vertexData;
  }
  inline gcs::VertexData<gc::Vector3> toVertexData(EigenVectorX3D &&vector) {
    gcs::VertexData<gc::Vector3> vertexData(mesh);
    gc::EigenMap<double, 3>(vertexData) = vector;
    return vertexData;
  }

  /**
   * @brief Find the corresponding scalar vector (vertexData) by rowwise
   * (vertexwise) projecting matrix (vector vertexData) onto angle-weighted
   * normal vector
   */
  inline gcs::VertexData<double>
  ontoNormal(gcs::VertexData<gc::Vector3> &vector) {
    gcs::VertexData<double> vertexData(mesh);
    vertexData.raw() =
        rowwiseDotProduct(gc::EigenMap<double, 3>(vector),
                          gc::EigenMap<double, 3>(vpg.vertexNormals));
    return vertexData;
  }
  inline gcs::VertexData<double>
  ontoNormal(gcs::VertexData<gc::Vector3> &&vector) {
    gcs::VertexData<double> vertexData(mesh);
    vertexData.raw() =
        rowwiseDotProduct(gc::EigenMap<double, 3>(vector),
                          gc::EigenMap<double, 3>(vpg.vertexNormals));
    return vertexData;
  }
  inline EigenVectorX1D ontoNormal(EigenVectorX3D &vector) {
    return rowwiseDotProduct(vector,
                             gc::EigenMap<double, 3>(vpg.vertexNormals));
  }
  inline EigenVectorX1D ontoNormal(EigenVectorX3D &&vector) {
    return rowwiseDotProduct(vector,
                             gc::EigenMap<double, 3>(vpg.vertexNormals));
  }
  inline double ontoNormal(gc::Vector3 &vector, gc::Vertex &v) {
    return dot(vector, vpg.vertexNormals[v]);
  }
  inline double ontoNormal(gc::Vector3 &&vector, gc::Vertex &v) {
    return dot(vector, vpg.vertexNormals[v]);
  }

  /**
   * @brief Find the corresponding matrix (vector vertexData) by rowwise
   * (vertexwise) appending angle-weighted normal vector to the scalar vector
   * (vertexData)
   */
  inline gcs::VertexData<gc::Vector3>
  addNormal(gcs::VertexData<double> &vector) {
    gcs::VertexData<gc::Vector3> vertexData(mesh);
    gc::EigenMap<double, 3>(vertexData) = rowwiseScaling(
        vector.raw(), gc::EigenMap<double, 3>(vpg.vertexNormals));
    return vertexData;
  }
  inline gcs::VertexData<gc::Vector3>
  addNormal(gcs::VertexData<double> &&vector) {
    gcs::VertexData<gc::Vector3> vertexData(mesh);
    gc::EigenMap<double, 3>(vertexData) = rowwiseScaling(
        vector.raw(), gc::EigenMap<double, 3>(vpg.vertexNormals));
    return vertexData;
  }

  inline EigenVectorX3D addNormal(EigenVectorX1D &vector) {
    return rowwiseScaling(vector, gc::EigenMap<double, 3>(vpg.vertexNormals));
  }
  inline EigenVectorX3D addNormal(EigenVectorX1D &&vector) {
    return rowwiseScaling(vector, gc::EigenMap<double, 3>(vpg.vertexNormals));
  }

  inline gc::Vector3 addNormal(double &vector, gc::Vertex &v) {
    return vector * vpg.vertexNormals[v];
  }
  inline gc::Vector3 addNormal(double &&vector, gc::Vertex &v) {
    return vector * vpg.vertexNormals[v];
  }

  /**
   * @brief Project the vector onto tangent plane by removing the
   * angle-weighted normal component
   */
  inline EigenVectorX3D toTangent(gcs::VertexData<gc::Vector3> &vector) {
    return gc::EigenMap<double, 3>(vector) -
           rowwiseScaling(
               rowwiseDotProduct(gc::EigenMap<double, 3>(vector),
                                 gc::EigenMap<double, 3>(vpg.vertexNormals)),
               gc::EigenMap<double, 3>(vpg.vertexNormals));
  }
  inline EigenVectorX3D toTangent(gcs::VertexData<gc::Vector3> &&vector) {
    return gc::EigenMap<double, 3>(vector) -
           rowwiseScaling(
               rowwiseDotProduct(gc::EigenMap<double, 3>(vector),
                                 gc::EigenMap<double, 3>(vpg.vertexNormals)),
               gc::EigenMap<double, 3>(vpg.vertexNormals));
  }

  inline gc::Vector3 toTangent(gc::Vector3 &vector, gc::Vertex &v) {
    return vector -
           gc::dot(vector, vpg.vertexNormals[v]) * vpg.vertexNormals[v];
  }
  inline gc::Vector3 toTangent(gc::Vector3 &&vector, gc::Vertex &v) {
    return vector -
           gc::dot(vector, vpg.vertexNormals[v]) * vpg.vertexNormals[v];
  }

  /**
   * @brief Find the masked force
   */
  inline gcs::VertexData<gc::Vector3>
  mask(gcs::VertexData<gc::Vector3> &vector) {
    gcs::VertexData<gc::Vector3> vertexData(mesh);
    gc::EigenMap<double, 3>(vertexData).array() =
        gc::EigenMap<double, 3>(vector).array() *
        gc::EigenMap<double, 3>(forceMask).array();
    return vertexData;
  }
  inline gcs::VertexData<gc::Vector3>
  mask(gcs::VertexData<gc::Vector3> &&vector) {
    gcs::VertexData<gc::Vector3> vertexData(mesh);
    gc::EigenMap<double, 3>(vertexData).array() =
        gc::EigenMap<double, 3>(vector).array() *
        gc::EigenMap<double, 3>(forceMask).array();
    return vertexData;
  }
  inline EigenVectorX3D mask(EigenVectorX3D &vector) {
    return vector.array() * gc::EigenMap<double, 3>(forceMask).array();
  }
  inline EigenVectorX3D mask(EigenVectorX3D &&vector) {
    return vector.array() * gc::EigenMap<double, 3>(forceMask).array();
  }
  inline gc::Vector3 mask(gc::Vector3 &vector, gc::Vertex &v) {
    return gc::Vector3{vector.x * forceMask[v].x, vector.y * forceMask[v].y,
                       vector.z * forceMask[v].z};
  }
  inline gc::Vector3 mask(gc::Vector3 &&vector, gc::Vertex &v) {
    return gc::Vector3{vector.x * forceMask[v].x, vector.y * forceMask[v].y,
                       vector.z * forceMask[v].z};
  }
};

struct Parameters {
  /// Bending modulus
  double Kb;
  /// linear constant of bending modulus vs protein density
  double Kbc;
  /// Spontaneous curvature
  double H0;
  /// radius of heterogenous domain
  std::vector<double> r_heter;
  /// Global stretching modulus
  double Ksg;
  /// Area reservior
  double A_res;
  /// Vertex shifting constant
  double Kst;
  /// Local stretching modulus
  double Ksl;
  /// Edge spring constant
  double Kse;
  /// pressure-volume modulus
  double Kv;
  /// volume reservoir
  double V_res;
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
  /// interior point parameter for protein density
  double lambdaPhi = 1e-10;
  /// sharpness of tanh transition
  double sharpness = 20;
  /// type of relation between H0 and protein density
  std::string relation = "linear";
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

struct Options {
  /// Whether or not do vertex shift
  bool isVertexShift = false;
  /// Whether or not consider protein binding
  bool isProteinVariation = false;
  /// Whether adopt reduced volume parametrization
  bool isReducedVolume = false;
  /// Whether adopt constant osmotic pressure
  bool isConstantOsmoticPressure = false;
  /// Whether adopt constant surface tension
  bool isConstantSurfaceTension = false;
  /// Whether precribe hetergenous membrane by geodesic distance
  bool isHeterogeneous = false;
  /// Whether edge flip
  bool isEdgeFlip = false;
  /// Whether split edge
  bool isSplitEdge = false;
  /// Whether collapse edge
  bool isCollapseEdge = false;
  /// Whether need reference mesh
  bool isRefMesh = false;
  /// Whether floating "the" vertex
  bool isFloatVertex = false;
  /// Whether Laplacian mean curvature
  bool isLaplacianMeanCurvature = false;
  /// Boundary condition: roller, pin, fixed, none
  std::string boundaryConditionType = "none";
  /// Whether open boundary mesh
  bool isOpenMesh = false;
  /// Whether compute geodesic distance
  bool isComputeGeodesics = false;
};

class DLL_PUBLIC System {
public:
  /// Parameters
  Parameters P;
  /// Options
  Options O;
  /// Mesh mutator;
  MeshMutator meshMutator;

  /// Cached mesh of interest
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  /// Embedding and other geometric details
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  /// reference embedding geometry
  std::unique_ptr<gcs::VertexPositionGeometry> refVpg;
  /// Energy
  Energy E;
  /// Time
  double time;

  /// Forces of the system
  Forces F;

  /// Mean target area per face
  double meanTargetFaceArea;
  /// Target total face area
  double refSurfaceArea;
  /// Maximal volume
  double refVolume;
  /// Mean target area per face
  double meanTargetEdgeLength;
  /// Target edge cross length ratio
  gcs::EdgeData<double> targetLcrs;
  /// Reference edge Length of reference mesh
  gcs::EdgeData<double> refEdgeLengths;
  /// Reference face Area of reference mesh
  gcs::FaceData<double> refFaceAreas;

  /// Cached geodesic distance
  gcs::VertexData<double> geodesicDistanceFromPtInd;
  /// V-E distribution matrix
  Eigen::SparseMatrix<double> D;

  /// L1 mechanical error norm
  double L1ErrorNorm;
  /// L1 chemical error norm
  double L1ChemErrorNorm;
  /// surface area
  double surfaceArea;
  /// Volume
  double volume;
  /// Cached vertex positions from the previous step
  gcs::VertexData<gc::Vector3> pastPositions;
  /// Cached protein surface density
  gcs::VertexData<double> proteinDensity;
  /// Spontaneous curvature gradient of the mesh
  gcs::FaceData<gc::Vector3> proteinDensityGradient;
  /// Cached vertex velocity
  gcs::VertexData<gc::Vector3> vel;
  /// Spontaneous curvature of the mesh
  gcs::VertexData<double> H0;
  /// Bending rigidity of the membrane
  gcs::VertexData<double> Kb;
  /// Random number engine
  pcg32 rng;
  std::normal_distribution<double> normal_dist;
  /// "the vertex"
  gcs::SurfacePoint thePoint;
  // "the vertex" tracker
  gcs::VertexData<bool> thePointTracker;
  // is Smooth
  bool isSmooth;
  // is
  gcs::VertexData<bool> smoothingMask;

  // ==========================================================
  // =============        Constructors           ==============
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
         Parameters &p, Options &o)
      : System(readMeshes(topologyMatrix, vertexMatrix, refVertexMatrix, nSub),
               p, o) {
    // Check confliciting parameters and options
    checkParametersAndOptions();

    // Initialize reference values
    initConstants();

    // Process the mesh by regularization and mutation
    processMesh();

    /// compute nonconstant values during simulation
    updateVertexPositions();
  };

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
  System(std::string inputMesh, std::string refMesh, size_t nSub,
         bool isContinue, Parameters &p, Options &o)
      : System(readMeshes(inputMesh, refMesh, nSub), p, o) {
    // Check confliciting parameters and options
    checkParametersAndOptions();

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
    processMesh();

    /// compute nonconstant values during simulation
    updateVertexPositions();
  };

#ifdef MEM3DG_WITH_NETCDF
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
         Parameters &p, Options &o)
      : System(readTrajFile(trajFile, startingFrame, nSub), p, o) {
    // Check confliciting parameters and options
    checkParametersAndOptions();

    // Initialize reference values
    initConstants();

    // Map continuation variables
    if (isContinue) {
      mapContinuationVariables(trajFile, startingFrame);
    }

    // Process the mesh by regularization and mutation
    processMesh();

    /// compute nonconstant values during simulation
    updateVertexPositions();
  };
#endif

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
         Parameters &p, Options &o)
      : System(std::move(std::get<0>(meshVpgTuple)),
               std::move(std::get<1>(meshVpgTuple)),
               std::move(std::get<2>(meshVpgTuple)), p, o){};

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
         Options &o)
      : mesh(std::move(ptrmesh_)), vpg(std::move(ptrvpg_)),
        refVpg(std::move(ptrrefVpg_)), P(p), O(o), time(0),
        E({0, 0, 0, 0, 0, 0, 0, 0, 0}), F(*mesh, *vpg),
        proteinDensity(*mesh, 0.5), targetLcrs(*mesh), refEdgeLengths(*mesh),
        refFaceAreas(*mesh), D(), geodesicDistanceFromPtInd(*mesh, 0),
        thePointTracker(*mesh, false), pastPositions(*mesh, {0, 0, 0}),
        vel(*mesh, {0, 0, 0}), H0(*mesh),
        proteinDensityGradient(*mesh, {0, 0, 0}), Kb(*mesh), isSmooth(true),
        smoothingMask(*mesh, false) {

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

  /**
   * @brief Map the continuation variables
   *
   */
  void mapContinuationVariables(std::string plyFile);

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
  void checkParametersAndOptions();

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
   * @brief Mesh processing by regularization or mutation
   */
  void processMesh();

  /**
   * @brief Update the vertex position and recompute cached values
   * (all quantities that characterizes the current energy state)
   * Careful: 1. when using eigenMap: memory address may change after update!!
   * Careful: 2. choosing to update geodesics and spatial properties may lead to
   * failing in backtrack!!
   */
  void updateVertexPositions(bool isUpdateGeodesics = false);

  // ==========================================================
  // ================        Pressure        ==================
  // ==========================================================
  /**
   * @brief Compute bending force component of the system
   */
  EigenVectorX1D computeBendingForce();

  /**
   * @brief Compute capillary force component of the system
   */
  EigenVectorX1D computeCapillaryForce();

  /**
   * @brief Compute osmotic force component of the system
   */
  EigenVectorX1D computeOsmoticForce();

  /**
   * @brief Compute fundamental three forces at the same time
   */
  EigenVectorX3D computeVectorForces();

  /**
   * @brief Compute chemical potential of the system
   */
  EigenVectorX1D computeChemicalPotential();

  /**
   * @brief Compute line tension force component of the system
   */
  EigenVectorX1D computeLineCapillaryForce();

  /**
   * @brief Compute external force component of the system
   */
  EigenVectorX1D computeExternalForce();

  /**
   * @brief Compute all forces of the system
   */
  void computePhysicalForces();

  /**
   * @brief Compute DPD forces of the system
   */
  std::tuple<EigenVectorX3D, EigenVectorX3D> computeDPDForces(double dt);

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
   * @brief Compute adsorption energy
   */
  void computeAdsorptionEnergy();

  /**
   * @brief Compute Dirichlet energy
   */
  void computeDirichletEnergy();

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
  bool edgeFlip();

  /**
   * @brief Get regularization pressure component of the system
   */
  bool growMesh();

  // ==========================================================
  // =============          Helpers             ===============
  // ==========================================================
  /**
   * @brief Get length cross ratio of the mesh
   */

  void computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg,
                               gcs::EdgeData<double> &targetLcrs);
  double computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg,
                                 gcs::Edge &e) const;
  double computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg,
                                 gcs::Edge &&e) const;

  /**
   * @brief Get gradient of quantities on face
   */
  void computeGradient(gcs::VertexData<double> &quantities,
                       gcs::FaceData<gc::Vector3> &gradient);

  /**
   * @brief Find "the" vertex
   */
  void findThePoint(gcs::VertexPositionGeometry &vpg,
                    gcs::VertexData<double> &geodesicDistance,
                    double range = 1e10);

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
   * @brief global smoothing after mutation of the mesh
   */
  void globalSmoothing(gcs::VertexData<bool> &smoothingMask, double tol = 1e-6,
                       double stepSize = 1);
};
} // namespace mem3dg
