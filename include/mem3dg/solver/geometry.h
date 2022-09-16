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

class DLL_PUBLIC Geometry {
public:
  /// Cached mesh of interest
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  /// Embedding and other geometric details
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  /// surface area
  double surfaceArea;
  /// Volume
  double volume;
  /// reservoir area
  double reservoirArea;
  /// reservoir volume
  double reservoirVolume;

  /// Cached geodesic distance
  gcs::VertexData<double> geodesicDistance;
  /// if has boundary
  bool isOpenMesh;
  /// defined notable vertex of the mesh
  gcs::VertexData<bool> notableVertex;

  // ==========================================================
  // =============        Constructors           ==============
  // ==========================================================

  // =======================================
  // =======       Matrices         ========
  // =======================================
  Geometry(EigenVectorX3sr &topologyMatrix, EigenVectorX3dr &vertexMatrix)
      : Geometry(readMatrices(topologyMatrix, vertexMatrix)){};

  // =======================================
  // =======       Mesh Files       ========
  // =======================================
  Geometry(std::string inputMesh) : Geometry(readMeshFile(inputMesh)){};

  // =======================================
  // =======       NetCDF Files     ========
  // =======================================

  // =======================================
  // =======       Tuple            ========
  // =======================================
  Geometry(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                      std::unique_ptr<gcs::VertexPositionGeometry>>
               meshVpgTuple)
      : Geometry(std::move(std::get<0>(meshVpgTuple)),
                 std::move(std::get<1>(meshVpgTuple))){};

  // =======================================
  // =======    Geometry Central    ========
  // =======================================

  Geometry(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
           std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_)
      : mesh(std::move(ptrmesh_)), vpg(std::move(ptrvpg_)) {
    geodesicDistance = gcs::VertexData<double>(*mesh, 0);
    notableVertex = gc::VertexData<bool>(*mesh, false);

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
  ~Geometry() {
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

  /**
   * @brief Construct a tuple of unique_ptrs from mesh and refMesh path
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMeshFile(std::string inputMesh);

  /**
   * @brief Save RichData to .ply file
   *
   */
  void saveGeometry(std::string PathToSave);

#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Construct a tuple of unique_ptrs from netcdf path
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readTrajFile(std::string trajFile, int startingFrame);
#endif

  // ==========================================================
  // ================     init.cpp           ==================
  // ==========================================================

  /**
   * @brief Initialize system
   *
   */
  void initialize(bool ifMute = false);

  /**
   * @brief Update the vertex position and recompute cached values
   * (all quantities that characterizes the current energy state)
   * Careful: 1. when using eigenMap: memory address may change after update!!
   */
  void updateConfigurations();

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
  computeVertexSchlafliLaplacianMeanCurvatureVectors(
      gcs::VertexData<double> &spontaneousCurvature);

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
  // =============          misc.cpp            ===============
  // ==========================================================
  /**
   * @brief infer the target surface area of the system
   */
  double inferTargetSurfaceArea();

  /**
   * @brief update cache of geodesicDistance
   */
  EigenVectorX1d computeGeodesicDistance();

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
};
} // namespace solver
} // namespace mem3dg
