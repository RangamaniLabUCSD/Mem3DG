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

#include <geometrycentral/surface/halfedge_element_types.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/heat_method_distance.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector2.h>
#include <geometrycentral/utilities/vector3.h>

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

#ifdef MEM3DG_WITH_NETCDF
#include "mem3dg/solver/mutable_trajfile.h"
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
  /// Embedding and other geometric details of reference mesh
  std::unique_ptr<gcs::VertexPositionGeometry> refVpg;
  /// reference length cross ratio
  gcs::HalfedgeData<double> refLcrs;
  /// surface area
  double surfaceArea;
  /// Volume
  double volume;

  /// Cached geodesic distance
  gcs::VertexData<double> geodesicDistance;
  /// defined notable vertex of the mesh
  gcs::VertexData<bool> notableVertex;

  // =======================================
  // =======       Matrices         ========
  // =======================================
  /**
   * @brief Construct Geometry
   * @param faceMatrix mesh face matrix
   * @param vertexMatrix mesh vertex matrix
   * @param notableVertex notable vertex data
   * @return geometry instance
   *
   */
  Geometry(EigenVectorX3sr &faceMatrix, EigenVectorX3dr &vertexMatrix,
           Eigen::Matrix<bool, Eigen::Dynamic, 1> &notableVertex_)
      : Geometry(readMatrices(faceMatrix, vertexMatrix), notableVertex_){};

  /**
   * @brief Construct Geometry
   * @param faceMatrix mesh face matrix
   * @param vertexMatrix mesh vertex matrix
   * @return geometry instance
   *
   */
  Geometry(EigenVectorX3sr &faceMatrix, EigenVectorX3dr &vertexMatrix)
      : Geometry(readMatrices(faceMatrix, vertexMatrix)){};

  // =======================================
  // =======       Mesh Files       ========
  // =======================================
  /**
   * @brief Construct Geometry
   * @param inputMesh file for initial mesh
   * @param notableVertex_ notable vertex data
   * @return geometry instance
   */
  Geometry(std::string inputMesh,
           Eigen::Matrix<bool, Eigen::Dynamic, 1> &notableVertex_)
      : Geometry(readMeshFile(inputMesh), notableVertex_){};

  /**
   * @brief Construct Geometry
   * @param inputMesh file for initial mesh
   * @return geometry instance
   */
  Geometry(std::string inputMesh) : Geometry(readMeshFile(inputMesh)){};

#ifdef MEM3DG_WITH_NETCDF
  // =======================================
  // =======       NetCDF Files     ========
  // =======================================
  /**
   * @brief Construct Geometry
   * @param trajFile NetCDF trajectory file
   * @param startingFrame starting frame
   * @return geometry instance
   */
  Geometry(std::string trajFile, int startingFrame)
      : Geometry(readTrajFile(trajFile, startingFrame)){};
#endif

private:
  // =======================================
  // =======       Tuple            ========
  // =======================================
  Geometry(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                      std::unique_ptr<gcs::VertexPositionGeometry>>
               tuple)
      : Geometry(std::move(std::get<0>(tuple)),
                 std::move(std::get<1>(tuple))){};

  Geometry(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                      std::unique_ptr<gcs::VertexPositionGeometry>>
               tuple,
           Eigen::Matrix<bool, Eigen::Dynamic, 1> &notableVertex_)
      : Geometry(std::move(std::get<0>(tuple)), std::move(std::get<1>(tuple)),
                 notableVertex_){};

  Geometry(std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                      std::unique_ptr<gcs::VertexPositionGeometry>,
                      Eigen::Matrix<bool, Eigen::Dynamic, 1>>
               tuple)
      : Geometry(std::move(std::get<0>(tuple)), std::move(std::get<1>(tuple)),
                 std::get<2>(tuple)){};

public:
  // =======================================
  // =======    Geometry Central    ========
  // =======================================
  /**
   * @brief Construct Geometry
   * @param ptrmesh_ manifold surface mesh
   * @param ptrvpg_ vertex position geometry
   * @param notableVertex_ notable vertex data
   * @return geometry instance
   */
  Geometry(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
           std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_,
           Eigen::Matrix<bool, Eigen::Dynamic, 1> &notableVertex_)
      : Geometry(std::move(ptrmesh_), std::move(ptrvpg_)) {
    notableVertex.raw() = notableVertex_;
    // auto &nv = notableVertex.raw();
    // if (std::count(nv.begin(), nv.end(), true) == 0)
    //   mem3dg_runtime_error("notableVertex cannot be all false!");
    computeGeodesicDistance();
  }

  /**
   * @brief Construct Geometry
   * @param ptrmesh_ manifold surface mesh
   * @param ptrvpg_ vertex position geometry
   * @return geometry instance
   */
  Geometry(std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrmesh_,
           std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg_)
      : mesh(std::move(ptrmesh_)), vpg(std::move(ptrvpg_)),
        geodesicDistance{*mesh, 0}, notableVertex{*mesh, false},
        refLcrs{*mesh, 0} {
    updateReferenceMeshFromVPG();

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

    // notableVertex[0] = true;
    volume = getMeshVolume(*mesh, *vpg, true);
    surfaceArea = vpg->faceAreas.raw().sum();
    updateReferenceConfigurations();
    if (!mesh->hasBoundary() && mesh->genus() != 0) {
      mem3dg_runtime_error(
          "Do not support closed mesh with nonzero number of genus!")
    }
  }

  /**
   * @brief Destroy the Geometry
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
   * @brief Construct a tuple of unique_ptrs from face matrix and vertex
   * position matrix
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMatrices(EigenVectorX3sr &faceVertexMatrix,
               EigenVectorX3dr &vertexPositionMatrix);
  // std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
  //            std::unique_ptr<gcs::VertexPositionGeometry>,
  //            std::unique_ptr<gcs::VertexPositionGeometry>>
  // readMatrices(EigenVectorX3sr &faceVertexMatrix,
  //              EigenVectorX3dr &vertexPositionMatrix,
  //              EigenVectorX3dr &refVertexPositionMatrix);
  /**
   * @brief Construct a tuple of unique_ptrs from mesh and refMesh path
   *
   */
  std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
             std::unique_ptr<gcs::VertexPositionGeometry>>
  readMeshFile(std::string inputMesh);

  // std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
  //            std::unique_ptr<gcs::VertexPositionGeometry>,
  //            std::unique_ptr<gcs::VertexPositionGeometry>>
  // readMeshFile(std::string inputMesh, std::string referenceMesh);

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
             std::unique_ptr<gcs::VertexPositionGeometry>,
             Eigen::Matrix<bool, Eigen::Dynamic, 1>>
  readTrajFile(std::string trajFile, int startingFrame);
#endif

  // ==========================================================
  // ================     init.cpp           ==================
  // ==========================================================

  inline void updateReferenceMeshFromVPG() {
    refVpg = vpg->copy();
    refVpg->requireEdgeLengths();
    refVpg->requireFaceAreas();
  }

  /**
   * @brief Update the vertex position and recompute cached values
   * (all quantities that characterizes the current energy state)
   * Careful: 1. when using eigenMap: memory address may change after update!!
   */
  void updateConfigurations();

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
  gc::Vector3 computeCornerAngleVariation(gcs::Halfedge he, gcs::Vertex v);
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
