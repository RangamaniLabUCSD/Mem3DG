/*
 * Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
 *
 * Copyright 2020- The Mem3DG Authors
 * and the project initiators Cuncheng Zhu, Christopher T. Lee, and
 * Padmini Rangamani.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Please help us support Mem3DG development by citing the research
 * papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/
 * for more information.
 */

#include "mem3dg/solver/geometry.h"
#include "mem3dg/macros.h"
#include <geometrycentral/surface/surface_mesh.h>

#ifdef MEM3DG_WITH_NETCDF
#include <netcdf>
#endif

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {
namespace solver {

void Geometry::saveGeometry(std::string PathToSave) {
  gcs::writeSurfaceMesh(*mesh, *vpg, PathToSave);
}

#ifdef MEM3DG_WITH_NETCDF
namespace nc = ::netCDF;

/**
 * @brief
 *
 * @param trajFile
 * @param startingFrame
 * @return std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
 * std::unique_ptr<gcs::VertexPositionGeometry>,
 * Eigen::Matrix<bool, Eigen::Dynamic, 1>>
 */
std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>,
           Eigen::Matrix<bool, Eigen::Dynamic, 1>>
Geometry::readTrajFile(std::string trajFile, int startingFrame) {
  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh_;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg_;

  MutableTrajFile fd = MutableTrajFile::openReadOnly(trajFile);
  fd.getNcFrame(startingFrame);
  std::tie(mesh_, vpg_) = gcs::makeManifoldSurfaceMeshAndGeometry(
      fd.getCoords(startingFrame), fd.getTopology(startingFrame));

  Eigen::Matrix<bool, Eigen::Dynamic, 1> notableVertex_here;
  try {
    notableVertex_here = fd.getNotableVertex(startingFrame);
  } catch (const nc::exceptions::NcException &e) {
    notableVertex_here =
        Eigen::Matrix<bool, Eigen::Dynamic, 1>::Zero(mesh_->nVertices());
    mem3dg_runtime_warning("Trajfile has no attribute `notable vertex`");
  }

  return std::make_tuple(std::move(mesh_), std::move(vpg_), notableVertex_here);
}
#endif

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
Geometry::readMeshFile(std::string inputMesh) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh_;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg_;

  // Load input mesh and geometry
  std::tie(mesh_, vpg_) = gcs::readManifoldSurfaceMesh(inputMesh);

  return std::make_tuple(std::move(mesh_), std::move(vpg_));
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
Geometry::readMatrices(EigenVectorX3sr &faceVertexMatrix,
                       EigenVectorX3dr &vertexPositionMatrix) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh_;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg_;

  // Load input mesh and geometry
  std::tie(mesh_, vpg_) = gcs::makeManifoldSurfaceMeshAndGeometry(
      vertexPositionMatrix, faceVertexMatrix);

  return std::make_tuple(std::move(mesh_), std::move(vpg_));
}

double Geometry::computeLengthCrossRatio(gcs::VertexPositionGeometry &vpg_,
                                         gcs::Halfedge &he) const {
  if (he.edge().isBoundary()) {
    return 1;
  } else {
    gcs::Edge lj = he.next().edge();
    gcs::Edge ki = he.twin().next().edge();
    gcs::Edge il = he.next().next().edge();
    gcs::Edge jk = he.twin().next().next().edge();
    return vpg_.edgeLengths[il] * vpg_.edgeLengths[jk] / vpg_.edgeLengths[ki] /
           vpg_.edgeLengths[lj];
  }
}

void Geometry::computeFaceTangentialDerivative(
    gcs::VertexData<double> &quantities, gcs::FaceData<gc::Vector3> &gradient) {
  if ((quantities.raw().array() == quantities.raw()[0]).all()) {
    gradient.fill({0, 0, 0});
  } else {
    for (gcs::Face f : mesh->faces()) {
      gc::Vector3 normal = vpg->faceNormals[f];
      gc::Vector3 gradientVec{0, 0, 0};
      for (gcs::Halfedge he : f.adjacentHalfedges()) {
        gradientVec += quantities[he.next().tipVertex()] *
                       gc::cross(normal, vecFromHalfedge(he, *vpg));
      }
      gradient[f] = gradientVec / 2 / vpg->faceAreas[f];
    }
  }
}

EigenVectorX1d Geometry::computeGeodesicDistance() {
  // assert(std::count(notableVertex.raw().begin(), notableVertex.raw().end(),
  //                   true) != 0);
  geodesicDistance.fill(std::numeric_limits<double>::max());

  gcs::HeatMethodDistanceSolver heatSolver(*vpg);
  // gc::Vertex centerVertex;
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    if (notableVertex[i]) {
      geodesicDistance.raw() = geodesicDistance.raw().cwiseMin(
          heatSolver.computeDistance(mesh->vertex(i)).raw());
    }
  }
  return geodesicDistance.raw();
}

void Geometry::updateReferenceConfigurations() {
  // refVpg->requireEdgeLengths();
  // refVpg->requireFaceAreas();
  for (std::size_t i = 0; i < mesh->nHalfedges(); ++i) {
    gcs::Halfedge he{mesh->halfedge(i)};
    refLcrs[he] = computeLengthCrossRatio(*refVpg, he);
  }
}

void Geometry::updateConfigurations() {
  // refresh cached quantities after regularization
  vpg->refreshQuantities();
  volume = getMeshVolume(*mesh, *vpg, true);
  surfaceArea = vpg->faceAreas.raw().sum();
}

} // namespace solver
} // namespace mem3dg
