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

#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include "mem3dg/macros.h"
#include "mem3dg/type_utilities.h"

namespace mem3dg {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

/**
 * @brief process vertex and face matrix soup from .ply file
 *
 * @param plyName path to the .ply file
 * @return tuple of processed face topology matrix and vertex position matrix
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
processSoup(std::string &plyName);

/**
 * @brief read vertex and face matrix from .ply file
 * @param plyName path to the .ply file
 * @return tuple of processed face topology matrix and vertex position matrix
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
getFaceAndVertexMatrix(std::string &plyName);

/**
 * @brief retrieve all richData element name from .ply file. Namely the list of
 * the places where data live in, such as vertex, edge or face.
 * @param plyName PLY file to read
 * @return list of all element names
 */
DLL_PUBLIC std::vector<std::string> getDataElementName(std::string &plyName);

/**
 * @brief
 *
 * @param plyName
 * @param elementName
 * @return
 */
DLL_PUBLIC std::vector<std::string>
getDataPropertyName(std::string &plyName, std::string &elementName);
/**
 * @brief
 *
 * @param plyName
 * @param elementName
 * @param vertexProperties
 * @return
 */
DLL_PUBLIC Eigen::Matrix<double, Eigen::Dynamic, 1>
getData(std::string &plyName, std::string &elementName,
        std::string &vertexProperties);

/**
 * @brief Construct an hexagon mesh in PolygonSoup form
 *
 * @param R         Radius of the the hexagon
 * @param nSub      Number of subdivision
 */
DLL_PUBLIC std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                      std::unique_ptr<gcs::VertexPositionGeometry>>
hexagon(double R, int nSub);

/**
 * @brief Construct an hexagon mesh in PolygonSoup form
 *
 * @param R         Radius of the the hexagon
 * @param nSub      Number of subdivision
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
getHexagonMatrix(double R, int nSub);

/**
 * @brief Construct an icosphere mesh in PolygonSoup form
 *
 * @param R         Radius of the the cylinder
 * @param nR        Number of element radially
 * @param nh        Number of element axially
 */
DLL_PUBLIC std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                      std::unique_ptr<gcs::VertexPositionGeometry>>
cylinder(double R, int nR, int nh);

/**
 * @brief Construct an icosphere mesh in PolygonSoup form
 *
 * @param R         Radius of the the cylinder
 * @param nR        Number of element radially
 * @param nh        Number of element axially
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
getCylinderMatrix(double R, int nR, int nh, double freq = 1, double amp = 0);

/**
 * @brief Construct an icosphere mesh in PolygonSoup form
 *
 * @param R         Radius of the icosphere
 * @param nSub      Iterations of quadrisections to perform
 */
DLL_PUBLIC std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                      std::unique_ptr<gcs::VertexPositionGeometry>>
icosphere(double R, int nSub);

/**
 * @brief Construct an icosphere mesh
 *
 * @param R         Radius of the icosphere
 * @param nSub      Iterations of quadrisections to perform
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
getIcosphereMatrix(double R, int nSub);

/**
 * @brief Hard code a tetrahedron in PolygonSoup form
 *
 */
DLL_PUBLIC
std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
tetrahedron();
/**
 * @brief Hard code a tetrahedron in PolygonSoup form
 *
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
getTetrahedronMatrix();

/**
 * @brief Hard code a diamond in PolygonSoup form
 *
 */
DLL_PUBLIC
std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
diamond(double dihedral);

/**
 * @brief Hard code a diamond in PolygonSoup form
 *
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
getDiamondMatrix(double dihedral);

/**
 * @brief Subdivide a manifold mesh
 *
 * @param mesh Reference to the pointer to the mesh object
 * @param vpg  Reference to the pointer to the geometry object
 * @param nSub Iterations of quadrisections to perform
 */
DLL_PUBLIC void
linearSubdivide(std::unique_ptr<gcs::ManifoldSurfaceMesh> &mesh,
                std::unique_ptr<gcs::VertexPositionGeometry> &vpg,
                std::size_t nSub);

/**
 * @brief Subdivide a mesh in Polygon Soup form
 *
 * @param faces
 * @param coords
 * @param nSub
 * @return
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
linearSubdivide(Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &faces,
                Eigen::Matrix<double, Eigen::Dynamic, 3> &coords,
                std::size_t nSub);

/**
 * @brief Subdivide a manifold mesh in Loop scheme
 *
 * @param mesh Reference to the pointer to the mesh object
 * @param vpg  Reference to the pointer to the geometry object
 * @param nSub Iterations of quadrisections to perform
 */
DLL_PUBLIC void
loopSubdivide(std::unique_ptr<gcs::ManifoldSurfaceMesh> &ptrMesh,
              std::unique_ptr<gcs::VertexPositionGeometry> &ptrVpg,
              std::size_t nSub);

/**
 * @brief Perform loop subdivision on a Polygon Soup mesh
 *
 * @param faces
 * @param coords
 * @param nSub
 * @return
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
loopSubdivide(Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &faces,
              Eigen::Matrix<double, Eigen::Dynamic, 3> &coords,
              std::size_t nSub);

/**
 * @brief Find the face index and the barycentric coordinate of the face surface
 * point closest to a embedded coordinate in Euclidean distance.
 *
 * @param vpg vertex position geometry
 * @param embeddedCoordinate the target embedded coordinate expressed in 3D
 * cartesian coordinate
 * @param accountedCoordinate array to identify the accounted coordinate. For
 * example, it finds closest vertex in the x-y plane when set to be {true, true,
 * false}. Defaults to {true, true, true}
 * @param filter Limit the scope of search to a subset of vertices, Defaults to
 * all true by overloading
 * @return tuple of face index and the barycentric coordinate
 */
DLL_PUBLIC std::tuple<std::size_t, std::array<double, 3>>
getFaceSurfacePointClosestToEmbeddedCoordinate(
    gcs::VertexPositionGeometry &vpg,
    const std::array<double, 3> &embeddedCoordinate_,
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> &filter,
    const std::array<bool, 3> &accountedCoordinate = std::array<bool, 3>{
        true, true, true});

/**
 * @brief Find the face index and the barycentric coordinate of the face surface
 * point closest to a embedded coordinate in Euclidean distance.
 *
 * @param faceMatrix face topology matrix of the mesh (F x 3)
 * @param vertexMatrix vertex position matrix of the mesh (V x 3)
 * @param embeddedCoordinate the target embedded coordinate expressed in 3D
 * cartesian coordinate
 * @param filter Limit the scope of search to a subset of vertices, Defaults to
 * all true by overloading
 * @param accountedCoordinate array to identify the accounted coordinate. For
 * example, it finds closest vertex in the x-y plane when set to be {true, true,
 * false}. Defaults to {true, true, true}
 *
 * @return tuple of face index and the barycentric coordinate
 */
DLL_PUBLIC std::tuple<std::size_t, std::array<double, 3>>
getFaceSurfacePointClosestToEmbeddedCoordinate(
    const EigenVectorX3sr &faceMatrix, const EigenVectorX3dr &vertexMatrix,
    const std::array<double, 3> &embeddedCoordinate,
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> &filter,
    const std::array<bool, 3> &accountedCoordinate = std::array<bool, 3>{
        true, true, true});
DLL_PUBLIC std::tuple<std::size_t, std::array<double, 3>>
getFaceSurfacePointClosestToEmbeddedCoordinate(
    const EigenVectorX3sr &faceMatrix, const EigenVectorX3dr &vertexMatrix,
    const std::array<double, 3> &embeddedCoordinate,
    const std::array<bool, 3> &accountedCoordinate = std::array<bool, 3>{
        true, true, true});
/**
 * @brief Find the index of the vertex closest to a embedded coordinate in
 * Euclidean distance.
 *
 * @param vertexMatrix vertex position matrix of the mesh (N x 3)
 * @param embeddedCoordinate the target embedded coordinate expressed in 3D
 * cartesian coordinate
 * @param accountedCoordinate array to identify the accounted coordinate.
 * For example, it finds closest vertex in the x-y plane when set to be
 * {true, true, false}. Defaults to {true, true, true}
 * @param filter Limit the scope of search to a subset of vertices, Defaults
 * to all true by overloading
 *
 * @return vertex index
 */
DLL_PUBLIC std::size_t getVertexClosestToEmbeddedCoordinate(
    const EigenVectorX3dr &vertexMatrix,
    const std::array<double, 3> &embeddedCoordinate,
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> &filter,
    const std::array<bool, 3> &accountedCoordinate = std::array<bool, 3>{
        true, true, true});
DLL_PUBLIC std::size_t getVertexClosestToEmbeddedCoordinate(
    const EigenVectorX3dr &vertexMatrix,
    const std::array<double, 3> &embeddedCoordinate,
    const std::array<bool, 3> &accountedCoordinate = std::array<bool, 3>{
        true, true, true});

/**
 * @brief Find the index of the vertex furthest away from the boundaries
 *
 * @param mesh surface mesh
 * @param vpg vertex position geometry
 *
 * @return vertex index
 */
DLL_PUBLIC std::size_t
getVertexFurthestFromBoundary(gc::ManifoldSurfaceMesh &mesh,
                              gc::VertexPositionGeometry &vpg);
DLL_PUBLIC std::size_t
getVertexFurthestFromBoundary(const EigenVectorX3sr &faceMatrix,
                              const EigenVectorX3dr &vertexMatrix);

} // namespace mem3dg
