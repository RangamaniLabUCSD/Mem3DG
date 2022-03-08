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

namespace mem3dg {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

/**
 * @brief read vertex and face matrix from .ply file
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
processSoup(std::string &plyName);

/**
 * @brief read vertex and face matrix from .ply file
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
readMesh(std::string &plyName);

/**
 * @brief Read data in the format of matrix from .ply file
 *
 * @param plyName   PLY file to read
 * @return Vector of
 */
DLL_PUBLIC std::vector<std::string> readData(std::string &plyName);

/**
 * @brief
 *
 * @param plyName
 * @param elementName
 * @return
 */
DLL_PUBLIC std::vector<std::string> readData(std::string &plyName,
                                             std::string &elementName);
/**
 * @brief
 *
 * @param plyName
 * @param elementName
 * @param vertexProperties
 * @return
 */
DLL_PUBLIC Eigen::Matrix<double, Eigen::Dynamic, 1>
readData(std::string &plyName, std::string &elementName,
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
DLL_PUBLIC void subdivide(std::unique_ptr<gcs::ManifoldSurfaceMesh> &mesh,
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
subdivide(Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &faces,
          Eigen::Matrix<double, Eigen::Dynamic, 3> &coords, std::size_t nSub);

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

} // namespace mem3dg
