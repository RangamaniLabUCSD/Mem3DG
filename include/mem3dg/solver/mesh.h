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

#include "geometrycentral/utilities/vector3.h"
#include "macros.h"
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

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
DLL_PUBLIC std::tuple<Eigen::Matrix<size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
getCylinderMatrix(double R, int nR, int nh, double freq = 1, double amp = 0);

/**
 * @brief Construct an icosphere mesh in PolygonSoup form
 *
 * @param n         Iterations of quadrisections to perform
 * @param R         Radius of the icosphere
 */
DLL_PUBLIC std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
                      std::unique_ptr<gcs::VertexPositionGeometry>>
icosphere(int n, double R);

/**
 * @brief Construct an icosphere mesh
 *
 * @param n         Iterations of quadrisections to perform
 * @param R         Radius of the icosphere
 */
DLL_PUBLIC std::tuple<Eigen::Matrix<size_t, Eigen::Dynamic, 3>,
                      Eigen::Matrix<double, Eigen::Dynamic, 3>>
getIcosphereMatrix(int n, double R);

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
DLL_PUBLIC std::tuple<Eigen::Matrix<size_t, Eigen::Dynamic, 3>,
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
DLL_PUBLIC std::tuple<Eigen::Matrix<size_t, Eigen::Dynamic, 3>,
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

} // namespace mem3dg
