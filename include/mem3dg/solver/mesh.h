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
 * @param coords    Reference to vector of vertex coordinates
 * @param polygons  Polygons of the mesh
 * @param n         Iterations of quadrisections to perform
 */
DLL_PUBLIC void icosphere(std::vector<gc::Vector3> &coords,
                          std::vector<std::vector<std::size_t>> &polygons,
                          int n, double R = 1);

/**
 * @brief Load a reference geometry object on existing mesh object
 *
 * @param ptrMesh    Reference to the unique pointer to the mesh object
 * @param ptrRefVpg  Pointer to the reference geometry object     
 * @param coords     Coordinates Eigen matrix V x 3     
 */
DLL_PUBLIC void
loadRefMesh(std::unique_ptr<gcs::ManifoldSurfaceMesh> &ptrMesh,
            gcs::VertexPositionGeometry *&ptrRefVpg,
            Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> coords);

/**
 * @brief Hard code a tetrahedron in PolygonSoup form
 *
 * @param coords    Reference to vector of vertex coordinates
 * @param polygons  Polygons of the mesh
 */
DLL_PUBLIC
void tetrahedron(std::vector<gc::Vector3> &coords,
                 std::vector<std::vector<std::size_t>> &polygons);

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

/**
 * @brief generate a icosphere in .ply format to path
 *
 * @param nSub Iterations of quadrisections to perform
 * @param path String of output path
 * @param R    Radius of the sphere
 */
DLL_PUBLIC int genIcosphere(size_t nSub, std::string path, double R);

} // end namespace ddgsolver
