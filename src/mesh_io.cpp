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

#include <cassert>
#include <cmath>
#include <iostream>

#include "mem3dg/mem3dg"
#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_centers.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/surface_mesh_factories.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "igl/cylinder.h"
#include "igl/loop.h"
#include <math.h>
#include <stdexcept>
#include <tuple>
namespace mem3dg {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void loopSubdivide(std::unique_ptr<gcs::ManifoldSurfaceMesh> &ptrMesh,
                   std::unique_ptr<gcs::VertexPositionGeometry> &ptrVpg,
                   std::size_t nSub) {
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> coords;
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor> faces;

  igl::loop(gc::EigenMap<double, 3>(ptrVpg->inputVertexPositions),
            ptrMesh->getFaceVertexMatrix<std::size_t>(), coords, faces, nSub);

  std::tie(ptrMesh, ptrVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(coords, faces);
}

std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
           Eigen::Matrix<double, Eigen::Dynamic, 3>>
loopSubdivide(Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &faces,
              Eigen::Matrix<double, Eigen::Dynamic, 3> &coords,
              std::size_t nSub) {

  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor> newFaces;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> newCoords;

  igl::loop(coords, faces, newCoords, newFaces, nSub);

  return std::tie(newFaces, newCoords);
}

void linearSubdivide(std::unique_ptr<gcs::ManifoldSurfaceMesh> &mesh,
                     std::unique_ptr<gcs::VertexPositionGeometry> &vpg,
                     std::size_t nSub) {
  for (std::size_t iter = 0; iter < nSub; ++iter) {
    gcs::VertexData<bool> isOrigVert(*mesh, true);
    gcs::EdgeData<bool> isOrigEdge(*mesh, true);
    std::vector<gcs::Edge> toFlip;

    for (gcs::Edge e : mesh->edges()) { // loop over all edges
      if (!isOrigEdge[e])
        continue; // don't keep processing new edges

      // gather both vertices incident on the edge, and their positions
      gcs::Vertex oldA = e.halfedge().tipVertex();
      gcs::Vertex oldB = e.halfedge().tailVertex();
      gc::Vector3 oldAPos = vpg->inputVertexPositions[oldA];
      gc::Vector3 oldBPos = vpg->inputVertexPositions[oldB];

      // split the edge
      gcs::Vertex newV = mesh->splitEdgeTriangular(e).vertex();
      isOrigVert[newV] = false;

      // position the new vertex
      gc::Vector3 newPos = 0.5 * (oldAPos + oldBPos);
      vpg->inputVertexPositions[newV] = newPos;

      // iterate through the edges incident on the new vertex
      for (gcs::Edge e : newV.adjacentEdges()) {
        isOrigEdge[e] = false;                    // mark the new edges
        gcs::Vertex otherV = e.otherVertex(newV); // other side of edge

        // if this is a new edge between an old an new vertex, save for flipping
        if (isOrigVert[otherV] && otherV != oldA && otherV != oldB) {
          toFlip.push_back(e);
        }
      }
    }

    for (gcs::Edge e : toFlip) {
      mesh->flip(e);
    }
  }
  mesh->compress();
  std::tie(mesh, vpg) = gcs::makeManifoldSurfaceMeshAndGeometry(
      gc::EigenMap<double, 3, Eigen::RowMajor>(vpg->inputVertexPositions),
      mesh->getFaceVertexMatrix<std::size_t>());
}

std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
           Eigen::Matrix<double, Eigen::Dynamic, 3>>
linearSubdivide(Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> &faces,
                Eigen::Matrix<double, Eigen::Dynamic, 3> &coords,
                std::size_t nSub) {
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> newCoords;
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor> newFaces;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(coords, faces);
  linearSubdivide(ptrMesh, ptrVpg, nSub);
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> meshMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  meshMatrix = ptrMesh->getFaceVertexMatrix<std::size_t>();
  vertexMatrix = gc::EigenMap<double, 3>(ptrVpg->inputVertexPositions);
  return std::tie(meshMatrix, vertexMatrix);
}

std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
           Eigen::Matrix<double, Eigen::Dynamic, 3>>
getCylinderMatrix(double R, int nR, int nh, double freq, double amp) {
  if (nR < 3 || nR < 2) {
    mem3dg_runtime_error("getCylinderMatrix: nR > 2 and nh > 1!");
  }
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> coords;
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor> faces;
  coords.resize(nR * nh, 3);
  faces.resize(2 * nR * (nh - 1), 3);

  double sideLength = 2 * R * sin(3.141592654 / nR);
  double totalHeight = (nh - 1) * sideLength / 2 * sqrt(3);

  int faceIndex = 0;
  for (int clock = 0; clock < nR; clock++) {
    for (int floor = 0; floor < nh; floor++) {
      double z = floor * sideLength / 2 * sqrt(3);
      double localR = R + amp * sin(freq * 2 * constants::PI / totalHeight * z);
      double x = localR * cos(2 * constants::PI * clock / nR +
                              constants::PI / nR * floor);
      double y = localR * sin(2 * constants::PI * clock / nR +
                              constants::PI / nR * floor);
      coords(clock + floor * nR, 0) = x;
      coords(clock + floor * nR, 1) = y;
      coords(clock + floor * nR, 2) = z;
      if (floor > 0) {
        faces(faceIndex, 0) = ((clock + 0) % nR) + (floor - 1) * nR;
        faces(faceIndex, 1) = ((clock + 1) % nR) + (floor - 1) * nR;
        faces(faceIndex, 2) = ((clock + 0) % nR) + (floor + 0) * nR;
        faceIndex++;
        faces(faceIndex, 0) = ((clock + 1) % nR) + (floor - 1) * nR;
        faces(faceIndex, 1) = ((clock + 1) % nR) + (floor + 0) * nR;
        faces(faceIndex, 2) = ((clock + 0) % nR) + (floor + 0) * nR;
        faceIndex++;
      }
    }
  }

  if (faceIndex != faces.rows()) {
    mem3dg_runtime_error("getCylinderMatrix: faceIndex not agree!");
  }
  return std::tie(faces, coords);
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
cylinder(double R, int nR, int nh) {
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> coords;
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor> faces;
  std::tie(faces, coords) = getCylinderMatrix(R, nR, nh);
  return gcs::makeManifoldSurfaceMeshAndGeometry(coords, faces);
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
icosphere(double R, int nSub) {

  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;

  // Initialize vertex coordinates
  static const double t = (1.0 + sqrt(5.0)) / 2.0;
  auto makeNormedVertex = [](double x, double y, double z) -> gc::Vector3 {
    return gc::Vector3{std::move(x), std::move(y), std::move(z)}.normalize();
  };

  coords.emplace_back(makeNormedVertex(-1, t, 0));
  coords.emplace_back(makeNormedVertex(1, t, 0));
  coords.emplace_back(makeNormedVertex(-1, -t, 0));
  coords.emplace_back(makeNormedVertex(1, -t, 0));
  coords.emplace_back(makeNormedVertex(0, -1, t));
  coords.emplace_back(makeNormedVertex(0, 1, t));
  coords.emplace_back(makeNormedVertex(0, -1, -t));
  coords.emplace_back(makeNormedVertex(0, 1, -t));
  coords.emplace_back(makeNormedVertex(t, 0, -1));
  coords.emplace_back(makeNormedVertex(t, 0, 1));
  coords.emplace_back(makeNormedVertex(-t, 0, -1));
  coords.emplace_back(makeNormedVertex(-t, 0, 1));

  // Initialize Faces
  polygons.emplace_back(std::vector<std::size_t>{0, 11, 5});
  polygons.emplace_back(std::vector<std::size_t>{0, 5, 1});
  polygons.emplace_back(std::vector<std::size_t>{0, 1, 7});
  polygons.emplace_back(std::vector<std::size_t>{0, 7, 10});
  polygons.emplace_back(std::vector<std::size_t>{0, 10, 11});
  polygons.emplace_back(std::vector<std::size_t>{1, 5, 9});
  polygons.emplace_back(std::vector<std::size_t>{5, 11, 4});
  polygons.emplace_back(std::vector<std::size_t>{11, 10, 2});
  polygons.emplace_back(std::vector<std::size_t>{10, 7, 6});
  polygons.emplace_back(std::vector<std::size_t>{7, 1, 8});
  polygons.emplace_back(std::vector<std::size_t>{3, 9, 4});
  polygons.emplace_back(std::vector<std::size_t>{3, 4, 2});
  polygons.emplace_back(std::vector<std::size_t>{3, 2, 6});
  polygons.emplace_back(std::vector<std::size_t>{3, 6, 8});
  polygons.emplace_back(std::vector<std::size_t>{3, 8, 9});
  polygons.emplace_back(std::vector<std::size_t>{4, 9, 5});
  polygons.emplace_back(std::vector<std::size_t>{2, 4, 11});
  polygons.emplace_back(std::vector<std::size_t>{6, 2, 10});
  polygons.emplace_back(std::vector<std::size_t>{8, 6, 7});
  polygons.emplace_back(std::vector<std::size_t>{9, 8, 1});

  auto getMidPoint = [&coords](std::size_t t1, std::size_t t2) -> std::size_t {
    coords.emplace_back(((coords[t1] + coords[t2]) / 2).normalize());
    return coords.size() - 1;
  };

  // Preallocate space
  std::size_t finalSize = polygons.size() * std::pow(4, nSub);
  std::vector<std::vector<std::size_t>> polygons_new;
  polygons_new.reserve(finalSize);
  polygons.reserve(finalSize);

  // Subdivide n times by quadrisection
  for (std::size_t iter = 0; iter < nSub; ++iter) {
    std::size_t sz = polygons.size();
    polygons_new.clear();
    for (std::size_t f = 0; f < sz; ++f) {
      std::size_t a = getMidPoint(polygons[f][0], polygons[f][1]);
      std::size_t b = getMidPoint(polygons[f][1], polygons[f][2]);
      std::size_t c = getMidPoint(polygons[f][2], polygons[f][0]);

      polygons_new.emplace_back(std::vector<std::size_t>{polygons[f][0], a, c});
      polygons_new.emplace_back(std::vector<std::size_t>{polygons[f][1], b, a});
      polygons_new.emplace_back(std::vector<std::size_t>{polygons[f][2], c, b});
      polygons_new.emplace_back(std::vector<std::size_t>{a, b, c});
    }
    std::swap(polygons, polygons_new);
  }

  // scale the icosphere
  if (R != 1) {
    for (std::size_t iter = 0; iter < coords.size(); ++iter) {
      coords[iter] *= R;
    }
  }

  gcs::SimplePolygonMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();
  return gcs::makeManifoldSurfaceMeshAndGeometry(soup.polygons,
                                                 soup.vertexCoordinates);
}

std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
           Eigen::Matrix<double, Eigen::Dynamic, 3>>
getIcosphereMatrix(double R, int nSub) {
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::tie(mesh, vpg) = icosphere(R, nSub);
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> meshMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  meshMatrix = mesh->getFaceVertexMatrix<std::size_t>();
  vertexMatrix = gc::EigenMap<double, 3>(vpg->inputVertexPositions);
  return std::tie(meshMatrix, vertexMatrix);
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
hexagon(double R, int nSub) {
  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;

  // Initialize vertex coordinates
  auto makeNormedVertex = [](double x, double y, double z) -> gc::Vector3 {
    return gc::Vector3{std::move(x), std::move(y), std::move(z)}.normalize();
  };

  for (std::size_t i = 0; i < 6; i++) {
    coords.emplace_back(gc::Vector3{R * cos(constants::PI / 3 * i),
                                    R * sin(constants::PI / 3 * i), 0});
  }
  coords.emplace_back(gc::Vector3{0, 0, 0});

  // Initialize Faces
  polygons.emplace_back(std::vector<std::size_t>{6, 0, 1});
  polygons.emplace_back(std::vector<std::size_t>{6, 1, 2});
  polygons.emplace_back(std::vector<std::size_t>{6, 2, 3});
  polygons.emplace_back(std::vector<std::size_t>{6, 3, 4});
  polygons.emplace_back(std::vector<std::size_t>{6, 4, 5});
  polygons.emplace_back(std::vector<std::size_t>{6, 5, 0});

  gcs::SimplePolygonMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::tie(mesh, vpg) = gcs::makeManifoldSurfaceMeshAndGeometry(
      soup.polygons, soup.vertexCoordinates);
  linearSubdivide(mesh, vpg, nSub);

  return std::make_tuple(std::move(mesh), std::move(vpg));
}

std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
           Eigen::Matrix<double, Eigen::Dynamic, 3>>
getHexagonMatrix(double R, int nSub) {
  // std::cout << "hello" << std::endl;
  // Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> meshMatrix;
  // Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::tie(mesh, vpg) = hexagon(R, nSub);
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> meshMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  meshMatrix = mesh->getFaceVertexMatrix<std::size_t>();
  vertexMatrix = gc::EigenMap<double, 3>(vpg->inputVertexPositions);
  return std::tie(meshMatrix, vertexMatrix);
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
tetrahedron() {
  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;

  // Initialize vertex coordinates
  auto makeNormedVertex = [](double x, double y, double z) -> gc::Vector3 {
    return gc::Vector3{std::move(x), std::move(y), std::move(z)}.normalize();
  };

  coords.emplace_back(makeNormedVertex(0, 0, 2));
  coords.emplace_back(makeNormedVertex(1.632993, -0.942809, -0.666667));
  coords.emplace_back(makeNormedVertex(0.000000, 1.885618, -0.666667));
  coords.emplace_back(makeNormedVertex(-1.632993, -0.942809, -0.666667));

  // Initialize Faces
  polygons.emplace_back(std::vector<std::size_t>{1, 0, 3});
  polygons.emplace_back(std::vector<std::size_t>{2, 0, 1});
  polygons.emplace_back(std::vector<std::size_t>{3, 0, 2});
  polygons.emplace_back(std::vector<std::size_t>{3, 2, 1});

  gcs::SimplePolygonMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();
  return gcs::makeManifoldSurfaceMeshAndGeometry(soup.polygons,
                                                 soup.vertexCoordinates);
}

std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
           Eigen::Matrix<double, Eigen::Dynamic, 3>>
getTetrahedronMatrix() {
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::tie(mesh, vpg) = tetrahedron();
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> meshMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  meshMatrix = mesh->getFaceVertexMatrix<std::size_t>();
  vertexMatrix = gc::EigenMap<double, 3>(vpg->inputVertexPositions);
  return std::tie(meshMatrix, vertexMatrix);
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
diamond(double dihedral) {
  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;

  // Initialize vertex coordinates
  auto makeVertex = [](double x, double y, double z) -> gc::Vector3 {
    return gc::Vector3{std::move(x), std::move(y), std::move(z)};
  };

  coords.emplace_back(makeVertex(0, 0.5, 0));
  coords.emplace_back(makeVertex(0, -0.5, 0));
  coords.emplace_back(makeVertex(cos(dihedral) * std::sqrt(3) / 2, 0,
                                 -sin(dihedral) * std::sqrt(3) / 2));
  coords.emplace_back(makeVertex(-std::sqrt(3) / 2, 0, 0));
  // Initialize Faces
  polygons.emplace_back(std::vector<std::size_t>{0, 1, 2});
  polygons.emplace_back(std::vector<std::size_t>{0, 3, 1});

  gcs::SimplePolygonMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();
  return gcs::makeManifoldSurfaceMeshAndGeometry(soup.polygons,
                                                 soup.vertexCoordinates);
}

std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
           Eigen::Matrix<double, Eigen::Dynamic, 3>>
getDiamondMatrix(double dihedral) {
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::tie(mesh, vpg) = diamond(dihedral);
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> meshMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  meshMatrix = mesh->getFaceVertexMatrix<std::size_t>();
  vertexMatrix = gc::EigenMap<double, 3>(vpg->inputVertexPositions);
  return std::tie(meshMatrix, vertexMatrix);
}

std::tuple<Eigen::Matrix<std::size_t, Eigen::Dynamic, 3>,
           Eigen::Matrix<double, Eigen::Dynamic, 3>>
getFaceAndVertexMatrix(std::string &plyName) {
  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3> meshMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrVpg) = gcs::readManifoldSurfaceMesh(plyName);
  meshMatrix = ptrMesh->getFaceVertexMatrix<std::size_t>();
  // Copy assignment?
  vertexMatrix = gc::EigenMap<double, 3>(ptrVpg->inputVertexPositions);
  return std::tie(meshMatrix, vertexMatrix);
}

Eigen::Matrix<double, Eigen::Dynamic, 1> getData(std::string &plyName,
                                                 std::string &elementName,
                                                 std::string &propertyName) {
  // Declare pointers to mesh, geometry and richdata objects
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
      gcs::RichSurfaceMeshData::readMeshAndData(plyName);
  // return ptrRichData->getVertexProperty<double>(propertyName).raw();
  // This is part of RichSurfaceMeshData::getElementProperty
  std::vector<double> rawData = ptrRichData->plyData.getElement(elementName)
                                    .getProperty<double>(propertyName);
  if (rawData.size() != ptrRichData->plyData.getElement(elementName).count) {
    mem3dg_runtime_error("Property " + propertyName +
                         " does not have size equal to number of " +
                         elementName);
  }

  Eigen::Matrix<double, Eigen::Dynamic, 1> property;
  property.resize(rawData.size(), 1);

  for (std::size_t i = 0; i < rawData.size(); i++) {
    property[i] = rawData[i];
  }

  return property;
}

std::vector<std::string> getDataPropertyName(std::string &plyName,
                                             std::string &elementName) {
  // Declare pointers to mesh, geometry and richdata objects
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
      gcs::RichSurfaceMeshData::readMeshAndData(plyName);
  return ptrRichData->plyData.getElement(elementName).getPropertyNames();
}

std::vector<std::string> getDataElementName(std::string &plyName) {
  // Declare pointers to mesh, geometry and richdata objects
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData;
  std::tie(ptrMesh, ptrRichData) =
      gcs::RichSurfaceMeshData::readMeshAndData(plyName);
  return ptrRichData->plyData.getElementNames();
}

std::tuple<Eigen::Matrix<size_t, Eigen::Dynamic, 3>,
           Eigen::Matrix<double, Eigen::Dynamic, 3>>
processSoup(std::string &plyName) {
  Eigen::Matrix<size_t, Eigen::Dynamic, 3> meshMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 3> vertexMatrix;
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  gcs::SimplePolygonMesh soup(plyName);
  soup.mergeIdenticalVertices();
  soup.stripFacesWithDuplicateVertices();
  soup.triangulate();
  std::tie(ptrMesh, ptrVpg) = gcs::makeManifoldSurfaceMeshAndGeometry(
      soup.polygons, soup.vertexCoordinates);
  meshMatrix = ptrMesh->getFaceVertexMatrix<size_t>();
  vertexMatrix = gc::EigenMap<double, 3>(ptrVpg->inputVertexPositions);
  return std::tie(meshMatrix, vertexMatrix);
}

std::size_t getVertexClosestToEmbeddedCoordinate(
    const EigenVectorX3dr &vertexMatrix,
    const std::array<double, 3> &embeddedCoordinate,
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> &filter,
    const std::array<bool, 3> &accountedCoordinate) {
  if (filter.rows() != vertexMatrix.rows())
    mem3dg_runtime_error(
        "number of rows of filter does not match with the vertexMatrix!");

  std::size_t dimension = 0;
  dimension = std::accumulate(accountedCoordinate.begin(),
                              accountedCoordinate.end(), dimension);
  if (dimension < 1)
    mem3dg_runtime_error(
        "number of accounted coordinate can not be less than 1!");
  std::size_t closestVertex;
  double shortestDistanceSq = std::numeric_limits<double>::max();
  for (std::size_t v = 0; v < vertexMatrix.rows(); ++v) {
    if (!filter[v])
      continue;
    double distanceSq = 0.0;
    for (std::size_t i = 0; i < accountedCoordinate.size(); ++i) {
      if (accountedCoordinate[i])
        distanceSq += pow(vertexMatrix.row(v)[i] - embeddedCoordinate[i], 2);
    }
    if (distanceSq < shortestDistanceSq) {
      shortestDistanceSq = distanceSq;
      closestVertex = v;
    }
  }
  return closestVertex;
}

std::size_t getVertexClosestToEmbeddedCoordinate(
    const EigenVectorX3dr &vertexMatrix,
    const std::array<double, 3> &embeddedCoordinate,
    const std::array<bool, 3> &accountedCoordinate) {
  Eigen::Matrix<bool, Eigen::Dynamic, 1> filter;
  filter.setConstant(vertexMatrix.rows(), true);
  return getVertexClosestToEmbeddedCoordinate(vertexMatrix, embeddedCoordinate,
                                              filter, accountedCoordinate);
}

std::tuple<std::size_t, std::array<double, 3>>
getFaceSurfacePointClosestToEmbeddedCoordinate(
    gcs::VertexPositionGeometry &vpg,
    const std::array<double, 3> &embeddedCoordinate_,
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> &filter,
    const std::array<bool, 3> &accountedCoordinate) {
  std::size_t dimension = 0;
  dimension = std::accumulate(accountedCoordinate.begin(),
                              accountedCoordinate.end(), dimension);
  if (dimension < 1)
    mem3dg_runtime_error(
        "number of accounted coordinate can not be less than 1!");

  // initialize embedded point and the closest vertex
  gc::Vector3 embeddedCoordinate{embeddedCoordinate_[0], embeddedCoordinate_[1],
                                 embeddedCoordinate_[2]};
  gcs::Vertex closestVertex =
      vpg.mesh.vertex(getVertexClosestToEmbeddedCoordinate(
          toMatrix(vpg.vertexPositions), embeddedCoordinate_, filter,
          accountedCoordinate));
  std::size_t face;
  std::array<double, 3> coord;

  if (dimension == 1) {
    face = closestVertex.halfedge().face().getIndex();
    gc::Vector3 corner = correspondBarycentricCoordinates(
        gc::Vector3{1, 0, 0}, closestVertex.halfedge());
    coord = {corner.x, corner.y, corner.z};
    return std::tie(face, coord);
  }

  gc::Vector3 projectionNormal{0, 0, 0};
  std::array<std::size_t, 2> directedPlane;
  if (dimension == 2) {
    for (std::size_t i = 0; i < accountedCoordinate.size(); ++i) {
      if (!accountedCoordinate[i]) {
        projectionNormal[i] = 1;
        directedPlane = {i % 3, (i + 1) % 3};
      }
    }
  }

  double shortestDistance = std::numeric_limits<double>::max();
  // loop through faces at the neighborhood
  for (gcs::Halfedge he : closestVertex.outgoingHalfedges()) {
    if (!he.isInterior())
      continue;

    if (dimension == 3) {
      projectionNormal = vpg.faceNormal(he.face());
      directedPlane = {0, 1};
      double max = abs(projectionNormal.z);
      for (std::size_t i = 1; i < 3; ++i) {
        if (abs(projectionNormal[i]) > max) {
          max = abs(projectionNormal[i]);
          directedPlane = {i, (i + 1) % 3};
        }
      }
    }

    // project the embedded point
    gc::Vector3 projectedEmbeddedCoordinate =
        embeddedCoordinate -
        gc::dot(embeddedCoordinate - vpg.vertexPositions[closestVertex],
                projectionNormal) *
            projectionNormal;

    // determine the best behaving choice to do inverse barycentric map
    gc::Vector2
        faceVertexOnDirectedPlane1 =
            gc::Vector2{vpg.vertexPositions[he.vertex()][directedPlane[0]],
                        vpg.vertexPositions[he.vertex()][directedPlane[1]]},
        faceVertexOnDirectedPlane2 =
            gc::Vector2{
                vpg.vertexPositions[he.next().vertex()][directedPlane[0]],
                vpg.vertexPositions[he.next().vertex()][directedPlane[1]]},
        faceVertexOnDirectedPlane3 =
            gc::Vector2{vpg.vertexPositions[he.next().next().vertex()]
                                           [directedPlane[0]],
                        vpg.vertexPositions[he.next().next().vertex()]
                                           [directedPlane[1]]},
        projectedEmbeddedCoordinateOnDirectedPlane =
            gc::Vector2{projectedEmbeddedCoordinate[directedPlane[0]],
                        projectedEmbeddedCoordinate[directedPlane[1]]};
    gc::Vector3 candidateCoord = cartesianToBarycentric(
        faceVertexOnDirectedPlane1, faceVertexOnDirectedPlane2,
        faceVertexOnDirectedPlane3, projectedEmbeddedCoordinateOnDirectedPlane);

    // find the closest surface point
    candidateCoord = gc::componentwiseMax(candidateCoord, gc::Vector3{0, 0, 0});
    candidateCoord /= gc::sum(candidateCoord);
    candidateCoord = correspondBarycentricCoordinates(candidateCoord, he);
    gcs::SurfacePoint candidate(he.face(), candidateCoord);
    double distance =
        (embeddedCoordinate - candidate.interpolate(vpg.vertexPositions))
            .norm();
    if (distance < shortestDistance) {
      face = he.face().getIndex();
      coord = {
          candidateCoord.x, candidateCoord.y,
          candidateCoord.z}; // note that this is barycentric, not cartesian
      shortestDistance = distance;
    }
  }
  return std::tie(face, coord);
}

DLL_PUBLIC std::tuple<std::size_t, std::array<double, 3>>
getFaceSurfacePointClosestToEmbeddedCoordinate(
    const EigenVectorX3sr &faceMatrix, const EigenVectorX3dr &vertexMatrix,
    const std::array<double, 3> &embeddedCoordinate,
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> &filter,
    const std::array<bool, 3> &accountedCoordinate) {
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::tie(mesh, vpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(vertexMatrix, faceMatrix);

  return getFaceSurfacePointClosestToEmbeddedCoordinate(
      *vpg, embeddedCoordinate, filter, accountedCoordinate);
}
DLL_PUBLIC std::tuple<std::size_t, std::array<double, 3>>
getFaceSurfacePointClosestToEmbeddedCoordinate(
    const EigenVectorX3sr &faceMatrix, const EigenVectorX3dr &vertexMatrix,
    const std::array<double, 3> &embeddedCoordinate,
    const std::array<bool, 3> &accountedCoordinate) {
  Eigen::Matrix<bool, Eigen::Dynamic, 1> filter;
  filter.setConstant(vertexMatrix.rows(), true);
  return getFaceSurfacePointClosestToEmbeddedCoordinate(
      faceMatrix, vertexMatrix, embeddedCoordinate, filter,
      accountedCoordinate);
}

DLL_PUBLIC std::size_t
getVertexFurthestFromBoundary(gc::ManifoldSurfaceMesh &mesh,
                              gc::VertexPositionGeometry &vpg) {
  std::vector<gcs::Vertex> boundaryLoop;
  for (gcs::BoundaryLoop bl : mesh.boundaryLoops()) {
    for (gcs::Vertex v : bl.adjacentVertices()) {
      boundaryLoop.push_back(v);
    }
  }
  // if (notableVertex.raw().cast<int>().sum() != 1)
  //   mem3dg_runtime_error("the number of found open mesh center is not 1");
  return gcs::findCenter(mesh, vpg, boundaryLoop, 2).nearestVertex().getIndex();
}

DLL_PUBLIC std::size_t
getVertexFurthestFromBoundary(const EigenVectorX3sr &faceMatrix,
                              const EigenVectorX3dr &vertexMatrix) {
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::tie(mesh, vpg) =
      gcs::makeManifoldSurfaceMeshAndGeometry(vertexMatrix, faceMatrix);
  return getVertexFurthestFromBoundary(*mesh, *vpg);
}

} // namespace mem3dg
