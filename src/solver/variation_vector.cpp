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

#include "mem3dg/solver/geometry.h"

namespace mem3dg {
namespace solver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

gc::Vector3 Geometry::computeCornerAngleVariation(gcs::Corner c,
                                                  gcs::Vertex v) {
  gcs::Halfedge he = c.halfedge();
  gc::Vector3 n = vpg->faceNormals[c.face()];
  gc::Vector3 ej = vecFromHalfedge(he, *vpg);
  gc::Vector3 ei = vecFromHalfedge(he.next(), *vpg);
  gc::Vector3 ek = vecFromHalfedge(he.next().next(), *vpg);
  if (c.vertex() == v) { // vi
    gc::Vector3 grad_anglek = -gc::cross(n, ej).normalize() / gc::norm(ej);
    gc::Vector3 grad_anglej = -gc::cross(n, ek).normalize() / gc::norm(ek);
    return -(grad_anglek + grad_anglej);
  } else if (he.next().vertex() == v) { // vk
    return -gc::cross(n, ej).normalize() / gc::norm(ej);
  } else if (he.next().next().vertex() == v) { // vj
    return -gc::cross(n, ek).normalize() / gc::norm(ek);
  } else {
    mem3dg_runtime_error("Unexpected combination of corner and vertex!");
    return gc::Vector3{0, 0, 0};
  }
}

gc::Vector3 Geometry::computeCornerAngleVariation(gcs::Halfedge he,
                                                  gcs::Vertex v) {
  if (he.isInterior()) {
    return computeCornerAngleVariation(he.corner(), v);
  } else {
    return gc::Vector3{0, 0, 0};
  }
}

gc::Vector3 Geometry::computeDihedralAngleVariation(gcs::Halfedge he,
                                                    gcs::Vertex v) {
  double l = vpg->edgeLengths[he.edge()];
  if (he.edge().isBoundary()) {
    return gc::Vector3{0, 0, 0};
  } else if (he.vertex() == v) {
    return (vpg->halfedgeCotanWeights[he.next().next()] *
                vpg->faceNormals[he.face()] +
            vpg->halfedgeCotanWeights[he.twin().next()] *
                vpg->faceNormals[he.twin().face()]) /
           l;
  } else if (he.next().vertex() == v) {
    return (vpg->halfedgeCotanWeights[he.twin().next().next()] *
                vpg->faceNormals[he.twin().face()] +
            vpg->halfedgeCotanWeights[he.next()] *
                vpg->faceNormals[he.face()]) /
           l;
  } else if (he.next().next().vertex() == v) {
    return (-(vpg->halfedgeCotanWeights[he.next().next()] +
              vpg->halfedgeCotanWeights[he.next()]) *
            vpg->faceNormals[he.face()]) /
           l;
  } else {
    mem3dg_runtime_error("Unexpected combination of halfedge and vertex!");
    return gc::Vector3{0, 0, 0};
  }
}

std::tuple<gc::Vector3, gc::Vector3>
Geometry::computeHalfedgeSchlafliVector(gcs::VertexPositionGeometry &vpg,
                                        gc::Halfedge &he) {
  // this is a raw and less readable implementation than using
  // computeDihedralAngleVariation()
  std::size_t fID = he.face().getIndex();
  std::size_t heID_twin = he.twin().getIndex();
  std::size_t fID_he_twin = he.twin().face().getIndex();
  std::size_t heID_twin_next = he.twin().next().getIndex();
  std::size_t heID_he_next_next = he.next().next().getIndex();
  gc::Vertex vj = he.tipVertex();
  bool boundaryVertex = he.vertex().isBoundary();
  bool boundaryEdge = he.edge().isBoundary();
  bool interiorHalfedge = he.isInterior();
  bool interiorTwinHalfedge = he.twin().isInterior();
  gc::Vector3 schlafliVec1{0, 0, 0};
  gc::Vector3 schlafliVec2{0, 0, 0};
  if (!boundaryEdge) {
    schlafliVec1 =
        vpg.halfedgeCotanWeights[heID_he_next_next] * vpg.faceNormals[fID] +
        vpg.halfedgeCotanWeights[heID_twin_next] * vpg.faceNormals[fID_he_twin];
  }
  if (boundaryVertex && boundaryEdge) {
    schlafliVec2 = interiorHalfedge
                       ? (-(vpg.halfedgeCotanWeights[he] +
                            vpg.halfedgeCotanWeights[heID_he_next_next]) *
                          vpg.faceNormals[fID])
                       : (-(vpg.halfedgeCotanWeights[heID_twin] +
                            vpg.halfedgeCotanWeights[heID_twin_next]) *
                          vpg.faceNormals[fID_he_twin]);
  } else if (!boundaryVertex && vj.isBoundary()) {
    schlafliVec2 =
        vpg.halfedgeCotanWeights[heID_he_next_next] * vpg.faceNormals[fID] +
        vpg.halfedgeCotanWeights[heID_twin_next] * vpg.faceNormals[fID_he_twin];

    if (!he.next().edge().isBoundary())
      schlafliVec2 -= (vpg.halfedgeCotanWeights[he] +
                       vpg.halfedgeCotanWeights[heID_he_next_next]) *
                      vpg.faceNormals[fID];

    if (!he.twin().next().next().edge().isBoundary())
      schlafliVec2 -= (vpg.halfedgeCotanWeights[heID_twin] +
                       vpg.halfedgeCotanWeights[heID_twin_next]) *
                      vpg.faceNormals[fID_he_twin];
  } else {
    schlafliVec2 =
        -(vpg.halfedgeCotanWeights[he] * vpg.faceNormals[fID] +
          vpg.halfedgeCotanWeights[heID_twin] * vpg.faceNormals[fID_he_twin]);
  }
  return std::make_tuple(schlafliVec1, schlafliVec2);
}

gc::Vector3 Geometry::computeHalfedgeGaussianCurvatureVector(
    gcs::VertexPositionGeometry &vpg, gc::Halfedge &he) {
  gc::Vector3 gaussVec{0, 0, 0};
  if (!he.edge().isBoundary()) {
    // gc::Vector3 eji{} = -vecFromHalfedge(he, *vpg);
    gaussVec = 0.5 * vpg.edgeDihedralAngles[he.edge()] *
               (-vecFromHalfedge(he, vpg)).unit();
  }
  return gaussVec;
}

gc::Vector3
Geometry::computeHalfedgeMeanCurvatureVector(gcs::VertexPositionGeometry &vpg,
                                             gc::Halfedge &he) {
  std::size_t fID = he.face().getIndex();
  std::size_t fID_he_twin = he.twin().face().getIndex();
  bool interiorHalfedge = he.isInterior();
  bool interiorTwinHalfedge = he.twin().isInterior();
  gc::Vector3 areaGrad{0, 0, 0};
  if (interiorHalfedge) {
    areaGrad +=
        0.25 * gc::cross(vpg.faceNormals[fID], vecFromHalfedge(he.next(), vpg));
  }
  if (interiorTwinHalfedge) {
    areaGrad += 0.25 * gc::cross(vpg.faceNormals[fID_he_twin],
                                 vecFromHalfedge(he.twin().next().next(), vpg));
  }
  return areaGrad / 2;
}

gc::Vector3
Geometry::computeHalfedgeVolumeVariationVector(gcs::VertexPositionGeometry &vpg,
                                               gc::Halfedge &he) {

  // Note: the missing contribution from faces only contributes to z -
  // axis forces
  // volGrad = vpg->faceNormals[fID] * vpg->faceAreas[fID] /
  // 3;
  std::size_t fID = he.face().getIndex();
  bool interiorHalfedge = he.isInterior();
  gc::Vector3 volGrad{0, 0, 0};
  if (interiorHalfedge) {
    volGrad = vpg.faceNormals[fID] * vpg.faceAreas[fID] / 3;
  }
  return volGrad;
}

gc::VertexData<gc::Vector3>
Geometry::computeVertexSchlafliLaplacianMeanCurvatureVectors(
    gcs::VertexData<double> &spontaneousCurvature) {
  mesh->compress();
  gc::VertexData<gc::Vector3> vector(*mesh, {0, 0, 0});
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex v{mesh->vertex(i)};
    double Hi = vpg->vertexMeanCurvatures[i] / vpg->vertexDualAreas[i];
    double H0i = spontaneousCurvature[i];
    for (gc::Halfedge he : v.outgoingHalfedges()) {
      std::size_t i_vj = he.tipVertex().getIndex();
      double Hj = vpg->vertexMeanCurvatures[i_vj] / vpg->vertexDualAreas[i_vj];
      double H0j = spontaneousCurvature[i_vj];
      gc::Vector3 vec1;
      gc::Vector3 vec2;
      std::tie(vec1, vec2) = computeHalfedgeSchlafliVector(*vpg, he);
      vector[v] += (Hi - H0i) * vec1 + (Hj - H0j) * vec2;
    }
  }
  return vector;
}

gc::VertexData<gc::Vector3> Geometry::computeVertexGaussianCurvatureVectors() {
  return halfedgeVectorToVertexVector(*mesh, *vpg,
                                      computeHalfedgeGaussianCurvatureVector);
}

gc::VertexData<gc::Vector3> Geometry::computeVertexMeanCurvatureVectors() {
  return halfedgeVectorToVertexVector(*mesh, *vpg,
                                      computeHalfedgeMeanCurvatureVector);
}

gc::VertexData<gc::Vector3> Geometry::computeVertexVolumeVariationVectors() {
  return halfedgeVectorToVertexVector(*mesh, *vpg,
                                      computeHalfedgeVolumeVariationVector);
}

gcs::VertexData<gc::Vector3> Geometry::halfedgeVectorToVertexVector(
    gcs::ManifoldSurfaceMesh &mesh, gcs::VertexPositionGeometry &vpg,
    std::function<gc::Vector3(gcs::VertexPositionGeometry &vpg, gc::Halfedge &)>
        computeHalfedgeVariationalVector) {
  mesh.compress();
  gc::VertexData<gc::Vector3> vector(mesh, {0, 0, 0});
  for (std::size_t i = 0; i < mesh.nVertices(); ++i) {
    gc::Vertex v{mesh.vertex(i)};
    for (gc::Halfedge he : v.outgoingHalfedges()) {
      vector[v] += computeHalfedgeVariationalVector(vpg, he);
    }
  }
  return vector;
}

gc::Vector3
Geometry::computeHalfedgeSquaredIntegratedDerivativeNormVariationVector(
    const gcs::VertexData<double> &quantities, const gcs::Halfedge &he) {
  if (!he.isInterior()) {
    throw std::runtime_error(
        "computeGradientNormGradient: halfedge is not interior!");
  }

  // quantities
  double qj = quantities[he.next().next().vertex()];
  double qi = quantities[he.vertex()];
  double qk = quantities[he.next().vertex()];

  if (qj == qi && qj == qk) {
    return gc::Vector3({0, 0, 0});
  } else {
    // Edge and normal vector
    gc::Vector3 n = vpg->faceNormals[he.face()];
    gc::Vector3 ej = vecFromHalfedge(he, *vpg);
    gc::Vector3 ei = vecFromHalfedge(he.next(), *vpg);
    gc::Vector3 ek = vecFromHalfedge(he.next().next(), *vpg);

    // exterior angle of triangles (angles formed by e_perp)
    double anglek = gc::angle(ej, ei);
    double anglej = gc::angle(ei, ek);
    double anglei = gc::angle(ek, ej);

    // gradient of edge length wrt he.vertex()
    gc::Vector3 grad_ejnorm = -ej.normalize();
    gc::Vector3 grad_eknorm = ek.normalize();

    // gradient of exterior angle wrt he.vertex()
    gc::Vector3 grad_anglek =
        -computeCornerAngleVariation(he.next().corner(), he.vertex());
    gc::Vector3 grad_anglej =
        -computeCornerAngleVariation(he.next().next().corner(), he.vertex());
    gc::Vector3 grad_anglei =
        -computeCornerAngleVariation(he.corner(), he.vertex());
    // gc::Vector3 grad_anglek = gc::cross(n, ej).normalize() / gc::norm(ej);
    // gc::Vector3 grad_anglej = gc::cross(n, ek).normalize() / gc::norm(ek);
    // gc::Vector3 grad_anglei = -(grad_anglek + grad_anglej);

    // chain rule
    gc::Vector3 grad_cosanglek = -sin(anglek) * grad_anglek;
    gc::Vector3 grad_cosanglei = -sin(anglei) * grad_anglei;
    gc::Vector3 grad_cosanglej = -sin(anglej) * grad_anglej;

    // g = qj * ej_perp +  qi * ei_perp +  qk * ek_perp
    // gradient of |g|^2
    return 2 * qj * qj * gc::norm(ej) * grad_ejnorm +
           2 * qk * qk * gc::norm(ek) * grad_eknorm +
           2 * qj * qi * gc::norm(ei) *
               (grad_ejnorm * cos(anglek) + gc::norm(ej) * grad_cosanglek) +
           2 * qi * qk * gc::norm(ei) *
               (grad_eknorm * cos(anglej) + gc::norm(ek) * grad_cosanglej) +
           2 * qj * qk *
               (grad_ejnorm * gc::norm(ek) * cos(anglei) +
                gc::norm(ej) * grad_eknorm * cos(anglei) +
                gc::norm(ej) * gc::norm(ek) * grad_cosanglei);
  }
}
} // namespace solver
} // namespace mem3dg
