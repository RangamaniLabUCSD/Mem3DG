
#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include "ddgsolver/force.h"
#include "ddgsolver/meshops.h"
#include "ddgsolver/util.h"

namespace ddgsolver {
  namespace gc = ::geometrycentral;
  namespace gcs = ::geometrycentral::surface;

  /// helper function for calculating cotangent
  inline double cotan(double phi) { return (1 / tan(phi)); }

  /// helper function for calculaing d(H at vertex v) / d(x at vertex v)
  gc::Vector3 getdHdx(gcs::VertexPositionGeometry& vpg, gcs::Vertex v, gcs::VertexData<double>& H) {
    // setZero and alias geometry quantities for convenience
    gc::Vector3 dHdx({ 0.0,0.0,0.0 });
    gc::Vector3 dAvdx({ 0.0,0.0,0.0 });
    const gcs::CornerData<double>& cornerAngles = vpg.cornerAngles;
    const gcs::EdgeData<double>& dihedralAngles = vpg.edgeDihedralAngles;
    const gcs::FaceData<gc::Vector3>& face_n = vpg.faceNormals;
    const gcs::FaceData<double>& face_a = vpg.faceAreas;

    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      gcs::Halfedge base_he = he.next();
      gc::Vector3 edgeGradient = -vecFromHalfedge(he, vpg).normalize();
      gc::Vector3 base_vec = vecFromHalfedge(base_he, vpg);
      gc::Vector3 dAdx = -gc::cross(base_vec, face_n[he.face()]);

      double phi1 = cornerAngles[he.twin().next().next().corner()];
      double phi2 = cornerAngles[he.next().next().corner()];
      gc::Vector3 dTdx = (cotan(phi1) * face_n[he.twin().face()]
        + cotan(phi2) * face_n[he.face()]) / vpg.edgeLengths[he.edge()];
      double H_component = vpg.edgeLengths[he.edge()] * dihedralAngles[he.edge()] / (4 * vpg.vertexDualAreas[v]);
      dHdx += (edgeGradient * dihedralAngles[he.edge()] + vpg.edgeLengths[he.edge()] * dTdx)
        / (4 * vpg.vertexDualAreas[v]) - H_component / vpg.vertexDualAreas[v] * (dAdx / 3);
      H[v] += H_component;
      dAvdx += dAdx / 3;
    }
    return dHdx;
  }

  /// helper function for calculaing d(H at the neighbors of vertex v) / d(x at vertex v)


  void Force::getConservativeForces() {

    /// A. PRESSURE FORCES
    pressureForces.fill({ 0.0, 0.0, 0.0 });
    volume = 0;
    double face_volume;
    gcs::FaceData<int> sign_of_volume(mesh);
    for (gcs::Face f : mesh.faces()) {
      face_volume = signedVolumeFromFace(f, vpg);
      volume += face_volume;
      if (face_volume < 0) {
        sign_of_volume[f] = -1;
      }
      else {
        sign_of_volume[f] = 1;
      }
    }
    std::cout << "total volume:  " << volume / maxVolume / Vt << std::endl;

    /// B. STRETCHING FORCES
    stretchingForces.fill({ 0.0,0.0,0.0 });
    const gcs::FaceData<gc::Vector3>& face_n = vpg.faceNormals;
    const gcs::FaceData<double>& face_a = vpg.faceAreas;
    auto faceArea_e = EigenMap(vpg.faceAreas);
    surfaceArea = faceArea_e.sum();
    std::cout << "area: " << surfaceArea / initialSurfaceArea << std::endl;

    /// C. BENDING FORCES
    bendingForces.fill({ 0.0,0.0,0.0 });

    /// D. LOOPING VERTICES
    for (gcs::Vertex v : mesh.vertices()) {
      bendingForces[v] +=

      for (gcs::Halfedge he : v.outgoingHalfedges()) {
        // Pressure forces
        gcs::Halfedge base_he = he.next();
        gc::Vector3 p1 = vpg.inputVertexPositions[base_he.vertex()];
        gc::Vector3 p2 = vpg.inputVertexPositions[base_he.next().vertex()];
        gc::Vector3 dVdx = 0.5 * gc::cross(p1, p2) / 3.0;
        assert(gc::dot(dVdx, vpg.inputVertexPositions[v] - p1) *
          sign_of_volume[he.face()] >
          0);
        pressureForces[v] +=
          -0.5 * Kv * (volume - maxVolume * Vt) / (maxVolume * Vt) * dVdx;

        // Stretching forces
        gc::Vector3 edgeGradient = -vecFromHalfedge(he, vpg).normalize();
        gc::Vector3 base_vec = vecFromHalfedge(base_he, vpg);
        gc::Vector3 dAdx = -gc::cross(base_vec, face_n[he.face()]);
        assert((gc::dot(dAdx, vecFromHalfedge(he, vpg))) < 0);
        if (Ksl != 0) {
          stretchingForces[v] += -2 * Ksl * dAdx *
            (face_a[base_he.face()] - initialFaceAreas[base_he.face()]) /
            initialFaceAreas[base_he.face()];
        }
        if (Ksg != 0) {
          stretchingForces[v] +=
            -2 * Ksg * dAdx * (surfaceArea - initialSurfaceArea) / initialSurfaceArea;
        }
        if (Kse != 0) {
          stretchingForces[v] += -Kse * edgeGradient *
            (vpg.edgeLengths[he.edge()] - targetEdgeLength[he.edge()]) / targetEdgeLength[he.edge()];
        }

    }

  }
} // end namespace ddgsolver
