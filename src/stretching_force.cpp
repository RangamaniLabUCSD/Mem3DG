
#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/IterativeLinearSolvers>

#include "ddgsolver/force.h"
#include "ddgsolver/util.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

template <typename T>
void log(gcs::FaceData<T> face_a, gcs::HalfedgeMesh &mesh, std::string name) {
  for (gcs::Face f : mesh.faces()) {
    std::cout << name << face_a[f] << std::endl;
  }
}

void Force::getStretchingForces() {
  stretchingForces.fill({ 0.0,0.0,0.0 });
  const gcs::FaceData<gc::Vector3> &face_n = vpg.faceNormals;
  // log(face_n, mesh,"face normal");

  const gcs::FaceData<double> &face_a = vpg.faceAreas;
  // log(face_a, mesh, "faceArea");

  // gcs::VertexData<size_t> &v_ind = vpg.vertexIndices;

  /*Eigen::Matrix<double, Eigen::Dynamic, 3> local_force;
  local_force.setZero(mesh.nVertices(), 3);*/

  /*Eigen::Matrix<double, Eigen::Dynamic, 3> global_force;
  global_force.setZero(mesh.nVertices(), 3);*/


  auto faceArea_e = EigenMap(vpg.faceAreas);
  surfaceArea = faceArea_e.sum();
  std::cout << "area: " << surfaceArea / initialSurfaceArea << std::endl;

  for (gcs::Vertex v : mesh.vertices()) {
    gc::Vector3 localForce{0.0, 0.0, 0.0};
    gc::Vector3 globalForce{0.0, 0.0, 0.0};
    gc::Vector3 edgeForce{ 0.0, 0.0, 0.0 };

    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      gc::Vector3 edgeGradient = -vecFromHalfedge(he, vpg).normalize();

      gcs::Halfedge base_he = he.next();
      gc::Vector3 base_vec = vecFromHalfedge(base_he, vpg);
      gc::Vector3 gradient = -gc::cross(base_vec, face_n[he.face()]);
      assert((gc::dot(gradient, vecFromHalfedge(he, vpg))) < 0);
      
      if(Ksl != 0){
        localForce += -2 * Ksl * gradient *
            (face_a[base_he.face()] - initialFaceAreas[base_he.face()]) /
            initialFaceAreas[base_he.face()];
      }
      
      if (Ksg != 0) {
        globalForce +=
          -2 * Ksg * gradient * (surfaceArea - initialSurfaceArea) / initialSurfaceArea;
      }
      
      if (Kse != 0) {
        edgeForce += -Kse * edgeGradient * 
        (vpg.edgeLengths[he.edge()] - targetEdgeLength[he.edge()]) / targetEdgeLength[he.edge()];
      }
  
    }
    stretchingForces[v] = localForce + globalForce + edgeForce;

  }

}
} // end namespace ddgsolver
