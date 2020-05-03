
#include <assert.h>
#include <iostream>
#include <cmath>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/IterativeLinearSolvers>

#include "ddgsolver/force.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

template <typename T>
void log(gcs::FaceData<T> face_a, gcs::HalfedgeMesh &mesh, std::string name) {
  for (gcs::Face f : mesh.faces()) {
    std::cout << name << face_a[f] << std::endl;
  }
}

void Force::getStretchingForces(double Ksl, double Ksg) {
  const gcs::FaceData<gc::Vector3> &face_n = vpg.faceNormals;
  // log(face_n, mesh,"face normal");

  const gcs::FaceData<double> &face_a = vpg.faceAreas;
  // log(face_a, mesh, "faceArea");

  //gcs::VertexData<size_t> &v_ind = vpg.vertexIndices;

  /*Eigen::Matrix<double, Eigen::Dynamic, 3> local_force;
  local_force.setZero(mesh.nVertices(), 3);*/


  /*Eigen::Matrix<double, Eigen::Dynamic, 3> global_force;
  global_force.setZero(mesh.nVertices(), 3);*/

  double total_area = 0.0;

  for (gcs::Face f : mesh.faces()) {
    total_area += face_a[f];
  }

  for (gcs::Vertex v : mesh.vertices()) {
    gc::Vector3 localForce{ 0.0,0.0,0.0 };
    gc::Vector3 globalForce{ 0.0,0.0,0.0 };

    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      gcs::Halfedge base_he = he.next();
      gc::Vector3 base_vec = -vecFromHalfedge(base_he, vpg);
      // std::cout << "base vector" << base_vec << std::endl;
      gc::Vector3 gradient = gc::cross(base_vec, face_n[he.face()]);
      assert((gc::dot(gradient, vecFromHalfedge(he, vpg))) < 0);
      // std::cout << "gradient" << gradient << std::endl;
      // auto force_v = force.row(v_ind[v]);
      // Eigen::Map<Eigen::Matrix<double, 1, 3>> force_v (&gradient.x, 3);
      localForce += -2 * Ksl * gradient
            * (face_a[base_he.face()] - targetFaceAreas[base_he.face()]) 
            / targetFaceAreas[base_he.face()];
      globalForce += -2 * Ksg * gradient
            * (total_area - targetSurfaceArea) / total_area;
      stretchingForces[v] = localForce + globalForce;
      // force.row(v_ind[v]) << gradient.x, gradient.y, gradient.z;
    }
  }
  
  // std::cout << "force" << force << std::endl;
}