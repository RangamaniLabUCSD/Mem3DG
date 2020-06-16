
#include <cassert>
#include <cmath>
#include <iostream>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include "ddgsolver/force.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getVelocityFromPastPosition(double& timeStep) {
  for (gcs::Vertex v : mesh.vertices()) {
    vertexVelocity[v] = (vpg.inputVertexPositions[v] - pastPositions[v]) / timeStep;
    //std::cout << vertexVelocity[v] << std::endl;
    // std::cout << "I am here" << pastPositions[v] << std::endl;
  }
}

void Force::getDampingForces() {
  // Compute approximate vertex positions
  // TODO: this can be computed on vertex position update and cached to prevent
  // the sequential loop.
  dampingForces.fill({ 0.0,0.0,0.0 });

  for (gcs::Vertex v : mesh.vertices()) {
    for (gcs::Vertex v_adj : v.adjacentVertices()) {
      gc::Vector3 velo_diff = vertexVelocity[v] - vertexVelocity[v_adj];
      gc::Vector3 posi_diff_unit =
          (vpg.inputVertexPositions[v] - vpg.inputVertexPositions[v_adj])
              .normalize();
      dampingForces[v] += -gamma * (gc::dot(velo_diff, posi_diff_unit) * posi_diff_unit);
    }
  }
}
} // end namespace ddgsolver
