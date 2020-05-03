
#include <iostream>
#include <cmath>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include "ddgsolver/force.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getDampingForces(double gamma)
{
  // Compute approximate vertex positions
  // TODO: this can be computed on vertex position update and cached to prevent
  // the sequential loop.
  gcs::VertexData<gc::Vector3> velocity(mesh);
  for (gcs::Vertex v : mesh.vertices()) {
    velocity[v] = (vpg.inputVertexPositions[v] - pastPositions[v]) / timestep;
    // std::cout << "I am here" << pastPositions[v] << std::endl;
  }

  for (gcs::Vertex v : mesh.vertices()) {
    for (gcs::Vertex v_adj : v.adjacentVertices()) {
      gc::Vector3 velo_diff = velocity[v] - velocity[v_adj];
      gc::Vector3 posi_diff_unit =
          (vpg.inputVertexPositions[v] - vpg.inputVertexPositions[v_adj])
              .normalize();
      dampingForces[v] += (gc::dot(velo_diff, posi_diff_unit) * posi_diff_unit);
    }
  }
}