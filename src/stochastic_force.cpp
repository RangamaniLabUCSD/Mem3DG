
#include <iostream>
#include <math.h>
#include <random>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include <pcg_random.hpp>

#include "ddgsolver/force.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::stochastic_force(double sigma) {
  gcs::VertexData<size_t> &v_ind = vpg.vertexIndices;

  gcs::EdgeData<double> random_var(mesh);
  std::default_random_engine random_generator;
  std::normal_distribution<double> normal_dist(0, sigma);

  for (gcs::Edge e : mesh.edges()) {
    random_var[e] = normal_dist(random_generator);
  }
  for (gcs::Vertex v : mesh.vertices()) {
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      gc::Vector3 posi_diff_unit = vec_from_halfedge(he, vpg).normalize();
      for (size_t i = 0; i < 3; i++) {
        stochasticForces(v_ind[v], i) +=
            (random_var[he.edge()] * posi_diff_unit)[i];
      }
    }
  }
}