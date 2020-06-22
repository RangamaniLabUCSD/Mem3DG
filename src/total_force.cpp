
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include <pcg_random.hpp>

#include "ddgsolver/force.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::updateForces() {
  std::cerr << "Function updateForces() is incomplete!"  << std::endl;

  // DPD Forces
  auto dampingForces_e = ddgsolver::EigenMap<double, 3>(dampingForces);
  auto stochasticForces_e = ddgsolver::EigenMap<double, 3>(stochasticForces);
  dampingForces_e.setZero();
  stochasticForces_e.setZero();

  // std::default_random_engine random_generator;
  gcs::EdgeData<double> random_var(mesh);
  std::normal_distribution<double> normal_dist(0, sigma);

  for (gcs::Edge e : mesh.edges()) {
    random_var[e] = normal_dist(rng);
  }

  // alias positions
  const auto &pos = vpg.inputVertexPositions;

  for (gcs::Vertex v : mesh.vertices()) {
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      auto v_adj = he.next().vertex();

      gc::Vector3 dVel = vel[v] - vel[v_adj];
      gc::Vector3 dPos_n = (pos[v] - pos[v_adj]).normalize();
      dampingForces[v] += -gamma * (gc::dot(dVel, dPos_n) * dPos_n);
      stochasticForces[v] += random_var[he.edge()] * dPos_n;
    }
  }
}
}
