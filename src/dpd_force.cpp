
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

void Force::getVelocityFromPastPosition(double dt) {
  auto vel_e = EigenMap<double, 3>(vel);
  auto pos_e = EigenMap<double, 3>(vpg.inputVertexPositions);
  auto pastpos_e = EigenMap<double, 3>(pastPositions);
  vel_e = (pos_e - pastpos_e) / dt;
}

void Force::getDPDForces() {
  // Reset forces to zero
  auto dampingForces_e = EigenMap<double, 3>(dampingForces);
  auto stochasticForces_e = EigenMap<double, 3>(stochasticForces);
  dampingForces_e.setZero();
  stochasticForces_e.setZero();

  // alias positions
  const auto &pos = vpg.inputVertexPositions;

  // std::default_random_engine random_generator;
  // gcs::EdgeData<double> random_var(mesh);
  std::normal_distribution<double> normal_dist(0, sigma);

  for (gcs::Edge e : mesh.edges()) {
    gcs::Halfedge he = e.halfedge();
    gcs::Vertex v1 = he.vertex();
    gcs::Vertex v2 = he.next().vertex();

    gc::Vector3 dVel12 = vel[v1] - vel[v2];
    gc::Vector3 dPos12_n = (pos[v1] - pos[v2]).normalize();

    // gc::Vector3 dVel21 = vel[v2] - vel[v1];
    // gc::Vector3 dPos21_n = (pos[v2] - pos[v1]).normalize();

    // std::cout << -gamma * (gc::dot(dVel12, dPos12_n) * dPos12_n)
    //           << " == " << -gamma * (gc::dot(-dVel12, -dPos12_n) * -dPos12_n)
    //           << " == " << -gamma * (gc::dot(dVel21, dPos21_n) * dPos21_n)
    //           << std::endl;
    
    gc::Vector3 df = gamma * (gc::dot(dVel12, dPos12_n) * dPos12_n);
    dampingForces[v1] -= df;
    dampingForces[v2] += df;

    double noise = normal_dist(rng);
    stochasticForces[v1] += noise * dPos12_n;
    stochasticForces[v2] -= noise * dPos12_n;
  }
}

void Force::pcg_test() {
  // Generate a normal distribution around that mean
  std::normal_distribution<> normal_dist(0, 2);

  // Make a copy of the RNG state to use later
  pcg32 rng_checkpoint = rng;

  // Produce histogram
  std::map<int, int> hist;
  for (int n = 0; n < 10000; ++n) {
    ++hist[std::round(normal_dist(rng))];
  }
  std::cout << "Normal distribution around " << 0 << ":\n";
  for (auto p : hist) {
    std::cout << std::fixed << std::setprecision(1) << std::setw(2) << p.first
              << ' ' << std::string(p.second / 30, '*') << '\n';
  }

  // Produce information about RNG usage
  std::cout << "Required " << (rng - rng_checkpoint) << " random numbers.\n";
}

} // end namespace ddgsolver
