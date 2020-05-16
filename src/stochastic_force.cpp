
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

void Force::getStochasticForces(double &sigma) {
  stochasticForces.fill({ 0,0,0 });
  gcs::EdgeData<double> random_var(mesh);
  std::default_random_engine random_generator;
  std::normal_distribution<double> normal_dist(0, sigma);

  for (gcs::Edge e : mesh.edges()) {
    random_var[e] = normal_dist(random_generator);
  }
  for (gcs::Vertex v : mesh.vertices()) {
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      gc::Vector3 posi_diff_unit = -vecFromHalfedge(he, vpg).normalize();
      stochasticForces[v] += random_var[he.edge()] * posi_diff_unit;
    }
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
